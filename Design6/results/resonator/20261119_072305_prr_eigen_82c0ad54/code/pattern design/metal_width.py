"""Minimum metal-width rule shared by QComponents, layout scripts and GDS checks.

Rule (skill/chip_drc.md, C12a): no metal narrower than MIN_METAL_WIDTH_UM on the exported mask,
including the ground strips that remain between two etches (pocket / flux-line slot,
pocket / charge-line slot). Junction metal (leads, arrays, SQUID electrodes) is drawn by the
junction process and is checked separately; pass only base-metal geometry here.

Two levels use the same rule:

* QComponents call ``require_strip`` for every ground strip their own options define, so an
  option set that would leave a strip below the minimum raises at ``make()`` time.
* ``thin_metal`` measures arbitrary metal (a component's local ground, a whole chip, a GDS)
  with a morphological opening; ``code/general/check_metal_width.py`` applies it to GDS files.
"""
from shapely import maximum_inscribed_circle
from shapely.geometry import MultiPolygon, Point, Polygon
from shapely.ops import unary_union

MIN_METAL_WIDTH_UM = 4.0


def min_metal_width(units_per_um=1e-3):
    """Rule value in design units (qiskit-metal designs use mm: 1 um = 1e-3)."""
    return MIN_METAL_WIDTH_UM * units_per_um


class MetalWidthError(ValueError):
    """An option set or a mask leaves metal narrower than the rule."""


def require_strip(owner, label, width, minimum):
    """Raise if a designed ground/metal strip (same units as ``minimum``) is below the rule."""
    if width + 1e-9 * max(1.0, abs(minimum)) < minimum:
        raise MetalWidthError(
            f'{owner}: {label} is {width:.6g}, below the minimum metal width {minimum:.6g} '
            f'(same units). Change the options so the strip is at least the minimum.')


def _parts(geom):
    if geom.is_empty:
        return []
    return list(geom.geoms) if hasattr(geom, 'geoms') else [geom]


def _convex_vertices(polygon):
    """Convex (outward) vertices of every ring of a polygon, as Points."""
    out = []
    rings = [(polygon.exterior, True)] + [(r, False) for r in polygon.interiors]
    for ring, outer in rings:
        coords = list(ring.coords)[:-1]
        n = len(coords)
        if n < 3:
            continue
        ccw = ring.is_ccw
        for i in range(n):
            (x0, y0), (x1, y1), (x2, y2) = coords[i - 1], coords[i], coords[(i + 1) % n]
            cross = (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1)
            # metal lies left of a CCW exterior / right of a CCW hole ring
            convex = cross > 0 if ccw == outer else cross < 0
            if convex:
                out.append(Point(x1, y1))
    return out


def thin_metal(metal, min_width, *, tolerance=None, resolution=16):
    """Return [(region, width_estimate)] where ``metal`` is narrower than ``min_width``.

    metal: shapely (Multi)Polygon in any length unit; min_width in the same unit.
    Uses an opening with radius (min_width - tolerance) / 2. The opening also clips every
    convex corner; a residual that lies within one radius of a single convex vertex is such
    a corner (any obtuse or right-angle corner, any rounded corner) and is not reported.
    Acute spikes and every neck/strip narrower than ``min_width`` are reported. Limitation:
    a sliver thinner than 0.2 * min_width that touches a convex corner is treated as that
    corner's residual (zero-area slivers are rule B14, checked separately).
    width_estimate = diameter of the largest circle inside the residual.
    """
    if tolerance is None:
        tolerance = 1e-3 * min_width
    r = (min_width - tolerance) / 2
    metal = unary_union(metal) if isinstance(metal, (list, tuple)) else metal
    if metal.is_empty:
        return []
    opened = metal.buffer(-r, quad_segs=resolution).buffer(r, quad_segs=resolution)
    residual = metal.difference(opened)
    polygons = [p for p in _parts(metal) if isinstance(p, Polygon)]
    found = []
    reach = r * 1.02 + tolerance
    for piece in _parts(residual):
        if piece.area <= (tolerance * tolerance):
            continue
        circle = maximum_inscribed_circle(piece, tolerance / 10)
        width = 2 * circle.length if circle.length > 0 else 0.0
        near = piece.buffer(tolerance)
        corner = False
        for polygon in polygons:
            if not polygon.intersects(near):
                continue
            for v in _convex_vertices(polygon):
                # right-angle corner residual: inscribed diameter 0.343 r -> 0.17 min_width
                if near.contains(v) and (piece.hausdorff_distance(v) <= reach or width < 0.2 * min_width):
                    corner = True
                    break
            if corner:
                break
        if not corner:
            found.append((piece, width))
    return found


def local_ground(window, etches, metal=()):
    """Ground metal inside ``window``: window minus all etches, plus positive metal on top."""
    ground = window.difference(unary_union(list(etches))) if etches else window
    if metal:
        ground = unary_union([ground, *metal])
    return ground


def check_component(component, *, window_margin_um=50.0, ignore_names=()):
    """Raise MetalWidthError if a placed QComponent leaves metal below the rule.

    The component's base metal is its non-subtract polys (junction table excluded) plus the
    ground plane inside its bounding box grown by ``window_margin_um``, with its own subtract
    polys removed. The margin edge itself is excluded from the result (window clipping)."""
    design = component.design
    unit = design.parse_value('1um')
    table = design.qgeometry.tables['poly']
    rows = table[table.component == component.id]
    rows = rows[~rows.name.isin(ignore_names)]
    if rows.empty:
        return []
    etches = [g for g in rows[rows.subtract].geometry]
    metal = [g for g in rows[~rows.subtract].geometry]
    bounds = unary_union(list(rows.geometry)).bounds
    m = window_margin_um * unit
    window = Polygon.from_bounds(bounds[0] - m, bounds[1] - m, bounds[2] + m, bounds[3] + m)
    ground = local_ground(window, etches, metal)
    inner = window.buffer(-2 * MIN_METAL_WIDTH_UM * unit)
    found = [(g, w) for g, w in thin_metal(ground, min_metal_width(unit)) if g.intersects(inner)]
    if found:
        lines = [f'  width {w / unit:.3f} um at x {g.centroid.x / unit:.1f}, y {g.centroid.y / unit:.1f} um'
                 for g, w in found]
        raise MetalWidthError(f'{component.name}: metal below {MIN_METAL_WIDTH_UM} um\n' + '\n'.join(lines))
    return found
