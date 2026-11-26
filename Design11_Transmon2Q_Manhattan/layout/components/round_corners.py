# -*- coding: utf-8 -*-
r"""In-memory 90-degree corner rounding helper for Design11_Transmon2Q_Manhattan QComponents.

Fillets exact 90-degree vertices with circular arcs (default R = 5 um, auto-scaled on short edges)
while leaving external CPW pin butt joints and non-90-degree Manhattan tapered/45-degree leads untouched.
"""
import numpy as np
from shapely.geometry import MultiPolygon, Point, Polygon, box
from qiskit_metal import draw


def _clean(geom, min_area=1e-10):
    parts = [g for g in getattr(geom, 'geoms', [geom]) if g.area > min_area]
    return parts[0] if len(parts) == 1 else MultiPolygon(parts)


def _round_ring(ring, radius, preserve=None, n_pts=7):
    pts = np.asarray(ring.coords, dtype=float)[:-1]
    if len(pts) < 4:
        return list(ring.coords)
    inc = np.roll(pts, 1, axis=0) - pts
    out = np.roll(pts, -1, axis=0) - pts
    b = np.linalg.norm(inc, axis=1)
    a = np.linalg.norm(out, axis=1)
    valid = (b > 1e-9) & (a > 1e-9)
    inc_u = np.zeros_like(inc)
    out_u = np.zeros_like(out)
    inc_u[valid] = inc[valid] / b[valid, None]
    out_u[valid] = out[valid] / a[valid, None]
    is_right = valid & (np.abs(np.sum(inc_u * out_u, axis=1)) < 1e-6)
    short_edge = np.minimum(b, a)
    setback = np.where(short_edge >= 2.0 * radius - 1e-9, radius, 0.45 * short_edge)
    angles = np.linspace(0.0, np.pi / 2.0, n_pts)
    res = []
    for pt, u, v, r, ok in zip(pts, inc_u, out_u, setback, is_right):
        if ok and r >= 5e-5:
            if preserve is not None and preserve.contains(Point(pt[0], pt[1])):
                res.append(pt)
                continue
            arc = pt + r * ((1.0 - np.sin(angles[:, None])) * u + (1.0 - np.cos(angles[:, None])) * v)
            res.extend(arc)
        else:
            res.append(pt)
    return res


def _round_geom(geom, radius, preserve=None, n_pts=7):
    if radius <= 0 or geom is None or geom.is_empty:
        return geom
    g = geom.simplify(0.0, preserve_topology=True)
    if isinstance(g, Polygon):
        out = Polygon(_round_ring(g.exterior, radius, preserve, n_pts),
                      [_round_ring(r, radius, preserve, n_pts) for r in g.interiors])
        return _clean(out) if out.is_valid else geom
    if isinstance(g, MultiPolygon):
        parts = [_round_geom(p, radius, preserve, n_pts) for p in g.geoms]
        return _clean(draw.union(*parts))
    return geom


def _pin_preserve_zone(comp, pad_mm=0.025):
    boxes = []
    for pin in comp.pins.values():
        x, y = [float(v) for v in pin.middle]
        boxes.append(box(x - pad_mm, y - pad_mm, x + pad_mm, y + pad_mm))
    return draw.union(*boxes) if boxes else None


def round_component_polys(comp, radius=5e-3, skip_names=('squid_ground', 'jj_pad_ground'), n_pts=7):
    """Round 90-degree corners of polygons belonging to ``comp`` except ``skip_names`` and CPW pins."""
    if radius <= 0:
        return
    pres = _pin_preserve_zone(comp)
    tbl = comp.design.qgeometry.tables['poly']
    for idx in tbl.index[tbl['component'] == comp.id]:
        name = str(tbl.at[idx, 'name'])
        if any(s in name for s in skip_names):
            continue
        tbl.at[idx, 'geometry'] = _round_geom(tbl.at[idx, 'geometry'], radius, preserve=pres, n_pts=n_pts)


def square_merged_etch_corners(comp, poly_name, raw_geom, other_etch, pad=5e-3):
    """Undo the rounding of ``poly_name`` (subtract polygon of ``comp``) where it meets ``other_etch``.

    ``raw_geom`` is the polygon before ``round_component_polys``. Inside ``other_etch`` buffered by ``pad`` the
    unrounded shape is restored (union), so two merging etches meet without a thin ground cusp between a fillet
    and the other etch edge. Elsewhere the rounded corners are kept.
    """
    if raw_geom is None or other_etch is None:
        return
    tbl = comp.design.qgeometry.tables['poly']
    idx = tbl.index[(tbl['component'] == comp.id) & (tbl['name'] == poly_name)]
    for i in idx:
        g = tbl.at[i, 'geometry']
        tbl.at[i, 'geometry'] = _clean(g.union(raw_geom.intersection(other_etch.buffer(pad, join_style=2))))

