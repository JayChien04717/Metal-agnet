# -*- coding: utf-8 -*-
r"""QGeometry-level keep-out clearance checks (Stage 1 layout DRC).

Collects keep-out polygons from components implementing ``keepouts()`` (markers,
coupons; see ``skill/qcomponents.md`` and ``skill/layout_structure.md``) and
measures clearance in micrometres to every other component's QGeometry polygons
and buffered CPW paths before GDS export.
"""
import math


def collect_keepouts(design, extra_keepouts=None):
    """Return ``[(owner_name, region_name, polygon_mm), ...]`` from ``design``.

    Every component exposing a callable ``keepouts()`` returning
    ``{region_name: shapely.Polygon}`` in design units (mm) is included.
    Optional ``extra_keepouts`` maps ``{label: polygon_mm}`` with owner ``None``.
    """
    items = []
    for comp_name, comp in design.components.items():
        fn = getattr(comp, 'keepouts', None)
        if callable(fn):
            for region_name, poly in fn().items():
                items.append((comp_name, str(region_name), poly))
    if extra_keepouts:
        for label, poly in extra_keepouts.items():
            items.append((None, str(label), poly))
    return items


def _component_geometries(design, ignore_names=()):
    """Map ``comp_name -> [shapely geometries in mm]`` across poly and path tables."""
    ignored = set(ignore_names)
    id_to_name = {c.id: n for n, c in design.components.items() if n not in ignored}
    by_comp = {}
    for table_name in ('poly', 'path'):
        table = design.qgeometry.tables.get(table_name)
        if table is None or len(table) == 0:
            continue
        for row in table.itertuples():
            comp_name = id_to_name.get(row.component)
            if comp_name is None:
                continue
            geom = row.geometry
            if geom is None or geom.is_empty:
                continue
            if table_name == 'path':
                w = getattr(row, 'width', 0.0)
                if w is not None and not math.isnan(float(w)) and float(w) > 0.0:
                    geom = geom.buffer(float(w) / 2.0, cap_style=2)
            by_comp.setdefault(comp_name, []).append(geom)
    return by_comp


def standard_site_marker_keepouts(sites=('NTHU', 'TSRI'), margin_um=50.0):
    """Return ``{label: Polygon_mm}`` for the fixed NTHU TiAu and/or TSRI A12 chip-marker openings.

    Used to enforce >= 100 um clearance from *all* e-beam site markers on every
    Design, even when only one site's marker component is placed on the mask.
    """
    import shapely
    from shapely.affinity import scale
    from shapely.geometry import box

    out = {}
    site_set = {s.upper() for s in sites}
    if 'NTHU' in site_set:
        from nthu_tiau_marker import chip_marker_sites, site_opening
        for name, polys in chip_marker_sites().items():
            sq_um = site_opening(polys, margin=margin_um)
            out[f'NTHU_TiAu:{name}'] = scale(sq_um, 1e-3, 1e-3, origin=(0, 0))
    if 'TSRI' in site_set:
        from tsri_a12_mark import a12_groups
        for name, polys in a12_groups().items():
            x0, y0, x1, y1 = shapely.union_all(polys).bounds
            sq_um = box(x0 - margin_um, y0 - margin_um, x1 + margin_um, y1 + margin_um)
            out[f'TSRI_A12:{name}'] = scale(sq_um, 1e-3, 1e-3, origin=(0, 0))
    return out


def check_keepouts(design, min_um=100.0, ignore_components=(), extra_keepouts=None,
                   include_site_markers=('NTHU', 'TSRI'), raise_on_violation=True):
    """Check clearance from all keep-out regions to every other component in ``design``.

    Parameters
    ----------
    design : qiskit_metal.designs.DesignPlanar
        Populated Qiskit Metal design (coordinates in mm).
    min_um : float
        Minimum allowed clearance in micrometres (default: 100.0 um per
        ``skill/chip_drc.md`` rule B11a and ``skill/fabrication_output.md``).
    ignore_components : iterable of str
        Component names excluded from the clearance check.
    extra_keepouts : dict[str, shapely.Polygon] | None
        Optional extra keep-out polygons in mm not attached to a component.
    include_site_markers : iterable of str
        E-beam site markers always checked in addition to placed components
        (default: ``('NTHU', 'TSRI')``). Pass ``()`` to check only placed components.
    raise_on_violation : bool
        If True (default), raise ``ValueError`` when any clearance is below ``min_um``.

    Returns
    -------
    dict
        ``{'min_um', 'passed', 'worst_um', 'violations', 'clearances'}``.
    """
    merged_extra = {}
    if include_site_markers:
        merged_extra.update(standard_site_marker_keepouts(sites=include_site_markers))
    if extra_keepouts:
        merged_extra.update(extra_keepouts)

    keepouts = collect_keepouts(design, extra_keepouts=merged_extra or None)
    # Deduplicate when a site marker component is already placed in the design:
    # ignore the marker component itself when checking its own site keep-outs.
    marker_owners = {c_name for c_name, c in design.components.items()
                     if type(c).__name__ in ('TSRIA12Mark', 'TiAuChipMarker', 'TiAuWaferMarker')}
    by_comp = _component_geometries(design, ignore_names=set(ignore_components) | marker_owners)

    clearances = {}
    violations = []
    worst_um = None

    for owner, region, poly in keepouts:
        key = f'{owner}:{region}' if owner is not None else region
        best_um = None
        nearest = None
        for comp_name, geoms in by_comp.items():
            if comp_name == owner:
                continue
            for geom in geoms:
                d_um = poly.distance(geom) * 1e3
                if best_um is None or d_um < best_um:
                    best_um = d_um
                    nearest = comp_name
        if best_um is None:
            clearances[key] = {'nearest': None, 'distance_um': None}
            continue
        rounded = round(float(best_um), 2)
        clearances[key] = {'nearest': nearest, 'distance_um': rounded}
        if worst_um is None or best_um < worst_um:
            worst_um = best_um
        if best_um < min_um:
            violations.append({
                'keepout': key,
                'component': nearest,
                'distance_um': rounded,
            })

    report = {
        'min_um': float(min_um),
        'passed': len(violations) == 0,
        'worst_um': None if worst_um is None else round(float(worst_um), 2),
        'violations': violations,
        'clearances': clearances,
    }
    if raise_on_violation and violations:
        summary = ', '.join(
            f"{v['keepout']} vs {v['component']} = {v['distance_um']} um"
            for v in violations
        )
        raise ValueError(f'Keep-out clearance below {min_um:g} um: {summary}')
    return report


def check_gds_keepouts(gds_path, min_um=100.0, sites=('NTHU', 'TSRI'), margin_um=50.0,
                       cheese_max_um=5.2, layer=1, raise_on_violation=False):
    """Check a 10 mm single-die GDS (layer 1 metal) against all site marker openings.

    Ignores small ground-cheese holes (<= ``cheese_max_um``) and the site marker
    openings themselves when measuring clearance to active circuitry / routes / pads.
    """
    from pathlib import Path
    import gdstk
    import shapely
    from shapely.affinity import scale
    from shapely.geometry import Polygon, box

    lib = gdstk.read_gds(str(gds_path))
    top, = lib.top_level()
    metal = shapely.union_all([
        Polygon(p.points).buffer(0) for p in top.get_polygons() if p.layer == layer
    ])
    raw_etch = box(-5000.0, -5000.0, 5000.0, 5000.0).difference(metal)
    active = [
        e for e in shapely.get_parts(raw_etch)
        if e.area > 1.0 and not e.buffer(-cheese_max_um / 2.0).is_empty
    ]

    site_polys_um = {
        k: scale(v, 1e3, 1e3, origin=(0, 0))
        for k, v in standard_site_marker_keepouts(sites=sites, margin_um=margin_um).items()
    }
    # Exclude already-cut marker openings themselves from active circuitry
    all_openings = shapely.union_all(list(site_polys_um.values())) if site_polys_um else Polygon()
    circuitry_parts = [e.difference(all_openings.buffer(0.1)) for e in active]
    circuitry_parts = [p for p in circuitry_parts if not p.is_empty and p.area > 0.5]
    circuitry = shapely.union_all(circuitry_parts) if circuitry_parts else Polygon()

    clearances = {}
    violations = []
    worst_um = None
    for key, sq in site_polys_um.items():
        d_um = round(float(sq.distance(circuitry)), 2) if not circuitry.is_empty else None
        clearances[key] = {'distance_um': d_um}
        if d_um is not None:
            if worst_um is None or d_um < worst_um:
                worst_um = d_um
            if d_um < min_um:
                violations.append({'keepout': key, 'distance_um': d_um})

    report = {
        'gds': str(Path(gds_path)),
        'sites': [s.upper() for s in sites],
        'min_um': float(min_um),
        'passed': len(violations) == 0,
        'worst_um': worst_um,
        'violations': violations,
        'clearances': clearances,
    }
    if raise_on_violation and violations:
        summary = ', '.join(f"{v['keepout']} = {v['distance_um']} um" for v in violations)
        raise ValueError(f'{gds_path}: marker clearance below {min_um:g} um: {summary}')
    return report


def _cli():
    import argparse
    import json
    from pathlib import Path
    import sys

    root = Path(__file__).resolve().parents[2]
    for p in (root / 'code' / 'general', root / 'pattern design' / 'marker'):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))

    ap = argparse.ArgumentParser(description='Check single-die GDS clearance against NTHU TiAu and TSRI A12 markers.')
    ap.add_argument('gds', nargs='+', help='One or more single-die GDS files')
    ap.add_argument('--min-um', type=float, default=100.0, help='Minimum clearance in um (default: 100)')
    ap.add_argument('--sites', default='NTHU,TSRI', help='Comma-separated sites to check (default: NTHU,TSRI)')
    args = ap.parse_args()
    sites = tuple(s.strip() for s in args.sites.split(',') if s.strip())

    failed = False
    results = []
    for path in args.gds:
        rep = check_gds_keepouts(path, min_um=args.min_um, sites=sites, raise_on_violation=False)
        results.append(rep)
        if not rep['passed']:
            failed = True
    print(json.dumps(results if len(results) > 1 else results[0], indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == '__main__':
    _cli()
