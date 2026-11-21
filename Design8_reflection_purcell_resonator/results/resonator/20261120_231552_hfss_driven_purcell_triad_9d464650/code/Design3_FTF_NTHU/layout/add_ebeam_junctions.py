"""Place the e-beam junction drawings on a Design3_FTF_NTHU wafer fab mask (draft).

Reads an existing fab GDS (never modified), adds FluxoniumJunction at FX1 / FX2
and two single TransmonJunction drawings at the TC SQUID hooks on separate
layers, checks the placement against the base metal (layer 1/0) and writes a
new GDS plus a site figure. Coordinates are in um. Run with the qmetal Python:

    python -X utf8 Design3_FTF_NTHU/layout/add_ebeam_junctions.py

The base-metal geometry is copied unchanged; no solver is run.
"""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import sys

import gdstk
from shapely import set_precision, union_all
from shapely.affinity import scale
from shapely.geometry import LineString, Polygon, box

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'pattern design'), str(ROOT / 'pattern design' / 'archive')]
from qiskit_metal import designs  # noqa: E402
from fluxonium_junction_layout3 import FluxoniumJunction  # noqa: E402
from transmon_junction_layout3 import TransmonJunction  # noqa: E402

PARENT = ROOT / 'Design3_FTF_NTHU/gds/Design3_2026-11-18_wafer_r11_fab.gds'
OUTPUT = ROOT / 'Design3_FTF_NTHU/gds/Design3_2026-11-19_wafer_r11_ebeamjj_fab.gds'
FIGURE = ROOT / 'Design3_FTF_NTHU/fig/Design3_2026-11-19_wafer_r11_ebeamjj_fab_sites.png'
BASE = (1, 0)
# Junction drawing layers in the output (source drawing layers 1-4 -> 11-14).
LAYERS = dict(junction=11, recess=12, array=13, array_bar=14)
OVERLAP = 1.0   # um past each hook bottom, as the source contacts overlap its base metal


def record_path(path):
    """Use Design-relative paths; allow scratch exports on another Windows drive."""
    try:
        return os.path.relpath(path, ROOT / 'Design3_FTF_NTHU').replace('\\', '/')
    except ValueError:
        return path.resolve().as_posix()

# Hook geometry measured on the r11 fab mask (um). FX stems end at +/-15 um
# from the site centre and the 2 um hooks reach +/-19 um. The TC island-side
# hooks (top) reach y = -1054.5 at x = -15 / 15; the ground-side hooks
# (bottom) are offset by +2 um in x and reach y = -1067.5.
FX_SITES = {'FX1': dict(x=-340.0, y=-800.0, mirror=False),
            'FX2': dict(x=340.0, y=-800.0, mirror=True)}
FX_HOOK_BOTTOM = 19.0
TC_SITES = {'TC_L': -15.0, 'TC_R': 15.0}
TC_HOOK_TOP_Y, TC_HOOK_BOTTOM_Y = -1054.5, -1067.5
TC_PAD_X, TC_LEAD_X = -4.6665, -2.5585          # TransmonJunction local pin x
TC_PAD_END, TC_LEAD_END = 3.0, -3.183            # local ends at contact_extension = 0


def placements():
    """Return [(name, class, options)] in um for every junction drawing."""
    fx_ext = FX_HOOK_BOTTOM + OVERLAP - 16.0
    rows = [(name, FluxoniumJunction, dict(pos_x=f'{s["x"]}um', pos_y=f'{s["y"]}um',
                                            mirror=s['mirror'], contact_extension=f'{fx_ext}um'))
            for name, s in FX_SITES.items()]
    top, bottom = TC_HOOK_TOP_Y + OVERLAP, TC_HOOK_BOTTOM_Y - OVERLAP
    tc_ext = ((top - bottom) - (TC_PAD_END - TC_LEAD_END)) / 2
    oy = top - TC_PAD_END - tc_ext
    for name, hook_x in TC_SITES.items():
        rows.append((name, TransmonJunction, dict(pos_x=f'{hook_x - TC_PAD_X:.4f}um', pos_y=f'{oy:.4f}um',
                                                  squid=False, contact_extension=f'{tc_ext:.4f}um')))
    return rows


def build():
    design = designs.DesignPlanar()
    layer_options = dict(layer=LAYERS['junction'], layer_recess=LAYERS['recess'])
    for name, cls, options in placements():
        extra = dict(layer_array=LAYERS['array'], layer_array_bar=LAYERS['array_bar']) \
            if cls is FluxoniumJunction else {}
        c = cls(design, name, options={**options, **layer_options, **extra})
        assert c.status == 'good', (name, c.status)
    table = design.qgeometry.tables['poly']
    names = {c.id: c.name for c in design.components.values()}
    shapes = []
    for row in table.itertuples():
        component = names[row.component]
        # Snap once before checks/export so figures and checks describe the
        # actual 1 nm GDS grid, including TC's nominal half-nm placements.
        geom = set_precision(scale(row.geometry, 1000, 1000, origin=(0, 0)), 0.001)
        assert geom.is_valid and not geom.is_empty and geom.area > 0
        shapes.append((component, row.name, int(row.layer), geom))
    return design, shapes


def base_metal(polygons, window):
    # Clip the perforated chip in GDS first; avoid repeatedly buffering the
    # complete ground polygon (including all remote perforations) in Shapely.
    hits = [p for p in polygons
            if box(*p.bounding_box()[0], *p.bounding_box()[1]).intersects(window)]
    x0, y0, x1, y1 = window.bounds
    clipped = gdstk.boolean(hits, gdstk.rectangle((x0, y0), (x1, y1)), 'and', precision=0.001)
    return union_all([Polygon(p.points).buffer(0) for p in clipped])


def check(shapes, base):
    """Contacts overlap base metal by about OVERLAP; everything else stays in the etch."""
    report = {}
    clearances = {}
    noncontact_count = 0
    local_metal = {component: base_metal(base, union_all(
        [g for c, _, _, g in shapes if c == component]).envelope.buffer(30))
        for component in {c for c, _, _, _ in shapes}}
    for component, name, layer, geom in shapes:
        metal = local_metal[component]
        touch = geom.intersection(metal)
        # Contacts (and, for FX, the contact-side recess strips 5-8, which also
        # overlap the base metal in the source) may touch metal at their outer end.
        is_contact = (component.startswith('FX') and (
            name.startswith('contact_') or (layer == LAYERS['recess'] and name in
                                            {f'reference_{i}' for i in (5, 6, 7, 8)}))) or \
            (component.startswith('TC') and name.split('_')[0] in ('pad', 'lead'))
        if not is_contact:
            assert not geom.intersects(metal), (component, name, touch.area, 'touches base metal')
            clearances[component] = min(clearances.get(component, float('inf')), geom.distance(metal))
            noncontact_count += 1
            continue
        gx0, gy0, gx1, gy1 = geom.bounds
        site_y = FX_SITES[component]['y'] if component.startswith('FX') else \
            (TC_HOOK_TOP_Y + TC_HOOK_BOTTOM_Y) / 2
        top = gy1 - site_y > site_y - gy0
        # Overlap along the shape's own centre line (the mask rounds the hook
        # bottom, so it is slightly longer off the hook centre).
        cx = (gx0 + gx1) / 2
        line = metal.intersection(LineString([(cx, gy0), (cx, gy1)]))
        length = line.length
        assert 0.85 * OVERLAP < length < OVERLAP + 0.4, (component, name, length)
        lo, hi = line.bounds[1], line.bounds[3]
        assert abs(hi - gy1) < 1e-6 if top else abs(lo - gy0) < 1e-6, (component, name, 'not at outer end')
        # Anything else touching metal is a side sliver along the hook wall.
        end = box(gx0 - 1, gy1 - OVERLAP - 1, gx1 + 1, gy1 + 1) if top else \
            box(gx0 - 1, gy0 - 1, gx1 + 1, gy0 + OVERLAP + 1)
        side = touch.difference(end)
        side_width = 0.0 if side.is_empty else max(
            p.bounds[2] - p.bounds[0] for p in getattr(side, 'geoms', [side]))
        side_allowed = component.startswith('FX') and layer == LAYERS['recess'] and \
            name in {f'reference_{i}' for i in (5, 6, 7, 8)}
        assert side.is_empty or (side_allowed and side_width <= 0.005 + 1e-9), \
            (component, name, 'cuts into hook wall', side_width)
        report[f'{component}/{name}'] = dict(
            centreline_overlap_um=round(length, 3), overlap_area_um2=round(touch.area, 4),
            side_sliver_width_um=round(side_width, 4),
            clearance_to_hook_wall_um=round(geom.difference(end).distance(metal.difference(end)), 3)
            if side.is_empty else 0.0)
    return dict(contact_overlaps=report, noncontact_shapes_checked=noncontact_count,
                noncontact_min_clearance_um={k: round(v, 6) for k, v in sorted(clearances.items())})


def plot(shapes, base, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    windows = {'FX1 (array to -x)': (-410, -825, -325, -775), 'FX1 junction': (-344, -806, -330, -790),
               'FX2 (array to +x, mirrored)': (325, -825, 410, -775), 'FX2 junction': (330, -806, 344, -790),
               'TC SQUID (single TransmonJunction x2)': (-25, -1072, 25, -1050),
               'TC_L junction': (-19, -1069, -9, -1053)}
    colors = {LAYERS['junction']: 'crimson', LAYERS['recess']: 'purple',
              LAYERS['array']: 'seagreen', LAYERS['array_bar']: 'orange'}
    fig, axes = plt.subplots(3, 2, figsize=(13, 15))
    for ax, (title, w) in zip(axes.flat, windows.items()):
        W = box(*w)
        metal = base_metal(base, W)
        for g in getattr(metal, 'geoms', [metal]):
            if g.is_empty:
                continue
            ax.fill(*g.exterior.xy, color='0.75', lw=0)
            for h in g.interiors:
                ax.fill(*h.xy, color='white', lw=0)
        for _, _, layer, geom in shapes:
            if geom.intersects(W):
                ax.fill(*geom.exterior.xy, color=colors[layer], alpha=.45, lw=.6, ec=colors[layer])
        ax.set_xlim(w[0], w[2]); ax.set_ylim(w[1], w[3]); ax.set_aspect('equal')
        ax.set_title(title, fontsize=10); ax.grid(alpha=.25)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, alpha=a) for c, a in
               (('0.75', 1), ('crimson', .45), ('purple', .45), ('seagreen', .45), ('orange', .45))]
    fig.legend(handles, ['base metal 1/0 (r11 fab)', f'junction {LAYERS["junction"]}/0',
                         f'recess {LAYERS["recess"]}/0', f'array {LAYERS["array"]}/0',
                         f'array bar {LAYERS["array_bar"]}/0'], loc='lower center', ncol=5)
    fig.suptitle('Design3_FTF_NTHU r11 fab + e-beam junction drawings (draft GDS check, um)')
    fig.tight_layout(rect=(0, .03, 1, .97)); fig.savefig(path, dpi=110); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--figure', type=Path, default=FIGURE)
    parser.add_argument('--record', type=Path, help='Optional JSON export record; refuses overwrite')
    args = parser.parse_args()
    targets = [args.output, args.figure] + ([args.record] if args.record else [])
    assert len({p.resolve() for p in targets}) == len(targets), 'output paths must differ'
    assert PARENT.resolve() not in {p.resolve() for p in targets}, 'parent is immutable'
    for target in targets:
        if target.exists():
            raise FileExistsError(f'Refusing to overwrite {target}; choose a new name')
    library = gdstk.read_gds(str(PARENT))
    assert (library.unit, library.precision) == (1e-6, 1e-9)
    parent_top, = library.top_level()
    base = parent_top.get_polygons()
    assert {(p.layer, p.datatype) for p in base} == {BASE}, 'parent must hold only 1/0'
    design, shapes = build()
    report = check(shapes, base)

    out = gdstk.Library(unit=1e-6, precision=1e-9)
    cell = out.new_cell(parent_top.name)
    cell.add(*base)
    for _, _, layer, geom in shapes:
        cell.add(gdstk.Polygon(list(geom.exterior.coords)[:-1], layer=layer))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Preserve the parent's existing polygon partition (up to 8184 vertices).
    out.write_gds(str(args.output), max_points=8190)

    # Read back: base metal identical, junction layers match the components.
    top, = gdstk.read_gds(str(args.output)).top_level()
    polys = top.get_polygons()
    new_base = [p for p in polys if (p.layer, p.datatype) == BASE]
    assert len(new_base) == len(base)
    # Identical polygon partitions and coordinates establish zero XOR exactly,
    # without a costly global union of the chip's ground and perforations.
    assert all(a.points.shape == b.points.shape and (a.points == b.points).all()
               for a, b in zip(base, new_base)), 'base metal changed'
    overlay_xor = {}
    for layer in LAYERS.values():
        # Compare at the file's 1 nm grid. TC placement includes half-nm
        # coordinates; serialization rounds them, rather than preserving floats.
        expected_polys = [gdstk.Polygon(list(g.exterior.coords)[:-1], layer=layer)
                          for _, _, lay, g in shapes if lay == layer]
        actual_polys = [p for p in polys if (p.layer, p.datatype) == (layer, 0)]
        difference = gdstk.boolean(expected_polys, actual_polys, 'xor', precision=0.001)
        overlay_xor[f'{layer}/0'] = sum(abs(p.area()) for p in difference)
        assert overlay_xor[f'{layer}/0'] < 1e-6, (layer, overlay_xor)
        expected = union_all([g for _, _, lay, g in shapes if lay == layer])
        got = union_all([Polygon(p.points) for p in polys if (p.layer, p.datatype) == (layer, 0)])
        assert expected.hausdorff_distance(got) <= 0.001, layer
    counts = {f'{l}/0': sum(1 for p in polys if (p.layer, p.datatype) == (l, 0)) for l in LAYERS.values()}
    readback_overlay = [p for p in polys if p.layer in LAYERS.values()]
    assert len(readback_overlay) == len(shapes)
    readback_shapes = []
    for (component, name, layer, expected), actual in zip(shapes, readback_overlay):
        assert (actual.layer, actual.datatype) == (layer, 0)
        geom = Polygon(actual.points)
        assert geom.is_valid and not geom.is_empty and geom.area > 0
        assert geom.hausdorff_distance(expected) < 1e-9
        readback_shapes.append((component, name, layer, geom))
    report = check(readback_shapes, new_base)

    args.figure.parent.mkdir(parents=True, exist_ok=True)
    plot(readback_shapes, new_base, args.figure)
    record = dict(
        status='draft; junction drawings not fabrication or EM validated',
        parent_fabrication_gds=str(PARENT.relative_to(ROOT / 'Design3_FTF_NTHU')).replace('\\', '/'),
        parent_sha256=sha256(PARENT.read_bytes()).hexdigest(),
        fabrication_gds=record_path(args.output),
        fabrication_gds_sha256=sha256(args.output.read_bytes()).hexdigest(),
        figure=record_path(args.figure),
        layer_map={**{k: f'{v}/0' for k, v in LAYERS.items()}, 'base_metal': '1/0'},
        contact_overlap_um=OVERLAP,
        placements=[dict(name=n, component=c.__name__, options=o) for n, c, o in placements()],
        parameter_snapshot={name: dict(c.options) for name, c in design.components.items()},
        polygons_per_layer=counts,
        checks=report,
        component_sha256={f: sha256((ROOT / 'pattern design/archive' / f).read_bytes()).hexdigest()
                          for f in ('fluxonium_junction_layout3.py', 'transmon_junction_layout3.py')},
        export_script='layout/add_ebeam_junctions.py',
        export_script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        base_metal_xor_um2=0.0,
        base_metal_check='exact polygon count, order and vertex equality after GDS readback',
        junction_layer_xor_um2=overlay_xor,
        gds_grid_nm=1,
        grid_policy='component polygons snapped to 1 nm before placement checks and export; base unchanged',
        em_rerun=False)
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(record, indent=1))


if __name__ == '__main__':
    main()
