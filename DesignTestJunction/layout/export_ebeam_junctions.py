"""Overlay the e-beam junction drawings (tapered hook contacts) on Design3 / Design6 fab masks.

Draft development export. Reads a parent fabrication GDS (never modified),
measures every Dolan hook mouth on its base metal (1/0), places
TransmonJunction (one per SQUID junction, ``squid=False``) and
FluxoniumJunction with ``hook_contact=True`` on layers 11-14, checks the
placement, writes a new GDS + site figure + JSON record into DesignTestJunction.
Coordinates are in um. Run from the repository root with the qmetal Python:

    python -X utf8 DesignTestJunction/layout/export_ebeam_junctions.py --design Design3
    python -X utf8 DesignTestJunction/layout/export_ebeam_junctions.py --design Design6

No solver is run; the base-metal geometry is copied unchanged.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import gdstk
from shapely import set_precision, union_all
from shapely.affinity import scale
from shapely.geometry import Polygon, box

DESIGN = Path(__file__).resolve().parents[1]
ROOT = DESIGN.parent
sys.path[:0] = [str(ROOT / 'pattern design' / 'ebeam_junction')]
from qiskit_metal import designs  # noqa: E402
from fluxonium_junction import FluxoniumJunction  # noqa: E402
from transmon_junction import TransmonJunction  # noqa: E402

BASE = (1, 0)
LAYERS = dict(junction=11, recess=12, array=13, array_bar=14)
REVISION = '2026-11-19_r3'
HOOK_W, HOOK_D, HOOK_DX = 2.0, 4.0, 2.0   # fab notch size and ground-side x offset (my_tmon_fab)

# Parents and nominal sites. SQUID sites: (name, island-hook x, nominal mouth midpoint y);
# the mouths themselves are measured on the parent metal. FX: (name, x, y, mirror).
SITES = {
    'Design3': dict(
        parent='Design3/gds/Design3_2026-11-18_wafer_r11_fab.gds',
        squid=[('TC_J1', -15.0, -1061.0), ('TC_J2', 15.0, -1061.0)],
        fluxonium=[('FX1', -340.0, -800.0, False), ('FX2', 340.0, -800.0, True)]),
    'Design6': dict(
        parent='Design6/gds/Design6_2026-11-19_wafer_prr_r1_txt02_fab.gds',
        squid=[('Q1_J1', -465.0, -203.0), ('Q1_J2', -435.0, -203.0),
               ('Q2_J1', 435.0, -203.0), ('Q2_J2', 465.0, -203.0),
               ('TC_J1', -15.0, -198.0), ('TC_J2', 15.0, -198.0)],
        fluxonium=[]),
}


def local_metal(polygons, window):
    hits = [p for p in polygons
            if box(*p.bounding_box()[0], *p.bounding_box()[1]).intersects(window)]
    x0, y0, x1, y1 = window.bounds
    clipped = gdstk.boolean(hits, gdstk.rectangle((x0, y0), (x1, y1)), 'and', precision=0.001)
    return union_all([Polygon(p.points).buffer(0) for p in clipped])


def notch_open(metal, x, mouth, sign):
    """2 x 4 um notch at x opening from the mouth edge (sign +1 up, -1 down), closed at its end."""
    # The notch end is rounded (radius 1 um), so probe the straight part only.
    inner = box(x - .9, min(mouth + sign * .2, mouth + sign * 3.0), x + .9, max(mouth + sign * .2, mouth + sign * 3.0))
    end = box(x - .4, min(mouth + sign * 4.2, mouth + sign * 4.6), x + .4, max(mouth + sign * 4.2, mouth + sign * 4.6))
    return metal.intersection(inner).area < 1e-3 and metal.intersection(end).area > .3


def measure(base, spec):
    """Return the placements [(name, class, options, info)] from measured hook mouths."""
    rows = []
    for name, x, ym in spec['squid']:
        metal = local_metal(base, box(x - 12, ym - 15, x + 14, ym + 15))
        top = metal.intersection(box(x - 3, ym, x + 3, ym + 12)).bounds[1]
        bottom = metal.intersection(box(x + HOOK_DX - 3, ym - 12, x + HOOK_DX + 3, ym)).bounds[3]
        assert notch_open(metal, x, top, +1), (name, 'island hook not found')
        assert notch_open(metal, x + HOOK_DX, bottom, -1), (name, 'ground hook not found')
        gap, mid = top - bottom, (top + bottom) / 2
        assert 3 < gap < 8, (name, gap)
        options = dict(pos_x=f'{x - TransmonJunction.PAD_CENTRE_X:.4f}um',
                       pos_y=f'{mid - TransmonJunction.HOOK_CENTRE_Y:.4f}um',
                       squid=False, hook_contact=True, hook_gap=f'{gap:.4f}um',
                       hook_width=f'{HOOK_W}um', hook_depth=f'{HOOK_D}um', hook_offset=f'{HOOK_DX}um')
        rows.append((name, TransmonJunction, options,
                     dict(island_hook=[x, round(top, 4)], ground_hook=[x + HOOK_DX, round(bottom, 4)])))
    for name, x, y, mirror in spec['fluxonium']:
        metal = local_metal(base, box(x - 12, y - 30, x + 12, y + 30))
        top = metal.intersection(box(x - 5, y + 5, x + 5, y + 30)).bounds[1]
        bottom = metal.intersection(box(x - 5, y - 30, x + 5, y - 5)).bounds[3]
        assert notch_open(metal, x, top, +1) and notch_open(metal, x, bottom, -1), (name, 'FX hooks')
        assert abs((top - y) - (y - bottom)) < .01, (name, top, bottom)
        options = dict(pos_x=f'{x}um', pos_y=f'{y}um', mirror=mirror, hook_contact=True,
                       hook_mouth=f'{top - y:.4f}um', hook_width=f'{HOOK_W}um', hook_depth=f'{HOOK_D}um')
        rows.append((name, FluxoniumJunction, options,
                     dict(top_hook=[x, round(top, 4)], bottom_hook=[x, round(bottom, 4)])))
    return rows


def build(rows):
    design = designs.DesignPlanar()
    for name, cls, options, _ in rows:
        extra = dict(layer_array=LAYERS['array'], layer_array_bar=LAYERS['array_bar']) \
            if cls is FluxoniumJunction else {}
        c = cls(design, name, options={**options, 'layer': LAYERS['junction'],
                                       'layer_recess': LAYERS['recess'], **extra})
        assert c.status == 'good', (name, c.status)
    table = design.qgeometry.tables['poly']
    names = {c.id: c.name for c in design.components.values()}
    shapes = []
    for row in table.itertuples():
        geom = set_precision(scale(row.geometry, 1000, 1000, origin=(0, 0)), 0.001)
        assert geom.is_valid and not geom.is_empty and geom.area > 0 and geom.geom_type == 'Polygon'
        shapes.append((names[row.component], row.name, int(row.layer), geom))
    return design, shapes


def is_contact(component, name, layer):
    if layer == LAYERS['junction'] and name.split('_')[0] in ('pad', 'lead'):
        return True    # TransmonJunction tapered contacts (neck + taper on the junction layer)
    return (layer == LAYERS['array'] and name.startswith('contact_')) or \
        (layer == LAYERS['recess'] and name in {f'reference_{i}' for i in (5, 6, 7, 8)})


def check(shapes, base, rows):
    """Contacts fill their own hook only; every other shape stays off the base metal."""
    info = {name: i for name, _, _, i in rows}
    report, clearances = {}, {}
    metals = {c: local_metal(base, union_all([g for cc, _, _, g in shapes if cc == c]).envelope.buffer(20))
              for c in info}
    for component, name, layer, geom in shapes:
        metal = metals[component]
        if not is_contact(component, name, layer):
            assert not geom.intersects(metal), (component, name, 'touches base metal')
            clearances[component] = min(clearances.get(component, 1e9), geom.distance(metal))
            continue
        i = info[component]
        if 'island_hook' in i:
            (hx, hy), sign = (i['island_hook'], 1) if name.startswith('pad') else (i['ground_hook'], -1)
        else:
            top = geom.centroid.y > (i['top_hook'][1] + i['bottom_hook'][1]) / 2
            (hx, hy), sign = (i['top_hook'], 1) if top else (i['bottom_hook'], -1)
        # The metal it touches must be on its own side of the junction gap.
        touch = geom.intersection(metal)
        assert not touch.is_empty, (component, name, 'does not reach base metal')
        assert (touch.bounds[1] >= hy - 1e-6) if sign > 0 else (touch.bounds[3] <= hy + 1e-6), \
            (component, name, 'touches the opposite electrode')
        entry = dict(metal_overlap_um2=round(touch.area, 4))
        if layer != LAYERS['recess']:
            notch = box(hx - HOOK_W / 2, min(hy, hy + sign * HOOK_D), hx + HOOK_W / 2, max(hy, hy + sign * HOOK_D))
            uncovered = notch.difference(metal).difference(geom).area
            assert uncovered < 1e-3, (component, name, 'hook not filled', uncovered)
            entry['hook_uncovered_um2'] = round(uncovered, 6)
            gap_side = box(hx - 10, min(hy, hy - sign * 10), hx + 10, max(hy, hy - sign * 10))
            entry['width_at_mouth_um'] = round(geom.intersection(
                box(hx - 10, hy - 1e-3, hx + 10, hy + 1e-3)).bounds[2] - geom.intersection(
                box(hx - 10, hy - 1e-3, hx + 10, hy + 1e-3)).bounds[0], 4)
            entry['outside_hook_clearance_um'] = round(geom.intersection(gap_side).distance(
                metal.difference(box(hx - 10, min(hy, hy + sign * 10), hx + 10, max(hy, hy + sign * 10)))), 4)
        report[f'{component}/{name}'] = entry
    return dict(contacts=report, noncontact_min_clearance_um={k: round(v, 4) for k, v in sorted(clearances.items())})


def plot(shapes, base, rows, title, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = {LAYERS['junction']: 'crimson', LAYERS['recess']: 'goldenrod',
              LAYERS['array']: 'royalblue', LAYERS['array_bar']: 'darkturquoise'}
    windows = []
    for name, cls, options, i in rows:
        if cls is TransmonJunction:
            (x, top), (_, bot) = i['island_hook'], i['ground_hook']
            windows.append((name, (x - 6, bot - 6, x + 8, top + 6)))
        else:
            x, y = float(options['pos_x'][:-2]), float(options['pos_y'][:-2])
            windows.append((name, (x - 12, y - 23, x + 12, y + 23)))
    n = len(windows)
    cols = min(n, 3)
    fig, axes = plt.subplots((n + cols - 1) // cols, cols, figsize=(5 * cols, 5.4 * ((n + cols - 1) // cols)),
                             squeeze=False)
    for ax in axes.flat[n:]:
        ax.axis('off')
    for ax, (name, w) in zip(axes.flat, windows):
        W = box(*w)
        metal = local_metal(base, W)
        for g in getattr(metal, 'geoms', [metal]):
            if g.is_empty:
                continue
            ax.fill(*g.exterior.xy, color='0.72', lw=0)
            for h in g.interiors:
                ax.fill(*h.xy, color='white', lw=0)
        for _, _, layer, geom in shapes:
            if geom.intersects(W):
                ax.fill(*geom.exterior.xy, color=colors[layer], alpha=.5, lw=.5, ec=colors[layer])
        ax.set_xlim(w[0], w[2]); ax.set_ylim(w[1], w[3]); ax.set_aspect('equal')
        ax.set_title(name, fontsize=10); ax.grid(alpha=.25)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, alpha=a) for c, a in
               (('0.72', 1),) + tuple((colors[l], .5) for l in LAYERS.values())]
    fig.legend(handles, ['base metal 1/0'] + [f'{k} {v}/0' for k, v in LAYERS.items()],
               loc='lower center', ncol=5)
    fig.suptitle(title)
    fig.tight_layout(rect=(0, .04, 1, .96)); fig.savefig(path, dpi=110); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--design', choices=sorted(SITES), required=True)
    args = parser.parse_args()
    spec = SITES[args.design]
    parent = ROOT / spec['parent']
    stem = f'{DESIGN.name}_{REVISION}_{args.design.lower()}_ebeamjj'
    gds, fig = DESIGN / 'gds' / f'{stem}_fab.gds', DESIGN / 'fig' / f'{stem}_sites.png'
    record_path = DESIGN / 'parameters' / f'{stem}_record.json'
    for target in (gds, fig, record_path):
        if target.exists():
            raise FileExistsError(f'Refusing to overwrite {target}; choose a new REVISION')

    library = gdstk.read_gds(str(parent))
    assert (library.unit, library.precision) == (1e-6, 1e-9)
    parent_top, = library.top_level()
    base = parent_top.get_polygons()
    assert {(p.layer, p.datatype) for p in base} == {BASE}, 'parent must hold only 1/0'
    rows = measure(base, spec)
    design, shapes = build(rows)
    check(shapes, base, rows)

    out = gdstk.Library(unit=1e-6, precision=1e-9)
    cell = out.new_cell(parent_top.name)
    cell.add(*base)
    for _, _, layer, geom in shapes:
        cell.add(gdstk.Polygon(list(geom.exterior.coords)[:-1], layer=layer))
    out.write_gds(str(gds), max_points=8190)

    # Read back: base metal identical, overlay matches the components, checks rerun.
    top, = gdstk.read_gds(str(gds)).top_level()
    polys = top.get_polygons()
    new_base = [p for p in polys if (p.layer, p.datatype) == BASE]
    assert len(new_base) == len(base) and all(
        a.points.shape == b.points.shape and (a.points == b.points).all() for a, b in zip(base, new_base)), \
        'base metal changed'
    overlay = [p for p in polys if p.layer in LAYERS.values()]
    assert len(overlay) == len(shapes)
    readback = []
    for (component, name, layer, expected), actual in zip(shapes, overlay):
        geom = Polygon(actual.points)
        assert (actual.layer, actual.datatype) == (layer, 0) and geom.hausdorff_distance(expected) < 1e-6
        readback.append((component, name, layer, geom))
    report = check(readback, new_base, rows)
    plot(readback, new_base, rows, f'{args.design} fab + e-beam junctions, hook_contact=True (draft, um)', fig)

    record = dict(
        status='draft development export; junction drawings not fabrication or EM validated',
        design=args.design, revision=REVISION,
        parent_fabrication_gds=spec['parent'], parent_sha256=sha256(parent.read_bytes()).hexdigest(),
        fabrication_gds=gds.relative_to(DESIGN).as_posix(), fabrication_gds_sha256=sha256(gds.read_bytes()).hexdigest(),
        figure=fig.relative_to(DESIGN).as_posix(),
        layer_map={**{k: f'{v}/0' for k, v in LAYERS.items()}, 'base_metal': '1/0'},
        contact_layers={'TransmonJunction': f'{LAYERS["junction"]}/0 (neck + taper merged)',
                        'FluxoniumJunction': f'{LAYERS["array"]}/0 contacts, {LAYERS["recess"]}/0 flanking strips'},
        measured_hooks={name: i for name, _, _, i in rows},
        placements=[dict(name=n, component=c.__name__, options=o) for n, c, o, _ in rows],
        parameter_snapshot={name: {k: str(v) for k, v in c.options.items()} for name, c in design.components.items()},
        polygons_per_layer={f'{l}/0': sum(1 for p in overlay if p.layer == l) for l in LAYERS.values()},
        checks=report,
        component_sha256={f: sha256((ROOT / 'pattern design/ebeam_junction' / f).read_bytes()).hexdigest()
                          for f in ('fluxonium_junction.py', 'transmon_junction.py')},
        export_script='layout/export_ebeam_junctions.py',
        export_script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        base_metal_check='exact polygon count, order and vertex equality after GDS readback',
        junction_areas='all junctions use the component defaults (junction_width 0.25um, finger_width 0.1um); '
                       'asymmetric-SQUID or EJ-specific areas are not mapped',
        em_rerun=False)
    record_path.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(gds=record['fabrication_gds'], figure=record['figure'], checks=report), indent=1))


if __name__ == '__main__':
    main()
