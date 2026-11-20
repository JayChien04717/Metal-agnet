"""Design3 r11 fab mask + NCKU-style e-beam junctions (DesignTestJunction r4, draft).

Reads the Design3 r11 fabrication GDS (never modified), measures the Dolan
hook mouths on its base metal (1/0), places the NCKU EBL standard components
``FluxoniumJunction`` (FX1 side='left', FX2 side='right') and
``TransmonJunction`` (TC SQUID, JJ2 = JJ1 translated by 30 um), checks the
placement and writes a new GDS, a site figure and a JSON record into
DesignTestJunction.

Junction sizes come from code/qubit/jj_resistance_model.py (NCKU W80):
FX small JJ EJ/h = 4 GHz (Design3 input); TC EJ1 = EJ2 = 11.528 GHz, the
symmetric NCKU Design3_3 assignment chosen by the user on 2026-11-20 (Design3's
own CZ model uses 16.401 / 6.655 GHz).

Layers: base metal stays 1/0 with unchanged coordinates. The NCKU junction
layers 1/2/3/4 are written as 11/12/13/14 (small-JJ main, shadow, array main,
array shadow), so the repository's 1/0 figure tools still apply. NCKU's own
files put the base on 5/0 and the junctions on 1-4.

Run from the repository root with the qmetal Python:

    python -X utf8 DesignTestJunction/layout/export_ncku_ebeam_junctions.py

No solver is run. The script refuses to overwrite its outputs; bump REVISION.
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
sys.path[:0] = [str(ROOT / 'pattern design' / 'ebeam_junction'), str(ROOT / 'code' / 'qubit')]
from qiskit_metal import designs  # noqa: E402
from fluxonium_junction import FluxoniumJunction  # noqa: E402
from transmon_junction import TransmonJunction  # noqa: E402
import jj_resistance_model as jjm  # noqa: E402

REVISION = '2026-11-20_r4'
PARENT = 'Design3/gds/Design3_2026-11-18_wafer_r11_fab.gds'
BASE = (1, 0)
# NCKU layer -> output layer
LAYERS = dict(jj=11, jj_shadow=12, array=13, array_shadow=14)
HOOK_W, HOOK_D = 2.0, 4.0           # Dolan notch in the base metal (my_tmon_fab / FX stems)
TC_HOOK_DX = 2.0                    # ground hook x offset from the island hook (Design3 r11)
TC_SPACING = 30.0

EJ = dict(FX=4.0, TC_JJ1=11.528, TC_JJ2=11.528)     # GHz per JJ
SITES = dict(TC=(0.0, -1061.0), FX1=(-340.0, -800.0, 'left'), FX2=(340.0, -800.0, 'right'))


def local_metal(polygons, window):
    hits = [p for p in polygons
            if box(*p.bounding_box()[0], *p.bounding_box()[1]).intersects(window)]
    x0, y0, x1, y1 = window.bounds
    clipped = gdstk.boolean(hits, gdstk.rectangle((x0, y0), (x1, y1)), 'and', precision=0.001)
    return union_all([Polygon(p.points).buffer(0) for p in clipped])


def notch_open(metal, x, mouth, sign):
    """2 x 4 um notch at x opening from the mouth (sign +1 up, -1 down), closed at its end."""
    lo, hi = sorted((mouth + sign * .2, mouth + sign * 3.0))
    inner = box(x - .9, lo, x + .9, hi)
    lo, hi = sorted((mouth + sign * 4.2, mouth + sign * 4.6))
    end = box(x - .4, lo, x + .4, hi)
    return metal.intersection(inner).area < 1e-3 and metal.intersection(end).area > .3


def measure(base):
    """Measured hooks and the component placements [(name, class, options, hooks)]."""
    rows, hooks = [], {}
    sx, sy = SITES['TC']
    for tag, x in (('jj1', sx - TC_SPACING / 2), ('jj2', sx + TC_SPACING / 2)):
        metal = local_metal(base, box(x - 12, sy - 15, x + 14, sy + 15))
        top = metal.intersection(box(x - 3, sy, x + 3, sy + 12)).bounds[1]
        bottom = metal.intersection(box(x + TC_HOOK_DX - 3, sy - 12, x + TC_HOOK_DX + 3, sy)).bounds[3]
        assert notch_open(metal, x, top, +1), (tag, 'island hook not found')
        assert notch_open(metal, x + TC_HOOK_DX, bottom, -1), (tag, 'ground hook not found')
        # NCKU frame: island mouth at +2.5, ground mouth at -2.5 from the origin.
        assert abs(top - (sy + 2.5)) < 1e-6 and abs(bottom - (sy - 2.5)) < 1e-6, (tag, top, bottom)
        hooks[f'TC_{tag}'] = dict(island=[x, round(top, 4), +1], ground=[x + TC_HOOK_DX, round(bottom, 4), -1])
    tc = dict(pos_x=f'{sx}um', pos_y=f'{sy}um', squid_spacing=f'{TC_SPACING}um',
              upper_contact_body_width='4um', lower_contact_offset=f'{TC_HOOK_DX}um',
              lower_contact_lead_offset='1um',
              **jjm.squid_options(EJ['TC_JJ1'], EJ['TC_JJ2']))
    rows.append(('TC', TransmonJunction, tc, {k: v for k, v in hooks.items()}))
    fx = jjm.size_from_ej(EJ['FX']).options()
    for name in ('FX1', 'FX2'):
        x, y, side = SITES[name]
        metal = local_metal(base, box(x - 12, y - 30, x + 12, y + 30))
        top = metal.intersection(box(x - 5, y + 5, x + 5, y + 30)).bounds[1]
        bottom = metal.intersection(box(x - 5, y - 30, x + 5, y - 5)).bounds[3]
        assert notch_open(metal, x, top, +1) and notch_open(metal, x, bottom, -1), (name, 'FX hooks')
        assert abs((top - y) - (y - bottom)) < .01, (name, top, bottom)
        end = top - y + HOOK_D + 1.0          # body ends 1 um past the hook bottom (NCKU: 20)
        options = dict(pos_x=f'{x}um', pos_y=f'{y}um', side=side, array_contact_body_width='4um',
                       array_contact_end=f'{end:.4f}um', **fx)
        rows.append((name, FluxoniumJunction, options,
                     {f'{name}_top': dict(island=[x, round(top, 4), +1]),
                      f'{name}_bottom': dict(ground=[x, round(bottom, 4), -1])}))
    return rows


def build(rows):
    design = designs.DesignPlanar()
    layer_options = {TransmonJunction: dict(layer_jj=LAYERS['jj'], layer_jj_shadow=LAYERS['jj_shadow']),
                     FluxoniumJunction: dict(layer_jj=LAYERS['jj'], layer_jj_shadow=LAYERS['jj_shadow'],
                                             layer_array=LAYERS['array'],
                                             layer_array_shadow=LAYERS['array_shadow'])}
    for name, cls, options, _ in rows:
        c = cls(design, name, options={**options, **layer_options[cls]})
        assert c.status == 'good', (name, c.status)
    names = {c.id: c.name for c in design.components.values()}
    shapes = []
    for row in design.qgeometry.tables['poly'].itertuples():
        geom = set_precision(scale(row.geometry, 1000, 1000, origin=(0, 0)), 0.001)
        assert geom.is_valid and not geom.is_empty and geom.area > 0 and geom.geom_type == 'Polygon'
        shapes.append((names[row.component], row.name, int(row.layer), geom))
    return design, shapes


def contact_side(component, name, layer):
    """+1 / -1 for shapes allowed to reach into a hook (island / ground side), else None."""
    if component == 'TC' and layer == LAYERS['jj']:
        if name.startswith('upper_contact'):
            return +1
        if name.startswith('lower_contact'):
            return -1
    if component.startswith('FX'):
        if layer == LAYERS['array'] and name in ('contact_top', 'contact_bottom'):
            return +1 if name.endswith('top') else -1
        if layer == LAYERS['jj_shadow'] and name.startswith('contact_shadow_'):
            return +1 if '_top_' in name else -1
    return None


def check(shapes, base, rows):
    """Hooks filled by their own contact only; other shapes off the base metal; no main overlaps."""
    hooks = {c: h for c, _, _, h in rows}
    metals = {c: local_metal(base, union_all([g for cc, _, _, g in shapes if cc == c]).envelope.buffer(20))
              for c in hooks}
    contacts, clearance = {}, {}
    for component, name, layer, geom in shapes:
        metal = metals[component]
        side = contact_side(component, name, layer)
        if side is None:
            assert not geom.intersects(metal), (component, name, 'touches base metal')
            clearance[component] = min(clearance.get(component, 1e9), geom.distance(metal))
            continue
        # The hook this contact serves: the nearest hook on its side.
        cands = [h[k] for h in hooks[component].values() for k in ('island', 'ground')
                 if k in h and h[k][2] == side]
        hx, hy, _ = min(cands, key=lambda h: abs(h[0] - geom.centroid.x))
        touch = geom.intersection(metal)
        assert not touch.is_empty, (component, name, 'does not reach base metal')
        ok = touch.bounds[1] >= hy - 1e-6 if side > 0 else touch.bounds[3] <= hy + 1e-6
        assert ok, (component, name, 'touches the opposite electrode')
        entry = dict(metal_overlap_um2=round(touch.area, 4))
        if layer != LAYERS['jj_shadow']:
            lo, hi = sorted((hy, hy + side * HOOK_D))
            notch = box(hx - HOOK_W / 2, lo, hx + HOOK_W / 2, hi)
            uncovered = notch.difference(metal).difference(geom).area
            assert uncovered < 1e-3, (component, name, 'hook not filled', uncovered)
            entry.update(hook=[hx, hy], hook_uncovered_um2=round(uncovered, 6),
                         past_hook_bottom_um=round(side * ((geom.bounds[3] if side > 0 else geom.bounds[1])
                                                           - (hy + side * HOOK_D)), 4))
        contacts[f'{component}/{name}'] = entry

    def layer_union(component, layer):
        return union_all([g for c, _, l, g in shapes if c == component and l == layer])
    overlaps = {}
    for component in hooks:
        main_jj, shadow = layer_union(component, LAYERS['jj']), layer_union(component, LAYERS['jj_shadow'])
        main_array = layer_union(component, LAYERS['array'])
        o = dict(jj_main_vs_shadow_um2=round(main_jj.intersection(shadow).area, 6),
                 jj_main_vs_array_main_um2=round(main_jj.intersection(main_array).area, 6),
                 array_main_vs_jj_shadow_um2=round(main_array.intersection(shadow).area, 6))
        assert o['jj_main_vs_shadow_um2'] == 0 and o['jj_main_vs_array_main_um2'] == 0, (component, o)
        # NCKU's inner return shadow overlaps the return link by 10 nm (0.059 um^2); allow only that.
        assert o['array_main_vs_jj_shadow_um2'] < 0.06, (component, o)
        overlaps[component] = o
    return dict(contacts=contacts, main_shadow_overlaps=overlaps,
                noncontact_min_clearance_um={k: round(v, 4) for k, v in sorted(clearance.items())})


def plot(shapes, base, rows, title, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = {LAYERS['jj']: '#d62728', LAYERS['jj_shadow']: '#9467bd',
              LAYERS['array']: '#1f77b4', LAYERS['array_shadow']: '#17becf'}
    sx, sy = SITES['TC']
    windows = [('TC SQUID (JJ1 | JJ2)', (sx - 21, sy - 10, sx + 23, sy + 10)),
               ('TC JJ1 core', (sx - 16.2, sy - 1.8, sx - 10.6, sy + 2.2))]
    for name in ('FX1', 'FX2'):
        x, y, side = SITES[name]
        windows.append((f'{name} ({side})', (x - 58, y - 23, x + 12, y + 23) if side == 'left'
                        else (x - 12, y - 23, x + 58, y + 23)))
        jx = x + 4.5 if side == 'left' else x - 6.5
        jy = y + 2.8 if side == 'left' else y - 2.9
        windows.append((f'{name} small JJ', (jx - 1.8, jy - 1.6, jx + 1.8, jy + 1.4)))
    fig, axes = plt.subplots(3, 2, figsize=(13, 15))
    for ax, (name, w) in zip(axes.flat, windows):
        W = box(*w)
        metal = local_metal(base, W)
        for g in getattr(metal, 'geoms', [metal]):
            if g.is_empty:
                continue
            ax.fill(*g.exterior.xy, color='0.75', lw=0)
            for h in g.interiors:
                ax.fill(*h.xy, color='white', lw=0)
        for layer in (LAYERS['array_shadow'], LAYERS['array'], LAYERS['jj_shadow'], LAYERS['jj']):
            for _, _, l, geom in shapes:
                if l == layer and geom.intersects(W):
                    ax.fill(*geom.exterior.xy, color=colors[l], alpha=.35 if l in (12, 14) else .7,
                            lw=.5, ec=colors[l])
        ax.set_xlim(w[0], w[2]); ax.set_ylim(w[1], w[3]); ax.set_aspect('equal')
        ax.set_title(name, fontsize=10); ax.grid(alpha=.25)
        ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)')
    labels = ['base metal 1/0', 'small-JJ main 11/0', 'shadow 12/0', 'array main 13/0', 'array shadow 14/0']
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, alpha=.7) for c in
               ('0.75', colors[11], colors[12], colors[13], colors[14])]
    fig.legend(handles, labels, loc='lower center', ncol=5)
    fig.suptitle(title)
    fig.tight_layout(rect=(0, .03, 1, .97)); fig.savefig(path, dpi=120); plt.close(fig)


def main():
    argparse.ArgumentParser(description=__doc__.split('\n')[0]).parse_args()
    parent = ROOT / PARENT
    stem = f'{DESIGN.name}_{REVISION}_design3_nckujj'
    gds = DESIGN / 'gds' / f'{stem}_fab.gds'
    fig = DESIGN / 'notes' / 'fig' / f'{stem}_sites.png'
    record_path = DESIGN / 'parameters' / f'{stem}_record.json'
    for target in (gds, fig, record_path):
        if target.exists():
            raise FileExistsError(f'Refusing to overwrite {target}; choose a new REVISION')
    fig.parent.mkdir(parents=True, exist_ok=True)

    library = gdstk.read_gds(str(parent))
    assert (library.unit, library.precision) == (1e-6, 1e-9)
    parent_top, = library.top_level()
    base = parent_top.get_polygons()
    assert {(p.layer, p.datatype) for p in base} == {BASE}, 'parent must hold only 1/0'
    rows = measure(base)
    design, shapes = build(rows)
    check(shapes, base, rows)

    out = gdstk.Library(unit=1e-6, precision=1e-9)
    cell = out.new_cell(parent_top.name)
    cell.add(*base)
    for _, _, layer, geom in shapes:
        cell.add(gdstk.Polygon(list(geom.exterior.coords)[:-1], layer=layer))
    out.write_gds(str(gds), max_points=8190)

    # Read back: base metal identical, junction shapes identical, checks rerun on the file.
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
    plot(readback, new_base, rows, 'Design3 r11 fab + NCKU-style e-beam junctions (draft, um)', fig)

    sizes = {k: {kk: round(vv, 6) for kk, vv in vars(jjm.size_from_ej(e)).items()} for k, e in EJ.items()}
    record = dict(
        status='draft development export; junction drawings not fabrication, evaporation or EM validated',
        design='Design3', revision=REVISION,
        parent_fabrication_gds=PARENT, parent_sha256=sha256(parent.read_bytes()).hexdigest(),
        fabrication_gds=gds.relative_to(DESIGN).as_posix(),
        fabrication_gds_sha256=sha256(gds.read_bytes()).hexdigest(),
        figure=fig.relative_to(DESIGN).as_posix(),
        layer_map={'base_metal': '1/0 (NCKU 5/0)', 'small_jj_main': '11/0 (NCKU 1/0)',
                   'shadow': '12/0 (NCKU 2/0)', 'array_main': '13/0 (NCKU 3/0)',
                   'array_shadow': '14/0 (NCKU 4/0)'},
        junction_energy_inputs=dict(
            values_ghz=EJ,
            FX='Design3 fluxonium small-JJ EJ/h = 4 GHz (Design3 README, CZ model input)',
            TC='symmetric 11.528 / 11.528 GHz per JJ (NCKU Design3_3 W80 row; user choice 2026-11-20). '
               'Design3 CZ model branches 16.401 / 6.655 GHz are NOT used; the 9.746 GHz half-flux EJ is lost.',
            model='code/qubit/jj_resistance_model.py (NCKU W80: Delta 180 ueV, RA 625 Ohm um^2, +30/+70 nm); '
                  'estimates, not measured'),
        junction_sizes=sizes,
        measured_hooks={k: v for _, _, _, h in rows for k, v in h.items()},
        placements=[dict(name=n, component=c.__name__, options=o) for n, c, o, _ in rows],
        parameter_snapshot={name: {k: str(v) for k, v in c.options.items()}
                            for name, c in design.components.items()},
        polygons_per_layer={f'{l}/0': sum(1 for p in overlay if p.layer == l) for l in LAYERS.values()},
        checks=report,
        component_sha256={f: sha256((ROOT / 'pattern design/ebeam_junction' / f).read_bytes()).hexdigest()
                          for f in ('fluxonium_junction.py', 'transmon_junction.py')},
        model_sha256=sha256((ROOT / 'code/qubit/jj_resistance_model.py').read_bytes()).hexdigest(),
        export_script='layout/export_ncku_ebeam_junctions.py',
        export_script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        base_metal_check='exact polygon count, order and vertex equality after GDS readback',
        not_covered=['Design3 corner test junctions (no e-beam JJ added)', 'evaporation simulation',
                     'HFSS/Q3D with the drawn junctions', 'process qualification'],
        em_rerun=False)
    record_path.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(gds=record['fabrication_gds'], figure=record['figure'],
                          sha256=record['fabrication_gds_sha256'],
                          clearance=report['noncontact_min_clearance_um'],
                          overlaps=report['main_shadow_overlaps']), indent=1))


if __name__ == '__main__':
    main()
