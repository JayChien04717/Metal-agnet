"""Stage 1 layout check of FluxoniumDelft (pattern design/fluxonium/fluxonium_delft.py).

Builds one FluxoniumDelft at the Design2_FTF_Andersen Q1 position, runs the component checks
(4 um base-metal width, pins outward, junction ends on metal), compares the junction-free base
metal with Q1 of the lab mask (Design2 r11 reference GDS; the lab Manhattan stubs and the
coupler crescent are expected differences), exports a clean simulation GDS (junctions as layer
50/0 annotations) and renders the figures in the plot_design_figs.py style.
No fab GDS, e-beam layers or EM runs.

Usage (repository root, qmetal):
    python -X utf8 -B code/qubit/preview_fluxonium_delft.py --out <new folder>
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'code' / 'general'))
import project_paths  # noqa: E402
project_paths.configure(__file__)

import gdstk  # noqa: E402
import numpy as np  # noqa: E402
from qiskit_metal import designs  # noqa: E402
from shapely.affinity import translate  # noqa: E402
from shapely.geometry import Point, Polygon, box  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

from fluxonium_delft import FluxoniumDelft  # noqa: E402
from metal_width import check_component  # noqa: E402
from gds_tools import export_gds  # noqa: E402
import plot_design_figs as figs  # noqa: E402

LAB_GDS = ROOT / 'Design2_FTF_Andersen' / 'gds' / 'Design2_2026-11-18_wafer_r11_reference.gds'
Q1 = (-485.0, -1000.0)          # um, Design2 Q1 centre in the lab mask
# Compare window (local um): the qubit, its three lines; leaves out the chip-level parts
# (feedline / TC route beyond the pins) and the TC coupler crescent next to pad_bot.
COMPARE = box(-490, -440, 440, 500).difference(box(50, -275, 450, -95))
EXPECTED = unary_union([box(-23, 64, -7.5, 68), box(-23, -68, -7.5, -64)])   # lab Manhattan arm stubs


def build(**options):
    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.size.size_x = '4mm'
    design.chips.main.size.size_y = '4mm'
    FluxoniumDelft(design, 'FX', options=dict(pos_x=f'{Q1[0]}um', pos_y=f'{Q1[1]}um', **options))
    return design


def check(design, name='FX'):
    comp = design.components[name]
    out = []
    centre = np.array([design.parse_value(comp.options.pos_x), design.parse_value(comp.options.pos_y)])
    for pn, pin in comp.pins.items():
        outward = np.dot(pin.normal, pin.middle - centre) > 0
        out.append(f'pin {pn}: middle {np.round(pin.middle * 1e3, 3)} um, normal {np.round(pin.normal, 3)}, '
                   f'width {pin.width * 1e3:.1f} um, outward {outward}')
        assert outward, pn
    poly = comp.qgeometry_table('poly')
    metal = unary_union(list(poly[~poly['subtract']].geometry)).buffer(design.parse_value('1nm'))
    for _, row in comp.qgeometry_table('junction').iterrows():
        ends = [Point(row.geometry.coords[0]), Point(row.geometry.coords[-1])]
        ok = all(metal.intersects(e) for e in ends)
        out.append(f'junction {row["name"]}: both ends on metal {ok}')
        assert ok, row['name']
    return out


def export(design, path):
    """Clean sim GDS: metal on 1, junction lines as 50/0 rectangles (um)."""
    from shapely.affinity import scale
    jt = design.qgeometry.tables['junction']
    k = 1 / design.parse_value('1um')
    marks = [scale(r.geometry, k, k, origin=(0, 0)).buffer(r.width * k / 2, cap_style=2)
             for _, r in jt.iterrows()]
    export_gds(design, str(path))
    lib = gdstk.read_gds(str(path))
    top, = lib.top_level()
    for m in marks:
        top.add(gdstk.Polygon(list(m.exterior.coords)[:-1], layer=50, datatype=0))
    lib.write_gds(str(path))


def gds_metal(path, layer=1):
    lib = gdstk.read_gds(str(path))
    top = max(lib.top_level(), key=lambda c: len(c.polygons))
    polys = [Polygon(p.points).buffer(0) for p in top.flatten().polygons if p.layer == layer]
    return translate(unary_union(polys), -Q1[0], -Q1[1])


def compare(component_gds):
    new = gds_metal(component_gds).intersection(COMPARE)
    lab = gds_metal(LAB_GDS).intersection(COMPARE)
    only_lab, only_new = lab.difference(new), new.difference(lab)
    unexpected_lab = only_lab.difference(EXPECTED.buffer(0.5))
    lines = [f'compare with lab Q1 ({LAB_GDS.name}), window {COMPARE.bounds} um minus coupler crescent:',
             f'  metal only in lab: {only_lab.area:.2f} um^2 '
             f'(Manhattan stubs {only_lab.intersection(EXPECTED.buffer(0.5)).area:.2f}, other {unexpected_lab.area:.2f})',
             f'  metal only in FluxoniumDelft: {only_new.area:.2f} um^2']
    return lines, only_lab, only_new


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)

    base = build(jj_metal='False')
    check_component(base.components['FX'])              # 4 um rule: base metal only
    lines = check(build())
    lines.append('check_component (4 um base-metal width): passed')
    mirrored = build(mirror='True', orientation='90')
    lines += [f'mirror + 90 deg: {s}' for s in check(mirrored)]

    mask_gds = out / 'fluxonium_delft_mask_check.gds'    # junction-free, for the lab comparison
    export_gds(base, str(mask_gds))
    cmp_lines, only_lab, only_new = compare(mask_gds)
    lines += cmp_lines

    gds = out / 'fluxonium_delft_sim.gds'
    export(build(), gds)
    metal, junction, _ = figs.load(gds)
    rev = 'draft 2026-11-22'
    x, y = Q1
    figs.render_part(out / 'fluxonium_delft.png', f'FluxoniumDelft · {rev} · whole qubit', metal, junction,
                     figs.fit([x - 520, x + 470, y - 470, y + 540]))
    figs.render_part(out / 'junction.png', f'FluxoniumDelft · {rev} · JJ + array + flux tip', metal, junction,
                     figs.fit([x - 100, x + 40, y - 50, y + 50]))
    figs.render_part(out / 'flux_launch.png', f'FluxoniumDelft · {rev} · flux-line launch', metal, junction,
                     figs.fit([x - 510, x - 400, y - 60, y + 70]))
    figs.render_part(out / 'charge_line.png', f'FluxoniumDelft · {rev} · charge line', metal, junction,
                     figs.fit([x - 110, x + 50, y - 460, y - 250]))
    lab_metal, lab_junction, _ = figs.load(LAB_GDS)
    figs.render_part(out / 'lab_q1_reference.png', 'Design2 lab Q1 · 2026-11-18_wafer_r11 · reference',
                     lab_metal, lab_junction, figs.fit([x - 520, x + 470, y - 470, y + 540]))
    (out / 'checks.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
