"""Draw the Design2 or Design3_3 junction model over its standard GDS."""
import argparse
import json
from pathlib import Path
import sys

import gdstk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def get_design(name):
    folder = ROOT / name
    if name == 'Design2':
        sys.path.insert(0, str(folder / 'pattern design'))
        from layout3_full import build_linear_design
        model = build_linear_design(folder / 'gds/layout3 ftf.gds')
        return folder / 'gds/ftf_layout.gds', model, ('fx_left_top', 'transmon', 'fx_right_top')
    sys.path.insert(0, str(folder))
    from layout import build_design
    return (folder / 'gds/Design3_3_2026-11-15_nocheese.gds', build_design(for_process=True),
            ('FX1', 'TC', 'FX2'))


def design33_cz_parameters():
    folder = ROOT / 'Design3_3'
    sys.path.insert(0, str(folder / 'gate_model'))
    from validate_gate import extracted
    import verify as cz

    electrical, _ = extracted(folder / 'results/arms6p7_18p9_cmat_fF.csv')
    ec_tc = electrical['Ec_TC [GHz]']
    idle = cz.coupler_EJ(ec_tc, cz.TC_IDLE_GHZ)
    at_cz = cz.coupler_EJ(ec_tc, cz.TC_CZ_GHZ)
    return cz.FT, idle, at_cz, cz.TC_IDLE_GHZ, cz.TC_CZ_GHZ


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('design', choices=('Design2', 'Design3_3'))
    args = parser.parse_args()
    mask, model, components = get_design(args.design)
    design2 = args.design == 'Design2'
    if design2:
        estimates = json.loads(
            (ROOT / 'Design2/simulation/results/layout3_full/transition_estimates.json')
            .read_text(encoding='utf-8')
        )
    else:
        fx, idle, at_cz, f_idle, f_cz = design33_cz_parameters()
        tc_junctions = ((idle + at_cz) / 2, (idle - at_cz) / 2)
    top, = gdstk.read_gds(mask).top_level()
    metal = [p for p in top.get_polygons() if p.layer == 1]
    table = model.qgeometry.tables['junction']

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4), layout='constrained')
    for ax, component, label in zip(axes, components, ('FX1', 'TC SQUID', 'FX2')):
        rows = list(table.loc[table.component.eq(model.components[component].id)].itertuples())
        segments = [np.asarray(row.geometry.coords) * 1000 for row in rows]
        points = np.concatenate(segments)
        x, y = (points.min(axis=0) + points.max(axis=0)) / 2
        xl, yl = (x - 43, x + 43), (y - 38, y + 38)
        visible = []
        for p in metal:
            (x0, y0), (x1, y1) = p.bounding_box()
            if x1 > xl[0] and x0 < xl[1] and y1 > yl[0] and y0 < yl[1]:
                visible.append(p.points)
        ax.set_facecolor('#555f6c')
        ax.add_collection(PolyCollection(visible, facecolor='#e4e7eb', edgecolor='none'))
        ax.set(xlim=xl, ylim=yl, xlabel='x (µm)', ylabel='y (µm)', title=label, aspect='equal')
        for i, (row, line) in enumerate(zip(rows, segments)):
            is_array = row.name.startswith('array')
            color = '#4bb4dc' if is_array else '#f0802a'
            ax.plot(*line.T, lw=3, color=color, solid_capstyle='butt')
            if not is_array and not design2 and label == 'TC SQUID':
                left = i == 0
                ax.annotate(f'CZ model branch {chr(65+i)}\nEJ/h = {tc_junctions[i]:.3f} GHz',
                            xy=line.mean(axis=0), xytext=(x + (-40 if left else 40), y - 34),
                            ha='left' if left else 'right', va='bottom', fontsize=9,
                            bbox=dict(boxstyle='round', fc='white', ec=color),
                            arrowprops=dict(arrowstyle='->', color=color, lw=1.3))
        if label.startswith('FX'):
            if design2:
                ej, el = estimates['FX_EJ_GHz'], estimates['FX_EL_GHz']
                bias = estimates['FX_flux_in_Phi0']
            else:
                ej, el = fx[:2] if label == 'FX1' else fx[2:]
                bias = 0.5
            note = f'Fluxonium Hamiltonian at {bias:g} Φ₀\nEJ/h = {ej:g} GHz\nEL/h = {el:g} GHz'
        elif design2:
            note = 'SQUID junction sites'
        else:
            note = (f'CZ SQUID effective EJ/h\n{idle:.3f} GHz at zero flux · f01 {f_idle:g} GHz\n'
                    f'{at_cz:.3f} GHz at half flux · f01 {f_cz:g} GHz')
        ax.text(.04, .95, note, va='top', transform=ax.transAxes, fontsize=9,
                bbox=dict(boxstyle='round', fc='white', ec='0.6'))

    model_title = 'Junction model' if design2 else 'CZ junction parameters'
    fig.suptitle(f'{args.design} · {model_title} on standard base metal', fontsize=15)
    caption = 'Orange: junction sites   ·   Blue: array legs'
    if not design2:
        caption += '   ·   TC model branch energies are shown without a fabrication-side assignment'
    fig.supxlabel(caption, fontsize=11)
    out = mask.parents[1] / 'fig/junction_parameters.png'
    fig.savefig(out, dpi=150)
    print(out)


if __name__ == '__main__':
    main()
