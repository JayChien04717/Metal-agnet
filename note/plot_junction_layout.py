"""Draw the Design2 or Design3_3 junction model over its standard GDS."""
import argparse
from pathlib import Path
import sys

import gdstk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
from scipy.constants import e, h

ROOT = Path(__file__).resolve().parents[1]
PHI_0 = h / (2 * e)


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


def junction_params(inductance):
    ic = PHI_0 / (2 * np.pi * inductance)
    rn = np.pi * 180e-6 / (2 * ic)  # Ambegaokar–Baratoff with Al Δ/e = 180 µV
    return ic * PHI_0 / (2 * np.pi * h) / 1e9, ic / 1e-9, rn / 1e3


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('design', choices=('Design2', 'Design3_3'))
    args = parser.parse_args()
    mask, model, components = get_design(args.design)
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
            if is_array:
                continue
            lj = float(row.hfss_inductance.removesuffix('nH'))
            ej, ic, rn = junction_params(lj * 1e-9)
            xoff = -40 if len(rows) == 2 and i == 0 else 40
            ax.annotate(f'JJ {i+1}\nL = {lj:.2f} nH\nEJ/h = {ej:.2f} GHz\nIc = {ic:.1f} nA\nRn ≈ {rn:.1f} kΩ',
                        xy=line.mean(axis=0), xytext=(x + xoff, y - 34),
                        ha='left' if xoff < 0 else 'right', va='bottom', fontsize=8.5,
                        bbox=dict(boxstyle='round', fc='white', ec=color),
                        arrowprops=dict(arrowstyle='->', color=color, lw=1.3))
        if label.startswith('FX'):
            array_note = 'Array model  2 × 100 nH'
            if args.design == 'Design3_3':
                array_note += '\nCZ model  EJ/h = 4 GHz, EL/h = 0.6 GHz'
            ax.text(.04, .95, array_note, va='top', transform=ax.transAxes, fontsize=9,
                    bbox=dict(boxstyle='round', fc='white', ec='#4bb4dc'))
        elif args.design == 'Design3_3':
            ax.text(.04, .95, 'CZ model EJ/h (GHz)\n16.401 idle · 6.655 at half flux', va='top',
                    transform=ax.transAxes, fontsize=9, bbox=dict(boxstyle='round', fc='white', ec='0.6'))

    fig.suptitle(f'{args.design} · Junction model on standard base metal', fontsize=15)
    fig.supxlabel('Orange: lumped JJ   ·   Blue: lumped array legs   ·   HFSS lumped L values; Rn estimated for Al Δ = 180 µeV')
    out = mask.parents[1] / 'fig/junction_parameters.png'
    fig.savefig(out, dpi=150)
    print(out)


if __name__ == '__main__':
    main()
