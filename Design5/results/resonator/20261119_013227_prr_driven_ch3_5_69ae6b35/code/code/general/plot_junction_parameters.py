"""Annotated parameter maps for FluxoniumJunction and TransmonJunction.

Highlights where junction_gap, junction_width and finger_width act (default
options). Run with the repository qmetal Python; writes the two PNG files
next to the component sources unless --out-dir is given.
"""
import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from qiskit_metal import designs

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'pattern design' / g) for g in ('', 'ebeam_junction')]
from fluxonium_junction import FluxoniumJunction
from transmon_junction import TransmonJunction

UM = 0.001
GAP = '#f5a300'


def geometry(cls, **options):
    d = designs.DesignPlanar()
    c = cls(d, 'c', options=options)
    rows = d.qgeometry.tables['poly']
    geom = {(int(r.layer), r['name']): r.geometry for _, r in rows.iterrows()}
    return c, {k: [v/UM for v in g.bounds] for k, g in geom.items()}, geom


def strip(a, b):
    """Rectangle filling the clearance between two axis-aligned boxes (um)."""
    ax0, ay0, ax1, ay1 = a
    bx0, by0, bx1, by1 = b
    if ax1 <= bx0 or bx1 <= ax0:
        x0, x1 = (ax1, bx0) if ax1 <= bx0 else (bx1, ax0)
        y0, y1 = max(ay0, by0), min(ay1, by1)
    else:
        y0, y1 = (ay1, by0) if ay1 <= by0 else (by1, ay0)
        x0, x1 = max(ax0, bx0), min(ax1, bx1)
    return x0, y0, x1, y1


def draw_layers(ax, geom):
    for (layer, _), g in geom.items():
        x, y = g.exterior.xy
        col = {1: 'crimson', 2: 'purple'}.get(layer, 'seagreen')
        ax.fill([v/UM for v in x], [v/UM for v in y], alpha=0.28, color=col, ec=col, lw=0.8)


def mark(ax, rect, color, label=None, pad=0.004):
    x0, y0, x1, y1 = rect
    ax.add_patch(Rectangle((x0-pad, y0-pad), x1-x0+2*pad, y1-y0+2*pad, fill=True,
                           fc=color, ec=color, alpha=0.75, lw=1.5, zorder=5))
    if label:
        ax.annotate(label, ((x0+x1)/2, (y0+y1)/2), fontsize=8, color=color,
                    xytext=(8, 8), textcoords='offset points', zorder=6, weight='bold')


def dim(ax, p0, p1, text, color, offset):
    ax.annotate('', p0, p1, arrowprops=dict(arrowstyle='<->', color=color, lw=1.6), zorder=7)
    ax.annotate(text, ((p0[0]+p1[0])/2, (p0[1]+p1[1])/2), xytext=offset,
                textcoords='offset points', color=color, fontsize=10, weight='bold', zorder=7,
                bbox=dict(fc='white', ec=color, lw=0.8, alpha=0.9))


def legend(ax):
    handles = [Rectangle((0, 0), 1, 1, fc=c, alpha=a, ec=c) for c, a in
               (('crimson', .3), ('purple', .3), (GAP, .75))]
    ax.legend(handles, ['layer (red)', 'layer_recess (purple)', 'junction_gap'], loc='upper left',
              bbox_to_anchor=(1.01, 1), fontsize=8)


def fluxonium(path):
    c, b, geom = geometry(FluxoniumJunction)
    fig, ax = plt.subplots(figsize=(11, 6.5))
    draw_layers(ax, geom)
    for metal, recess in ((0, 0), (1, 0), (1, 1), (1, 9), (2, 11), (2, 9), (2, 10), (2, 3), (2, 4)):
        mark(ax, strip(b[1, f'reference_{metal}'], b[2, f'reference_{recess}']), GAP)
    x0, y0, x1, y1 = b[1, 'reference_2']
    dim(ax, (4.6, y0), (4.6, y1), 'junction_width = 0.25 um\n(about centre)', 'crimson', (-75, 45))
    x0, y0, x1, y1 = b[1, 'reference_1']
    dim(ax, (x0, 2.25), (x1, 2.25), 'finger_width = 0.1 um\n(grows to the left)', 'darkred', (-60, -75))
    ax.text(2.75, 3.65, 'junction_gap = 0.075 um at all 9 orange strips\n'
            '(pad landing, finger sides and tip, bridge, both strips at the contact)',
            color='#b06f00', fontsize=9)
    ax.set_xlim(2.6, 6.6); ax.set_ylim(1.6, 4.0); ax.set_aspect('equal')
    ax.set_xlabel('x (um, local)'); ax.set_ylabel('y (um, local)')
    ax.set_title('FluxoniumJunction: main junction parameters (defaults)')
    legend(ax)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def transmon(path):
    c, b, geom = geometry(TransmonJunction)
    fig, ax = plt.subplots(figsize=(11, 6.5))
    draw_layers(ax, geom)
    for side in ('L',):
        for metal, recess in (('pad', 'landing'), ('finger', 'landing'), ('finger', 'below_bar'),
                              ('finger', 'bridge'), ('bar', 'below_bar'), ('bar', 'bridge'),
                              ('finger', 'tip'), ('bar', 'above_bar')):
            mark(ax, strip(b[1, f'{metal}_{side}'], b[2, f'{recess}_{side}']), GAP)
    x0, y0, x1, y1 = b[1, 'bar_L']
    dim(ax, (-3.0, y0), (-3.0, y1), 'junction_width = 0.25 um\n(bottom fixed on lead,\ngrows upward)', 'crimson', (15, 25))
    x0, y0, x1, y1 = b[1, 'finger_L']
    dim(ax, (x0, -0.3), (x1, -0.3), 'finger_width = 0.1 um\n(grows outward, -x)', 'darkred', (-60, -75))
    ax.text(-5.15, 1.75, 'junction_gap = 0.075 um at all 8 orange strips\n'
            '(left junction shown; the right one is its mirror about x = 0)',
            color='#b06f00', fontsize=9)
    ax.set_xlim(-5.2, -2.0); ax.set_ylim(-1.1, 2.1); ax.set_aspect('equal')
    ax.set_xlabel('x (um, local)'); ax.set_ylabel('y (um, local)')
    ax.set_title('TransmonJunction (squid=True): junction parameters (defaults)')
    legend(ax)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path)
    args = parser.parse_args()
    out = {'fluxonium': ROOT / 'pattern design/ebeam_junction/fluxonium_junction_parameters.png',
           'transmon': ROOT / 'pattern design/ebeam_junction/transmon_junction_parameters.png'}
    if args.out_dir:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        out = {k: args.out_dir / v.name for k, v in out.items()}
    fluxonium(out['fluxonium'])
    transmon(out['transmon'])
    for v in out.values():
        print(v)


if __name__ == '__main__':
    main()
