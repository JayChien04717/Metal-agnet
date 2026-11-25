"""Annotated parameter maps for FluxoniumJunctionDolan / TransmonJunctionDolan.

Draws the default geometry and labels the option
names. Run with the qmetal Python; writes two PNG files next to the component
sources unless --out-dir is given:

    python -I -X utf8 code/general/plot_ebeam_junction_parameters.py
"""
import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from qiskit_metal import designs

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'pattern design' / 'ebeam_junction')]
from fluxonium_junction_dolan import FluxoniumJunctionDolan  # noqa: E402
from transmon_junction_dolan import TransmonJunctionDolan  # noqa: E402

COLORS = {1: ('#d62728', 'small-JJ main dose (1/0)'), 2: ('#9467bd', 'shadow (2/0)'),
          3: ('#1f77b4', 'array main dose (3/0)'), 4: ('#17becf', 'array shadow (4/0)')}


def shapes(cls, **options):
    d = designs.DesignPlanar()
    cls(d, 'c', options=options)
    rows = d.qgeometry.tables['poly']
    return [(int(r.layer), [(x * 1000, y * 1000) for x, y in r.geometry.exterior.coords])
            for _, r in rows.iterrows()]


def draw(ax, geom, xlim, ylim, title):
    for layer in (4, 3, 2, 1):
        for l, pts in geom:
            if l == layer:
                face, _ = COLORS[l]
                ax.add_patch(MplPolygon(pts, closed=True, fc=face,
                                        alpha=0.35 if l in (2, 4) else 0.75, ec=face, lw=0.6))
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('x (um)')
    ax.set_ylabel('y (um)')
    ax.grid(alpha=0.2)


def dim(ax, a, b, text, offset=(0, 0), color='k', fs=7):
    ax.annotate('', xy=b, xytext=a, arrowprops=dict(arrowstyle='<->', color=color, lw=0.8))
    ax.text((a[0] + b[0]) / 2 + offset[0], (a[1] + b[1]) / 2 + offset[1], text,
            fontsize=fs, color=color, ha='center', va='center',
            bbox=dict(fc='white', ec='none', alpha=0.8, pad=0.5))


def legend(fig):
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, alpha=0.6) for c, _ in COLORS.values()]
    fig.legend(handles, [t for _, t in COLORS.values()], loc='lower center', ncol=4, fontsize=8)


def fluxonium(out):
    left, right = shapes(FluxoniumJunctionDolan), shapes(FluxoniumJunctionDolan, side='right')
    fig, axs = plt.subplots(2, 2, figsize=(13, 10), gridspec_kw=dict(height_ratios=[1, 1.15]))
    ax = axs[0, 0]
    draw(ax, left, (-60, 13), (-22, 22), "side='left' overview")
    dim(ax, (-0.501, 10.5), (1.499, 10.5), 'array_contact_lead_width\n(funnel mouth, 2)', (0, 2.5))
    dim(ax, (-2, 15), (2, 15), 'array_contact_body_width (4)', (-9, 0))
    dim(ax, (-55.108, -7), (-53.108, -7), 'array_return_link_width (2)', (0, -1.6))
    dim(ax, (-20, 6.2), (-18.799, 6.2), 'pitch', (0, 1.1))
    ax.text(-30, 0, 'num_array = 43 per arm\narray_element_length 0.998, array_bridge_gap 0.202\n'
            'array_row_width 2', fontsize=7, ha='center')
    dim(ax, (6.728, -5), (8.728, -5), 'jj_return_lead_width (2)', (2, -1.2))
    for y, t in ((6, 'array_contact_taper_start'), (8, 'array_contact_taper_end'),
                 (20, 'array_contact_end')):
        ax.axhline(y, xmin=0.84, xmax=0.9, color='gray', lw=0.6)
        ax.text(5.5, y, f'y = {y}: {t}', fontsize=6.5, va='center')
    ax = axs[0, 1]
    draw(ax, right, (-13, 60), (-22, 22), "side='right' : rotated 180 deg, JJ re-routed")
    ax.text(-6, -9, 'finger stays LEFT of\nthe cross electrode', fontsize=7, ha='center')
    ax = axs[1, 0]
    draw(ax, left, (3.4, 7.2), (1.6, 3.95), "side='left': small JJ")
    dim(ax, (4.402, 2.2), (4.482, 2.2), 'jj_finger_width\n(80 nm)', (-0.35, 0))
    dim(ax, (4.482, 2.65), (4.73, 2.65), 'jj_bridge_gap\n(0.248)', (0.05, 0.38), fs=6.5)
    dim(ax, (6.0, 2.604), (6.0, 2.696), 'jj_cross_electrode_width (92 nm)', (0, -0.28))
    dim(ax, (5.6, 2.529), (5.6, 2.604), 'shadow_clearance', (0.42, -0.1), fs=6)
    ax.text(5.7, 3.35, 'fixed: finger right edge, electrode top edge', fontsize=6.5, ha='center')
    ax = axs[1, 1]
    draw(ax, right, (-9.4, -3.8), (-4.2, -1.6), "side='right': small JJ")
    dim(ax, (-6.71, -2.4), (-6.63, -2.4), 'jj_finger_width', (-0.55, 0))
    dim(ax, (-6.63, -2.75), (-6.382, -2.75), 'jj_bridge_gap', (0.1, 0.25), fs=6.5)
    dim(ax, (-5.2, -3.092), (-5.2, -3.0), 'jj_cross_electrode_width', (0, -0.3))
    legend(fig)
    fig.suptitle('FluxoniumJunctionDolan (default options)', fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(out, dpi=170)
    plt.close(fig)


def transmon(out):
    pair = shapes(TransmonJunctionDolan, jj2_cross_electrode_width='157nm', jj1_cross_electrode_width='384nm')
    one = shapes(TransmonJunctionDolan)
    fig, axs = plt.subplots(1, 2, figsize=(13, 6.2), gridspec_kw=dict(width_ratios=[1.25, 1]))
    ax = axs[0]
    draw(ax, pair, (-19, 23), (-9, 9.5), 'SQUID (squid=True), JJ2 = JJ1 translated by squid_spacing\n'
         '(example: Design6_TXT Q1, jj1 / jj2 = 384 / 157 nm)')
    dim(ax, (-15, 8.4), (15, 8.4), 'squid_spacing (30)', (0, 0.6))
    dim(ax, (-17, 6.6), (-13, 6.6), 'upper_contact_body_width\n(4; 7.6 rounded bay)', (0, -1.6), fs=6.5)
    dim(ax, (-13, -8.2), (-9, -8.2), 'lower_contact_body_width (4)', (0, -0.55), fs=6.5)
    dim(ax, (-15, -4.5), (-11, -4.5), 'lower_contact_\noffset (4)', (0, 0.9), fs=6.5)
    dim(ax, (-12, 0.9), (-10, 0.9), 'lower_contact_lead_width (2)', (4.6, 0), fs=6.5)
    ax.text(-12.4, 2.6, 'JJ1', ha='left', fontsize=10)
    ax.text(17.6, 2.6, 'JJ2', ha='left', fontsize=10)
    ax = axs[1]
    draw(ax, one, (-15.8, -11.2), (-1.5, 1.8), 'one junction (JJ1, default 271 nm)')
    dim(ax, (-14.742, 1.0), (-14.662, 1.0), 'jj_finger_width (80 nm)', (-0.9, 0))
    dim(ax, (-14.662, -0.15), (-14.413, -0.15), 'jj_bridge_gap (0.249)', (0.25, 0.16), fs=6.5)
    dim(ax, (-12.6, -0.575), (-12.6, -0.304), 'jj1_cross_electrode_width', (0, 0.45))
    dim(ax, (-13.9, -0.65), (-13.9, -0.575), 'shadow_clearance', (0.7, -0.2), fs=6)
    ax.text(-13.3, 1.55, 'upper contact body: y = 1.5 .. 7.5', fontsize=6.5, ha='center')
    legend(fig)
    fig.suptitle('TransmonJunctionDolan (default options)', fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(out, dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--out-dir', type=Path, default=ROOT / 'pattern design' / 'ebeam_junction')
    args = parser.parse_args()
    fluxonium(args.out_dir / 'fluxonium_junction_dolan_parameters.png')
    transmon(args.out_dir / 'transmon_junction_dolan_parameters.png')
    print('wrote', args.out_dir)


if __name__ == '__main__':
    main()
