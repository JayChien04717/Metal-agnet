# -*- coding: utf-8 -*-
r"""Compute CZ pulse crosstalk / distortion metrics and render L-vs-T SQUID layout comparison figure."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as MplPoly
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DESIGN_DIR = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / 'code' / 'general'),
    str(DESIGN_DIR / 'layout'),
    str(DESIGN_DIR / 'simulation' / 'flux'),
]
import project_paths
project_paths.configure(__file__)
import hfss_mutual_d12 as HM


def draw_design(ax, d, xlim, ylim, title):
    tbl = d.qgeometry.tables['poly']
    ax.set_facecolor('#2b3a4a')
    etch_patches, metal_patches = [], []
    for _, r in tbl.iterrows():
        g = r['geometry']
        geoms = list(getattr(g, 'geoms', [g]))
        for poly in geoms:
            if poly.is_empty:
                continue
            coords = np.asarray(poly.exterior.coords) * 1e3
            if r['subtract']:
                etch_patches.append(MplPoly(coords, closed=True))
            else:
                metal_patches.append(MplPoly(coords, closed=True))
    ax.add_collection(PatchCollection(etch_patches, facecolor='#f4f4f6', edgecolor='none', zorder=2))
    ax.add_collection(PatchCollection(metal_patches, facecolor='#d4a017', edgecolor='#6b4f00', linewidth=0.35, zorder=3))
    jtbl = d.qgeometry.tables['junction']
    for _, r in jtbl.iterrows():
        c = np.asarray(r['geometry'].coords) * 1e3
        ax.plot(c[:, 0], c[:, 1], color='crimson', lw=1.6, zorder=4)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('x (um)')
    ax.set_ylabel('y (um)')


def main():
    flux_dir = DESIGN_DIR / 'results' / 'flux'
    tags = [
        ('L_r4_p12', 'L-r4 (30x42)'),
        ('T_r4loop_p8', 'T + r4 (30x42)'),
        ('T_c36_24_p8', 'T + 24x30'),
        ('T_c36_30_p8', 'T + 30x30 (sq)'),
        ('T_c36_36_p8', 'T + 36x30'),
        ('T_c32_36_bh4_p8', 'T + 36x28 (bh4)'),
    ]
    print('=== Distortion & CZ-pulse Crosstalk Analysis ===')
    for tag, lbl in tags:
        p = sorted(flux_dir.glob(f'*hfss_mutual_{tag}_*'))[-1]
        s = json.loads((p / 'summaries' / 'summary.json').read_text(encoding='utf-8'))
        sp = s['M_spread_pct']
        lf = s['L_fl_pH']
        xt = s['crosstalk_pct']
        dphi_q1_from_tc = 0.224 * xt['Q1<-TC'] / 100.0
        dphi_tc_from_q1 = 0.116 * xt['TC<-Q1'] / 100.0
        dphi_q2_total = 0.224 * xt['Q2<-TC'] / 100.0 + 0.116 * xt['Q2<-Q1'] / 100.0
        print(
            f"{lbl:16s} | spread={sp['Q1']:.3f}/{sp['TC']:.3f}/{sp['Q2']:.3f}% "
            f"| L_fl={lf['Q1']:.1f}/{lf['TC']:.1f}/{lf['Q2']:.1f} pH "
            f"| CZ err (mPhi0): Q1<-TC={dphi_q1_from_tc*1e3:+.2f}, TC<-Q1={dphi_tc_from_q1*1e3:+.2f}, Q2_tot={dphi_q2_total*1e3:+.2f}"
        )

    out_dir = sorted(flux_dir.glob('*T_screen_agg_*'))[-1] / 'figures'
    d_L = HM.build_core_design({})
    d_T30 = HM.build_core_design(dict(
        fbl=dict(fbl_type='T', fbl_x='35um', fbl_length='120um', fbl_sep='10.5um'),
        squid=dict(jj_gap='36um', squid_spacing='30um', squid_bar_height='6um'),
    ))
    d_T24 = HM.build_core_design(dict(
        fbl=dict(fbl_type='T', fbl_x='32um', fbl_length='110um', fbl_sep='10.5um'),
        squid=dict(jj_gap='36um', squid_spacing='24um', squid_bar_height='6um'),
    ))
    d_T36 = HM.build_core_design(dict(
        fbl=dict(fbl_type='T', fbl_x='38um', fbl_length='130um', fbl_sep='10.5um'),
        squid=dict(jj_gap='36um', squid_spacing='36um', squid_bar_height='6um'),
    ))

    fig, axes = plt.subplots(2, 3, figsize=(15, 9.5))
    draw_design(axes[0, 0], d_L,   (-620, 620), (-560, 240), 'Core overview: L-r4 baseline')
    draw_design(axes[0, 1], d_T30, (-620, 620), (-560, 240), 'Core overview: T + 30x30 (sq, recommended)')
    draw_design(axes[0, 2], d_T30, (-120, 150), (-480, -250), 'TC SQUID + T-line zoom (+x arm toward Q2)')
    draw_design(axes[1, 0], d_L,   (-490, -310), (-300, -120), 'Q1 SQUID zoom: L-r4 (30x42)')
    draw_design(axes[1, 1], d_T30, (-520, -320), (-290, -120), 'Q1 SQUID zoom: T + 30x30 (jj_gap=36, sp=30)')
    draw_design(axes[1, 2], d_T24, (-510, -330), (-290, -120), 'Q1 SQUID zoom: T + 24x30 (jj_gap=36, sp=24)')
    fig.tight_layout()
    fig_path = out_dir / 'layout_compare_L_vs_T.png'
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)
    print('SAVED_FIG=' + str(fig_path))


if __name__ == '__main__':
    main()
