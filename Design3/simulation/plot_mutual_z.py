"""Impedance response of the flux line -> fluxonium loop extraction (reading guide).

Reads results/runs/<tag>_Z.csv written by hfss_mutual_fx.py:
    broadband (default mutual_fx_broadband, 0.01-6 GHz) and
    low-f (default mutual_fx_p15, 10-100 MHz, converged mesh)
and writes ../fig/mutual/mutual_fx_impedance.png.

Usage:  python plot_mutual_z.py [broadband_tag] [lowf_tag]
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, 'results', 'runs')
FIG = os.path.join(os.path.dirname(HERE), 'fig', 'mutual')

BLUE, ORANGE = '#2a78d6', '#eb6834'      # HFSS data / model
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'


def load(tag):
    d = pd.read_csv(os.path.join(RUNS, f'{tag}_Z.csv'))
    # runs before the f_min port fix: loop and flux columns swapped if L_loop is pH-sized
    if d['L_loop_nH'].iloc[0] < 1:
        d['L_loop_nH'], d['L_fl_pH'] = d['L_fl_pH'] / 1e3, d['L_loop_nH'] * 1e3
    return d


def style(ax):
    ax.grid(True, color=GRID, lw=0.8)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2)


def main(bb_tag='mutual_fx_broadband', lf_tag='mutual_fx_p15'):
    bb, lf = load(bb_tag), load(lf_tag)
    f = bb['f_MHz'].to_numpy() / 1e3                     # GHz
    imz = bb['ImZ21_ohm'].to_numpy()
    M0 = lf['M_pH'].iloc[0]                              # quasi-static M, converged mesh
    # loop resonance: where Im Z22 (= w L_loop) changes sign
    fr = f[np.argmax(bb['L_loop_nH'].to_numpy() < 0)]
    w = 2 * np.pi * f * 1e9

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))

    # (a) Im Z21 vs f: straight line w*M at low f, pole at the loop resonance
    ax = axs[0]
    ax.loglog(f, np.abs(imz), color=BLUE, lw=2, label='HFSS |Im Z21|')
    ax.loglog(f, w * M0 * 1e-12, color=ORANGE, lw=2, ls='--', label=f'ωM, M = {M0:.3f} pH')
    ax.axvline(fr, color=INK2, lw=1, ls=':')
    ax.text(fr, ax.get_ylim()[0] * 3, f'  loop resonance\n  {fr:.2f} GHz', color=INK2, fontsize=9)
    ax.set_xlabel('f (GHz)', color=INK)
    ax.set_ylabel('|Im Z21| (Ω)', color=INK)
    ax.set_title('(a) Im Z21: slope = M at low f', color=INK, loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    style(ax)

    # (b) M_eff = Im Z21 / w: flat plateau = M, then the resonance
    ax = axs[1]
    meff = bb['M_pH'].to_numpy()
    ax.semilogx(f, meff, color=BLUE, lw=2, label='HFSS Im Z21 / ω')
    ax.semilogx(f[f < fr * 0.98], M0 / (1 - (f[f < fr * 0.98] / fr) ** 2), color=ORANGE, lw=2, ls='--',
                label='M / (1 − f²/f_r²), inductive only')
    ax.axhline(M0, color=INK2, lw=1, ls=':')
    ax.axvspan(f[0], 0.1, color=BLUE, alpha=0.08, lw=0)
    ax.text(0.012, M0 * 1.15, 'read M here\n(plateau)', color=INK2, fontsize=9)
    ax.set_ylim(-4 * M0, 12 * M0)
    ax.set_xlabel('f (GHz)', color=INK)
    ax.set_ylabel('Im Z21 / ω (pH)', color=INK)
    ax.set_title('(b) effective M vs frequency', color=INK, loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc='upper left', bbox_to_anchor=(0, 0.88))
    style(ax)

    # (c) the plateau, converged mesh
    ax = axs[2]
    fl = lf['f_MHz'].to_numpy()
    ax.plot(fl, lf['M_pH'], color=BLUE, lw=2, marker='o', ms=8, label='M = Im Z21 / ω')
    ax.set_ylim(M0 * 0.97, M0 * 1.04)
    ax.set_xlabel('f (MHz)', color=INK)
    ax.set_ylabel('M (pH)', color=INK)
    ax.set_title(f'(c) 10–100 MHz: M = {M0:.3f} pH, I_Φ0 = {2.0678e-15 / (M0 * 1e-12) * 1e3:.2f} mA',
                 color=INK, loc='left', fontsize=11)
    ax.text(fl[-1], lf['M_pH'].iloc[-1], f"  +{(lf['M_pH'].iloc[-1] / M0 - 1) * 100:.1f} %",
            color=INK2, fontsize=9, va='center')
    ax.text(0.02, 0.04, f"L_loop = Im Z22/ω = {lf['L_loop_nH'].iloc[0]:.1f} nH\n"
                        f"L_fl = Im Z11/ω = {lf['L_fl_pH'].iloc[0]:.1f} pH\n"
                        f"|Z21 − Z12| / |Z21| < {lf['Z21_Z12_mismatch'].max():.0e}",
            transform=ax.transAxes, color=INK2, fontsize=9)
    style(ax)

    fig.tight_layout()
    out = os.path.join(FIG, 'mutual_fx_impedance.png')
    fig.savefig(out, dpi=130, facecolor='#fcfcfb')
    print(out, 'f_r =', fr, 'GHz')


if __name__ == '__main__':
    main(*sys.argv[1:3])
