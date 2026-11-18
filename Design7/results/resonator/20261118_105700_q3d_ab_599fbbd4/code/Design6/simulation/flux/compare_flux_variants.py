"""Compare the Design6 flux-line variants (results of flux_metrics.py / flux_xtalk_matrix.py / flux_square_spectro.py).

A = L (layout: sep 10 um, bar 100 um), B = L short (sep 9, 60 um), C = T (sep 9, arms 2 x 50 um),
D = U hairpin (sep 9, legs 40 um apart, return leg shorted 120 um below).
Output: results/2026-11-14/flux_variants_compare.{png,json} and a printed table.
Usage (qmetal or scqenv): python driven/compare_flux_variants.py
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import json, os
import numpy as np

SIM = _S('Design6/simulation')
R = _S('Design6/results/general/legacy/2026-11-14')
VARIANTS = [('A', 'L sep10/100 (layout)', 'flux_loc_fl10', ''), ('B', 'L short sep9/60', 'flux_loc_Lshort', 'Lshort'),
            ('C', 'T sep9, 2x50', 'flux_loc_T', 'T'), ('D', 'U hairpin, drop 120', 'flux_loc_U', 'U'),
            ('B2', 'L sep9/45, fbl_x 12 (layout 11-15)', 'flux_loc_B2', 'B2')]
NAMES = ['Q1', 'TC', 'Q2']


def load(run, sq_tag):
    d = _relocate(os.path.join(R, run))
    try:
        m = json.load(open(_relocate(os.path.join(d, 'metrics.json')))); x = json.load(open(_relocate(os.path.join(d, 'xtalk_matrix.json'))))
    except FileNotFoundError:
        return None
    sq = _relocate(os.path.join(R, 'flux_pulse', f'square_spectro_Q1{"_" + sq_tag if sq_tag else ""}.json'))
    if not os.path.exists(sq) and sq_tag:
        sq = _relocate(os.path.join(R, run, f'square_spectro_Q1_{sq_tag}.json'))
    s = json.load(open(sq)) if os.path.exists(sq) else None
    return dict(m=m, x=x, s=s)


def main():
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    data = {k: (lab, load(run, t)) for k, lab, run, t in VARIANTS}
    data = {k: v for k, v in data.items() if v[1] is not None}
    rows = []
    for k, (lab, d) in data.items():
        m, x = d['m'], d['x']
        C = np.array(x['C_normalised_pct'])
        off = [abs(C[i, j]) for i in range(3) for j in range(3) if i != j]
        rows.append(dict(variant=k, label=lab, M_pH=m['M_self_pH'], I_phi0_mA=m['I_phi0_mA'],
                         ringup_max_err_30ns_pct={n: m['step_first_30ns'][n]['max_err_edge_to_30ns_pct'] for n in NAMES},
                         xtalk_pct=C.tolist(), xtalk_max_pct=max(off), xtalk_rms_pct=float(np.sqrt(np.mean(np.square(off)))),
                         diag_phi0_per_V=np.diag(np.array(x['C_chip'])).tolist()))
    json.dump(rows, open(_S('Design6/results/flux/2026-11-14/flux_variants_compare.json'), 'w'), indent=1)
    print(f"{'var':4s}{'M Q1/TC/Q2 (pH)':24s}{'I_phi0 (mA)':20s}{'xtalk max/rms (%)':20s}{'ring-up 30ns max err (%) Q1/TC/Q2'}")
    for r in rows:
        print(f"{r['variant']:4s}" + '/'.join(f"{r['M_pH'][n]:.2f}" for n in NAMES).ljust(24)
              + '/'.join(f"{r['I_phi0_mA'][n]:.2f}" for n in NAMES).ljust(20)
              + f"{r['xtalk_max_pct']:.1f} / {r['xtalk_rms_pct']:.1f}".ljust(20)
              + '/'.join(f"{r['ringup_max_err_30ns_pct'][n]:.3f}" for n in NAMES))
    fig, ax = plt.subplots(1, len(rows), figsize=(4.6 * len(rows), 4.4))
    ax = np.atleast_1d(ax)
    lim = max(r['xtalk_max_pct'] for r in rows)
    for a, r in zip(ax, rows):
        C = np.array(r['xtalk_pct']); Cd = C.copy(); np.fill_diagonal(Cd, np.nan)
        im = a.imshow(Cd, cmap='RdBu_r', vmin=-lim, vmax=lim)
        for i in range(3):
            for j in range(3):
                a.text(j, i, '100' if i == j else f'{C[i, j]:.1f}', ha='center', va='center', fontsize=11)
        a.set_xticks(range(3)); a.set_xticklabels([f'line {n}' for n in NAMES], fontsize=8)
        a.set_yticks(range(3)); a.set_yticklabels([f'loop {n}' for n in NAMES], fontsize=8)
        a.set_title(f"{r['variant']}: {r['label']}\nM " + '/'.join(f"{r['M_pH'][n]:.2f}" for n in NAMES) + ' pH', fontsize=9)
    fig.colorbar(im, ax=ax.tolist(), fraction=0.02, label='C_ij / C_jj (%)')
    fig.savefig(_S('Design6/results/flux/2026-11-14/flux_variants_compare.png'), dpi=120, bbox_inches='tight'); plt.close(fig)
    print('saved', _S('Design6/results/flux/2026-11-14/flux_variants_compare.png'))


if __name__ == '__main__':
    main()
