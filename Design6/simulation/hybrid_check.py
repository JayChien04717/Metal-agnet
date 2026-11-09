"""Validate the hybrid EM-circuit model: HFSS block (hfss_block6.py) + ideal lines + Q3D caps (raw, no fitted corrections)
against full-channel HFSS runs of the same geometry.

Usage: python hybrid_check.py <block tag> <run tag>[:dLr] [<run tag>[:dLr] ...]
  e.g. python hybrid_check.py blk_T5 c6_T5w c6_T5r:0.1      (c6_T5r = resonator +0.1 mm -> Lr_rest + 0.1)
Output: results/hybrid_<block>.png, dips printed (hybrid / ideal model / HFSS).
"""
import os, sys, json, copy, glob, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results'); sys.path.insert(0, HERE)
import theory_abcd as T
from scipy.signal import find_peaks


def dips(f, s21, prom=3):
    a = 20 * np.log10(np.abs(s21)); i, _ = find_peaks(-a, prominence=prom)
    return [(round(float(f[k]), 4), round(float(a[k]), 1)) for k in i]


def run_data(tag):
    d = [x for x in glob.glob(os.path.join(RES, f'*_{tag}')) if os.path.isdir(x)][-1]
    P = json.load(open(os.path.join(d, f'{tag}_theory_params.json')))
    for k in ('cin', 'cp', 'cj'): P[k] = tuple(P[k])
    S = pd.read_csv(os.path.join(d, f'{tag}_S.csv'), index_col=0).map(complex)
    return P, S.index.values.astype(float), S.S21.values


def hybrid_P(P, blk_tag, dLr=0.0, dLf=0.0):
    meta = json.load(open(os.path.join(RES, 'blocks', f'{blk_tag}_meta.json')))
    Q = copy.deepcopy(P)
    Q['block'] = T.load_block(os.path.join(RES, 'blocks', f'{blk_tag}_S3.csv'))
    Q['Lf_rest'] = meta['Lf_rest'] + dLf; Q['Lr_rest'] = meta['Lr_rest'] + dLr
    return Q


if __name__ == '__main__':
    blk = sys.argv[1]; runs = sys.argv[2:] or ['c6_T5w', 'c6_T5r:0.1']
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(len(runs), 1, figsize=(10, 3.8 * len(runs)), squeeze=False)
    for ax, r in zip(axs[:, 0], runs):
        tag, dl = (r.split(':') + ['0'])[:2]; dl = float(dl)
        P, f, s21 = run_data(tag)
        Ph = hybrid_P(P, blk, dLr=dl if tag != 'c6_T5w' else 0.0)
        # the archived params already contain the run's own L_res; the block meta was taken from the base layout,
        # so the resonator extension of the run is added to Lr_rest explicitly (dl)
        _, sh = T.s21_theory(f, Ph); _, si = T.s21_theory(f, P)
        print(f'{tag:8s} HFSS {dips(f, s21)}\n         hybrid {dips(f, sh)}\n         ideal  {dips(f, si)}')
        ax.plot(f, 20 * np.log10(np.abs(s21)), label='HFSS full channel')
        ax.plot(f, 20 * np.log10(np.abs(sh)).clip(-90), '--', label=f'hybrid (block {blk} + lines + Q3D C)')
        ax.plot(f, 20 * np.log10(np.abs(si)).clip(-90), ':', label='ideal model (Q3D C)')
        ax.set_title(f'{tag}  (resonator {dl:+.2f} mm)'); ax.set_ylabel('S21 (dB)'); ax.grid(alpha=.3); ax.legend(fontsize=8)
    axs[-1, 0].set_xlabel('f (GHz)'); fig.tight_layout()
    out = os.path.join(RES, f'hybrid_{blk}.png'); fig.savefig(out, dpi=120); print('saved', out)
