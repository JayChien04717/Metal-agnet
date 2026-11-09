"""Overlay HFSS S21 / S11 of runs in results/<time>_<tag>/ with their theory curves and tabulate capacitances.
python compare_runs.py [tag ...]  ->  results/compare_S.png, results/runs_summary.csv"""
import os, sys, json, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.signal import find_peaks
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results')
runs = sorted(glob.glob(os.path.join(RES, '*', '*_meta.json')))
if len(sys.argv) > 1:
    runs = [r for r in runs if os.path.basename(r)[:-10] in sys.argv[1:]]
fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True); rows = []
for mp in runs:
    d = os.path.dirname(mp); t = os.path.basename(mp)[:-10]
    m = json.load(open(mp)); db = pd.read_csv(os.path.join(d, f'{t}_dB.csv'), index_col=0)
    lab = f"{os.path.basename(d)}{' (filter only)' if m['filter_only'] else ''}"
    l, = a1.plot(db.index, db.S21_dB, label=lab, lw=1); a2.plot(db.index, db.S11_dB, lw=1, color=l.get_color())
    tp = os.path.join(d, f'{t}_theory_dB.csv')
    if os.path.exists(tp):
        th = pd.read_csv(tp, index_col=0); a1.plot(th.index, th.S21_dB.clip(lower=-90), '--', lw=.8, color=l.get_color())
    idx, pr = find_peaks(-db.S21_dB.values, prominence=0.5); c = m['caps']
    rows.append(dict(run=os.path.basename(d), filter_only=m['filter_only'], C_cin_fF=c.get('C_cin_fF'), C_cp_fF=c.get('C_cp_fF'),
                     C_cj_fF=c.get('C_cj_fF'), dips_GHz=';'.join(f'{db.index[i]:.4f}({db.S21_dB.values[i]:.0f}dB)' for i in idx),
                     L_filter=m['L_FILTER'], L_res=m['L_RES_TOTAL'], note=m['note']))
a1.set_ylabel('S21 (dB)  (solid HFSS, dashed theory)'); a2.set_ylabel('S11 (dB)'); a2.set_xlabel('f (GHz)')
for a in (a1, a2): a.grid(alpha=.3)
if runs: a1.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(RES, 'compare_S.png'), dpi=120)
df = pd.DataFrame(rows); df.to_csv(os.path.join(RES, 'runs_summary.csv'), index=False); print(df.to_string())
