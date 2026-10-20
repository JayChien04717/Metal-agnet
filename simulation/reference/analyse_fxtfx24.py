"""Ec / J / g of one FTF set of the FXTFX24-W7-0x-V.3.1 design (Q1 - Q3 - Q2), from the
Q3D Maxwell matrix saved in that notebook (cell 33), with the same method as our design
(note/cap_matrix.py). Writes results/reference/FXTFX24_cell33_cmat_fF.csv; if the notebook is not
found (it lives outside this repo), that saved matrix is used instead.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.constants import e, h

HERE = os.path.dirname(os.path.abspath(__file__))     # simulation/reference
SIM = os.path.dirname(HERE)                           # simulation
ROOT = os.path.dirname(SIM)
sys.path.insert(0, os.path.join(ROOT, 'note'))
import cap_matrix as cm  # noqa: E402

NB = os.path.join(ROOT, '..', 'FXTFX24-W7-0x-V.3.1-ALL_Sims.ipynb')


def load_cell33():
    nb = json.load(open(NB, encoding='utf-8'))
    txt = None
    for o in nb['cells'][33]['outputs']:
        t = o.get('data', {}).get('text/plain')
        if t:
            txt = ''.join(t)
    frames = []
    for b in [b for b in txt.split('\n\n') if b.strip()]:
        lines = [l for l in b.split('\n') if l.strip()]
        cols = lines[0].replace(chr(92), '').split()      # drop pandas line-continuation '\'
        rows = {}
        for l in lines[1:]:
            p = l.split()
            rows[p[0]] = [float(x) for x in p[1:]]
        frames.append(pd.DataFrame(rows, index=cols).T)
    C = pd.concat(frames, axis=1)
    return C.loc[C.columns, C.columns]


if __name__ == '__main__':
    saved = os.path.join(SIM, 'results', 'reference', 'FXTFX24_cell33_cmat_fF.csv')
    if os.path.exists(NB):
        C = load_cell33()
        C.to_csv(saved)
    else:                       # notebook not next to this repo: use the matrix saved earlier
        C = pd.read_csv(saved, index_col=0)
    print('nets:', list(C.columns))

    nodes = ['pad_bot_Q1', 'pad_top_Q1', 'tmon_0_Q3', 'pad_top_Q2', 'pad_bot_Q2']
    print('\n--- Ec / J, 5 nodes (readout pads, flux and charge lines grounded) ---')
    cm.analyse(C, nodes=nodes)
    print('C_sigma [fF]', {k: round(v, 2) for k, v in cm.total_capacitances(C, nodes).items()})

    print('\n--- f-r: beta with the readout pad as node; notebook resonators cavity1 6.109 GHz (Q1),'
          ' cavity2 6.475 GHz (Q2); Z0 = 45 Ohm ---')
    cm.readout_g(C, ('pad_bot_Q1', 'pad_top_Q1'), 'readout_padNline_Q1', Z0=45.0, f_r=6.109e9)
    cm.readout_g(C, ('pad_top_Q2', 'pad_bot_Q2'), 'readout_padNline_Q2', Z0=45.0, f_r=6.475e9)

    print('\n--- charge lines (f01 0.5 GHz, |n01| 0.1 assumed, as for our design) ---')
    for q, pads in (('Q1', ('pad_bot_Q1', 'pad_top_Q1')), ('Q2', ('pad_top_Q2', 'pad_bot_Q2'))):
        r = cm.charge_line(C, pads, f'charge_line_{q}', 0.5e9, 0.1)
        print(q, ', '.join(f'{k} = {v:.4g}' for k, v in r.items()))

    n7 = nodes + ['readout_padNline_Q1', 'readout_padNline_Q2']
    Cm = cm.reduced_matrix(C, n7)
    Mi = np.linalg.inv(cm.transform_matrix(7))
    Ci = np.linalg.inv(Mi.T @ Cm @ Mi)
    k = e * np.sqrt(np.pi / 2 / 50 / h) / 2 / np.pi / 1e6
    print(f'\nnotebook cell-36 g formula on this matrix: g1 = {abs(k * Ci[1, 5]):.1f} MHz, '
          f'g2 = {abs(k * Ci[4, 6]):.1f} MHz')
    print('--- 7 nodes (readout pads as floating nodes, as in notebook cell 36) ---')
    cm.analyse(C, nodes=n7)
