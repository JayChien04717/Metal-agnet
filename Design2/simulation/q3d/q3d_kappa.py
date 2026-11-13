"""Design2 resonator linewidth kappa from Q3D: feedline TL + the two lambda/2 resonators (lab
style: each starts at an OpenToGround 32 um left of the feedline and runs parallel to it for
`start_straight`). The qubits are not rendered; the resonators' qubit ends and the feedline
ends are open-terminated.

C_1 = -C[resonator net, feedline net], then (lab formula, FXTFX24 cell 42)

    Q_c = pi / (Z_0 * Z_r * (C_1 * omega)**2),   kappa = f / Q_c

Usage:  python q3d_kappa.py TAG ['{"Q1": {"start_straight": 0.8}, ...}']
"""
import json
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # simulation/
sys.path[:0] = [os.path.join(HERE, 'q3d'), HERE, os.path.join(os.path.dirname(HERE), 'pattern design')]
from qiskit_metal.analyses.quantization import LOManalysis  # noqa: E402
from ftf_layout import build_design, RESONATORS, CPW_Z0  # noqa: E402
from q3d_ftf import SETUP, RESULTS  # noqa: E402

# kappa ~ C_1^2: converged setup (as Design1 / 3 KAPPA_SETUP)
KAPPA_SETUP = dict(SETUP, max_passes=12, percent_error=0.2)


def kappa_mhz(C_1, f_ghz, Z0=CPW_Z0, Zr=CPW_Z0):
    omega = 2 * np.pi * f_ghz * 1e9
    Q_c = np.pi / (Z0 * Zr * (C_1 * omega) ** 2)
    return f_ghz * 1e9 / Q_c * 1e-6, Q_c


def target_C1(f_ghz, kappa_mhz_, Z0=CPW_Z0, Zr=CPW_Z0):
    Q_c = f_ghz * 1e3 / kappa_mhz_
    return np.sqrt(np.pi / (Q_c * Z0 * Zr)) / (2 * np.pi * f_ghz * 1e9)


def run(tag='kappa', resonators=None):
    d = build_design(for_process=True, resonators=resonators)
    names = list(RESONATORS)
    comps = ['TL'] + [f'otg_{q}' for q in names] + [f'res_{q}' for q in names]
    opens = [('TL', 'start'), ('TL', 'end')] + [(f'res_{q}', 'end') for q in names]
    lom = LOManalysis(d, 'q3d')
    for k, v in KAPPA_SETUP.items():
        lom.sim.setup[k] = v
    t0 = time.time()
    lom.sim.run(components=comps, open_terminations=opens)
    cmat = lom.sim.capacitance_matrix
    cmat.to_csv(os.path.join(RESULTS, f'{tag}_cmat_fF.csv'))
    print(f'=== {tag} ({(time.time() - t0) / 60:.1f} min) ===')
    print(cmat.round(3).to_string())
    fl = [n for n in cmat.index if 'TL' in n][0]
    rows = []
    for q in names:
        r = dict(RESONATORS[q], **(resonators or {}).get(q, {}))
        rn = [n for n in cmat.index if f'res_{q}' in n or f'otg_{q}' in n][0]
        C1 = -cmat.loc[rn, fl] * 1e-15
        kap, Qc = kappa_mhz(C1, r['f'])
        rows.append(dict(resonator=q, net=rn, f_GHz=r['f'], start_straight_mm=r['start_straight'],
                         C_1_fF=C1 * 1e15, C_1_target_fF=target_C1(r['f'], r['kappa']) * 1e15,
                         Q_c=Qc, kappa_MHz=kap, kappa_target_MHz=r['kappa']))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, f'{tag}.csv'), index=False)
    print('\nfeedline net:', fl)
    print(df.drop(columns=['net']).round(3).to_string(index=False))
    try:
        lom.sim.close()
    except Exception:
        pass
    return df


if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'kappa', json.loads(sys.argv[2]) if len(sys.argv) > 2 else None)
