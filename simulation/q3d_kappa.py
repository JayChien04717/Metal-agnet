"""Resonator linewidth kappa from the Q3D capacitance matrix of the full readout
circuit: feedline (with the three coupled-line tees) + the three resonators.

C_1 = -C[resonator net, feedline net] from the Maxwell matrix, then the lab formula
(FXTFX24 cell 42):

    omega = 2 pi f
    Q_c   = pi / (Z_0 * Z_r * (C_1 * omega)**2)
    kappa = f / Q_c

The qubits are not rendered; the resonators' qubit ends and the feedline ends at the
launch pads are open-terminated.

Usage:  python q3d_kappa.py [tag]
"""
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from qiskit_metal.analyses.quantization import LOManalysis  # noqa: E402
from ftf_layout import build_design, RESONATORS, CPW_Z0  # noqa: E402
from q3d_ftf import SETUP, RESULTS  # noqa: E402

# kappa ~ C_1^2 and this structure is large: 6 passes leave C_1 noisy (+-0.5 fF),
# so use a converged setup here
KAPPA_SETUP = dict(SETUP, max_passes=12, percent_error=0.2)

Z_0 = Z_r = CPW_Z0          # feedline and resonator use the same CPW


def kappa_mhz(C_1, f_ghz, Z0=Z_0, Zr=Z_r):
    omega = 2 * np.pi * f_ghz * 1e9
    Q_c = np.pi / (Z0 * Zr * (C_1 * omega) ** 2)
    return f_ghz * 1e9 / Q_c * 1e-6, Q_c


def net_of(cmat, *keys):
    for n in cmat.index:
        if any(k in n for k in keys):
            return n
    raise KeyError(f'{keys} not in {list(cmat.index)}')


def run(tag='kappa'):
    d = build_design(for_process=True)
    names = list(RESONATORS)
    feed = ['feed_in', 'feed_out'] + [n for n in d.components if n.startswith('feed_') and n not in ('feed_in', 'feed_out')]
    comps = feed + [f'tee_{k}' for k in names] + [f'res_{k}' for k in names]
    opens = [('feed_in', 'start'), ('feed_out', 'end')] + [(f'res_{k}', 'end') for k in names]
    lom = LOManalysis(d, 'q3d')
    for k, v in KAPPA_SETUP.items():
        lom.sim.setup[k] = v
    t0 = time.time()
    lom.sim.run(components=comps, open_terminations=opens)
    cmat = lom.sim.capacitance_matrix
    cmat.to_csv(os.path.join(RESULTS, f'{tag}_cmat_fF.csv'))
    print(f'=== {tag} ({(time.time() - t0) / 60:.1f} min) ===')
    print(cmat.round(3).to_string())

    fl = net_of(cmat, 'prime_cpw', 'feed')
    rows = []
    for k, r in RESONATORS.items():
        rn = net_of(cmat, f'res_{k}', f'second_cpw_tee_{k}')
        C1 = -cmat.loc[rn, fl] * 1e-15
        kap, Qc = kappa_mhz(C1, r['f'])
        rows.append(dict(resonator=k, net=rn, f_GHz=r['f'], C_1_fF=C1 * 1e15,
                         C_1_target_fF=r['Cc'], Q_c=Qc, kappa_MHz=kap, kappa_target_MHz=r['kappa'],
                         tee_coupling_length_mm=r['coupling_length']))
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
    run(sys.argv[1] if len(sys.argv) > 1 else 'kappa')
