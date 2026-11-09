"""Q3D capacitance of one interdigital coupling capacitor of the Purcell filter.

Renders 'cap_in' (taper + comb + leads, ends open) of purcell_layout.build_design().
C_c = -C[comb_a, comb_b] from the Maxwell matrix; then the lab formula gives the
port Q_c = pi / (Z0 Zr (C_c w)^2) and kappa_F = f0 / (Q_c / 2) for two equal ports.

Usage:  python q3d_cap.py [tag] [overlap_um]
"""
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # simulation/ (this script: simulation/<kind>/)
import sys as _sys; _sys.path[:0] = [os.path.join(HERE, _k) for _k in ('driven', 'eigenmode', 'q3d')] + [HERE]
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from qiskit_metal.analyses.quantization import LOManalysis  # noqa: E402
import purcell_layout as pl  # noqa: E402

RESULTS = os.path.join(HERE, 'results')
os.makedirs(RESULTS, exist_ok=True)
SETUP = dict(max_passes=12, min_passes=2, percent_error=0.2, freq_ghz=6.7,
             solution_order='Medium', auto_increase_solution_order=False)


def run(tag='cap', overlap_um=None):
    d = pl.build_design(overlap_um=overlap_um)
    lom = LOManalysis(d, 'q3d')
    for k, v in SETUP.items():
        lom.sim.setup[k] = v
    t0 = time.time()
    lom.sim.run(name=tag, components=['cap_in'], open_terminations=[('cap_in', 'a'), ('cap_in', 'b')])
    cmat = lom.sim.capacitance_matrix
    cmat.to_csv(os.path.join(RESULTS, f'{tag}_cmat_fF.csv'))
    print(f'=== {tag} ({(time.time() - t0) / 60:.1f} min) ===')
    print(cmat.round(3).to_string())
    names = list(cmat.index)
    ground = [n for n in names if 'ground' in n.lower()]
    sig = [n for n in names if n not in ground]
    Cc = -cmat.loc[sig[0], sig[1]]                       # fF
    w = 2 * np.pi * pl.F0 * 1e9
    Qc = np.pi / (pl.CPW_Z0 ** 2 * (Cc * 1e-15 * w) ** 2)
    kappa = pl.F0 * 1e3 / (Qc / 2)
    print(f'\nCc = {Cc:.2f} fF (target {pl.CC:.1f}, Alley {pl.CAP_FF:.1f}); '
          f'Q_c = {Qc:.1f}; kappa_F = {kappa:.1f} MHz (target {pl.KAPPA_F}); '
          f'C to ground: {[round(cmat.loc[s, s] + cmat.loc[s, sig[1 - i]], 2) for i, s in enumerate(sig)]} fF')
    try:
        lom.sim.close()
    except Exception:
        pass
    return Cc, kappa


if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'cap',
        float(sys.argv[2]) if len(sys.argv) > 2 else None)
