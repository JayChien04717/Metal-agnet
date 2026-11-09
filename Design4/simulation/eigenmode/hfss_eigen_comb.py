"""HFSS eigenmode frequency of each hanger resonator of design3_purcell (tee + resonator,
open at the far end; the feedline section of the tee is open at both ends).

Usage:  python hfss_eigen_comb.py [FX1 TC FX2]
"""
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # simulation/ (this script: simulation/<kind>/)
import sys as _sys; _sys.path[:0] = [os.path.join(HERE, _k) for _k in ('driven', 'eigenmode', 'q3d')] + [HERE]
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim  # noqa: E402
import design3_purcell as dp  # noqa: E402

RESULTS = os.path.join(HERE, 'results')
SETUP = dict(min_freq_ghz=4.0, n_modes=2, max_delta_f=0.05, max_passes=12, min_passes=3,
             min_converged=2, pct_refinement=30, basis_order=1)


def run(names):
    rows = []
    for name in names:
        d = dp.build_design()
        sim = EigenmodeSim(d, 'hfss')
        for k, v in SETUP.items():
            sim.setup[k] = v
        t0 = time.time()
        sim.run(name=f'eig_{name}', components=[f'res_{name}', f'tee_{name}'],
                open_terminations=[(f'tee_{name}', 'prime_start'), (f'tee_{name}', 'prime_end'),
                                   (f'res_{name}', 'end')])
        conv = sim.convergence_f
        last = conv.iloc[-1]
        freqs = sorted(float(v) for v in last.values if pd.notna(v))
        target = dp.RESONATORS[name]['f']
        f_res = min(freqs, key=lambda f: abs(f - target))
        rows.append(dict(resonator=name, f_target_GHz=target, f_hfss_GHz=f_res,
                         error_MHz=(f_res - target) * 1e3, modes=str([round(f, 4) for f in freqs]),
                         passes=len(conv), minutes=round((time.time() - t0) / 60, 1)))
        print(f'{name}: modes {freqs} GHz -> {f_res:.4f} GHz (design {target}) after {len(conv)} passes')
        try:
            sim.close()
        except Exception:
            pass
    pd.DataFrame(rows).to_csv(os.path.join(RESULTS, 'eigen_comb.csv'), index=False)


if __name__ == '__main__':
    run(sys.argv[1:] or list(dp.RESONATORS))
