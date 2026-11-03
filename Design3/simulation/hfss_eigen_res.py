"""HFSS eigenmode frequency of each readout resonator.

For every resonator one HFSS eigenmode design renders: the resonator route, its
coupled-line tee (the feedline section is open at both ends, so its own modes are far
above 10 GHz), and its qubit (simulation version: lumped inductors, so the readout pad
loads the resonator as in the real device).

The qubits sit far below the resonators (fluxonium: lumped 200 nH, ~3 GHz;
coupler: Lj = 50 nH per junction, ~2.3 GHz). The resonator mode is therefore only
dispersively pulled, by about g^2 / Delta ~ +7-10 MHz.

Usage:  python hfss_eigen_res.py [FX1 TC FX2]
"""
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim  # noqa: E402
from ftf_layout import build_design, RESONATORS  # noqa: E402
from q3d_ftf import RESULTS  # noqa: E402

# 17 passes: the 15-pass FX1 run had settled to ~1 MHz per pass; 25 passes pushed
# hf3d to >14 GB of RAM (30 % mesh refinement per pass)
SETUP = dict(min_freq_ghz=4.0, n_modes=2, max_delta_f=0.05, max_passes=17, min_passes=3,
             min_converged=2, pct_refinement=30, basis_order=1)
TC_LJ = '50nH'


def run(names):
    rows = []
    for name in names:
        d = build_design(for_process=True, tmon_options=dict(Lj=TC_LJ))
        sim = EigenmodeSim(d, 'hfss')
        for k, v in SETUP.items():
            sim.setup[k] = v
        comps = [f'res_{name}', f'tee_{name}', name]
        opens = [(f'tee_{name}', 'prime_start'), (f'tee_{name}', 'prime_end')]
        t0 = time.time()
        sim.run(components=comps, open_terminations=opens)
        conv = sim.convergence_f
        conv.to_csv(os.path.join(RESULTS, f'eigen_{name}_convergence.csv'))
        last = conv.iloc[-1]
        freqs = sorted(float(v) for v in last.values if pd.notna(v))
        target = RESONATORS[name]['f']
        f_res = min(freqs, key=lambda f: abs(f - target))
        rows.append(dict(resonator=name, f_target_GHz=target, f_hfss_GHz=f_res,
                         error_MHz=(f_res - target) * 1e3, all_modes_GHz=str([round(f, 4) for f in freqs]),
                         passes=len(conv), minutes=round((time.time() - t0) / 60, 1),
                         length_mm=round(RESONATORS[name]['length'], 3)))
        print(f'{name}: modes {freqs} GHz -> resonator {f_res:.4f} GHz (target {target}) '
              f'after {len(conv)} passes, {(time.time() - t0) / 60:.1f} min')
        try:
            sim.close()
        except Exception:
            pass
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, 'eigen_resonators.csv'), index=False)
    print('\n', df.to_string(index=False))
    return df


if __name__ == '__main__':
    run(sys.argv[1:] or list(RESONATORS))
