"""Plain open readout line (NO Purcell filter, NO interdigital caps): Design3's three tees on a
straight feedline with the three hanger resonators, 50 Ohm lumped ports on the two outer tee
prime pins. Discrete windows around each resonance; the hanger dips must show here.

Usage: python hfss_plain.py [f1 f2 ...]   (GHz window centres; default 6.5575 6.6917 6.7849)
"""
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from qiskit_metal.analyses.simulation import ScatteringImpedanceSim  # noqa: E402
import design3_purcell as dp  # noqa: E402

RESULTS = os.path.join(HERE, 'results')
HALF_MHZ, STEP_MHZ = 15.0, 0.3
COMPS = ([f'tee_{k}' for k in dp.RESONATORS] + [f'res_{k}' for k in dp.RESONATORS]
         + ['feed_FX1_TC', 'feed_TC_FX2'])


def run(f0):
    d = dp.build_design()
    sim = ScatteringImpedanceSim(d, 'hfss')
    for k, v in dict(freq_ghz=f0, max_delta_s=0.01, max_passes=10, min_passes=2, min_converged=1,
                     pct_refinement=30, basis_order=1).items():
        sim.setup[k] = v
    sim.setup.sweep_setup.update(name='Sweep', start_ghz=f0 - HALF_MHZ * 1e-3, stop_ghz=f0 + HALF_MHZ * 1e-3,
                                 count=int(round(2 * HALF_MHZ / STEP_MHZ)) + 1, type='Discrete',
                                 save_fields=False)
    tag = f'plain_{f0:.4f}'
    t0 = time.time()
    sim.run(name=tag, components=COMPS,
            open_terminations=[(f'res_{k}', 'end') for k in dp.RESONATORS],
            port_list=[('tee_FX1', 'prime_start', 50), ('tee_FX2', 'prime_end', 50)])
    _, _, S = sim.renderer.get_params(['S11', 'S21'])
    S.index.name = 'f_GHz'
    S.to_csv(os.path.join(RESULTS, f'{tag}_S.csv'))
    f = np.asarray(S.index, float)
    s21 = 20 * np.log10(np.abs(S['S21'].values))
    i = int(np.argmin(s21))
    print(f'=== {tag} ({(time.time() - t0) / 60:.1f} min) ===')
    print(f'WINDOW {f0:.4f}: min S21 {s21[i]:.2f} dB at {f[i]:.4f} GHz, window max {s21.max():.2f} dB, '
          f'depth {s21.max() - s21[i]:.2f} dB')
    try:
        sim.close()
    except Exception:
        pass


if __name__ == '__main__':
    for a in (sys.argv[1:] or ['6.5575', '6.6917', '6.7849']):
        run(float(a))
