"""Loaded eigenmode of ONE tee+resonator on a plain line with 50 Ohm ports (no Purcell filter).
Gives the real f and kappa (Q) of the hanger: tells where / how wide the dip must be, or that the
resonator is not coupled (Q huge). Usage: python hfss_eigen_plain1.py [FX1|TC|FX2]
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
import design3_purcell as dp  # noqa: E402

RESULTS = os.path.join(HERE, 'results')
SETUP = dict(min_freq_ghz=6.0, n_modes=3, max_delta_f=0.01, max_passes=17, min_passes=8,
             min_converged=3, pct_refinement=30, basis_order=1)


def run(name):
    d = dp.build_design()
    sim = EigenmodeSim(d, 'hfss')
    for k, v in SETUP.items():
        sim.setup[k] = v
    t0 = time.time()
    sim.run(name=f'eigplain_{name}', components=[f'tee_{name}', f'res_{name}'],
            open_terminations=[(f'res_{name}', 'end')],
            port_list=[(f'tee_{name}', 'prime_start', 50), (f'tee_{name}', 'prime_end', 50)])
    print(f'=== {name} ({(time.time() - t0) / 60:.1f} min) ===')
    print(sim.convergence_f.round(4).to_string())
    freqs, kap = sim.renderer.pinfo.design.solutions.eigenmodes()
    print('RAW freqs', freqs, 'kappa/2pi (GHz)', kap)
    rows = [dict(f_GHz=f, kappa_MHz=(kap[i] * 1e3 if kap else float('nan')),
                 Q=(f / kap[i] if kap and kap[i] > 0 else float('inf'))) for i, f in enumerate(freqs)]
    df = pd.DataFrame(rows)
    print(df.to_string())
    df.to_csv(os.path.join(RESULTS, f'eigplain_{name}_modes.csv'), index=False)
    try:
        sim.close()
    except Exception:
        pass


if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'FX1')
