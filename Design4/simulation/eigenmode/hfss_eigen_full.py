"""HFSS eigenmodes of the whole Design3 + Purcell filter chain (no ports, all line ends open).

The modes in 6-7 GHz are the filter mode and the three hybridised readout resonators; their
frequencies are where the narrow S21 dips of the driven sweep sit (hfss_combined.py).

Usage:  python hfss_eigen_full.py [tag]
"""
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # simulation/ (this script: simulation/<kind>/)
import sys as _sys; _sys.path[:0] = [os.path.join(HERE, _k) for _k in ('driven', 'eigenmode', 'q3d')] + [HERE]
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim  # noqa: E402
import design3_purcell as dp  # noqa: E402
from hfss_combined import COMPS  # noqa: E402

RESULTS = os.path.join(HERE, 'results')
SETUP = dict(min_freq_ghz=6.3, n_modes=6, max_delta_f=0.01, max_passes=17, min_passes=8,
             min_converged=3, pct_refinement=30, basis_order=1)


def run(tag='eig_full'):
    d = dp.build_design()
    sim = EigenmodeSim(d, 'hfss')
    for k, v in SETUP.items():
        sim.setup[k] = v
    t0 = time.time()
    sim.run(name=tag, components=COMPS,
            open_terminations=[('feed_in', 'start'), ('out_line', 'end')]
            + [(f'res_{k}', 'end') for k in dp.RESONATORS])
    conv = sim.convergence_f
    last = conv.iloc[-1]
    freqs = sorted(float(v) for v in last.values if pd.notna(v))
    print(f'=== {tag} ({(time.time() - t0) / 60:.1f} min, {len(conv)} passes) ===')
    print('CONV'); print(conv.round(4).to_string())
    print('modes (GHz):', [round(f, 4) for f in freqs])
    print('resonator design:', {k: r['f'] for k, r in dp.RESONATORS.items()})
    pd.Series(freqs, name='f_GHz').to_csv(os.path.join(RESULTS, f'{tag}_modes.csv'), index=False)
    try:
        sim.close()
    except Exception:
        pass


if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'eig_full')
