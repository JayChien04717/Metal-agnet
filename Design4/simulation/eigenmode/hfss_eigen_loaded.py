"""HFSS eigenmodes of the full chain with 50 Ohm lumped ports on both ends (loaded modes).
Prints f, kappa/2pi (= 2*Im) and Q of each mode -> which mode is the Purcell filter (low Q) and
which are the resonators (high Q), at the frequencies the ports actually pull them to.
Usage: python hfss_eigen_loaded.py [tag]
"""
import os, sys, time, warnings
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


def run(tag='eig_loaded'):
    d = dp.build_design()
    sim = EigenmodeSim(d, 'hfss')
    for k, v in SETUP.items():
        sim.setup[k] = v
    t0 = time.time()
    sim.run(name=tag, components=COMPS,
            open_terminations=[(f'res_{k}', 'end') for k in dp.RESONATORS],
            port_list=[('feed_in', 'start', 50), ('out_line', 'end', 50)])
    print(f'=== {tag} ({(time.time() - t0) / 60:.1f} min, {len(sim.convergence_f)} passes) ===')
    print(sim.convergence_f.round(4).to_string())
    sol = sim.renderer.pinfo.design.solutions
    freqs, kap = sol.eigenmodes()
    rows = []
    for i, f in enumerate(freqs):
        k = kap[i] if kap else float('nan')
        rows.append(dict(f_GHz=f, kappa_MHz=k * 1e3, Q=(f / k if kap and k > 0 else float('inf'))))
    df = pd.DataFrame(rows)
    print(df.to_string())
    df.to_csv(os.path.join(RESULTS, f'{tag}_modes.csv'), index=False)
    try:
        sim.close()
    except Exception:
        pass


if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'eig_loaded')
