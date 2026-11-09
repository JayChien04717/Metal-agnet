"""HFSS driven-modal S21 of Design3's readout line + Purcell filter (no qubits).

Renders feedline, tees, the three hanger resonators (open ends), the two interdigital caps,
the filter meander and the output line of design3_purcell.build_design(). 50 Ohm lumped
ports on the outer ends of feed_in and out_line (launch pads not rendered).
Looks for the Purcell pass-band peak and the resonator dips inside it.

Usage:  python hfss_combined.py [tag]
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
from qiskit_metal.analyses.simulation import ScatteringImpedanceSim  # noqa: E402
import design3_purcell as dp  # noqa: E402

RESULTS = os.path.join(HERE, 'results')
os.makedirs(RESULTS, exist_ok=True)
SETUP = dict(freq_ghz=6.7, max_delta_s=0.02, max_passes=8, min_passes=2, min_converged=1,
             pct_refinement=30, basis_order=1)
SWEEP = dict(name='Sweep', start_ghz=6.5, stop_ghz=6.9, count=201, type='Discrete',
             save_fields=False)
COMPS = (['feed_in', 'feed_out', 'out_line', 'cap_in', 'cap_out', 'filter_in']
         + [f'tee_{k}' for k in dp.RESONATORS] + [f'res_{k}' for k in dp.RESONATORS]
         + ['feed_FX1_TC', 'feed_TC_FX2'])


def run(tag='comb_v1'):
    d = dp.build_design()
    sim = ScatteringImpedanceSim(d, 'hfss')
    for k, v in SETUP.items():
        sim.setup[k] = v
    for k, v in SWEEP.items():
        sim.setup.sweep_setup[k] = v
    t0 = time.time()
    sim.run(name=tag, components=COMPS,
            open_terminations=[(f'res_{k}', 'end') for k in dp.RESONATORS],
            port_list=[('feed_in', 'start', 50), ('out_line', 'end', 50)])
    freqs, _, S = sim.renderer.get_params(['S11', 'S21'])
    S = S.copy()
    S.index.name = 'f_GHz'
    S.to_csv(os.path.join(RESULTS, f'{tag}_S.csv'))
    f = np.asarray(S.index, float)
    s21 = 20 * np.log10(np.abs(S['S21'].values))
    i0 = int(np.argmax(s21))
    print(f'=== {tag} ({(time.time() - t0) / 60:.1f} min) ===')
    print(f'Purcell peak: {f[i0]:.4f} GHz, S21 = {s21[i0]:.2f} dB')
    from scipy.signal import find_peaks
    # dips = peaks of -S21 on top of the smooth pass band (prominence > 1 dB)
    idx, pr = find_peaks(-s21, prominence=1.0)
    for i, p in zip(idx, pr['prominences']):
        print(f'  dip at {f[i]:.4f} GHz: S21 = {s21[i]:.1f} dB (depth {p:.1f} dB)')
    print('resonator targets:', {k: r['f'] for k, r in dp.RESONATORS.items()})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(f, s21, label='S21')
    ax.plot(f, 20 * np.log10(np.abs(S['S11'].values)), label='S11')
    for k, r in dp.RESONATORS.items():
        ax.axvline(r['f'], color='gray', ls=':')
    ax.set_xlabel('f (GHz)'); ax.set_ylabel('dB'); ax.legend(); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(RESULTS, f'{tag}_S.png'), dpi=120)
    try:
        sim.close()
    except Exception:
        pass


if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'comb_v1')
