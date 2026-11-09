"""Driven S21 with FINE adaptive mesh + INTERPOLATING sweep, kappa = 3 MHz design.

 plain    : three tees + resonators on a plain line, ports at the outer tee prime pins (no Purcell)
 combined : full Design3 + Purcell chain (hfss_combined.run)

Usage: python hfss_interp.py plain|combined
"""
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # simulation/ (this script: simulation/<kind>/)
import sys as _sys; _sys.path[:0] = [os.path.join(HERE, _k) for _k in ('driven', 'eigenmode', 'q3d')] + [HERE]
sys.path.insert(0, HERE)
import hfss_combined as hc  # noqa: E402
import design3_purcell as dp  # noqa: E402
from qiskit_metal.analyses.simulation import ScatteringImpedanceSim  # noqa: E402

# fine mesh: tight delta S, many passes, 30 % refinement per pass, 2nd order basis
FINE = dict(freq_ghz=6.7, max_delta_s=0.003, max_passes=15, min_passes=6, min_converged=2,
            pct_refinement=30, basis_order=2)
SWEEP = dict(name='Sweep', start_ghz=6.4, stop_ghz=7.0, count=1201, type='Interpolating', save_fields=False)


def plain(tag='interp_plain'):
    d = dp.build_design()
    sim = ScatteringImpedanceSim(d, 'hfss')
    for k, v in FINE.items():
        sim.setup[k] = v
    sim.setup.sweep_setup.update(SWEEP)
    comps = ([f'tee_{k}' for k in dp.RESONATORS] + [f'res_{k}' for k in dp.RESONATORS]
             + ['feed_FX1_TC', 'feed_TC_FX2'])
    t0 = time.time()
    sim.run(name=tag, components=comps, open_terminations=[(f'res_{k}', 'end') for k in dp.RESONATORS],
            port_list=[('tee_FX1', 'prime_start', 50), ('tee_FX2', 'prime_end', 50)])
    _, _, S = sim.renderer.get_params(['S11', 'S21'])
    S.index.name = 'f_GHz'
    S.to_csv(os.path.join(hc.RESULTS, f'{tag}_S.csv'))
    f = np.asarray(S.index, float)
    s21 = 20 * np.log10(np.abs(S['S21'].values))
    from scipy.signal import find_peaks
    idx, pr = find_peaks(-s21, prominence=0.5)
    print(f'=== {tag} ({(time.time() - t0) / 60:.1f} min) ===')
    for i, p in zip(idx, pr['prominences']):
        print(f'  dip at {f[i]:.4f} GHz: S21 = {s21[i]:.2f} dB (depth {p:.2f} dB)')
    print('min S21', s21.min().round(2), 'max', s21.max().round(2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(f, s21, label='S21')
    ax.plot(f, 20 * np.log10(np.abs(S['S11'].values)), label='S11')
    for r in dp.RESONATORS.values():
        ax.axvline(r['f'], color='gray', ls=':')
    ax.set_xlabel('f (GHz)'); ax.set_ylabel('dB'); ax.legend(); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(hc.RESULTS, f'{tag}_S.png'), dpi=120)
    try:
        sim.close()
    except Exception:
        pass


def combined(tag='interp_comb'):
    hc.SETUP.update(FINE)
    hc.SWEEP.update(SWEEP)
    hc.run(tag)


if __name__ == '__main__':
    {'plain': plain, 'combined': combined}[sys.argv[1] if len(sys.argv) > 1 else 'plain']()
