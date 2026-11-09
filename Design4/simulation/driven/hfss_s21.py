"""HFSS driven-modal S-parameters of the transmission Purcell filter.

Renders cap_in + filter + cap_out of purcell_layout.build_design(); 50 Ohm lumped ports on
the outer CPW pins of the two caps (cap_in.a, cap_out.b), so the launch pads and long
feed lines are left out (they are matched 50.9 Ohm CPW). Reports f0 (peak |S21|), the
3 dB bandwidth (= kappa_F), peak |S21| and |S11| at f0.

Usage:  python hfss_s21.py [tag] [overlap_um]
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
import purcell_layout as pl  # noqa: E402

RESULTS = os.path.join(HERE, 'results')
os.makedirs(RESULTS, exist_ok=True)
SETUP = dict(freq_ghz=6.7, max_delta_s=0.02, max_passes=10, min_passes=2, min_converged=1,
             pct_refinement=30, basis_order=1)
SWEEP = dict(name='Sweep', start_ghz=5.7, stop_ghz=7.7, count=201, type='Interpolating',
             save_fields=False)


def run(tag='s21', overlap_um=None):
    d = pl.build_design(overlap_um=overlap_um)
    sim = ScatteringImpedanceSim(d, 'hfss')
    for k, v in SETUP.items():
        sim.setup[k] = v
    for k, v in SWEEP.items():
        sim.setup.sweep_setup[k] = v
    t0 = time.time()
    sim.run(name=tag, components=['cap_in', 'filter', 'cap_out'],
            open_terminations=[],
            port_list=[('cap_in', 'a', 50), ('cap_out', 'b', 50)])
    freqs, _, S = sim.renderer.get_params(['S11', 'S21'])
    S = S.copy()
    S.index.name = 'f_GHz'
    S.to_csv(os.path.join(RESULTS, f'{tag}_S.csv'))
    f = np.asarray(S.index, float)
    s21 = 20 * np.log10(np.abs(S['S21'].values))
    s11 = 20 * np.log10(np.abs(S['S11'].values))
    i0 = int(np.argmax(s21))
    above = np.where(s21 >= s21[i0] - 3)[0]
    bw = (f[above[-1]] - f[above[0]]) * 1e3
    print(f'=== {tag} ({(time.time() - t0) / 60:.1f} min) ===')
    print(f'f0 = {f[i0]:.4f} GHz (target {pl.F0}), peak S21 = {s21[i0]:.2f} dB, '
          f'S11 at f0 = {s11[i0]:.1f} dB, 3 dB bandwidth = {bw:.0f} MHz (target {pl.KAPPA_F})')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(f, s21, label='S21'); ax.plot(f, s11, label='S11')
    ax.axvline(pl.F0, color='gray', ls=':'); ax.set_ylim(-60, 2)
    ax.set_xlabel('f (GHz)'); ax.set_ylabel('dB'); ax.legend(); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(RESULTS, f'{tag}_S.png'), dpi=120)
    try:
        sim.close()
    except Exception:
        pass
    return f[i0], bw


if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 's21',
        float(sys.argv[2]) if len(sys.argv) > 2 else None)
