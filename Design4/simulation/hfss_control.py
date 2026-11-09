"""Control experiment: one tee + hanger resonator (FX1) on a plain 50 Ohm line, ports on the two
ends of the tee's prime line. A standard hanger must show an S21 dip of ~1-2 MHz width at the
resonator frequency. If it does not, the HFSS setup (ports / tee / open end) is the problem,
not the Purcell-filter chain.

Usage:  python hfss_control.py [name] [f_center_GHz]
"""
import os
import sys
import warnings

warnings.filterwarnings('ignore')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from qiskit_metal.analyses.simulation import ScatteringImpedanceSim  # noqa: E402
import design3_purcell as dp  # noqa: E402

RESULTS = os.path.join(HERE, 'results')


def run(name='FX1', fc=6.5598, half_mhz=10.0, step_mhz=0.25, sweep_type='Discrete', tag=None):
    d = dp.build_design()
    sim = ScatteringImpedanceSim(d, 'hfss')
    for k, v in dict(freq_ghz=fc, max_delta_s=0.02, max_passes=8, min_passes=2, min_converged=1,
                     pct_refinement=30, basis_order=1).items():
        sim.setup[k] = v
    sim.setup.sweep_setup.update(name='Sweep', start_ghz=fc - half_mhz * 1e-3, stop_ghz=fc + half_mhz * 1e-3,
                                 count=int(round(2 * half_mhz / step_mhz)) + 1, type=sweep_type, save_fields=False)
    sim.run(name=f'ctrl_{tag or name}', components=[f'tee_{name}', f'res_{name}'],
            open_terminations=[(f'res_{name}', 'end')],
            port_list=[(f'tee_{name}', 'prime_start', 50), (f'tee_{name}', 'prime_end', 50)])
    _, _, S = sim.renderer.get_params(['S11', 'S21'])
    S.index.name = 'f_GHz'
    S.to_csv(os.path.join(RESULTS, f'ctrl_{tag or name}_S.csv'))
    f = np.asarray(S.index, float)
    s21 = 20 * np.log10(np.abs(S['S21'].values))
    i = int(np.argmin(s21))
    print(f'CONTROL {name}: S21 min {s21[i]:.2f} dB at {f[i]:.4f} GHz, window max {s21.max():.2f} dB, '
          f'depth {s21.max() - s21[i]:.2f} dB')
    try:
        sim.close()
    except Exception:
        pass


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'wide':
        # wide Interpolating sweep 6.2-6.9 GHz: a simple hanger resonance (1-2 MHz) is resolved reliably
        run('FX1', 6.55, half_mhz=350.0, step_mhz=0.5, sweep_type='Interpolating', tag='FX1_wide')
    elif len(sys.argv) > 1 and sys.argv[1] == 'coarse':
        # Discrete 6.50-6.62 GHz, 1 MHz step (dip width ~1.2 MHz from Q3D C_1 = 11.3 fF)
        run('FX1', 6.56, half_mhz=60.0, step_mhz=1.0, sweep_type='Discrete', tag='FX1_coarse')
    else:
        run(sys.argv[1] if len(sys.argv) > 1 else 'FX1', float(sys.argv[2]) if len(sys.argv) > 2 else 6.5598)
