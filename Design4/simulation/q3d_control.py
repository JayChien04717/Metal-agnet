"""Q3D coupling check of one tee + resonator (control of hfss_control.py): C_1 = -C[resonator, prime line]
and the lab formula kappa = f / Q_c, Q_c = pi / (Z0 Zr (C_1 w)^2).  Design3 got C_1 ~ 12 fF.

Usage:  python q3d_control.py [name]
"""
import os
import sys
import warnings

warnings.filterwarnings('ignore')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from qiskit_metal.analyses.quantization import LOManalysis  # noqa: E402
import design3_purcell as dp  # noqa: E402
import purcell_layout as pl  # noqa: E402

RESULTS = os.path.join(HERE, 'results')
SETUP = dict(max_passes=12, min_passes=2, percent_error=0.2, freq_ghz=6.7,
             solution_order='Medium', auto_increase_solution_order=False)


def run(name='FX1'):
    d = dp.build_design()
    lom = LOManalysis(d, 'q3d')
    for k, v in SETUP.items():
        lom.sim.setup[k] = v
    lom.sim.run(name=f'q3dctrl_{name}', components=[f'tee_{name}', f'res_{name}'],
                open_terminations=[(f'res_{name}', 'end'), (f'tee_{name}', 'prime_start'),
                                   (f'tee_{name}', 'prime_end')])
    cm = lom.sim.capacitance_matrix
    cm.to_csv(os.path.join(RESULTS, f'q3dctrl_{name}_cmat_fF.csv'))
    print(cm.round(3).to_string())
    nets = [n for n in cm.index if 'ground' not in n.lower()]
    print('NETS', nets)
    f = dp.RESONATORS[name]['f']
    for a in nets:
        for b in nets:
            if a < b:
                C1 = -cm.loc[a, b]
                w = 2 * np.pi * f * 1e9
                Qc = np.pi / (pl.CPW_Z0 ** 2 * (C1 * 1e-15 * w) ** 2)
                print(f'C({a},{b}) = {C1:.3f} fF -> Q_c {Qc:.0f}, kappa {f * 1e3 / Qc:.2f} MHz '
                      f'(design target kappa {dp.RESONATORS[name]["kappa"]} MHz, Cc {dp.RESONATORS[name]["Cc"]:.1f} fF)')
    try:
        lom.sim.close()
    except Exception:
        pass


if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'FX1')
