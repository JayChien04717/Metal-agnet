"""Q3D of one CoupledLineTee (10/6 um CPW, 3 um coupling space) vs coupling length -> Cc(L) fit and
the length for kappa = 3 MHz (lab formula Qc = pi / (Z0 Zr (Cc w)^2)). Uses Design3's tee_cc();
results are written to Design4 only.

Usage: python q3d_tee_sweep.py 1.0 1.3 1.6
"""
import os
import sys
import warnings

warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # simulation/ (this script: simulation/<kind>/)
import sys as _sys; _sys.path[:0] = [os.path.join(HERE, _k) for _k in ('driven', 'eigenmode', 'q3d')] + [HERE]
D3 = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'Design3')
sys.path.insert(0, os.path.join(D3, 'pattern design'))
sys.path.insert(0, os.path.join(D3, 'simulation', 'q3d'))
from q3d_tee import tee_cc  # noqa: E402
import ftf_layout as f3  # noqa: E402

RESULTS = os.path.join(HERE, 'results')

if __name__ == '__main__':
    rows = []
    for L in [float(x) for x in sys.argv[1:]] or [1.0, 1.3, 1.6]:
        cc, cmat = tee_cc(L)
        cmat.to_csv(os.path.join(RESULTS, f'tee3_cl{L:.3f}_cmat_fF.csv'))
        rows.append(dict(coupling_length_mm=L, Cc_fF=cc))
        print(f'coupling_length = {L} mm -> Cc = {cc:.3f} fF', flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, 'tee3_sweep.csv'), index=False)
    p = np.polyfit(df.coupling_length_mm, df.Cc_fF, 1)
    print(df.to_string(index=False))
    print(f'linear fit: Cc = {p[0]:.3f} fF/mm * L + {p[1]:.3f} fF')
    for k, f in (('FX1', 6.6), ('TC', 6.7), ('FX2', 6.8)):
        cc = f3.coupling_cap(f, 3.0)
        print(f'{k}: kappa 3 MHz needs Cc {cc:.2f} fF -> coupling length {(cc - p[1]) / p[0]:.3f} mm')
