"""Q3D of one CoupledLineTee alone: feedline-resonator coupling capacitance Cc vs
coupling length, so the tee lengths can be set for the target kappa.

Kappa from Cc uses the lab formula (FXTFX24 cell 42): Qc = pi / (Z0 Zr (Cc w)^2).

Usage:  python q3d_tee.py 0.3 0.6 0.9        (coupling lengths in mm)
"""
import os
import sys
import warnings

warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
from qiskit_metal import designs
from qiskit_metal.qlibrary.couplers.coupled_line_tee import CoupledLineTee
from qiskit_metal.analyses.quantization import LOManalysis

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from ftf_layout import TEE_DOWN, FILLET, CPW_WIDTH, CPW_GAP  # noqa: E402
from q3d_ftf import SETUP, RESULTS  # noqa: E402


def tee_cc(coupling_length_mm):
    d = designs.DesignPlanar()
    d.overwrite_enabled = True
    d.variables['cpw_width'] = CPW_WIDTH
    d.variables['cpw_gap'] = CPW_GAP
    CoupledLineTee(d, 'tee', options=dict(
        pos_x='0mm', pos_y='0mm', prime_width='cpw_width', prime_gap='cpw_gap',
        second_width='cpw_width', second_gap='cpw_gap', coupling_space='3um',
        coupling_length=f'{coupling_length_mm}mm', down_length=f'{TEE_DOWN}mm',
        fillet=FILLET, open_termination=True))
    lom = LOManalysis(d, 'q3d')
    for k, v in SETUP.items():
        lom.sim.setup[k] = v
    lom.sim.run(components=['tee'], open_terminations=[('tee', 'prime_start'), ('tee', 'prime_end'),
                                                       ('tee', 'second_end')])
    cmat = lom.sim.capacitance_matrix
    a = [n for n in cmat.index if n.startswith('prime_cpw')][0]
    b = [n for n in cmat.index if n.startswith('second_cpw')][0]
    return -cmat.loc[a, b], cmat


if __name__ == '__main__':
    lengths = [float(x) for x in sys.argv[1:]] or [0.3, 0.6, 0.9]
    rows = []
    for L in lengths:
        cc, cmat = tee_cc(L)
        cmat.to_csv(os.path.join(RESULTS, f'tee_cl{L:.3f}_cmat_fF.csv'))
        rows.append(dict(coupling_length_mm=L, Cc_fF=cc))
        print(f'coupling_length = {L} mm -> Cc = {cc:.3f} fF')
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, 'tee_sweep.csv'), index=False)
    p = np.polyfit(df.coupling_length_mm, df.Cc_fF, 1)
    print(df.to_string(index=False))
    print(f'linear fit: Cc = {p[0]:.3f} fF/mm * L + {p[1]:.3f} fF')
