"""Q3D sweep of the FX readout-pad width (both fluxoniums, current gaps) to bring the
readout g from ~160 MHz to ~80 MHz (chi ~ kappa / 2, Design1/notes/zz_method.md).

Usage:  python sweep_readout_pad.py [width_um[:gap_um] ...]   (default 40 50 60; gap = FX1 gap, FX2 +1.45 um)
"""
import os
import sys

import pandas as pd
from qiskit_metal import Dict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # simulation/ (q3d_ftf.py)

from q3d_ftf import Q3DRunner, RESULTS

GAPS = {'FX1': '19.75um', 'FX2': '21.2um'}          # current layout (FX_PER_QUBIT)


def ro(w, q, gap=None):
    sep = GAPS[q] if gap is None else f'{gap + (1.45 if q == "FX2" else 0):g}um'
    return dict(readout_line_options=Dict(make_rol=True, pad_width=f'{w}um', pad_height='40um', pad_sep=sep))


def main(cases):
    r = Q3DRunner()
    rows = []
    try:
        for w, gap in cases:
            tag = f'ropad_w{w}' + ('' if gap is None else f'_g{gap:g}')
            res = r.run(tag, fx1_options=ro(w, 'FX1', gap), fx2_options=ro(w, 'FX2', gap), with_resonators=True)
            rows.append({'pad_w_um': w, 'gap_um': gap, 'Ec_FX1': res['Ec_FX1 [GHz]'], 'Ec_FX2': res['Ec_FX2 [GHz]'],
                         'Ec_TC': res['Ec_TC [GHz]'], 'J_1c': res['J_1c [MHz]'], 'J_12': res['J_12 [MHz]'],
                         'beta_FX1': res['readout_FX1']['beta = Cc/Csigma'],
                         'g_FX1': res['readout_FX1']['g/2pi [MHz]'], 'g_FX2': res['readout_FX2']['g/2pi [MHz]'],
                         'C_line_FX1_fF': res['readout_FX1']['C_line [fF]']})
    finally:
        r.close()
        if rows:
            df = pd.DataFrame(rows)
            df.to_csv(os.path.join(RESULTS, 'sweep_ropad.csv' if len(sys.argv) < 2 else 'sweep_ropad_gap.csv'), index=False)
            print('\n', df.round(4).to_string())


if __name__ == '__main__':
    cs = [(int(a.split(':')[0]), float(a.split(':')[1]) if ':' in a else None) for a in sys.argv[1:]]
    main(cs or [(40, None), (50, None), (60, None)])
