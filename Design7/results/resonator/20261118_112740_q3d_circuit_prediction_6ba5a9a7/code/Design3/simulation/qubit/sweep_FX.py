"""Q3D tuning of the fluxonium (all three resonators rendered):
Ec_FX -> 1.3 GHz (pad size), readout g -> 150 MHz (FX1 readout pad gap),
and the charge-line coupling.

The defaults of fx() / fx1() are the current layout (ftf_layout.py FX_OPTIONS, FX_PER_QUBIT);
add a set below to sweep around them.

Usage:  python sweep_FX.py <set name>
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import os
import sys

import pandas as pd

HERE = _S('Design3/simulation/qubit')
sys.path.insert(0, _I('Design3/simulation/qubit'))          # simulation/ (q3d_ftf.py)

from q3d_ftf import Q3DRunner, RESULTS


def fx(pad_w='125um', pad_h='50um'):
    """shared FX pad size (both fluxoniums)"""
    return dict(pad_width=pad_w, pad_height=pad_h)


def fx1(ro_sep='19.75um', ro_w='100um', ro_h='40um'):
    """FX1 readout pad (sets g_FX1)"""
    return dict(readout_line_options=dict(make_rol=True, pad_width=ro_w, pad_height=ro_h, pad_sep=ro_sep))


SETS = {
    # readout gap of FX1 around the current 19.75 um
    'fx1_gap': {
        'fx1_gap18': (fx(), fx1(ro_sep='18um')),
        'fx1_gap19.75': (fx(), fx1()),
        'fx1_gap22': (fx(), fx1(ro_sep='22um')),
    },
}


def main(name):
    r = Q3DRunner()
    rows = []
    try:
        for tag, (f, f1) in SETS[name].items():
            res = r.run(tag, fx_options=f, fx1_options=f1, with_resonators=True)
            rows.append({'tag': tag, 'Ec_FX1': res['Ec_FX1 [GHz]'], 'Ec_FX2': res['Ec_FX2 [GHz]'],
                         'Ec_TC': res['Ec_TC [GHz]'],
                         'g_FX1': res['readout_FX1']['g/2pi [MHz]'], 'g_FX2': res['readout_FX2']['g/2pi [MHz]'],
                         'beta_FX1': res['readout_FX1']['beta = Cc/Csigma'],
                         'g_TC': res['readout_TC']['g_transmon/2pi [MHz]'],
                         'C_cl_FX1': res['charge_FX1']['C_line-qubit [fF]'],
                         'T1_cl_FX1_ms': res['charge_FX1']['T1 limit [us]'] / 1e3})
    finally:
        r.close()
        if rows:
            df = pd.DataFrame(rows)
            df.to_csv(_relocate(os.path.join(RESULTS, f'sweep_{name}.csv')), index=False)
            print('\n', df.round(4).to_string())


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'fx1_gap')
