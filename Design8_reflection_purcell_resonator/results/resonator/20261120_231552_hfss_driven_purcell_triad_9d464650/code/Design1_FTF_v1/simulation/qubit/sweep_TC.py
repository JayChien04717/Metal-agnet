"""Q3D tuning of the transmon coupler with all three resonators rendered:
Ec_TC -> 0.2 GHz (stem length), g_TC -> 150 MHz (readout pad width),
charge line -> T1 limit ~ 1 ms (charge-line tip position).

The defaults of tc() are the current layout (ftf_layout.py TMON_OPTIONS).

Usage:  python sweep_TC.py <set name>
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

HERE = _S('Design1_FTF_v1/simulation/qubit')
sys.path.insert(0, _I('Design1_FTF_v1/simulation/qubit'))          # simulation/ (q3d_ftf.py)

from q3d_ftf import Q3DRunner, RESULTS


def tc(stem='158um', cl_sep='250um', ro_width='273um', cl_x='120um'):
    return dict(pad_head_length=stem,
                charge_line_options=dict(make_cl=True, cl_x=cl_x, cl_sep=cl_sep),
                readout_line_options=dict(make_rol=True, ro_width=ro_width, ro_height='20um', ro_sep='5um'))


SETS = {
    # charge line right of the stem at x = +120 um: tip distance below the bar etch
    # (latest sweep, results/latest/sweep_tc5.csv; cl_sep = 250 um chosen)
    'tc5': {
        'tc5_x120_cl250': tc(cl_sep='250um'),
        'tc5_x120_cl210': tc(cl_sep='210um'),
        'tc5_x120_cl180': tc(cl_sep='180um'),
    },
}


def main(name):
    r = Q3DRunner()
    rows = []
    try:
        for tag, tm in SETS[name].items():
            res = r.run(tag, tmon_options=tm, with_resonators=True)
            rows.append({'tag': tag, 'Ec_FX1': res['Ec_FX1 [GHz]'], 'Ec_TC': res['Ec_TC [GHz]'],
                         'J_1c': res['J_1c [MHz]'], 'g_TC_transmon': res['readout_TC']['g_transmon/2pi [MHz]'],
                         'g_FX1': res['readout_FX1']['g/2pi [MHz]'],
                         'C_cl_TC': res['charge_TC']['C_line-qubit [fF]'],
                         'T1_cl_TC_us': res['charge_TC']['T1 limit [us]']})
    finally:
        r.close()
        if rows:
            df = pd.DataFrame(rows)
            df.to_csv(_relocate(os.path.join(RESULTS, f'sweep_{name}.csv')), index=False)
            print('\n', df.round(3).to_string())


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'tc5')
