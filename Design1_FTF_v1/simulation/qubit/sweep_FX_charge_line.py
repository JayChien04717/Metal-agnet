"""Q3D sweep of the fluxonium charge-line tip position (cl_sep: distance of the open
end below the pocket edge; negative = the stub reaches into the pocket).

Current layout: cl_sep = 15 um -> 0.45 fF, T1 limit 2.0 ms
(f01 = 0.5 GHz, |n01| = 0.1 assumed, as in cap_matrix.charge_line).

Usage:  python sweep_FX_charge_line.py 25 15 5      (cl_sep values in um)
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

if __name__ == '__main__':
    seps = [float(x) for x in sys.argv[1:]] or [25, 15, 5]
    r = Q3DRunner()
    rows = []
    try:
        for s in seps:
            tag = f'fxcl_sep{s:g}'
            fx = dict(charge_line_options=dict(make_cl=True, cl_sep=f'{s:g}um', cl_length='100um'))
            res = r.run(tag, fx_options=fx, with_resonators=True)
            c = res['charge_FX1']
            rows.append({'cl_sep_um': s, 'C_cl [fF]': c['C_line-qubit [fF]'], 'beta': c['beta'],
                         'T1_limit_ms': c['T1 limit [us]'] / 1e3,
                         'P_10MHz_dBm': c['P for 10 MHz Rabi [dBm at chip]'],
                         'Ec_FX1': res['Ec_FX1 [GHz]']})
    finally:
        r.close()
        if rows:
            df = pd.DataFrame(rows)
            df.to_csv(_relocate(os.path.join(RESULTS, 'sweep_fx_charge.csv')), index=False)
            print('\n', df.round(4).to_string())
