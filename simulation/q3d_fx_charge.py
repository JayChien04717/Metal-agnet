"""Q3D sweep of the fluxonium charge-line tip position (cl_sep: distance of the open
end below the pocket edge; negative = the stub reaches into the pocket).

Target: ~3x the v3 coupling (0.18 fF, beta 0.0047) with the T1 limit ~1 ms
(f01 = 0.5 GHz, |n01| = 0.1 assumed, as in cap_matrix.charge_line).

Usage:  python q3d_fx_charge.py 5 -10 -30      (cl_sep values in um)
"""
import os
import sys

import pandas as pd
from q3d_ftf import Q3DRunner, RESULTS

if __name__ == '__main__':
    seps = [float(x) for x in sys.argv[1:]] or [5, -10, -30]
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
            df.to_csv(os.path.join(RESULTS, 'sweep_fx_charge.csv'), index=False)
            print('\n', df.round(4).to_string())
