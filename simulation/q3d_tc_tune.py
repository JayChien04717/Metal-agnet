"""Q3D tuning of the transmon coupler with all three resonators rendered:
Ec_TC -> 0.2 GHz (stem length), g_TC -> 150 MHz (readout pad width),
charge line -> T1 limit >~ 1 ms (charge-line tip position).

Usage:  python q3d_tc_tune.py <set name>
"""
import os
import sys

import pandas as pd
from q3d_ftf import Q3DRunner, RESULTS


def tc(stem, cl_sep, ro_width='300um'):
    return dict(pad_head_length=stem,
                charge_line_options=dict(make_cl=True, cl_x='-100um', cl_sep=cl_sep),
                readout_line_options=dict(make_rol=True, ro_width=ro_width, ro_height='20um', ro_sep='5um'))


SETS = {
    'tc1': {
        'tc1_stem150_cl250': tc('150um', '250um'),
        'tc1_stem120_cl250_ro270': tc('120um', '250um', '270um'),
    },
    'tc2': {
        'tc2_stem150_ro265': tc('150um', '250um', '265um'),
        'tc2_stem160_ro260': tc('160um', '250um', '260um'),
    },
    'tc3': {
        'tc3_stem170_ro272': tc('170um', '250um', '272um'),
        'tc3_stem170_ro275': tc('170um', '250um', '275um'),
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
            df.to_csv(os.path.join(RESULTS, f'sweep_{name}.csv'), index=False)
            print('\n', df.round(3).to_string())


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'tc1')
