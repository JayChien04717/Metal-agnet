"""Q3D sweep over FX pad / coupler geometry (one Ansys session).

Each variant: tag -> (fx_options, tmon_options[, fx_x in mm]). Results land in results/,
and a table of each set is written to results/sweep_<set>.csv.

Usage:  python q3d_sweep.py <name of a VARIANT SET below>
"""
import os
import sys

import pandas as pd

from q3d_ftf import Q3DRunner, RESULTS

ro_small = dict(make_rol=True, pad_width='120um', pad_height='30um', pad_sep='40um')
fx130 = dict(pad_width='130um', pad_height='50um', pocket_width='300um', lead_x='-40um')
tc_small = dict(pad_head_length='250um', pad_arm_width='15um',
                charge_line_options=dict(make_cl=True, cl_x='-100um'))

VARIANT_SETS = {
    # step 1: sensitivity of Ec to pad size, readout pad and coupler size
    'step1': {
        's1_pad160x60': (dict(pad_width='160um', pad_height='60um'), {}),
        's1_pad120x50': (dict(pad_width='120um', pad_height='50um'), {}),
        's1_pad120x50_ro': (dict(pad_width='120um', pad_height='50um', readout_line_options=ro_small), {}),
        's1_tc_short': ({}, dict(pad_head_length='250um', pad_arm_width='15um')),
    },
    # step 2: smaller FX pocket and FX closer together -> shorter coupler bar
    'step2': {
        's2_x035': (fx130, tc_small, 0.35),
        's2_x030': (fx130, tc_small, 0.30),
    },
    # step 3: interpolated FX_X for Ec_TC = 0.2 GHz
    'step3': {
        's3_x034': (fx130, tc_small, 0.34),
    },
}


def main(set_name):
    runner = Q3DRunner()
    rows = []
    try:
        for tag, v in VARIANT_SETS[set_name].items():
            fx, tm, fx_x = (v + (None,))[:3]
            res = runner.run(tag, fx_options=fx, tmon_options=tm, fx_x=fx_x)
            cs = res.pop('C_sigma [fF]')
            rows.append(dict(tag=tag, fx=str(fx), tmon=str(tm), fx_x=fx_x, **res,
                             **{f'Csig_{k} [fF]': v for k, v in cs.items()}))
    finally:
        runner.close()
        if rows:
            df = pd.DataFrame(rows)
            df.to_csv(os.path.join(RESULTS, f'sweep_{set_name}.csv'), index=False)
            print('\n', df.drop(columns=['fx', 'tmon']).round(3).to_string())


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'step1')
