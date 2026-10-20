"""Q3D tuning of the fluxonium (all three resonators rendered):
Ec_FX -> 1.3 GHz (pad size), readout g -> 150 MHz (readout pad height / gap),
and the charge-line coupling.

Usage:  python sweep_FX.py <set name>
Past runs: results/history/6_v4_pocket_readout
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # simulation/ (q3d_ftf.py)

from q3d_ftf import Q3DRunner, RESULTS


def fx(pad_w='130um', pad_h='50um', ro_h='40um', ro_sep='30um', ro_w='100um'):
    return dict(pad_width=pad_w, pad_height=pad_h,
                readout_line_options=dict(make_rol=True, pad_width=ro_w, pad_height=ro_h, pad_sep=ro_sep))


SETS = {
    'fx1': {
        'fx1_sep15': fx(ro_sep='15um'),
        'fx1_sep15_pad120': fx(pad_w='120um', ro_sep='15um'),
        'fx1_sep20_roh60': fx(ro_h='60um', ro_sep='20um'),
    },
}


def main(name):
    r = Q3DRunner()
    rows = []
    try:
        for tag, f in SETS[name].items():
            res = r.run(tag, fx_options=f, with_resonators=True)
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
            df.to_csv(os.path.join(RESULTS, f'sweep_{name}.csv'), index=False)
            print('\n', df.round(4).to_string())


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'fx1')
