"""Save a run record next to its S-parameters: S11/S21 in dB + the capacitances (latest Q3D csv of cin/cp/cj)
and the main layout numbers, so runs can be compared (see compare_runs.py).
Used by hfss_chan6.py:  run_record.save(tag, d6, filter_only, note)"""
import os, json, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, 'results')


def caps_from_q3d():
    out = {}
    for n in ('cin', 'cp', 'cj'):
        p = os.path.join(RES, f'q3d6_{n}_cmat_fF.csv')
        if os.path.exists(p):
            cm = pd.read_csv(p, index_col=0)
            sig = [c for c in cm.index if 'ground' not in c.lower()]
            out[f'C_{n}_fF'] = round(float(-cm.loc[sig[0], sig[1]]), 3)
            out[f'C_{n}_file_time'] = time.strftime('%Y-%m-%d %H:%M', time.localtime(os.path.getmtime(p)))
    return out


def save(tag, d6, filter_only=False, note=''):
    S = pd.read_csv(os.path.join(RES, f'{tag}_S.csv'), index_col=0)
    s11 = S['S11'].map(complex); s21 = S['S21'].map(complex)
    db = pd.DataFrame({'S11_dB': 20 * np.log10(np.abs(s11)), 'S21_dB': 20 * np.log10(np.abs(s21))}, index=S.index)
    db.index.name = 'f_GHz'; db.to_csv(os.path.join(RES, f'{tag}_dB.csv'))
    meta = dict(tag=tag, date=time.strftime('%Y-%m-%d %H:%M'), filter_only=bool(filter_only), note=note,
                caps=caps_from_q3d(), ENC=d6.ENC, L_FILTER=round(d6.L_FILTER, 4), L_RES_TOTAL=round(d6.L_RES_TOTAL, 4),
                CJ=dict(head_w_um=d6.CJ_HW, head_h_um=d6.CJ_HH, gap_um=d6.CJ_HG), MSP=[d6.MSP_F, d6.MSP_R],
                XM_YM=[d6.XM, d6.YM], XMF_YMF=[d6.XMF, d6.YMF])
    json.dump(meta, open(os.path.join(RES, f'{tag}_meta.json'), 'w'), indent=1, default=str)
    return meta
