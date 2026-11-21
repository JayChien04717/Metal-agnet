"""Re-examination of the Wallraff thesis Table 2.4 reproduction and its transfer to Design6_TXT (re-analysis, no EM rerun).

1. Closed-form J (code/general/purcell_j.py) for the five Table 2.4 channels vs the table's simulated / measured J
   and vs the old full-circuit reproduction (Design6_TXT/simulation/wallraff_model.py, git 2eb2bef, deleted in fafa3ac;
   its output JSON is passed in as an input snapshot).
2. Per-channel errors of that old reproduction in omega_p, omega_r, omega_p - omega_r and kappa_p.
3. Closed-form J for Design6_TXT r9 (+0.1 mm, no claw) with the Q3D hammer, vs readout_model and HFSS; the lumped C_j that
   would be needed to reproduce the HFSS J.
Wallraff CPW: v = 1.147e8 m/s (old global calibration), Z0 = 50 ohm. C_AB = 0."""
from argparse import ArgumentParser
from pathlib import Path
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code/general'))
from run_record import create_run, finish_run  # noqa: E402
from purcell_j import j_closed_form            # noqa: E402

Q = ['Q36', 'Q35', 'Q34', 'Q33', 'Q32']
TAB = dict(Cpf=[48.57, 45.47, 37.08, 45.10, 37.55], Cpfg=[98.80, 93.33, 78.11, 92.56, 79.00],
           Cpr=[2.50, 2.77, 2.42, 2.5, 2.66], Cqr=[12.26, 17.95, 12.26, 17.95, 12.26],
           Cqrg=[17.61, 27.03, 17.61, 27.04, 17.61], Cqg=[105.36, 99.35, 105.36, 99.35, 105.36],
           Lpf=[1200] * 5, Lpfg=[1920.45, 2268.05, 2202.36, 1844.14, 2296.28],
           Lrr=[1303.6, 1300, 1303.6, 1300, 1303.6], Lrrg=[2489.65, 2749.33, 2603.73, 2316.18, 2720.55])
SIM = dict(wp=[7250, 6690, 7050, 7440, 6850], wr=[7238, 6671, 7039, 7417, 6839], J=[15.84, 15.04, 15.80, 15.55, 15.41],
           kp=[35.37, 40.98, 52.07, 63.16, 59.91])
MEAS = dict(wp=[7258, 6705, 7056, 7451, 6865], wr=[7227, 6655, 7018, 7394, 6821], J=[20.54, 18.20, 20.59, 16.27, 17.63],
            kp=[55.58, 28.43, 53.61, 47.01, 60.41])


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('old_json', type=Path, help='wallraff_model.json from git 2eb2bef')
    ap.add_argument('poles', type=Path, help='Design6_TXT readout_poles run (+0.1 mm)')
    ap.add_argument('model', type=Path, help='Design6_TXT readout_model run')
    a = ap.parse_args()
    old_p, po, mo = (p.resolve() for p in (a.old_json, a.poles, a.model))
    folder = create_run('Design6_TXT', 'resonator', 'wallraff_check', entrypoint=__file__,
                        inputs=[old_p, po / 'summary.json', mo / 'run.json', mo / 'model.json'], upstream=[str(po), str(mo)],
                        parameters=dict(EM_rerun=False, wallraff_v_m_s=1.147e8, wallraff_Z0=50.0, C_AB_fF=0.0,
                                        old_reproduction='Design6_TXT/simulation/wallraff_model.py at git 2eb2bef'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        old = json.loads(old_p.read_text())
        cl = 1e12 / (50 * 1.147e8)          # fF/mm
        wall = {}
        for i, q in enumerate(Q):
            end_p = TAB['Cpf'][i] + TAB['Cpfg'][i]                       # feed node ~ AC ground (25 ohm << 1/wC_pf)
            end_r = TAB['Cqrg'][i] + 1 / (1 / TAB['Cqr'][i] + 1 / TAB['Cqg'][i])
            f = (SIM['wp'][i] + SIM['wr'][i]) / 2
            cf = j_closed_form(f, TAB['Cpr'][i], cl, TAB['Lpf'][i] / 1000, (TAB['Lpf'][i] + TAB['Lpfg'][i]) / 1000, end_p,
                               TAB['Lrr'][i] / 1000, (TAB['Lrr'][i] + TAB['Lrrg'][i]) / 1000, end_r)
            o = old['table'][q]
            wall[q] = dict(J_closed_form=cf['J_MHz'], J_old_circuit=o['J'], J_table_sim=SIM['J'][i], J_table_meas=MEAS['J'][i],
                           J_meas_over_sim=MEAS['J'][i] / SIM['J'][i],
                           wp_minus_wr=dict(table_sim=SIM['wp'][i] - SIM['wr'][i], table_meas=MEAS['wp'][i] - MEAS['wr'][i],
                                            old_circuit=o['wp'] - o['wr']),
                           old_error_MHz=dict(wp=o['wp'] - SIM['wp'][i], wr=o['wr'] - SIM['wr'][i]),
                           kappa_p=dict(table_sim=SIM['kp'][i], table_meas=MEAS['kp'][i], old_circuit=o['kp'],
                                        old_rel_error=o['kp'] / SIM['kp'][i] - 1),
                           closed_form_detail=cf)
        mrun = json.loads((mo / 'run.json').read_text(encoding='utf-8'))['parameters']
        mjs = json.loads((mo / 'model.json').read_text(encoding='utf-8'))['prediction_current']['bare']
        hp = json.loads((po / 'summary.json').read_text(encoding='utf-8'))['channels']
        geo = mrun['geometry']; c6 = 1e12 / (geo['z0_ohm'] * geo['velocity_m_s'])
        cj = mrun['hammer_deembedded_fF'][0]; cp = mrun['cp_fF']
        d6 = {}
        for q in ('Q1', 'Q2'):
            Lf, Lr = mrun['lengths_mm'][q]; h = hp[q]['wide_only']
            f = (h['bare_fp_GHz'] + h['bare_fr_GHz']) / 2 * 1000
            cf = j_closed_form(f, cj, c6, 1.20, Lf, cp[0] + cp[2], 1.30, Lr, 0.0)
            d6[q] = dict(J_closed_form=cf['J_MHz'], J_readout_model=mjs[q]['J_MHz'], J_hfss=h['J_MHz'],
                         hfss_over_closed_form=h['J_MHz'] / cf['J_MHz'],
                         equivalent_lumped_Cj_fF=cj * h['J_MHz'] / cf['J_MHz'], q3d_Cj_fF=cj, closed_form_detail=cf)
        out = dict(wallraff=wall, design6=d6, units='MHz, fF, mm',
                   findings=['closed-form J reproduces Table 2.4 J_sim within 5 %; the Wallraff measured J is 5-30 % above J_sim',
                             'old full-circuit reproduction: J 4-10 % low, but wp - wr is 4.1-6.2x the table value '
                             '(+45 to +108 vs +11 to +23 MHz), wr 20-29 MHz low and wp 5-64 MHz high, kappa_p -55 to +45 %',
                             'Design6_TXT: closed form equals readout_model; HFSS J needs a lumped C_j of ~2.2 fF, 3.3x the Q3D hammer'])
        (folder / 'check.json').write_text(json.dumps(out, indent=1))
        finish_run(folder, status='completed', summary=out)
        for q, r in wall.items():
            print(f"{q}: J cf {r['J_closed_form']:.2f} old {r['J_old_circuit']:.2f} sim {r['J_table_sim']} meas {r['J_table_meas']} | "
                  f"wp-wr table {r['wp_minus_wr']['table_sim']:+d} old {r['wp_minus_wr']['old_circuit']:+.1f} | "
                  f"kp table {r['kappa_p']['table_sim']} old {r['kappa_p']['old_circuit']:.1f}")
        for q, r in d6.items():
            print(f"Design6_TXT {q}: J cf {r['J_closed_form']:.2f} model {r['J_readout_model']:.2f} HFSS {r['J_hfss']:.2f} "
                  f"-> C_j equiv {r['equivalent_lumped_Cj_fF']:.2f} fF ({r['hfss_over_closed_form']:.2f}x)")
    except Exception:
        import traceback
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()}); raise


if __name__ == '__main__':
    main()
