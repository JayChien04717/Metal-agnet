# -*- coding: utf-8 -*-
r"""Q3D-only channel solve (no HFSS offsets): filter / resonator lengths and hammer head_gap (2026-11-26).

Inputs (all Q3D):
  --check    results/resonator/<q3d_model_check run>  in-situ C_j, end loads, tap_r of the CURRENT build.py channels
                                                      (full-fixture Q3D) and the calibrated line (straight-CPW Q3D)
  --coupons  existing hammer coupon runs at the same spacing / head_height (shape of C_j(head_gap) only; the in-situ
             C_j sets the scale: C_j_insitu(g) = C_j_insitu(g0) * C_coupon(g) / C_coupon(g0)). No new coupon runs.
Targets: f_p = f_r = PURCELL_CHANNELS f_GHz, J = sqrt(kappa_r kappa_p) / 2 (kappa_r 3, kappa_p 30 MHz).
kappa_p: the C_p plates are kept (driven S21 gives 28.6 / 30.6 MHz at the r16 plates).
Model as in q3d_model_check.py: tapped_lambda4_f / j_hammer with the Q3D-calibrated line; the stub load at the tap
is C_j, the remaining in-situ excess C sits at the open end (unchanged when only the shorted arm length changes).
Known model - eigenmode residuals at r16 (eigen pass 2, approximate): f_p -29 / +11 MHz, f_r -15 / +10 MHz (Q1 / Q2),
not applied here (pure Q3D prediction).

    python -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/resonator/q3d_length_solve.py --check <run> --coupons <run>[,<run>]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator'),
                str(ROOT / 'Design11_Transmon2Q_Manhattan' / 'layout')]

import numpy as np

import build as D11
import purcell_readout_model as PRM
from run_record import create_run, finish_run

DESIGN = 'Design11_Transmon2Q_Manhattan'
KAPPA_R, KAPPA_P = 3.0, 30.0
TAP_F = 1.20


def coupon_curve(runs, head_height_um):
    """(gap_um, C_coupon_fF) from existing cj_g* coupon runs with the build.py spacing and this head height."""
    pts = {}
    for run in runs:
        prm = json.loads((run / 'run.json').read_text(encoding='utf-8')).get('parameters', {})
        hm = prm.get('hammer', {})
        if float(str(hm.get('line_spacing', '0um')).replace('um', '')) != D11.HAMMER_LINE_SPACING_UM:
            continue
        if float(str(hm.get('head_height', '0um')).replace('um', '')) != head_height_um:
            continue
        for f in (run / 'raw').glob('cj_g*/caps.json'):
            r = json.loads(f.read_text(encoding='utf-8'))
            pts[float(f.parent.name[4:].replace('p', '.'))] = float(r['mutual_fF'])
    g = np.array(sorted(pts))
    return g, np.array([pts[x] for x in g])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', required=True)
    ap.add_argument('--coupons', required=True)
    a = ap.parse_args()
    chk_run = ROOT / a.check
    chk = json.loads((chk_run / 'summaries' / 'summary.json').read_text(encoding='utf-8'))
    cpw = chk['line']
    J_t = float(np.sqrt(KAPPA_R * KAPPA_P) / 2)
    stub = D11.HAMMER_LINE_SPACING_UM / 2000.0
    runs = [ROOT / c.strip() for c in a.coupons.split(',') if c.strip()]
    out = dict(targets=dict(J_MHz=J_t, kappa_r_MHz=KAPPA_R, kappa_p_MHz=KAPPA_P), line=cpw, channels={})
    for q, c in chk['channels'].items():
        ch = D11.PURCELL_CHANNELS[q]
        f = ch['f_GHz']
        g, cc = coupon_curve(runs, ch['head_height_um'])
        logc = lambda x: float(np.interp(x, g, np.log(cc)))                       # log-linear in gap
        gap_for = lambda c_rel: float(np.interp(np.log(c_rel), np.log(cc[::-1]), g[::-1]))
        g0, cj0 = ch['head_gap_um'], c['C_j_insitu_fF']
        scale = cj0 / np.exp(logc(g0))                                          # in-situ / coupon at the current gap
        end_f, end_r, tap_r = c['end_load_filter_fF'], c['end_load_res_fF'], c['tap_r_mm']
        cj = cj0
        for _ in range(50):
            Lf = PRM.length_for_f(f, TAP_F, end_f, stub, cj, cpw=cpw)
            Lr = PRM.length_for_f(f, tap_r, end_r, stub, cj, cpw=cpw)
            cj_new = J_t / PRM.j_hammer(f, 1.0, Lf, TAP_F, end_f, Lr, tap_r, end_r, cpw=cpw)
            if abs(cj_new - cj) < 1e-7:
                break
            cj = cj_new
        gap = gap_for(cj / scale)
        out['channels'][q] = dict(
            f_GHz=f, L_filter_mm=float(Lf), L_res_mm=float(Lr), head_gap_um=gap, head_height_um=ch['head_height_um'],
            C_j_target_insitu_fF=float(cj), C_j_r16_insitu_fF=cj0, insitu_over_coupon=float(scale),
            gap_in_coupon_range=bool(g.min() <= gap <= g.max()), cp_plate_um=ch['cp_plate_um'],
            J_model_MHz=PRM.j_hammer(f, cj, Lf, TAP_F, end_f, Lr, tap_r, end_r, cpw=cpw),
            check_f_p_GHz=PRM.tapped_lambda4_f(Lf, TAP_F, end_f, stub, cj, cpw=cpw),
            check_f_r_GHz=PRM.tapped_lambda4_f(Lr, tap_r, end_r, stub, cj, cpw=cpw),
            r16=dict(L_filter_mm=ch['L_filter'], L_res_mm=ch['L_res'], head_gap_um=g0),
            delta_um=dict(L_filter=(Lf - ch['L_filter']) * 1e3, L_res=(Lr - ch['L_res']) * 1e3, head_gap=gap - g0))
    folder = create_run(DESIGN, 'resonator', 'q3d_length_solve', entrypoint=__file__,
                        inputs=[ROOT / DESIGN / 'layout' / 'build.py'],
                        parameters=dict(check_run=a.check, coupon_runs=a.coupons.split(','), targets=out['targets']))
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    (folder / 'summaries' / 'summary.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    finish_run(folder, status='completed', summary=out)
    print('RUN_FOLDER=' + str(folder))
    print(json.dumps(out['channels'], indent=1))


if __name__ == '__main__':
    main()
