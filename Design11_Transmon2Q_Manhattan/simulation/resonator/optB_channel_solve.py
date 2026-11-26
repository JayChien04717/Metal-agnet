# -*- coding: utf-8 -*-
r"""Readout option B channel solve (circuit model only, no EM): hammer gap, C_p plate width and lengths.

    python -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/resonator/optB_channel_solve.py \
        --core <results/qubit/<q3d_core run>> --coupons <run>[,<run> ...] [--kappa-p 30] [--anchor '{"Q1": {...}}' --anchor-run <run>]

Targets (option B, user-selected 2026-11-25, E_C fixed): kappa_r = |2 chi| = 3.0 MHz, kappa_p = --kappa-p
(20 MHz until r9, 30 MHz from 2026-11-26) -> J = sqrt(kappa_r kappa_p) / 2 at f_p = f_r = f_target.
Inputs: Q3D core run (claw end load = readout-net self C with Q1, Q2, TC islands floating, Kron-reduced),
Q3D coupon run (hammer C_j and head excess vs head_gap; EnclosingCap C_p mutual / plate ground vs plate_width),
full-feed kappa_p model (feed_kappa_model.py) for the plate width that gives kappa_p = 20 MHz.
Resonator tap = resonator_tap + resonator_open_extension (PurcellReadout docstring), filter tap 1.20 mm,
hammer stub = line_spacing / 2, stub load = C_j + head excess (skill/purcell_resonator_design.md section 1).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator'),
                str(ROOT / 'Design11_Transmon2Q_Manhattan' / 'layout'),
                str(Path(__file__).resolve().parent)]
import project_paths
project_paths.configure(__file__)

import numpy as np
import pandas as pd
from scipy.optimize import brentq

import build as D11
import feed_kappa_model as FK
import purcell_readout_model as PRM
import purcell_resonator_pair_sim as S
from run_record import create_run, finish_run

DESIGN = 'Design11_Transmon2Q_Manhattan'
KAPPA_R, KAPPA_P = 3.0, 20.0
J_TARGET = float(np.sqrt(KAPPA_R * KAPPA_P) / 2)
TAP_F, STUB = 1.20, 0.200


def claw_end_loads(core_run: Path):
    cm = pd.read_csv(core_run / 'raw' / 'cmat_fF.csv', index_col=0)
    names = [n for n in cm.index if 'ground' not in n.lower()]
    isl = [n for n in names if n.lower().startswith(('cross_', 'tmon_'))]
    M = cm.loc[names, names].to_numpy(float)
    out = {}
    for q in ('Q1', 'Q2'):
        ro = next(n for n in names if 'readout' in n.lower() and q.lower() in n.lower())
        out[q] = float(S.kron_floating(M, [names.index(ro)], [names.index(n) for n in isl])[0, 0])
    return out, isl


def coupon_tables(cp_runs):
    """No parallel-line reference coupon (user decision 2026-11-25): C_j = -Maxwell(filter, resonator) of the hammer
    coupon directly; head excess = net ground C - analytic straight 10/6 CPW of the same main-line length
    (2 x coupon_lead, two etched open ends) - C' x stub (s/2), i.e. relative to the line model length_for_f uses.
    Several coupon runs (same coupon geometry and lead) are merged; a later run overrides the same point."""
    cp_runs = [cp_runs] if isinstance(cp_runs, Path) else list(cp_runs)
    raw, leads = {}, set()
    for run in cp_runs:
        raw.update({f.parent.name: json.loads(f.read_text(encoding='utf-8')) for f in (run / 'raw').glob('*/caps.json')})
        leads.add(json.loads((run / 'run.json').read_text(encoding='utf-8')).get('parameters', {}).get('coupon_lead_um', 400.0))
    if len(leads) != 1:
        raise RuntimeError(f'coupon runs use different coupon_lead_um {sorted(leads)}')
    cl = PRM.line('10/6')['c_line']
    lead_mm = leads.pop() / 1000
    straight = cl * 2 * lead_mm + 2 * PRM.open_end_c('10/6')
    cj = []
    for k, r in raw.items():
        if k.startswith('cj_g'):
            ex = np.mean([r['ground_fF'][a] - straight - cl * STUB for a in ('filter', 'resonator')])
            cj.append((float(k[4:].replace('p', '.')), r['mutual_fF'], float(ex)))
    cp = sorted((float(k[4:].replace('p', '.')), r['mutual_fF'], r['ground_fF']['plate'], r['ground_fF']['frame'])
                for k, r in raw.items() if k.startswith('cp_w'))
    return np.array(sorted(cj)), np.array(cp)


def main():
    global KAPPA_P, J_TARGET
    ap = argparse.ArgumentParser()
    ap.add_argument('--core', required=True)
    ap.add_argument('--coupons', required=True, help='q3d_coupons run(s), comma-separated (merged)')
    ap.add_argument('--overrides', default='{}', help='build_design overrides of the core geometry (for the extension)')
    ap.add_argument('--kappa-p', type=float, default=KAPPA_P, help='target kappa_p (MHz); J = sqrt(kappa_r kappa_p) / 2')
    ap.add_argument('--anchor', default='',
                    help='JSON {Q1: {f_p_GHz, f_r_GHz}, Q2: ...}: HFSS bare frequencies of the CURRENT build.py channels; '
                         'the model-HFSS offset of each line is removed from the new lengths (note/purcell_filter_method.md: '
                         'anchor trims to HFSS)')
    ap.add_argument('--anchor-run', default='', help='HFSS run the anchor values come from (recorded only)')
    ap.add_argument('--kp-correction', default='',
                    help='JSON {Q1: r, Q2: r}: HFSS-fitted kappa_p / feed-model kappa_p of the CURRENT build.py plates '
                         '(fit_s21_purcell.py); the feed model is then solved for kappa_p / r')
    a = ap.parse_args()
    KAPPA_P = float(a.kappa_p)
    J_TARGET = float(np.sqrt(KAPPA_R * KAPPA_P) / 2)
    anchor = json.loads(a.anchor) if a.anchor else {}
    kp_corr = json.loads(a.kp_correction) if a.kp_correction else {}
    core, cpr = ROOT / a.core, [ROOT / c.strip() for c in a.coupons.split(',') if c.strip()]
    ce, isl = claw_end_loads(core)
    cj_t, cp_t = coupon_tables(cpr)
    gap_for = lambda c: float(np.exp(np.interp(np.log(c), np.log(cj_t[::-1, 1]), np.log(cj_t[::-1, 0]))))
    ex_for = lambda g: float(np.interp(g, cj_t[:, 0], cj_t[:, 2]))
    if not (cj_t[:, 1].min() <= 1.0 and len(cp_t) >= 2):
        raise RuntimeError('coupon tables incomplete')

    # plate width for kappa_p = 20 MHz from the full-feed model (linear in width between coupon points)
    d = D11.build_design(variant='sim', overrides=json.loads(a.overrides), channels=True, check_keepout_um=None)
    L = FK.feed_lengths(d)
    w_cin = 140.0
    i_cin = int(np.argmin(abs(cp_t[:, 0] - w_cin)))
    c_in = dict(series=cp_t[i_cin, 1], ga=cp_t[i_cin, 3], gb=cp_t[i_cin, 2])
    out = dict(targets=dict(kappa_r_MHz=KAPPA_R, two_chi_MHz=-KAPPA_R, kappa_p_MHz=KAPPA_P, J_MHz=J_TARGET),
               inputs=dict(core_run=str(a.core), coupon_runs=a.coupons.split(','), floating_islands=isl,
                           claw_end_fF=ce, cj_table=cj_t.tolist(), cp_table=cp_t.tolist(), feed_lengths_mm=L),
               channels={})

    def kp_model(i, w):
        mut = np.interp(w, cp_t[:, 0], cp_t[:, 1])
        frame = np.interp(w, cp_t[:, 0], cp_t[:, 3])
        f0 = list(D11.PURCELL_CHANNELS.values())[i]['f_GHz'] * 1e9
        z = FK.node_z(f0, L, c_in, frame, FK.V, L['stub'])[i]
        return PRM.kappa_hanger(f0 / 1e9, mut) * z.real / (FK.Z0 / 2)

    for i, q in enumerate(('Q1', 'Q2')):
        ch = D11.PURCELL_CHANNELS[q]
        f = ch['f_GHz']
        prr = d.components['prr_' + q].options
        tap_r = (d.parse_value(prr.resonator_tap) + d.parse_value(prr.resonator_open_extension))
        kp_target = KAPPA_P / float(kp_corr.get(q, 1.0))          # feed-model target that gives KAPPA_P in HFSS
        if not kp_model(i, cp_t[0, 0]) <= kp_target <= kp_model(i, cp_t[-1, 0]):
            raise RuntimeError(f'{q}: model kappa_p {kp_target:.1f} MHz outside the coupon plate-width table '
                               f'({cp_t[0, 0]:g}-{cp_t[-1, 0]:g} um); run wider C_p coupons')
        w = brentq(lambda w: kp_model(i, w) - kp_target, cp_t[0, 0], cp_t[-1, 0])
        # HFSS anchor: model frequency of the current geometry vs HFSS -> per-line offset
        offs = dict(filter=0.0, resonator=0.0)
        if q in anchor:
            w0 = ch['cp_plate_um']
            ce_p0 = float(np.interp(w0, cp_t[:, 0], cp_t[:, 1]) + np.interp(w0, cp_t[:, 0], cp_t[:, 2]))
            g0 = ch['head_gap_um']
            st0 = float(np.exp(np.interp(np.log(g0), np.log(cj_t[:, 0]), np.log(cj_t[:, 1])))) + ex_for(g0)
            fp0 = PRM.tapped_lambda4_f(ch['L_filter'], TAP_F, ce_p0, STUB, st0)
            fr0 = PRM.tapped_lambda4_f(ch['L_res'], tap_r, ce[q], STUB, st0)
            offs = dict(filter=anchor[q]['f_p_GHz'] - fp0, resonator=anchor[q]['f_r_GHz'] - fr0,
                        model_current_GHz=dict(filter=fp0, resonator=fr0), hfss_GHz=anchor[q])
        cend_p = float(np.interp(w, cp_t[:, 0], cp_t[:, 1]) + np.interp(w, cp_t[:, 0], cp_t[:, 2]))
        cj, gap = 0.6, 42.0
        for _ in range(30):                       # J, lengths and the head load depend on each other
            stub_c = cj + ex_for(gap)
            Lf = PRM.length_for_f(f - offs['filter'], TAP_F, cend_p, STUB, stub_c)
            Lr = PRM.length_for_f(f - offs['resonator'], tap_r, ce[q], STUB, stub_c)
            cj_new = J_TARGET / PRM.j_hammer(f, 1.0, Lf, TAP_F, cend_p, Lr, tap_r, ce[q])
            gap = gap_for(cj_new)
            if abs(cj_new - cj) < 1e-6:
                break
            cj = cj_new
        out['channels'][q] = dict(
            f_GHz=f, plate_width_um=float(w), C_p_end_load_fF=cend_p, kappa_p_feed_model_MHz=float(kp_model(i, w)), kappa_p_hfss_correction=float(kp_corr.get(q, 1.0)),
            kappa_p_expected_hfss_MHz=float(kp_model(i, w) * kp_corr.get(q, 1.0)),
            C_j_fF=float(cj), head_gap_um=float(gap), head_excess_fF=ex_for(gap),
            claw_end_fF=ce[q], tap_r_mm=float(tap_r), tap_f_mm=TAP_F,
            L_filter_mm=float(Lf), L_res_mm=float(Lr),
            J_MHz=PRM.j_hammer(f, cj, Lf, TAP_F, cend_p, Lr, tap_r, ce[q]),
            kappa_r_MHz=PRM.kappa_readout(J_TARGET, KAPPA_P, 0.0),
            gap_in_coupon_range=bool(cj_t[:, 0].min() <= gap <= cj_t[:, 0].max()),
            width_in_coupon_range=bool(cp_t[:, 0].min() <= w <= cp_t[:, 0].max()),
            hfss_anchor_offset_GHz=offs)

    folder = create_run(DESIGN, 'resonator', 'optB_channel_solve', entrypoint=__file__,
                        inputs=[ROOT / DESIGN / 'layout/build.py', core / 'raw' / 'cmat_fF.csv'],
                        parameters=dict(core=a.core, coupons=a.coupons, overrides=json.loads(a.overrides),
                                        targets=out['targets'], anchor=anchor, anchor_run=a.anchor_run, kp_correction=kp_corr))
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    (folder / 'summaries' / 'summary.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    finish_run(folder, status='completed', summary=out)
    print('RUN_FOLDER=' + str(folder))
    print(json.dumps(out['channels'], indent=2))


if __name__ == '__main__':
    main()
