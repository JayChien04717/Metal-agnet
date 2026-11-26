# -*- coding: utf-8 -*-
r"""Step D (2026-11-26): predict f_p, f_r and J from Q3D only and compare with the HFSS complex fit.

Inputs
  --line    results/resonator/<q3d_cpw_line run>       C' and eps_eff of the 10/6 CPW (step B) -> v, Z0
  --insitu  results/resonator/<q3d_channel_insitu run>  in-situ Maxwell matrix per channel (step C)
  --fit     results/resonator/<fit_s21 run>              HFSS f_p, f_r, J (no wire bonds on the channel lines)
Model (purcell_readout_model with the Q3D-calibrated line instead of the analytic one):
  * nets: filter (PurcellReadout filter + C_p plate), resonator (PurcellReadout resonator + extension + claw), cross
    (floating), and the ground plane (which includes the C_p frame: the feed is a low-impedance node);
  * C_j = -C(filter, resonator) in situ;
  * line ground C: filter / resonator row sums to ground; the cross is Kron-reduced (floating island);
  * end loads: C_excess = C_ground - C' (L + s/2) - C_far, where C_far is the fringe of the Q3D-only etched gap at the
    opened short end (open_end_fF of the line run). The stub load at the tap is C_j (+ head_excess_fF, default 0);
    the rest of C_excess sits at the open end (C_p plate / claw). Misplacing part of it between the open end and
    the tap changes f by only ~(1 - cos^2(k tap)) of that part;
  * f from tapped_lambda4_f, J from j_hammer, both with the calibrated line.
The HFSS - model offsets per line and the J ratio are reported; the goal is small, channel-independent offsets.

    python -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/resonator/q3d_model_check.py --line <run> --insitu <run> --fit <run>
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
import pandas as pd

import build as D11
import purcell_readout_model as PRM
from run_record import create_run, finish_run

DESIGN = 'Design11_Transmon2Q_Manhattan'
C0 = 299792458.0


def calibrated_line(line_run):
    s = json.loads((line_run / 'summaries' / 'summary.json').read_text(encoding='utf-8'))
    v = C0 / np.sqrt(s['eps_eff'])
    cl = s['c_line_fF_per_mm']
    return dict(width_um=10, gap_um=6, v=v, z0=1e15 / (v * cl * 1e3), c_line=cl,
                source=f'Q3D {line_run.name}'), s


def net(names, key, q=None):
    hit = [n for n in names if key in n.lower() and (q is None or n.endswith(f'_{q}'))]
    if len(hit) != 1:
        raise RuntimeError(f'net {key!r} not unique in {names}')
    return hit[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--line', required=True)
    ap.add_argument('--insitu', required=True)
    ap.add_argument('--fit', required=True)
    ap.add_argument('--far-end-fF', type=float, default=None, help='fringe of the opened short end (default: line run open_end_fF, floored at 0)')
    ap.add_argument('--head-excess-fF', type=float, default=0.0)
    a = ap.parse_args()
    line_run, ins_run, fit_run = ROOT / a.line, ROOT / a.insitu, ROOT / a.fit
    cpw, line_sum = calibrated_line(line_run)
    c_far = a.far_end_fF if a.far_end_fF is not None else max(0.0, line_sum['open_end_fF'])
    ins = json.loads((ins_run / 'summaries' / 'summary.json').read_text(encoding='utf-8'))
    fit = json.loads((fit_run / 'summaries' / 'summary.json').read_text(encoding='utf-8'))['channels']
    stub = D11.HAMMER_LINE_SPACING_UM / 2000.0
    out = dict(line=cpw, c_far_fF=c_far, head_excess_fF=a.head_excess_fF, channels={})
    for q, ch in D11.PURCELL_CHANNELS.items():
        if q not in ins:
            continue
        names = ins[q]['nets']
        M = pd.DataFrame(ins[q]['maxwell_fF'], index=names, columns=names)
        nf, nr, nx = net(names, 'filter', q), net(names, 'readout_connector', q), net(names, 'cross', q)
        cj = float(-M.loc[nf, nr])
        cqr = float(-M.loc[nr, nx])
        c_xf = float(-M.loc[nx, nf])
        # every other net (ground plane incl. the C_p frame, Q3D outer ground) counts as ground
        c_xg = float(M.loc[nx, nx]) - cqr - c_xf                               # cross to ground
        cf_g = float(M.loc[nf, nf]) - cj - c_xf                                # filter to ground (incl. C_p mutual)
        cr_g = float(M.loc[nr, nr]) - cj - cqr + cqr * c_xg / (cqr + c_xg)     # resonator to ground, cross floating
        tap_f = 1.20
        tap_r = float(D11_tap_r(q))       # resonator_tap + open extension, as built
        Lf, Lr = ch['L_filter'], ch['L_res']
        ex_f = cf_g - cpw['c_line'] * (Lf + stub) - c_far
        ex_r = cr_g - cpw['c_line'] * (Lr + stub) - c_far
        c_stub = cj + a.head_excess_fF
        end_f, end_r = ex_f - a.head_excess_fF, ex_r - a.head_excess_fF
        fp = PRM.tapped_lambda4_f(Lf, tap_f, end_f, stub, c_stub, cpw=cpw)
        fr = PRM.tapped_lambda4_f(Lr, tap_r, end_r, stub, c_stub, cpw=cpw)
        J = PRM.j_hammer(0.5 * (fp + fr), cj, Lf, tap_f, end_f, Lr, tap_r, end_r, cpw=cpw)
        h = fit.get(q, {})
        out['channels'][q] = dict(
            C_j_insitu_fF=cj, C_qr_fF=cqr, C_cross_ground_fF=c_xg, C_filter_ground_fF=cf_g, C_res_ground_kron_fF=cr_g,
            end_load_filter_fF=end_f, end_load_res_fF=end_r, tap_r_mm=tap_r,
            model=dict(f_p_GHz=fp, f_r_GHz=fr, J_MHz=J),
            hfss=dict(f_p_GHz=h.get('f_p_GHz'), f_r_GHz=h.get('f_r_GHz'), J_MHz=h.get('J_MHz'), kappa_p_MHz=h.get('kappa_p_MHz')),
            offset_MHz=dict(filter=(h['f_p_GHz'] - fp) * 1e3 if h else None, resonator=(h['f_r_GHz'] - fr) * 1e3 if h else None),
            J_ratio_hfss_over_model=(h['J_MHz'] / J) if h else None)
    folder = create_run(DESIGN, 'resonator', 'q3d_model_check', entrypoint=__file__,
                        inputs=[ROOT / DESIGN / 'layout' / 'build.py'],
                        parameters=dict(line_run=a.line, insitu_run=a.insitu, fit_run=a.fit, c_far_fF=c_far,
                                        head_excess_fF=a.head_excess_fF))
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    (folder / 'summaries' / 'summary.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    finish_run(folder, status='completed', summary=out)
    print('RUN_FOLDER=' + str(folder))
    print(json.dumps(out['channels'], indent=1))


def D11_tap_r(q):
    """Resonator tap distance used by build.py: PurcellReadout resonator_tap + resonator_open_extension."""
    d = D11.build_design(variant='sim', channels=True, check_keepout_um=None)
    o = d.components[f'prr_{q}'].options
    return d.parse_value(o.resonator_tap) + d.parse_value(o.resonator_open_extension)


if __name__ == '__main__':
    main()
