# -*- coding: utf-8 -*-
r"""Q3D coupons for Design11 PurcellReadout channels at the layout's actual hammer / C_p geometry.

    python -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/resonator/q3d_readout_coupons.py \
        --gaps 30,45,60,86 --plate-widths 140,180,220 --tag optB_sweep

Cases (8 passes / 0.3 %, local solves, skill/purcell_resonator_design.md section 3):
  cj_g<gap>  PurcellReadout draw_mode='coupon' with build.py hammer options (head 20 x 40 um, 10/6 CPW) and head_gap=<gap>
  (no parallel-line reference coupon: C_j = mutual(cj) directly, head excess against the analytic straight CPW;
   see optB_channel_solve.coupon_tables)
  cp_w<pw>   EnclosingCap (build.ENC_OPT with plate_width=<pw>) alone -> C_p mutual, plate / frame ground
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / 'code' / 'general'),
    str(ROOT / 'code' / 'resonator'),
    str(ROOT / 'Design11_Transmon2Q_Manhattan' / 'layout'),
]
import project_paths
project_paths.configure(__file__)

from qiskit_metal import designs
import build as D11
import hfss_dirs
import purcell_readout_model as PRM
import purcell_resonator_pair_sim as S
import save_aedt
from enclosing_cap import EnclosingCap
from purcell_readout import PurcellReadout
from round_corners import round_component_polys
from run_record import create_run, finish_run

DESIGN = 'Design11_Transmon2Q_Manhattan'
HAMMER = dict(line_spacing='400um', tee_half='100um', lead='100um', head_width='20um', head_height='40um',
              hammer_ground_gap='16um', corner_radius='5um')
COUPON_LEAD_UM = 400.0


def _small(sx, sy):
    d = designs.DesignPlanar()
    d.overwrite_enabled = True
    d.chips.main.material = 'silicon_eps11p45'
    d.chips.main.size.update(size_x=f'{sx}mm', size_y=f'{sy}mm', size_z='-500um', center_x='0mm', center_y='0mm')
    return d


def tag(v):
    """AEDT object / design names may not contain '.', so 38.5 -> '38p5'."""
    return f'{v:g}'.replace('.', 'p')


def untag(s):
    return float(s.replace('p', '.'))


def cases(gaps, widths, hammer):
    out = {}
    for g in gaps:
        d = _small(2.0, 2.0)
        PurcellReadout(d, f'cj_g{tag(g)}', options=dict(draw_mode='coupon', coupon_lead=f'{COUPON_LEAD_UM:g}um',
                                                       head_gap=f'{g:g}um', **hammer, **D11.CPW_RO))
        out[f'cj_g{tag(g)}'] = (d, [f'cj_g{tag(g)}'], [], ('filter', 'resonator'))
    for w in widths:
        d = _small(1.5, 1.2)
        EnclosingCap(d, 'cp', options=dict(D11.ENC_OPT, plate_width=f'{w:g}um', **D11.CPW_RO))
        round_component_polys(d.components['cp'], 5e-3, skip_names=())
        out[f'cp_w{tag(w)}'] = (d, ['cp'], [('cp', 'a'), ('cp', 'b')], ('frame|lead_a', 'plate|lead_b'))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gaps', default='30,45,60,86')
    ap.add_argument('--plate-widths', default='140,180,220')
    ap.add_argument('--head', default='', help='override head as W,H in um (default build.py 20,40)')
    ap.add_argument('--tag', default='coupons')
    a = ap.parse_args()
    gaps = [float(v) for v in a.gaps.split(',') if v]
    widths = [float(v) for v in a.plate_widths.split(',') if v]
    hammer = dict(HAMMER)
    if a.head:
        w, h = a.head.split(',')
        hammer.update(head_width=f'{w}um', head_height=f'{h}um')

    folder = create_run(DESIGN, 'resonator', f'q3d_coupons_{a.tag}', entrypoint=__file__,
                        inputs=[ROOT / DESIGN / 'layout/build.py'],
                        parameters=dict(coupon_lead_um=COUPON_LEAD_UM, gaps_um=gaps, plate_widths_um=widths, hammer=hammer, enc_opt=D11.ENC_OPT,
                                        cpw=D11.CPW_RO, max_passes=8, percent_error=0.3, freq_ghz=6.25))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        hfss_dirs.use_d_drive()
        save_aedt.use_project(ROOT / DESIGN / 'simulation' / 'aedt' / 'D11.aedt')
        res = {}
        for name, (d, comps, opens, keys) in cases(gaps, widths, hammer).items():
            res[name] = S.q3d_case(name, d, comps, opens, keys, folder, folder.name[-4:],
                                   max_passes=8, pct_error=0.3, freq_ghz=6.25)
        summary = dict(cj={}, cp={})
        for g in gaps:
            r = res[f'cj_g{tag(g)}']
            summary['cj'][f'{g:g}'] = dict(cj=r['mutual_fF'], ground_fF=r['ground_fF'])
        for w in widths:
            r = res[f'cp_w{tag(w)}']
            summary['cp'][f'{w:g}'] = dict(mutual_fF=r['mutual_fF'], ground_fF=r['ground_fF'],
                                           kappa_hanger_6p25_MHz=PRM.kappa_hanger(6.25, r['mutual_fF']))
        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(json.dumps(summary, indent=2, default=float), encoding='utf-8')
        finish_run(folder, status='completed', summary=summary)
        print(json.dumps(summary, indent=2, default=float), flush=True)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
