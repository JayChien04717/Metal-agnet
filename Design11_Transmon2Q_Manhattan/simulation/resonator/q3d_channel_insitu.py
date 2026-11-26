# -*- coding: utf-8 -*-
r"""In-situ Q3D of each PurcellReadout channel in the actual Design11 layout (step C, 2026-11-26).

Instead of the isolated parallel-line hammer coupon (not used, user 2026-11-26), each channel is rendered from
layout/build.py as built (200 um hammer spacing, current PURCELL_CHANNELS):
    prr_<Q>      PurcellReadout with short_ends='open' (etched gap at each meander end so filter and resonator stay
                 separate nets; note/purcell_filter_method.md "render shorted ends open")
    cp_<Q>       EnclosingCap C_p (its feed-side frame ends in ground: the feed is a low-impedance node)
    res_open_<Q> resonator open extension to the claw
    <Q>          Xmon with its readout claw (cross = floating island, flux line as ground, no junction metal)
No wire bonds (the fabricated chip has none on these lines); hammer straps are not rendered in Q3D (they are
ground-to-ground bridges, method note section 5).
Output: full Maxwell matrix per channel (raw/<case>/matrix.csv, caps.json); step D
(simulation/resonator/q3d_model_check.py) turns it into in-situ C_j, line ground C, C_qr and C_p loads.

    python -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/resonator/q3d_channel_insitu.py [--passes 4] [--mesh-um 10]
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator'),
                str(ROOT / 'Design11_Transmon2Q_Manhattan' / 'layout'),
                str(ROOT / 'Design11_Transmon2Q_Manhattan' / 'simulation' / 'qubit')]
import project_paths
project_paths.configure(__file__)

import build as D11
import hfss_dirs
import purcell_resonator_pair_sim as S
import save_aedt
from run_record import create_run, finish_run

DESIGN = 'Design11_Transmon2Q_Manhattan'
CORE_OV = dict(xmon=dict(squid=dict(contact_overlap='0um', shift_toward_arm='0um')),
               coupler=dict(squid_contact_overlap='0um', squid_shift_toward_stem='0um'))


def clip_ground_pads(d, names):
    """Same as q3d_core._clip_ground_pads_for_q3d: squid / jj pad grounds flush with the pocket for Q3D."""
    from shapely.ops import unary_union
    tbl = d.qgeometry.tables['poly']
    for cname in names:
        cid = d.components[cname].id
        sub_u = unary_union(list(tbl[(tbl['component'] == cid) & (tbl['subtract'])]['geometry']))
        mask = (tbl['component'] == cid) & (~tbl['subtract']) & (tbl['name'].str.contains('squid_ground|jj_pad_ground'))
        for idx in tbl.index[mask]:
            tbl.at[idx, 'geometry'] = tbl.at[idx, 'geometry'].intersection(sub_u)


def channel_design(qs):
    qs = [qs] if isinstance(qs, str) else list(qs)
    d = D11.build_design(variant='sim', overrides=CORE_OV, channels=True, flux_as_ground=True, check_keepout_um=None)
    for q in qs:
        prr = d.components[f'prr_{q}']
        prr.options.short_ends = 'open'
        prr.rebuild()
    clip_ground_pads(d, qs)
    return d


# all components of the HFSS driven fixture (hfss_readout_s21.build_driven_fixture), no ports
FULL_FIXTURE = ['feed_cin', 'cin', 'feed_in', 'feed_mid', 'feed_out', 'Q1', 'Q2'] + [
    f'{p}_{q}' for q in ('Q1', 'Q2') for p in ('n0', 'stub_p', 'cp', 'prr', 'res_open')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--passes', type=int, default=4)
    ap.add_argument('--pct', type=float, default=0.1)
    ap.add_argument('--mesh-um', type=float, default=10.0, help='length mesh on all metal (user 2026-11-26: 4 passes, 10 um)')
    ap.add_argument('--channels', default='Q1,Q2')
    ap.add_argument('--tag', default='r16')
    ap.add_argument('--full', action='store_true',
                    help='one case with the whole HFSS fixture (feed, C_in, both channels, both qubits), nets separate')
    a = ap.parse_args()
    chans = [c.strip() for c in a.channels.split(',') if c.strip()]
    folder = create_run(DESIGN, 'resonator', f'q3d_channel_insitu_{a.tag}', entrypoint=__file__,
                        inputs=[ROOT / DESIGN / 'layout' / 'build.py'],
                        parameters=dict(channels=chans, max_passes=a.passes, percent_error=a.pct, mesh_um=a.mesh_um,
                                        purcell_channels=D11.PURCELL_CHANNELS,
                                        hammer_line_spacing_um=D11.HAMMER_LINE_SPACING_UM,
                                        short_ends='open', wire_bonds='none', full_fixture=a.full,
                                        components=FULL_FIXTURE if a.full else 'per channel'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        hfss_dirs.use_d_drive()
        save_aedt.use_project(ROOT / DESIGN / 'simulation' / 'aedt' / 'D11.aedt')
        out = {}
        if a.full:
            d = channel_design(chans)
            r = S.q3d_case('full', d, FULL_FIXTURE, [], (), folder, folder.name[-4:],
                           max_passes=a.passes, pct_error=a.pct, freq_ghz=6.25, mesh_um=a.mesh_um)
            for q in chans:                        # same matrix for each channel; step D picks nets by suffix
                out[q] = dict(nets=r['maxwell_order'], maxwell_fF=r['maxwell_fF'], minutes=r['minutes'],
                              converged=r['solver_converged'], case='full')
        for q in ([] if a.full else chans):
            d = channel_design(q)
            comps = [q, f'prr_{q}', f'cp_{q}', f'res_open_{q}']
            opens = []
            r = S.q3d_case(f'ch_{q}', d, comps, opens, (), folder, folder.name[-4:],
                           max_passes=a.passes, pct_error=a.pct, freq_ghz=6.25,
                               mesh_um=a.mesh_um)
            out[q] = dict(nets=r['maxwell_order'], maxwell_fF=r['maxwell_fF'], minutes=r['minutes'],
                          converged=r['solver_converged'])
        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
        finish_run(folder, status='completed', summary={q: dict(nets=v['nets']) for q, v in out.items()})
        print(json.dumps({q: v['nets'] for q, v in out.items()}, indent=2), flush=True)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
