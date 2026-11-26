# -*- coding: utf-8 -*-
r"""HFSS Eigenmode verification of Purcell filter + readout resonator frequencies for Design11_Transmon2Q_Manhattan.

Renders the loaded readout fixture (feedline with C_in, both TapJunctions, both C_p, both PurcellReadout
modules with PEC straps, open extensions, and Q1 + Q2 with Manhattan pads & lumped Lj; flux lines as ground,
TC omitted to keep mesh focused on readout channels).

Solver defaults (skill/em_simulation.md):
  - HFSS eigenmode, max_passes = 6, basis_order = 1
  - Explicit 10 um length mesh on all metal objects (mesh_all)
  - Auto wire bonds every 350 um on etched paths + PEC straps at PurcellReadout.strap_points()
  - 50 Ohm ports at feed_cin start (P8 side) and feed_out start (P7 side)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / 'code' / 'general'),
    str(ROOT / 'code' / 'resonator'),
    str(ROOT / 'Design11_Transmon2Q_Manhattan' / 'layout'),
]
import project_paths
project_paths.configure(__file__)

import numpy as np
import build as D11
import hfss_dirs
import purcell_resonator_pair_sim as S
import save_aedt
from run_record import create_run, finish_run

CHANNELS = ('Q2', 'Q1')  # sorted by frequency: Q2 (6.15 GHz), Q1 (6.35 GHz)


def build_fixture():
    d = D11.build_design(variant='sim', channels=True, flux_as_ground=True, check_keepout_um=None)
    t = d.qgeometry.tables['path']
    t['hfss_wire_bonds'] = t['subtract'].astype(bool)
    comps = ['feed_cin', 'cin', 'feed_in', 'feed_mid', 'feed_out', 'Q1', 'Q2']
    for q in ('Q1', 'Q2'):
        comps += [f'{p}_{q}' for p in ('n0', 'stub_p', 'cp', 'prr', 'res_open')]
    ports = [('feed_cin', 'start', 50), ('feed_out', 'start', 50)]
    opens = []
    return d, comps, ports, opens


def analyse_modes(raw: dict) -> dict:
    """Pair the 4 complex readout/Purcell modes (two per channel) and invert to bare (f_p, f_r, J, kappa_p)."""
    freqs = np.asarray(raw['frequencies_GHz'], float)
    imags = np.asarray(raw['imag_GHz'], float)
    modes = [dict(idx=i, f_GHz=float(f), im_GHz=float(im), kappa_MHz=float(2.0 * abs(im) * 1e3))
             for i, (f, im) in enumerate(zip(freqs, imags))]
    # Filter modes in the readout band (5.6 - 6.9 GHz)
    in_band = [m for m in modes if 5.6 <= m['f_GHz'] <= 6.9]
    channels = {}
    if len(in_band) >= 4:
        # Sort by frequency; lowest pair -> Q2 (6.15 GHz), next pair -> Q1 (6.35 GHz)
        in_band.sort(key=lambda m: m['f_GHz'])
        for ch, pair in (('Q2', in_band[0:2]), ('Q1', in_band[2:4])):
            z1 = pair[0]['f_GHz'] + 1j * abs(pair[0]['im_GHz'])
            z2 = pair[1]['f_GHz'] + 1j * abs(pair[1]['im_GHz'])
            inv = S.invert(z1, z2)
            # Effective resonator linewidth kappa_r ≈ 4 J^2 / kappa_p (on resonance)
            kp = inv['kappa_p_MHz']
            J = abs(inv['J_MHz'])
            kr_est = (4.0 * J * J / kp) if kp > 0 else float('nan')
            inv['kappa_r_est_MHz'] = float(kr_est)
            inv['modes'] = pair
            inv['target_f_GHz'] = D11.PURCELL_CHANNELS[ch]['f_GHz']
            inv['df_p_MHz'] = (inv['f_p_GHz'] - inv['target_f_GHz']) * 1e3
            inv['df_r_MHz'] = (inv['f_r_GHz'] - inv['target_f_GHz']) * 1e3
            channels[ch] = inv
    return dict(all_modes=modes, in_band_count=len(in_band), channels=channels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--passes', type=int, default=6)
    ap.add_argument('--min-freq', type=float, default=5.4)
    ap.add_argument('--n-modes', type=int, default=6)
    ap.add_argument('--tag', default='ro_eig')
    ap.add_argument('--fine-mesh-um', type=float, default=0.0,
                    help='local length mesh (um) on the hammer tees, readout claws (qubit fork) and qubit crosses; 0 = off')
    a = ap.parse_args()

    d, comps, ports, opens = build_fixture()
    folder = create_run(
        'Design11_Transmon2Q_Manhattan', 'resonator', f'hfss_eigen_{a.tag}',
        entrypoint=__file__,
        inputs=[ROOT / 'Design11_Transmon2Q_Manhattan/design.json',
                ROOT / 'Design11_Transmon2Q_Manhattan/layout/build.py'],
        parameters=dict(
            components=comps, ports=ports, open_terminations=opens,
            max_passes=a.passes, min_freq_ghz=a.min_freq, n_modes=a.n_modes, fine_mesh_um=a.fine_mesh_um,
            basis_order=1, mesh_length=S.MESH, wire_bonds=S.WB,
            purcell_channels=D11.PURCELL_CHANNELS,
        ),
    )
    print('RUN_FOLDER=' + str(folder), flush=True)

    hfss_dirs.use_d_drive()
    project = ROOT / 'Design11_Transmon2Q_Manhattan' / 'simulation' / 'aedt' / 'D11.aedt'
    project.parent.mkdir(parents=True, exist_ok=True)
    save_aedt.use_project(project)

    fm = None
    if a.fine_mesh_um > 0:
        fm = {f'{a.fine_mesh_um:g}um': [f'{n}_{c}' for q in ('Q1', 'Q2')
                                         for n, c in (('filter_tee', f'prr_{q}'), ('resonator_tee', f'prr_{q}'),
                                                      ('readout_connector_arm', q), ('cross', q))]}
    raw = S.run_eigen(folder, d, comps, ports, opens, minf=a.min_freq, n_modes=a.n_modes, passes=a.passes, tag='E',
                      fine_mesh=fm)
    ana = analyse_modes(raw)
    summary = dict(raw=raw, analysis=ana)
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    (folder / 'summaries' / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    finish_run(folder, status='completed', summary=summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
