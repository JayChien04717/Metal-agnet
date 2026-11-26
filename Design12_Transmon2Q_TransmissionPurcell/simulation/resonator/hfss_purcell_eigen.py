# -*- coding: utf-8 -*-
r"""HFSS eigenmode frequency check of the Design12 shared lambda/4 Purcell filter + RR_Q1 / RR_Q2 (readout-only fixture).

Same fixture as hfss_purcell_s21.py (build_design(readout_only=True): no qubits / TC / control lines, resonator open
ends extended by the claw-equivalent length; 50 Ohm lumped port sheets at P7 and P8 so the filter mode is loaded).
Solver (skill/em_simulation.md): eigenmode, first-order basis, 10 um length mesh on every rendered metal object,
--passes adaptive passes (default 6), wire bonds every 350 um on feed_in / feed_out only.

Mode identification: the lowest-Q mode near f_PF is the filter; the two highest-Q modes are the resonators (lower
frequency -> lower-target resonator). With |delta| >> J the hybridisation shift J^2 / delta is a few MHz, so the
resonator modes are close to the bare f_r; the shift estimate from the r1 driven J is reported alongside.

Low-memory option (user 2026-11-26): --crop x0,x1,y0,y1 (mm) shrinks the rendered chip (substrate, ground and air
box) to a window around the readout fixture instead of the full 10 x 10 mm die; --modes 3 solves only the filter and
the two resonators. HFSS 2021 R2 eigenmode has no iterative solver, so memory is cut by volume and mode count.

Usage (repository root, qmetal):
    python -X utf8 -B Design12_Transmon2Q_TransmissionPurcell/simulation/resonator/hfss_purcell_eigen.py --passes 6 --tag r2
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code' / 'general'))
import project_paths
project_paths.configure(__file__)
DESIGN = project_paths.design_dir(__file__)
sys.path.insert(0, str(DESIGN / 'layout'))
sys.path.insert(0, str(Path(__file__).resolve().parent))


def classify(raw, targets, f_pf):
    modes = [dict(f_GHz=f, imag_GHz=i, Q=(f / (2 * i)) if i > 0 else float('inf'))
             for f, i in zip(raw['frequencies_GHz'], raw['imag_GHz'])]
    near = [m for m in modes if abs(m['f_GHz'] - f_pf) < 0.8]
    filt = min(near, key=lambda m: m['Q']) if near else None
    rest = sorted([m for m in near if m is not filt], key=lambda m: -m['Q'])[:len(targets)]
    rest = sorted(rest, key=lambda m: m['f_GHz'])
    names = sorted(targets, key=targets.get)
    out = dict(all_modes=modes)
    if filt:
        out['PF'] = dict(f_GHz=filt['f_GHz'], kappa_MHz=2 * filt['imag_GHz'] * 1e3, Q=filt['Q'], target_GHz=f_pf,
                         df_MHz=(filt['f_GHz'] - f_pf) * 1e3)
    for q, m in zip(names, rest):
        out[q] = dict(f_GHz=m['f_GHz'], kappa_MHz=2 * m['imag_GHz'] * 1e3, Q=m['Q'], target_GHz=targets[q],
                      df_MHz=(m['f_GHz'] - targets[q]) * 1e3)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--passes', type=int, default=6)
    ap.add_argument('--min-freq', type=float, default=5.5)
    ap.add_argument('--modes', type=int, default=5)
    ap.add_argument('--tag', default='r2')
    ap.add_argument('--crop', default='', help='x0,x1,y0,y1 in mm: render only this chip window (low memory)')
    a = ap.parse_args()
    from run_record import create_run, finish_run
    import hfss_dirs
    import save_aedt
    import purcell_resonator_pair_sim as S
    from hfss_purcell_s21 import build_fixture
    import build as B
    targets = {q: B.RR_CHANNELS[q]['f_GHz'] for q in ('Q1', 'Q2')}
    folder = create_run(
        DESIGN.name, 'resonator', f'purcell_eigen_{a.tag}', entrypoint=__file__,
        inputs=[DESIGN / 'design.json', DESIGN / 'layout' / 'build.py'],
        parameters=dict(fixture='readout_only (PF + RR_Q1 + RR_Q2, no qubits / TC / control lines)',
                        ports='50 Ohm lumped port sheets at P7 (feed_in start) and P8 (feed_out end)',
                        max_passes=a.passes, min_freq_GHz=a.min_freq, n_modes=a.modes, basis_order=1,
                        max_delta_f_pct=0.1, crop_mm=a.crop or 'full 10 x 10 mm chip', mesh_length=S.MESH, mesh_targets='every rendered metal object',
                        wire_bonds=dict(S.WB, on='feed_in, feed_out only'), geometry_revision=B.REVISION,
                        pf=B.PF, rr=B.RR_CHANNELS, hfss_anchor=B.HFSS_ANCHOR))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        hfss_dirs.use_d_drive()
        project = DESIGN / 'simulation' / 'aedt' / 'D12.aedt'
        save_aedt.use_project(project)
        d, _, comps, ports, opens = build_fixture()
        if a.crop:
            x0, x1, y0, y1 = map(float, a.crop.split(','))
            d.chips.main.size.update(size_x=f'{x1 - x0}mm', size_y=f'{y1 - y0}mm',
                                     center_x=f'{(x0 + x1) / 2}mm', center_y=f'{(y0 + y1) / 2}mm')
        raw = S.run_eigen(folder, d, comps, ports, opens, minf=a.min_freq, n_modes=a.modes, passes=a.passes, tag='E12')
        out = classify(raw, targets, B.PF['f_GHz'])
        summary = dict(hfss_eigen=raw, modes=out)
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
