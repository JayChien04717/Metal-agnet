# -*- coding: utf-8 -*-
r"""Q3D calibration of the 10/6 um CPW line constants used by purcell_readout_model (step B, 2026-11-26).

One straight CPW (no second line, no hammer) with etched open ends, on the Design11 substrate
(silicon_eps11p45, 500 um) and, for the vacuum reference, the same geometry with a vacuum substrate.
Two lengths cancel the open-end fringe:
    C'   = (C(L2) - C(L1)) / (L2 - L1)          per-unit-length C to ground
    C0'  = same with a vacuum substrate          eps_eff = C' / C0'
    v    = c / sqrt(eps_eff),  Z0 = 1 / (v C')   (quasi-TEM)
    C_open = (C(L1) - C' L1) / 2                 fringe C of one etched open end
Compared with the analytic constants in purcell_readout_model.CPW['10/6'] (conformal mapping, eps_eff 6.225).

    python -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/resonator/q3d_cpw_line.py [--lengths 1,3] [--passes 4] [--mesh-um 10]
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator'),
                str(ROOT / 'Design11_Transmon2Q_Manhattan' / 'layout')]
import project_paths
project_paths.configure(__file__)

import numpy as np
from qiskit_metal import Dict, designs, draw
from qiskit_metal.qlibrary.core import QComponent

import build as D11
import hfss_dirs
import purcell_readout_model as PRM
import purcell_resonator_pair_sim as S
import save_aedt
from run_record import create_run, finish_run

DESIGN = 'Design11_Transmon2Q_Manhattan'
C0 = 299792458.0


class OneLine(QComponent):
    """One straight CPW along y, length `length`, etched open ends (end gap = trace_gap); net 'line'."""
    default_options = Dict(length='1mm', trace_width='10um', trace_gap='6um')

    def make(self):
        p = self.p
        h = p.length / 2
        ln = draw.LineString([(0, -h), (0, h)])
        cut = draw.LineString([(0, -h - p.trace_gap), (0, h + p.trace_gap)])
        self.add_qgeometry('poly', {'line': ln.buffer(p.trace_width / 2, cap_style=2)}, layer=1)
        self.add_qgeometry('poly', {'line_cut': cut.buffer(p.trace_width / 2 + p.trace_gap, cap_style=2)},
                           layer=1, subtract=True)


def small(sx, sy, material):
    d = designs.DesignPlanar()
    d.overwrite_enabled = True
    d.chips.main.material = material
    d.chips.main.size.update(size_x=f'{sx}mm', size_y=f'{sy}mm', size_z='-500um', center_x='0mm', center_y='0mm')
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lengths', default='1,3', help='two line lengths in mm')
    ap.add_argument('--passes', type=int, default=4)
    ap.add_argument('--pct', type=float, default=0.1)
    ap.add_argument('--mesh-um', type=float, default=10.0, help='length mesh on all metal (user 2026-11-26: 4 passes, 10 um)')
    ap.add_argument('--tag', default='cpw106')
    a = ap.parse_args()
    L1, L2 = sorted(float(v) for v in a.lengths.split(','))
    mats = dict(si='silicon_eps11p45', vac='vacuum')
    folder = create_run(DESIGN, 'resonator', f'q3d_cpw_line_{a.tag}', entrypoint=__file__, inputs=[],
                        parameters=dict(lengths_mm=[L1, L2], max_passes=a.passes, percent_error=a.pct, mesh_um=a.mesh_um,
                                        cpw=D11.CPW_RO, materials=mats, chip_width_mm=1.5, chip_margin_mm=1.0))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        hfss_dirs.use_d_drive()
        save_aedt.use_project(ROOT / DESIGN / 'simulation' / 'aedt' / 'D11.aedt')
        caps = {}
        for mk, mat in mats.items():
            for L in (L1, L2):
                name = f'line_{mk}_L{L:g}'.replace('.', 'p')
                d = small(1.5, L + 1.0, mat)
                OneLine(d, 'line', options=dict(length=f'{L:g}mm', **D11.CPW_RO))
                r = S.q3d_case(name, d, ['line'], [], ('line',), folder, folder.name[-4:],
                               max_passes=a.passes, pct_error=a.pct, freq_ghz=6.25,
                               mesh_um=a.mesh_um)
                caps[(mk, L)] = r['ground_fF']['line']
        cp = (caps[('si', L2)] - caps[('si', L1)]) / (L2 - L1)
        c0 = (caps[('vac', L2)] - caps[('vac', L1)]) / (L2 - L1)
        eps = cp / c0
        v = C0 / np.sqrt(eps)
        z0 = 1e15 / (v * cp * 1e3)
        model = PRM.line('10/6')
        out = dict(c_line_fF_per_mm=cp, c_line_vacuum_fF_per_mm=c0, eps_eff=eps, v_m_s=v, z0_ohm=z0,
                   open_end_fF=(caps[('si', L1)] - cp * L1) / 2,
                   raw_ground_fF={f'{k[0]}_L{k[1]:g}': c for k, c in caps.items()},
                   analytic=dict(c_line_fF_per_mm=model['c_line'], eps_eff=(C0 / model['v']) ** 2, v_m_s=model['v'],
                                 z0_ohm=model['z0'], open_end_fF=PRM.open_end_c('10/6')),
                   ratio_q3d_over_analytic=dict(c_line=cp / model['c_line'], v=v / model['v'], z0=z0 / model['z0']))
        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
        finish_run(folder, status='completed', summary=out)
        print(json.dumps(out, indent=2), flush=True)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
