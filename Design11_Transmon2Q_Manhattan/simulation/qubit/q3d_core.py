# -*- coding: utf-8 -*-
r"""Core Q3D capacitance extraction for Design11_Transmon2Q_Manhattan (8 passes / 0.5%).

Includes Q1, Q2, TC (with Manhattan junction landing metal on layer 1), readout claws, and XY lines;
flux lines are drawn as ground (flux_as_ground=True). Extracts C_sigma, Ec, g_n(Q1-TC, Q2-TC, Q1-Q2),
and readout coupling g_qr.
"""
from argparse import ArgumentParser
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [
    str(ROOT / 'code' / 'general'),
    str(ROOT / 'Design11_Transmon2Q_Manhattan' / 'layout'),
]
import project_paths
project_paths.configure(__file__)

import build as D11
from run_record import create_run, finish_run

COMPS = ['Q1', 'Q2', 'TC', 'xy_Q1', 'xy_Q2', 'xy_end_Q1', 'xy_end_Q2']
OPENS = [('Q1', 'readout'), ('Q2', 'readout'), ('xy_Q1', 'start'), ('xy_Q2', 'start')]
SETUP = dict(min_passes=2, max_passes=8, percent_error=0.5, freq_ghz=5.0,
             solution_order='Medium', auto_increase_solution_order=False)
E_CHARGE = 1.602176634e-19
H_PLANCK = 6.62607015e-34


def analyse_cmat(cmat_fF):
    """Extract Ec [GHz], C_sigma [fF], g_n [MHz], and g_qr [MHz] from Q3D Maxwell capacitance matrix."""
    cols = list(cmat_fF.columns)

    def find_net(substrs):
        for c in cols:
            cl = c.lower()
            if all(s.lower() in cl for s in substrs):
                return c
        raise KeyError(f'Net matching {substrs} not found in {cols}')

    n_q1 = find_net(['cross', 'q1'])
    n_q2 = find_net(['cross', 'q2'])
    n_tc = find_net(['tc', 'pad']) if any('tc' in c.lower() and 'pad' in c.lower() for c in cols) else find_net(['tc'])
    n_ro1 = find_net(['readout', 'q1'])
    n_ro2 = find_net(['readout', 'q2'])

    # 3x3 island Maxwell capacitance submatrix (all other nets grounded)
    islands = [n_q1, n_q2, n_tc]
    C_3 = cmat_fF.loc[islands, islands].to_numpy(dtype=float) * 1e-15  # Farads
    Cinv = np.linalg.inv(C_3)

    ec_ghz = {k: float(E_CHARGE**2 * Cinv[i, i] / (2.0 * H_PLANCK) / 1e9) for i, k in enumerate(('Q1', 'Q2', 'TC'))}
    c_sig = {k: float(1.0 / Cinv[i, i] * 1e15) for i, k in enumerate(('Q1', 'Q2', 'TC'))}
    g_n_mhz = {
        'Q1-TC': float(4.0 * E_CHARGE**2 * Cinv[0, 2] / H_PLANCK / 1e6),
        'Q2-TC': float(4.0 * E_CHARGE**2 * Cinv[1, 2] / H_PLANCK / 1e6),
        'Q1-Q2': float(4.0 * E_CHARGE**2 * Cinv[0, 1] / H_PLANCK / 1e6),
    }
    # Readout coupling g_qr (lambda/4 50-ohm resonator estimate: C_r = pi / (4 * omega_r * Z0))
    f_q = {'Q1': 4.739e9, 'Q2': 4.342e9}
    f_r = {'Q1': 6.35e9, 'Q2': 6.15e9}
    ro_nets = {'Q1': (n_q1, n_ro1, 0), 'Q2': (n_q2, n_ro2, 1)}
    readout = {}
    for k, (nq, nro, idx) in ro_nets.items():
        c_qr_fF = float(-cmat_fF.loc[nq, nro])
        c_r_fF = float(np.pi / (4.0 * 2.0 * np.pi * f_r[k] * 50.0) * 1e15)
        g_qr_mhz = float(0.5 * (c_qr_fF / np.sqrt(c_sig[k] * c_r_fF)) * np.sqrt(f_q[k] * f_r[k]) / 1e6)
        readout[k] = dict(C_qr_fF=c_qr_fF, C_r_lambda4_fF=c_r_fF, g_qr_MHz=g_qr_mhz)

    return dict(
        nets=dict(Q1=n_q1, Q2=n_q2, TC=n_tc, RO1=n_ro1, RO2=n_ro2),
        Ec_GHz=ec_ghz, C_sigma_fF=c_sig, g_n_MHz=g_n_mhz, readout=readout,
    )


def _clip_ground_pads_for_q3d(d):
    """Clip squid_ground / jj_pad_ground to each component's subtract pocket so Q3D sees flush contact, zero overlap."""
    from shapely.ops import unary_union
    tbl = d.qgeometry.tables['poly']
    for cname in ('Q1', 'Q2', 'TC'):
        if cname not in d.components:
            continue
        cid = d.components[cname].id
        subs = list(tbl[(tbl['component'] == cid) & (tbl['subtract'])]['geometry'])
        if not subs:
            continue
        sub_u = unary_union(subs)
        mask = (tbl['component'] == cid) & (~tbl['subtract']) & (tbl['name'].str.contains('squid_ground|jj_pad_ground'))
        for idx in tbl.index[mask]:
            tbl.at[idx, 'geometry'] = tbl.at[idx, 'geometry'].intersection(sub_u)


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--max-passes', type=int, default=8)
    ap.add_argument('--percent-error', type=float, default=0.5)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--overrides', default='{}',
                    help='JSON merged into build_design overrides, e.g. {"ro_claw_Q1": {"claw_sep": "12um"}, "xmon_Q1": {...}}')
    ap.add_argument('--tag', default='r1')
    a = ap.parse_args()
    extra = json.loads(a.overrides)

    ov = dict(
        xmon=dict(squid=dict(contact_overlap='0um', shift_toward_arm='0um')),
        coupler=dict(squid_contact_overlap='0um', squid_shift_toward_stem='0um'),
    )
    for k, v in extra.items():
        ov[k] = {**ov.get(k, {}), **v} if isinstance(v, dict) else v
    d = D11.build_design(variant='sim', overrides=ov, channels=False, flux_as_ground=True, check_keepout_um=100.0)
    _clip_ground_pads_for_q3d(d)
    setup = dict(SETUP, max_passes=a.max_passes, percent_error=a.percent_error)
    folder = create_run(
        'Design11_Transmon2Q_Manhattan', 'qubit', f'q3d_core_{a.tag}', entrypoint=__file__,
        inputs=[ROOT / 'Design11_Transmon2Q_Manhattan/design.json',
                ROOT / 'Design11_Transmon2Q_Manhattan/layout/build.py',
                ROOT / 'Design11_Transmon2Q_Manhattan/layout/components/round_corners.py'],
        parameters=dict(components=COMPS, open_terminations=OPENS, Q3D=setup, dry_run=a.dry_run, overrides=extra),
    )
    print('RUN_FOLDER=' + str(folder), flush=True)
    if a.dry_run:
        finish_run(folder, status='completed', summary=dict(dry_run=True))
        return

    from qiskit_metal.analyses.quantization import LOManalysis
    import hfss_dirs, save_aedt

    hfss_dirs.use_d_drive()
    # Short project/design names keep AEDT .aedtresults paths well under Windows 260-char MAX_PATH.
    project = ROOT / 'Design11_Transmon2Q_Manhattan' / 'simulation' / 'aedt' / 'D11.aedt'
    project.parent.mkdir(parents=True, exist_ok=True)
    save_aedt.use_project(project)

    lom = LOManalysis(d, 'q3d')
    lom.sim.setup.update(setup)
    t0 = time.time()
    try:
        lom.sim.run(name='D11_' + folder.name[-4:], components=COMPS, open_terminations=OPENS)
    except Exception:
        try:
            oDesk = lom.sim.renderer._pinfo.desktop._desktop
            print('=== AEDT MESSAGES ===', flush=True)
            for msg in oDesk.GetMessages('', '', 0):
                print(' ', msg, flush=True)
        except Exception as e2:
            print('Could not read AEDT messages:', e2, flush=True)
        raise
    cmat = lom.sim.capacitance_matrix
    cmat.to_csv(folder / 'raw/cmat_fF.csv')
    res = analyse_cmat(cmat)
    res['elapsed_s'] = round(time.time() - t0, 1)
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    (folder / 'summaries/summary.json').write_text(json.dumps(res, indent=2), encoding='utf-8')
    finish_run(folder, status='completed', summary=res)
    print(json.dumps(res, indent=2), flush=True)


if __name__ == '__main__':
    main()
