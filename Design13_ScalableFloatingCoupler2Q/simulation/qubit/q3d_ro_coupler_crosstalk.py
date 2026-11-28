# -*- coding: utf-8 -*-
r"""Q3D check of readout-arc <-> coupler-C-arc parasitic coupling (Design13 r6 circular design language).

Model (sim geometry, flux lines drawn as ground, JJ landing metal omitted):
- Q0: D5-like circular qubit at (0, 0) with readout arc at NE (open-terminated pin -> own net),
  JJ at SW; the XY stub is not drawn (it is grounded at the qubit frequency).
- 4 C-arc floating couplers around Q0: TC_E (ori 0), TC_N (90), TC_W (180), TC_S (270);
  pad_l of every coupler is the arc next to Q0, pad_r the arc next to the neighbour qubit.
- 4 neighbour qubits Q_E / Q_N / Q_W / Q_S (no readout arc) so the far arcs are loaded realistically.

Quantities (skill/cap_matrix_analysis.md):
- direct mutual capacitances C(RO, pad_l_k), C(RO, Q0), C(Q0, pad_l_k)
- coupler mode = difference coordinate of (pad_l, pad_r); Ec of Q0 and couplers from the inverse matrix
- linear coupling g_ij = 0.5 * C_ij,eff / sqrt(C_i C_j) * sqrt(f_i f_j), with C_ij,eff = -C~_ij from the
  transformed Maxwell matrix and C_r = pi / (4 w_r Z0) for the lambda/4 readout resonator
  (the readout node is the arc pad only, so g_qr / g_rc are pad-level estimates; their RATIO is the robust number).
"""
from __future__ import annotations

from argparse import ArgumentParser
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DESIGN_DIR = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'aedt'), str(DESIGN_DIR / 'layout')]
import project_paths
project_paths.configure(__file__)

import build as D13
from run_record import create_run, finish_run

TCS =[('TC_E', 1.0, 0.0, '0'), ('TC_N', 0.0, 1.0, '90'), ('TC_W', -1.0, 0.0, '180'), ('TC_S', 0.0, -1.0, '270')]
NEIGH = [('Q_E', 2.0, 0.0), ('Q_N', 0.0, 2.0), ('Q_W', -2.0, 0.0), ('Q_S', 0.0, -2.0)]
SETUP = dict(min_passes=2, max_passes=8, percent_error=0.5, freq_ghz=5.0,
             solution_order='Medium', auto_increase_solution_order=False)
E = 1.602176634e-19
H = 6.62607015e-34
F_Q = 4.65e9      # D5 target
F_C = 6.2e9       # coupler idle target (6.0-6.4 GHz)
F_R = 7.17e9      # D5 readout target
Z0 = 50.0


def build_model(q_over=None, c_over=None):
    from qiskit_metal import designs
    d = designs.DesignPlanar()
    d.overwrite_enabled = True
    d.chips.main.size.update(dict(size_x='6mm', size_y='6mm', size_z='-500um', center_x='0mm', center_y='0mm'))
    d.chips.main.material = 'silicon_eps11p45'
    q_opt = dict(D13.QUBIT_OPT, jj_metal='False', **(q_over or {}))
    D13.ScalablePadQubit(d, 'Q0', options=dict(**D13._pos(0.0, 0.0), make_readout_pad='True', ro_exit_mode='diag',
                                               ro_corner='NE', jj_corner='SW', **q_opt))
    for nm, x, y in NEIGH:
        D13.ScalablePadQubit(d, nm, options=dict(**D13._pos(x, y), make_readout_pad='False', ro_corner='NE',
                                                 jj_corner='SW', **q_opt))
    c_opt = dict(D13.COUPLER_OPT, jj_metal='False', flux_as_ground='True', **(c_over or {}))
    for nm, x, y, ori in TCS:
        D13.FloatingTunableCoupler(d, nm, options=dict(**D13._pos(x, y), orientation=ori, **c_opt))
    return d


def clip_ground_stems(d):
    """Make ground-connected metal Q3D-safe.

    - JJ ground stems: clip to the owning pocket (flush contact with ground, no overlap).
    - Readout shield spokes: they are unetched ground, so subtract them from EVERY pocket (they also cross the
      neighbouring coupler-arc etch) and drop the overlay polygon -> identical metal, no extra net."""
    from shapely.ops import unary_union
    tbl = d.qgeometry.tables['poly']
    spokes = []
    for cname in ['Q0'] + [n for n, *_ in NEIGH]:
        cid = d.components[cname].id
        own = unary_union(list(tbl[(tbl['component'] == cid) & (tbl['subtract'])]['geometry']))
        for idx in tbl.index[(tbl['component'] == cid) & (~tbl['subtract']) & (tbl['name'] == 'jj_gnd_stem')]:
            tbl.at[idx, 'geometry'] = tbl.at[idx, 'geometry'].intersection(own)
        sp_idx = list(tbl.index[(tbl['component'] == cid) & (~tbl['subtract']) & (tbl['name'].str.startswith('gnd_spoke'))])
        spokes += [tbl.at[i, 'geometry'] for i in sp_idx]
        tbl.drop(index=sp_idx, inplace=True)
    if spokes:
        sp = unary_union(spokes)
        for idx in tbl.index[(tbl['subtract']) & (tbl['layer'].astype(int) == 1)]:
            g = tbl.at[idx, 'geometry']
            if g.intersects(sp):
                tbl.at[idx, 'geometry'] = g.difference(sp)


def cut_diag_endcap(d, comp='Q0', pin='readout'):
    """Open-terminate a (45-deg) pin by cutting a rotated endcap into the component pocket.

    qiskit-metal's Q3D open_terminations only draws axis-aligned endcaps, which leave the corners of a
    diagonal stem touching ground (the readout pad then merges into ground_main_plane)."""
    from shapely.geometry import Polygon
    pd_ = d.components[comp].pins[pin]
    w, g = float(d.parse_value(pd_['width'])), float(d.parse_value(pd_['gap']))
    mid, nrm = np.asarray(pd_['middle'], float), np.asarray(pd_['normal'], float)
    tng = np.array([-nrm[1], nrm[0]])
    c0, c1 = mid - 0.001 * nrm, mid + g * nrm
    half = 0.5 * w + g
    cap = Polygon([tuple(c0 + half * tng), tuple(c1 + half * tng), tuple(c1 - half * tng), tuple(c0 - half * tng)])
    tbl = d.qgeometry.tables['poly']
    cid = d.components[comp].id
    idx = tbl.index[(tbl['component'] == cid) & (tbl['subtract']) & (tbl['name'] == 'pocket')][0]
    tbl.at[idx, 'geometry'] = tbl.at[idx, 'geometry'].union(cap)


def analyse(cmat):
    cols = list(cmat.columns)

    def net(*keys):
        for c in cols:
            cl = c.lower()
            if all(k.lower() in cl for k in keys):
                return c
        raise KeyError(f'{keys} not in {cols}')

    n_q0 = net('pad', 'q0')
    n_ro = net('readout', 'q0')
    pads = {tc: (net('pad_l', tc), net('pad_r', tc)) for tc, *_ in TCS}
    n_nb = {nm: net('pad', nm) for nm, *_ in NEIGH}
    nodes = [n_q0, n_ro] + [p for tc in pads for p in pads[tc]] + list(n_nb.values())
    C = cmat.loc[nodes, nodes].to_numpy(float) * 1e-15
    n = len(nodes)
    M = np.eye(n)
    for k, tc in enumerate(pads):        # couplers: rows (sum, diff) for (pad_l, pad_r)
        i = 2 + 2 * k
        M[i, i], M[i, i + 1], M[i + 1, i], M[i + 1, i + 1] = 1.0, 1.0, 1.0, -1.0
    Mi = np.linalg.inv(M)
    Ct = Mi.T @ C @ Mi
    Ci = np.linalg.inv(Ct)
    idx = {'Q0': 0, 'RO': 1}
    for k, tc in enumerate(pads):
        idx[tc] = 2 + 2 * k + 1           # difference coordinate = coupler mode

    def csig(i):
        return 1.0 / Ci[i, i]

    C_r = np.pi / (4.0 * 2.0 * np.pi * F_R * Z0)
    out = dict(nets=dict(Q0=n_q0, RO=n_ro, pads=pads), C_r_lambda4_fF=C_r * 1e15)
    out['Ec_MHz'] = {k: float(E ** 2 * Ci[i, i] / (2 * H) / 1e6) for k, i in idx.items() if k != 'RO'}
    out['C_sigma_fF'] = {k: float(csig(i) * 1e15) for k, i in idx.items() if k != 'RO'}
    out['direct_mutual_fF'] = {
        'RO-Q0': float(-cmat.loc[n_ro, n_q0]),
        **{f'RO-{tc}.pad_l': float(-cmat.loc[n_ro, pads[tc][0]]) for tc in pads},
        **{f'RO-{tc}.pad_r': float(-cmat.loc[n_ro, pads[tc][1]]) for tc in pads},
        **{f'Q0-{tc}.pad_l': float(-cmat.loc[n_q0, pads[tc][0]]) for tc in pads},
    }

    def g_lin(i, j, fi, fj, Cj=None):
        c_eff = -Ct[i, j]
        Ci_s = csig(i)
        Cj_s = Cj if Cj is not None else csig(j)
        return float(0.5 * c_eff / np.sqrt(Ci_s * Cj_s) * np.sqrt(fi * fj) / 1e6)

    out['g_MHz'] = {'Q0-RO (g_qr)': g_lin(idx['Q0'], idx['RO'], F_Q, F_R, C_r)}
    for tc in pads:
        out['g_MHz'][f'{tc}-RO (g_cr)'] = g_lin(idx[tc], idx['RO'], F_C, F_R, C_r)
        out['g_MHz'][f'Q0-{tc} (g_qc)'] = g_lin(idx['Q0'], idx[tc], F_Q, F_C)
    g_qr = out['g_MHz']['Q0-RO (g_qr)']
    out['ratio_g_cr_over_g_qr'] = {tc: float(out['g_MHz'][f'{tc}-RO (g_cr)'] / g_qr) for tc in pads}
    # dispersive (Jaynes-Cummings, linear) shift of the resonator from each coupler: g^2 / Delta
    out['resonator_pull_from_coupler_MHz'] = {
        tc: float(out['g_MHz'][f'{tc}-RO (g_cr)'] ** 2 / ((F_R - F_C) / 1e6)) for tc in pads}
    # transmon dispersive shift (alpha ~ -Ec): chi = g^2 alpha / (Delta (Delta + alpha)), n_crit = (Delta / 2g)^2
    ec_q = out['Ec_MHz']['Q0']
    d_qr = (F_Q - F_R) / 1e6
    chi = g_qr ** 2 * (-ec_q) / (d_qr * (d_qr - ec_q))
    out['readout'] = dict(g_qr_MHz=g_qr, Delta_MHz=d_qr, alpha_est_MHz=-ec_q, chi_MHz=float(chi),
                          two_chi_MHz=float(2 * chi), n_crit=float((d_qr / (2 * g_qr)) ** 2))
    # sanity checks (skill/cap_matrix_analysis.md section 7)
    Cf = cmat.loc[nodes, nodes].to_numpy(float)
    out['sanity'] = dict(
        max_asym_rel=float(np.max(np.abs(Cf - Cf.T)) / np.max(np.abs(Cf))),
        diag_positive=bool(np.all(np.diag(Cf) > 0)),
        offdiag_nonpositive=bool(np.all(Cf[~np.eye(n, dtype=bool)] <= 1e-6)),
        min_row_sum_fF=float(np.min(Cf.sum(axis=1))),
    )
    return out


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--tag', default='r6')
    ap.add_argument('--max-passes', type=int, default=SETUP['max_passes'])
    ap.add_argument('--qubit-overrides', default='{}')
    ap.add_argument('--coupler-overrides', default='{}')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()
    q_over, c_over = json.loads(a.qubit_overrides), json.loads(a.coupler_overrides)
    setup = dict(SETUP, max_passes=a.max_passes)

    d = build_model(q_over, c_over)
    clip_ground_stems(d)
    cut_diag_endcap(d, 'Q0', 'readout')
    comps = ['Q0'] + [n for n, *_ in NEIGH] + [n for n, *_ in TCS]
    opens = []   # readout pin opened by cut_diag_endcap (renderer endcaps are axis-aligned only)
    folder = create_run(
        'Design13_ScalableFloatingCoupler2Q', 'qubit', f'q3d_ro_coupler_xtalk_{a.tag}', entrypoint=__file__,
        inputs=[DESIGN_DIR / 'design.json', DESIGN_DIR / 'layout' / 'build.py'],
        parameters=dict(revision=D13.REVISION, components=comps, open_terminations=opens, Q3D=setup,
                        qubit_overrides=q_over, coupler_overrides=c_over, dry_run=a.dry_run,
                        f_targets_GHz=dict(qubit=F_Q / 1e9, coupler=F_C / 1e9, readout=F_R / 1e9),
                        note='Stage 2 check requested by user: readout arc <-> coupler C-arc parasitic coupling'),
    )
    print('RUN_FOLDER=' + str(folder), flush=True)
    if a.dry_run:
        finish_run(folder, status='completed', summary=dict(dry_run=True))
        return

    from qiskit_metal.analyses.quantization import LOManalysis
    import hfss_dirs
    import save_aedt

    hfss_dirs.use_d_drive()
    project = DESIGN_DIR / 'simulation' / 'aedt' / 'D13Q.aedt'
    project.parent.mkdir(parents=True, exist_ok=True)
    save_aedt.use_project(project)

    lom = LOManalysis(d, 'q3d')
    lom.sim.setup.update(setup)
    t0 = time.time()
    try:
        lom.sim.run(name='X_' + folder.name[-4:], components=comps, open_terminations=opens)
    except Exception:
        try:
            desk = lom.sim.renderer._pinfo.desktop._desktop
            for msg in desk.GetMessages('', '', 0):
                print('  AEDT:', msg, flush=True)
        except Exception as exc:
            print('Could not read AEDT messages:', exc, flush=True)
        finish_run(folder, status='failed', summary=dict(error='Q3D run failed'))
        raise

    cmat = lom.sim.capacitance_matrix
    (folder / 'raw').mkdir(parents=True, exist_ok=True)
    cmat.to_csv(folder / 'raw' / 'cmat_fF.csv')
    conv = None
    try:
        conv = lom.sim.get_convergence()
    except Exception:
        pass
    if conv is not None:
        try:
            conv.to_csv(folder / 'raw' / 'convergence.csv')
        except Exception:
            pass
    res = analyse(cmat)
    res['elapsed_s'] = round(time.time() - t0, 1)
    res['Q3D_setup'] = setup
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    (folder / 'summaries' / 'summary.json').write_text(json.dumps(res, indent=2, default=str), encoding='utf-8')
    finish_run(folder, status='completed', validation='single run; pass-count convergence not separately verified',
               summary=res)
    print(json.dumps(res, indent=2, default=str), flush=True)


if __name__ == '__main__':
    main()
