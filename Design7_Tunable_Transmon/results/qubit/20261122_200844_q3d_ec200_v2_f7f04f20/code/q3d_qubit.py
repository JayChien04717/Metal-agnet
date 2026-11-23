"""Single-qubit Q3D capacitance extraction for Design7_Tunable_Transmon.

Extracts Ec, C_sigma, readout beta, g (lambda/4 resonator), readout open-end load C_end,
and exact transmon Hamiltonian parameters (f01, alpha, 2*chi) for a symmetric SQUID.

Usage (from repository root, qmetal Python, AEDT running):
    python -X utf8 -B Design7_Tunable_Transmon/simulation/qubit/q3d_qubit.py [--tag baseline]
"""
from argparse import ArgumentParser
from pathlib import Path
import json
import math
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

REPO = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(REPO / 'code' / 'general')]
import project_paths  # noqa: E402

project_paths.configure(__file__)
DESIGN = project_paths.design_dir(__file__)
sys.path[:0] = [
    str(DESIGN / 'layout'),
    str(DESIGN / 'layout' / 'components'),
    *project_paths.pattern_design_paths(),
]

from qiskit_metal import Dict, designs  # noqa: E402
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround  # noqa: E402
from tunable_transmon_pocket import TunableTransmonPocket  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

E_CHARGE = 1.602176634e-19
H_PLANCK = 6.62607015e-34
HBAR = H_PLANCK / (2.0 * math.pi)
PHI0 = H_PLANCK / (2.0 * E_CHARGE)


def solve_transmon_spectrum(ej_ghz, ec_ghz, n_cut=25):
    """Exact charge-basis diagonalization of H = 4 Ec n^2 - Ej cos(phi)."""
    ns = np.arange(-n_cut, n_cut + 1, dtype=float)
    dim = len(ns)
    H = np.diag(4.0 * ec_ghz * ns ** 2)
    off = -0.5 * ej_ghz * np.ones(dim - 1)
    H += np.diag(off, 1) + np.diag(off, -1)
    evals, evecs = np.linalg.eigh(H)
    f01 = float(evals[1] - evals[0])
    f12 = float(evals[2] - evals[1])
    alpha = f12 - f01
    n_op = np.diag(ns)
    n01 = float(abs(evecs[:, 0].conj() @ n_op @ evecs[:, 1]))
    return dict(f01_ghz=f01, f12_ghz=f12, alpha_ghz=alpha, n01=n01)


def solve_ej_for_f01(target_f01_ghz, ec_ghz):
    """Bisection for EJ(0) [GHz] giving exact transmon f01 = target_f01_ghz."""
    lo, hi = 1.0, 40.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if solve_transmon_spectrum(mid, ec_ghz)['f01_ghz'] < target_f01_ghz:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def analyze_cmat(cmat_df, f_q_target_ghz=4.0, f_r_list_ghz=(6.3, 6.5, 6.7, 6.9), z0=50.9):
    """Extract Ec, C_sigma, beta, g (lambda/4), 2*chi, and readout pad self-capacitance."""
    cols = list(cmat_df.columns)
    pad_top = [c for c in cols if 'pad_top' in c][0]
    pad_bot = [c for c in cols if 'pad_bot' in c][0]
    ro_pad = [c for c in cols if 'readout_connector_pad' in c][0]

    nodes = [pad_top, pad_bot, ro_pad]
    C = cmat_df.loc[nodes, nodes].to_numpy(dtype=float) * 1e-15  # F

    # Sum / difference transform on (pad_top, pad_bot); readout pad unchanged
    M = np.array([
        [1.0,  1.0, 0.0],
        [1.0, -1.0, 0.0],
        [0.0,  0.0, 1.0],
    ])
    Mi = np.linalg.inv(M)
    Ct = Mi.T @ C @ Mi
    Ci = np.linalg.inv(Ct)

    # Index 1 is the differential qubit mode
    c_sigma_f = 1.0 / Ci[1, 1]
    ec_ghz = (E_CHARGE ** 2 / (2.0 * H_PLANCK) * Ci[1, 1]) * 1e-9
    beta = float(abs(Ci[1, 2] / Ci[2, 2]))

    # Readout pad loading capacitance to ground + qubit (Maxwell C_rr)
    c_ro_self_fF = float(cmat_df.loc[ro_pad, ro_pad])
    c_ro_top_fF = float(-cmat_df.loc[ro_pad, pad_top])
    c_ro_bot_fF = float(-cmat_df.loc[ro_pad, pad_bot])

    # Symmetric SQUID EJ(0) for f01 = f_q_target_ghz
    ej_sum_ghz = solve_ej_for_f01(f_q_target_ghz, ec_ghz)
    spec = solve_transmon_spectrum(ej_sum_ghz, ec_ghz)
    ej_branch_ghz = ej_sum_ghz / 2.0
    # Lj = Phi0^2 / (4 pi^2 h EJ)
    lj_branch_nH = (PHI0 ** 2 / (4.0 * math.pi ** 2 * H_PLANCK * (ej_branch_ghz * 1e9))) * 1e9

    readout = {}
    for fr in f_r_list_ghz:
        wr = 2.0 * math.pi * fr * 1e9
        v_zpf_qtr = wr * math.sqrt(2.0 * HBAR * z0 / math.pi)  # lambda/4 open-end V_zpf
        g_mhz = (2.0 * E_CHARGE * beta * v_zpf_qtr * spec['n01'] / H_PLANCK) * 1e-6
        delta_ghz = spec['f01_ghz'] - fr
        # Dispersive shift 2*chi (full |0>-to-|1> frequency separation of resonator)
        two_chi_mhz = 2.0 * (g_mhz ** 2 / 1e3) * (spec['alpha_ghz']) / (delta_ghz * (delta_ghz + spec['alpha_ghz']))
        readout[f'{fr:.2f}GHz'] = dict(
            f_r_ghz=fr,
            v_zpf_uV=v_zpf_qtr * 1e6,
            g_mhz=g_mhz,
            two_chi_mhz=two_chi_mhz,
        )

    # Sanity checks
    sym_err = float(np.max(np.abs(C - C.T)) / np.max(np.abs(C)))
    row_sums_fF = {n: float(cmat_df.loc[n].sum()) for n in nodes}

    return dict(
        nodes=nodes,
        Ec_GHz=ec_ghz,
        C_sigma_fF=c_sigma_f * 1e15,
        beta=beta,
        C_ro_self_fF=c_ro_self_fF,
        C_ro_top_fF=c_ro_top_fF,
        C_ro_bot_fF=c_ro_bot_fF,
        symmetric_squid=dict(
            f01_target_GHz=f_q_target_ghz,
            EJ_sum_GHz=ej_sum_ghz,
            EJ_branch_GHz=ej_branch_ghz,
            Lj_branch_nH=lj_branch_nH,
            f01_GHz=spec['f01_ghz'],
            alpha_GHz=spec['alpha_ghz'],
            n01=spec['n01'],
        ),
        readout_lambda4=readout,
        sanity=dict(symmetry_rel_err=sym_err, row_sums_fF=row_sums_fF),
    )


def build_single_qubit(tp_overrides=None):
    geo = json.loads((DESIGN / 'parameters' / 'geometry.json').read_text(encoding='utf-8'))
    tp = dict(geo['transmon_pocket'])
    if tp_overrides:
        tp.update(tp_overrides)
    cpw = geo['cpw']

    d = designs.DesignPlanar(overwrite_enabled=True)
    d.chips.main.size.size_x = '2.5mm'
    d.chips.main.size.size_y = '2.5mm'
    d.chips.main.size.size_z = '-500um'
    d.chips.main.material = 'silicon'

    def _um(v):
        return f'{v:g}um'

    q1 = TunableTransmonPocket(d, 'Q1', options=Dict(
        chip='main', pos_x='0mm', pos_y='0mm', orientation='0',
        pad_width=_um(tp['pad_width_um']),
        pad_height=_um(tp['pad_height_um']),
        pad_gap=_um(tp['pad_gap_um']),
        pad_corner_radius=_um(tp['pad_corner_radius_um']),
        pocket_width=_um(tp['pocket_width_um']),
        pocket_height=_um(tp['pocket_height_um']),
        pocket_corner_radius=_um(tp['pocket_corner_radius_um']),
        squid_width=_um(tp['squid_width_um']),
        squid_loop_width=_um(tp['squid_loop_width_um']),
        jj_cut=_um(tp['jj_cut_um']),
        jj_width=_um(tp['jj_width_um']),
        jj_metal='True',
        fab_contacts='False',
        connection_pads=Dict(
            readout=Dict(
                pad_gap=_um(tp['readout_pad_gap_um']),
                pad_width=_um(tp['readout_pad_width_um']),
                pad_height=_um(tp['readout_pad_height_um']),
                cpw_width=_um(cpw['trace_width_um']),
                cpw_gap=_um(cpw['trace_gap_um']),
                cpw_extend=_um(tp['readout_cpw_extend_um']),
                loc_W='0', loc_H='+1',
            )
        ),
    ))
    pin_xy = q1.pins['readout'].middle
    OpenToGround(d, 'ro_open', options=Dict(
        chip='main',
        pos_x=f'{pin_xy[0]:.9f}mm',
        pos_y=f'{pin_xy[1]:.9f}mm',
        orientation='90',
        width=_um(cpw['trace_width_um']),
        gap=_um(cpw['trace_gap_um']),
        termination_gap=_um(cpw['trace_gap_um']),
    ))
    return d, tp


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--tag', default='q3d_q1')
    ap.add_argument('--max-passes', type=int, default=8)
    ap.add_argument('--percent-error', type=float, default=0.5)
    ap.add_argument('--pad-width', type=float, default=None)
    ap.add_argument('--pad-height', type=float, default=None)
    ap.add_argument('--pad-gap', type=float, default=None)
    ap.add_argument('--pocket-width', type=float, default=None)
    ap.add_argument('--pocket-height', type=float, default=None)
    ap.add_argument('--ro-width', type=float, default=None)
    ap.add_argument('--ro-gap', type=float, default=None)
    a = ap.parse_args()

    ov = {}
    for k, v in (
        ('pad_width_um', a.pad_width),
        ('pad_height_um', a.pad_height),
        ('pad_gap_um', a.pad_gap),
        ('pocket_width_um', a.pocket_width),
        ('pocket_height_um', a.pocket_height),
        ('readout_pad_width_um', a.ro_width),
        ('readout_pad_gap_um', a.ro_gap),
    ):
        if v is not None:
            ov[k] = v

    d, tp = build_single_qubit(ov)
    setup = dict(freq_ghz=5.0, max_passes=a.max_passes, min_passes=2,
                 percent_error=a.percent_error, solution_order='Medium',
                 auto_increase_solution_order=False)

    folder = create_run(
        DESIGN.name, 'qubit', a.tag, entrypoint=__file__,
        inputs=[DESIGN / 'design.json', DESIGN / 'parameters' / 'geometry.json',
                DESIGN / 'layout' / 'components' / 'tunable_transmon_pocket.py'],
        parameters=dict(transmon_pocket=tp, q3d_setup=setup, components=['Q1'],
                        open_terminations=[('Q1', 'readout')]),
    )
    print('RUN_FOLDER=' + str(folder), flush=True)

    try:
        from qiskit_metal.analyses.quantization import LOManalysis
        import save_aedt
        import hfss_dirs

        temp = hfss_dirs.use_d_drive()
        # Use short project/design name (D7.aedt) so AEDT .aedtresults paths stay below Windows 260-char MAX_PATH
        proj = DESIGN / 'simulation' / 'aedt' / 'D7.aedt'
        save_aedt.use_project(str(proj))

        lom = LOManalysis(d, 'q3d')
        lom.sim.setup.update(setup)
        t0 = time.time()
        dname = f'Q_{folder.name[-4:]}'
        try:
            lom.sim.run(name=dname, components=['Q1', 'ro_open'],
                        open_terminations=[('Q1', 'readout')])
        except Exception:
            try:
                desk = save_aedt._desktop()._desktop
                for lvl in (2, 1, 0):
                    msgs = desk.GetMessages(Path(proj).stem, dname + '_q3d', lvl)
                    if msgs:
                        print(f'AEDT messages (lvl={lvl}):\n' + '\n'.join(msgs), flush=True)
                        break
            except Exception as m_exc:
                print('Could not read AEDT messages:', m_exc, flush=True)
            raise
        elapsed = time.time() - t0

        cmat = lom.sim.capacitance_matrix
        try:
            save_aedt._desktop().get_active_project()._project.Save()
        except Exception:
            pass
        cmat.to_csv(folder / 'raw' / 'cmat_fF.csv')
        res = analyze_cmat(cmat)
        res['elapsed_s'] = round(elapsed, 1)
        res['q3d_setup'] = setup
        res['transmon_pocket'] = tp
        res['solver_temp'] = temp['after']

        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(
            json.dumps(res, indent=2) + '\n', encoding='utf-8'
        )
        finish_run(folder, status='completed', summary=res)

        print(f"\n=== Q3D Result ({elapsed:.1f} s) ===")
        print(f"  Ec = {res['Ec_GHz']*1e3:.2f} MHz  (C_sigma = {res['C_sigma_fF']:.2f} fF)")
        print(f"  beta = {res['beta']:.5f}  (C_ro_self = {res['C_ro_self_fF']:.2f} fF, "
              f"C_ro_top = {res['C_ro_top_fF']:.2f} fF, C_ro_bot = {res['C_ro_bot_fF']:.2f} fF)")
        sq = res['symmetric_squid']
        print(f"  Symmetric SQUID (f01=4.00 GHz): EJ_sum={sq['EJ_sum_GHz']:.3f} GHz, "
              f"EJ_branch={sq['EJ_branch_GHz']:.3f} GHz, Lj_branch={sq['Lj_branch_nH']:.2f} nH, "
              f"alpha={sq['alpha_GHz']*1e3:.1f} MHz")
        for k, r in res['readout_lambda4'].items():
            print(f"  {k}: g = {r['g_mhz']:.2f} MHz, 2*chi = {r['two_chi_mhz']:+.3f} MHz")
    except Exception as exc:
        finish_run(folder, status='failed', summary=dict(error=repr(exc)))
        raise


if __name__ == '__main__':
    main()
