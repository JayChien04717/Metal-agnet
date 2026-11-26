# -*- coding: utf-8 -*-
r"""HFSS DrivenModal mutual inductance (M) and 3x3 crosstalk extraction for Design12.

Method (skill/mutual_inductance_analysis.md):
  - DrivenModal Z-parameters, discrete sweep 10-100 MHz (10 points), adaptive at 100 MHz.
  - Flux-line ports: 50 Ohm lumped ports at the flux_bias_line pins.
  - SQUID ports: rect_jj1 as PORT ONLY (50 Ohm, no lumped L in parallel -> open-loop EMF).
    rect_jj2 keeps its lumped Lj (Q1: 34.79 nH, Q2: 41.11 nH, TC: 14.03 nH), which carries
    no current with rect_jj1 open and makes every JJ port unambiguously identifiable by
    Im(Z_ii) / omega (TC ~14 nH < Q1 ~35 nH < Q2 ~41 nH).
  - M_ij = Im(Z_ji) / (2 * pi * f) quoted at f_min = 10 MHz (checks reciprocity Z21=Z12,
    lossless Re Z21 ~ 0, flat plateau over 10-100 MHz, and L_geo = Im(Z_jj)/w - L_j2).

Modes:
  --mode matrix : 6-port fixture rendering Q1, TC, Q2 (full 3x3 M and crosstalk matrix).
  --mode single : 2-port fixture rendering one device (Q1, Q2, or TC) for fast screening.

Environment: qmetal (C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parents[3]
DESIGN_DIR = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / 'code' / 'general'),
    str(DESIGN_DIR / 'layout'),
    str(DESIGN_DIR / 'layout' / 'components'),
]
import project_paths
project_paths.configure(__file__)

import numpy as np
import pandas as pd
from qiskit_metal import Dict, designs

import build as D12
from claw_coupler_manhattan import ClawCouplerManhattan
from fluxline import B2_FLUX_OPTIONS
from round_corners import round_component_polys, square_merged_etch_corners
from run_record import create_run, finish_run
from xmon_manhattan import XmonManhattan

PHI0 = 2.067833848e-15  # Wb
NAMES = ['Q1', 'TC', 'Q2']
LJ_OTHER_NH = dict(Q1=34.79, Q2=41.11, TC=14.03)  # rect_jj2 lumped Lj (nH)
CZ_FLUX_PHI0 = dict(Q1=0.116, TC=0.224, Q2=0.0)


def _merge_dict(base, extra):
    out = dict(base)
    for k, v in (extra or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge_dict(out[k], v)
        else:
            out[k] = v
    return out


def build_core_design(overrides=None, round_polys=True, measure_lgeo=False):
    """Build Design12 core (Q1, TC, Q2) matching build.py r4 with optional SQUID / flux-line overrides.

    Overrides keys:
      fbl        : dict merged into B2_FLUX_OPTIONS for Q1, Q2, TC
      fbl_Q1/Q2/TC: per-device flux-line overrides (merged after fbl)
      squid      : dict merged into Xmon squid options (jj_gap, squid_spacing, squid_bar_height, ...)
      squid_Q1/Q2: per-qubit squid overrides
      coupler    : dict merged into COUPLER_OPT (jj_gap, squid_spacing, squid_bar_height, ...)
      xmon / xmon_Q1 / xmon_Q2 : Xmon overrides
    """
    ov = overrides or {}
    xmon_opt = _merge_dict(D12.XMON_OPT, ov.get('xmon', {}))
    coupler_opt = _merge_dict(D12.COUPLER_OPT, ov.get('coupler', {}))
    qx = float(ov.get('qx_mm', D12.QX))

    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon_eps11p45'
    design.chips.main.size.update(dict(size_x='10mm', size_y='10mm', size_z='-500um', center_x='0mm', center_y='0mm'))
    design.variables.update(cpw_width=D12.CPW_RO['trace_width'], cpw_gap=D12.CPW_RO['trace_gap'])

    base_fbl = getattr(D12, 'R5_FLUX_OPTIONS', B2_FLUX_OPTIONS)
    base_sq = getattr(D12, 'R5_SQUID_OPTIONS', {})
    arms = {}
    for name, sx in (('Q1', -1), ('Q2', +1)):
        fbl_dict = _merge_dict(base_fbl, ov.get('fbl', {}))
        fbl_dict = _merge_dict(fbl_dict, ov.get('fbl_' + name, {}))
        fbl_dict['make_fbl'] = True
        if 'mirror' not in ov.get('fbl_' + name, {}):
            fbl_dict['mirror'] = str(sx < 0)
        fbl = Dict(**fbl_dict)

        q_ro = _merge_dict(D12.RO_CLAW, D12.RO_CLAW_BY_QUBIT[name])
        q_ro = _merge_dict(q_ro, ov.get('ro_claw', {}))
        q_ro = _merge_dict(q_ro, ov.get('ro_claw_' + name, {}))

        q_opt = _merge_dict(xmon_opt, D12.XMON_BY_QUBIT[name])
        q_opt = _merge_dict(q_opt, ov.get('xmon_' + name, {}))
        sq_base = _merge_dict(dict(XmonManhattan.default_options.squid), base_sq)
        sq_dict = _merge_dict(sq_base, q_opt.get('squid', {}))
        sq_dict = _merge_dict(sq_dict, ov.get('squid', {}))
        sq_dict = _merge_dict(sq_dict, ov.get('squid_' + name, {}))
        q_opt['squid'] = Dict(**sq_dict)

        arms[name] = design.parse_value(q_opt['cross_length'])
        XmonManhattan(design, name, options=dict(
            **D12._pos(sx * qx, D12.YQ), **q_opt,
            Lj_a='1e-6nH' if measure_lgeo else D12.LJ[name + 'a'],
            Lj_b='1e-6nH' if measure_lgeo else D12.LJ[name + 'b'],
            jj_metal='True', flux_as_ground='False',
            flux_bias_line_options=fbl,
            connection_pads=dict(readout=dict(q_ro, connector_location='90')),
        ))
        if round_polys:
            raw_etcher = D12._poly_geom(design, name, 'readout_connector_etcher')
            round_component_polys(design.components[name], 5e-3)
            square_merged_etch_corners(design.components[name], 'readout_connector_etcher', raw_etcher,
                                       D12._poly_geom(design, name, 'cross_etch'), pad=5e-3)

    tc_fbl = _merge_dict(base_fbl, ov.get('fbl', {}))
    tc_fbl = _merge_dict(tc_fbl, ov.get('fbl_TC', {}))
    tc_fbl['make_fbl'] = True
    coupler_opt['flux_bias_line_options'] = Dict(**tc_fbl)
    # Propagate common squid keys (jj_gap, squid_spacing, squid_bar_height) if not explicitly overridden in coupler
    for k in ('jj_gap', 'squid_spacing', 'squid_bar_height'):
        if k in ov.get('squid', {}) and k not in ov.get('coupler', {}):
            coupler_opt[k] = ov['squid'][k]

    claw_sep = design.parse_value(coupler_opt['claw_sep'])
    L_bar = (2.0 * qx - arms['Q1'] - arms['Q2']) - 2.0 * claw_sep
    x_tc = 0.5 * (arms['Q1'] - arms['Q2'])
    ClawCouplerManhattan(design, 'TC', options=dict(
        **D12._pos(x_tc, D12.YQ), pad_arm_length=f'{L_bar * 1e3:.4f}um',
        qubit_arm_width=xmon_opt['cross_width'],
        Lj='1e-6nH' if measure_lgeo else D12.LJ['TC'],
        jj_metal='True', flux_as_ground='False', **coupler_opt,
    ))
    if round_polys:
        round_component_polys(design.components['TC'], 5e-3)
    return design


def _loop_metrics(design, dev_name):
    """Extract geometric loop dimensions, area, perimeter, and Manhattan site coordinates."""
    c = design.components[dev_name]
    p = c.p
    if dev_name in ('Q1', 'Q2'):
        jg_um = float(p.squid.jj_gap) * 1e3
        sp_um = float(p.squid.squid_spacing) * 1e3
        bh_um = float(p.squid.squid_bar_height) * 1e3
    else:
        jg_um = float(p.jj_gap) * 1e3
        sp_um = float(p.squid_spacing) * 1e3
        bh_um = float(p.squid_bar_height) * 1e3
    # Inner crossbar-to-crossbar height = 2*y_root - 1um = jg - bh - 1um (y_root = jg/2 - bh/2)
    h_inner_um = jg_um - bh_um - 1.0
    # Centerline loop height = jg (from island crossbar center to ground crossbar center)
    h_center_um = jg_um
    # Inner width between Manhattan lead tips/trunks (approx sp - lead_tip_width at center, sp - lead_base_width at root)
    w_center_um = sp_um
    area_center_um2 = w_center_um * h_center_um
    perim_center_um = 2.0 * (w_center_um + h_center_um)
    sites = c.manhattan_sites()
    return dict(
        jj_gap_um=jg_um, squid_spacing_um=sp_um, squid_bar_height_um=bh_um,
        h_inner_um=h_inner_um, w_center_um=w_center_um, h_center_um=h_center_um,
        area_center_um2=area_center_um2, perim_center_um=perim_center_um,
        manhattan_sites=sites,
    )


def analyse_matrix(Z, freqs, devices, conv_df=None):
    """Identify ports by Im(Z_ii)/w and extract M_ij, crosstalk, L_loop, L_geo, sanity checks."""
    n_dev = len(devices)
    n_ports = 2 * n_dev
    w = 2.0 * np.pi * freqs
    L0_nH = np.array([Z[0, i, i].imag / w[0] * 1e9 for i in range(n_ports)])

    # Identify JJ ports: the n_dev ports with largest Im(Z_ii)/w, matched to devices sorted by LJ_OTHER_NH
    jj_indices = sorted(np.argsort(L0_nH)[-n_dev:], key=lambda idx: L0_nH[idx])
    devs_by_lj = sorted(devices, key=lambda d: LJ_OTHER_NH[d])
    jj_map = {dev: int(idx) for dev, idx in zip(devs_by_lj, jj_indices)}

    # Identify flux ports: remaining ports, matched by maximum |Im(Z[jj, fl])|
    fl_pool = [i for i in range(n_ports) if i not in jj_map.values()]
    fl_map = {}
    used = set()
    for dev in devices:
        cand = [i for i in fl_pool if i not in used]
        best = max(cand, key=lambda i: abs(Z[0, jj_map[dev], i].imag))
        fl_map[dev] = int(best)
        used.add(best)

    # Build M[loop][line] vs f
    M_pH = {}
    M_fmin = {}
    M_spread_pct = {}
    recip_err = {}
    re_over_im = {}
    L_loop_nH = {}
    L_geo_pH = {}
    L_fl_pH = {}

    for di in devices:
        ji = jj_map[di]
        fi = fl_map[di]
        L_loop_nH[di] = float(Z[0, ji, ji].imag / w[0] * 1e9)
        L_geo_pH[di] = float((Z[0, ji, ji].imag / w[0] * 1e9 - LJ_OTHER_NH[di]) * 1e3)
        L_fl_pH[di] = float(Z[0, fi, fi].imag / w[0] * 1e12)
        M_pH[di] = {}
        M_fmin[di] = {}
        for dj in devices:
            fj = fl_map[dj]
            z_ji = Z[:, ji, fj]
            z_ij = Z[:, fj, ji]
            zm = 0.5 * (z_ji + z_ij)
            m_arr = zm.imag / w * 1e12
            M_pH[di][dj] = m_arr
            M_fmin[di][dj] = float(m_arr[0])
            if di == dj:
                M_spread_pct[di] = float(np.ptp(m_arr) / abs(m_arr[0]) * 100.0)
                recip_err[di] = float(np.max(np.abs(z_ji - z_ij) / np.maximum(np.abs(zm), 1e-18)))
                re_over_im[di] = float(abs(zm[0].real) / max(abs(zm[0].imag), 1e-18))

    I_phi0_mA = {d: float(PHI0 / abs(M_fmin[d][d] * 1e-12) * 1e3) for d in devices}
    I_0p5_mA = {d: 0.5 * I_phi0_mA[d] for d in devices}
    I_cz_mA = {d: CZ_FLUX_PHI0[d] * I_phi0_mA[d] for d in devices if CZ_FLUX_PHI0[d] > 0}

    crosstalk_pct = {}
    for di in devices:
        for dj in devices:
            if di != dj:
                # M_ij / M_jj: flux in loop i per unit flux applied to loop j by line j
                crosstalk_pct[f'{di}<-{dj}'] = float(100.0 * M_fmin[di][dj] / M_fmin[dj][dj])

    conv_info = {}
    if conv_df is not None and len(conv_df) > 0:
        last = conv_df.iloc[-1]
        conv_info = dict(
            passes=int(len(conv_df)),
            columns=list(conv_df.columns),
            last_row={str(k): (float(v) if isinstance(v, (int, float, np.number)) else str(v)) for k, v in last.items()},
        )

    return dict(
        devices=devices,
        port_map=dict(jj={k: v + 1 for k, v in jj_map.items()}, flux={k: v + 1 for k, v in fl_map.items()}),
        L0_nH_by_port=[float(x) for x in L0_nH],
        M_fmin_pH=M_fmin,
        M_self_pH={d: M_fmin[d][d] for d in devices},
        M_spread_pct=M_spread_pct,
        I_phi0_mA=I_phi0_mA,
        I_0p5_mA=I_0p5_mA,
        I_cz_mA=I_cz_mA,
        L_loop_nH=L_loop_nH,
        L_geo_pH=L_geo_pH,
        L_fl_pH=L_fl_pH,
        crosstalk_pct=crosstalk_pct,
        reciprocity_rel_err=recip_err,
        re_over_im_Z21=re_over_im,
        f_MHz=[float(x / 1e6) for x in freqs],
        convergence=conv_info,
    ), M_pH


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tag', required=True)
    ap.add_argument('--mode', choices=('matrix', 'single'), default='matrix')
    ap.add_argument('--device', choices=('Q1', 'TC', 'Q2'), default='Q1', help='Device for --mode single')
    ap.add_argument('--min-passes', type=int, default=12)
    ap.add_argument('--max-passes', type=int, default=12)
    ap.add_argument('--min-converged', type=int, default=3)
    ap.add_argument('--max-delta-s', type=float, default=0.0002)
    ap.add_argument('--npts', type=int, default=10, help='Discrete sweep points 10-100 MHz')
    ap.add_argument('--mesh', default='', help='Optional length-based mesh seed on SQUID/flux line (e.g. 6um)')
    ap.add_argument('--overrides', default='{}', help='JSON overrides for build_core_design')
    ap.add_argument('--role', default='baseline', help='screening | baseline | convergence_check')
    ap.add_argument('--no-round', action='store_true', help='Skip round_component_polys (test HFSS subtract)')
    ap.add_argument('--measure-lgeo', action='store_true', help='Set rect_jj2 Lj=1e-6 nH to read pure geometric L_loop')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    ov = json.loads(a.overrides)
    d = build_core_design(ov, round_polys=not a.no_round, measure_lgeo=a.measure_lgeo)
    devices = NAMES if a.mode == 'matrix' else [a.device]
    comps = list(devices)
    opens = [(dev, 'readout') for dev in devices if dev in ('Q1', 'Q2')]
    if a.measure_lgeo:
        opens += [(dev, 'flux_bias_line') for dev in devices]
        ports = []
    else:
        ports = [(dev, 'flux_bias_line', 50) for dev in devices]
    jj_ports = [(dev, 'rect_jj1', 50, False) for dev in devices]
    loop_geom = {dev: _loop_metrics(d, dev) for dev in devices}

    setup_dict = dict(
        freq_ghz=0.1,
        max_delta_s=a.max_delta_s,
        max_passes=a.max_passes,
        min_passes=a.min_passes,
        min_converged=a.min_converged,
        pct_refinement=30,
        basis_order=1,
        sweep_start_ghz=0.01,
        sweep_stop_ghz=0.10,
        sweep_count=a.npts,
        sweep_type='Discrete',
        mesh_seed=a.mesh or None,
    )

    run_dir = create_run(
        'Design12_Transmon2Q_TransmissionPurcell',
        'flux',
        f'hfss_mutual_{a.tag}',
        entrypoint=__file__,
        inputs=[DESIGN_DIR / 'design.json', DESIGN_DIR / 'layout' / 'build.py'],
        parameters=dict(
            mode=a.mode, devices=devices, role=a.role,
            components=comps, open_terminations=opens,
            port_list=ports, jj_to_port=jj_ports,
            overrides=ov, setup=setup_dict, loop_geom=loop_geom,
        ),
    )
    print('RUN_FOLDER=' + str(run_dir), flush=True)
    if a.dry_run:
        finish_run(run_dir, status='completed', summary=dict(dry_run=True, loop_geom=loop_geom))
        print(json.dumps(loop_geom, indent=2))
        return

    from qiskit_metal.analyses.simulation.scattering_impedance import ScatteringImpedanceSim
    import hfss_dirs
    import save_aedt

    hfss_dirs.use_d_drive()
    # Keep project & design names very short so AEDT .aedtresults paths stay well under Windows 260-char MAX_PATH.
    project = DESIGN_DIR / 'simulation' / 'aedt' / 'D12F.aedt'
    project.parent.mkdir(parents=True, exist_ok=True)
    save_aedt.use_project(project)

    fine = []
    if a.mesh:
        tbl = d.qgeometry.tables['poly']
        for dev in devices:
            cid = d.components[dev].id
            for r in ('flux_bias_line', 'squid_ground'):
                if len(tbl[(tbl['component'] == cid) & (tbl['name'] == r)]) > 0:
                    fine.append(f'{r}_{dev}')

    class FluxSim(ScatteringImpedanceSim):
        def _analyze(self):
            if fine and a.mesh:
                self.renderer.modeler.mesh_length('squid_flux_mesh', fine, MaxLength=a.mesh)
                print(f'Seeded mesh {a.mesh} on {fine}', flush=True)
            super()._analyze()

    sim = FluxSim(d, 'hfss')
    for k, v in dict(
        freq_ghz=0.1,
        max_delta_s=a.max_delta_s,
        max_passes=a.max_passes,
        min_passes=a.min_passes,
        min_converged=a.min_converged,
        pct_refinement=30,
        basis_order=1,
        sweep_setup=Dict(
            name='Sweep', start_ghz=0.01, stop_ghz=0.10, count=a.npts,
            step_ghz=None, type='Discrete', save_fields=False,
        ),
    ).items():
        sim.setup[k] = v

    t0 = time.time()
    des_name = f'F_{run_dir.name[-4:]}'
    try:
        sim.run_sim(
            name=des_name,
            components=comps,
            open_terminations=opens,
            port_list=ports,
            jj_to_port=jj_ports,
            box_plus_buffer=True,
        )
    except Exception:
        try:
            oDesk = sim.renderer._pinfo.desktop._desktop
            print('=== AEDT MESSAGES ===', flush=True)
            for msg in oDesk.GetMessages('', '', 0):
                print(' ', msg, flush=True)
        except Exception as e2:
            print('Could not read AEDT messages:', e2, flush=True)
        finish_run(run_dir, status='failed', summary=dict(error='HFSS Analyze failed'))
        raise

    if a.measure_lgeo:
        freqs, curves, _ = sim.renderer.get_params(['Z11'])
        freqs = np.asarray(freqs, dtype=float)
        z11 = np.asarray(curves[0], dtype=complex)
        w = 2.0 * np.pi * freqs
        l_geo_pH = float(z11[0].imag / w[0] * 1e12)
        summary = dict(
            tag=a.tag, role=a.role, mode='measure_lgeo', devices=devices,
            L_loop_nH={devices[0]: l_geo_pH * 1e-3}, L_geo_pH={devices[0]: l_geo_pH},
            L_geo_spread_pct=float(np.ptp(z11.imag / w) / abs(z11[0].imag / w[0]) * 100.0),
            elapsed_s=round(time.time() - t0, 1), loop_geom=loop_geom, overrides=ov, setup=setup_dict,
        )
        (run_dir / 'summaries').mkdir(parents=True, exist_ok=True)
        (run_dir / 'summaries' / 'summary.json').write_text(json.dumps(summary, indent=2, default=float), encoding='utf-8')
        finish_run(run_dir, status='completed', summary=summary)
        print(json.dumps(summary, indent=2), flush=True)
        return

    n_ports = 2 * len(devices)
    pnames = [f'Z{i}{j}' for i in range(1, n_ports + 1) for j in range(1, n_ports + 1)]
    freqs, curves, _ = sim.renderer.get_params(pnames)
    freqs = np.asarray(freqs, dtype=float)
    Z = np.zeros((len(freqs), n_ports, n_ports), dtype=complex)
    for nm, cv in zip(pnames, curves):
        Z[:, int(nm[1]) - 1, int(nm[2]) - 1] = np.asarray(cv, dtype=complex)

    raw_dir = run_dir / 'raw'
    raw_dir.mkdir(parents=True, exist_ok=True)
    np.save(raw_dir / 'Z.npy', Z)
    np.save(raw_dir / 'f.npy', freqs)

    conv_df = None
    try:
        conv_df, _ = sim.renderer.pinfo.setup.get_convergence()
        conv_df.to_csv(raw_dir / 'convergence.csv')
    except Exception as ex:
        print('convergence table not saved:', repr(ex), flush=True)

    summary, M_pH = analyse_matrix(Z, freqs, devices, conv_df=conv_df)
    summary['elapsed_s'] = round(time.time() - t0, 1)
    summary['role'] = a.role
    summary['tag'] = a.tag
    summary['loop_geom'] = loop_geom
    summary['overrides'] = ov
    summary['setup'] = setup_dict

    # Save Z.csv
    rows = []
    w = 2.0 * np.pi * freqs
    for k, fk in enumerate(freqs):
        row = dict(f_MHz=fk / 1e6)
        for di in devices:
            for dj in devices:
                row[f'M_{di}<-{dj}_pH'] = float(M_pH[di][dj][k])
            ji = summary['port_map']['jj'][di] - 1
            fi = summary['port_map']['flux'][di] - 1
            row[f'L_loop_{di}_nH'] = float(Z[k, ji, ji].imag / w[k] * 1e9)
            row[f'L_fl_{di}_pH'] = float(Z[k, fi, fi].imag / w[k] * 1e12)
            row[f'ReZ_{di}_ohm'] = float(Z[k, ji, fi].real)
            row[f'ImZ_{di}_ohm'] = float(Z[k, ji, fi].imag)
        rows.append(row)
    pd.DataFrame(rows).to_csv(raw_dir / 'Z.csv', index=False)

    (run_dir / 'summaries').mkdir(parents=True, exist_ok=True)
    (run_dir / 'summaries' / 'summary.json').write_text(json.dumps(summary, indent=2, default=float), encoding='utf-8')

    # Plot M(f) and crosstalk
    (run_dir / 'figures').mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2 if len(devices) > 1 else 1, figsize=(11 if len(devices) > 1 else 6, 4))
    ax0 = ax[0] if len(devices) > 1 else ax
    for dev in devices:
        ax0.plot(freqs / 1e6, M_pH[dev][dev], 'o-', label=f'{dev}: {summary["M_self_pH"][dev]:.3f} pH')
    ax0.set_xlabel('f (MHz)')
    ax0.set_ylabel('M (pH)')
    ax0.set_title(f'Self mutual inductance ({a.tag})')
    ax0.grid(True, alpha=0.3)
    ax0.legend()
    if len(devices) > 1:
        for di in devices:
            for dj in devices:
                if di != dj:
                    xt = 100.0 * M_pH[di][dj] / M_pH[dj][dj]
                    ax[1].plot(freqs / 1e6, xt, 'o-', label=f'{di}<-{dj}: {xt[0]:+.3f}%')
        ax[1].set_xlabel('f (MHz)')
        ax[1].set_ylabel('M_ij / M_jj (%)')
        ax[1].set_title(f'Flux crosstalk ({a.tag})')
        ax[1].grid(True, alpha=0.3)
        ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(run_dir / 'figures' / 'mutual_z.png', dpi=140)
    plt.close(fig)

    finish_run(run_dir, status='completed', summary=summary)
    print(json.dumps({k: summary[k] for k in (
        'tag', 'role', 'elapsed_s', 'M_self_pH', 'M_spread_pct', 'I_phi0_mA',
        'I_0p5_mA', 'I_cz_mA', 'L_loop_nH', 'L_geo_pH', 'L_fl_pH',
        'crosstalk_pct', 'reciprocity_rel_err', 're_over_im_Z21', 'convergence',
    )}, indent=2), flush=True)


if __name__ == '__main__':
    main()
