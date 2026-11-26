# -*- coding: utf-8 -*-
r"""HFSS Driven S21 + Q3D coupon comparison for Design11_Transmon2Q_Manhattan PurcellReadout channels.

Usage:
  python -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/resonator/hfss_readout_s21.py --mode driven
  python -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/resonator/hfss_readout_s21.py --mode q3d
  python -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/resonator/hfss_readout_s21.py --mode both

Solver defaults (skill/em_simulation.md):
  - HFSS Driven S21: max_passes=4, basis_order=1, explicit 10 um mesh (mesh_all),
    auto wire bonds every 350 um + PEC straps at PurcellReadout.strap_points().
    One wide Fast sweep (1 MHz step) only; no fine sweep (locate f_p / f_r with eigenmode).
  - Q3D coupons: max_passes=8, percent_error=0.3%, freq=6.25 GHz (cj_Q1, cj_Q2, ref, cp).
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

import numpy as np
from qiskit_metal import designs
import build as D11
import hfss_dirs
import purcell_readout_model as PRM
import purcell_resonator_pair_sim as S
import save_aedt
from enclosing_cap import EnclosingCap
from purcell_readout import PurcellReadout
from run_record import create_run, finish_run

DESIGN = 'Design11_Transmon2Q_Manhattan'
CHANNELS = ('Q2', 'Q1')  # Q2: 6.15 GHz, Q1: 6.35 GHz


def _small_design(sx: float, sy: float):
    d = designs.DesignPlanar()
    d.overwrite_enabled = True
    d.chips.main.material = 'silicon_eps11p45'
    d.chips.main.size.update(size_x=f'{sx}mm', size_y=f'{sy}mm', size_z='-500um', center_x='0mm', center_y='0mm')
    return d


def build_driven_fixture():
    d = D11.build_design(variant='sim', channels=True, flux_as_ground=True, check_keepout_um=None)
    t = d.qgeometry.tables['path']
    t['hfss_wire_bonds'] = t['subtract'].astype(bool)
    comps = ['feed_cin', 'cin', 'feed_in', 'feed_mid', 'feed_out', 'Q1', 'Q2']
    for q in ('Q1', 'Q2'):
        comps += [f'{p}_{q}' for p in ('n0', 'stub_p', 'cp', 'prr', 'res_open')]
    ports = [('feed_cin', 'start', 50), ('feed_out', 'start', 50)]
    opens = []
    return d, comps, ports, opens


def build_q3d_coupons():
    """Build local Q3D coupon designs for cj_Q1, cj_Q2, ref (straight 16/9 CPW pair), and cp (EnclosingCap)."""
    cases = {}
    for q in ('Q1', 'Q2'):
        ch = D11.PURCELL_CHANNELS[q]
        d = _small_design(2.0, 2.0)
        PurcellReadout(d, f'cj_{q}', options=dict(
            draw_mode='coupon', coupon_lead='400um',
            line_spacing='400um', tee_half='100um', lead='100um',
            head_width='25um', head_height='65um', head_gap=f"{ch['head_gap_um']:.1f}um",
            hammer_ground_gap='16um', corner_radius='5um', **D11.CPW_RO,
        ))
        cases[f'cj_{q}'] = (d, [f'cj_{q}'], [], ('filter', 'resonator'))

    d_ref = _small_design(2.0, 2.0)
    S.TwoLines(d_ref, 'ref', options=dict(spacing='400um', half_length='400um', **D11.CPW_RO))
    cases['ref'] = (d_ref, ['ref'], [], ('line_f', 'line_r'))

    d_cp = _small_design(1.5, 1.0)
    EnclosingCap(d_cp, 'cp', options=dict(**D11.ENC_OPT, **D11.CPW_RO))
    D11._round_component_polys(d_cp.components['cp'], 5e-3, skip_names=())
    cases['cp'] = (d_cp, ['cp'], [('cp', 'a'), ('cp', 'b')], ('frame|lead_a', 'plate|lead_b'))
    return cases


def analyse_s21(folder: Path, raw: dict) -> dict:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    f, s = S.load_s(folder / 'raw' / raw.get('s_file', 'S_fine.csv'))
    targets = {q: D11.PURCELL_CHANNELS[q]['f_GHz'] for q in CHANNELS}
    windows = S.channel_windows(f, s, targets)

    out = {}
    fig, ax = plt.subplots(figsize=(9.5, 4.5))
    ax.plot(f, 20 * np.log10(np.maximum(abs(s), 1e-12)), 'k-', lw=1.8, label='HFSS Driven |S21|')
    for q in CHANNELS:
        w = windows[q]
        try:
            r = S.s21_channel(f, s, w)
            ax.plot(r['f'], 20 * np.log10(np.maximum(abs(r['fit']), 1e-12)), '--', lw=1.2, label=f'{q} 2-pole fit')
            poles = r['poles']
            J = abs(poles['J_MHz'])
            kp = poles['kappa_p_MHz']
            kr = (4.0 * J * J / kp) if kp > 0 else float('nan')
            out[q] = dict(
                window_GHz=r['window_GHz'],
                zeros=r['zeros'],
                poles=poles,
                kappa_r_est_MHz=float(kr),
                target_f_GHz=targets[q],
                df_p_MHz=(poles['f_p_GHz'] - targets[q]) * 1e3,
                df_r_MHz=(poles['f_r_GHz'] - targets[q]) * 1e3,
                extra_dips=r['extra_dips'],
            )
        except Exception as e:
            out[q] = dict(window_GHz=list(w), error=str(e))
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('|S21| (dB)')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title(f"Design11_Transmon2Q_Manhattan Driven S21 ({raw['passes']} passes, dS={raw['last_delta_s']:.3g})")
    fig.tight_layout()
    (folder / 'figures').mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / 'figures' / 's21_channels.png', dpi=150)
    plt.close(fig)
    return out


def analyse_q3d_coupons(q3d_res: dict) -> dict:
    """Compute C_j, head_excess, C_p and model-predicted (f_p, f_r, J, kappa_p, kappa_r) from Q3D coupons."""
    ln = PRM.line('16/9')
    cl = ln['c_line']
    cp_c = q3d_res['cp']
    cp_mut = float(cp_c['mutual_fF'])
    cp_plate = float(cp_c['ground_fF']['plate'])
    cp_frame = float(cp_c['ground_fF']['frame'])

    # Load adopted core Q3D C_qr for claw load estimate
    dj = json.loads((ROOT / DESIGN / 'design.json').read_text(encoding='utf-8'))
    c_qr_map = {'Q1': 9.93, 'Q2': 11.03}
    try:
        q3d_path = ROOT / dj['adopted_runs']['q3d_core']['path'] / 'summaries' / 'summary.json'
        if q3d_path.exists():
            qs = json.loads(q3d_path.read_text(encoding='utf-8'))
            c_qr_map = {q: float(qs['readout'][q]['C_qr_fF']) for q in ('Q1', 'Q2')}
    except Exception:
        pass

    derived = dict(cp_mutual_fF=cp_mut, cp_plate_ground_fF=cp_plate, cp_frame_ground_fF=cp_frame, channels={})
    for q in ('Q1', 'Q2'):
        hl = S.hammer_loads(q3d_res[f'cj_{q}'], q3d_res['ref'], cl, 0.400)
        ch = D11.PURCELL_CHANNELS[q]
        # Circuit-model prediction of J and kappa_p using extracted C_j and C_p
        w0 = 2.0 * np.pi * ch['f_GHz'] * 1e9
        Z0 = ln['z0']
        L_f = ch['L_filter']
        L_r = ch['L_res']
        tap_f = 1.20
        tap_r = 1.30
        c_line = cl * 1e-12  # F/m? PRM c_line is fF/mm -> *1e-12 F/m
        # Use PRM helper if available, else closed-form tap-coupled lambda/4 formulas:
        # C_eq = (pi / 4) / (w0 * Z0) for lambda/4 open-short resonator
        C_eq_fF = (np.pi / 4.0) / (w0 * Z0) * 1e15
        sin_f = np.sin(np.pi / 2.0 * (L_f - tap_f) / L_f)
        sin_r = np.sin(np.pi / 2.0 * (L_r - tap_r) / L_r)
        J_MHz = 0.5 * (hl['cj'] / C_eq_fF) * sin_f * sin_r * ch['f_GHz'] * 1e3
        # Purcell linewidth kappa_p = (w0^2 * Z0 * C_p^2) / C_eq / (2*pi) (single-sided 50 Ohm / 2 = 25 Ohm load on feedline)
        Z_feed = 25.0  # matched 50 Ohm line on both sides of TapJunction
        kp_MHz = (w0 ** 2 * Z_feed * (cp_mut * 1e-15) ** 2) / (C_eq_fF * 1e-15) / (2.0 * np.pi) / 1e6
        kr_MHz = 4.0 * J_MHz ** 2 / kp_MHz if kp_MHz > 0 else float('nan')
        derived['channels'][q] = dict(
            head_gap_um=ch['head_gap_um'],
            C_j_fF=float(hl['cj']),
            head_excess_fF=float(hl['head_excess']),
            C_qr_fF=c_qr_map[q],
            J_q3d_model_MHz=float(J_MHz),
            kappa_p_q3d_model_MHz=float(kp_MHz),
            kappa_r_q3d_model_MHz=float(kr_MHz),
            target_C_j_fF=dj['targets']['purcell_readout'][q]['C_j_fF'],
            target_J_MHz=dj['targets']['purcell_readout'][q]['J_MHz'],
            target_kappa_p_MHz=dj['targets']['purcell_readout'][q]['kappa_p_MHz'],
            target_kappa_r_MHz=dj['targets']['purcell_readout'][q]['kappa_r_MHz'],
        )
    return derived


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['driven', 'q3d', 'both'], default='both')
    ap.add_argument('--passes', type=int, default=4)
    ap.add_argument('--band', default='5.75,6.75')
    ap.add_argument('--tag', default='r5_rcorner')
    ap.add_argument('--fine-mesh-um', type=float, default=0.0,
                    help='local length mesh (um) on the hammer tees, readout claws (qubit fork) and qubit crosses; 0 = off')
    a = ap.parse_args()
    lo, hi = map(float, a.band.split(','))

    folder = create_run(
        DESIGN, 'resonator', f'purcell_s21_{a.mode}_{a.tag}',
        entrypoint=__file__,
        inputs=[ROOT / DESIGN / 'design.json', ROOT / DESIGN / 'layout/build.py'],
        parameters=dict(
            mode=a.mode, max_passes=a.passes, band_GHz=[lo, hi], fine_mesh_um=a.fine_mesh_um,
            basis_order=1, mesh_length=S.MESH, wire_bonds=S.WB,
            purcell_channels=D11.PURCELL_CHANNELS,
        ),
    )
    print('RUN_FOLDER=' + str(folder), flush=True)

    try:
        hfss_dirs.use_d_drive()
        project = ROOT / DESIGN / 'simulation' / 'aedt' / 'D11.aedt'
        project.parent.mkdir(parents=True, exist_ok=True)
        save_aedt.use_project(project)

        summary = {}
        if a.mode in ('driven', 'both'):
            d, comps, ports, opens = build_driven_fixture()
            fm = None
            if a.fine_mesh_um > 0:
                objs = [f'{n}_{c}' for q in ('Q1', 'Q2')
                        for n, c in (('filter_tee', f'prr_{q}'), ('resonator_tee', f'prr_{q}'),
                                     ('readout_connector_arm', q), ('cross', q))]
                fm = {f'{a.fine_mesh_um:g}um': objs}
            raw = S.run_driven(
                folder, d, comps, ports, opens,
                f0=0.5 * (lo + hi), band=(lo, hi), passes=a.passes,
                fine=False, wide_step_ghz=0.001, tag='S11', fine_mesh=fm,
            )
            s21_ana = analyse_s21(folder, raw)
            summary['hfss_driven'] = dict(raw=raw, channels=s21_ana)

        if a.mode in ('q3d', 'both'):
            cases = build_q3d_coupons()
            q3d_res = {}
            for cname, (dc, ccomps, copens, ckeys) in cases.items():
                q3d_res[cname] = S.q3d_case(
                    cname, dc, ccomps, copens, ckeys, folder, folder.name[-4:],
                    max_passes=8, pct_error=0.3, freq_ghz=6.25,
                )
            summary['q3d_coupons'] = dict(derived=analyse_q3d_coupons(q3d_res))

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
