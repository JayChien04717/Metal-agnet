# -*- coding: utf-8 -*-
"""Design10_CPW_BPF_5p5_6p5GHz HFSS Driven verification (2.0 - 7.5 GHz) for the
single-port (1-port) 3-stage 5.5-6.5 GHz reflection CPW bandpass filter (>40 dB in 2-4 GHz).

Ports:
    - Physical Port 1: ('feed_in', 'start', 50)   -- 50 Ohm single reflection input port
    - Probe Port 2:    ('R3', 'end', 50000)       -- non-perturbing 50 kOhm voltage probe across
      the 6 um OpenToGround gap at the far end of R3, analytically terminated in 50 kOhm
      (Gamma_2 = 0.998002) to recover exact unloaded 1-port reflection S11 (phase & group delay)
      and normalized open-end transfer |V_open / V_inc| (5.5-6.5 GHz passband & 2-4 GHz stopband).

Solver settings:
    - HFSS DrivenModal: max_passes = 3, basis_order = 1, freq_ghz = 6.0 GHz,
      max_delta_s = 0.02, pct_refinement = 20, explicit 10 um mesh on all 10 signal objects,
      Fast sweep 2.0 - 7.5 GHz (2201 pts, 2.5 MHz step).

Usage (from repository root, qmetal Python):
    python -X utf8 -I -B Design10_CPW_BPF_5p5_6p5GHz/simulation/resonator/hfss_s21.py
"""

import argparse
import io
import json
import math
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general')]
import project_paths  # noqa: E402

project_paths.configure(__file__)
DESIGN = project_paths.design_dir(__file__)
sys.path[:0] = [
    str(DESIGN / 'layout'),
    *project_paths.pattern_design_paths(),
]

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from qiskit_metal import Dict  # noqa: E402

import build as B  # noqa: E402
from run_record import _long_path, create_run, finish_run  # noqa: E402

MESH_LENGTH = '10um'
DEFAULT_PASSES = 3
STRAP_DIAMETER_MM = 0.010
Z_PROBE_OHM = 50000.0
PORTS = [('feed_in', 'start', 50), ('R3', 'end', int(Z_PROBE_OHM))]
COMPS = ['feed_in', 'T1', 'R1', 'T2', 'R2', 'T3', 'R3', 'S1', 'S2', 'S3']


def _write(path, data):
    lp = _long_path(path)
    _long_path(Path(path).parent).mkdir(parents=True, exist_ok=True)
    lp.write_text(json.dumps(data, indent=2, default=float) + '\n', encoding='utf-8')


def _parse_complex_col(series):
    if np.iscomplexobj(series):
        return series.to_numpy(complex)
    return series.map(lambda v: complex(str(v).replace('i', 'j').strip())).to_numpy(complex)


def _terminate_open_probe(s_df, z_probe_ohm=Z_PROBE_OHM, z0_ohm=50.0):
    s11 = _parse_complex_col(s_df['S11'])
    s21 = _parse_complex_col(s_df['S21'])
    s22 = _parse_complex_col(s_df['S22'])
    g2 = (z_probe_ohm - z0_ohm) / (z_probe_ohm + z0_ohm)
    denom = 1.0 - g2 * s22
    s11_open = s11 + (g2 * (s21 ** 2)) / denom
    s21_probe = (math.sqrt(1.0 - g2 ** 2) * s21) / denom
    return s11_open, s21_probe


def build_verified_straps(design, geo):
    from shapely.geometry import Point
    from shapely.ops import unary_union

    cpw = geo['cpw']
    w_mm = cpw['trace_width_um'] / 1000.0
    g_mm = cpw['trace_gap_um'] / 1000.0
    span_mm = (w_mm + 2.0 * g_mm) + 0.024

    y_main = float(geo['y_main_mm'])
    x_nodes = [float(x) for x in geo['stubs']['x_mm']]

    straps = []

    def add(name, x_mm, y_mm, angle_deg):
        straps.append(dict(name=name, x_mm=float(x_mm), y_mm=float(y_mm),
                           span_mm=float(span_mm), angle_deg=float(angle_deg)))

    # 1. Input feed (LP_IN.tie=-3.22 -> T1.prime_start=-2.30)
    add('strap_feed_0', -2.90, y_main, 90.0)
    add('strap_feed_1', -2.55, y_main, 90.0)

    # 2. All 3 arms around T1, T2, T3
    for idx, xk in enumerate(x_nodes, start=1):
        add(f'strap_T{idx}_L', xk - 0.13, y_main, 90.0)
        add(f'strap_T{idx}_R', xk + 0.13, y_main, 90.0)
        add(f'strap_T{idx}_S', xk, y_main + 0.18, 0.0)

    # 3. Straight shorted stubs S1, S2, S3
    for rep in design.metadata['stubs']:
        xk = float(rep['x_mm'])
        y_end = float(rep['y_short_mm'])
        span_y = y_end - y_main
        if span_y > 1.5:
            for k, frac in enumerate((0.35, 0.62, 0.85)):
                add(f"strap_{rep['name']}_{k}", xk, y_main + frac * span_y, 0.0)
        elif span_y > 0.6:
            add(f"strap_{rep['name']}_mid", xk, y_main + 0.65 * span_y, 0.0)

    # 4. Resonator U-loops R1, R2, R3
    for rep in design.metadata['resonators']:
        x_l = float(rep['x_leg_l_mm'])
        x_r = float(rep['x_leg_r_mm'])
        x_m = 0.5 * (x_l + x_r)
        y_bot = float(rep['y_loop_mm'])
        n_v = 4 if (y_main - y_bot) > 2.2 else 2
        for k, yv in enumerate(np.linspace(y_main - 0.35, y_bot + 0.35, n_v)):
            add(f"strap_{rep['name']}_L{k}", x_l, float(yv), 0.0)
            add(f"strap_{rep['name']}_R{k}", x_r, float(yv), 0.0)
        if (x_r - x_l) > 0.9:
            add(f"strap_{rep['name']}_B0", x_m - 0.28, y_bot, 90.0)
            add(f"strap_{rep['name']}_B1", x_m + 0.28, y_bot, 90.0)
        else:
            add(f"strap_{rep['name']}_B0", x_m, y_bot, 90.0)

    sig_geoms = []
    for _, row in design.qgeometry.tables['poly'].iterrows():
        if not bool(row['subtract']):
            sig_geoms.append(row['geometry'])
    for _, row in design.qgeometry.tables['path'].iterrows():
        if not bool(row['subtract']):
            fline = B._fillet_linestring(row['geometry'], float(row['fillet']) if row['fillet'] else 0.0)
            sig_geoms.append(fline.buffer(float(row['width']) / 2.0, cap_style=2))
    all_signal = unary_union(sig_geoms)

    foot_r = STRAP_DIAMETER_MM / 2.0
    worst_um = 1e9
    for st in straps:
        a = math.radians(st['angle_deg'])
        dx, dy = 0.5 * st['span_mm'] * math.cos(a), 0.5 * st['span_mm'] * math.sin(a)
        for sgn in (-1.0, +1.0):
            fx, fy = st['x_mm'] + sgn * dx, st['y_mm'] + sgn * dy
            d_mm = float(Point(fx, fy).buffer(foot_r).distance(all_signal))
            worst_um = min(worst_um, d_mm * 1e3)
            if d_mm < 0.008:
                raise RuntimeError(f"Wirebond {st['name']} too close: {d_mm * 1e3:.2f} um")
    print(f'Verified {len(straps)} PEC ground wirebonds (min foot-to-signal = {worst_um:.2f} um >= 8 um)', flush=True)
    return straps


def draw_straps(renderer, straps):
    from pyEPR.ansys import parse_units

    names = []
    for st in straps:
        a = np.radians(st['angle_deg'])
        names.append(renderer.modeler.draw_wirebond(
            pos=parse_units(np.array([st['x_mm'], st['y_mm']])),
            ori=np.array([np.cos(a), np.sin(a)]),
            width=parse_units(f"{st['span_mm']}mm"),
            height=parse_units(f"{st['span_mm']}mm"),
            z=0,
            wire_diameter=f'{STRAP_DIAMETER_MM}mm',
            NumSides=6,
            name='g_wb',
            material='pec',
            solve_inside=False,
        ))
    return names


def mesh_all_metal(renderer, design, mesh_length):
    id2n = {c.id: n for n, c in design.components.items()}
    objs = []
    for kind in ('path', 'poly'):
        for _, r in design.qgeometry.tables[kind].iterrows():
            comp = id2n[int(r['component'])]
            if not bool(r['subtract']) and comp in COMPS:
                objs.append(f"{r['name']}_{comp}")
    renderer.modeler.mesh_length('metal_mesh', objs, MaxLength=mesh_length)
    print(f'Applied mesh_length={mesh_length} on {len(objs)} signal metal objects', flush=True)
    return objs


def analyze_and_plot(freqs_ghz, s11_open, s21_probe, conv_df, mesh_length, out_paths):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    s21_db = 20.0 * np.log10(np.maximum(np.abs(s21_probe), 1e-12))
    s21_norm = s21_db - float(np.max(s21_db))
    s11_db = 20.0 * np.log10(np.maximum(np.abs(s11_open), 1e-12))
    ph_rad = np.unwrap(np.angle(s11_open))
    tau_ns = np.abs(np.gradient(ph_rad) / (2.0 * np.pi * np.gradient(freqs_ghz) * 1e9)) * 1e9

    pb = (freqs_ghz >= 5.5) & (freqs_ghz <= 6.5)
    sb = (freqs_ghz >= 2.0) & (freqs_ghz <= 4.0)
    idx3 = np.where(s21_norm >= -3.0)[0]
    f3_lo, f3_hi = float(freqs_ghz[idx3[0]]), float(freqs_ghz[idx3[-1]])

    def at_f(ft):
        return float(s21_norm[np.argmin(np.abs(freqs_ghz - ft))])

    last_ds = float(conv_df['Max Mag. Delta S'].dropna().iloc[-1])
    tets = int(conv_df['Solved Elements'].dropna().iloc[-1])

    metrics = dict(
        band_3dB_GHz=[round(f3_lo, 4), round(f3_hi, 4)],
        f0_center_3dB_GHz=round(0.5 * (f3_lo + f3_hi), 4),
        pb_5p5_6p5_min_norm_dB=round(float(np.min(s21_norm[pb])), 3),
        sb_2_4_max_norm_dB=round(float(np.max(s21_norm[sb])), 2),
        sb_2_4_meets_40dB=bool(np.max(s21_norm[sb]) <= -40.0),
        norm_at_2GHz_dB=round(at_f(2.0), 2),
        norm_at_3GHz_dB=round(at_f(3.0), 2),
        norm_at_4GHz_dB=round(at_f(4.0), 2),
        norm_at_5GHz_dB=round(at_f(5.0), 2),
        norm_at_6GHz_dB=round(at_f(6.0), 2),
        norm_at_7GHz_dB=round(at_f(7.0), 2),
        tau_g_5p5_6p5_median_ns=round(float(np.median(tau_ns[pb])), 3),
        tau_g_3GHz_ns=round(float(tau_ns[np.argmin(np.abs(freqs_ghz - 3.0))]), 4),
        passes=int(len(conv_df)),
        solved_elements=tets,
        last_delta_s=round(last_ds, 6),
        mesh_length=mesh_length,
    )

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(9.2, 9.2), sharex=True)

    ax0.axvspan(2.0, 4.0, color='#fadbd8', alpha=0.45, label='Stopband target: 2–4 GHz (≤ -40 dB)')
    ax0.axvspan(5.5, 6.5, color='#d6eaf8', alpha=0.50, label='Passband target: 5.5–6.5 GHz')
    ax0.plot(freqs_ghz, s21_norm, color='#1b4f72', lw=2.0, label='HFSS normalized |S21_probe| (50 kΩ open-end probe)')
    ax0.axhline(-40.0, color='#922b21', lw=0.9, ls=':')
    ax0.axhline(-3.0, color='#555555', lw=0.8, ls=':')
    ax0.set_ylim(-75, 2)
    ax0.set_ylabel('Normalized |S21_probe| (dB)')
    ax0.set_title(
        f"{DESIGN.name} · 1-Port 3-Stage Reflection BPF · HFSS Driven Only "
        f"({metrics['passes']} passes, mesh={mesh_length}, ΔS={last_ds:.4f}, {tets:,} tets)\n"
        f"-3 dB: {f3_lo:.2f}–{f3_hi:.2f} GHz | 5.5–6.5 GHz min={metrics['pb_5p5_6p5_min_norm_dB']:.2f} dB | "
        f"2–4 GHz max={metrics['sb_2_4_max_norm_dB']:.1f} dB",
        fontsize=9.5,
        weight='bold',
    )
    ax0.grid(True, alpha=0.3)
    ax0.legend(loc='lower right', fontsize=8.3, frameon=True)

    ax1.axvspan(2.0, 4.0, color='#fadbd8', alpha=0.35)
    ax1.axvspan(5.5, 6.5, color='#d6eaf8', alpha=0.45)
    ax1.plot(freqs_ghz, s11_db, color='#c0392b', lw=1.8, label='HFSS |S11| (dB, left)')
    ax1.set_ylim(-3.0, 0.5)
    ax1.set_ylabel('HFSS |S11| (dB)', color='#c0392b')
    ax1.tick_params(axis='y', labelcolor='#c0392b')
    ax1.grid(True, alpha=0.3)
    ax1r = ax1.twinx()
    ax1r.plot(freqs_ghz, np.degrees(ph_rad), color='#7d3c98', lw=1.5, ls='--', label='HFSS ∠S11 unwrapped (deg, right)')
    ax1r.set_ylabel('∠S11 (deg)', color='#7d3c98')
    ax1r.tick_params(axis='y', labelcolor='#7d3c98')
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax1r.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc='lower left', fontsize=8.3, frameon=True)
    ax1.set_title('1-Port Reflection S11: Magnitude & Unwrapped Phase', fontsize=9.5, weight='bold')

    ax2.axvspan(2.0, 4.0, color='#fadbd8', alpha=0.35)
    ax2.axvspan(5.5, 6.5, color='#d6eaf8', alpha=0.45)
    ax2.plot(freqs_ghz, tau_ns, color='#1e8449', lw=1.8, label='HFSS 1-port reflection group delay τ_g(S11)')
    ax2.set_yscale('log')
    ax2.set_xlim(freqs_ghz.min(), freqs_ghz.max())
    ax2.set_xlabel('Frequency (GHz)')
    ax2.set_ylabel('Group Delay τ_g(S11) (ns)')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend(loc='upper left', fontsize=8.3, frameon=True)

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=220, facecolor='white')
    plt.close(fig)
    png_bytes = buf.getvalue()
    for p in out_paths:
        _long_path(Path(p).parent).mkdir(parents=True, exist_ok=True)
        _long_path(p).write_bytes(png_bytes)
    return metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--passes', type=int, default=DEFAULT_PASSES)
    ap.add_argument('--mesh', type=str, default=MESH_LENGTH)
    ap.add_argument('--fmin', type=float, default=2.0)
    ap.add_argument('--fmax', type=float, default=7.5)
    ap.add_argument('--count', type=int, default=2201)
    args = ap.parse_args()

    geo = B.load_config()
    design = B.build_design(wire_bonds=False)
    straps = build_verified_straps(design, geo)

    params = dict(
        solver='hfss_driven',
        max_passes=args.passes,
        max_delta_s=0.02,
        mesh_length=args.mesh,
        basis_order=1,
        freq_ghz=6.0,
        sweep=dict(start_ghz=args.fmin, stop_ghz=args.fmax, count=args.count, type='Fast'),
        components=COMPS,
        ports=PORTS,
        strap_count=len(straps),
        geometry_revision=geo['geometry_revision'],
    )

    folder = create_run(
        DESIGN.name,
        'resonator',
        'hfss_driven_1port_3stage',
        entrypoint=__file__,
        inputs=[DESIGN / 'parameters' / 'geometry.json'],
        parameters=params,
    )
    print('RUN_FOLDER=' + str(folder), flush=True)

    try:
        from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
        import hfss_dirs
        import save_aedt

        hfss_dirs.use_d_drive()
        short_aedt = Path('D:/HFSS_temp/d10_1p.aedt')
        Path(str(short_aedt) + '.lock').unlink(missing_ok=True)
        save_aedt.use_project(str(short_aedt))

        rec = json.loads(_long_path(folder / 'run.json').read_text(encoding='utf-8'))
        rec['environment']['solver_version'] = save_aedt._desktop().get_version()
        rec['geometry_revision'] = geo['geometry_revision']
        _write(folder / 'run.json', rec)

        mesh_info = {}

        class DrivenSim(ScatteringImpedanceSim):
            def _analyze(self):
                mesh_info['straps'] = draw_straps(self.renderer, straps)
                mesh_info['mesh_objects'] = mesh_all_metal(self.renderer, design, args.mesh)
                super()._analyze()

        sim = DrivenSim(design, 'hfss')
        sim.setup.update(
            freq_ghz=6.0,
            max_delta_s=0.02,
            max_passes=args.passes,
            min_passes=1,
            min_converged=1,
            pct_refinement=20,
            basis_order=1,
            vars=Dict(),
        )
        sim.setup.sweep_setup.update(
            name='Fast_2_7p5GHz',
            start_ghz=args.fmin,
            stop_ghz=args.fmax,
            count=args.count,
            type='Fast',
            save_fields=False,
        )

        t0 = time.time()
        sim.run(name=f'D10_1P_{folder.name[-8:]}', components=COMPS, open_terminations=[], port_list=PORTS)

        _, _, s_df = sim.renderer.get_params(['S11', 'S21', 'S22'])
        s_df.to_csv(_long_path(folder / 'raw' / 'S_params_raw50.csv'))
        conv = sim.renderer.pinfo.setup.get_convergence()[0]
        conv.to_csv(_long_path(folder / 'raw' / 'hfss_convergence.csv'))
        sim.renderer.pinfo.project._project.Save()

        freqs = s_df.index.to_numpy(float)
        s11_open, s21_probe = _terminate_open_probe(s_df)

        s21_db = 20.0 * np.log10(np.maximum(np.abs(s21_probe), 1e-12))
        s21_norm = s21_db - float(np.max(s21_db))
        s11_db = 20.0 * np.log10(np.maximum(np.abs(s11_open), 1e-12))
        out_tbl = pd.DataFrame(dict(
            freq_GHz=freqs,
            S21_probe_norm_dB=np.round(s21_norm, 4),
            S21_probe_dB=np.round(s21_db, 4),
            S11_dB=np.round(s11_db, 4),
            S11_deg=np.round(np.degrees(np.angle(s11_open)), 3),
        ))
        out_tbl.to_csv(_long_path(folder / 'summaries' / 'hfss_1port_2_7p5GHz.csv'), index=False)

        run_fig = folder / 'figures' / 'hfss_1port_s21_s11.png'
        notes_fig = DESIGN / 'notes' / 'fig' / 'hfss_1port_s21_s11.png'
        metrics = analyze_and_plot(freqs, s11_open, s21_probe, conv, args.mesh, [run_fig, notes_fig])
        metrics['strap_count'] = len(mesh_info.get('straps', []))
        metrics['minutes'] = round((time.time() - t0) / 60.0, 2)
        metrics['aedt_design'] = sim.renderer.pinfo.design.name

        summary = dict(metrics=metrics, figure=str(notes_fig.relative_to(ROOT).as_posix()))
        print(json.dumps(metrics, indent=2), flush=True)
        finish_run(folder, status='completed', summary=summary)
    except Exception:
        _long_path(folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
