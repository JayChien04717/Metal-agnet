# -*- coding: utf-8 -*-
"""Design9_CPW_Bandpass_Filter HFSS Driven S21/S11 verification for the
4-8 GHz 4-stub quarter-wave short-circuited CPW bandpass filter (10/6 um CPW).

Solver configuration (Design9 user override recorded in run.json & README.md):
    - HFSS DrivenModal: max_passes = 3 (user requested 3 passes after explicit mesh),
      first-order basis, freq_ghz = 6.0 GHz, max_delta_s = 0.02, pct_refinement = 20.
    - Explicit length-based mesh: MaxLength = 10 um on all 13 signal metal objects
      (feed_in, T1..T4, L12, L23, L34, feed_out, S1..S4).
    - Broadband Interpolating sweep: 0.5 - 13.5 GHz (651 pts, 20 MHz step).
    - Ground wire bonds: geometry-verified explicit PEC straps across every T-junction
      arm, stub leg/top, connecting U-loop leg/bottom, and input/output feeds
      (all feet verified >= 8 um clear of signal metal).
    - Ports:
        Port 1: ('feed_in', 'start', 50)
        Port 2: ('feed_out', 'end', 50)

Usage (from repository root, qmetal Python):
    python -X utf8 -I -B Design9_CPW_Bandpass_Filter/simulation/resonator/hfss_s21.py
    python -X utf8 -I -B Design9_CPW_Bandpass_Filter/simulation/resonator/hfss_s21.py --reanalyze <run_dir>
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
PORTS = [('feed_in', 'start', 50), ('feed_out', 'end', 50)]
COMPS = [
    'feed_in', 'T1', 'L12', 'T2', 'L23', 'T3', 'L34', 'T4', 'feed_out',
    'S1', 'S2', 'S3', 'S4',
]


def _write(path, data):
    lp = _long_path(path)
    _long_path(Path(path).parent).mkdir(parents=True, exist_ok=True)
    lp.write_text(json.dumps(data, indent=2, default=float) + '\n', encoding='utf-8')


def _parse_complex_col(series):
    if np.iscomplexobj(series):
        return series.to_numpy(complex)
    return series.map(lambda v: complex(str(v).replace('i', 'j').strip())).to_numpy(complex)


def build_verified_straps(design, geo):
    """Build explicit PEC ground wirebond placements across every T-junction, stub, loop,
    and feed line, mathematically verifying every landing foot is >= 8 um clear of signal metal.
    """
    from shapely.geometry import Point
    from shapely.ops import unary_union

    cpw = geo['cpw']
    w_mm = cpw['trace_width_um'] / 1000.0
    g_mm = cpw['trace_gap_um'] / 1000.0
    single_span_mm = (w_mm + 2.0 * g_mm) + 0.024  # 46 um -> feet 12 um outside trench

    x_nodes = [float(x) for x in geo['stubs']['x_mm']]
    lead_off = geo['connecting_lines']['lead_offset_x_um'] / 1000.0

    straps = []

    def add_strap(name, x_mm, y_mm, angle_deg):
        straps.append(dict(
            name=name,
            x_mm=float(x_mm),
            y_mm=float(y_mm),
            span_mm=float(single_span_mm),
            angle_deg=float(angle_deg),
        ))

    # 1. Input and output 50-ohm feeds (horizontal CPW -> vertical strap angle=90)
    add_strap('strap_feed_in_0', -3.85, 0.0, 90.0)
    add_strap('strap_feed_in_1', -3.45, 0.0, 90.0)
    add_strap('strap_feed_out_0', +3.45, 0.0, 90.0)
    add_strap('strap_feed_out_1', +3.85, 0.0, 90.0)

    # 2. All 3 arms around each T-junction T1..T4 (critical for slotline suppression)
    for idx, xk in enumerate(x_nodes, start=1):
        add_strap(f'strap_T{idx}_left', xk - 0.16, 0.0, 90.0)
        add_strap(f'strap_T{idx}_right', xk + 0.16, 0.0, 90.0)
        add_strap(f'strap_T{idx}_stub_root', xk, +0.18, 0.0)

    # 3. Along all 4 folded short-circuited stubs S1..S4
    for rep in design.metadata['stubs']:
        name = rep['name']
        x_tap = float(rep['x_tap_mm'])
        x_sh = float(rep['x_short_mm'])
        y_top = float(rep['y_top_mm'])
        # Left vertical leg (connected to TapJunction)
        for k, yv in enumerate((0.80, 1.45, 2.10)):
            add_strap(f'strap_{name}_legL_{k}', x_tap, yv, 0.0)
        # Top horizontal bar midpoint
        add_strap(f'strap_{name}_top', 0.5 * (x_tap + x_sh), y_top, 90.0)
        # Right vertical leg (ending in ShortToGround at y = 0.35 mm)
        for k, yv in enumerate((0.55, 1.25, 1.95)):
            add_strap(f'strap_{name}_legR_{k}', x_sh, yv, 0.0)

    # 4. Along all 3 connecting quarter-wave U-loops L12, L23, L34
    for idx, rep in enumerate(design.metadata['connecting_lines']):
        name = rep['name']
        x_l = x_nodes[idx] + lead_off
        x_r = x_nodes[idx + 1] - lead_off
        x_m = 0.5 * (x_l + x_r)
        y_bot = float(rep['y_loop_mm'])
        for k, yv in enumerate((-0.45, -1.10)):
            add_strap(f'strap_{name}_legL_{k}', x_l, yv, 0.0)
            add_strap(f'strap_{name}_legR_{k}', x_r, yv, 0.0)
        add_strap(f'strap_{name}_bot_0', x_m - 0.30, y_bot, 90.0)
        add_strap(f'strap_{name}_bot_1', x_m + 0.30, y_bot, 90.0)

    # Mathematically verify every landing foot is >= 8 um clear of all signal metal
    sig_geoms = []
    for _, row in design.qgeometry.tables['poly'].iterrows():
        if not bool(row['subtract']):
            sig_geoms.append(row['geometry'])
    for _, row in design.qgeometry.tables['path'].iterrows():
        if not bool(row['subtract']):
            fline = B._fillet_linestring(row['geometry'], float(row['fillet']) if row['fillet'] else 0.0)
            sig_geoms.append(fline.buffer(float(row['width']) / 2.0, cap_style=2))
    all_signal = unary_union(sig_geoms)

    foot_radius_mm = STRAP_DIAMETER_MM / 2.0
    min_clear_mm = 0.008
    worst_um = 1e9
    for st in straps:
        a = math.radians(st['angle_deg'])
        dx = 0.5 * st['span_mm'] * math.cos(a)
        dy = 0.5 * st['span_mm'] * math.sin(a)
        for sgn in (-1.0, +1.0):
            fx = st['x_mm'] + sgn * dx
            fy = st['y_mm'] + sgn * dy
            dist_mm = float(Point(fx, fy).buffer(foot_radius_mm).distance(all_signal))
            worst_um = min(worst_um, dist_mm * 1e3)
            if dist_mm < min_clear_mm:
                raise RuntimeError(
                    f"Wirebond foot of {st['name']} at ({fx:.3f}, {fy:.3f}) mm is only "
                    f"{dist_mm * 1e3:.2f} um from signal metal (< {min_clear_mm * 1e3:.1f} um)"
                )
    print(f'Verified {len(straps)} PEC ground wirebonds (min foot-to-signal clearance = {worst_um:.2f} um >= 8 um)',
          flush=True)
    return straps


def plot_wirebond_map(design, straps, out_png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    metal, _, _, chip_win = B.design_to_metal_um(design)
    w_c, h_c = chip_win[1] - chip_win[0], chip_win[3] - chip_win[2]
    fig, ax = plt.subplots(figsize=(9.0, 9.0 * h_c / w_c))
    B._draw_metal(ax, metal, chip_win)
    for st in straps:
        a = math.radians(st['angle_deg'])
        dx = 0.5 * st['span_mm'] * 1e3 * math.cos(a)
        dy = 0.5 * st['span_mm'] * 1e3 * math.sin(a)
        x, y = st['x_mm'] * 1e3, st['y_mm'] * 1e3
        ax.plot([x - dx, x + dx], [y - dy, y + dy], color='#d7322b', lw=1.8, zorder=6)
        ax.plot(x, y, 'o', color='#d7322b', ms=2.5, zorder=6)
    B._scale_bar(ax, chip_win)
    fig.suptitle(f'{DESIGN.name} · Verified PEC ground wirebonds ({len(straps)} straps)',
                 x=0.02, ha='left', fontsize=10.5, weight='bold')
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=180, facecolor='white')
    plt.close(fig)
    _long_path(Path(out_png).parent).mkdir(parents=True, exist_ok=True)
    _long_path(out_png).write_bytes(buf.getvalue())


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


def mesh_all_metal(renderer, design, comps, mesh_length):
    id2n = {c.id: n for n, c in design.components.items()}
    objs = []
    for kind in ('path', 'poly'):
        for _, r in design.qgeometry.tables[kind].iterrows():
            comp = id2n[int(r['component'])]
            if not bool(r['subtract']) and comp in comps:
                objs.append(f"{r['name']}_{comp}")
    renderer.modeler.mesh_length('metal_mesh', objs, MaxLength=mesh_length)
    print(f'Applied explicit mesh_length={mesh_length} on {len(objs)} signal metal objects: {objs}', flush=True)
    return objs


def analytic_abcd_s(freqs_ghz, geo):
    z0 = float(geo['cpw']['assumed_z0_ohm'])
    f0 = float(geo['filter_target']['f0_GHz'])
    stub_scales = [float(v) for v in geo['stubs']['electrical_scale']]
    conn_scales = [float(v) for v in geo['connecting_lines']['electrical_scale']]
    theta = (np.pi / 2.0) * (freqs_ghz / f0)

    def abcd_tline(z, th):
        c, s = np.cos(th), np.sin(th)
        return np.array([[c, 1j * z * s], [1j * s / z, c]])

    def abcd_sc(z, th):
        t = np.tan(th)
        y = np.where(np.abs(t) < 1e-12, -1j * 1e12, -1j / (z * t))
        return np.array([[np.ones_like(th), np.zeros_like(th)], [y, np.ones_like(th)]])

    def mm(a, b):
        return np.array([
            [a[0, 0] * b[0, 0] + a[0, 1] * b[1, 0], a[0, 0] * b[0, 1] + a[0, 1] * b[1, 1]],
            [a[1, 0] * b[0, 0] + a[1, 1] * b[1, 0], a[1, 0] * b[0, 1] + a[1, 1] * b[1, 1]],
        ])

    m = abcd_sc(z0, theta * stub_scales[0])
    for cs, ss in zip(conn_scales, stub_scales[1:]):
        m = mm(m, abcd_tline(z0, theta * cs))
        m = mm(m, abcd_sc(z0, theta * ss))

    a, b, c, d = m[0, 0], m[0, 1], m[1, 0], m[1, 1]
    den = a + b / z0 + c * z0 + d
    s21 = 2.0 / den
    s11 = (a + b / z0 - c * z0 - d) / den
    return s11, s21


def analyze_s_params(freqs_ghz, s11, s21, s22, conv_df, geo, mesh_length):
    s21_db = 20.0 * np.log10(np.maximum(np.abs(s21), 1e-12))
    s11_db = 20.0 * np.log10(np.maximum(np.abs(s11), 1e-12))
    s22_db = 20.0 * np.log10(np.maximum(np.abs(s22), 1e-12))

    # Main passband region around 3-9 GHz
    main_mask = (freqs_ghz >= 2.0) & (freqs_ghz <= 10.5)
    f_m = freqs_ghz[main_mask]
    s21_m = s21_db[main_mask]

    # Peak transmission and -3 dB / -1 dB contiguous band around peak
    pk_val = float(np.max(s21_m))
    idx_3 = np.where(s21_m >= pk_val - 3.0)[0]
    f3_lo, f3_hi = float(f_m[idx_3[0]]), float(f_m[idx_3[-1]])
    idx_1 = np.where(s21_m >= pk_val - 1.0)[0]
    f1_lo, f1_hi = float(f_m[idx_1[0]]), float(f_m[idx_1[-1]])
    f0_hfss = 0.5 * (f3_lo + f3_hi)

    pb_mask = (freqs_ghz >= 4.0) & (freqs_ghz <= 8.0)
    def at_f(f_target):
        return float(s21_db[np.argmin(np.abs(freqs_ghz - f_target))])

    # Upper transmission zero near 2*f0
    hi_mask = (freqs_ghz >= 9.5) & (freqs_ghz <= 13.5)
    f_zero_hi = float(freqs_ghz[hi_mask][np.argmin(s21_db[hi_mask])])
    s21_zero_hi = float(np.min(s21_db[hi_mask]))

    # Effective phase velocity calibrated from f0_hfss
    v_ph_assumed = float(geo['cpw']['assumed_v_ph_m_s'])
    f0_target = float(geo['filter_target']['f0_GHz'])
    v_ph_eff = v_ph_assumed * (f0_hfss / f0_target)

    last_delta_s = float(conv_df['Max Mag. Delta S'].dropna().iloc[-1])
    solved_elements = int(conv_df['Solved Elements'].dropna().iloc[-1])

    return dict(
        f0_target_GHz=f0_target,
        f0_hfss_center_3dB_GHz=round(f0_hfss, 4),
        band_3dB_GHz=[round(f3_lo, 4), round(f3_hi, 4)],
        bw_3dB_GHz=round(f3_hi - f3_lo, 4),
        band_1dB_GHz=[round(f1_lo, 4), round(f1_hi, 4)],
        peak_s21_dB=round(pk_val, 3),
        pb_4_8_min_s21_dB=round(float(np.min(s21_db[pb_mask])), 3),
        pb_4_8_max_s11_dB=round(float(np.max(s11_db[pb_mask])), 3),
        pb_4_8_max_s22_dB=round(float(np.max(s22_db[pb_mask])), 3),
        s21_at_2GHz_dB=round(at_f(2.0), 2),
        s21_at_3GHz_dB=round(at_f(3.0), 2),
        s21_at_4GHz_dB=round(at_f(4.0), 2),
        s21_at_6GHz_dB=round(at_f(6.0), 2),
        s21_at_8GHz_dB=round(at_f(8.0), 2),
        s21_at_9GHz_dB=round(at_f(9.0), 2),
        s21_at_10GHz_dB=round(at_f(10.0), 2),
        upper_zero_GHz=round(f_zero_hi, 4),
        upper_zero_s21_dB=round(s21_zero_hi, 2),
        effective_v_ph_m_s=round(v_ph_eff, 1),
        passes=int(len(conv_df)),
        solved_elements=solved_elements,
        last_delta_s=round(last_delta_s, 6),
        mesh_length=mesh_length,
    )


def plot_s_params(freqs_ghz, s11, s21, s22, metrics, geo, out_paths):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    s21_db = 20.0 * np.log10(np.maximum(np.abs(s21), 1e-12))
    s11_db = 20.0 * np.log10(np.maximum(np.abs(s11), 1e-12))
    s22_db = 20.0 * np.log10(np.maximum(np.abs(s22), 1e-12))
    s11_mod, s21_mod = analytic_abcd_s(freqs_ghz, geo)
    s21_mod_db = 20.0 * np.log10(np.maximum(np.abs(s21_mod), 1e-12))
    s11_mod_db = 20.0 * np.log10(np.maximum(np.abs(s11_mod), 1e-12))

    # Group delay of S21 in passband
    phase_s21 = np.unwrap(np.angle(s21))
    df_hz = np.gradient(freqs_ghz) * 1e9
    tau_ns = -np.gradient(phase_s21) / (2.0 * np.pi * df_hz) * 1e9

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9.2, 7.6), gridspec_kw={'height_ratios': [2.2, 1.0]})

    ax0.axvspan(4.0, 8.0, color='#d6eaf8', alpha=0.50, label='Target passband (4–8 GHz)')
    ax0.plot(freqs_ghz, s21_db, color='#1b4f72', lw=2.0, label='HFSS |S21|')
    ax0.plot(freqs_ghz, s11_db, color='#c0392b', lw=1.6, label='HFSS |S11|')
    ax0.plot(freqs_ghz, s22_db, color='#e67e22', lw=1.1, ls=':', label='HFSS |S22|')
    ax0.plot(freqs_ghz, s21_mod_db, color='#5dade2', lw=1.2, ls='--', label='Analytic ABCD |S21|')
    ax0.plot(freqs_ghz, s11_mod_db, color='#f1948a', lw=1.0, ls='--', label='Analytic ABCD |S11|')
    ax0.axhline(-3.0, color='#555555', lw=0.8, ls=':')
    ax0.axhline(-15.0, color='#888888', lw=0.8, ls=':')
    ax0.set_xlim(freqs_ghz.min(), freqs_ghz.max())
    ax0.set_ylim(-55, 2)
    ax0.set_ylabel('Magnitude (dB)')
    b3 = metrics['band_3dB_GHz']
    ax0.set_title(
        f"{DESIGN.name} · HFSS Driven S21/S11 ({metrics['passes']} passes, mesh={metrics['mesh_length']}, "
        f"ΔS={metrics['last_delta_s']:.4f}, {metrics['solved_elements']:,} tetras)\n"
        f"HFSS -3 dB: {b3[0]:.2f}–{b3[1]:.2f} GHz (f0={metrics['f0_hfss_center_3dB_GHz']:.2f} GHz) | "
        f"4–8 GHz: min |S21|={metrics['pb_4_8_min_s21_dB']:.2f} dB, max |S11|={metrics['pb_4_8_max_s11_dB']:.2f} dB",
        fontsize=9.8,
        weight='bold',
    )
    ax0.grid(True, alpha=0.3)
    ax0.legend(loc='lower center', ncol=3, fontsize=8.2, frameon=True)

    # Bottom panel: passband zoom (3–9 GHz) for S21, S11 and group delay
    ax1.axvspan(4.0, 8.0, color='#d6eaf8', alpha=0.50)
    ax1.plot(freqs_ghz, s21_db, color='#1b4f72', lw=1.8, label='HFSS |S21| (left)')
    ax1.set_xlim(3.0, 9.0)
    ax1.set_ylim(-6.0, 0.5)
    ax1.set_xlabel('Frequency (GHz)')
    ax1.set_ylabel('|S21| zoom (dB)', color='#1b4f72')
    ax1.tick_params(axis='y', labelcolor='#1b4f72')
    ax1.grid(True, alpha=0.3)

    ax1r = ax1.twinx()
    ax1r.plot(freqs_ghz, tau_ns, color='#1e8449', lw=1.4, ls='-.', label='Group delay τ_g (right)')
    ax1r.set_ylim(0, max(1.0, min(3.0, float(np.percentile(tau_ns[(freqs_ghz >= 4.0) & (freqs_ghz <= 8.0)], 98)) * 1.5)))
    ax1r.set_ylabel('Group delay τ_g (ns)', color='#1e8449')
    ax1r.tick_params(axis='y', labelcolor='#1e8449')
    ax1.set_title('Passband Zoom (3–9 GHz): Insertion Loss & Group Delay', fontsize=9.5, weight='bold')

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=200, facecolor='white')
    plt.close(fig)
    png_bytes = buf.getvalue()
    for p in out_paths:
        _long_path(Path(p).parent).mkdir(parents=True, exist_ok=True)
        _long_path(p).write_bytes(png_bytes)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--passes', type=int, default=DEFAULT_PASSES, help='HFSS driven adaptive passes (default 3)')
    ap.add_argument('--mesh', type=str, default=MESH_LENGTH, help='Explicit metal length mesh (default 10um)')
    ap.add_argument('--fmin', type=float, default=0.5, help='Sweep start GHz (default 0.5)')
    ap.add_argument('--fmax', type=float, default=13.5, help='Sweep stop GHz (default 13.5)')
    ap.add_argument('--count', type=int, default=651, help='Sweep point count (default 651 -> 20 MHz step)')
    ap.add_argument('--reanalyze', type=Path, default=None, help='Reanalyze an existing raw HFSS run folder')
    args = ap.parse_args()

    geo = B.load_config()
    if args.reanalyze is not None:
        raw_run = args.reanalyze.resolve()
        s_df = pd.read_csv(_long_path(raw_run / 'raw' / 'S_params.csv'), index_col=0)
        conv = pd.read_csv(_long_path(raw_run / 'raw' / 'hfss_convergence.csv'))
        freqs = s_df.index.to_numpy(float)
        s11 = _parse_complex_col(s_df['S11'])
        s21 = _parse_complex_col(s_df['S21'])
        s22 = _parse_complex_col(s_df['S22'])
        metrics = analyze_s_params(freqs, s11, s21, s22, conv, geo, args.mesh)
        notes_fig = DESIGN / 'notes' / 'fig' / 'hfss_s21.png'
        plot_s_params(freqs, s11, s21, s22, metrics, geo, [raw_run / 'figures' / 'hfss_s21.png', notes_fig])
        print(json.dumps(metrics, indent=2), flush=True)
        return

    design = B.build_design(wire_bonds=False)
    straps = build_verified_straps(design, geo)

    params = dict(
        solver='hfss_driven',
        max_passes=args.passes,
        solver_override_reason='User explicitly requested 3 adaptive passes after setting explicit length-based mesh.',
        max_delta_s=0.02,
        mesh_length=args.mesh,
        basis_order=1,
        freq_ghz=float(geo['filter_target']['f0_GHz']),
        sweep=dict(start_ghz=args.fmin, stop_ghz=args.fmax, count=args.count, type='Interpolating'),
        components=COMPS,
        ports=PORTS,
        strap_count=len(straps),
        straps=straps,
        geometry_revision=geo['geometry_revision'],
    )

    folder = create_run(
        DESIGN.name,
        'resonator',
        'hfss_driven_s21',
        entrypoint=__file__,
        inputs=[DESIGN / 'parameters' / 'geometry.json'],
        parameters=params,
    )
    print('RUN_FOLDER=' + str(folder), flush=True)
    plot_wirebond_map(design, straps, folder / 'figures' / 'wirebond_map.png')

    try:
        from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
        import hfss_dirs
        import save_aedt

        hfss_dirs.use_d_drive()
        short_aedt = Path('D:/HFSS_temp/d9bpf.aedt')
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
                mesh_info['mesh_objects'] = mesh_all_metal(self.renderer, design, COMPS, args.mesh)
                super()._analyze()

        sim = DrivenSim(design, 'hfss')
        f0 = float(geo['filter_target']['f0_GHz'])
        sim.setup.update(
            freq_ghz=f0,
            max_delta_s=0.02,
            max_passes=args.passes,
            min_passes=1,
            min_converged=1,
            pct_refinement=20,
            basis_order=1,
            vars=Dict(),
        )
        sim.setup.sweep_setup.update(
            name='Broadband',
            start_ghz=args.fmin,
            stop_ghz=args.fmax,
            count=args.count,
            type='Interpolating',
            save_fields=False,
        )

        t0 = time.time()
        try:
            sim.run(
                name=f'D9BPF_{folder.name[-8:]}',
                components=COMPS,
                open_terminations=[],
                port_list=PORTS,
            )
        except Exception:
            try:
                pinfo = sim.renderer.pinfo
                msgs = pinfo.desktop._desktop.GetMessages(pinfo.project_name, pinfo.design_name, 0)
                print('AEDT_MESSAGES:\n' + '\n'.join(msgs), flush=True)
            except Exception as e_msg:
                print('Could not read AEDT messages:', e_msg, flush=True)
            raise

        _, _, s_df = sim.renderer.get_params(['S11', 'S21', 'S22'])
        s_df.to_csv(_long_path(folder / 'raw' / 'S_params.csv'))
        conv = sim.renderer.pinfo.setup.get_convergence()[0]
        conv.to_csv(_long_path(folder / 'raw' / 'hfss_convergence.csv'))
        sim.renderer.pinfo.project._project.Save()

        freqs = s_df.index.to_numpy(float)
        s11 = _parse_complex_col(s_df['S11'])
        s21 = _parse_complex_col(s_df['S21'])
        s22 = _parse_complex_col(s_df['S22'])

        metrics = analyze_s_params(freqs, s11, s21, s22, conv, geo, args.mesh)
        metrics['mesh_object_count'] = len(mesh_info.get('mesh_objects', []))
        metrics['mesh_objects'] = mesh_info.get('mesh_objects', [])
        metrics['strap_count'] = len(mesh_info.get('straps', []))
        metrics['minutes'] = round((time.time() - t0) / 60.0, 2)
        metrics['aedt_design'] = sim.renderer.pinfo.design.name

        run_fig = folder / 'figures' / 'hfss_s21.png'
        notes_fig = DESIGN / 'notes' / 'fig' / 'hfss_s21.png'
        plot_s_params(freqs, s11, s21, s22, metrics, geo, [run_fig, notes_fig])

        summary = dict(
            metrics=metrics,
            figure=str(notes_fig.relative_to(ROOT).as_posix()),
            validation=(
                f"HFSS Driven {metrics['passes']}-pass run (user override: max_passes={args.passes}) with explicit "
                f"{args.mesh} length mesh on {metrics['mesh_object_count']} signal metal objects and "
                f"{metrics['strap_count']} geometry-verified PEC ground straps; "
                f"final Max Mag. Delta S = {metrics['last_delta_s']:.5f} (configured target 0.02, "
                f"met={metrics['last_delta_s'] <= 0.02})."
            ),
        )
        print(json.dumps(metrics, indent=2), flush=True)
        finish_run(folder, status='completed', summary=summary)
    except Exception:
        _long_path(folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
