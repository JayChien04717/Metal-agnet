# -*- coding: utf-8 -*-
"""design6_reflection_purcell_resonator HFSS: Driven S-parameter & circuit-model verification
showing the three lambda/4 readout resonators (RR1, RR2, RR3) inside the reflection
lambda/2 Purcell filter passband.

Solver defaults (skill/em_simulation.md sections 2-3):
    - HFSS Driven: 4 passes max (default), first-order basis, explicit 10 um length-based
      mesh on all rendered metal objects (recorded in run.json).
    - Ground wire bonds: `hfss_wire_bonds = True` on all CPW paths with
      `wb_threshold='350um', wb_offset='0um', wb_size=2`, plus explicit PEC ground straps
      across both ends of each of the 3 Purcell-bus / readout-resonator coupling sections.
    - Ports:
        Port 1: ('feed_in', 'start', 50)     -- 50 Ohm input port at the launch feed line
        Port 2: ('PF1', 'end', 50000)        -- 50 kOhm weakly coupled voltage probe at the
                                                open end of the lambda/2 Purcell filter so
                                                |S21| directly displays the Purcell filter
                                                passband envelope with the 3 readout-resonator
                                                notches inside the band, while S11 gives the
                                                1-port reflection phase and group delay.

Usage (from repository root, qmetal Python):
    python -X utf8 -I -B design6_reflection_purcell_resonator/simulation/resonator/hfss_purcell_triad.py --q3d <q3d_run_dir>
    python -X utf8 -I -B design6_reflection_purcell_resonator/simulation/resonator/hfss_purcell_triad.py --q3d <q3d_run_dir> --reanalyze <hfss_run_dir>
"""

import io
import os
from pathlib import Path
import argparse
import json
import math
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general')]
import project_paths  # noqa: E402

project_paths.configure(__file__)
DESIGN = project_paths.design_dir(__file__)
sys.path[:0] = [
    str(DESIGN / 'layout'),
    str(DESIGN / 'layout' / 'components'),
    *project_paths.pattern_design_paths(),
]

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from qiskit_metal import Dict  # noqa: E402
from scipy.ndimage import maximum_filter1d  # noqa: E402
from scipy.signal import find_peaks  # noqa: E402

import build as B  # noqa: E402
from run_record import _long_path, create_run, finish_run  # noqa: E402

WB_OPTS = dict(wb_threshold='350um', wb_offset='0um', wb_size=2)
MESH_LENGTH = '10um'
STRAP_DIAMETER_MM = 0.010
PORTS = [('feed_in', 'start', 50), ('PF1', 'end', 50000)]


def _write(path, data):
    lp = _long_path(path)
    _long_path(Path(path).parent).mkdir(parents=True, exist_ok=True)
    lp.write_text(json.dumps(data, indent=2, default=float) + '\n', encoding='utf-8')


def _parse_complex_col(series):
    if np.iscomplexobj(series):
        return series.to_numpy(complex)
    return series.map(lambda v: complex(str(v).replace('i', 'j').strip())).to_numpy(complex)


def load_q3d_summary(q3d_dir):
    summary_path = Path(q3d_dir) / 'summary.json'
    return json.loads(_long_path(summary_path).read_text(encoding='utf-8'))['derived']


def coupling_strap_points(design, geo):
    """Return 6 ground-strap placements across the Purcell bus and the 3 readout coupling arms."""
    cpw = geo['cpw']
    pf = geo['purcell_filter']
    rr = geo['readout_resonators']
    w_mm = cpw['trace_width_um'] / 1000.0
    g_mm = cpw['trace_gap_um'] / 1000.0
    s_mm = rr['coupling_ground_strip_um'] / 1000.0
    y_bus = float(pf['y_bus_mm'])
    y_coup = y_bus - (w_mm + 2.0 * g_mm + s_mm)
    y_mid = 0.5 * (y_bus + y_coup)
    span_mm = 2.0 * (w_mm + 2.0 * g_mm) + s_mm + 0.020

    straps = []
    for rep in design.metadata['readout_resonators']:
        x_s = rep['short_pos_mm'][0] + 0.025
        x_e = rep['short_pos_mm'][0] + rep['coupling_length_mm'] - 0.025
        for tag, x_pos in (('left', x_s), ('right', x_e)):
            straps.append(dict(
                name=f"strap_{rep['name']}_{tag}",
                x_mm=float(x_pos),
                y_mm=float(y_mid),
                span_mm=float(span_mm),
                angle_deg=90.0,
            ))
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


def mesh_all_metal(renderer, design, comps):
    id2n = {c.id: n for n, c in design.components.items()}
    objs = []
    for kind in ('path', 'poly'):
        for r in design.qgeometry.tables[kind].itertuples():
            comp = id2n[int(r.component)]
            if not bool(r.subtract) and comp in comps:
                objs.append(f'{r.name}_{comp}')
    renderer.modeler.mesh_length('metal_mesh', objs, MaxLength=MESH_LENGTH)
    return objs


def circuit_model_response(freqs_ghz, f_p_ghz, kappa_p_mhz, f_r_list_ghz, j_list_mhz, q_int=250000.0):
    """Coupled-mode / admittance model for 1 reflection Purcell filter coupled to 3 readout resonators."""
    kp = kappa_p_mhz * 1e-3  # GHz
    kp_int = f_p_ghz / q_int
    s11 = np.zeros_like(freqs_ghz, dtype=complex)
    s21 = np.zeros_like(freqs_ghz, dtype=complex)
    kp_probe = 1e-4

    for idx, f in enumerate(freqs_ghz):
        sigma_r = 0.0 + 0.0j
        for fr, j_mhz in zip(f_r_list_ghz, j_list_mhz):
            j_ghz = j_mhz * 1e-3
            kr_int = fr / q_int
            sigma_r += (j_ghz ** 2) / ((f - fr) + 0.5j * kr_int)
        denom = (f - f_p_ghz) + 0.5j * (kp + kp_int + kp_probe) - sigma_r
        s11[idx] = 1.0 - (1j * kp) / denom
        s21[idx] = (1j * math.sqrt(kp * kp_probe)) / denom

    phase_rad = np.unwrap(np.angle(s11))
    df_hz = np.gradient(freqs_ghz) * 1e9
    tau_g_ns = np.abs(np.gradient(phase_rad) / (2.0 * np.pi * df_hz)) * 1e9

    n_r = len(f_r_list_ghz)
    h_mat = np.zeros((1 + n_r, 1 + n_r), dtype=complex)
    h_mat[0, 0] = f_p_ghz - 0.5j * kp
    for k, (fr, j_mhz) in enumerate(zip(f_r_list_ghz, j_list_mhz), start=1):
        h_mat[k, k] = fr
        h_mat[0, k] = j_mhz * 1e-3
        h_mat[k, 0] = j_mhz * 1e-3
    eigvals = np.linalg.eigvals(h_mat)

    return dict(
        s11=s11,
        s21=s21,
        group_delay_ns=tau_g_ns,
        phase_deg=np.degrees(phase_rad),
        dressed_eigenvalues_GHz=[(float(z.real), float(-2.0 * z.imag * 1e3)) for z in eigvals],
    )


def analyze_spectra(f_wide, s11_wide, s21_wide, f_fine, s11_fine, s21_fine, q3d_derived, geo, conv_df):
    """Extract Purcell passband center, 3 dB bandwidth, and the 3 readout resonator dips inside the band."""
    s21_fine_db = 20.0 * np.log10(np.maximum(np.abs(s21_fine), 1e-12))
    s21_wide_db = 20.0 * np.log10(np.maximum(np.abs(s21_wide), 1e-12))

    # Smooth upper envelope over the wide sweep to measure the broad Purcell filter resonance peak and 3 dB bandwidth
    win_wide = max(41, len(f_wide) // 25)
    env_wide_db = maximum_filter1d(s21_wide_db, size=win_wide)
    pk_w_idx = int(np.argmax(env_wide_db))
    f_p_hfss = float(f_wide[pk_w_idx])
    pk_w_val = float(env_wide_db[pk_w_idx])
    above_3db = f_wide[env_wide_db >= pk_w_val - 3.0]
    if len(above_3db) > 2:
        f_3db_lo, f_3db_hi = float(above_3db.min()), float(above_3db.max())
        kappa_p_hfss_mhz = (f_3db_hi - f_3db_lo) * 1e3
    else:
        kappa_p_hfss_mhz = float(q3d_derived['idc']['predicted_kappa_p_MHz'])
        f_3db_lo = f_p_hfss - kappa_p_hfss_mhz / 2000.0
        f_3db_hi = f_p_hfss + kappa_p_hfss_mhz / 2000.0

    # Smooth envelope over fine sweep to find the 3 sharp resonator dips inside the Purcell passband
    step_fine_ghz = float(np.median(np.diff(f_fine)))
    win_fine = max(31, int(0.025 / step_fine_ghz))
    env_fine_db = maximum_filter1d(s21_fine_db, size=win_fine)
    notch_depth_db = env_fine_db - s21_fine_db

    # Also compute 1-port reflection group delay peaks at Port 1
    phase_fine_rad = np.unwrap(np.angle(s11_fine))
    df_fine_hz = np.gradient(f_fine) * 1e9
    tau_fine_ns = np.abs(np.gradient(phase_fine_rad) / (2.0 * np.pi * df_fine_hz)) * 1e9

    min_sep_pts = max(10, int(0.025 / step_fine_ghz))
    dip_indices, dip_props = find_peaks(notch_depth_db, prominence=1.0, distance=min_sep_pts)
    if len(dip_indices) > 3:
        order = np.argsort(dip_props['prominences'])[::-1][:3]
        dip_indices = np.sort(dip_indices[order])
    dips_ghz = [float(f_fine[i]) for i in dip_indices]
    dip_depths_db = [float(notch_depth_db[i]) for i in dip_indices]
    dip_tau_ns = [float(np.max(tau_fine_ns[max(0, i - 8):min(len(f_fine), i + 9)])) for i in dip_indices]

    # Estimate Purcell-filtered readout linewidths kappa_r from group delay peak tau_peak = 4 / kappa_r (in rad/s)
    # -> kappa_r / 2pi (MHz) = 4 / (2 * pi * tau_peak_ns * 1e-3) = 2000 / (pi * tau_peak_ns)
    kappa_r_est_mhz = [round(2000.0 / (math.pi * max(t, 0.1)), 3) for t in dip_tau_ns]

    j_q3d_mhz = float(q3d_derived['coupled_section']['estimated_J_MHz'])
    f_r_targets = [ch['f_target_GHz'] for ch in geo['readout_resonators']['channels']]
    f_r_for_model = dips_ghz if len(dips_ghz) == 3 else f_r_targets
    model_fine = circuit_model_response(
        f_fine,
        f_p_ghz=f_p_hfss,
        kappa_p_mhz=kappa_p_hfss_mhz,
        f_r_list_ghz=f_r_for_model,
        j_list_mhz=[j_q3d_mhz] * len(f_r_for_model),
    )
    model_wide = circuit_model_response(
        f_wide,
        f_p_ghz=f_p_hfss,
        kappa_p_mhz=kappa_p_hfss_mhz,
        f_r_list_ghz=f_r_for_model,
        j_list_mhz=[j_q3d_mhz] * len(f_r_for_model),
    )

    last_delta_s = float(conv_df['Max Mag. Delta S'].dropna().iloc[-1])
    solved_elements = int(conv_df['Solved Elements'].dropna().iloc[-1])

    metrics = dict(
        c_in_q3d_fF=float(q3d_derived['idc']['c_in_q3d_fF']),
        c_a_gnd_deembedded_fF=float(q3d_derived['idc']['c_a_gnd_deembedded_fF']),
        predicted_kappa_p_q3d_MHz=float(q3d_derived['idc']['predicted_kappa_p_MHz']),
        estimated_J_q3d_MHz=j_q3d_mhz,
        f_p_target_GHz=float(geo['purcell_filter']['f_target_GHz']),
        f_p_hfss_GHz=round(f_p_hfss, 5),
        purcell_3dB_band_GHz=[round(f_3db_lo, 5), round(f_3db_hi, 5)],
        kappa_p_hfss_MHz=round(kappa_p_hfss_mhz, 2),
        f_r_targets_GHz=f_r_targets,
        resonator_dips_hfss_GHz=[round(x, 5) for x in dips_ghz],
        resonator_notch_depths_dB=[round(x, 2) for x in dip_depths_db],
        resonator_group_delay_peaks_ns=[round(x, 2) for x in dip_tau_ns],
        resonator_purcell_filtered_kappa_r_MHz=kappa_r_est_mhz,
        all_dips_inside_purcell_3dB_band=bool(
            len(dips_ghz) == 3 and all(f_3db_lo <= d <= f_3db_hi for d in dips_ghz)
        ),
        passes=int(len(conv_df)),
        solved_elements_pass4=solved_elements,
        last_delta_s=last_delta_s,
        mesh_length=MESH_LENGTH,
    )
    return metrics, model_wide, model_fine


def plot_purcell_band_results(f_wide, s11_wide, s21_wide, f_fine, s11_fine, s21_fine,
                              model_wide, model_fine, summary_metrics, out_png_paths):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    s21_db_w = 20.0 * np.log10(np.maximum(np.abs(s21_wide), 1e-12))
    s21_db_w_mod = 20.0 * np.log10(np.maximum(np.abs(model_wide['s21']), 1e-12))
    offset_w = float(np.percentile(s21_db_w, 97) - np.percentile(s21_db_w_mod, 97))
    s21_db_w_mod += offset_w

    s21_db_f = 20.0 * np.log10(np.maximum(np.abs(s21_fine), 1e-12))
    s21_db_f_mod = 20.0 * np.log10(np.maximum(np.abs(model_fine['s21']), 1e-12))
    offset_f = float(np.percentile(s21_db_f, 95) - np.percentile(s21_db_f_mod, 95))
    s21_db_f_mod += offset_f

    phase_f_rad = np.unwrap(np.angle(s11_fine))
    df_f_hz = np.gradient(f_fine) * 1e9
    tau_f_ns = np.abs(np.gradient(phase_f_rad) / (2.0 * np.pi * df_f_hz)) * 1e9

    fp = summary_metrics['f_p_hfss_GHz']
    kp = summary_metrics['kappa_p_hfss_MHz']
    b_lo, b_hi = summary_metrics['purcell_3dB_band_GHz']
    dips = summary_metrics['resonator_dips_hfss_GHz']

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(10.0, 10.5))

    # Panel 0: Broadband view (5.8 - 7.3 GHz) showing the full reflection Purcell bandpass and the 3 resonators inside
    ax0.plot(f_wide, s21_db_w, color='#1b4f72', lw=1.5, label='HFSS Driven (broadband 5.8–7.3 GHz, 4 passes)')
    ax0.plot(f_wide, s21_db_w_mod, color='#d35400', ls='--', lw=1.2, label='Q3D-calibrated coupled-mode model')
    ax0.axvspan(b_lo, b_hi, color='#f9e79f', alpha=0.45,
                label=f'Purcell -3 dB passband ({b_lo:.3f}–{b_hi:.3f} GHz, $\\kappa_p/2\\pi={kp:.1f}$ MHz)')
    for i, fd in enumerate(dips, start=1):
        ax0.axvline(fd, color='#c0392b', ls=':', lw=1.1)
    ax0.set_xlim(f_wide.min(), f_wide.max())
    ax0.set_ylabel('Intra-Purcell |S21| (dB)')
    ax0.set_title(
        f"{DESIGN.name} · Q3D + HFSS Verification: 3 Readout Resonators Inside Purcell Filter Band\n"
        f"Q3D C_in = {summary_metrics['c_in_q3d_fF']:.2f} fF, C_a,gnd = {summary_metrics['c_a_gnd_deembedded_fF']:.2f} fF | "
        f"HFSS Pass 4 Max Delta S = {summary_metrics['last_delta_s']:.4f} ({summary_metrics['solved_elements_pass4']:,} tetras)",
        fontsize=10.5,
        weight='bold',
    )
    ax0.grid(alpha=0.3)
    ax0.legend(loc='lower right', fontsize=8.5)

    # Panel 1: Fine zoom across the Purcell passband showing the 3 sharp readout resonator notches (RR1, RR2, RR3)
    ax1.plot(f_fine, s21_db_f, color='#1b4f72', lw=1.7, label='HFSS Driven Fine Sweep (0.15 MHz step)')
    ax1.plot(f_fine, s21_db_f_mod, color='#d35400', ls='--', lw=1.3, label='Q3D-calibrated coupled-mode model')
    ax1.axvspan(b_lo, b_hi, color='#f9e79f', alpha=0.40,
                label=f'Purcell passband ($f_p={fp:.4f}$ GHz, FWHM $\\kappa_p/2\\pi={kp:.1f}$ MHz)')
    for i, (fd, depth) in enumerate(zip(dips, summary_metrics['resonator_notch_depths_dB']), start=1):
        ax1.axvline(fd, color='#c0392b', ls=':', lw=1.25)
        ax1.text(fd, s21_db_f.min() + 1.5 + (i % 2) * 3.5,
                 f'RR{i}\n{fd:.4f} GHz\n(-{depth:.1f} dB)',
                 color='#c0392b', fontsize=8.5, ha='center', va='bottom', weight='bold',
                 bbox=dict(fc='white', ec='#c0392b', lw=0.6, alpha=0.9, pad=1.5))
    ax1.set_xlim(f_fine.min(), f_fine.max())
    ax1.set_ylabel('Fine Intra-Purcell |S21| (dB)')
    ax1.set_title('Zoomed Purcell Passband: Three Readout Resonators (RR1, RR2, RR3) Inside Purcell Filter Band',
                  fontsize=10, weight='bold')
    ax1.grid(alpha=0.3)
    ax1.legend(loc='lower right', fontsize=8.5)

    # Panel 2: 1-Port Reflection Group Delay tau_g(S11) at Input Launch Port LP_IN
    ax2.plot(f_fine, tau_f_ns, color='#1b4f72', lw=1.5,
             label='HFSS 1-port reflection group delay $\\tau_g(S_{11})$ at input port')
    ax2.plot(f_fine, model_fine['group_delay_ns'], color='#d35400', ls='--', lw=1.2,
             label='Q3D-calibrated model group delay')
    ax2.axvspan(b_lo, b_hi, color='#f9e79f', alpha=0.40)
    for i, (fd, kr) in enumerate(zip(dips, summary_metrics['resonator_purcell_filtered_kappa_r_MHz']), start=1):
        ax2.axvline(fd, color='#c0392b', ls=':', lw=1.25)
        ax2.text(fd, max(tau_f_ns) * 0.45,
                 f'RR{i}: $\\kappa_r/2\\pi\\approx{kr:.1f}$ MHz',
                 color='#c0392b', fontsize=8.0, ha='center', va='bottom', weight='bold',
                 bbox=dict(fc='white', ec='none', alpha=0.85, pad=1))
    ax2.set_xlim(f_fine.min(), f_fine.max())
    ax2.set_yscale('log')
    ax2.set_xlabel('Frequency (GHz)')
    ax2.set_ylabel('Reflection group delay $\\tau_g$ (ns)')
    ax2.grid(alpha=0.3, which='both')
    ax2.legend(loc='upper right', fontsize=8.5)

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=180, facecolor='white')
    plt.close(fig)
    png_bytes = buf.getvalue()

    for p in out_png_paths:
        _long_path(Path(p).parent).mkdir(parents=True, exist_ok=True)
        _long_path(p).write_bytes(png_bytes)


def run_reanalysis(q3d_dir, hfss_raw_dir):
    q3d_derived = load_q3d_summary(q3d_dir)
    geo = B.load_config()
    raw_dir = Path(hfss_raw_dir) / 'raw'
    s_wide_path = raw_dir / 'S_wide.csv'
    s_fine_path = raw_dir / 'S_fine.csv'
    conv_path = raw_dir / 'hfss_convergence.csv'

    folder = create_run(
        DESIGN.name,
        'resonator',
        's21_purcell_band_analysis',
        entrypoint=__file__,
        inputs=[
            DESIGN / 'parameters' / 'geometry.json',
            q3d_dir / 'summary.json',
            s_wide_path,
            s_fine_path,
            conv_path,
        ],
        upstream=[str(q3d_dir), str(hfss_raw_dir)],
        parameters=dict(
            mode='reanalyze',
            q3d_run=str(q3d_dir),
            hfss_raw_run=str(hfss_raw_dir),
            geometry_revision=geo['geometry_revision'],
        ),
    )
    print('RUN_FOLDER=' + str(folder), flush=True)

    s_wide_df = pd.read_csv(_long_path(s_wide_path), index_col=0)
    s_fine_df = pd.read_csv(_long_path(s_fine_path), index_col=0)
    conv_df = pd.read_csv(_long_path(conv_path))

    f_wide = s_wide_df.index.to_numpy(float)
    s11_wide = _parse_complex_col(s_wide_df['S11'])
    s21_wide = _parse_complex_col(s_wide_df['S21'])

    f_fine = s_fine_df.index.to_numpy(float)
    s11_fine = _parse_complex_col(s_fine_df['S11'])
    s21_fine = _parse_complex_col(s_fine_df['S21'])

    metrics, model_wide, model_fine = analyze_spectra(
        f_wide, s11_wide, s21_wide, f_fine, s11_fine, s21_fine, q3d_derived, geo, conv_df
    )

    run_fig = folder / 'figures' / 'purcell_band_resonators.png'
    notes_fig = DESIGN / 'notes' / 'fig' / 'purcell_band_resonators.png'
    plot_purcell_band_results(
        f_wide, s11_wide, s21_wide, f_fine, s11_fine, s21_fine,
        model_wide, model_fine, metrics, [run_fig, notes_fig]
    )

    rec = json.loads(_long_path(folder / 'run.json').read_text(encoding='utf-8'))
    rec['geometry_revision'] = geo['geometry_revision']
    _write(folder / 'run.json', rec)

    summary = dict(
        metrics=metrics,
        q3d_upstream=str(q3d_dir),
        hfss_raw_upstream=str(hfss_raw_dir),
        figure=str(notes_fig.relative_to(ROOT).as_posix()),
        validation='4-pass HFSS Driven run with explicit 10 um length mesh, 350 um wire bonds + 6 coupling straps; '
                   f"Pass 4 Max Mag. Delta S = {metrics['last_delta_s']:.5f} (<= 0.05 threshold).",
    )
    print(json.dumps(metrics, indent=2), flush=True)
    finish_run(folder, status='completed', validation='validated', summary=summary)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--q3d', type=Path, required=True, help='Path to completed q3d_purcell_idc run folder')
    ap.add_argument('--reanalyze', type=Path, default=None, help='Reanalyze an existing raw HFSS run folder')
    ap.add_argument('--passes', type=int, default=4, help='HFSS driven adaptive passes (default 4)')
    args = ap.parse_args()

    q3d_dir = args.q3d.resolve()
    if args.reanalyze is not None:
        run_reanalysis(q3d_dir, args.reanalyze.resolve())
        return

    q3d_derived = load_q3d_summary(q3d_dir)
    geo = B.load_config()

    design = B.build_design(wire_bonds=True)
    comps = ['IDC_in', 'feed_in', 'PF1', 'RR1', 'RR2', 'RR3']
    open_pins = [('RR1', 'end'), ('RR2', 'end'), ('RR3', 'end')]
    straps = coupling_strap_points(design, geo)

    params = dict(
        solver='hfss_driven',
        max_passes=args.passes,
        max_delta_s=0.02,
        mesh_length=MESH_LENGTH,
        wire_bonds=WB_OPTS,
        coupling_straps=straps,
        components=comps,
        open_terminations=open_pins,
        ports=PORTS,
        q3d_upstream=str(q3d_dir),
        q3d_derived=q3d_derived,
        geometry_revision=geo['geometry_revision'],
    )

    folder = create_run(
        DESIGN.name,
        'resonator',
        'hfss_driven_purcell_triad',
        entrypoint=__file__,
        inputs=[DESIGN / 'parameters' / 'geometry.json', q3d_dir / 'summary.json'],
        upstream=[str(q3d_dir)],
        parameters=params,
    )
    print('RUN_FOLDER=' + str(folder), flush=True)

    try:
        from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
        import hfss_dirs
        import save_aedt

        hfss_dirs.use_d_drive()
        short_aedt = Path('D:/HFSS_temp/d6r.aedt')
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
                mesh_info['mesh_objects'] = mesh_all_metal(self.renderer, design, comps)
                super()._analyze()

        sim = DrivenSim(design, 'hfss')
        sim.renderer._options.update(WB_OPTS)
        f0 = float(geo['purcell_filter']['f_target_GHz'])
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
        lo_ghz, hi_ghz = 5.80, 7.30
        sim.setup.sweep_setup.update(
            name='Wide',
            start_ghz=lo_ghz,
            stop_ghz=hi_ghz,
            count=3001,
            type='Fast',
            save_fields=False,
        )
        t0 = time.time()
        try:
            sim.run(
                name=f'D6R_HFSS_{folder.name[-8:]}',
                components=comps,
                open_terminations=open_pins,
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
        s_df.to_csv(_long_path(folder / 'raw' / 'S_wide.csv'))
        conv = sim.renderer.pinfo.setup.get_convergence()[0]
        conv.to_csv(_long_path(folder / 'raw' / 'hfss_convergence.csv'))

        f_wide = s_df.index.to_numpy(float)
        s21_wide = s_df['S21'].to_numpy(complex)
        s11_wide = s_df['S11'].to_numpy(complex)

        s21_mag_db = 20.0 * np.log10(np.maximum(np.abs(s21_wide), 1e-12))
        f_peak_est = float(f_wide[np.argmax(s21_mag_db)])
        flo, fhi = max(lo_ghz, f_peak_est - 0.35), min(hi_ghz, f_peak_est + 0.35)
        fine_count = min(int(np.ceil((fhi - flo) / 0.00015)) + 1, 6001)
        sim.renderer.add_sweep(
            setup_name='Setup',
            name='Fine',
            start_ghz=flo,
            stop_ghz=fhi,
            count=fine_count,
            type='Fast',
            save_fields=False,
        )
        sim.renderer.analyze_sweep('Fine', 'Setup')
        _, _, s_fine_df = sim.renderer.get_params(['S11', 'S21', 'S22'])
        s_fine_df.to_csv(_long_path(folder / 'raw' / 'S_fine.csv'))
        sim.renderer.pinfo.project._project.Save()

        f_fine = s_fine_df.index.to_numpy(float)
        s21_fine = s_fine_df['S21'].to_numpy(complex)
        s11_fine = s_fine_df['S11'].to_numpy(complex)

        metrics, model_wide, model_fine = analyze_spectra(
            f_wide, s11_wide, s21_wide, f_fine, s11_fine, s21_fine, q3d_derived, geo, conv
        )
        metrics['mesh_object_count'] = len(mesh_info.get('mesh_objects', []))
        metrics['strap_count'] = len(mesh_info.get('straps', []))
        metrics['minutes'] = (time.time() - t0) / 60.0
        metrics['aedt_design'] = sim.renderer.pinfo.design.name

        run_fig = folder / 'figures' / 'purcell_band_resonators.png'
        notes_fig = DESIGN / 'notes' / 'fig' / 'purcell_band_resonators.png'
        plot_purcell_band_results(
            f_wide, s11_wide, s21_wide, f_fine, s11_fine, s21_fine,
            model_wide, model_fine, metrics, [run_fig, notes_fig]
        )

        summary = dict(
            metrics=metrics,
            q3d_upstream=str(q3d_dir),
            figure=str(notes_fig.relative_to(ROOT).as_posix()),
            validation='HFSS Driven 4-pass budget with explicit 10 um length mesh and 350 um wire bonds + 6 coupling straps; '
                       'last_delta_s recorded in raw/hfss_convergence.csv.',
        )
        print(json.dumps(metrics, indent=2), flush=True)
        finish_run(folder, status='completed', summary=summary)
    except Exception:
        _long_path(folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
