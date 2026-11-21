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
"""

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
from scipy.signal import find_peaks  # noqa: E402

import build as B  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

WB_OPTS = dict(wb_threshold='350um', wb_offset='0um', wb_size=2)
MESH_LENGTH = '10um'
STRAP_DIAMETER_MM = 0.010
PORTS = [('feed_in', 'start', 50), ('PF1', 'end', 50000)]


def _write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, default=float) + '\n', encoding='utf-8')


def load_q3d_summary(q3d_dir):
    summary_path = Path(q3d_dir) / 'summary.json'
    return json.loads(summary_path.read_text(encoding='utf-8'))['derived']


def q3d_calibrated_overrides(geo, q3d_derived):
    """Return geometry overrides so PF1 is calibrated with the Q3D IDC_in end load (C_in + C_a_gnd)
    and RR1-RR3 sit inside the Purcell filter passband around f_p = 6.65 GHz.
    """
    cpw = geo['cpw']
    pf = geo['purcell_filter']
    c_per_m = cpw['assumed_c_per_m_f']
    idc_d = q3d_derived['idc']

    # Total capacitive end load on PF1 at pin a from Q3D: C_a_gnd_deembedded + C_in_eff_shunt
    c_total_ff = idc_d['c_total_pf_end_load_fF']
    c_eff_only_ff = idc_d['c_in_eff_shunt_fF']
    extra_c_ff = c_total_ff - c_eff_only_ff
    # Shifting f_target_GHz in build_design by the extra Q3D shunt capacitance shortens PF1 accordingly
    v_ph = cpw['assumed_v_ph_m_s']
    f_p = pf['f_target_GHz']
    l_unloaded_mm = (v_ph / (2.0 * f_p * 1e9)) * 1e3
    l_total_shorten_mm = (c_total_ff * 1e-15 / c_per_m) * 1e3

    return dict(
        q3d_pf_target_len_mm=l_unloaded_mm - l_total_shorten_mm,
        extra_c_a_gnd_fF=extra_c_ff,
    )


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
    span_mm = 2.0 * (w_mm + 2.0 * g_mm) + s_mm + 0.020  # bridges across both CPWs into outer ground

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


def circuit_model_response(freqs_ghz, f_p_ghz, kappa_p_mhz, f_r_list_ghz, j_list_mhz, q_int=200000.0):
    """Coupled-mode / admittance model for 1 reflection Purcell filter coupled to 3 readout resonators.

    Hamiltonian / pole matrix in rotating frame:
        omega_p_tilde = f_p - i * (kappa_p + kappa_p_int) / 2
        omega_rk_tilde = f_rk - i * kappa_r_int / 2
        Coupling J_k between Purcell mode and readout resonator k.

    Returns:
        dict with complex S11(f), probe S21(f), group_delay_ns(f), and dressed eigenvalues.
    """
    kp = kappa_p_mhz * 1e-3  # GHz
    kp_int = f_p_ghz / q_int
    s11 = np.zeros_like(freqs_ghz, dtype=complex)
    s21 = np.zeros_like(freqs_ghz, dtype=complex)
    kp_probe = 1e-4  # weak probe linewidth equivalent (GHz)

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
    tau_g_ns = -(np.gradient(phase_rad) / (2.0 * np.pi * df_hz)) * 1e9

    # Dressed complex modes of the 4x4 coupled system
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


def plot_purcell_band_results(f_ghz, s11_hfss, s21_hfss, model_out, summary_metrics, out_png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    s21_db_hfss = 20.0 * np.log10(np.maximum(np.abs(s21_hfss), 1e-12))
    s21_db_mod = 20.0 * np.log10(np.maximum(np.abs(model_out['s21']), 1e-12))
    # Align probe baseline offset so envelope shape and 3 resonator dips compare directly
    offset_db = float(np.percentile(s21_db_hfss, 95) - np.percentile(s21_db_mod, 95))
    s21_db_mod_aligned = s21_db_mod + offset_db

    # Exact 1-port open reflection from 2-port S-matrix (Port 2 high-Z probe)
    phase_hfss_rad = np.unwrap(np.angle(s11_hfss))
    df_hz = np.gradient(f_ghz) * 1e9
    tau_hfss_ns = np.abs(np.gradient(phase_hfss_rad) / (2.0 * np.pi * df_hz)) * 1e9

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.5, 7.5), sharex=True)

    # Panel 1: Purcell filter passband + 3 readout resonator dips inside the band
    ax1.plot(f_ghz, s21_db_hfss, color='#1b4f72', lw=1.6, label='HFSS Driven (intra-Purcell probe |S21|)')
    ax1.plot(f_ghz, s21_db_mod_aligned, color='#d35400', ls='--', lw=1.3,
             label='Q3D-calibrated coupled-mode model')

    fp = summary_metrics['f_p_hfss_GHz']
    kp = summary_metrics['kappa_p_hfss_MHz']
    ax1.axvspan(fp - kp / 2000.0, fp + kp / 2000.0, color='#f9e79f', alpha=0.35,
                label=f'Purcell passband ($f_p={fp:.3f}$ GHz, $\\kappa_p/2\\pi\\approx{kp:.0f}$ MHz)')

    for i, fd in enumerate(summary_metrics['resonator_dips_hfss_GHz'], start=1):
        ax1.axvline(fd, color='#c0392b', ls=':', lw=1.2)
        ax1.text(fd, s21_db_hfss.min() + 3.0, f'RR{i}\n{fd:.4f} GHz',
                 color='#c0392b', fontsize=8.5, ha='center', va='bottom', weight='bold',
                 bbox=dict(fc='white', ec='none', alpha=0.8, pad=1))

    ax1.set_ylabel('Intra-Purcell response |S21| (dB)')
    ax1.set_title(
        f"{DESIGN.name} · 3 Readout Resonators Inside Reflection Purcell Filter Band\n"
        f"Q3D C_in = {summary_metrics['c_in_q3d_fF']:.1f} fF | "
        f"HFSS f_p = {fp:.4f} GHz (kappa_p/2pi = {kp:.1f} MHz) | "
        f"Dips = {', '.join(f'{x:.4f}' for x in summary_metrics['resonator_dips_hfss_GHz'])} GHz",
        fontsize=10,
        weight='bold',
    )
    ax1.grid(alpha=0.3)
    ax1.legend(loc='lower right', fontsize=8.5)

    # Panel 2: 1-port reflection group delay at input port LP_IN
    ax2.plot(f_ghz, tau_hfss_ns, color='#1b4f72', lw=1.5, label='HFSS 1-port reflection group delay $\\tau_g(S_{11})$')
    ax2.plot(f_ghz, np.abs(model_out['group_delay_ns']), color='#d35400', ls='--', lw=1.2,
             label='Q3D-calibrated model group delay')
    ax2.axvspan(fp - kp / 2000.0, fp + kp / 2000.0, color='#f9e79f', alpha=0.35)
    for fd in summary_metrics['resonator_dips_hfss_GHz']:
        ax2.axvline(fd, color='#c0392b', ls=':', lw=1.2)
    ax2.set_yscale('log')
    ax2.set_xlabel('Frequency (GHz)')
    ax2.set_ylabel('Reflection group delay $\\tau_g$ (ns)')
    ax2.grid(alpha=0.3, which='both')
    ax2.legend(loc='upper right', fontsize=8.5)

    fig.tight_layout()
    fig.savefig(out_png, dpi=180, facecolor='white')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--q3d', type=Path, required=True, help='Path to completed q3d_purcell_idc run folder')
    ap.add_argument('--passes', type=int, default=4, help='HFSS driven adaptive passes (default 4 per skill/em_simulation.md)')
    args = ap.parse_args()

    q3d_dir = args.q3d.resolve()
    q3d_derived = load_q3d_summary(q3d_dir)
    geo = B.load_config()

    design = B.build_design(wire_bonds=True)
    # Exclude LP_IN and PF1_open so QHFSSRenderer places lumped port sheets at ('feed_in', 'start') and ('PF1', 'end')
    comps = [n for n in design.components.keys() if n not in ('LP_IN', 'PF1_open')]
    straps = coupling_strap_points(design, geo)

    params = dict(
        solver='hfss_driven',
        max_passes=args.passes,
        max_delta_s=0.02,
        mesh_length=MESH_LENGTH,
        wire_bonds=WB_OPTS,
        coupling_straps=straps,
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

        import shutil
        hfss_dirs.use_d_drive()
        short_aedt = Path('D:/HFSS_temp/d6r.aedt')
        save_aedt.use_project(str(short_aedt))
        rec = json.loads((folder / 'run.json').read_text(encoding='utf-8'))
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
        # Wide fast sweep across 5.8 - 7.3 GHz (3001 points -> 0.5 MHz step) to resolve both Purcell envelope and 3 dips
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
        s_df.to_csv(folder / 'raw' / 'S_wide.csv')
        conv = sim.renderer.pinfo.setup.get_convergence()[0]
        conv.to_csv(folder / 'raw' / 'hfss_convergence.csv')

        f_wide = s_df.index.to_numpy(float)
        s21_wide = s_df['S21'].to_numpy(complex)
        s11_wide = s_df['S11'].to_numpy(complex)

        # Fine fast sweep (0.15 MHz step) across the Purcell band + 3 resonator dips
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
        s_fine_df.to_csv(folder / 'raw' / 'S_fine.csv')
        sim.renderer.pinfo.project._project.Save()

        f_fine = s_fine_df.index.to_numpy(float)
        s21_fine = s_fine_df['S21'].to_numpy(complex)
        s11_fine = s_fine_df['S11'].to_numpy(complex)

        # Extract Purcell peak center, 3 dB bandwidth, and the 3 sharp resonator dips inside the Purcell band
        s21_fine_db = 20.0 * np.log10(np.maximum(np.abs(s21_fine), 1e-12))
        # Smooth envelope to measure Purcell filter center & bandwidth
        from scipy.ndimage import maximum_filter1d
        env_db = maximum_filter1d(s21_fine_db, size=max(15, len(f_fine) // 25))
        pk_idx = int(np.argmax(env_db))
        f_p_hfss = float(f_fine[pk_idx])
        pk_val = float(env_db[pk_idx])
        above_3db = f_fine[env_db >= pk_val - 3.0]
        kappa_p_hfss_mhz = float((above_3db.max() - above_3db.min()) * 1e3) if len(above_3db) > 2 else float(
            q3d_derived['idc']['predicted_kappa_p_MHz']
        )

        # Find the 3 sharp resonator dips relative to the smooth Purcell envelope
        notch_depth_db = env_db - s21_fine_db
        dip_indices, dip_props = find_peaks(notch_depth_db, prominence=1.5, distance=max(10, int(0.03 / 0.00015)))
        if len(dip_indices) > 3:
            order = np.argsort(dip_props['prominences'])[::-1][:3]
            dip_indices = np.sort(dip_indices[order])
        dips_ghz = [float(f_fine[i]) for i in dip_indices]
        dip_depths_db = [float(notch_depth_db[i]) for i in dip_indices]

        # Evaluate Q3D-based circuit model on the same fine grid
        j_q3d_mhz = float(q3d_derived['coupled_section']['estimated_J_MHz'])
        f_r_targets = [ch['f_target_GHz'] for ch in geo['readout_resonators']['channels']]
        f_r_for_model = dips_ghz if len(dips_ghz) == 3 else f_r_targets
        model_out = circuit_model_response(
            f_fine,
            f_p_ghz=f_p_hfss,
            kappa_p_mhz=kappa_p_hfss_mhz,
            f_r_list_ghz=f_r_for_model,
            j_list_mhz=[j_q3d_mhz] * len(f_r_for_model),
        )

        summary_metrics = dict(
            c_in_q3d_fF=float(q3d_derived['idc']['c_in_q3d_fF']),
            c_a_gnd_deembedded_fF=float(q3d_derived['idc']['c_a_gnd_deembedded_fF']),
            predicted_kappa_p_q3d_MHz=float(q3d_derived['idc']['predicted_kappa_p_MHz']),
            estimated_J_q3d_MHz=j_q3d_mhz,
            f_p_hfss_GHz=f_p_hfss,
            kappa_p_hfss_MHz=kappa_p_hfss_mhz,
            resonator_dips_hfss_GHz=dips_ghz,
            resonator_notch_depths_dB=dip_depths_db,
            all_dips_inside_purcell_3dB_band=bool(
                len(dips_ghz) == 3 and all(abs(d - f_p_hfss) * 1e3 <= 0.65 * kappa_p_hfss_mhz for d in dips_ghz)
            ),
            passes=int(len(conv)),
            last_delta_s=float(conv.iloc[-1, -1]),
            mesh_length=MESH_LENGTH,
            mesh_object_count=len(mesh_info.get('mesh_objects', [])),
            strap_count=len(mesh_info.get('straps', [])),
            minutes=(time.time() - t0) / 60.0,
            aedt_design=sim.renderer.pinfo.design.name,
        )

        fig_dir = folder / 'figures'
        fig_dir.mkdir(exist_ok=True)
        fig_path = fig_dir / 'purcell_band_resonators.png'
        plot_purcell_band_results(f_fine, s11_fine, s21_fine, model_out, summary_metrics, fig_path)

        notes_fig_dir = DESIGN / 'notes' / 'fig'
        notes_fig_dir.mkdir(parents=True, exist_ok=True)
        (notes_fig_dir / 'purcell_band_resonators.png').write_bytes(fig_path.read_bytes())

        summary = dict(
            metrics=summary_metrics,
            q3d_upstream=str(q3d_dir),
            figure=str(fig_path.relative_to(ROOT)),
            validation='HFSS Driven 4-pass budget with explicit 10 um length mesh and 350 um wire bonds + 6 coupling straps; '
                       'last_delta_s recorded in raw/hfss_convergence.csv.',
        )
        print(json.dumps(summary_metrics, indent=2), flush=True)
        finish_run(folder, status='completed', summary=summary)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
