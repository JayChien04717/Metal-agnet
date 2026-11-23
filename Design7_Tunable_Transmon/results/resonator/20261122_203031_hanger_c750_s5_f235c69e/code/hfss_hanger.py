"""Single-channel HFSS driven S21 simulation for Design7_Tunable_Transmon lambda/4 hanger resonator.

Follows skill/em_simulation.md defaults:
- 4 adaptive passes, max_delta_s = 0.01, first-order basis
- Explicit 10 um length-based mesh on feedline and resonator traces
- Ground wire bonds on every etched path (wb_threshold = '350um')
- Broadband sweep to locate dip + fine sweep around dip to extract f_r, Q_c, and kappa/2pi

Usage (from repository root, qmetal Python, AEDT running):
    python -X utf8 -B Design7_Tunable_Transmon/simulation/resonator/hfss_hanger.py \
        --tag hanger_base --f-target 6.5 --coup-len 300 --coup-strip 8
"""
from argparse import ArgumentParser
from collections import OrderedDict
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
sys.path[:0] = [str(DESIGN / 'layout'), *project_paths.pattern_design_paths()]

from qiskit_metal import Dict, designs  # noqa: E402
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround  # noqa: E402
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight  # noqa: E402
import build as d7  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402


def read_s2p(path):
    """Read a 2-port Touchstone file (MA / RI / DB) and return (freq_ghz, s11, s21)."""
    fmt = 'MA'
    freq_scale = 1.0
    rows = []
    with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('!'):
                continue
            if s.startswith('#'):
                toks = s.upper().split()
                for t in ('MA', 'RI', 'DB'):
                    if t in toks:
                        fmt = t
                if 'HZ' in toks and 'GHZ' not in toks and 'MHZ' not in toks:
                    freq_scale = 1e-9
                elif 'MHZ' in toks:
                    freq_scale = 1e-3
                continue
            vals = [float(x) for x in s.split()]
            if len(vals) >= 9:
                rows.append(vals[:9])
    arr = np.asarray(rows, float)
    freq = arr[:, 0] * freq_scale

    def pair(a, b):
        if fmt == 'RI':
            return a + 1j * b
        if fmt == 'MA':
            return a * np.exp(1j * np.deg2rad(b))
        # DB
        return (10.0 ** (a / 20.0)) * np.exp(1j * np.deg2rad(b))

    return freq, pair(arr[:, 1], arr[:, 2]), pair(arr[:, 3], arr[:, 4])


def export_sweep_s2p(design_obj, setup_name, sweep_name, out_path):
    sol = design_obj._design.GetModule('Solutions')
    out_str = str(Path(out_path).resolve())
    sol.ExportNetworkData('', [f'{setup_name}:{sweep_name}'], 3, out_str, ['All'],
                          True, 50, 'S', -1, 0, 15, True, True, False)
    return read_s2p(out_str)


def fit_notch_s21(freq_ghz, s21_complex):
    """Extract f_r, Q_c, and kappa/2pi from complex S21(f) of a lossless hanger resonator.

    For a lossless side-coupled hanger resonator on a matched feedline:
        S21(f) = A * exp(i*(phi0 + tau*(f - f0))) * (1 - (Q_l / |Q_c|) * exp(i*phi) / (1 + 2i Q_l (f/f_r - 1)))
    For lossless PEC (Q_i -> inf, Q_l = Q_c), |S21|^2 has a Lorentzian dip to ~0 with FWHM = f_r / Q_c = kappa/2pi.
    """
    from scipy.optimize import least_squares

    f = np.asarray(freq_ghz, float)
    s = np.asarray(s21_complex, complex)
    mag = np.abs(s)
    idx_min = int(np.argmin(mag))
    f0_init = float(f[idx_min])

    # Estimate FWHM from half-power points (-3 dB relative to baseline)
    base = float(np.median(mag))
    half_lvl = base / math.sqrt(2.0)
    below = np.where(mag < half_lvl)[0]
    if len(below) >= 2:
        fwhm_init = max(float(f[below[-1]] - f[below[0]]), 1e-4)
    else:
        fwhm_init = 0.001

    # Fit complex notch model in window around dip (+- 10 * fwhm or full fine sweep)
    def model(params, freq):
        fr, ql, qc_mag, phi, a_re, a_im, tau = params
        bg = (a_re + 1j * a_im) * np.exp(1j * tau * (freq - fr))
        denom = 1.0 + 2j * ql * (freq / fr - 1.0)
        return bg * (1.0 - (ql / qc_mag) * np.exp(1j * phi) / denom)

    ql_init = f0_init / fwhm_init
    bg_init = s[0] if abs(f[0] - f0_init) > abs(f[-1] - f0_init) else s[-1]
    p0 = [f0_init, ql_init, ql_init, 0.0, bg_init.real, bg_init.imag, 0.0]
    lb = [f.min(), 100.0, 100.0, -math.pi, -2.0, -2.0, -200.0]
    ub = [f.max(), 5e6, 5e6, math.pi, 2.0, 2.0, 200.0]

    def resid(p):
        diff = model(p, f) - s
        return np.concatenate([diff.real, diff.imag])

    sol = least_squares(resid, p0, bounds=(lb, ub), max_nfev=4000)
    fr, ql, qc_mag, phi, a_re, a_im, tau = [float(x) for x in sol.x]
    qc = qc_mag / max(math.cos(phi), 0.1)
    kappa_mhz = (fr / qc) * 1e3
    kappa_l_mhz = (fr / ql) * 1e3
    rmse = float(np.sqrt(np.mean(np.abs(model(sol.x, f) - s) ** 2)))
    return dict(
        f_r_GHz=fr,
        Q_l=ql,
        Q_c=qc,
        phi_rad=phi,
        kappa_c_MHz=kappa_mhz,
        kappa_l_MHz=kappa_l_mhz,
        f_min_raw_GHz=f0_init,
        s21_min_dB=float(20.0 * np.log10(max(mag[idx_min], 1e-12))),
        fit_rmse=rmse,
    )


def build_hanger_fixture(f_target_ghz, coup_len_um, coup_strip_um, v_ph=1.2016e8):
    geo, _ = d7.load_config()
    geo = json.loads(json.dumps(geo))
    geo['hanger_resonator']['coupling_length_um'] = float(coup_len_um)
    geo['hanger_resonator']['coupling_ground_strip_um'] = float(coup_strip_um)

    d = designs.DesignPlanar(overwrite_enabled=True)
    d.chips.main.size.size_x = '3.2mm'
    d.chips.main.size.size_y = '2.8mm'
    d.chips.main.size.size_z = '-500um'
    d.chips.main.material = 'silicon'

    cpw_cfg = geo['cpw']
    cpw = dict(trace_width=f"{cpw_cfg['trace_width_um']:g}um",
               trace_gap=f"{cpw_cfg['trace_gap_um']:g}um")

    y_feed = -0.8
    OpenToGround(d, 'feed_in', options=dict(
        pos_x='-1.35mm', pos_y=f'{y_feed}mm', orientation='180',
        width=cpw['trace_width'], gap=cpw['trace_gap'], termination_gap='6um',
    ))
    OpenToGround(d, 'feed_out', options=dict(
        pos_x='1.35mm', pos_y=f'{y_feed}mm', orientation='0',
        width=cpw['trace_width'], gap=cpw['trace_gap'], termination_gap='6um',
    ))
    RouteStraight(d, 'feedline', options=Dict(
        chip='main',
        pin_inputs=Dict(
            start_pin=Dict(component='feed_in', pin='open'),
            end_pin=Dict(component='feed_out', pin='open'),
        ),
        **cpw,
    ))

    l_target = d7.quarter_wave_mm(f_target_ghz, v_ph=v_ph, load_shorten_mm=0.0)
    rep = d7.add_hanger_resonator(d, 'R_test', 0.0, y_feed, +1, l_target, geo, end_pin=None)

    # Enable ground wire bonds on all CPW paths (skill/em_simulation.md section 3)
    d.qgeometry.tables['path']['hfss_wire_bonds'] = True
    return d, rep, geo


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--tag', default='hfss_hanger')
    ap.add_argument('--f-target', type=float, default=6.5)
    ap.add_argument('--coup-len', type=float, default=300.0)
    ap.add_argument('--coup-strip', type=float, default=8.0)
    ap.add_argument('--max-passes', type=int, default=4)
    ap.add_argument('--mesh-um', type=float, default=10.0)
    ap.add_argument('--f-start', type=float, default=6.0)
    ap.add_argument('--f-stop', type=float, default=7.0)
    a = ap.parse_args()

    d, rep, geo = build_hanger_fixture(a.f_target, a.coup_len, a.coup_strip)
    folder = create_run(
        DESIGN.name, 'resonator', a.tag, entrypoint=__file__,
        inputs=[DESIGN / 'design.json', DESIGN / 'parameters' / 'geometry.json',
                DESIGN / 'layout' / 'build.py'],
        parameters=dict(
            f_target_GHz=a.f_target,
            coupling_length_um=a.coup_len,
            coupling_ground_strip_um=a.coup_strip,
            resonator_length_mm=rep['actual_filleted_mm'],
            hfss_max_passes=a.max_passes,
            hfss_max_delta_s=0.01,
            explicit_mesh_length_um=a.mesh_um,
            mesh_targets=['trace_feedline', 'trace_R_test'],
            wire_bonds=dict(wb_threshold='350um', wb_offset='0um', wb_size=2),
            broad_sweep_GHz=[a.f_start, a.f_stop],
        ),
    )
    print('RUN_FOLDER=' + str(folder), flush=True)

    try:
        from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
        import save_aedt
        import hfss_dirs

        temp = hfss_dirs.use_d_drive()
        # Use short project/design name (D7.aedt) so AEDT .aedtresults paths stay below Windows 260-char MAX_PATH
        proj = DESIGN / 'simulation' / 'aedt' / 'D7.aedt'
        save_aedt.use_project(str(proj))

        sim = ScatteringImpedanceSim(d, 'hfss')
        hfss = sim.renderer
        hfss.start()
        hfss.options['wb_threshold'] = '350um'
        hfss.options['wb_offset'] = '0um'
        hfss.options['wb_size'] = 1
        hfss.options['x_buffer_width_mm'] = 0.2
        hfss.options['y_buffer_width_mm'] = 0.2

        dname = f'R_{folder.name[-4:]}'
        hfss.activate_ansys_design(dname, 'drivenmodal')
        hfss.render_design(
            selection=['feedline', 'R_test'],
            open_pins=[('R_test', 'end')],
            port_list=[('feedline', 'start', 50), ('feedline', 'end', 50)],
            box_plus_buffer=False,
        )

        # Explicit 10 um length-based mesh on CPW traces (skill/em_simulation.md section 2)
        pinfo = hfss.pinfo
        oModeler = pinfo.design.modeler
        all_objs = list(oModeler._modeler.GetObjectsInGroup('Sheets'))
        mesh_objs = [o for o in all_objs if 'feedline' in o or 'R_test' in o]
        mesh_objs = [o for o in mesh_objs if not o.startswith('wire')]
        print('Mesh objects:', mesh_objs, flush=True)

        setup = pinfo.design.create_dm_setup(
            freq_ghz=a.f_target,
            name='Setup1',
            max_delta_s=0.01,
            max_passes=a.max_passes,
            min_passes=2,
            min_converged=1,
            pct_refinement=30,
            basis_order=1,
        )
        if mesh_objs:
            oModeler.mesh_length('cpw_mesh', mesh_objs, MaxLength=f'{a.mesh_um:g}um')

        t0 = time.time()
        # 1. Broadband fast sweep to locate the resonance dip
        setup.insert_sweep(
            a.f_start, a.f_stop, count=1001, name='SweepBroad', type='Fast', save_fields=False
        )
        pinfo.design._design.Analyze('Setup1 : SweepBroad')
        freqs_b, s11_b, s21_b = export_sweep_s2p(
            pinfo.design, 'Setup1', 'SweepBroad', folder / 'raw' / 'broad.s2p'
        )
        f_dip = float(freqs_b[int(np.argmin(np.abs(s21_b)))])
        print(f'Broad sweep dip at {f_dip:.5f} GHz '
              f'(min |S21| = {20*np.log10(max(np.min(np.abs(s21_b)), 1e-12)):.2f} dB, '
              f'med |S21| = {20*np.log10(max(np.median(np.abs(s21_b)), 1e-12)):.2f} dB)', flush=True)

        pd.DataFrame({'freq_ghz': freqs_b, 's21_re': s21_b.real, 's21_im': s21_b.imag}).to_csv(
            folder / 'raw' / 'S21_broad.csv', index=False
        )

        # 2. Fine sweep around dip (+- 15 MHz, 601 points -> 0.05 MHz step)
        f_lo, f_hi = round(f_dip - 0.015, 5), round(f_dip + 0.015, 5)
        setup.insert_sweep(
            f_lo, f_hi, count=601, name='SweepFine', type='Fast', save_fields=False
        )
        pinfo.design._design.Analyze('Setup1 : SweepFine')
        freqs_f, s11_f, s21_f = export_sweep_s2p(
            pinfo.design, 'Setup1', 'SweepFine', folder / 'raw' / 'fine.s2p'
        )
        elapsed = time.time() - t0
        pd.DataFrame({'freq_ghz': freqs_f, 's21_re': s21_f.real, 's21_im': s21_f.imag}).to_csv(
            folder / 'raw' / 'S21_fine.csv', index=False
        )

        fit = fit_notch_s21(freqs_f, s21_f)
        v_eff = 4.0 * (fit['f_r_GHz'] * 1e9) * (rep['actual_filleted_mm'] * 1e-3)
        summary = dict(
            fit=fit,
            v_eff_m_s=v_eff,
            resonator_length_mm=rep['actual_filleted_mm'],
            coupling_length_um=a.coup_len,
            coupling_ground_strip_um=a.coup_strip,
            elapsed_s=round(elapsed, 1),
            solver_temp=temp['after'],
            mesh_objects=mesh_objs,
        )
        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(
            json.dumps(summary, indent=2) + '\n', encoding='utf-8'
        )
        finish_run(folder, status='completed', summary=summary)

        print(f"\n=== HFSS Hanger S21 Result ({elapsed:.1f} s) ===")
        print(f"  f_r     = {fit['f_r_GHz']:.5f} GHz  (target {a.f_target:.3f} GHz, v_eff = {v_eff:.4e} m/s)")
        print(f"  Q_c     = {fit['Q_c']:.1f}  (Q_l = {fit['Q_l']:.1f})")
        print(f"  kappa/2pi = {fit['kappa_c_MHz']:.3f} MHz  (target 0.500 MHz)")
    except Exception as exc:
        finish_run(folder, status='failed', summary=dict(error=repr(exc)))
        raise
    finally:
        try:
            import save_aedt
            desk = save_aedt._desktop()
            desk.get_active_project()._project.Save()
            desk._desktop.CloseProject('D7')
        except Exception:
            pass


if __name__ == '__main__':
    main()
