# -*- coding: utf-8 -*-
"""design6_reflection_purcell_resonator Q3D: extract the 3D capacitance matrix of
the hexagon-tapered InterdigitalCap (`IDC_in`) and the parallel Purcell-bus /
readout-resonator coupled CPW section.

Solver settings (skill/em_simulation.md section 2):
    Q3D max_passes = 8, min_passes = 2, percent_error = 0.3, freq_ghz = 6.65.

Cases:
    idc   -- InterdigitalCap (`IDC_in`) coupon alone with open pins a (filter side) and b (feed side).
             Extracts series C_in = -C_ab, filter-side shunt C_a_gnd, and feed-side shunt C_b_gnd.
    coup  -- Parallel coupled CPW coupon (length = coupling_length_um = 380 um, 10/6 um CPW,
             6 um ground strip -> 28 um center-to-center spacing) with 200 um perpendicular leads.
    ref   -- Single isolated CPW line of the same total length for line-capacitance reference.

Usage (from repository root, qmetal Python):
    python -X utf8 -I -B design6_reflection_purcell_resonator/simulation/resonator/q3d_purcell_idc.py
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
from qiskit_metal import Dict, designs, draw  # noqa: E402
from qiskit_metal.qlibrary.core import QComponent  # noqa: E402

from interdigital_cap import InterdigitalCap, alley_cap_fF  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

MAX_PASSES = 8
PCT_ERROR = 0.3


class CoupledCpwCoupon(QComponent):
    """Two parallel horizontal CPW lines of length ``length`` separated by ``ground_strip``
    (center-to-center spacing = trace_width + 2*trace_gap + ground_strip), with etched open ends.
    """

    component_metadata = Dict(short_name='ccpw', _qgeometry_table_poly='True')
    default_options = Dict(
        length='380um',
        trace_width='10um',
        trace_gap='6um',
        ground_strip='6um',
        layer='1',
    )

    def make(self):
        p = self.p
        w, g, s, l_len = p.trace_width, p.trace_gap, p.ground_strip, p.length
        dy = (w + 2.0 * g + s) / 2.0
        for nm, y in (('bus_line', +dy), ('res_line', -dy)):
            metal = draw.rectangle(l_len, w, 0.0, y)
            pocket = draw.rectangle(l_len + 2.0 * g, w + 2.0 * g, 0.0, y)
            self.add_qgeometry('poly', {nm: metal}, layer=p.layer)
            self.add_qgeometry('poly', {f'{nm}_cut': pocket}, layer=p.layer, subtract=True)


def load_geo():
    return json.loads((DESIGN / 'parameters' / 'geometry.json').read_text(encoding='utf-8'))


def _write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, default=float) + '\n', encoding='utf-8')


def new_design(size_x_mm, size_y_mm, material='silicon'):
    d = designs.DesignPlanar(overwrite_enabled=True)
    d.chips.main.material = material
    d.chips.main.size.update(
        size_x=f'{size_x_mm}mm',
        size_y=f'{size_y_mm}mm',
        size_z='-500um',
        center_x='0mm',
        center_y='0mm',
    )
    return d


def case_design(case, geo):
    cpw = geo['cpw']
    idc = geo['interdigital_cap']
    rr = geo['readout_resonators']
    tw = f"{cpw['trace_width_um']:g}um"
    tg = f"{cpw['trace_gap_um']:g}um"

    if case == 'idc':
        d = new_design(1.2, 0.8, geo['chip']['material'])
        InterdigitalCap(d, 'IDC_in', options=Dict(
            pos_x='0mm',
            pos_y='0mm',
            orientation='0',
            finger_count=str(idc['finger_count']),
            finger_width=f"{idc['finger_width_um']:g}um",
            finger_gap=f"{idc['finger_gap_um']:g}um",
            overlap=f"{idc['overlap_um']:g}um",
            end_gap=f"{idc['end_gap_um']:g}um",
            bus_width=f"{idc['bus_width_um']:g}um",
            taper_length=f"{idc['taper_length_um']:g}um",
            lead_length=f"{idc['lead_length_um']:g}um",
            ground_gap=f"{idc['ground_gap_um']:g}um",
            trace_width=tw,
            trace_gap=tg,
            open_pins='True',
        ))
        return d, ['IDC_in'], [], ('comb_a', 'comb_b')

    if case == 'coup':
        d = new_design(1.0, 0.6, geo['chip']['material'])
        CoupledCpwCoupon(d, 'coup', options=Dict(
            length=f"{rr['coupling_length_um']:g}um",
            trace_width=tw,
            trace_gap=tg,
            ground_strip=f"{rr['coupling_ground_strip_um']:g}um",
        ))
        return d, ['coup'], [], ('bus_line', 'res_line')

    raise ValueError(f'Unknown Q3D case: {case}')


def solve_case(case, geo, folder, run_tag):
    from qiskit_metal.analyses.quantization import LOManalysis

    d, comps, opens, keys = case_design(case, geo)
    out = folder / 'raw' / case
    out.mkdir(parents=True, exist_ok=False)

    lom = LOManalysis(d, 'q3d')
    lom.sim.renderer.options['x_buffer_width_mm'] = 0.0
    lom.sim.renderer.options['y_buffer_width_mm'] = 0.0
    lom.sim.setup.update(
        max_passes=MAX_PASSES,
        min_passes=2,
        percent_error=PCT_ERROR,
        freq_ghz=geo['purcell_filter']['f_target_GHz'],
        solution_order='Medium',
        auto_increase_solution_order=False,
    )
    t0 = time.time()
    lom.sim.run(
        name=f'D6R_Q3D_{case}_{run_tag}',
        components=comps,
        open_terminations=opens,
        box_plus_buffer=False,
    )
    cm = lom.sim.capacitance_matrix
    cm.to_csv(out / 'matrix.csv')

    names = [n for n in cm.index if 'ground' not in n.lower()]
    signals = []
    for key in keys:
        hit = [n for n in names if any(k in n.lower() for k in key.split('|'))]
        if len(hit) != 1:
            raise RuntimeError(f'{case}: net for {key!r} not unique in {names}')
        signals.append(hit[0])

    sub = cm.loc[signals, signals].to_numpy(float)
    passes = []
    for k, m in lom.sim.capacitance_all_passes.items():
        frame = m if isinstance(m, pd.DataFrame) else pd.DataFrame(m, index=cm.index, columns=cm.columns)
        frame.to_csv(out / f'pass_{k}.csv')
        passes.append(frame.loc[signals, signals].to_numpy(float))

    delta = None
    if len(passes) > 1:
        delta = float(np.max(np.abs(passes[-1] - passes[-2]) / np.maximum(np.abs(passes[-1]), 1e-30)))

    short_keys = tuple(k.split('|')[0] for k in keys)
    res = dict(
        case=case,
        signals=dict(zip(short_keys, signals)),
        mutual_fF=float(-sub[0, 1]),
        ground_fF=dict(zip(short_keys, sub.sum(axis=1).tolist())),
        self_fF=dict(zip(short_keys, np.diag(sub).tolist())),
        maxwell_fF=sub.tolist(),
        max_passes=MAX_PASSES,
        percent_error=PCT_ERROR,
        passes_saved=len(passes),
        last_pass_max_relative_change=delta,
        solver_converged=bool(lom.sim.is_converged),
        minutes=(time.time() - t0) / 60.0,
        aedt_design=lom.sim.renderer.pinfo.design.name,
    )
    _write(out / 'caps.json', res)
    lom.sim.renderer.pinfo.project._project.Save()
    print(json.dumps(res), flush=True)
    return res


def derive_metrics(results, geo):
    cpw = geo['cpw']
    idc = geo['interdigital_cap']
    pf = geo['purcell_filter']
    rr = geo['readout_resonators']

    z0 = cpw['assumed_z0_ohm']
    v_ph = cpw['assumed_v_ph_m_s']
    c_per_m = cpw['assumed_c_per_m_f']
    c_line_ff_per_mm = c_per_m * 1e12  # 163.4 fF/mm
    f_p = pf['f_target_GHz']
    omega_p = 2.0 * math.pi * f_p * 1e9

    out = {}
    if 'idc' in results:
        r_idc = results['idc']
        c_in_q3d = r_idc['mutual_fF']
        # De-embed the short 45 um CPW lead included on each side of IDC_in
        lead_mm = idc['lead_length_um'] / 1000.0
        c_lead_ff = c_line_ff_per_mm * lead_mm
        c_a_gnd_raw = r_idc['ground_fF']['comb_a']
        c_b_gnd_raw = r_idc['ground_fF']['comb_b']
        c_a_gnd_excess = c_a_gnd_raw - c_lead_ff
        c_b_gnd_excess = c_b_gnd_raw - c_lead_ff

        c_alley = alley_cap_fF(
            idc['finger_count'],
            idc['overlap_um'],
            finger_width_um=idc['finger_width_um'],
            h_um=abs(geo['chip']['size_z_um']),
            eps_r=11.9,
        )
        # Single-port reflection lambda/2 Purcell filter loaded at input by (C_in, C_a_gnd) + 50 Ohm port
        w_cin_z0 = omega_p * (c_in_q3d * 1e-15) * z0
        c_in_eff_shunt_ff = c_in_q3d / (1.0 + w_cin_z0 ** 2)
        c_total_pf_end_ff = c_a_gnd_excess + c_in_eff_shunt_ff
        kappa_p_mhz = (2.0 * (f_p * 1e3) * (w_cin_z0 ** 2)) / (math.pi * (1.0 + w_cin_z0 ** 2))
        q_p = (f_p * 1e3) / kappa_p_mhz

        out['idc'] = dict(
            c_in_q3d_fF=c_in_q3d,
            c_in_alley_fF=c_alley,
            q3d_over_alley_ratio=c_in_q3d / c_alley,
            c_a_gnd_raw_fF=c_a_gnd_raw,
            c_b_gnd_raw_fF=c_b_gnd_raw,
            c_a_gnd_deembedded_fF=c_a_gnd_excess,
            c_b_gnd_deembedded_fF=c_b_gnd_excess,
            c_in_eff_shunt_fF=c_in_eff_shunt_ff,
            c_total_pf_end_load_fF=c_total_pf_end_ff,
            pf_end_load_shorten_mm=(c_total_pf_end_ff * 1e-15 / c_per_m) * 1e3,
            predicted_kappa_p_MHz=kappa_p_mhz,
            predicted_Q_p=q_p,
        )

    if 'coup' in results:
        r_c = results['coup']
        l_coup_mm = rr['coupling_length_um'] / 1000.0
        c_m_ff = r_c['mutual_fF']
        c11_ff = r_c['self_fF']['bus_line']
        c22_ff = r_c['self_fF']['res_line']
        k_c = c_m_ff / math.sqrt(c11_ff * c22_ff)
        # Distributed coupled-line coupling rate J between lambda/2 filter and lambda/4 resonator
        # With k_L approx k_C on homogeneous substrate, effective coupled-section J:
        theta_c = omega_p * (l_coup_mm * 1e-3) / v_ph
        j_est_mhz = (f_p * 1e3) * (math.sqrt(2.0) / math.pi) * k_c * math.sin(theta_c)
        out['coupled_section'] = dict(
            coupling_length_um=rr['coupling_length_um'],
            ground_strip_um=rr['coupling_ground_strip_um'],
            c_mutual_fF=c_m_ff,
            c_mutual_per_mm_fF=c_m_ff / l_coup_mm,
            c_bus_self_fF=c11_ff,
            c_res_self_fF=c22_ff,
            k_c=k_c,
            estimated_J_MHz=j_est_mhz,
        )

    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cases', default='idc,coup')
    args = ap.parse_args()
    cases = [c.strip() for c in args.cases.split(',') if c.strip()]

    geo = load_geo()
    params = dict(
        cases=cases,
        max_passes=MAX_PASSES,
        percent_error=PCT_ERROR,
        geometry_revision=geo['geometry_revision'],
        interdigital_cap=geo['interdigital_cap'],
        readout_coupling=dict(
            coupling_length_um=geo['readout_resonators']['coupling_length_um'],
            coupling_ground_strip_um=geo['readout_resonators']['coupling_ground_strip_um'],
        ),
    )
    folder = create_run(
        DESIGN.name,
        'resonator',
        'q3d_purcell_idc',
        entrypoint=__file__,
        inputs=[DESIGN / 'parameters' / 'geometry.json'],
        parameters=params,
    )
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        import hfss_dirs
        import save_aedt

        hfss_dirs.use_d_drive()
        save_aedt.use_project(str(DESIGN / 'simulation' / 'aedt' / 'd6r.aedt'))
        rec = json.loads((folder / 'run.json').read_text(encoding='utf-8'))
        rec['environment']['solver_version'] = save_aedt._desktop().get_version()
        rec['geometry_revision'] = geo['geometry_revision']
        _write(folder / 'run.json', rec)

        results = {}
        for c in cases:
            results[c] = solve_case(c, geo, folder, folder.name[-8:])

        derived = derive_metrics(results, geo)
        summary = dict(
            results=results,
            derived=derived,
            validation='Q3D 8-pass budget; solver_converged and last-pass relative change recorded per case.',
        )
        print(json.dumps(derived, indent=2), flush=True)
        finish_run(folder, status='completed', summary=summary)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
