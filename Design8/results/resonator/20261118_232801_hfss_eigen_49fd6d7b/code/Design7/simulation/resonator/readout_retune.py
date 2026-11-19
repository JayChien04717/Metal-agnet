"""Size the Design7 short-leg trim using the v3 assembled bare-frequency reference."""

from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from scipy.optimize import root

from readout_ab import L, check_port_planes, ground_straps, write
from readout_circuit import COUPONS, capacitor_values, create_run, finish_run, read
from readout_exchange import bare_parameters


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('sizing_run', type=Path)
    parser.add_argument('circuit_run', type=Path)
    parser.add_argument('s21_fit_run', type=Path)
    args = parser.parse_args()
    sizing, circuit, s21 = (path.resolve() for path in
                            (args.sizing_run, args.circuit_run, args.s21_fit_run))
    geometry_file = sizing / 'raw/geometry_15mhz.json'
    sizing_summary = sizing / 'summary.json'
    circuit_manifest = circuit / 'run.json'
    s21_manifest = s21 / 'run.json'
    s21_summary = s21 / 'summary.json'
    geometry = read(geometry_file)
    base_manifest = read(circuit_manifest)
    s21_run = read(s21_manifest)
    if not geometry['geometry_revision'] == base_manifest['geometry_revision'] == s21_run['geometry_revision']:
        raise ValueError('The reference geometry and both analysis runs must use the same revision.')
    base = base_manifest['parameters']['circuit']
    fit = read(s21_summary)['fit']['inversion']
    target_f = read(sizing_summary)['prediction']['bare_GHz'][0]
    target_j = geometry['target']['J_over_2pi_MHz']
    ref_bare = np.array([fit['bare_fp_GHz'], fit['bare_fr_GHz']])
    seeds = (6.9, 6.9)
    old_height = geometry['hammer_um']['height']
    low, high = map(read, COUPONS)
    lo, hi = low['head_height_um'], high['head_height_um']
    cap_slope = (capacitor_values(high) - capacitor_values(low)) / (hi - lo)

    def model(trim_um):
        filter_trim, res_trim, height_trim = trim_um
        candidate = dict(
            base, filter_short_arm_mm=base['filter_short_arm_mm'] + filter_trim / 1000,
            res_short_arm_mm=base['res_short_arm_mm'] + res_trim / 1000,
            cj_fF=(np.array(base['cj_fF']) + height_trim * cap_slope).tolist(),
        )
        poles, _, exchange = bare_parameters(candidate, seeds)
        return np.array([pole.real for pole in poles]), exchange, candidate

    model_ref, model_ref_j, _ = model([0, 0, 0])

    def residual(trim_um):
        bare, exchange, _ = model(trim_um)
        return np.r_[1000 * (ref_bare + bare - model_ref - target_f), exchange - target_j]

    solution = root(residual, [20, -35, 0])
    height = old_height + solution.x[2]
    if not solution.success or np.linalg.norm(residual(solution.x)) > 1e-5 or not lo <= height <= hi:
        raise RuntimeError(f'Could not solve within the hammer coupon bracket: {solution.message}')
    model_bare, exchange, candidate_circuit = model(solution.x)
    base_design, old_lengths = L.build('opposing', geometry_file)
    geometry['geometry_revision'] = 'd7_readout_aligned_v4'
    geometry['channel']['filter_length_mm'] += solution.x[0] / 1000
    geometry['channel']['resonator_length_mm'] += solution.x[1] / 1000
    geometry['hammer_um']['height'] = height
    geometry['target']['bare_frequency_GHz'] = target_f
    folder = create_run(
        'Design7', 'resonator', 'retune_v3', entrypoint=__file__,
        inputs=[geometry_file, sizing_summary, circuit_manifest, s21_manifest, s21_summary, *COUPONS],
        upstream=[str(path) for path in (sizing, circuit, s21, *(path.parent for path in COUPONS))],
        parameters=dict(EM_rerun=False, geometry_revision=geometry['geometry_revision'],
                        target_bare_GHz=target_f, target_J_over_2pi_MHz=target_j,
                        model_root_seeds_GHz=seeds, circuit=candidate_circuit),
    )
    manifest = read(folder / 'run.json')
    manifest['geometry_revision'] = geometry['geometry_revision']
    write(folder / 'run.json', manifest)
    candidate_file = folder / 'raw/geometry_retune.json'
    write(candidate_file, geometry)
    design, lengths = L.build('opposing', candidate_file)
    layout_checks = L.checks(design, lengths)
    check_port_planes(design)
    straps = ground_straps(design, 'opposing')
    ends = ('filter_tail_short', 'res_tail_short')
    old_ends, new_ends = ([d.components[name].pins['end'].middle for name in ends]
                         for d in (base_design, design))
    if not np.allclose(np.subtract(new_ends, old_ends),
                       np.column_stack((np.zeros(2), solution.x[:2] / 1000)), atol=1e-8):
        raise RuntimeError('Short-end displacement differs from the specified straight-leg trim.')
    probe = np.eye(3)
    slopes = [500 * (model(step)[0][i] - model(-step)[0][i]) for i, step in enumerate(probe[:2])]
    summary = dict(
        order=['filter', 'resonator'],
        v3_assembled_bare_GHz=ref_bare, v3_assembled_conditional_J_MHz=fit['J_MHz'],
        v3_circuit_bare_GHz=model_ref, v3_circuit_J_MHz=model_ref_j,
        short_leg_trim_um=solution.x[:2], local_short_leg_slope_MHz_per_um=slopes,
        old_drawn_lengths_mm=old_lengths, new_drawn_lengths_mm=lengths,
        old_short_end_mm=old_ends, new_short_end_mm=new_ends,
        old_hammer_span_um=old_height, new_hammer_span_um=height,
        predicted_J_over_2pi_MHz=exchange, unanchored_model_bare_GHz=model_bare,
        anchored_bare_GHz=ref_bare + model_bare - model_ref, layout_checks=layout_checks,
        ground_strap_count=len(straps), port_plane_checks='passed',
        assumptions=[
            'Use the v3 assembled HFSS bare frequencies as local offsets on the circuit frequency shifts.',
            'Size exchange with the circuit. The strongly detuned v3 HFSS exchange is conditional.',
            'Apply the Design4 hammer coupon slopes to the v3 capacitances and hold the nonlocal remainder fixed.',
        ],
    )
    write(folder / 'summary.json', summary)
    finish_run(folder, status='completed', validation='layout_checked_pending_em')
    print(folder, flush=True)


if __name__ == '__main__':
    main()
