"""Size the opposing 15 MHz draft from the recorded circuit and hammer coupons."""

from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from scipy.optimize import root

from readout_ab import L, write
from readout_circuit import COUPONS, LOCAL, capacitor_values, create_run, finish_run, read
from readout_exchange import bare_parameters, estimate

def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('comparison_run', type=Path)
    args = ap.parse_args()
    source = args.comparison_run.resolve() / 'run.json'
    original = read(source)['parameters']['circuit']
    geometry = read(L.SETTINGS)
    local = LOCAL / 'summary_out.json'
    remainder = original['cj_fF'][0] - read(local)['hammer']['mutual_fF']
    first, last = [read(path) for path in COUPONS]
    lo, hi = first['head_height_um'], last['head_height_um']
    low_caps, high_caps = map(capacitor_values, (first, last))
    target = 15.0
    seeds = (6.9, 6.9)

    def candidate_circuit(length, height):
        caps = low_caps + (height - lo) / (hi - lo) * (high_caps - low_caps)
        caps[0] += remainder
        return dict(original, filter_short_arm_mm=length, cj_fF=caps.tolist())

    def residual(candidate):
        bare, _, exchange = bare_parameters(candidate_circuit(*candidate), seeds)
        return [1000 * (bare[1].real - bare[0].real), exchange - target]

    fit = root(residual, [2.81, 74])
    if not fit.success or np.linalg.norm(residual(fit.x)) > 1e-4 or not lo <= fit.x[1] <= hi:
        raise RuntimeError(f'Could not size the candidate within the coupon bracket: {fit.message}')
    length, height = map(float, fit.x)
    delta = length - original['filter_short_arm_mm']
    geometry['geometry_revision'] = 'd7_readout_15mhz_v3'
    geometry['status'] = 'draft_pending_em'
    geometry['variants'] = {'opposing': geometry['variants']['opposing']}
    geometry['channel']['filter_length_mm'] += delta
    geometry['hammer_um']['height'] = height
    geometry['target'] = dict(J_over_2pi_MHz=target, bare_detuning_MHz=0)
    prediction = estimate(candidate_circuit(length, height), seeds)

    folder = create_run(
        'Design7', 'resonator', 'size_15mhz_candidate', entrypoint=__file__,
        inputs=[L.SETTINGS, source, local, *COUPONS], upstream=[str(source.parent), *(str(p.parent) for p in COUPONS)],
        parameters=dict(EM_rerun=False, target_J_over_2pi_MHz=target, target_bare_detuning_MHz=0,
                        solver_root_seeds_GHz=seeds, assumed_full_minus_local_mutual_fF=remainder))
    record = read(folder / 'run.json')
    record['geometry_revision'] = geometry['geometry_revision']
    write(folder / 'run.json', record)
    candidate = folder / 'raw/geometry_15mhz.json'
    write(candidate, geometry)
    design, measured = L.build('opposing', candidate)
    geometry_checks = L.checks(design, measured)
    endpoint = design.components['filter_tail_short'].pins['end'].middle.tolist()
    summary = dict(prediction=prediction, layout_checks=geometry_checks, hammer_height_um=height,
                   filter_short_arm_mm=length, filter_short_trim_um=1000 * delta,
                   filter_short_end_mm=endpoint, assumed_full_mutual_fF=candidate_circuit(length, height)['cj_fF'][0],
                   status='model_sized_geometry_pending_q3d_and_hfss')
    finish_run(folder, status='completed', summary=summary)
    print(folder, flush=True)
    print(summary, flush=True)


if __name__ == '__main__':
    main()
