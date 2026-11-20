"""Audit the Design7 bare-frequency discrepancy without rerunning EM."""

from argparse import ArgumentParser
from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from readout_ab import write
from readout_circuit import CASE, create_run, finish_run, read
from readout_exchange import hammer_admittances, pole
from readout_fit import fit_spectrum


def poles_by_trace(path, picks):
    frame = pd.read_csv(path, index_col=0)
    f = frame.index.to_numpy(float)
    results = {}
    for label in ('S11', 'S21', 'S12', 'S22'):
        voltage = frame[label].map(complex).to_numpy()
        results[label], _ = fit_spectrum(f, voltage, picks)
    return results


def audit(q3d, hfss, model, fit, comparison):
    p = read(model / 'run.json')['parameters']['circuit']
    model_result = read(model / 'summary.json')
    hfss_fit = read(fit / 'summary.json')['fit']
    inv = hfss_fit['inversion']
    fitted_bare = {'filter': inv['bare_fp_GHz'], 'resonator': inv['bare_fr_GHz']}
    seeds = model_result['prediction']['bare_GHz']
    picks = [(f, width / 2) for f, width in
             zip(hfss_fit['pole_frequencies_GHz'], hfss_fit['fwhm_MHz'])]

    def bare_res(values):
        return pole(lambda f: hammer_admittances(f, values)[1], seeds[1]).real

    target = fitted_bare['resonator']
    length = p['res_short_arm_mm']
    extra_mm = brentq(lambda delta: bare_res(dict(p, res_short_arm_mm=length + delta)) - target, -.2, .2)
    extra_fF = brentq(lambda delta: bare_res(dict(p, cj_fF=[*p['cj_fF'][:2], p['cj_fF'][2] + delta]))
                      - target, -25, 25)
    velocity_scale = brentq(lambda scale: bare_res(dict(p, velocity_m_s=p['velocity_m_s'] * scale))
                            - target, .96, 1.04)
    velocity_model = dict(p, velocity_m_s=p['velocity_m_s'] * velocity_scale)
    filter_after_velocity = pole(lambda f: hammer_admittances(f, velocity_model)[0], seeds[0]).real

    c_per_mm = 1e12 / (p['z0_ohm'] * p['velocity_m_s'])
    static_ground = dict(
        feed=c_per_mm * sum(p[k] for k in ('feed_in_mm', 'feed_out_mm', 'feed_branch_mm')) + p['cp_fF'][1],
        filter=c_per_mm * sum(p[k] for k in ('filter_open_arm_mm', 'filter_short_arm_mm',
                                             'filter_hammer_arm_mm')) + p['cp_fF'][2] + p['cj_fF'][1],
        resonator=c_per_mm * sum(p[k] for k in ('res_open_arm_mm', 'res_short_arm_mm',
                                                'res_hammer_arm_mm')) + p['cj_fF'][2])
    q3d_ground = model_result['q3d']['ground_fF']
    saved_geometry = [(run / CASE / 'geometry.json').read_bytes() for run in (q3d, hfss)]
    if saved_geometry[0] != saved_geometry[1]:
        raise ValueError('Q3D and HFSS geometry snapshots differ')
    report = read(comparison / 'summary.json')['s21']['pole_fits']
    return dict(
        shared_geometry_sha256=sha256(saved_geometry[0]).hexdigest(),
        drawn_resonator_length_mm=read(hfss / CASE / 'boundary.json')['layout']['lengths_mm']['resonator'],
        circuit_resonator_cpw_route_length_mm=p['res_open_arm_mm'] + p['res_short_arm_mm'],
        circuit_exact_bare_GHz=seeds, hfss_inverted_bare_GHz=fitted_bare,
        hfss_narrow_pole_to_inverted_resonator_MHz=1000 * (target - hfss_fit['pole_frequencies_GHz'][0]),
        hfss_fit_by_trace=poles_by_trace(hfss / CASE / 'S_fine.csv', picks),
        equivalent_resonator_short_extension_um=1000 * extra_mm,
        equivalent_extra_hammer_resonator_shunt_fF=extra_fF,
        velocity_scale_to_match_resonator=velocity_scale,
        filter_bare_GHz_after_matching_velocity=filter_after_velocity,
        circuit_C_per_mm_fF=c_per_mm, circuit_static_ground_fF=static_ground,
        q3d_static_ground_fF=q3d_ground,
        q3d_minus_circuit_static_ground_fF={k: q3d_ground[k] - v for k, v in static_ground.items()},
        compared_pole_frequencies_GHz={k: v['pole_frequencies_GHz'] for k, v in report.items()},
        note='Length, shunt and velocity changes are separate equivalent perturbations, not identified causes.')


def main():
    ap = ArgumentParser(description=__doc__)
    for name in ('q3d', 'hfss', 'model', 'fit', 'comparison'):
        ap.add_argument(name, type=Path)
    args = ap.parse_args()
    runs = [getattr(args, k).resolve() for k in ('q3d', 'hfss', 'model', 'fit', 'comparison')]
    revision, = set(read(p / 'run.json')['geometry_revision'] for p in runs)
    result = audit(*runs)
    q3d, hfss, model, fit, comparison = runs
    folder = create_run(
        'Design7', 'resonator', 'bare_frequency_audit', entrypoint=__file__,
        inputs=[*(p / 'run.json' for p in runs), *(p / 'summary.json' for p in (model, fit, comparison)),
                q3d / CASE / 'geometry.json', hfss / CASE / 'geometry.json',
                hfss / CASE / 'boundary.json', hfss / CASE / 'S_fine.csv'],
        upstream=[str(p) for p in runs], parameters=dict(EM_rerun=False, geometry_revision=revision))
    manifest = read(folder / 'run.json')
    manifest['geometry_revision'] = revision
    write(folder / 'run.json', manifest)
    finish_run(folder, status='completed', summary=result)
    print(folder)
    for key, value in result.items():
        if key == 'hfss_fit_by_trace':
            value = {trace: item['inversion']['bare_fr_GHz'] for trace, item in value.items()}
        print(key, value)


if __name__ == '__main__':
    main()
