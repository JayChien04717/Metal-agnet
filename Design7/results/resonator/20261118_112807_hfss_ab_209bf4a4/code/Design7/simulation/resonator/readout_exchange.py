"""Estimate model exchange from the two driving-point admittances at the hammer."""

from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from scipy.optimize import root

from readout_circuit import read, create_run, finish_run
from readout_ab import write


def hammer_admittances(f_ghz, p):
    omega = 2 * np.pi * f_ghz * 1e9
    y0 = 1 / p['z0_ohm']

    def tangent(length):
        return np.tan(omega * length * 1e-3 / p['velocity_m_s'])

    def line(load_y, length):
        t = tangent(length)
        return y0 * (load_y + 1j * y0 * t) / (y0 + 1j * load_y * t)

    def short(length):
        return -1j * y0 / tangent(length)

    def cap(fF):
        return 1j * omega * fF * 1e-15

    cp, cp_feed, cp_filter = map(cap, p['cp_fF'])
    cj, cj_filter, cj_res = map(cap, p['cj_fF'])
    port = 1 / p['port_impedance_ohm']
    feed = line(port, p['feed_in_mm']) + line(port, p['feed_out_mm'])
    feed_load = line(feed, p['feed_branch_mm']) + cp_feed
    cp_load = cp_filter + cp * feed_load / (cp + feed_load)
    filter_tee = short(p['filter_short_arm_mm']) + line(cp_load, p['filter_open_arm_mm'])
    res_tee = short(p['res_short_arm_mm']) + 1j * y0 * tangent(p['res_open_arm_mm'])
    return (line(filter_tee, p['filter_hammer_arm_mm']) + cj + cj_filter,
            line(res_tee, p['res_hammer_arm_mm']) + cj + cj_res, cj)


def pole(admittance, seed):
    def residual(parts):
        y = admittance(complex(*parts))
        return [y.real, y.imag]

    answer = root(residual, [seed, 0], tol=1e-10)
    if np.linalg.norm(residual(answer.x)) > 1e-10:
        raise RuntimeError(f'Could not find mode near {seed} GHz: {answer.message}')
    return complex(*answer.x)


def bare_parameters(p, seeds):
    bare = [pole(lambda f: hammer_admittances(f, p)[i], seed) for i, seed in enumerate(seeds)]
    capacitances = []
    step_ghz = .000002
    for index, mode in enumerate(bare):
        plus = hammer_admittances(mode.real + step_ghz, p)[index]
        minus = hammer_admittances(mode.real - step_ghz, p)[index]
        # C_eff is half the derivative of susceptance with angular frequency.
        capacitances.append(float((plus - minus).imag / (8 * np.pi * step_ghz * 1e9)))
    exchange = (1000 * np.sqrt(bare[0].real * bare[1].real) * p['cj_fF'][0] * 1e-15 /
                (2 * np.sqrt(np.prod(capacitances))))
    return bare, capacitances, float(exchange)


def estimate(p, seeds):
    bare, capacitances, exchange = bare_parameters(p, seeds)

    def determinant(f):
        yp, yr, cross = hammer_admittances(f, p)
        return yp * yr - cross ** 2

    bare_f = [z.real for z in bare]
    center = np.mean(bare_f)
    split = np.hypot((bare_f[1] - bare_f[0]) / 2, exchange / 1000)
    hybrid = [pole(determinant, seed) for seed in (center - split, center + split)]
    return dict(J_over_2pi_MHz=float(exchange), bare_GHz=[z.real for z in bare],
                detuning_MHz=1000 * (bare[1].real - bare[0].real),
                hammer_plane_effective_capacitance_fF=[c * 1e15 for c in capacitances],
                hybrid_GHz=[z.real for z in hybrid],
                hybrid_FWHM_kHz=[2e6 * z.imag for z in hybrid],
                method='Linearize the bare susceptances at each hammer plane; retain Cj in the diagonal load.',
                limitation='Circuit estimate with the full static mutual placed at the hammer and zero intrinsic loss.')


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('comparison_run', type=Path)
    args = parser.parse_args()
    source = args.comparison_run.resolve()
    upstream = read(source / 'run.json')
    p = upstream['parameters']['circuit']
    seeds = (5.84, 6.9)
    result = estimate(p, seeds)
    folder = create_run('Design7', 'resonator', 'q3d_exchange_estimate', entrypoint=__file__,
                        inputs=[source / 'run.json'], upstream=[str(source)],
                        parameters=dict(root_seeds_GHz=seeds, EM_rerun=False))
    manifest = read(folder / 'run.json')
    manifest['geometry_revision'] = upstream['geometry_revision']
    write(folder / 'run.json', manifest)
    finish_run(folder, status='completed', summary=result)
    print(folder)
    print(result)


if __name__ == '__main__':
    main()
