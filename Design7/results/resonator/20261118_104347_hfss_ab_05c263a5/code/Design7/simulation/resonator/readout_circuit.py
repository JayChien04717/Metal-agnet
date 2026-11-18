"""Predict the Design7 transmission using Q3D values and unfitted CPW lines."""

from argparse import ArgumentParser
import json
from pathlib import Path
import sys
import traceback

import numpy as np
import pandas as pd
from scipy.special import ellipk

HOME = Path(__file__).resolve().parents[2]
ROOT = HOME.parent
sys.path[:0] = [str(HOME / 'layout'), str(ROOT / 'code/general'),
                str(ROOT / 'Design6/simulation/resonator')]
import design7_layout as L
from readout_poles import load
from run_record import create_run, finish_run

LOCAL = ROOT / 'Design6/results/resonator/20261116_221901_readout_caps_0e896b6c'
CASE = Path('raw/opposing_full')


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def inputs(q3d_run):
    geometry = read(q3d_run / CASE / 'geometry.json')
    boundary = read(q3d_run / CASE / 'boundary.json')
    if boundary['substrate'] != L.P.d5.SI_MATERIAL:
        raise ValueError('The current CPW material differs from the recorded Q3D fixture')
    design, lengths = L.build('opposing')
    if not all(np.isclose(lengths[name], length) for name, length in
               boundary['layout']['lengths_mm'].items()):
        raise ValueError('Circuit line lengths differ from the Q3D fixture')

    def match_local(component, folder, ignore=()):
        local = read(LOCAL / 'raw' / folder / 'geometry.json')['options'][component]
        prior = {k: design.parse_value(v) for k, v in local.items() if str(v).endswith('um')}
        current = design.components[component].options
        for key, expected in prior.items():
            if key.startswith('pos_') or key in ignore:
                continue
            if not np.isclose(design.parse_value(current[key]), expected):
                raise ValueError(f'{component}.{key} differs from the selected Q3D coupon')

    match_local('cp', 'cp_coupon')
    match_local('cj', 'cj_55um', ignore=('lead_length',))
    local = read(LOCAL / 'summary_out.json')
    case = read(q3d_run / CASE / 'result.json')
    mutual = case['mutual_fF']
    cp = [mutual['feed-filter'], *local['cp']['ground_fF']]
    cj = [mutual['filter-resonator'], *local['hammer_excess_ground_fF']]
    opts = geometry['options']
    jf_y, jr_y, cj_y = [design.parse_value(opts[name]['pos_y']) for name in ('jf', 'jr', 'cj')]
    n0 = design.components['n0'].p
    ftee, rtee = (design.components[name].p for name in ('jf', 'jr'))
    pl = lambda name: L.P.d5.f3.path_length(design, name)
    w, g = (design.parse_value(L.P.CPW[key]) for key in ('trace_width', 'trace_gap'))
    effective_eps = (1 + L.P.d5.SI_EPS) / 2
    k2 = (w / (w + 2 * g)) ** 2
    return dict(cp_fF=cp, cj_fF=cj, feed_res_mutual_omitted_fF=mutual['feed-resonator'],
                source_cp_local_converged=local['cp']['solver_converged'],
                z0_ohm=float(30 * np.pi / np.sqrt(effective_eps) * ellipk(1 - k2) / ellipk(k2)),
                velocity_m_s=float(299792458 / np.sqrt(effective_eps)),
                feed_in_mm=pl('feed_in') + n0.prime_length / 2,
                feed_out_mm=pl('feed_out') + n0.prime_length / 2,
                feed_branch_mm=pl('stub_p') + n0.stub_length,
                filter_open_arm_mm=pl('feed_jf') + ftee.prime_length / 2,
                res_open_arm_mm=pl('res_open_run') + rtee.prime_length / 2,
                filter_short_arm_mm=pl('filter_tail') + ftee.stub_length,
                res_short_arm_mm=pl('res_tail') + rtee.stub_length,
                filter_hammer_arm_mm=abs(jf_y - cj_y), res_hammer_arm_mm=abs(jr_y - cj_y),
                port_impedance_ohm=50)


def s21(f_ghz, p):
    """Cascade standard lossless TL and pi sections over the as-drawn branches."""
    omega = 2 * np.pi * np.asarray(f_ghz) * 1e9
    z0, v = p['z0_ohm'], p['velocity_m_s']
    y0 = 1 / z0

    def phase(length):
        return omega / v * length * 1e-3

    def line(load_y, length):
        t = np.tan(phase(length))
        return y0 * (load_y + 1j * y0 * t) / (y0 + 1j * load_y * t)

    def shunt(c_fF):
        return 1j * omega * c_fF * 1e-15

    def pi(load_y, cap):
        mutual, ground_in, ground_out = map(shunt, cap)
        y = load_y + ground_out
        return ground_in + mutual * y / (mutual + y)

    def short(length):
        return -1j * y0 / np.tan(phase(length))

    yr = short(p['res_short_arm_mm']) + 1j * y0 * np.tan(phase(p['res_open_arm_mm']))
    yr = line(pi(line(yr, p['res_hammer_arm_mm']), p['cj_fF']), p['filter_hammer_arm_mm'])
    yf = line(yr + short(p['filter_short_arm_mm']), p['filter_open_arm_mm'])
    yb = line(pi(yf, p['cp_fF']), p['feed_branch_mm'])

    def abcd(length):
        c, s = np.cos(phase(length)), np.sin(phase(length))
        return np.array([[c, 1j * z0 * s], [1j * y0 * s, c]]).transpose(2, 0, 1)

    tee = np.broadcast_to(np.eye(2, dtype=complex), (len(omega), 2, 2)).copy()
    tee[:, 1, 0] = yb
    matrix = abcd(p['feed_in_mm']) @ tee @ abcd(p['feed_out_mm'])
    z = p['port_impedance_ohm']
    return 2 / (matrix[:, 0, 0] + matrix[:, 0, 1] / z + z * matrix[:, 1, 0] + matrix[:, 1, 1])


def compare(hfss_run, folder, params):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    f, measured = load(hfss_run / CASE / 'S_fine.csv')
    predicted = s21(f, params)
    pd.DataFrame(dict(f_GHz=f, model_S21=predicted, hfss_S21=measured)).to_csv(
        folder / 'raw/s21_overlay.csv', index=False)
    fig, axs = plt.subplots(2, 1, figsize=(9, 6), sharex=True, layout='constrained')
    for label, series in (('Q3D + TL, no fit', predicted), ('HFSS driven', measured)):
        axs[0].plot(f, 20 * np.log10(np.maximum(abs(series), 1e-12)), label=label)
        axs[1].plot(f, np.rad2deg(np.unwrap(np.angle(series))), label=label)
    axs[0].set_ylabel('|S21| (dB)')
    axs[1].set_ylabel('Phase (degrees)')
    axs[1].set_xlabel('Frequency (GHz)')
    for ax in axs:
        ax.grid(alpha=.25)
        ax.legend()
    fig.suptitle('Design7 double-U fixture')
    fig.savefig(folder / 's21_q3d_vs_hfss.png', dpi=165)
    plt.close(fig)
    return dict(model_broad_dip_GHz=float(f[np.argmin(abs(predicted))]),
                hfss_broad_dip_GHz=float(f[np.argmin(abs(measured))]),
                complex_rms=float(np.sqrt(np.mean(abs(predicted - measured) ** 2))),
                magnitude_rms=float(np.sqrt(np.mean((abs(predicted) - abs(measured)) ** 2))),
                note='Complex RMS uses the raw port reference phase, without phase or length alignment.')


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('q3d_run', type=Path)
    parser.add_argument('hfss_run', type=Path)
    a = parser.parse_args()
    q3d_run, hfss_run = a.q3d_run.resolve(), a.hfss_run.resolve()
    params = inputs(q3d_run)
    folder = create_run(
        'Design7', 'resonator', 'q3d_circuit_comparison', entrypoint=__file__,
        inputs=[q3d_run / 'run.json', q3d_run / CASE / 'result.json', q3d_run / CASE / 'geometry.json',
                q3d_run / CASE / 'boundary.json', hfss_run / 'run.json', hfss_run / CASE / 'S_fine.csv',
                LOCAL / 'summary_out.json', LOCAL / 'raw/cp_coupon/geometry.json',
                LOCAL / 'raw/cj_55um/geometry.json', L.SETTINGS],
        upstream=[str(q3d_run), str(hfss_run), str(LOCAL)],
        parameters=dict(geometry_revision='d7_readout_ab_v2', hfss_parameter_fit=False,
                        model='standard Design5/6 transmission lines and pi capacitors', circuit=params))
    print(f'RUN_FOLDER={folder}', flush=True)
    try:
        result = compare(hfss_run, folder, params)
    except Exception:
        finish_run(folder, status='failed', summary=dict(error=traceback.format_exc()))
        raise
    finish_run(folder, status='completed', summary=result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
