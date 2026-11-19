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
from readout_ab import PORTS
from readout_fit import fit_spectrum, peak_guesses
from readout_poles import load
from run_record import create_run, finish_run

LOCAL = ROOT / 'Design6/results/resonator/20261116_221901_readout_caps_0e896b6c'
CASE = Path('raw/opposing_full')
D4_CAPS = ROOT / 'Design4/results/resonator'
COUPONS = [
    D4_CAPS / '20261116_200156_h60_g40_w25_a55_q3d_33be5e01/summary.json',
    D4_CAPS / '20261116_200519_h80_g40_w25_a55_q3d_58545fab/summary.json',
]


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def q3d_capacitances(case_folder, pass_number):
    result = read(case_folder / 'result.json')
    names = result['net_objects']
    order = result['order']
    if pass_number is None:
        conv = pd.read_csv(case_folder / 'convergence.csv', index_col=0)
        delta_column = next(column for column in conv if 'delta' in column.lower())
        deltas = conv[delta_column].dropna()
        return dict(pass_number=int(conv.index[-1]), mutual_fF=result['mutual_fF'],
                    ground_fF=result['ground_fF'], order=order, maxwell_fF=result['maxwell_fF'],
                    adaptive_delta_percent=float(deltas.iloc[-1]) if len(deltas) else None,
                    filter_res_mutual_relative_change=result.get('last_filter_res_mutual_relative_change'))

    # Qiskit Metal converts each adaptive-pass export to farads.
    def selected_matrix(k):
        data = pd.read_csv(case_folder / f'pass_{k}.csv', index_col=0)
        labels = [names[net] for net in order]
        return data.loc[labels, labels].to_numpy(float) * 1e15

    raw = selected_matrix(pass_number)
    cm = (raw + raw.T) / 2
    if (not np.allclose(raw, raw.T, atol=.01) or np.any(np.triu(cm, 1) > 0)
            or np.any(np.linalg.eigvalsh(cm) <= 0) or np.any(cm.sum(axis=1) <= 0)):
        raise ValueError(f'Invalid Maxwell matrix on pass {pass_number}')
    conv = pd.read_csv(case_folder / 'convergence.csv', index_col=0)
    delta_column = next(column for column in conv if 'delta' in column.lower())
    previous = selected_matrix(pass_number - 1)
    drift = abs(previous[1, 2] / cm[1, 2] - 1)
    return dict(pass_number=pass_number, adaptive_delta_percent=float(conv.loc[pass_number, delta_column]),
                filter_res_mutual_relative_change=float(drift), order=order, maxwell_fF=cm.tolist(),
                ground_fF=dict(zip(order, cm.sum(axis=1))),
                mutual_fF={f'{a}-{b}': float(-cm[i, j]) for i, a in enumerate(order)
                           for j, b in enumerate(order) if j > i})


def geometry_settings(q3d_run):
    record = read(q3d_run / 'run.json')
    for entry in record['inputs']:
        snapshot = q3d_run / entry['snapshot']
        if snapshot.suffix != '.json':
            continue
        config = read(snapshot)
        if 'variants' in config and config.get('geometry_revision') == record['geometry_revision']:
            return snapshot
    if record['geometry_revision'] == read(L.SETTINGS)['geometry_revision']:
        return L.SETTINGS
    raise ValueError('Missing matching geometry configuration in the Q3D input snapshots')


def capacitor_values(coupon):
    c = coupon['capacitance']
    excess = np.asarray(c['ground_fF']) - coupon['reference']['ground_fF'][0] / 2
    return np.array([c['mutual_fF'], *excess])


def inputs(q3d_run, caps):
    geometry = read(q3d_run / CASE / 'geometry.json')
    boundary = read(q3d_run / CASE / 'boundary.json')
    if boundary['substrate'] != L.P.d5.SI_MATERIAL:
        raise ValueError('The current CPW material differs from the recorded Q3D fixture')
    settings = geometry_settings(q3d_run)
    design, lengths = L.build('opposing', settings)
    if not all(np.isclose(lengths[name], length) for name, length in
               boundary['layout']['lengths_mm'].items()):
        raise ValueError('Circuit line lengths differ from the Q3D fixture')

    def match_local(component, path, ignore=()):
        local = read(path)['options'][component]
        prior = {k: design.parse_value(v) for k, v in local.items() if str(v).endswith('um')}
        current = design.components[component].options
        for key, expected in prior.items():
            if key.startswith('pos_') or key in ignore:
                continue
            if not np.isclose(design.parse_value(current[key]), expected):
                raise ValueError(f'{component}.{key} differs from the selected Q3D coupon')

    local_record = LOCAL / 'summary_out.json'
    cp_geometry = LOCAL / 'raw/cp_coupon/geometry.json'
    cj_geometry = LOCAL / 'raw/cj_55um/geometry.json'
    match_local('cp', cp_geometry)
    local = read(local_record)
    sources = [settings, local_record, cp_geometry]
    height = design.components['cj'].p.head_height
    previous_height = design.parse_value(read(cj_geometry)['options']['cj']['head_height'])
    if np.isclose(height, previous_height):
        match_local('cj', cj_geometry, ignore=('lead_length',))
        cj_shunts = local['hammer_excess_ground_fF']
        sources.append(cj_geometry)
    else:
        coupons = [read(path) for path in COUPONS]
        lo, hi = [item['head_height_um'] / 1000 for item in coupons]
        if not lo <= height <= hi:
            raise ValueError('Hammer span is outside the reference coupon bracket')
        for path in COUPONS:
            shape = path.parent / 'raw/cj_55um/geometry.json'
            match_local('cj', shape, ignore=('lead_length', 'head_height'))
            sources.extend((path, shape))
        a, b = [capacitor_values(item)[1:] for item in coupons]
        cj_shunts = (a + (height - lo) / (hi - lo) * (b - a)).tolist()
    mutual = caps['mutual_fF']
    cp = [mutual['feed-filter'], *local['cp']['ground_fF']]
    cj = [mutual['filter-resonator'], *cj_shunts]
    opts = geometry['options']
    jf_y, jr_y, cj_y = [design.parse_value(opts[name]['pos_y']) for name in ('jf', 'jr', 'cj')]
    n0 = design.components['n0'].p
    ftee, rtee = (design.components[name].p for name in ('jf', 'jr'))
    pl = lambda name: L.P.d5.f3.path_length(design, name)
    w, g = (design.parse_value(L.P.CPW[key]) for key in ('trace_width', 'trace_gap'))
    effective_eps = (1 + L.P.d5.SI_EPS) / 2
    k2 = (w / (w + 2 * g)) ** 2
    normals = [design.components[component].pins[pin].normal for component, pin in PORTS]
    # HFSS's positive-axis port integration lines use opposite relative signs
    # at the two outward normals. This polarity follows the port definitions.
    hfss_polarity = int(np.prod([np.sign(n[np.argmax(abs(n))]) for n in normals]))
    return dict(cp_fF=cp, cj_fF=cj, feed_res_mutual_omitted_fF=mutual['feed-resonator'],
                source_cp_local_converged=local['cp']['solver_converged'],
                hfss_s21_polarity=hfss_polarity,
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
                port_impedance_ohm=50), sources


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

    f, hfss_raw = load(hfss_run / CASE / 'S_fine.csv')
    measured = params['hfss_s21_polarity'] * hfss_raw
    predicted = s21(f, params)
    pd.DataFrame(dict(f_GHz=f, model_S21=predicted, hfss_S21_raw=hfss_raw,
                      hfss_S21_signal_to_ground=measured)).to_csv(
        folder / 'raw/s21_overlay.csv', index=False)
    fig, axs = plt.subplots(3, 1, figsize=(9, 8), layout='constrained')
    measured_dip = float(f[np.argmin(abs(measured))])
    model_dip = float(f[np.argmin(abs(predicted))])
    dips = [model_dip, measured_dip]
    pole_fits = {}
    for label, series in (('model', predicted), ('hfss', measured)):
        guesses = peak_guesses(f, series)
        if guesses:
            dips.extend(center for center, _ in guesses)
            pole_fits[label], _ = fit_spectrum(f, series, guesses)
    for label, series in (('Q3D + TL, no fit', predicted), ('HFSS driven', measured)):
        for ax in axs[:2]:
            ax.plot(f, 20 * np.log10(np.maximum(abs(series), 1e-12)), label=label)
            ax.set_ylabel('|S21| (dB)')
        axs[2].plot(f, np.rad2deg(np.unwrap(np.angle(series))), label=label)
    axs[1].set_xlim(max(f[0], min(dips) - .035), min(f[-1], max(dips) + .035))
    axs[2].set_ylabel('Phase (degrees)')
    for ax in axs:
        ax.grid(alpha=.25)
        ax.legend()
        ax.set_xlabel('Frequency (GHz)')
    fig.suptitle('Design7 double-U fixture')
    fig.savefig(folder / 's21_q3d_vs_hfss.png', dpi=165)
    plt.close(fig)
    broad = {}
    for label, series, center in (('model', predicted, model_dip), ('hfss', measured, measured_dip)):
        band = abs(f - center) < .075
        guesses = peak_guesses(f[band], series[band])
        if len(guesses) == 1:
            fit, _ = fit_spectrum(f[band], series[band], guesses)
            broad[label] = dict(frequency_GHz=fit['pole_frequencies_GHz'][0], fwhm_MHz=fit['fwhm_MHz'][0])
    return dict(model_deepest_dip_GHz=model_dip, hfss_deepest_dip_GHz=measured_dip,
                model_minus_hfss_dip_MHz=1000 * (model_dip - measured_dip), broad_pole=broad,
                pole_fits=pole_fits, frequency_band_GHz=[float(f[0]), float(f[-1])],
                complex_rms=float(np.sqrt(np.mean(abs(predicted - measured) ** 2))),
                magnitude_rms=float(np.sqrt(np.mean((abs(predicted) - abs(measured)) ** 2))),
                note='Complex RMS uses port-sign conversion from the outward normals. No phase or length is fitted.')


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('q3d_run', type=Path)
    parser.add_argument('hfss_run', type=Path, nargs='?', help='Omit to freeze a prediction before HFSS')
    parser.add_argument('--q3d-pass', type=int, help='Use an exported adaptive pass in place of LastAdaptive')
    parser.add_argument('--frozen', type=Path, help='Reuse a pre-HFSS prediction run without rebuilding its inputs')
    a = parser.parse_args()
    q3d_run = a.q3d_run.resolve()
    hfss_run = a.hfss_run.resolve() if a.hfss_run else None
    if a.frozen and (hfss_run is None or a.q3d_pass is not None):
        parser.error('--frozen requires an HFSS run and uses the pass recorded in the prediction')
    revision = read(q3d_run / 'run.json')['geometry_revision']
    case_folder = q3d_run / CASE
    inputs_to_save = [q3d_run / 'run.json', case_folder / 'result.json', case_folder / 'geometry.json',
                      case_folder / 'boundary.json', case_folder / 'convergence.csv']
    upstream = [str(q3d_run)]
    if a.frozen:
        frozen = a.frozen.resolve()
        record = read(frozen / 'run.json')
        if str(q3d_run) not in record['upstream'] or record['geometry_revision'] != revision:
            raise ValueError('Frozen prediction uses a different Q3D result or geometry revision')
        params = record['parameters']['circuit']
        caps = read(frozen / 'summary.json')['q3d']
        inputs_to_save += [frozen / 'run.json', frozen / 'summary.json']
        upstream.append(str(frozen))
    else:
        caps = q3d_capacitances(case_folder, a.q3d_pass)
        params, sources = inputs(q3d_run, caps)
        inputs_to_save.extend(sources)
        upstream.append(str(LOCAL))
        if a.q3d_pass is not None:
            inputs_to_save.extend(case_folder / f'pass_{k}.csv' for k in (a.q3d_pass - 1, a.q3d_pass))
    if hfss_run is not None:
        inputs_to_save.extend(check_pair(q3d_run, hfss_run, revision))
        upstream.append(str(hfss_run))
    folder = create_run(
        'Design7', 'resonator', 'q3d_circuit_comparison' if hfss_run else 'q3d_circuit_prediction',
        entrypoint=__file__, inputs=inputs_to_save, upstream=upstream,
        parameters=dict(geometry_revision=revision, hfss_parameter_fit=False,
                        model='standard Design5/6 transmission lines and pi capacitors', circuit=params,
                        q3d_pass=caps['pass_number']))
    run_path = folder / 'run.json'
    record = read(run_path)
    record['geometry_revision'] = revision
    run_path.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(f'RUN_FOLDER={folder}', flush=True)
    try:
        if hfss_run:
            result = dict(q3d=caps, s21=compare(hfss_run, folder, params))
        else:
            from readout_exchange import estimate

            target = read(geometry_settings(q3d_run)).get('target', {}).get('bare_frequency_GHz')
            if target is not None:
                seeds = (target, target)
            else:
                seeds = (6.9, 6.9) if revision == 'd7_readout_15mhz_v3' else (5.84, 6.9)
            result = dict(q3d=caps, prediction=estimate(params, seeds))
    except Exception:
        finish_run(folder, status='failed', summary=dict(error=traceback.format_exc()))
        raise
    finish_run(folder, status='completed', summary=result)
    print(json.dumps(result, indent=2), flush=True)


def check_pair(q3d_run, hfss_run, revision):
    q3d_boundary, hfss_boundary = [read(run / CASE / 'boundary.json') for run in (q3d_run, hfss_run)]
    if q3d_boundary['domain'] != hfss_boundary['domain']:
        raise ValueError('Q3D and HFSS ground domains differ')
    hfss_geometry = hfss_run / CASE / 'geometry.json'
    if read(q3d_run / CASE / 'geometry.json') != read(hfss_geometry):
        raise ValueError('Q3D and HFSS signal geometry snapshots differ')
    if revision != read(hfss_run / 'run.json')['geometry_revision']:
        raise ValueError('Q3D and HFSS geometry revisions differ')
    return [hfss_run / 'run.json', hfss_run / CASE / 'S_fine.csv', hfss_run / CASE / 'boundary.json',
            hfss_geometry]


if __name__ == '__main__':
    main()
