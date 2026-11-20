"""Compare matched resonator-only eigenmodes and probe the endpoint field phase."""

from argparse import ArgumentParser
from pathlib import Path
import traceback

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from readout_ab import HOME, create_run, finish_run, save_aedt, write
from readout_circuit import read
from readout_eigen_audit import CASES, circuit_reference


def read_cases(source, params):
    results, inputs = {}, []
    length = params['res_open_arm_mm'] + params['res_short_arm_mm']
    for case in CASES:
        folder = source / 'raw' / case
        result = folder / 'result.json'
        by_mode, overall = folder / 'convergence_f.csv', folder / 'convergence.csv'
        if result.is_file() and 'resonator_GHz' in read(result):
            results[case] = read(result)
            inputs.append(result)
        elif by_mode.is_file() and overall.is_file():
            # The lossless-mode exporter may return None for linewidth. The
            # adaptive frequency table already contains the solved frequencies.
            modes = pd.read_csv(by_mode, index_col=0)
            conv = pd.read_csv(overall, index_col=0)
            f = modes.iloc[-1].to_numpy(float)
            index, = np.flatnonzero((6.5 < f) & (f < 7.5))
            deltas = conv.iloc[:, -1].dropna()
            ref = circuit_reference(case, params, length)
            evidence = dict(aedt_design=f'D7_{source.name[-8:]}_{case}_hfss', aedt_setup='Setup',
                            modes_GHz=f.tolist(), mode_index_1based=int(index) + 1, resonator_GHz=float(f[index]),
                            circuit_GHz=ref, length_mm=length, adaptive_passes=len(conv),
                            last_delta_f_percent=float(deltas.iloc[-1]) if len(deltas) else None,
                            hfss_minus_circuit_MHz=float(1000 * (f[index] - ref)),
                            frequency_source='Last row of saved AEDT frequency-vs-pass report')
            results[case] = evidence
            inputs.extend((by_mode, overall))
        if case in results and by_mode.is_file():
            f_mode = pd.read_csv(by_mode, index_col=0).iloc[:, results[case]['mode_index_1based'] - 1]
            if len(f_mode) > 1:
                results[case]['resonator_last_step_MHz'] = float(1000 * (f_mode.iloc[-1] - f_mode.iloc[-2]))
    return results, inputs


def transitions(results):
    ordered = {
        'counterelectrode_grounding': ('isolated', 'grounded_hammer'),
        'hammer_replacing_open_stub': ('plain_stub', 'grounded_hammer'),
        'tee_with_open_stub': ('no_tee', 'plain_stub'),
        'open_gap_9_to_45_um': ('no_tee', 'open_gap'),
        'short_contact_flare': ('no_tee', 'short_flare'),
    }
    out = {}
    for label, (reference, changed) in ordered.items():
        if not {reference, changed} <= results.keys():
            continue
        a, b = results[reference], results[changed]
        hfss = 1000 * (b['resonator_GHz'] - a['resonator_GHz'])
        ideal = 1000 * (b['circuit_GHz'] - a['circuit_GHz'])
        steps = [item.get('resonator_last_step_MHz') for item in (a, b)]
        out[label] = dict(reference=reference, changed=changed, hfss_shift_MHz=hfss,
                          circuit_shift_MHz=ideal, residual_shift_MHz=hfss - ideal,
                          sum_absolute_last_mode_steps_MHz=(
                              sum(map(abs, steps)) if all(step is not None for step in steps) else None))
    return out


def plot_comparison(source, folder, results, differences):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), layout='constrained')
    for case, result in results.items():
        path = source / 'raw' / case / 'convergence_f.csv'
        if not path.is_file():
            continue
        series = pd.read_csv(path, index_col=0).iloc[:, result['mode_index_1based'] - 1]
        residual = 1000 * (series - result['circuit_GHz'])
        axes[0].plot(series.index, residual, 'o-', label=case.replace('_', ' '))
    axes[0].set_xlabel('HFSS adaptive pass')
    axes[0].set_ylabel('HFSS − circuit (MHz)')
    axes[0].legend(ncol=2)
    labels = ['Ground counterelectrode', 'Hammer vs plain stub', 'Tee and open stub',
              'Open gap 9 → 45 µm', 'Short contact flare']
    all_keys = ('counterelectrode_grounding', 'hammer_replacing_open_stub', 'tee_with_open_stub',
                'open_gap_9_to_45_um', 'short_contact_flare')
    data = [(label, differences[key]['residual_shift_MHz'])
            for label, key in zip(labels, all_keys) if key in differences]
    if data:
        names, values = zip(*data)
        axes[1].barh(names, values, color=['#306f97' if value < 0 else '#bc5a35' for value in values])
        axes[1].invert_yaxis()
    axes[1].set_xlabel('Change in HFSS frequency minus change in circuit frequency (MHz)')
    for ax, axis in zip(axes, ('y', 'x')):
        ax.grid(axis=axis, alpha=.2)
    axes[0].axhline(0, color='0.5', linewidth=.7)
    axes[1].axvline(0, color='0.5', linewidth=.7)
    fig.suptitle('Design7 resonator-only eigenmode audit')
    fig.savefig(folder / 'eigenmode_comparison.png', dpi=165)
    plt.close(fig)


def export_gap_strip(reporter, solution, out, x, y, sign, side, z):
    low, high = sorted((y + sign * .075, y + sign * .825))
    left, right = sorted((x + side * .0085, x + side * .0165))
    start, stop = [f'{left}mm', f'{low}mm', f'{z}mm'], [f'{right}mm', f'{high}mm', f'{z}mm']
    reporter.CalcStack('clear')
    reporter.EnterQty('E')
    reporter.CalcOp('ScalarX')
    reporter.CalcOp('CmplxMag')
    reporter.ExportOnGrid(out.as_posix(), start, stop, ['0.0005mm', '0.01mm', '0.01mm'], solution,
                          ['Phase:=', '0deg'], True, 'Cartesian', ['0mm'] * 3, False)
    rows = []
    for line in out.read_text(encoding='utf-8', errors='replace').splitlines():
        if len(line.split()) != 4:
            continue
        try:
            rows.append([float(v) for v in line.split()])
        except ValueError:
            continue
    a = np.asarray(rows)
    if len(a) < 40:
        raise ValueError(f'Too few field samples in {out}')
    # AEDT releases can export points in m or in the model's mm units.
    grid_y = a[:, 1] * (1000 if np.ptp(a[:, 1]) < .005 else 1)
    d = (grid_y - y) * sign
    keep = np.isfinite(a[:, 3]) & (d >= .07) & (d <= .83)
    # Each uniform-width transverse strip gives a voltage proxy at its y coordinate.
    strip = pd.DataFrame({'d': np.round(d[keep], 7), 'ex': a[keep, 3]})
    mean = strip.groupby('d')['ex'].mean()
    return mean.index.to_numpy(float), mean.to_numpy(float)


def terminal_fields(source, result, folder):
    from pyEPR.ansys import parse_units

    desk = save_aedt._desktop()
    project, = [p for p in desk.get_projects() if p.name == HOME.name]
    design = project.get_design(result['aedt_design'])
    setup = design.get_setup(result['aedt_setup'])
    setup.get_solutions().set_mode(result['mode_index_1based'])
    reporter = design._fields_calc
    options = read(source / 'raw/no_tee/geometry.json')['options']
    straps = read(source / 'raw/no_tee/boundary.json')['straps']
    traces = []
    for name, component, sign in (('open', 'readout_open', 1), ('short', 'res_tail_short', -1)):
        coords = options[component]
        x, y = [float(parse_units(coords[k])) * 1000 for k in ('pos_x', 'pos_y')]
        for side in (-1, 1):
            path = folder / 'raw' / f'{name}_{side:+d}_field.fld'
            distance, field = export_gap_strip(reporter, str(setup.solution_name), path, x, y, sign, side, .002)
            keep = np.ones_like(distance, dtype=bool)
            for strap in straps:
                sx, sy = strap['pos']
                if abs(sx - x) < .005 and abs(strap['ori'][0]) > .9:
                    keep &= abs(distance - (sy - y) * sign) > .025
            distance, field = distance[keep], field[keep]
            if len(field) < 40 or np.max(field) <= 0:
                raise ValueError(f'Empty modal field strip at {name} end')
            traces.append((name, distance, field / np.max(field)))
    return fit_terminal_phase(traces, result)


def fit_terminal_phase(traces, result):
    fits = []
    for minimum in (.09, .15, .21):
        selected = [(k, d[d >= minimum], v[d >= minimum]) for k, d, v in traces]

        def residual(parameters):
            wave, open_pad, short_pad, *amplitudes = parameters
            values = []
            for (kind, d, observed), amplitude in zip(selected, amplitudes):
                curve = np.cos(wave * (d + open_pad)) if kind == 'open' else np.sin(wave * (d + short_pad))
                values.extend(amplitude * curve - observed)
            return values

        initial = [.375, .010, .005, 1, 1, 3, 3]
        answer = least_squares(residual, initial,
                               bounds=([.29, -.10, -.10, 0, 0, 0, 0], [.48, .10, .10, 20, 20, 20, 20]))
        wave, open_pad, short_pad, *_ = answer.x
        freq = result['resonator_GHz']
        velocity = 2 * np.pi * freq * 1e6 / wave
        fits.append(dict(minimum_distance_mm=minimum, wave_per_mm=float(wave),
                         velocity_m_s=float(velocity), open_extension_um=float(1000 * open_pad),
                         short_extension_um=float(1000 * short_pad),
                         open_shift_MHz=float(-1000 * freq * open_pad / result['length_mm']),
                         short_shift_MHz=float(-1000 * freq * short_pad / result['length_mm']),
                         normalized_rms=float(np.sqrt(np.mean(answer.fun ** 2))),
                         any_bound_active=bool(np.any(answer.active_mask))))
    return fits


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('source', type=Path)
    ap.add_argument('--fields', action='store_true', help='Sample stored eigenmode fields on the no-tee case')
    args = ap.parse_args()
    source = args.source.resolve()
    parent = read(source / 'run.json')
    model_manifest = Path(parent['upstream'][1]) / 'run.json'
    params = read(model_manifest)['parameters']['circuit']
    results, inputs = read_cases(source, params)
    if not results:
        raise ValueError('No solved eigenmode cases found')
    revision = parent['geometry_revision']
    inputs = [source / 'run.json', model_manifest, *inputs]
    if args.fields:
        inputs.extend((source / 'raw/no_tee/geometry.json', source / 'raw/no_tee/boundary.json'))
    folder = create_run('Design7', 'resonator', 'eigen_audit_fit', entrypoint=__file__,
                        inputs=inputs, upstream=[str(source)],
                        parameters=dict(EM_rerun=False, sample_stored_fields=args.fields, geometry_revision=revision))
    record = read(folder / 'run.json')
    record['geometry_revision'] = revision
    write(folder / 'run.json', record)
    summary = dict(cases=results, transitions=transitions(results),
                   note='The no-tee residual includes both endpoint phases, the bends and distributed CPW effects.')
    plot_comparison(source, folder, results, summary['transitions'])
    if args.fields:
        try:
            summary['terminal_field_fits'] = terminal_fields(source, results['no_tee'], folder)
        except Exception:
            summary['terminal_field_error'] = traceback.format_exc()
    finish_run(folder, status='completed', summary=summary)
    print(folder)
    print(summary)


if __name__ == '__main__':
    main()
