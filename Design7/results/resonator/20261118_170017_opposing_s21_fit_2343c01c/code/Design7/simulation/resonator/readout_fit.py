"""Summarize the complex HFSS S21 of the Design7 double-U fixture."""

from argparse import ArgumentParser
import json
from pathlib import Path
import sys
import traceback

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths

HOME = Path(__file__).resolve().parents[2]
ROOT = HOME.parent
sys.path[:0] = [str(ROOT / 'Design6/simulation/resonator'), str(ROOT / 'code/general')]
import readout_poles as poles
from run_record import create_run, finish_run


def db(s):
    return 20 * np.log10(np.maximum(np.abs(s), 1e-12))


def peak_guesses(f, s):
    depth = -db(s)
    peaks, properties = find_peaks(depth, prominence=.1)
    if len(peaks) == 0:
        return []
    # Peak widths in power avoid instability when the exact null is very deep.
    widths = peak_widths(1 - abs(s) ** 2, peaks)[0] * np.median(np.diff(f)) * 500
    best = sorted(np.argsort(properties['prominences'])[-2:], key=lambda i: peaks[i])
    return [(float(f[peaks[i]]), float(max(widths[i], .1))) for i in best]


def fit_spectrum(f, s, guesses):
    grid_mhz = 1000 * np.median(np.diff(f))
    # Retain the full fine grid near each resonance and thin the broadband background.
    use = np.arange(len(f)) % max(len(f) // 650, 1) == 0
    for center, half_width in guesses:
        use |= np.abs(f - center) < min(.025, (5 * half_width + 1) / 1000)
    limits = (max(f[0], guesses[0][0] - .09), min(f[-1], guesses[-1][0] + .09))
    fitted, model, rms = poles.fit(f[use], s[use], guesses, limits)
    fitted = sorted(fitted, key=lambda p: p.real)
    linewidths = [float(2000 * p.imag) for p in fitted]
    summary = dict(pole_frequencies_GHz=[float(p.real) for p in fitted], fwhm_MHz=linewidths,
                   window_GHz=limits, rms_complex_S21=rms, fast_grid_spacing_MHz=grid_mhz,
                   widths_resolved=bool(min(linewidths) > 2 * grid_mhz))
    if len(fitted) == 2:
        inversion = poles.invert(*fitted)
        near = abs(inversion['detuning_fr_minus_fp_MHz']) <= 2 * (inversion['J_MHz'] or 0)
        summary.update(inversion=inversion, near_resonance=bool(near))
    return summary, model


def analyze(case, out, manual_guesses):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    wide_f, wide_s = poles.load(case / 'S_wide.csv')
    fine_f, fine_s = poles.load(case / 'S_fine.csv')
    guess = manual_guesses or peak_guesses(fine_f, fine_s)
    conv = pd.read_csv(case / 'convergence.csv')
    deltas = conv['Max Mag. Delta S'].dropna()
    last_delta = float(deltas.iloc[-1]) if len(deltas) else None
    fit, model = fit_spectrum(fine_f, fine_s, guess) if guess else (None, None)
    result = dict(adaptive_passes=len(conv), last_delta_s=last_delta, delta_s_target=.01,
                  delta_s_target_met=last_delta is not None and last_delta <= .01,
                  minima_GHz=[center for center, _ in guess], fit=fit,
                  fit_assumptions='One lossy filter mode, one lossless readout mode, real J. Fast sweep is an interpolant.')
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), layout='constrained')
    axes[0].plot(wide_f, db(wide_s), color='#215a8e', label='HFSS wide')
    axes[1].plot(fine_f, db(fine_s), color='#215a8e', label='HFSS fine')
    if fit:
        m = (fine_f >= fit['window_GHz'][0]) & (fine_f <= fit['window_GHz'][1])
        axes[1].plot(fine_f[m], db(model(fine_f[m])), '--', color='#bd5328', label='Complex pole fit')
        axes[1].set_xlim(fit['window_GHz'])
        for center in fit['pole_frequencies_GHz']:
            axes[1].axvline(center, color='0.55', linewidth=.7, linestyle=':')
    for axis in axes:
        axis.set_ylabel('|S21| (dB)')
        axis.set_xlabel('Frequency (GHz)')
        axis.grid(alpha=.2)
        axis.legend(loc='lower right')
    fig.suptitle('Design7 double-U fixture — driven HFSS S21')
    fig.savefig(out / 's21.png', dpi=165)
    plt.close(fig)
    return result


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('hfss_run', type=Path)
    parser.add_argument('--guess', nargs=2, type=float, metavar=('CENTER_GHZ', 'HALF_WIDTH_MHZ'),
                        action='append')
    args = parser.parse_args()
    source = args.hfss_run.resolve()
    case = source / 'raw/opposing_full'
    upstream = json.loads((source / 'run.json').read_text(encoding='utf-8'))
    revision = upstream['geometry_revision']
    folder = create_run(
        'Design7', 'resonator', 'opposing_s21_fit', entrypoint=__file__,
        inputs=[source / 'run.json', *[case / f for f in ('S_wide.csv', 'S_fine.csv', 'convergence.csv')]],
        upstream=[str(source)], parameters=dict(geometry_revision=revision, manual_guesses=args.guess))
    run_path = folder / 'run.json'
    record = json.loads(run_path.read_text(encoding='utf-8'))
    record['geometry_revision'] = revision
    run_path.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(f'RUN_FOLDER={folder}', flush=True)
    try:
        summary = analyze(case, folder, sorted(args.guess) if args.guess else None)
    except Exception:
        finish_run(folder, status='failed', summary=dict(error=traceback.format_exc()))
        raise
    finish_run(folder, status='completed', summary=summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
