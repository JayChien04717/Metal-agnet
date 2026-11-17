"""Fit and compare both loaded r10 readout channels with and without the series capacitor."""
from argparse import ArgumentParser
import csv
from pathlib import Path
import traceback

import numpy as np
import pandas as pd

from compare_loaded_s21_r9 import db, fit_candidate, read
from compare_loaded_trim_r9 import join_grids, spectrum_paths
from readout_poles import invert, load
from run_record import create_run, finish_run

LABELS = {'baseline': '0 fF', 'candidate': '73 fF'}
CHANNELS = ('Q2', 'Q1')


def fit_spectrum(wide, fine, seeds, window):
    freq, signal = load(fine)
    mask = (freq >= window[0]) & (freq <= window[1])
    poles, model, rms = fit_candidate(freq[mask], signal[mask], seeds)
    wide_poles, _, wide_rms = fit_candidate(*load(wide), seeds)
    width = np.median(np.diff(freq)) * 2 * 1000
    channels = {}
    for index, channel in enumerate(CHANNELS):
        sl = slice(2 * index, 2 * index + 2)
        pair = poles[sl]
        channels[channel] = dict(**invert(*pair), centers_GHz=pair.real.tolist(),
                                 FWHM_MHz=(2000 * pair.imag).tolist(),
                                 resolved=bool(2000 * min(pair.imag) > width),
                                 wide_fit=invert(*wide_poles[sl]))
    data = pd.read_csv(fine, index_col=0)
    power_error = np.abs(signal)**2 + np.abs(data['S11'].map(complex))**2 - 1
    reciprocal = 'S12' in data
    if reciprocal:
        s12 = data['S12'].map(complex).to_numpy()
        s22 = data['S22'].map(complex).to_numpy()
        power_error = np.r_[power_error, np.abs(s12)**2 + np.abs(s22)**2 - 1]
    diagnostic = dict(complex_fit_rms=rms, wide_complex_fit_rms=wide_rms,
                      max_lossless_power_error=float(max(abs(power_error))),
                      max_reciprocity_error=float(max(abs(signal - s12))) if reciprocal else None)
    full = join_grids(load(wide), (freq, signal))
    return channels, diagnostic, full, (freq[mask], model(freq[mask]))


def plot(folder, data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10.6, 4.8), layout='constrained')
    for (key, (channels, _, raw, fitted)), colors in zip(data.items(), (('#2c5e97', '#68a6bc'), ('#a3452b', '#eea52a'))):
        ax.plot(raw[0], db(raw[1]), color=colors[0], lw=1.35, label=f'HFSS · {LABELS[key]}')
        ax.plot(fitted[0], db(fitted[1]), '--', color=colors[1], lw=1.15, label=f'Four-pole fit · {LABELS[key]}')
    ax.set(xlim=(6.55, 7.25), xlabel='Frequency (GHz)', ylabel='$|S_{21}|$ (dB)',
           title='Design6 loaded r10 readout · P7 to P8')
    ax.legend(ncol=2, loc='lower left')
    ax.grid(alpha=.24)
    fig.savefig(folder / 's21_cin_comparison.png', dpi=180)
    plt.close(fig)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('candidate', type=Path)
    args = parser.parse_args()
    candidate = args.candidate.resolve()
    manifest = read(candidate / 'run.json')
    if manifest['status'] != 'completed' or manifest['parameters']['dry_run']:
        parser.error('The capacitor HFSS solve must be complete')
    baseline = Path(manifest['upstream'][0])
    sources = dict(baseline=baseline, candidate=candidate)
    spectra = {key: spectrum_paths(path) for key, path in sources.items()}
    if any(wide == fine for wide, fine in spectra.values()):
        parser.error('Both runs must include a fine sweep')
    folder = create_run('Design6', 'resonator', 'compare_readout_cin', entrypoint=__file__,
                        upstream=[str(p) for p in sources.values()],
                        inputs=[p / 'run.json' for p in sources.values()] + [p for pair in spectra.values() for p in pair],
                        parameters=dict(candidate_only=True, EM_rerun=False, parent_geometry=manifest['parameters']['parent_geometry'],
                                        fit='complex S21, four poles and cubic baseline, common fine window',
                                        inversion='real J, one damped mode per channel, zero separate resonator loss',
                                        direction='P7 to P8 via reciprocity of the two-port fixture'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        grids = [load(pair[1])[0] for pair in spectra.values()]
        window = [max(f[0] for f in grids), min(f[-1] for f in grids)]
        seeds = dict(Q2={'poles_GHz': [(6.802, .012), (6.833, .004)]},
                     Q1={'poles_GHz': [(6.978, .008), (7.004, .010)]})
        data = {key: fit_spectrum(*pair, seeds, window) for key, pair in spectra.items()}
        summary = dict(fit_window_GHz=window, variants={})
        for key, (channels, diagnostic, raw, fitted) in data.items():
            summary['variants'][key] = dict(channels=channels, diagnostics=diagnostic,
                                             solver=read(sources[key] / 'summary.json'))
            pd.DataFrame({'f_GHz': raw[0], 'S21_real': raw[1].real, 'S21_imag': raw[1].imag,
                          'S21_dB': db(raw[1])}).to_csv(folder / f'raw/{key}_s21.csv', index=False)
            pd.DataFrame({'f_GHz': fitted[0], 'S21_fit_real': fitted[1].real, 'S21_fit_imag': fitted[1].imag,
                          'S21_fit_dB': db(fitted[1])}).to_csv(folder / f'raw/{key}_fit.csv', index=False)
        plot(folder, data)
        metrics = dict(bare_fp_GHz='GHz', bare_fr_GHz='GHz', J_MHz='MHz', kappa_p_MHz='MHz',
                       FWHM_low_MHz='MHz', FWHM_high_MHz='MHz')
        with (folder / 'fit_comparison.csv').open('w', encoding='utf-8', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['channel', 'quantity', 'unit', 'baseline', '73_fF', 'difference_MHz'])
            for q in CHANNELS:
                def values(key):
                    ch = data[key][0][q]
                    return dict(ch, FWHM_low_MHz=ch['FWHM_MHz'][0], FWHM_high_MHz=ch['FWHM_MHz'][1])

                before, after = map(values, LABELS)
                for field, unit in metrics.items():
                    old, new = before[field], after[field]
                    delta = (new - old) * (1000 if unit == 'GHz' else 1) if None not in (old, new) else None
                    writer.writerow([q, field, unit, old, new, delta])
        finish_run(folder, status='completed', summary=summary)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
