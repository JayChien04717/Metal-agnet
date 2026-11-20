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

CHANNELS = ('Q2', 'Q1')


def port_residues(freq, poles, data):
    center = freq.mean()
    x = (freq - center) * 1000
    p = (poles - center) * 1000
    basis = np.column_stack([*((x / 300)**n for n in range(4)), *(1 / (x - pole) for pole in p)])
    terms = [item for item in ('S11', 'S21', 'S22') if item in data]
    raw = data[terms].apply(lambda col: col.map(complex)).to_numpy()
    amplitudes = np.linalg.lstsq(basis, raw, rcond=None)[0][-len(poles):]
    residue = dict(zip(terms, amplitudes.T))
    # A reciprocal simple pole has R_ij proportional to v_i v_j. Port 1 is P8.
    p8 = abs(residue['S11'])**2 / (abs(residue['S11'])**2 + abs(residue['S21'])**2)
    rank_error = None
    if 'S22' in residue:
        product = residue['S11'] * residue['S22']
        rank_error = abs(product - residue['S21']**2) / np.maximum(abs(product), abs(residue['S21']**2))
    return p8, rank_error


def fit_spectrum(wide, fine, seeds, window):
    freq, signal = load(fine)
    mask = (freq >= window[0]) & (freq <= window[1])
    poles, model, rms = fit_candidate(freq[mask], signal[mask], seeds)
    wide_poles, _, wide_rms = fit_candidate(*load(wide), seeds)
    width = np.median(np.diff(freq)) * 2 * 1000
    data = pd.read_csv(fine, index_col=0)
    p8, rank_error = port_residues(freq[mask], poles, data.iloc[mask])
    channels = {}
    for index, channel in enumerate(CHANNELS):
        sl = slice(2 * index, 2 * index + 2)
        pair = poles[sl]
        channels[channel] = dict(**invert(*pair), centers_GHz=pair.real.tolist(),
                                 FWHM_MHz=(2000 * pair.imag).tolist(),
                                 P8_emission_share_from_residue=p8[sl].tolist(),
                                 residue_rank_error=rank_error[sl].tolist() if rank_error is not None else None,
                                 resolved=bool(2000 * min(pair.imag) > width),
                                 wide_fit=invert(*wide_poles[sl]))
    power_error = np.abs(signal)**2 + np.abs(data['S11'].map(complex))**2 - 1
    reciprocal = 'S12' in data
    if reciprocal:
        s12 = data['S12'].map(complex).to_numpy()
        s22 = data['S22'].map(complex).to_numpy()
        power_error = np.r_[power_error, np.abs(s12)**2 + np.abs(s22)**2 - 1]
    diagnostic = dict(complex_fit_rms=rms, wide_complex_fit_rms=wide_rms,
                      max_lossless_power_error=float(max(abs(power_error))),
                      max_reciprocity_error=float(max(abs(signal - s12))) if reciprocal else None,
                      S21_dB_at_6p9_GHz=float(np.interp(6.9, freq, db(signal))))
    full = join_grids(load(wide), (freq, signal))
    return channels, diagnostic, full, (freq[mask], model(freq[mask]))


def plot(folder, data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(11, 6.7), layout='constrained')
    grid = fig.add_gridspec(2, 2, height_ratios=(1.2, 1))
    full = fig.add_subplot(grid[0, :])
    details = [fig.add_subplot(grid[1, i]) for i in range(2)]
    labels = {'baseline': 'Continuous feed', 'candidate': '73 fF'}
    colors = {'baseline': ('#2c5e97', '#68a6bc'), 'candidate': ('#a3452b', '#eea52a')}
    if 'untrimmed' in data:
        labels.update(untrimmed='73 fF, original', candidate='73 fF, Q2 filter retuned')
        colors['untrimmed'] = ('#346456', '#6aac99')
    for ax in (full, *details):
        for key, (_, _, raw, fitted) in data.items():
            raw_color, fit_color = colors[key]
            ax.plot(raw[0], db(raw[1]), color=raw_color, lw=1.35, label=f'HFSS · {labels[key]}')
            ax.plot(fitted[0], db(fitted[1]), '--', color=fit_color, lw=1.15, label=f'Four-pole fit · {labels[key]}')
        ax.set(xlabel='Frequency (GHz)', ylabel='$|S_{21}|$ (dB)')
        ax.grid(alpha=.24)
    full.set(xlim=(6.55, 7.25), title='Design6_TXT loaded r10 readout · P7 to P8')
    full.legend(ncol=2, loc='lower left')
    for ax, q in zip(details, CHANNELS):
        centers = [center for ch, *_ in data.values() for center in ch[q]['centers_GHz']]
        ax.set(xlim=(min(centers) - .024, max(centers) + .024), title=q)
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
    if manifest['parameters'].get('q2_filter_add_um'):
        sources['untrimmed'] = Path(manifest['upstream'][1])
    spectra = {key: spectrum_paths(path) for key, path in sources.items()}
    if any(wide == fine for wide, fine in spectra.values()):
        parser.error('Both runs must include a fine sweep')
    folder = create_run('Design6_TXT', 'resonator', 'compare_readout_cin', entrypoint=__file__,
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
            header = ['channel', 'quantity', 'unit', 'baseline', '73_fF', 'difference_MHz']
            if 'untrimmed' in data:
                header += ['untrimmed_73_fF', 'retuned_minus_untrimmed_MHz']
            writer.writerow(header)
            for q in CHANNELS:
                def values(key):
                    ch = data[key][0][q]
                    return dict(ch, FWHM_low_MHz=ch['FWHM_MHz'][0], FWHM_high_MHz=ch['FWHM_MHz'][1])

                before, after = map(values, ('baseline', 'candidate'))
                untrimmed = values('untrimmed') if 'untrimmed' in data else None
                for field, unit in metrics.items():
                    old, new = before[field], after[field]
                    scale = 1000 if unit == 'GHz' else 1
                    delta = (new - old) * scale if None not in (old, new) else None
                    row = [q, field, unit, old, new, delta]
                    if untrimmed is not None:
                        untrimmed_value = untrimmed[field]
                        row.extend([untrimmed_value, (new - untrimmed_value) * scale
                                    if None not in (new, untrimmed_value) else None])
                    writer.writerow(row)
        finish_run(folder, status='completed', summary=summary)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
