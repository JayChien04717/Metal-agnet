"""Compare the loaded Design6_TXT r9 HFSS pole fits before and after tail trims."""
import argparse
import csv
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

from compare_loaded_s21_r9 import FREQ, METRICS, db, fit_candidate, read
from readout_poles import invert, load
from run_record import create_run, finish_run


def spectrum_paths(run):
    wide = run / 'raw/S_wide.csv'
    fine = run / 'raw/S_fine.csv'
    return wide, fine if fine.is_file() else wide


def reference_seeds(previous):
    seeds = {}
    for q, target in FREQ.items():
        values = previous[q]['candidate']
        g = values['J_MHz']['hfss'] / 1000
        gamma = values['kappa_p_MHz']['hfss'] / 2000
        pair = np.linalg.eigvals([[target + 1j * gamma, g], [g, target]])
        seeds[q] = {'poles_GHz': [(pole.real, pole.imag) for pole in pair]}
    return seeds


def join_grids(wide, fine):
    fw, sw = wide
    ff, sf = fine
    lo, hi = fw < ff[0], fw > ff[-1]
    return np.r_[fw[lo], ff, fw[hi]], np.r_[sw[lo], sf, sw[hi]]


def compare_channels(poles, wide_poles, previous, grid_step_mhz):
    channels = {}
    for i, q in enumerate(('Q2', 'Q1')):
        pair = poles[2*i:2*i+2]
        current = invert(*pair)
        coarse = invert(*wide_poles[2*i:2*i+2])
        prior = previous[q]
        quantities = {}
        for key in (*METRICS, 'detuning_fr_minus_fp_MHz'):
            old = (prior['loaded_detuning_MHz'] if key == 'detuning_fr_minus_fp_MHz'
                   else prior['candidate'][key]['hfss'])
            new = current[key]
            factor = 1000 if key.endswith('GHz') else 1
            quantities[key] = dict(previous=old, trimmed=new,
                                   trimmed_minus_previous_MHz=(None if new is None or old is None
                                                               else factor * (new - old)),
                                   fine_minus_wide_MHz=(None if new is None or coarse[key] is None
                                                        else factor * (new - coarse[key])))
        widths = (2000 * pair.imag).tolist()
        resonator_fraction = None
        if current['J_MHz'] is not None and current['J_MHz'] > 0:
            # The second row of the coupled-mode eigenproblem gives v_r/v_p.
            ratio = current['J_MHz'] / 1000 / (pair - current['bare_fr_GHz'])
            resonator_fraction = (abs(ratio)**2 / (1 + abs(ratio)**2)).tolist()
        channels[q] = dict(fit=quantities,
                           target_GHz=FREQ[q],
                           target_error_MHz={key: 1000 * (current[key] - FREQ[q]) for key in METRICS[:2]},
                           hybrid_centers_GHz=dict(previous=prior['hybrid_centers_GHz'], trimmed=pair.real.tolist()),
                           hybrid_FWHM_MHz=dict(previous=prior['hybrid_FWHM_MHz'], trimmed=widths),
                           estimated_resonator_fraction_low_high=resonator_fraction,
                           resolved_on_fit_grid=min(widths) > 2 * grid_step_mhz)
    return channels


def plot(path, previous, trimmed, fine, fitted, channels):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fp, sp = previous
    ft, st = trimmed
    ff, sf = fine
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.8), sharey='row', layout='constrained')
    for col, q in enumerate(('Q2', 'Q1')):
        full, detail = axes[:, col]
        target = FREQ[q]
        report = channels[q]
        old, new = (report['hybrid_centers_GHz'][item] for item in ('previous', 'trimmed'))
        full.plot(fp, db(sp), color='#bc6740', label='Previous loaded')
        full.plot(ft, db(st), color='#203755', label='Trimmed loaded')
        full.set(xlim=(min(old + new) - .025, max(old + new) + .025), title=f'{q} · Loaded HFSS')
        detail.plot(ff, db(sf), color='#203755', label='Trimmed loaded')
        detail.plot(ff, db(fitted(ff)), '--', color='#198272', label='Four-pole fit')
        for key, color, label in (('bare_fp_GHz', '#87518a', '$f_p$'), ('bare_fr_GHz', '#cf8c37', '$f_r$')):
            detail.axvline(report['fit'][key]['trimmed'], color=color, linestyle=':', label=label)
        detail.set(xlim=(min(new + [target]) - .025, max(new + [target]) + .025),
                   xlabel='Frequency (GHz)', title=f'{q} · Trimmed candidate pole fit')
    for ax in axes.flat:
        ax.grid(alpha=.25)
        ax.legend(fontsize=8, loc='lower right')
    for ax in axes[:, 0]:
        ax.set_ylabel('|S21| (dB)')
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trimmed_run', type=Path)
    parser.add_argument('previous_fit', type=Path)
    args = parser.parse_args()
    run, fit_run = args.trimmed_run.resolve(), args.previous_fit.resolve()
    manifest = read(run / 'run.json')
    if manifest['status'] != 'completed':
        parser.error('The trimmed HFSS run must be completed')
    previous_run = Path(manifest['upstream'][0]).resolve()
    if str(previous_run) not in read(fit_run / 'run.json')['upstream']:
        parser.error('The comparison fit does not correspond to the untrimmed loaded parent')
    wide, fine = spectrum_paths(run)
    old_wide, old_fine = spectrum_paths(previous_run)
    settings = manifest['parameters']
    folder = create_run('Design6_TXT', 'resonator', 'compare_loaded_trim_r9', entrypoint=__file__,
                        inputs=list(dict.fromkeys([wide, fine, old_wide, old_fine, run / 'run.json',
                                                  run / 'summary.json', fit_run / 'summary.json'])),
                        upstream=[str(run), str(previous_run), str(fit_run)],
                        parameters=dict(parent_geometry=settings['parent_geometry'], candidate_only=True,
                                        EM_rerun=False, target_GHz=FREQ, length_offsets_um=settings['length_offsets_um'],
                                        fit='Four complex poles on complex S21 with a cubic background',
                                        inversion='Two coupled modes per channel, one damped bare mode'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        previous = read(fit_run / 'summary.json')
        seeds = reference_seeds(previous['channels'])
        wide_data, fine_data = load(wide), load(fine)
        poles, fitted, rms = fit_candidate(*fine_data, seeds)
        wide_poles, _, wide_rms = fit_candidate(*wide_data, seeds)
        step = float(np.median(np.diff(fine_data[0])) * 1000)
        channels = compare_channels(poles, wide_poles, previous['channels'], step)
        ff, sf = fine_data
        s11 = pd.read_csv(fine)['S11'].map(complex).to_numpy()
        summary = dict(channels=channels, lengths_mm=settings['lengths_mm'],
                       length_offsets_um=settings['length_offsets_um'], trimmed_solver=read(run / 'summary.json'),
                       previous_solver=previous['loaded_solver'],
                       fit_grid='Fine' if fine != wide else 'Wide', fit_grid_step_MHz=step,
                       complex_rms=rms, wide_complex_rms=wide_rms,
                       max_lossless_two_port_power_error=float(max(abs(abs(sf)**2 + abs(s11)**2 - 1))),
                       width_definition='FWHM = kappa/(2 pi) for a pole with positive imaginary half-width; '
                                        'hybrid poles shown in increasing real frequency order')
        plot(folder / 's21_comparison.png', join_grids(load(old_wide), load(old_fine)),
             join_grids(wide_data, fine_data), fine_data, fitted, channels)
        np.savez(folder / 'raw/fitted_s21.npz', f_GHz=ff, S21=sf, S21_fit=fitted(ff), poles_GHz=poles)
        with (folder / 'fit_comparison.csv').open('w', encoding='utf-8', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['channel', 'parameter', 'unit', 'previous', 'trimmed', 'trimmed_minus_previous_MHz'])
            for q in FREQ:
                for key, values in channels[q]['fit'].items():
                    writer.writerow([q, key, 'GHz' if key.endswith('GHz') else 'MHz',
                                     values['previous'], values['trimmed'], values['trimmed_minus_previous_MHz']])
        finish_run(folder, status='completed', summary=summary)
        for q in FREQ:
            print(q, channels[q])
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
