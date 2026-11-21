"""Fit the loaded r9 candidate S21 and compare with variant G and the r9 baseline."""
from argparse import ArgumentParser
import json
from pathlib import Path
import traceback

import numpy as np
from scipy.signal import find_peaks

from readout_poles import ROOT, fit, invert, load
from run_record import create_run, finish_run

RESULTS = ROOT / 'Design6_TXT/results/resonator'
BASE = RESULTS / '20261117_071647_readout_s21_clawdl_cf91f657'
BASE_POLES = RESULTS / '20261117_113914_readout_poles_af4e9eaf'
BASE_MODEL = RESULTS / '20261117_133256_readout_model_q3d_07c856b2'
FREQ = dict(Q1=7.0, Q2=6.8)
METRICS = ('bare_fp_GHz', 'bare_fr_GHz', 'J_MHz', 'kappa_p_MHz')


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def db(s):
    return 20 * np.log10(np.maximum(abs(s), 1e-10))


def fit_candidate(freq, s21, seed_channels):
    target_window = (max(freq[0], 6.61), min(freq[-1], 7.18))
    idx, props = find_peaks(-db(s21), prominence=.4)
    centers = freq[idx].tolist() if len(idx) == 4 else []
    if not centers:
        for q, target in sorted(FREQ.items(), key=lambda item: item[1]):
            peaks = sorted((i for i in range(len(idx)) if abs(freq[idx[i]] - target) < .095),
                           key=lambda i: props['prominences'][i], reverse=True)[:2]
            seen = [freq[idx[i]] for i in peaks]
            predicted = [p[0] for p in seed_channels[q]['poles_GHz']]
            if len(seen) == 1:
                seen += [max(predicted, key=lambda p: abs(p - seen[0]))]
            centers.extend(sorted(seen if seen else predicted))
    if len(set(centers)) != 4:
        raise RuntimeError(f'Ambiguous pole seed locations {centers}')
    best = None
    for widths in ((7, 7, 8, 8), (3, 12, 3, 12), (12, 3, 12, 3), (.2, 9, .2, 9), (9, .2, 9, .2)):
        candidate = fit(freq, s21, list(zip(centers, widths)), target_window)
        poles = np.sort_complex(candidate[0])
        if min(np.diff(poles.real)) > .001 and (best is None or candidate[2] < best[2]):
            best = (poles, *candidate[1:])
    if best is None:
        raise RuntimeError(f'No distinct four-pole fit from {centers}')
    return best


def compare(hfss, circuit):
    return {k: dict(hfss=hfss[k], q3d_model=circuit[k],
                    hfss_minus_model=(None if hfss[k] is None else
                                      (hfss[k] - circuit[k]) * (1000 if k.endswith('GHz') else 1)))
            for k in METRICS}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('loaded', type=Path)
    parser.add_argument('length_run', type=Path)
    args = parser.parse_args()
    loaded, tuned = args.loaded.resolve(), args.length_run.resolve()
    wide, fine = [loaded / 'raw' / name for name in ('S_wide.csv', 'S_fine.csv')]
    sim = read(loaded / 'summary.json')
    baseline = read(BASE_MODEL / 'compare.json')['variants']['G']
    baseline_hfss = read(BASE_POLES / 'summary.json')['channels']
    candidate = read(tuned / 'summary.json')['tuned']
    spectra = [BASE_MODEL / 'raw/model_s21.npz', tuned / 'raw/candidate_s21.npz']
    folder = create_run('Design6_TXT', 'resonator', 'compare_loaded_s21_r9', entrypoint=__file__,
                        inputs=[wide, *([fine] if fine.exists() else []), loaded / 'summary.json',
                                tuned / 'summary.json', BASE / 'raw/S_wide.csv', BASE_POLES / 'summary.json',
                                BASE_MODEL / 'compare.json', *spectra],
                        upstream=[str(p) for p in (loaded, tuned, BASE, BASE_POLES, BASE_MODEL)],
                        parameters=dict(EM_rerun=False, fit='4 complex poles, complex S21, cubic background',
                                        inversion='one lossy mode per channel; zero separate resonator loss',
                                        new_hfss_geometry='candidate with both claws and qubits',
                                        new_model_geometry='candidate with estimated 0.100 mm claw correction'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        fw, sw = load(wide)
        ff, sf = load(fine) if fine.exists() else (fw, sw)
        poles, fitted, rms = fit_candidate(ff, sf, candidate)
        coarse, _, coarse_rms = fit_candidate(fw, sw, candidate)
        step_mhz = np.median(np.diff(ff)) * 1000
        import pandas as pd
        s11 = pd.read_csv(fine if fine.exists() else wide)['S11'].map(complex).to_numpy()
        power_error = float(np.max(np.abs(np.abs(sf)**2 + np.abs(s11)**2 - 1)))
        out = dict(loaded_solver=sim, fit_grid='Fine' if fine.exists() else 'Wide',
                   fit_grid_step_MHz=float(step_mhz), fit_rms=rms, wide_rms=coarse_rms,
                   max_lossless_two_port_power_error=power_error, channels={})
        for i, q in enumerate(('Q2', 'Q1')):
            pair = poles[2*i:2*i+2]
            primary, check = invert(*pair), invert(*coarse[2*i:2*i+2])
            widths = [2000 * p.imag for p in pair]
            out['channels'][q] = dict(
                candidate=compare(primary, candidate[q]),
                baseline=compare(baseline_hfss[q]['wide_only'], {k: baseline[q][k]['model'] for k in METRICS}),
                loaded_detuning_MHz=primary['detuning_fr_minus_fp_MHz'],
                hybrid_centers_GHz=[float(p.real) for p in pair],
                hybrid_FWHM_MHz=[float(w) for w in widths],
                resolved=bool(min(widths) > 2 * step_mhz),
                fine_minus_wide={k: ((primary[k] - check[k]) * (1000 if k.endswith('GHz') else 1)
                                      if primary[k] is not None and check[k] is not None else None) for k in METRICS})
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fb, sb = load(BASE / 'raw/S_wide.csv')
        with np.load(spectra[0]) as b, np.load(spectra[1]) as c:
            fig, axes = plt.subplots(2, 2, figsize=(11, 6), sharey=True)
            for column, q in enumerate(('Q2', 'Q1')):
                ax, bx = axes[:, column]
                ax.plot(fb, db(sb), color='#223755', label='HFSS')
                ax.plot(b['f_GHz'], db(b['G']), color='#c36a46', label='Variant G')
                f1 = baseline_hfss[q]['wide_only']['bare_fp_GHz']
                ax.set(xlim=(f1-.075, f1+.11), title=f'{q} · r9 baseline, readout only')
                bx.plot(ff, db(sf), color='#223755', label='HFSS loaded')
                bx.plot(ff, db(fitted(ff)), '--', color='#44857b', label='Four-pole fit')
                bx.plot(c['f_GHz'], db(c['S21']), color='#c36a46', label='Variant G + claw correction')
                fp = out['channels'][q]['candidate']['bare_fp_GHz']['hfss']
                fr = out['channels'][q]['candidate']['bare_fr_GHz']['hfss']
                for x in (fp, fr):
                    bx.axvline(x, alpha=.38, color='#223755', linewidth=.9)
                bx.set(xlim=(min(fp, fr, FREQ[q])-.060, max(fp, fr, FREQ[q])+.060),
                       xlabel='Frequency (GHz)', title=f'{q} · length candidate, loaded HFSS')
            for ax in axes.flat:
                ax.grid(alpha=.2)
                ax.legend(fontsize=8, loc='lower right')
            for ax in axes[:, 0]:
                ax.set_ylabel('|S21| (dB)')
            fig.tight_layout()
            fig.savefig(folder / 's21_comparison.png', dpi=170)
            plt.close(fig)
        np.savez(folder / 'raw/s21_fit.npz', f_GHz=ff, hfss=sf, fitted=fitted(ff))
        finish_run(folder, status='completed', summary=out)
        print(json.dumps(out['channels'], indent=2))
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
