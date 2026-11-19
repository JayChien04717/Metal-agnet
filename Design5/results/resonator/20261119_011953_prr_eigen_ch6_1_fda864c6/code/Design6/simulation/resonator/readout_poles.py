"""Four-pole rational fit of Design6 two-channel readout S21; per-channel conditional J inversion.

Re-analysis only (no EM rerun). Fits S21 = polynomial + sum_k r_k/(x - p_k) with complex poles
on the HFSS S_wide grid (fast-sweep interpolant) and checks it against S_fine. Each channel's pair
(Purcell-like broad pole, readout-like narrow pole) is inverted assuming one lossy Purcell mode,
real J, no intrinsic readout loss, and no cross-channel coupling."""
from argparse import ArgumentParser
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code/general'))
from run_record import create_run, finish_run


def load(path):
    d = pd.read_csv(path, index_col=0)
    return d.index.to_numpy(float), d['S21'].map(complex).to_numpy()


def fit(f, s, guesses, window, degree=3):
    m = (f >= window[0]) & (f <= window[1]); f, s = f[m], s[m]
    c = np.mean(window); x = (f - c) * 1000  # MHz
    p0 = (np.array([g[0] for g in guesses]) - c) * 1000
    w0 = np.log([g[1] for g in guesses])

    def basis(x, q):
        n = len(guesses); poles = q[:n] + 1j * np.exp(q[n:])
        return np.column_stack([*((x / 300) ** k for k in range(degree + 1)), *(1 / (x - p) for p in poles)])

    def res(q):
        b = basis(x, q); r = b @ np.linalg.lstsq(b, s, rcond=None)[0] - s
        return np.r_[r.real, r.imag]
    q = least_squares(res, np.r_[p0, w0], bounds=(np.r_[p0 - 20, w0 - 8], np.r_[p0 + 20, w0 + 4]), max_nfev=4000).x
    coeff = np.linalg.lstsq(basis(x, q), s, rcond=None)[0]
    n = len(guesses); poles = c + (q[:n] + 1j * np.exp(q[n:])) / 1000
    model = lambda ff: basis((ff - c) * 1000, q) @ coeff
    return poles, model, float(np.sqrt(np.mean(abs(model(f) - s) ** 2)))


def invert(a, b):
    """Two poles (GHz, Im = half-width) -> bare fp, fr, kappa_p, J (MHz) with one lossy mode."""
    a, b = a * 1000, b * 1000; s = a + b; p = a * b
    fr = p.imag / s.imag; fp = s.real - fr; j2 = fp * fr - p.real
    return dict(bare_fp_GHz=fp / 1000, bare_fr_GHz=fr / 1000, kappa_p_MHz=2 * s.imag,
                J_MHz=float(np.sqrt(j2)) if j2 > 0 else None, J_squared_MHz2=float(j2),
                detuning_fr_minus_fp_MHz=fr - fp)


def main():
    ap = ArgumentParser(description=__doc__); ap.add_argument('source', type=Path)
    ap.add_argument('--auto-guesses', action='store_true',
                    help='take the four pole guesses from the S_wide dips (resolved, hybridized modes)')
    args = ap.parse_args(); src = args.source.resolve()
    has_fine = (src / 'raw/S_fine.csv').is_file()
    folder = create_run('Design6', 'resonator', 'readout_poles', entrypoint=__file__,
                        inputs=[src / 'run.json', src / 'raw/S_wide.csv', *([src / 'raw/S_fine.csv'] if has_fine else [])],
                        upstream=[str(src)],
                        parameters=dict(EM_rerun=False, model='cubic background + 4 complex poles on S21',
                                        auto_guesses=args.auto_guesses, S_fine_available=has_fine,
                                        assumption='one lossy Purcell mode per channel, real J, zero intrinsic loss, no cross-channel coupling'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    fw, sw = load(src / 'raw/S_wide.csv')
    ff, sf = load(src / 'raw/S_fine.csv') if has_fine else (fw[:0], sw[:0])
    if args.auto_guesses:
        from scipy.signal import find_peaks, peak_widths
        m = (fw >= 6.6) & (fw <= 7.5); db = -20 * np.log10(abs(sw[m]))
        idx, _ = find_peaks(db, prominence=.3); hw = peak_widths(db, idx, rel_height=.5)[0] * (fw[1] - fw[0]) * 500
        assert len(idx) == 4, f'expected 4 dips, got {fw[m][idx]}'
        guesses = [(fw[m][i], max(h, .1)) for i, h in zip(idx, hw)]
    else:  # (centre GHz, half-width MHz) guesses read off the 872da39e S_wide/S_fine dips
        guesses = [(6.904, 8), (7.099, .1), (7.111, 8), (7.278, .1)]
    poles, model, rms = fit(fw, sw, guesses, (6.6, 7.5))
    fine_rms = float(np.sqrt(np.mean(abs(model(ff) - sf) ** 2))) if has_fine else None
    order = np.argsort(poles.real); poles = poles[order]
    pairs = {'Q2': (poles[0], poles[1]), 'Q1': (poles[2], poles[3])}  # Q2 is longer (target 6.80 GHz) so lower
    channels = {}
    for k, (broad, narrow) in pairs.items():
        ch = dict(wide_only=invert(broad, narrow))
        resolved = 2 * min(broad.imag, narrow.imag) * 1000 > 2 * (fw[1] - fw[0]) * 1000   # FWHM > 2 grid steps
        ch['wide_only']['note'] = ('both widths resolved on the wide grid' if resolved else
                                   'narrow pole from 1 MHz fast-sweep grid; width unresolved, J not usable')
        ch['wide_only']['widths_resolved'] = bool(resolved)
        if has_fine and ff[0] + .002 < narrow.real < ff[-1] - .002:  # refit the narrow pole on the 25 kHz fine grid
            local, _, lrms = fit(ff, sf, [(narrow.real, .2)], (narrow.real - .002, narrow.real + .002), degree=2)
            ch['fine_local'] = dict(invert(broad, local[0]), local_window_MHz=4, local_rms=lrms,
                                    kappa_r_MHz=2 * local[0].imag * 1000)
        elif not has_fine:
            ch['fine_local'] = 'no S_fine in source run'
        else:
            ch['fine_local'] = 'narrow pole outside S_fine; needs a fine sweep around it'
        channels[k] = ch
    out = dict(source=str(src), poles_GHz=[[p.real, p.imag * 1000] for p in poles], pole_units='[GHz, half-width MHz]',
               wide_rms=rms, fine_check_rms=fine_rms,
               channels=channels,
               pairing='lower two poles = Q2, upper two = Q1: Q2 lines are 3.2 % longer (filter 3.94469 vs 3.82263 mm; targets 6.80 vs 7.00 GHz). Not verified by field plots',
               limitation='Wide data is an HFSS fast-sweep interpolant. Large detuning makes J sensitive to the narrow-pole width.')
    (folder / 'poles.json').write_text(json.dumps(out, indent=1))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(10, 7))
    for a, (f, s, t) in zip(ax, [(fw, sw, 'S_wide'), (ff, sf, 'S_fine') if has_fine else (fw, sw, 'S_wide (zoom)')]):
        a.plot(f, 20 * np.log10(abs(s)), 'k', lw=1.2, label=f'HFSS {t}')
        a.plot(f, 20 * np.log10(abs(model(f))), 'r--', lw=1, label='4-pole fit')
        a.set_ylabel('|S21| (dB)'); a.grid(alpha=.3); a.legend()
    ax[0].set_xlim(6.6, 7.5); ax[1].set_xlabel('Frequency (GHz)')
    if not has_fine:
        ax[1].set_xlim(min(p.real for p in poles) - .03, max(p.real for p in poles) + .03)
    def label(k, ch):
        v = ch['fine_local'] if isinstance(ch['fine_local'], dict) else ch['wide_only']
        ok = isinstance(ch['fine_local'], dict) or ch['wide_only']['widths_resolved']
        j = f"J {v['J_MHz']:.1f} MHz" if ok and v['J_MHz'] else 'J unresolved'
        return f"{k}: fp {v['bare_fp_GHz']:.4f}, fr {v['bare_fr_GHz']:.4f} GHz, {j}"
    txt = '  '.join(label(k, ch) for k, ch in out['channels'].items())
    fig.suptitle('Design6 readout-only S21 pole fit\n' + txt, fontsize=9); fig.tight_layout()
    fig.savefig(folder / 'pole_fit.png', dpi=160)
    finish_run(folder, status='completed', summary=out)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
