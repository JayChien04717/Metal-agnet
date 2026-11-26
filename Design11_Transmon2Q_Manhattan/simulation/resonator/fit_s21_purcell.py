# -*- coding: utf-8 -*-
r"""Re-analysis: complex least-squares fit of a Purcell filter + readout resonator channel to HFSS driven S21.

Environment: scqenv (analysis only, no EM).
    & "C:\Users\QEL\miniforge-pypy3\envs\scqenv\python.exe" -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/resonator/fit_s21_purcell.py \
        --run Design11_Transmon2Q_Manhattan/results/resonator/<purcell_s21_driven run> [--run ...]

Model (per channel window, frequencies in GHz, linewidths in MHz):
    S21(f) = (b0 + b1 (f - fc)) * [1 - eta e^{i phi} (kappa_p/2) / D(f)],
    D(f)   = kappa_p/2 + i (f - f_p) + J^2 / (kappa_ri/2 + i (f - f_r)),
a hanger filter (total linewidth kappa_p, external-coupling ratio / feed-mismatch amplitude eta, Fano phase phi) with
the readout resonator coupled only through the filter (J). The complex linear background b0 + b1 (f - fc) absorbs
the C_in standing wave of the feed. kappa_ri (resonator loss other than through the filter) is fixed to 0.001 MHz
(HFSS has no material loss here). Several starts (resonator above / below the filter, J 2-10 MHz) are fitted and
the lowest residual kept; f_r is also seeded at the narrowest feature (largest single-point deviation). Derived: kappa_r,eff = 4 J^2 kappa_p / (kappa_p^2 + 4 Delta^2), Delta = f_r - f_p.
Limitations: the 1 MHz sweep samples the narrow resonator feature with 1-3 points, so f_r and J carry
roughly +-1 MHz / +-1 MHz uncertainty (reported 1-sigma values come from the Jacobian and are optimistic).
The parent run is not modified; results go to a new results/resonator/<run>_fit_s21_purcell folder.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general')]

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from run_record import create_run, finish_run

DESIGN = 'Design11_Transmon2Q_Manhattan'
WINDOWS = dict(Q2=(6.08, 6.21), Q1=(6.29, 6.41))
TARGET = dict(Q1=6.35, Q2=6.15)
KAPPA_RI = 0.001          # MHz


def model(p, f):
    fp, fr, J, kp, eta, phi, b0r, b0i, b1r, b1i, fc = p
    df_p, df_r = (f - fp) * 1e3, (f - fr) * 1e3                      # MHz
    D = kp / 2 + 1j * df_p + J ** 2 / (KAPPA_RI / 2 + 1j * df_r)
    return (b0r + 1j * b0i + (b1r + 1j * b1i) * (f - fc)) * (1 - eta * np.exp(1j * phi) * (kp / 2) / D)


def fit_channel(f, s, q):
    lo, hi = WINDOWS[q]
    k = (f >= lo) & (f <= hi)
    f, s = f[k], s[k]
    fc = 0.5 * (lo + hi)
    edge = np.r_[0:5, len(f) - 5:len(f)]
    b1, b0 = np.polyfit(f[edge] - fc, s[edge], 1)
    m = abs(s)
    f_dip = float(f[np.argmin(m)])
    # narrow resonator feature: largest deviation of a point from a local quadratic through its neighbours
    dev = np.zeros(len(f))
    for i in range(3, len(f) - 3):
        idx = np.r_[i - 3:i - 1, i + 2:i + 4]
        c = np.polyfit(f[idx] - f[i], s[idx], 2)
        dev[i] = abs(s[i] - np.polyval(c, 0)) / abs(s[i])
    f_narrow = float(f[np.argmax(dev)])
    fr_seeds = [f_dip + sd * dr * 1e-3 for sd in (-1, 1) for dr in (3.0, 8.0, 15.0, 30.0)]
    fr_seeds += [f_narrow + x * 1e-3 for x in (-1.0, 0.0, 1.0)]
    best = None
    for fr0 in fr_seeds:
        for J0 in (2.0, 4.0, 6.0, 9.0):
            for fp0 in (f_dip, f_dip + 5e-3, f_dip - 5e-3):
                for kp0 in (15.0, 30.0, 50.0):
                    p0 = [fp0, fr0, J0, kp0, 1.0, 0.0, b0.real, b0.imag, b1.real, b1.imag]
                    res = lambda p: np.r_[(model(list(p) + [fc], f) - s).real, (model(list(p) + [fc], f) - s).imag]
                    lb = [lo, lo, 0.1, 1.0, 0.05, -np.pi, -np.inf, -np.inf, -np.inf, -np.inf]
                    ub = [hi, hi, 40.0, 300.0, 5.0, np.pi, np.inf, np.inf, np.inf, np.inf]
                    try:
                        r = least_squares(res, p0, bounds=(lb, ub), x_scale='jac', max_nfev=4000)
                    except Exception:
                        continue
                    if best is None or r.cost < best.cost:
                        best = r
    p = best.x
    dof = max(1, 2 * len(f) - len(p))
    s2 = 2 * best.cost / dof
    try:
        cov = np.linalg.inv(best.jac.T @ best.jac) * s2
        err = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        err = np.full(len(p), np.nan)
    fp, fr, J, kp, eta, phi = p[:6]
    D = (fr - fp) * 1e3
    kr = 4 * J ** 2 * kp / (kp ** 2 + 4 * D ** 2)
    rms = float(np.sqrt(np.mean(abs(model(list(p) + [fc], f) - s) ** 2)) / np.mean(abs(s)))
    out = dict(window_GHz=[lo, hi], n_points=int(len(f)), target_GHz=TARGET[q],
               f_p_GHz=float(fp), f_r_GHz=float(fr), J_MHz=float(J), kappa_p_MHz=float(kp),
               eta=float(eta), phi_rad=float(phi), sigma=dict(f_p_MHz=float(err[0] * 1e3), f_r_MHz=float(err[1] * 1e3),
                                                             J_MHz=float(err[2]), kappa_p_MHz=float(err[3])),
               detuning_fr_minus_fp_MHz=float(D), df_p_MHz=float((fp - TARGET[q]) * 1e3),
               df_r_MHz=float((fr - TARGET[q]) * 1e3), kappa_r_eff_MHz=float(kr), rel_rms_residual=rms)
    return out, (f, s, list(p) + [fc])


def plot(fits, png, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(13, 8))
    for col, q in enumerate(('Q2', 'Q1')):
        o, (f, s, p) = fits[q]
        ff = np.linspace(f.min(), f.max(), 6000)
        sm = model(p, ff)
        a = axs[0, col]
        a.plot(f, 20 * np.log10(abs(s)), 'o', ms=3, color='#d9363e', label='HFSS (1 MHz)')
        a.plot(ff, 20 * np.log10(abs(sm)), '-', lw=1.2, color='k', label='fit')
        for x, ls, lab in ((o['f_p_GHz'], '-', 'f_p'), (o['f_r_GHz'], ':', 'f_r')):
            a.axvline(x, color='#3d8be8', ls=ls, lw=1, label=f'{lab} {x:.4f}')
        a.axvline(o['target_GHz'], color='k', ls='--', lw=0.8)
        a.axvspan(o['target_GHz'] - 0.0075, o['target_GHz'] + 0.0075, color='green', alpha=0.08)
        a.set_title(f"{q}: f_p {o['f_p_GHz']:.4f}, f_r {o['f_r_GHz']:.4f} GHz, J {o['J_MHz']:.2f}, "
                    f"κp {o['kappa_p_MHz']:.1f}, κr,eff {o['kappa_r_eff_MHz']:.2f} MHz", fontsize=9)
        a.set_ylabel('|S21| (dB)'); a.grid(alpha=0.3); a.legend(fontsize=7)
        b = axs[1, col]
        b.plot(f, np.degrees(np.angle(s)), 'o', ms=3, color='#d9363e')
        b.plot(ff, np.degrees(np.angle(sm)), '-', lw=1.0, color='k')
        b.set_xlabel('f (GHz)'); b.set_ylabel('arg S21 (deg, wrapped)'); b.grid(alpha=0.3)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(); fig.savefig(png, dpi=120); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='append', required=True, help='parent purcell_s21_driven run (repeatable)')
    a = ap.parse_args()
    for parent in a.run:
        pr = (ROOT / parent).resolve()
        csv = pr / 'raw' / 'S_wide.csv'
        d = pd.read_csv(csv)
        f, s = d.iloc[:, 0].values.astype(float), d['S21'].apply(complex).values
        tag = pr.name.split('_purcell_s21_driven_')[-1].rsplit('_', 1)[0]
        folder = create_run(DESIGN, 'resonator', f'fit_s21_{tag}', entrypoint=__file__, inputs=[csv],
                            upstream=[str(pr.relative_to(ROOT))],
                            parameters=dict(parent_run=pr.name, windows_GHz=WINDOWS, kappa_ri_MHz=KAPPA_RI))
        fits = {q: fit_channel(f, s, q) for q in ('Q1', 'Q2')}
        out = dict(parent_run=pr.name, channels={q: v[0] for q, v in fits.items()})
        (folder / 'figures').mkdir(parents=True, exist_ok=True)
        plot(fits, folder / 'figures' / 's21_fit.png', f'{pr.name}: complex Purcell-channel fit')
        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
        finish_run(folder, status='completed', summary=out)
        print('RUN_FOLDER=' + str(folder))
        print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
