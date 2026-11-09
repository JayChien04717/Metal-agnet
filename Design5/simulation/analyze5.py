"""Analyse one single-channel Design5 HFSS run and append it to Design5/run_record.md. Only (Jay 2026-11-09):
  1. capacitances: Q3D (series C_in, C_p, C_j; source 'full' = whole-channel Q3D, or 'part') and fitted to HFSS,
  2. HFSS S21 and S11 (results/<time>_<tag>/<tag>_S.csv),
  3. one plot: raw HFSS S21 / S11, theory with the Q3D C (no fit), theory fitted by C
     (filter / resonator length corrections, series C scales, v scale) -> how far the C model is from HFSS.

Usage: python analyze5.py <tag> [note]       (needs results/<tag>_S.csv with complex S11, S21; hfss_chan5.py calls it)
       env RECORD_MD = other md file (tests), NO_ARCHIVE=1 keeps the files in results/
"""
import os, sys, json, glob, copy, shutil, time, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.optimize import least_squares, differential_evolution
from scipy.signal import find_peaks
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results'); sys.path.insert(0, HERE)
import theory5 as T
RECORD = os.environ.get('RECORD_MD') or os.path.join(os.path.dirname(HERE), 'run_record.md')

FIT_NAMES = ['dLf', 'dLr', 's_in', 's_p', 's_j', 'kv']
FIT_BOUNDS = [(-0.4, 0.4), (-0.4, 0.4), (0.7, 1.3), (0.7, 1.3), (0.3, 3.0), (0.95, 1.05)]
PRIOR = dict(s_in=0.15, s_p=0.15, s_j=0.5, kv=0.02)


def apply(P0, x):
    p = dict(zip(FIT_NAMES, x)); P = copy.deepcopy(P0)
    P['L_filt'] += p['dLf']; P['L_res'] += p['dLr']
    P['cin'] = (P0['cin'][0] * p['s_in'],) + tuple(P0['cin'][1:])      # only the series C is scaled
    P['cp'] = (P0['cp'][0] * p['s_p'],) + tuple(P0['cp'][1:])
    P['cj'] = (P0['cj'][0] * p['s_j'],) + tuple(P0['cj'][1:])
    P['v'] = P0['v'] * p['kv']
    return P


def theory_fit(P0, f, s11, s21):
    def res(x):
        t11, t21 = T.s21_theory(f, apply(P0, x))
        r = np.r_[(t21 - s21).real, (t21 - s21).imag, 0.3 * (t11 - s11).real, 0.3 * (t11 - s11).imag] / np.sqrt(len(f))
        p = dict(zip(FIT_NAMES, x))
        return np.r_[30 * r, 0.02 * np.array([(p[k] - 1) / s for k, s in PRIOR.items()])]
    de = differential_evolution(lambda x: float(np.sum(res(x) ** 2)), FIT_BOUNDS, seed=1, maxiter=120, popsize=18, tol=1e-8, polish=False)
    lo, hi = zip(*FIT_BOUNDS)
    sol = least_squares(res, de.x, bounds=(lo, hi), x_scale=[0.02, 0.02, .05, .05, .1, .005])
    return dict(zip(FIT_NAMES, map(float, sol.x))), apply(P0, sol.x)


def dips(f, s, prom=3):
    a = 20 * np.log10(np.abs(s)); i, _ = find_peaks(-a, prominence=prom, distance=max(1, int(0.002 / (f[1] - f[0]))))
    return [round(float(f[k]), 4) for k in i]


def analyze(tag, note=''):
    S = pd.read_csv(os.path.join(RES, f'{tag}_S.csv'), index_col=0).map(complex)
    f = S.index.values.astype(float); s11, s21 = S.S11.values, S.S21.values
    meta_p = os.path.join(RES, f'{tag}_meta.json'); meta = json.load(open(meta_p)) if os.path.exists(meta_p) else {}
    P0 = T.params_from_layout()
    xfit, Pf = theory_fit(P0, f, s11, s21)
    t11, t21 = T.s21_theory(f, P0); u11, u21 = T.s21_theory(f, Pf)
    rms = lambda t: float(np.sqrt(np.mean(np.abs(t - s21) ** 2)))
    out = dict(tag=tag, date=time.strftime('%Y-%m-%d %H:%M'), note=note or meta.get('note', ''), hfss=meta.get('hfss', {}),
               cap_source=P0['cap_source'], caps_q3d=dict(cin=P0['cin'], cp=P0['cp'], cj=P0['cj']),
               caps_fit=dict(cin=Pf['cin'], cp=Pf['cp'], cj=Pf['cj']), fit_params=xfit,
               dips=dict(hfss=dips(f, s21), q3d=dips(f, t21), fit=dips(f, u21)), rms_S21=dict(q3d=rms(t21), fit=rms(u21)))
    # ---- one figure: S21 and S11, raw HFSS / Q3D C / fitted ----
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    db = lambda s: 20 * np.log10(np.abs(s))
    fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    sty = [('HFSS (raw)', 'k', '-', 1.6), (f'theory, Q3D C ({P0["cap_source"]}, no fit)', 'C0', ':', 1.3), ('theory fitted by C', 'C3', '--', 1.1)]
    for a, name, curves in ((ax[0], 'S21', (s21, t21, u21)), (ax[1], 'S11', (s11, t11, u11))):
        for (lab, c, ls, lw), s in zip(sty, curves):
            a.plot(f, db(s).clip(-90), color=c, ls=ls, lw=lw, label=lab)
        a.set_ylabel(f'{name} (dB)'); a.legend(fontsize=8); a.grid(alpha=.3)
    ax[0].set_title(f'{tag}  {out["note"]}', fontsize=9); ax[1].set_xlabel('f (GHz)')
    fig.tight_layout()
    # ---- archive ----
    import runstore                                       # results/<date>/<HHMM>_<tag>/, nothing overwritten
    run_dir = RES if os.environ.get('NO_ARCHIVE') else runstore.new_run_dir(tag, RES)
    if run_dir != RES:
        for p in glob.glob(os.path.join(RES, f'{tag}_*')) + glob.glob(os.path.join(RES, f'{tag}.log')):
            try:
                shutil.move(p, os.path.join(run_dir, os.path.basename(p)))
            except (PermissionError, OSError):
                shutil.copy(p, os.path.join(run_dir, os.path.basename(p)))
        for p in glob.glob(os.path.join(RES, 'q3d5_*_caps.json')):
            shutil.copy(p, os.path.join(run_dir, f'{tag}_' + os.path.basename(p)))
    st = runstore.stamp()
    png = os.path.join(run_dir, f'{tag}_compare_{st}.png'); fig.savefig(png, dpi=120); plt.close(fig)
    json.dump(out, open(os.path.join(run_dir, f'{tag}_compare_{st}.json'), 'w'), indent=1, default=float)
    try:
        rel = os.path.relpath(png, os.path.dirname(RECORD))
    except ValueError:                                   # record on another drive (tests)
        rel = png
    record_md(out, rel.replace('\\', '/'))
    return out, run_dir


def record_md(o, png_rel):
    c, cf, p = o['caps_q3d'], o['caps_fit'], o['fit_params']; h = o.get('hfss', {})
    pct = lambda a, b: f'{100 * (b / a - 1):+.1f} %'
    L = [f'\n## {o["tag"]}  ({o["date"]})\n', f'{o["note"]}\n' if o['note'] else '',
         f'HFSS: {", ".join(f"{k} {v}" for k, v in h.items())}\n' if h else '',
         f'![{o["tag"]}]({png_rel.replace(" ", "%20")})\n',
         f'**1. Capacitances** (series, fF; Q3D source `{o["cap_source"]}`)\n',
         '| | C_in | C_p | C_j |', '|---|---|---|---|',
         f'| Q3D | {c["cin"][0]:.2f} | {c["cp"][0]:.2f} | {c["cj"][0]:.2f} |',
         f'| fitted to HFSS S21 | {cf["cin"][0]:.2f} | {cf["cp"][0]:.2f} | {cf["cj"][0]:.2f} |',
         f'| difference | {pct(c["cin"][0], cf["cin"][0])} | {pct(c["cp"][0], cf["cp"][0])} | {pct(c["cj"][0], cf["cj"][0])} |',
         '', f'Fitted ΔL filter {p["dLf"]:+.4f} mm, ΔL resonator {p["dLr"]:+.4f} mm, v scale {p["kv"]:.4f}; '
             f'rms |S21 − HFSS|: Q3D model {o["rms_S21"]["q3d"]:.4f}, fitted {o["rms_S21"]["fit"]:.4f}\n',
         f'**2. S21 dips** (GHz): HFSS {o["dips"]["hfss"]}; theory, Q3D C {o["dips"]["q3d"]}; theory fitted {o["dips"]["fit"]}\n']
    if not os.path.exists(RECORD):
        open(RECORD, 'w', encoding='utf-8').write('# Design5 run record\n')
    open(RECORD, 'a', encoding='utf-8').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    o, d = analyze(sys.argv[1], ' '.join(sys.argv[2:]))
    print(json.dumps(dict(caps_q3d=o['caps_q3d'], caps_fit=o['caps_fit'], dips=o['dips'], rms=o['rms_S21']), indent=1, default=float))
    print('run folder', d, '\nrecord', RECORD)
