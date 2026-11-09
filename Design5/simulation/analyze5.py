"""Analyse one Design5 HFSS run and append it to Design5/run_record.md.

For every run it records
  1. the capacitances used (series C_in, C_p, C_j and the pi shunts; source 'full' = whole-channel Q3D, or 'part'),
  2. HFSS S21 and S11 (dips, peak, plot),
  3. theory with the Q3D C (no fitting)                       -> omega_p, omega_r, J, kappa_p from the circuit,
  4. theory FIT by C: lengths + capacitor scales fitted to the complex S21 / S11 -> fitted C, then omega_p.. from it,
  5. DIRECT fit of S21 with eq. 2.137 (+ background)          -> omega_p, omega_r, J, kappa_p,
and a comparison table of 3-5.

Usage: python analyze5.py <tag> [note]       (needs results/<tag>_S.csv with complex S11, S21; hfss_chan5.py calls it)
       env RECORD_MD = other md file (tests), NO_ARCHIVE=1 keeps the files in results/
"""
import os, sys, json, glob, copy, shutil, time, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.optimize import least_squares, differential_evolution
from scipy.signal import find_peaks
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results'); sys.path.insert(0, HERE)
import theory5 as T, modes5 as M
RECORD = os.environ.get('RECORD_MD') or os.path.join(os.path.dirname(HERE), 'run_record.md')
F_CENTER = 6.71

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


def model_modes(P):
    """In-channel (C_in on the P1 side of n0, as the HFSS / measured S21 sees it) + matched-line kappa for reference."""
    r = M.extract(lambda z: T.branch_Y(z, P), P['Z0'], F_CENTER, span=1.2, Yenv=lambda z: T.n0_env_Y(z, P))
    rm = M.extract(lambda z: T.branch_Y(z, P), P['Z0'], F_CENTER, span=1.2)
    if r is not None and rm is not None:
        r['kp_matched'], r['J_matched'] = rm['kp'], rm['J']
    return r


def dips(f, s, prom=3):
    a = 20 * np.log10(np.abs(s)); i, _ = find_peaks(-a, prominence=prom, distance=max(1, int(0.002 / (f[1] - f[0]))))
    return [(round(float(f[k]), 4), round(float(a[k]), 1)) for k in i]


def fmt(r, k, nd):
    return '—' if (r is None or r.get(k) is None) else f'{r[k]:.{nd}f}'


def analyze(tag, note=''):
    S = pd.read_csv(os.path.join(RES, f'{tag}_S.csv'), index_col=0).map(complex)
    f = S.index.values.astype(float); s11, s21 = S.S11.values, S.S21.values
    meta_p = os.path.join(RES, f'{tag}_meta.json'); meta = json.load(open(meta_p)) if os.path.exists(meta_p) else {}
    P0 = T.params_from_layout()
    r_raw = model_modes(P0)
    xfit, Pf = theory_fit(P0, f, s11, s21); r_fit = model_modes(Pf)
    dfit, xd = M.direct_fit(f, s21)
    _, t21 = T.s21_theory(f, P0); t11f, t21f = T.s21_theory(f, Pf); t11, _ = T.s21_theory(f, P0)
    rms = lambda t: float(np.sqrt(np.mean(np.abs(t - s21) ** 2)))
    out = dict(tag=tag, date=time.strftime('%Y-%m-%d %H:%M'), note=note or meta.get('note', ''), hfss=meta.get('hfss', {}),
               cap_source=P0['cap_source'], caps_q3d=dict(cin=P0['cin'], cp=P0['cp'], cj=P0['cj']),
               caps_fit=dict(cin=Pf['cin'], cp=Pf['cp'], cj=Pf['cj']), fit_params=xfit,
               lengths=dict(L_filt=P0['L_filt'], L_res=P0['L_res'], L_filt_fit=Pf['L_filt'], L_res_fit=Pf['L_res']),
               dips_hfss=dips(f, s21), dips_theory=dips(f, t21), dips_fit=dips(f, t21f),
               modes=dict(theory_q3d=r_raw, theory_fit=r_fit, direct_fit=dfit),
               rms_S21=dict(theory_q3d=rms(t21), theory_fit=rms(t21f)))
    # ---- plot ----
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    db = lambda s: 20 * np.log10(np.abs(s))
    fig, ax = plt.subplots(3, 1, figsize=(10, 11))
    ax[0].plot(f, db(s21), 'k', lw=1.6, label='HFSS S21')
    ax[0].plot(f, db(t21).clip(-90), ':', label=f'theory, Q3D C ({P0["cap_source"]})')
    ax[0].plot(f, db(t21f).clip(-90), '--', label='theory fit by C')
    ax[0].set_ylabel('S21 (dB)'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3); ax[0].set_title(f'{tag}  {out["note"]}')
    ax[1].plot(f, db(s11), 'k', lw=1.6, label='HFSS S11'); ax[1].plot(f, db(t11), ':', label='theory, Q3D C')
    ax[1].plot(f, db(t11f), '--', label='theory fit by C'); ax[1].set_ylabel('S11 (dB)'); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    w = (f >= dfit['window'][0]) & (f <= dfit['window'][1])
    ax[2].plot(f[w], db(s21[w]), 'k', lw=1.6, label='HFSS S21')
    ax[2].plot(f[w], db(t21f[w]).clip(-90), '--', label='theory fit by C')
    ax[2].plot(f[w], db(np.conj(M.io_model(f[w], xd))).clip(-90), '-.', label='direct fit eq. 2.137')
    ax[2].set_xlabel('f (GHz)'); ax[2].set_ylabel('S21 (dB)'); ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)
    fig.tight_layout()
    # ---- archive ----
    stamp = time.strftime('%Y%m%d_%H%M')
    run_dir = RES if os.environ.get('NO_ARCHIVE') else os.path.join(RES, f'{stamp}_{tag}')
    os.makedirs(run_dir, exist_ok=True)
    if run_dir != RES:
        for p in glob.glob(os.path.join(RES, f'{tag}_*')):
            shutil.move(p, os.path.join(run_dir, os.path.basename(p)))
        for p in glob.glob(os.path.join(RES, 'q3d5_*_caps.json')) + glob.glob(os.path.join(RES, 'q3d5_part_*_cmat_fF.csv')):
            shutil.copy(p, os.path.join(run_dir, f'{tag}_' + os.path.basename(p)))
    png = os.path.join(run_dir, f'{tag}_analysis.png'); fig.savefig(png, dpi=120); plt.close(fig)
    json.dump(out, open(os.path.join(run_dir, f'{tag}_analysis.json'), 'w'), indent=1, default=float)
    try:
        rel = os.path.relpath(png, os.path.dirname(RECORD))
    except ValueError:                                   # record on another drive (tests)
        rel = png
    record_md(out, rel.replace('\\', '/'))
    return out, run_dir


def record_md(o, png_rel):
    c, cf = o['caps_q3d'], o['caps_fit']; R = o['modes']
    h = o.get('hfss', {})
    L = [f'\n## {o["tag"]}  ({o["date"]})\n', f'{o["note"]}\n' if o['note'] else '',
         f'HFSS: {", ".join(f"{k} {v}" for k, v in h.items())}\n' if h else '',
         f'![{o["tag"]}]({png_rel.replace(" ", "%20")})\n',
         f'**Capacitances** (fF; Q3D source: `{o["cap_source"]}`)\n',
         '| | C_in | C_p | C_j |', '|---|---|---|---|',
         f'| series, Q3D | {c["cin"][0]:.2f} | {c["cp"][0]:.2f} | {c["cj"][0]:.2f} |',
         f'| series, theory fit | {cf["cin"][0]:.2f} | {cf["cp"][0]:.2f} | {cf["cj"][0]:.2f} |',
         f'| shunt n0/line side, filter/node side | {c["cin"][1]:.1f} / {c["cin"][2]:.1f} | {c["cp"][1]:.1f} / {c["cp"][2]:.1f} | {c["cj"][1]:.1f} / {c["cj"][2]:.1f} |',
         '', f'Fit by C: dL_filter {o["fit_params"]["dLf"]:+.4f} mm, dL_res {o["fit_params"]["dLr"]:+.4f} mm, '
             f'v scale {o["fit_params"]["kv"]:.4f}; rms |ΔS21| Q3D {o["rms_S21"]["theory_q3d"]:.4f} → fit {o["rms_S21"]["theory_fit"]:.4f}\n',
         f'S21 dips: HFSS {o["dips_hfss"]}; theory (Q3D C) {o["dips_theory"]}; theory fit {o["dips_fit"]}\n',
         '| | ω_p/2π (GHz) | ω_r/2π (GHz) | J/2π (MHz) | κ_p/2π (MHz) |', '|---|---|---|---|---|']
    for k, lab in (('theory_q3d', 'theory, Q3D C'), ('theory_fit', 'theory fit by C'), ('direct_fit', 'direct fit eq. 2.137')):
        r = R[k]
        L.append(f'| {lab} | {fmt(r, "wp", 4)} | {fmt(r, "wr", 4)} | {fmt(r, "J", 2)} | {fmt(r, "kp", 2)} |')
    km = [f'{lab} {R[k]["kp_matched"]:.1f}' for k, lab in (('theory_q3d', 'Q3D'), ('theory_fit', 'fit'))
          if R[k] and 'kp_matched' in R[k]]
    if km:
        L.append(f'\nκ_p of the theory rows is in-channel (C_in on the P1 side of n0, like the S21 fit); '
                 f'κ_p on a matched line (Wallraff Table 2.4 definition): {", ".join(km)} MHz')
    d = R['direct_fit']
    L.append(f'\nDirect fit: window {d["window"][0]:.3f}–{d["window"][1]:.3f} GHz, rms/|S21| {d["rms_rel"]:.3f}, '
             f'A {d["A"]:.3f}, φ {d["phi"]:+.3f} rad\n')
    if not os.path.exists(RECORD):
        open(RECORD, 'w', encoding='utf-8').write('# Design5 run record\n')
    open(RECORD, 'a', encoding='utf-8').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    o, d = analyze(sys.argv[1], ' '.join(sys.argv[2:]))
    print(json.dumps(o['modes'], indent=1, default=float)); print('run folder', d, '\nrecord', RECORD)
