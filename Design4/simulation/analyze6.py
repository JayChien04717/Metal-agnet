"""Analysis of a Design6 verification run: do the Q3D capacitances give the couplings that HFSS S21 shows?

Per run (appended to Design6/run_record.md):
  1. capacitances: Q3D (this geometry, q3d6.py) and fitted to the HFSS S21,
  2. HFSS S21 / S11 (<tag>_S.csv),
  3. couplings: from the Q3D C (no fit) vs from the fitted model
       filter : omega_p, kappa_p          (Purcell <-> feed line, C_p)
       full   : omega_r, omega_p, J, kappa_p   (J: resonator <-> Purcell, C_j)
  4. one plot: raw HFSS S21 / S11, model with the Q3D C, fitted model.
Fit (complex S21 + 0.3 |S11|, differential evolution then least squares; free phase phi and delay tau absorb the
reference planes):
  filter : dL_filter, C_p scale
  full   : dL_res, C_j scale; dL_filter and C_p fixed from the latest filter analysis (already verified there)

Usage: python analyze6.py <tag>        (results/<tag>_S.csv just written, or the newest results/<date>/<HHMM>_<tag>/)
"""
import os, sys, json, glob, shutil, time, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.optimize import least_squares, differential_evolution
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results')
HFSS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE); sys.path.insert(1, os.path.join(os.path.dirname(HERE), 'pattern design'))
sys.path.insert(2, os.path.join(HFSS, 'Design5', 'simulation'))
import verify_model as VM, runstore
RECORD = os.environ.get('RECORD_MD') or os.path.join(os.path.dirname(HERE), 'run_record.md')


def find_run(tag):
    p = os.path.join(RES, f'{tag}_S.csv')
    if os.path.exists(p):
        return p, None
    hits = [d for d in runstore.find_run_dirs(tag, RES) if os.path.exists(os.path.join(d, f'{tag}_S.csv'))]
    if not hits:
        raise FileNotFoundError(tag)
    return os.path.join(hits[-1], f'{tag}_S.csv'), hits[-1]


def q3d_params(variant):
    """Model parameters with the Design6 Q3D series C (falls back to the Design5 values if q3d6 has not run)."""
    P = VM.params(); src = 'Design5 whole-channel Q3D'
    pf = os.path.join(RES, 'q3d6_full_caps.json'); pv = os.path.join(RES, f'q3d6_{variant}_caps.json')
    if os.path.exists(pf):
        q = json.load(open(pf)); P['cp'] = (q['series']['cp'],) + tuple(P['cp'][1:]); P['cj'] = (q['series']['cj'],) + tuple(P['cj'][1:])
        src = f'Design6 Q3D ({q["run_dir"]})'
    if variant == 'filter' and os.path.exists(pv):          # floating hammer plate: its own C_j and plate-to-ground C
        q = json.load(open(pv)); P['cj'] = (q['series']['cj'], P['cj'][1], q['ground_fF']['resonator'])
        src += f' + {q["run_dir"]}'
    P['cap_source'] = src
    return P


def feed_lengths(variant):
    import verify_layout as V
    d = V.build(variant); pl = lambda n: V.f3.path_length(d, n)
    return pl('feed_in') + 0.15, pl('feed_out') + 0.15


def model(f, P, res, lin, lout, phi=0.0, tau=0.0):
    w = 2 * np.pi * f * 1e9; b = w / P['v']
    Y = VM.branch_Y(f, P, res); den = 2 + P['Z0'] * Y; ph = np.exp(1j * (phi - 2 * np.pi * f * tau))
    return np.exp(-2j * b * lin * 1e-3) * (-P['Z0'] * Y / den), ph * np.exp(-1j * b * (lin + lout) * 1e-3) * 2 / den


def with_x(P, names, x):
    Q = dict(P); p = dict(zip(names, x))
    if 'dLf' in p: Q['L_filt'] = P['L_filt'] + p['dLf']
    if 'dLr' in p: Q['L_res'] = P['L_res'] + p['dLr']
    if 's_p' in p: Q['cp'] = (P['cp'][0] * p['s_p'],) + tuple(P['cp'][1:])
    if 's_j' in p: Q['cj'] = (P['cj'][0] * p['s_j'],) + tuple(P['cj'][1:])
    return Q, p.get('phi', 0.0), p.get('tau', 0.0)


def fit(f, s11, s21, P, res, lin, lout, names, bounds):
    """Stage A: physical parameters on |S21|, |S11| only (phase-free; dense grid over dLr x s_j for the full
    channel, where the two hybridised modes make the cost multi-modal), stage B: phase phi and delay tau,
    stage C: everything on complex S21 + |S11| by least squares."""
    phys = [n for n in names if n not in ('phi', 'tau')]; ip = [names.index(n) for n in phys]
    def r(x):
        Q, phi, tau = with_x(P, names, x); t11, t21 = model(f, Q, res, lin, lout, phi, tau)
        return np.r_[(t21 - s21).real, (t21 - s21).imag, 0.3 * (np.abs(t11) - np.abs(s11))]
    def r_mag(xp):
        Q, _, _ = with_x(P, phys, xp); t11, t21 = model(f, Q, res, lin, lout)
        return np.r_[np.abs(t21) - np.abs(s21), 0.3 * (np.abs(t11) - np.abs(s11))]
    pb = [bounds[i] for i in ip]
    if not phys:
        x0 = np.zeros(0)
    elif 'dLr' in phys and 's_j' in phys and len(phys) == 2:
        grid = [(a, b) for a in np.arange(pb[0][0], pb[0][1] + 1e-9, 0.004) for b in np.geomspace(pb[1][0], pb[1][1], 25)]
        costs = [float(np.sum(r_mag(np.array(g)) ** 2)) for g in grid]
        x0 = np.array(grid[int(np.argmin(costs))])
    else:
        x0 = differential_evolution(lambda xp: float(np.sum(r_mag(xp) ** 2)), pb, seed=1, maxiter=150, popsize=20, tol=1e-10, polish=False).x
    xp = least_squares(r_mag, x0, bounds=tuple(zip(*pb))).x if phys else x0
    x = np.zeros(len(names)); x[ip] = xp
    iph = [names.index(n) for n in ('phi', 'tau') if n in names]
    if iph:
        def r_ph(q):
            xx = x.copy(); xx[iph] = q; return r(xx)
        q = differential_evolution(lambda q: float(np.sum(r_ph(q) ** 2)), [bounds[i] for i in iph], seed=1, maxiter=100, popsize=20, polish=False).x
        x[iph] = q
    lo, hi = zip(*bounds)
    sol = least_squares(r, np.clip(x, lo, hi), bounds=(lo, hi))
    return dict(zip(names, map(float, sol.x)))


def couplings(P, res):
    pr = VM.predict(P)
    if not res:
        return pr['filter']
    return {k: (pr['full'][k] if pr['full'] else float('nan')) for k in ('wr', 'wp', 'J', 'kp')}


def analyze(tag):
    s_path, run_dir = find_run(tag)
    S = pd.read_csv(s_path, index_col=0).map(complex)
    f = S.index.values.astype(float); s11, s21 = S.S11.values, S.S21.values
    meta_p = os.path.join(os.path.dirname(s_path), f'{tag}_meta.json'); meta = json.load(open(meta_p)) if os.path.exists(meta_p) else {}
    variant = meta.get('variant', 'full' if 'full' in tag else 'filter'); res = variant == 'full'
    P = q3d_params(variant); lin, lout = feed_lengths(variant)
    fixed = {}
    if res:
        prev = sorted(glob.glob(os.path.join(RES, '*', '*', '*_compare_*.json')), key=os.path.getmtime)
        prev = [json.load(open(p)) for p in prev]; prev = [o for o in prev if o.get('variant') == 'filter']
        if prev:
            fixed = dict(dLf=prev[-1]['fit']['dLf'], s_p=prev[-1]['fit']['s_p'], from_run=prev[-1]['tag'])
            P = dict(P, L_filt=P['L_filt'] + fixed['dLf'], cp=(P['cp'][0] * fixed['s_p'],) + tuple(P['cp'][1:]))
            names, bounds = ['dLr', 's_j', 'phi', 'tau'], [(-0.3, 0.3), (0.3, 3.0), (-np.pi, np.pi), (-0.05, 0.05)]
        else:
            names, bounds = ['dLf', 'dLr', 's_p', 's_j', 'phi', 'tau'], [(-0.3, 0.3), (-0.3, 0.3), (0.6, 1.6), (0.3, 3.0), (-np.pi, np.pi), (-0.05, 0.05)]
    else:
        names, bounds = ['dLf', 's_p', 'phi', 'tau'], [(-0.3, 0.3), (0.6, 1.6), (-np.pi, np.pi), (-0.05, 0.05)]
    # model with the Q3D C: best phase / delay only (no physical parameter fitted)
    ph = fit(f, s11, s21, P, res, lin, lout, ['phi', 'tau'], [(-np.pi, np.pi), (-0.05, 0.05)])
    t11, t21 = model(f, P, res, lin, lout, ph['phi'], ph['tau'])
    x = fit(f, s11, s21, P, res, lin, lout, names, bounds)
    Pf, phi, tau = with_x(P, names, [x[n] for n in names]); u11, u21 = model(f, Pf, res, lin, lout, phi, tau)
    rms = lambda t: float(np.sqrt(np.mean(np.abs(t - s21) ** 2)))
    cq = dict(cp=P['cp'][0] / fixed.get('s_p', 1.0), cj=P['cj'][0])           # Q3D series C
    # couplings at the HFSS-fitted lengths (J and kappa_p depend on the frequencies / detuning too), once with the
    # Q3D C and once with the fitted C -> the difference is the capacitor model alone
    Pq = dict(Pf, cp=(cq['cp'],) + tuple(Pf['cp'][1:]), cj=(cq['cj'],) + tuple(Pf['cj'][1:]))
    out = dict(tag=tag, variant=variant, date=time.strftime('%Y-%m-%d %H:%M'), note=meta.get('note', ''), hfss=meta.get('hfss', {}),
               cap_source=P['cap_source'], fixed_from_filter=fixed,
               caps_q3d=cq, caps_fit=dict(cp=Pf['cp'][0], cj=Pf['cj'][0]), fit=x,
               couplings_design=couplings(dict(P, L_filt=P['L_filt'] - fixed.get('dLf', 0.0), cp=(cq['cp'],) + tuple(P['cp'][1:])), res),
               couplings_q3d=couplings(Pq, res), couplings_fit=couplings(Pf, res), rms_S21=dict(q3d=rms(t21), fit=rms(u21)))
    # ---- plot ----
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    db = lambda s: 20 * np.log10(np.abs(s))
    fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    sty = [('HFSS (raw)', 'k', '-', 1.6), ('model, Q3D C (no fit)', 'C0', ':', 1.4), ('model fitted by C', 'C3', '--', 1.1)]
    for a, nm, cs in ((ax[0], 'S21', (s21, t21, u21)), (ax[1], 'S11', (s11, t11, u11))):
        for (lab, c, ls, lw), s in zip(sty, cs):
            a.plot(f, db(s).clip(-90), color=c, ls=ls, lw=lw, label=lab)
        a.set_ylabel(f'{nm} (dB)'); a.grid(alpha=.3); a.legend(fontsize=8)
    cq, cf = out['couplings_q3d'], out['couplings_fit']
    txt = ', '.join(f'{k} {cq.get(k, float("nan")):.4g}→{cf.get(k, float("nan")):.4g}' for k in cq)
    ax[0].set_title(f'{tag} ({variant})  C_p {out["caps_q3d"]["cp"]:.2f}→{out["caps_fit"]["cp"]:.2f} fF, '
                    f'C_j {out["caps_q3d"]["cj"]:.2f}→{out["caps_fit"]["cj"]:.2f} fF (Q3D→fit)\n{txt}', fontsize=8)
    ax[1].set_xlabel('f (GHz)'); fig.tight_layout()
    # ---- store ----
    if run_dir is None:
        run_dir = runstore.new_run_dir(tag, RES)
        for p in glob.glob(os.path.join(RES, f'{tag}_*')) + glob.glob(os.path.join(RES, f'{tag}.log')):
            try: shutil.move(p, os.path.join(run_dir, os.path.basename(p)))
            except (PermissionError, OSError): shutil.copy(p, os.path.join(run_dir, os.path.basename(p)))
        for p in glob.glob(os.path.join(RES, 'q3d6_*_caps.json')):
            shutil.copy(p, os.path.join(run_dir, f'{tag}_' + os.path.basename(p)))
    st = runstore.stamp()
    png = os.path.join(run_dir, f'{tag}_compare_{st}.png'); fig.savefig(png, dpi=120); plt.close(fig)
    json.dump(out, open(os.path.join(run_dir, f'{tag}_compare_{st}.json'), 'w'), indent=1, default=float)
    try: rel = os.path.relpath(png, os.path.dirname(RECORD))
    except ValueError: rel = png
    record(out, rel.replace('\\', '/'))
    return out, run_dir


def record(o, png_rel):
    cq, cf = o['caps_q3d'], o['caps_fit']; h = o.get('hfss', {})
    pct = lambda a, b: f'{100 * (b / a - 1):+.1f} %'
    L = [f'\n## {o["tag"]} — {o["variant"]}  ({o["date"]})\n', f'{o["note"]}\n' if o['note'] else '',
         f'HFSS: {", ".join(f"{k} {v}" for k, v in h.items())}\n' if h else '',
         f'![{o["tag"]}]({png_rel.replace(" ", "%20")})\n',
         f'**Capacitances** (series, fF; Q3D: {o["cap_source"]})\n', '| | C_p | C_j |', '|---|---|---|',
         f'| Q3D | {cq["cp"]:.2f} | {cq["cj"]:.2f} |', f'| fitted to HFSS S21 | {cf["cp"]:.2f} | {cf["cj"]:.2f} |',
         f'| difference | {pct(cq["cp"], cf["cp"])} | {pct(cq["cj"], cf["cj"])} |', '']
    if o['fixed_from_filter']:
        L.append(f'C_p and the filter length fixed from `{o["fixed_from_filter"]["from_run"]}` (fit there).\n')
    names = dict(wr='ω_r (GHz)', wp='ω_p (GHz)', J='J (MHz)', kp='κ_p (MHz)')
    ks = list(o['couplings_q3d'])
    row = lambda lab, c: f'| {lab} | ' + ' | '.join(f'{c[k]:.4f}' if k in ('wr', 'wp') else f'{c[k]:.1f}' for k in ks) + ' |'
    L += ['**Couplings** (rows 2 and 3 at the HFSS-fitted lengths: their difference is the capacitor model alone)\n',
          '| | ' + ' | '.join(names[k] for k in ks) + ' |', '|---' * (len(ks) + 1) + '|',
          row('design (drawn lengths, Q3D C)', o['couplings_design']),
          row('fitted lengths, Q3D C', o['couplings_q3d']), row('fitted lengths, fitted C = HFSS', o['couplings_fit']),
          '', f'Fit: {", ".join(f"{k} {v:+.4f}" for k, v in o["fit"].items())}; rms |S21 − HFSS|: Q3D model {o["rms_S21"]["q3d"]:.4f}, fitted {o["rms_S21"]["fit"]:.4f}\n']
    if not os.path.exists(RECORD):
        open(RECORD, 'w', encoding='utf-8').write('# Design6 run record — capacitor model vs HFSS couplings\n\nMethod: `README.md`.\n')
    open(RECORD, 'a', encoding='utf-8').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    o, d = analyze(sys.argv[1])
    print(json.dumps({k: o[k] for k in ('caps_q3d', 'caps_fit', 'couplings_q3d', 'couplings_fit', 'rms_S21')}, indent=1, default=float))
    print('run folder', d)
