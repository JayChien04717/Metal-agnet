"""HFSS driven modal of the WHOLE Design5 chip (chip6_layout.py: 6 resonator + Purcell channels on one straight feed
line, ports in the middle of the left / right edges), then the comparison (only this, Jay 2026-11-09):
  1. capacitances: Q3D and fitted to HFSS,
  2. HFSS S21 / S11 (results/<time>_<tag>/<tag>_S.csv),
  3. one plot: raw HFSS S21 / S11, the chip circuit model with the Q3D C (no fit), the model fitted by C
     (12 lengths + C_in, C_p, C_j scales) -> how far the capacitance model is from HFSS.
Appended to Design5/run_record.md.

Usage: python hfss_chip6.py <tag> [f_start f_stop count type]     (env SOL_F, MAXP, BASIS, RUN_NOTE)
       python hfss_chip6.py --analyze <tag>                      (analysis only, results/<tag>_S.csv)
"""
import os, sys, json, time, glob, shutil, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.signal import find_peaks
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES = os.path.join(HERE, 'results')   # simulation/ (this script: simulation/<kind>/)
import sys as _sys; _sys.path[:0] = [os.path.join(HERE, _k) for _k in ('driven', 'eigenmode', 'q3d')] + [HERE]
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, HERE)
RECORD = os.environ.get('RECORD_MD') or os.path.join(os.path.dirname(HERE), 'run_record.md')


def run_hfss(tag, f0, f1, count, stype):
    from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
    import chip6_layout as C6, hfss_dirs
    temp_dirs = hfss_dirs.use_d_drive()      # HFSS solver scratch on D: (C: fills up -> COM error at Analyze)
    import save_aedt; save_aedt.use_project(save_aedt.project_file(HERE))   # one AEDT project per Design, never closed
    d = C6.build_chip()
    setup = dict(freq_ghz=float(os.environ.get('SOL_F', 6.80)), max_delta_s=0.003, max_passes=int(os.environ.get('MAXP', 14)),
                 min_passes=6, min_converged=2, pct_refinement=30, basis_order=int(os.environ.get('BASIS', 2)))
    sim = ScatteringImpedanceSim(d, 'hfss')
    for k, v in setup.items():
        sim.setup[k] = v
    sim.setup.sweep_setup.update(name='Sweep', start_ghz=f0, stop_ghz=f1, count=count, type=stype, save_fields=False)
    # everything except the launch pads and the pin-only ResEnds; ports on the feed-line ends at the pads
    comps = [n for n, c in d.components.items() if type(c).__name__ not in ('ResEnd', 'LaunchpadWirebond')]
    opens = [(f'res_a_ch{ch}', 'start') for ch, _, _ in C6.ORDER]          # resonator open ends (else shorted)
    nl = len(C6.ORDER)
    t0 = time.time()
    sim.run(name=save_aedt.design_name(tag), components=comps, open_terminations=opens, port_list=[('feed_in', 'start', 50), (f'fl_{nl}', 'end', 50)])
    conv_info = {}
    try:
        conv = sim.renderer.pinfo.setup.get_convergence()[0]
        conv.to_csv(os.path.join(RES, f'{tag}_convergence.csv')); print(conv.tail(4).to_string())
        conv_info = dict(passes=int(len(conv)), last_dS=float(conv.iloc[-1, -1]))
    except Exception as e:
        print('convergence not read:', repr(e)[:150])
    _, _, S = sim.renderer.get_params(['S11', 'S21'])
    S.index.name = 'f_GHz'; S.to_csv(os.path.join(RES, f'{tag}_S.csv'))
    minutes = round((time.time() - t0) / 60, 1)
    json.dump(dict(tag=tag, note=os.environ.get('RUN_NOTE', ''),
                   hfss=dict(sweep=f'{f0}-{f1} GHz, {count} pts, {stype}', **setup, **conv_info, minutes=minutes,
                             temp_dir=temp_dirs['after'])),
              open(os.path.join(RES, f'{tag}_meta.json'), 'w'), indent=1)
    print(f'=== {tag} HFSS done ({minutes} min) ===')
    try:                                               # keep the AEDT project + solved results with the run
        import save_aedt; save_aedt.snapshot(os.path.join(RES, f'{tag}_aedt'), tag)
    except Exception as e:
        print('AEDT project not saved:', repr(e)[:150])
    try: sim.close()
    except Exception: pass


FIT_G = ['s_in', 's_p', 's_j']                             # global: series C scales (same cap geometry on every channel)
G_BOUNDS = [(0.7, 1.3), (0.7, 1.3), (0.3, 3.0)]
PRIOR = dict(s_in=0.15, s_p=0.15, s_j=0.5)                 # soft priors toward the Q3D values
# no v scale: f ~ v / L, so a v scale is degenerate with the 12 length corrections (synthetic test: the fit
# drifted to v x 0.98 and missed every length); the lengths absorb a v error


def chip_model(tag_or_none=None):
    """Chip circuit model with the drawn geometry: returns (P0, chans [(ch, P)], segs, L)."""
    import theory5 as T, chip6_layout as C6, ftf_layout as f3, copy
    d = C6.build_chip(); P0 = T.params_from_layout()
    L = {c['ch']: c for c in json.load(open(C6.LEN_JSON))['channels']}
    pl = lambda n: f3.path_length(d, n); nl = len(C6.ORDER)
    P0['l_feed_in'] = pl('feed_in') + 0.04
    segs = [pl('fl_0') + 0.04 + 0.15] + [pl(f'fl_{k}') + 0.3 for k in range(1, nl)] + [pl(f'fl_{nl}') + 0.15]
    chans = []
    for ch, _, _ in C6.ORDER:
        P = copy.deepcopy(P0); P['L_filt'], P['L_res'] = L[ch]['L_filter'], L[ch]['L_res']; chans.append((ch, P))
    return P0, chans, segs, L


def apply_x(P0, chans, x):
    """x = [dLf_k, dLr_k for each channel in feed-line order] + [s_in, s_p, s_j] -> (P0', [P'])."""
    n = len(chans); g = dict(zip(FIT_G, x[2 * n:]))
    Q0 = dict(P0); Q0['cin'] = (P0['cin'][0] * g['s_in'],) + tuple(P0['cin'][1:])
    out = []
    for k, (_, P) in enumerate(chans):
        Q = dict(P)
        Q['L_filt'] = P['L_filt'] + x[2 * k]; Q['L_res'] = P['L_res'] + x[2 * k + 1]
        Q['cp'] = (P['cp'][0] * g['s_p'],) + tuple(P['cp'][1:]); Q['cj'] = (P['cj'][0] * g['s_j'],) + tuple(P['cj'][1:])
        out.append(Q)
    return Q0, out


def chip_fit(f, s11, s21, P0, chans, segs, L):
    """Theory S21 model fitted by C to the HFSS S21 (+ 0.3 S11).
    Rounds (x3): (a) each channel alone in its window, differential evolution on its dLf, dLr (others fixed);
                 (b) least squares on the shared C_in, C_p, C_j scales over the whole band.
    Final: least squares on all 12 lengths + 3 C scales together, soft priors on the C scales."""
    import chip6_theory as CT
    from scipy.optimize import least_squares, differential_evolution
    n = len(chans); x = np.r_[np.zeros(2 * n), np.ones(len(FIT_G))]
    def res(xx, m=slice(None)):
        Q0, Qs = apply_x(P0, chans, xx); t11, t21 = CT.chip_s21(f[m], Q0, Qs, segs, s11=True)
        return np.r_[(t21 - s21[m]).real, (t21 - s21[m]).imag, 0.3 * (t11 - s11[m]).real, 0.3 * (t11 - s11[m]).imag]
    def prior(xx):
        g = dict(zip(FIT_G, xx[2 * n:])); return 0.02 * np.array([(g[k] - 1) / s for k, s in PRIOR.items()])
    sc = 30 / np.sqrt(len(f))
    win = {k: (f >= L[ch]['f_r'] - 0.09) & (f <= L[ch]['f_r'] + 0.10) for k, (ch, _) in enumerate(chans)}
    for rnd in range(3):
        for k in range(n):
            def cost(p, k=k):
                xx = x.copy(); xx[2 * k:2 * k + 2] = p; return float(np.sum(res(xx, win[k]) ** 2))
            lim = 0.3 if rnd == 0 else 0.08
            b = [(x[2 * k + i] - lim, x[2 * k + i] + lim) for i in range(2)]
            x[2 * k:2 * k + 2] = differential_evolution(cost, b, seed=1, maxiter=40, popsize=12, tol=1e-10, polish=True).x
        def res_g(g):
            xx = x.copy(); xx[2 * n:] = g; return np.r_[res(xx) * sc, prior(xx)]
        x[2 * n:] = least_squares(res_g, x[2 * n:], bounds=tuple(zip(*G_BOUNDS))).x
    lo = [-0.5] * (2 * n) + [b[0] for b in G_BOUNDS]; hi = [0.5] * (2 * n) + [b[1] for b in G_BOUNDS]
    sol = least_squares(lambda xx: np.r_[res(xx) * sc, prior(xx)], np.clip(x, lo, hi), bounds=(lo, hi),
                        x_scale=[0.02] * (2 * n) + [0.05, 0.05, 0.1])
    return sol.x


def find_run(tag):
    """results/<tag>_S.csv of a run just finished, or the newest run folder (results/<date>/<HHMM>_<tag>/)."""
    import runstore
    p = os.path.join(RES, f'{tag}_S.csv')
    if os.path.exists(p):
        return p, None
    hits = [d for d in runstore.find_run_dirs(tag, RES) if os.path.exists(os.path.join(d, f'{tag}_S.csv'))]
    if not hits:
        raise FileNotFoundError(f'no S-parameters for {tag}')
    return os.path.join(hits[-1], f'{tag}_S.csv'), hits[-1]


def analyze_chip(tag):
    """Per run: (1) capacitances, (2) HFSS S21 / S11, (3) one plot: raw HFSS, theory with the Q3D C (no fit),
    theory fitted by C -> how far the capacitance model is from HFSS."""
    import chip6_theory as CT, chip6_layout as C6
    s_path, run_dir = find_run(tag)
    S = pd.read_csv(s_path, index_col=0).map(complex)
    f = S.index.values.astype(float); s11, s21 = S.S11.values, S.S21.values
    meta_p = os.path.join(os.path.dirname(s_path), f'{tag}_meta.json'); meta = json.load(open(meta_p)) if os.path.exists(meta_p) else {}
    P0, chans, segs, L = chip_model()
    t11, t21 = CT.chip_s21(f, P0, [P for _, P in chans], segs, s11=True)            # Q3D C, no fit
    x = chip_fit(f, s11, s21, P0, chans, segs, L)
    Q0, Qs = apply_x(P0, chans, x)
    u11, u21 = CT.chip_s21(f, Q0, Qs, segs, s11=True)                               # fitted by C
    n = len(chans); g = dict(zip(FIT_G, x[2 * n:]))
    db = lambda s: 20 * np.log10(np.abs(s))
    def dips(s, lo, hi):
        m = (f >= lo) & (f <= hi); a = db(s[m]); i, _ = find_peaks(-a, prominence=3, distance=max(1, int(0.002 / (f[1] - f[0]))))
        i = sorted(sorted(i, key=lambda k: a[k])[:2], key=lambda k: f[m][k])
        return [round(float(f[m][k]), 4) for k in i]
    rows = []
    for k, (ch, P) in enumerate(chans):
        fr = L[ch]['f_r']; lo, hi = fr - 0.09, fr + 0.10
        rows.append(dict(ch=ch, f_design=fr, dips_hfss=dips(s21, lo, hi), dips_q3d=dips(t21, lo, hi), dips_fit=dips(u21, lo, hi),
                         dL_filter=float(x[2 * k]), dL_res=float(x[2 * k + 1])))
    rows.sort(key=lambda r: r['ch'])
    rms = lambda t: float(np.sqrt(np.mean(np.abs(t - s21) ** 2)))
    out = dict(tag=tag, date=time.strftime('%Y-%m-%d %H:%M'), note=meta.get('note', ''), hfss=meta.get('hfss', {}),
               cap_source=P0['cap_source'], caps_q3d=dict(cin=P0['cin'], cp=P0['cp'], cj=P0['cj']),
               caps_fit=dict(cin=Q0['cin'], cp=Qs[0]['cp'], cj=Qs[0]['cj']), channels=rows,
               rms_S21=dict(q3d=rms(t21), fit=rms(u21)))
    # ---- one figure: S21 and S11 over the band, then S21 per channel ----
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(14, 14)); gs = fig.add_gridspec(4, 3)
    sty = [('HFSS (raw)', 'k', '-', 1.4), ('theory, Q3D C (no fit)', 'C0', ':', 1.2), ('theory fitted by C', 'C3', '--', 1.0)]
    for row, (name, curves) in enumerate((('S21', (s21, t21, u21)), ('S11', (s11, t11, u11)))):
        ax = fig.add_subplot(gs[row, :])
        for (lab, c, ls, lw), s in zip(sty, curves):
            ax.plot(f, db(s).clip(-90), color=c, ls=ls, lw=lw, label=lab)
        ax.set_ylabel(f'{name} (dB)'); ax.grid(alpha=.3); ax.legend(fontsize=8, loc='lower left')
        if row == 0:
            ax.set_title(f'{tag}  {out["note"]}', fontsize=9)
    for k, r in enumerate(rows):
        ax = fig.add_subplot(gs[2 + k // 3, k % 3]); fr = r['f_design']; m = (f >= fr - 0.09) & (f <= fr + 0.10)
        for (lab, c, ls, lw), s in zip(sty, (s21, t21, u21)):
            ax.plot(f[m], db(s[m]).clip(-90), color=c, ls=ls, lw=lw, label=lab)
        ax.set_title(f'ch{r["ch"]} S21 (design f_r {fr:.2f} GHz)', fontsize=9); ax.grid(alpha=.3)
        ax.set_xlabel('f (GHz)')
    fig.tight_layout()
    import runstore
    if run_dir is None:                                   # just finished: archive into results/<date>/<HHMM>_<tag>/
        run_dir = runstore.new_run_dir(tag, RES)
        for p in glob.glob(os.path.join(RES, f'{tag}_*')) + glob.glob(os.path.join(RES, f'{tag}.log')):
            try:
                shutil.move(p, os.path.join(run_dir, os.path.basename(p)))
            except (PermissionError, OSError):           # e.g. the run log still open by the shell redirection
                shutil.copy(p, os.path.join(run_dir, os.path.basename(p)))
        for p in [os.path.join(RES, 'q3d5_full_caps.json'), C6.LEN_JSON]:       # inputs used by this run
            if os.path.exists(p): shutil.copy(p, os.path.join(run_dir, f'{tag}_' + os.path.basename(p)))
    st = runstore.stamp()                                 # every (re-)analysis kept, time-stamped
    png = os.path.join(run_dir, f'{tag}_compare_{st}.png'); fig.savefig(png, dpi=110); plt.close(fig)
    json.dump(out, open(os.path.join(run_dir, f'{tag}_compare_{st}.json'), 'w'), indent=1, default=float)
    try:
        rel = os.path.relpath(png, os.path.dirname(RECORD))
    except ValueError:                                   # record on another drive (tests)
        rel = png
    record(out, rel.replace('\\', '/'))
    return out, run_dir


def record(o, png_rel):
    c, cf = o['caps_q3d'], o['caps_fit']; h = o.get('hfss', {})
    pct = lambda a, b: f'{100 * (b / a - 1):+.1f} %'
    Lr = [f'\n## {o["tag"]} — full chip, 6 channels  ({o["date"]})\n', f'{o["note"]}\n' if o['note'] else '',
          f'HFSS: {", ".join(f"{k} {v}" for k, v in h.items())}\n' if h else '',
          f'![{o["tag"]}]({png_rel.replace(" ", "%20")})\n',
          f'**1. Capacitances** (series, fF; Q3D source `{o["cap_source"]}`)\n',
          '| | C_in | C_p | C_j |', '|---|---|---|---|',
          f'| Q3D | {c["cin"][0]:.2f} | {c["cp"][0]:.2f} | {c["cj"][0]:.2f} |',
          f'| fitted to HFSS S21 | {cf["cin"][0]:.2f} | {cf["cp"][0]:.2f} | {cf["cj"][0]:.2f} |',
          f'| difference | {pct(c["cin"][0], cf["cin"][0])} | {pct(c["cp"][0], cf["cp"][0])} | {pct(c["cj"][0], cf["cj"][0])} |',
          '', f'rms |S21 − HFSS|: Q3D model {o["rms_S21"]["q3d"]:.4f}, fitted {o["rms_S21"]["fit"]:.4f}\n',
          '**2. S21 dips** (GHz) and fitted length corrections (mm)\n',
          '| ch | design f_r | HFSS | theory, Q3D C | theory fitted | ΔL filter | ΔL resonator |', '|---|---|---|---|---|---|---|']
    for r in o['channels']:
        j = lambda v: ', '.join(f'{x:.4f}' for x in v)
        Lr.append(f'| {r["ch"]} | {r["f_design"]:.2f} | {j(r["dips_hfss"])} | {j(r["dips_q3d"])} | {j(r["dips_fit"])} | '
                  f'{r["dL_filter"]:+.4f} | {r["dL_res"]:+.4f} |')
    open(RECORD, 'a', encoding='utf-8').write('\n'.join(Lr) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == '--analyze':
        o, rd = analyze_chip(sys.argv[2])
        print('caps Q3D', {k: round(v[0], 2) for k, v in o['caps_q3d'].items()}, '-> fit', {k: round(v[0], 2) for k, v in o['caps_fit'].items()}, 'rms', o['rms_S21'])
        for r in o['channels']:
            print(f'ch{r["ch"]}: HFSS {r["dips_hfss"]}, Q3D model {r["dips_q3d"]}, fit {r["dips_fit"]}, dLf {r["dL_filter"]:+.4f} dLr {r["dL_res"]:+.4f}')
        print('run folder', rd)
        sys.exit()
    tag = sys.argv[1]
    f0, f1 = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (6.1, 7.5)
    count = int(sys.argv[4]) if len(sys.argv) > 4 else int(round((f1 - f0) * 1000 / 0.25)) + 1
    stype = sys.argv[5] if len(sys.argv) > 5 else 'Fast'
    run_hfss(tag, f0, f1, count, stype)
    o, rd = analyze_chip(tag)
    for r in o['channels']:
        print(f'ch{r["ch"]} design {r["f_design"]:.2f}: HFSS dips {r["dips_hfss"]}, Q3D model {r["dips_q3d"]}, fit {r["dips_fit"]}')
    print('run folder', rd)
