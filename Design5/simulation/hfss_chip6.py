"""HFSS driven modal of the WHOLE Design5 chip (chip6_layout.py: 6 resonator + Purcell channels on one straight feed
line, ports in the middle of the left / right edges), then the chip analysis:
  * HFSS S21 / S11 vs the chip-level circuit model (chip6_theory.chip_s21, Q3D C, drawn lengths, no fitting),
  * dips per channel: HFSS vs model,
  * direct eq. 2.137 fit of S21 in a window around each channel -> omega_p, omega_r, J, kappa_p (in the chip),
  * model omega_p, omega_r, J, kappa_p of each channel (matched-line definition) for comparison,
appended to Design5/run_record.md; files in results/<time>_<tag>/.

Usage: python hfss_chip6.py <tag> [f_start f_stop count type]     (env SOL_F, MAXP, BASIS, RUN_NOTE)
       python hfss_chip6.py --analyze <tag>                      (analysis only, results/<tag>_S.csv)
"""
import os, sys, json, time, glob, shutil, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.signal import find_peaks
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results')
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, HERE)
RECORD = os.environ.get('RECORD_MD') or os.path.join(os.path.dirname(HERE), 'run_record.md')


def run_hfss(tag, f0, f1, count, stype):
    from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
    import chip6_layout as C6, hfss_dirs
    temp_dirs = hfss_dirs.use_d_drive()      # HFSS solver scratch on D: (C: fills up -> COM error at Analyze)
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
    sim.run(name=tag, components=comps, open_terminations=opens, port_list=[('feed_in', 'start', 50), (f'fl_{nl}', 'end', 50)])
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
    try: sim.close()
    except Exception: pass


def analyze_chip(tag):
    import theory5 as T, modes5 as M, chip6_theory as CT, chip6_layout as C6, ftf_layout as f3, copy
    S = pd.read_csv(os.path.join(RES, f'{tag}_S.csv'), index_col=0).map(complex)
    f = S.index.values.astype(float); s11, s21 = S.S11.values, S.S21.values
    meta_p = os.path.join(RES, f'{tag}_meta.json'); meta = json.load(open(meta_p)) if os.path.exists(meta_p) else {}
    # chip-level circuit model with the drawn geometry
    d = C6.build_chip(); P0 = T.params_from_layout()
    L = {c['ch']: c for c in json.load(open(C6.LEN_JSON))['channels']}
    pl = lambda n: f3.path_length(d, n)
    P0['l_feed_in'] = pl('feed_in') + 0.04
    nl = len(C6.ORDER)
    segs = [pl('fl_0') + 0.04 + 0.15] + [pl(f'fl_{k}') + 0.3 for k in range(1, nl)] + [pl(f'fl_{nl}') + 0.15]
    chans = []
    for ch, _, _ in C6.ORDER:
        P = copy.deepcopy(P0); P['L_filt'], P['L_res'] = L[ch]['L_filter'], L[ch]['L_res']; chans.append((ch, P))
    t21 = CT.chip_s21(f, P0, [P for _, P in chans], segs)
    db = lambda s: 20 * np.log10(np.abs(s))
    def dips(s, lo, hi):
        m = (f >= lo) & (f <= hi); a = db(s[m]); i, _ = find_peaks(-a, prominence=3, distance=max(1, int(0.002 / (f[1] - f[0]))))
        i = sorted(sorted(i, key=lambda k: a[k])[:2], key=lambda k: f[m][k])
        return [round(float(f[m][k]), 4) for k in i]
    rows = []
    for ch, P in sorted(chans):
        fr = L[ch]['f_r']; lo, hi = fr - 0.09, fr + 0.10
        rm = M.extract(lambda z: T.branch_Y(z, P), P['Z0'], fr + 0.0095, span=0.5)
        m = (f >= lo) & (f <= hi)
        try:
            dfit, _ = M.direct_fit(f[m], s21[m], window=0.09)
        except Exception as e:
            dfit = None
        rows.append(dict(ch=ch, f_design=fr, dips_hfss=dips(s21, lo, hi), dips_model=dips(t21, lo, hi), model=rm, direct=dfit))
    out = dict(tag=tag, date=time.strftime('%Y-%m-%d %H:%M'), note=meta.get('note', ''), hfss=meta.get('hfss', {}),
               cap_source=P0['cap_source'], caps=dict(cin=P0['cin'], cp=P0['cp'], cj=P0['cj']), channels=rows,
               rms_S21_model=float(np.sqrt(np.mean(np.abs(t21 - s21) ** 2))))
    # plots: full span + one zoom per channel
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(14, 12)); gs = fig.add_gridspec(3, 3)
    a0 = fig.add_subplot(gs[0, :])
    a0.plot(f, db(s21), 'k', lw=1.3, label='HFSS S21'); a0.plot(f, db(t21).clip(-90), '--', lw=1, label='chip model (Q3D C, no fit)')
    a0.plot(f, db(s11), color='0.6', lw=.8, label='HFSS S11')
    a0.set_title(f'{tag}  {out["note"]}'); a0.set_ylabel('dB'); a0.legend(fontsize=8); a0.grid(alpha=.3)
    for k, r in enumerate(rows):
        ax = fig.add_subplot(gs[1 + k // 3, k % 3]); fr = r['f_design']; m = (f >= fr - 0.09) & (f <= fr + 0.10)
        ax.plot(f[m], db(s21[m]), 'k', lw=1.3, label='HFSS'); ax.plot(f[m], db(t21[m]).clip(-90), '--', label='model')
        if r['direct']:
            xd = [r['direct'][n] for n in M.NAMES]
            ax.plot(f[m], db(np.conj(M.io_model(f[m], xd))).clip(-90), '-.', label='eq. 2.137 fit')
        ax.set_title(f'ch{r["ch"]}  design f_r {fr:.2f} GHz', fontsize=9); ax.grid(alpha=.3); ax.legend(fontsize=7)
    fig.tight_layout()
    stamp = time.strftime('%Y%m%d_%H%M'); run_dir = os.path.join(RES, f'{stamp}_{tag}'); os.makedirs(run_dir, exist_ok=True)
    for p in glob.glob(os.path.join(RES, f'{tag}_*')):
        try:
            shutil.move(p, os.path.join(run_dir, os.path.basename(p)))
        except (PermissionError, OSError):           # e.g. the run log still open by the shell redirection
            shutil.copy(p, os.path.join(run_dir, os.path.basename(p)))
    for p in [os.path.join(RES, 'q3d5_full_caps.json'), C6.LEN_JSON]:
        if os.path.exists(p): shutil.copy(p, os.path.join(run_dir, f'{tag}_' + os.path.basename(p)))
    png = os.path.join(run_dir, f'{tag}_chip_analysis.png'); fig.savefig(png, dpi=110); plt.close(fig)
    json.dump(out, open(os.path.join(run_dir, f'{tag}_chip_analysis.json'), 'w'), indent=1, default=float)
    try:
        rel = os.path.relpath(png, os.path.dirname(RECORD))
    except ValueError:                                   # record on another drive (tests)
        rel = png
    record(out, rel.replace('\\', '/'))
    return out, run_dir


def record(o, png_rel):
    f4 = lambda v: '—' if v is None else f'{v:.4f}'
    f1 = lambda v: '—' if v is None else f'{v:.1f}'
    h = o.get('hfss', {})
    Lr = [f'\n## {o["tag"]} — full chip, 6 channels  ({o["date"]})\n', f'{o["note"]}\n' if o['note'] else '',
          f'HFSS: {", ".join(f"{k} {v}" for k, v in h.items())}\n' if h else '',
          f'![{o["tag"]}]({png_rel.replace(" ", "%20")})\n',
          f'Capacitances (fF, `{o["cap_source"]}`): C_in {o["caps"]["cin"][0]:.2f}, C_p {o["caps"]["cp"][0]:.2f}, C_j {o["caps"]["cj"][0]:.2f}; '
          f'rms |ΔS21| HFSS vs chip model {o["rms_S21_model"]:.4f}\n',
          '| ch | design f_r | dips HFSS (GHz) | dips model (GHz) | model ω_p / ω_r (GHz) | model J / κ_p matched (MHz) | fit ω_p / ω_r (GHz) | fit J / κ_p (MHz) |',
          '|---|---|---|---|---|---|---|---|']
    for r in o['channels']:
        m, dfit = r['model'] or {}, r['direct'] or {}
        Lr.append(f'| {r["ch"]} | {r["f_design"]:.2f} | {", ".join(map(str, r["dips_hfss"]))} | {", ".join(map(str, r["dips_model"]))} | '
                  f'{f4(m.get("wp"))} / {f4(m.get("wr"))} | {f1(m.get("J"))} / {f1(m.get("kp"))} | '
                  f'{f4(dfit.get("wp"))} / {f4(dfit.get("wr"))} | {f1(dfit.get("J"))} / {f1(dfit.get("kp"))} |')
    Lr.append('\nThe fit κ_p is the in-chip value (C_in, standing waves); the model κ_p is the matched-line value.\n')
    open(RECORD, 'a', encoding='utf-8').write('\n'.join(Lr) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == '--analyze':
        o, rd = analyze_chip(sys.argv[2]); print(json.dumps(o['channels'], indent=1, default=float)[:3000]); print('run folder', rd)
        sys.exit()
    tag = sys.argv[1]
    f0, f1 = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (6.1, 7.5)
    count = int(sys.argv[4]) if len(sys.argv) > 4 else int(round((f1 - f0) * 1000 / 0.25)) + 1
    stype = sys.argv[5] if len(sys.argv) > 5 else 'Fast'
    run_hfss(tag, f0, f1, count, stype)
    o, rd = analyze_chip(tag)
    for r in o['channels']:
        print(f'ch{r["ch"]} design {r["f_design"]:.2f}: HFSS dips {r["dips_hfss"]}, model dips {r["dips_model"]}')
    print('run folder', rd)
