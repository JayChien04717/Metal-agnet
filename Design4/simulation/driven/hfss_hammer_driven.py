"""Driven S21 of a hammer design and the capacitor-model fit to it (no eigenmode involved).

Geometry: same build as hfss_hammer_straps.py (purcell_6500, straps on every CPW, 50 Ohm lumped ports feed_in / feed_out).
Driven setup: adaptive on S at --adapt-f (default 6.49 GHz, next to the resonator), length mesh seeds --mesh-line /
--mesh-cap (default 10 / 5 um, the project rule), first order, max --maxp (8) passes, Delta S --ds; Fast sweep
--f0..--f1 at --step MHz with fields saved (any frequency can be plotted in AEDT afterwards).
Fit (theory5 full channel: C_in, n0, stub - C_p - filter - C_j - resonator, + complex scale A e^{i phi} and delay tau):
  Q3D theory  C_p, C_j fixed at Q3D, only the two line lengths and scale/delay fitted
  S21 fit     C_p, C_j free as well
The circuit kappa_p and J of both come from the poles of the fitted circuit (2-mode inversion).
E field: |E| 1 um above the metal at the S21 minimum and at the resonator feature (E_<f>GHz.fld / .png).

Usage: python driven/hfss_hammer_driven.py --spacing 0.4 --gap 70 --trim=-11
       python driven/hfss_hammer_driven.py ... --fit-only           (refit the saved S.csv)
Output: results/<date>/hammer<spacing>_gap<gap>_s21_<trim>/ (S.csv, convergence.csv, meta.json, fit.json, s21_fit.png, E_*)
"""
import argparse, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np

SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DESIGN = os.path.dirname(SIM); HFSS = os.path.dirname(DESIGN)
sys.path[:0] = [os.path.join(SIM, 'eigenmode'), os.path.join(DESIGN, 'pattern design'), os.path.join(HFSS, 'Design5', 'simulation')]
import hfss_hammer_straps as E     # geometry build / straps / mesh ops only
import purcell_6500 as P


def tag(t):
    return f'm{-t:g}um' if t < 0 else f'p{t:g}um'


def run_driven(args, out):
    from qiskit_metal import Dict
    from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
    import save_aedt, hfss_dirs
    hfss_dirs.use_d_drive(); save_aedt.use_project(save_aedt.project_file(SIM))
    E.MESH_LINE, E.MESH_CAP = args.mesh_line, args.mesh_cap
    d, comps = E.build(args.trim, autobonds=False)

    class Drv(ScatteringImpedanceSim):
        def _analyze(self):
            names = P.draw_hammer_straps(self.renderer, d, which=args.straps)
            print(f'ground straps ({args.straps}): {len(names)} drawn', flush=True)
            E.mesh_ops(self.renderer, d, comps)
            super()._analyze()

    sim = Drv(d, 'hfss')
    for k, v in dict(freq_ghz=args.adapt_f, max_delta_s=args.ds, max_passes=args.maxp, min_passes=3, min_converged=2,
                     pct_refinement=20, basis_order=1, vars=Dict()).items():
        sim.setup[k] = v
    count = int(round((args.f1 - args.f0) * 1e3 / args.step)) + 1
    sim.setup.sweep_setup.update(name='Sweep', start_ghz=args.f0, stop_ghz=args.f1, count=count, type='Fast',
                                 save_fields=True)     # Jay 2026-11-12: keep the E field
    t0 = time.time()
    name = save_aedt.design_name(f'{E.prefix()}_s21_{tag(args.trim)}')
    sim.run(name=name, components=comps, open_terminations=[], port_list=E.PORTS)
    try:
        conv = sim.renderer.pinfo.setup.get_convergence()[0]; conv.to_csv(os.path.join(out, 'convergence.csv'))
        print(conv.to_string(), flush=True)
    except Exception as e:
        print('convergence not read:', repr(e)[:150])
    _, _, S = sim.renderer.get_params(['S11', 'S21'])
    S.index.name = 'f_GHz'; S.to_csv(os.path.join(out, 'S.csv'))
    json.dump(dict(design=name + '_hfss', trim_um=args.trim, spacing_mm=args.spacing, gap_um=args.gap, straps=args.straps,
                   mesh_seeds=[args.mesh_line, args.mesh_cap], adapt_f_GHz=args.adapt_f, max_delta_s=args.ds,
                   max_passes=args.maxp, basis_order=1, sweep=f'{args.f0}-{args.f1} GHz, {count} pts, Fast, fields saved',
                   minutes=round((time.time() - t0) / 60, 1)), open(os.path.join(out, 'meta.json'), 'w'), indent=1)
    print(f'driven done, {round((time.time() - t0) / 60, 1)} min', flush=True)
    try:
        f, s = load_S(out)
        export_fields(sim, d, out, args, key_freqs(f, s))
    except Exception as e:
        print('!!! E-field export failed (fields stay in the AEDT project):', repr(e)[:200], flush=True)
    try:
        sim.renderer.pinfo.project._project.Save()
    except Exception:
        pass


def load_S(out):
    import pandas as pd
    S = pd.read_csv(os.path.join(out, 'S.csv'), index_col=0)
    f = S.index.values.astype(float)
    col = [c for c in S.columns if 's21' in c.lower()][0]
    return f, np.array([complex(str(x).replace(' ', '').replace('i', 'j')) for x in S[col]])


def key_freqs(f, s):
    """(S21 minimum, resonator feature). Feature: largest deviation of complex S21 from its 4 MHz running median
    (the broad filter notch is smooth in linear S21); accepted only if > 10x the median deviation."""
    from scipy.ndimage import median_filter
    fmin = float(f[np.argmin(np.abs(s))])
    w = max(5, int(round(0.004 / np.median(np.diff(f)))) | 1)
    sm = median_filter(s.real, w, mode='nearest') + 1j * median_filter(s.imag, w, mode='nearest')
    dev = np.abs(s - sm)
    k = int(np.argmax(dev))
    return fmin, (float(f[k]) if dev[k] > 10 * np.median(dev) and dev[k] > 1e-3 else None)


def export_fields(sim, d, out, args, freqs):
    """|E| (complex magnitude, port 1 driven, 1 W) on a grid just above the chip at each frequency."""
    fr = sim.renderer.pinfo.design._design.GetModule('FieldsReporter')
    mm = lambda v: float(d.parse_value(v))
    c = d.chips.main.size
    x0, x1 = mm(c.center_x) - mm(c.size_x) / 2, mm(c.center_x) + mm(c.size_x) / 2
    y0, y1 = mm(c.center_y) - mm(c.size_y) / 2, mm(c.center_y) + mm(c.size_y) / 2
    st, z = args.field_step / 1e3, args.field_z / 1e3
    soln = f'{sim.sim_setup_name} : {sim.sweep_name}'
    for f in [x for x in freqs if x]:
        path = os.path.join(out, f'E_{f:.4f}GHz.fld').replace(os.sep, '/')
        fr.CalcStack('clear'); fr.EnterQty('E'); fr.CalcOp('CmplxMag'); fr.CalcOp('Mag')
        fr.ExportOnGrid(path, [f'{x0}mm', f'{y0}mm', f'{z}mm'], [f'{x1}mm', f'{y1}mm', f'{z}mm'],
                        [f'{st}mm', f'{st}mm', f'{st}mm'], soln, ['Freq:=', f'{f}GHz', 'Phase:=', '0deg'],
                        True, 'Cartesian', ['0mm', '0mm', '0mm'], False)
        plot_field(path, f)
        print(f'E field at {f:.5f} GHz: {path}', flush=True)


def plot_field(path, f):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    rows = []
    for line in open(path):
        p = line.split()
        if len(p) == 4:
            try:
                rows.append([float(v) for v in p])
            except ValueError:
                pass
    a = np.array(rows); a = a[np.isfinite(a[:, 3])]
    xs, ys = np.unique(a[:, 0]), np.unique(a[:, 1])
    img = np.full((len(ys), len(xs)), np.nan)
    img[np.searchsorted(ys, a[:, 1]), np.searchsorted(xs, a[:, 0])] = a[:, 3]
    k = 1e3 if xs.max() < 0.1 else 1
    v = img[np.isfinite(img) & (img > 0)]
    fig, ax = plt.subplots(figsize=(7, 6))
    m = ax.imshow(img, origin='lower', extent=[xs[0] * k, xs[-1] * k, ys[0] * k, ys[-1] * k], cmap='inferno',
                  norm=LogNorm(vmin=np.percentile(v, 5), vmax=v.max()))
    fig.colorbar(m, ax=ax, label='|E| (V/m), port 1 at 1 W')
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)'); ax.set_title(f'|E| just above the chip, {f:.4f} GHz', fontsize=10)
    fig.tight_layout(); fig.savefig(path[:-4] + '.png', dpi=130); plt.close(fig)


def fit(args, out):
    from scipy.optimize import least_squares
    import theory5, s21_model_vs_eigen as M
    from hfss_hammer_straps import point_J
    f, s = load_S(out)
    base, _ = P.tune_lengths()
    cpq, cjq = M.CP_Q3D, M.cj_q3d(args.gap)
    fmin, fres = key_freqs(f, s)
    fr0 = fres if fres else fmin - 0.013
    print(f'S21 minimum {fmin:.5f} GHz, resonator feature {fres}', flush=True)
    # starting lengths: Q3D circuit with bare f_p near the S21 minimum, f_r at the feature
    L0 = M.tune_bare(base, cpq, cjq, fmin + 0.008, fr0, (base['L_filt'], base['L_res']))
    model = lambda x, ff: theory5.s21_theory(ff, M.circuit(base, x[2], x[3], x[0], x[1]))[1]
    full = lambda x, ff: x[4] * np.exp(1j * (x[5] - 2 * np.pi * ff * x[6])) * model(x, ff)
    near = np.abs(f - fmin) < 0.04
    wgt = np.where(near, 3.0, 1.0)
    if fres:
        wgt = np.where(np.abs(f - fres) < 0.003, 10.0, wgt)     # the narrow resonator feature has few points
    off = ~near
    m0 = model([*L0, cpq, cjq], f)
    ph = np.unwrap(np.angle(s[off] / m0[off])); tau0 = -np.polyfit(2 * np.pi * f[off], ph, 1)[0]
    phi0 = float(np.angle(np.mean(s[off] / m0[off] * np.exp(2j * np.pi * f[off] * tau0))))
    A0 = float(np.median(np.abs(s[off] / m0[off])))
    xs = [1e-3, 1e-3, 1, 0.05, 0.01, 0.1, 1e-3]

    def solve(x0, free_c):
        idx = [0, 1, 4, 5, 6] + ([2, 3] if free_c else [])
        lo_all = [x0[0] - 0.08, x0[1] - 0.08, 20, 0.02, 0.2 * x0[4], x0[5] - np.pi, x0[6] - 0.5]
        hi_all = [x0[0] + 0.08, x0[1] + 0.08, 80, 4.0, 5 * x0[4], x0[5] + np.pi, x0[6] + 0.5]

        def r(y):
            x = list(x0)
            for i, v in zip(idx, y):
                x[i] = v
            e = wgt * (full(x, f) - s)
            return np.concatenate([e.real, e.imag])
        sl = least_squares(r, [x0[i] for i in idx], bounds=([lo_all[i] for i in idx], [hi_all[i] for i in idx]),
                           x_scale=[xs[i] for i in idx], diff_step=1e-6)
        x = list(x0)
        for i, v in zip(idx, sl.x):
            x[i] = v
        return x, sl.cost

    xq, _ = solve([*L0, cpq, cjq, A0, phi0, tau0], free_c=False)          # Q3D C, lengths + scale fitted
    best = None
    for cpi, cji in ((cpq, cjq), (cpq, 2 * cjq), (cpq * 1.15, cjq), (cpq * 0.9, 0.5 * cjq)):
        x, cost = solve([xq[0], xq[1], cpi, cji, *xq[4:]], free_c=True)
        if best is None or cost < best[1]:
            best = (x, cost)
    xf = best[0]

    def circ(x):
        p = M.circuit(base, x[2], x[3], x[0], x[1])
        pz = M.poles(p, [(fmin + 0.008, 40.0), (fr0, 1.0)])
        inv = point_J([z[0] for z in pz], [z[1] / 2e3 for z in pz])
        return p, pz, inv
    pq, pzq, iq = circ(xq); pf, pzf, iff = circ(xf)
    rms = lambda x: float(np.sqrt(np.mean(np.abs(full(x, f) - s) ** 2)))
    res = dict(s21_min_GHz=fmin, resonator_feature_GHz=fres,
               q3d=dict(C_fF=dict(cp=cpq, cj=cjq), L_mm=xq[:2], poles=pzq, kappa_p_MHz=iq['kappa_p_MHz'], J_MHz=iq['J_MHz'],
                        f_p_GHz=iq['f_p_GHz'], f_r_GHz=iq['f_r_GHz'], rms=rms(xq)),
               fit=dict(C_fF=dict(cp=xf[2], cj=xf[3]), L_mm=xf[:2], poles=pzf, kappa_p_MHz=iff['kappa_p_MHz'], J_MHz=iff['J_MHz'],
                        f_p_GHz=iff['f_p_GHz'], f_r_GHz=iff['f_r_GHz'], rms=rms(xf)),
               scale=dict(A=xf[4], phi=xf[5], tau_ns=xf[6]))
    json.dump(res, open(os.path.join(out, 'fit.json'), 'w'), indent=1, default=float)
    for k in ('q3d', 'fit'):
        r_ = res[k]
        print(f"{k:4s}: C_p {r_['C_fF']['cp']:.2f} fF, C_j {r_['C_fF']['cj']:.3f} fF -> kappa_p {r_['kappa_p_MHz']:.1f} MHz, "
              f"J {r_['J_MHz']:.2f} MHz, f_p {r_['f_p_GHz']:.5f}, f_r {r_['f_r_GHz']:.5f} GHz, rms {r_['rms']:.2e}", flush=True)

    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    db = lambda z: 20 * np.log10(np.abs(z))
    spans = [None, 40, 4] if fres else [None, 40]
    fig, axs = plt.subplots(1, len(spans), figsize=(6 * len(spans), 4.8))
    for ax, span in zip(axs, spans):
        ax.plot(f, db(s), color='k', lw=2.4, alpha=.4, label='HFSS driven S21')
        ax.plot(f, db(full(xf, f)), 'C3', lw=1.2,
                label=f"capacitor model fitted to S21: C_p {xf[2]:.2f}, C_j {xf[3]:.3f} fF\n"
                      f"  κ_p {iff['kappa_p_MHz']:.1f}, J {iff['J_MHz']:.2f} MHz")
        ax.plot(f, db(full(xq, f)), 'C0--', lw=1.1,
                label=f"Q3D capacitors: C_p {cpq:.2f}, C_j {cjq:.3f} fF (lengths fitted)\n"
                      f"  κ_p {iq['kappa_p_MHz']:.1f}, J {iq['J_MHz']:.2f} MHz")
        if span:
            c0 = fmin if span > 10 else fres
            ax.set_xlim(c0 - span / 1e3, c0 + span / 1e3)
        ax.set_xlabel('f (GHz)'); ax.set_ylabel('|S21| (dB)'); ax.grid(alpha=.3)
    axs[0].legend(fontsize=7, loc='lower left')
    fig.suptitle(f'{args.spacing * 1e3:.0f} µm spacing, hammer gap {args.gap:g} µm, trim {args.trim:+g} µm, bonds on every CPW: '
                 f'HFSS driven S21 vs capacitor model', fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(out, 's21_fit.png'), dpi=140)
    print('figure:', os.path.join(out, 's21_fit.png'))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--spacing', type=float, default=0.4)
    ap.add_argument('--gap', type=float, default=70)
    ap.add_argument('--trim', type=float, default=-11)
    ap.add_argument('--straps', default='cpw')
    ap.add_argument('--f0', type=float, default=6.40)
    ap.add_argument('--f1', type=float, default=6.60)
    ap.add_argument('--step', type=float, default=0.05, help='sweep step (MHz)')
    ap.add_argument('--adapt-f', type=float, default=6.49, help='adaptive solve frequency (GHz)')
    ap.add_argument('--ds', type=float, default=0.01, help='max Delta S')
    ap.add_argument('--maxp', type=int, default=8)
    ap.add_argument('--mesh-line', default='10um')
    ap.add_argument('--mesh-cap', default='5um')
    ap.add_argument('--field-step', type=float, default=10, help='E-field grid step (um)')
    ap.add_argument('--field-z', type=float, default=1, help='E-field grid height above the metal (um)')
    ap.add_argument('--out')
    ap.add_argument('--fit-only', action='store_true')
    args = ap.parse_args()
    E.SPACING_MM, E.GAP_UM, E.NO_HAMMER = args.spacing, args.gap, False
    out = args.out or os.path.join(SIM, 'results', time.strftime('%Y-%m-%d'),
                                   f'hammer{args.spacing * 1e3:.0f}_gap{args.gap:g}_s21_{tag(args.trim)}')
    os.makedirs(out, exist_ok=True)
    if not args.fit_only:
        run_driven(args, out)
    fit(args, out)


if __name__ == '__main__':
    main()
