"""S21 of a hammer eigenmode scan: circuit model with Q3D C vs the circuit fitted to the HFSS eigenmodes.

Circuit: theory5 ABCD channel (port - C_in - n0 - port, branch stub - C_p - filter - C_j - resonator), as the layout.
  Q3D-C circuit   C_p, C_j from Q3D; L_filt, L_res tuned so the bare f_p, f_r equal the HFSS per-point values.
  fitted circuit  per simulated trim, (L_filt, L_res, C_p, C_j) solved so the two complex circuit poles equal the two
                  HFSS eigenmodes (f and kappa of each): 4 unknowns, 4 numbers, no mode labelling.
Both circuits and HFSS are reduced to (f_p, kappa_p, f_r, J) with the same per-point 2-mode inversion (point_J).
S21 is the full channel (C_in in series, -16 dB baseline); the fitted curve is the HFSS-equivalent S21.

Usage: python eigenmode/s21_model_vs_eigen.py [--res results/2026-11-12/hammer400_gap70_cpwwb] [--gap 70]
Output in --res: s21_model_vs_eigen.png, s21_model_vs_eigen.json
"""
import argparse, json, os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np
from scipy.optimize import least_squares, root

SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DESIGN = os.path.dirname(SIM)
sys.path[:0] = [os.path.join(DESIGN, 'pattern design'), os.path.dirname(os.path.abspath(__file__))]
import purcell_6500 as P
import theory5, modes5
from hfss_hammer_straps import point_J

CP_Q3D = 39.76                      # 400 um spacing, Q3D (q3d400_gap_sweep), same at every gap
CJ_Q3D = {9: 2.202, 20: 1.291, 30: 0.916, 45: 0.619, 60: 0.464, 80: 0.356}


def cj_q3d(gap):
    g = np.array(sorted(CJ_Q3D), float); c = np.array([CJ_Q3D[k] for k in sorted(CJ_Q3D)])
    return float(np.exp(np.interp(np.log(gap), np.log(g), np.log(c))))


def circuit(base, cp, cj, Lf, Lr):
    p = dict(base); p['cp'] = (cp,) + tuple(base['cp'][1:]); p['cj'] = (cj,) + tuple(base['cj'][1:])
    p['L_filt'], p['L_res'] = Lf, Lr
    return p


def ext(p, fc):
    """modes5.extract in the real channel environment (C_in on the input side, matched output)."""
    return modes5.extract(lambda f: theory5.branch_Y(f, p), p['Z0'], fc, span=0.55,
                          Yenv=lambda f: theory5.n0_env_Y(f, p))


def poles(p, guesses):
    """The two complex channel poles (f GHz, kappa MHz), sorted by f: complex root search from the guesses
    (f GHz, kappa MHz) in the C_in environment; dip tracking (modes5.extract) if the two searches meet."""
    D = lambda z: p['Z0'] * (theory5.n0_env_Y(z / 1e3, p) + theory5.branch_Y(z / 1e3, p))
    z = [modes5.croot(D, f * 1e3 + 0.5j * k, h=0.005) for f, k in guesses]
    if abs(z[0] - z[1]) > 0.05:
        return sorted((c.real / 1e3, 2 * abs(c.imag)) for c in z)
    m = ext(p, float(np.mean([g[0] for g in guesses])))
    return sorted(m['modes']) if m is not None else sorted((c.real / 1e3, 2 * abs(c.imag)) for c in z)


def tune_bare(base, cp, cj, fp, fr, L0):
    """L_filt, L_res for bare f_p, f_r (zero of the filter / resonator reactance, modes5.extract definition)."""
    def r(L):
        m = ext(circuit(base, cp, cj, *L), 0.5 * (fp + fr))
        return np.array([1., 1.]) if m is None else np.array([m['wp'] - fp, m['wr'] - fr])
    s = root(r, L0)
    return tuple(s.x)


def s21(f, p):
    """Full channel with C_in (theory5.s21_theory)."""
    return theory5.s21_theory(f, p)[1]


def s21_cm(f, fp, kp, fr, J):
    dp, dr = (f - fp) * 1e3, (f - fr) * 1e3
    return 1 - (kp / 2) / (kp / 2 - 1j * dp + J ** 2 / (1e-9 - 1j * dr))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--res', default=os.path.join(SIM, 'results', '2026-11-12', 'hammer400_gap70_cpwwb'))
    ap.add_argument('--gap', type=float, default=70)
    ap.add_argument('--span', type=float, default=110, help='plot half-span (MHz)')
    a = ap.parse_args()
    S = json.load(open(os.path.join(a.res, 'summary.json')))
    base, _ = P.tune_lengths()
    cjq = cj_q3d(a.gap)
    rows = []
    L0 = (base['L_filt'], base['L_res'])
    for name in sorted(S['cases'], key=lambda n: S['cases'][n]['trim_um']):
        c = S['cases'][name]; h = S['per_point'][name]
        modes_h = sorted(zip(c['f_GHz'], c['kappa_MHz']))
        # Q3D-C circuit at the HFSS bare frequencies
        Lq = tune_bare(base, CP_Q3D, cjq, h['f_p_GHz'], h['f_r_GHz'], L0)
        q = circuit(base, CP_Q3D, cjq, *Lq)
        pq = poles(q, modes_h)
        iq = point_J([z[0] for z in pq], [z[1] / 2e3 for z in pq])

        # fitted circuit: poles == HFSS modes
        def r(x):
            pf = poles(circuit(base, x[2], x[3], x[0], x[1]), modes_h)
            return np.array([(pf[i][0] - modes_h[i][0]) * 1e3 for i in (0, 1)] +
                            [(pf[i][1] - modes_h[i][1]) for i in (0, 1)])
        starts = [(CP_Q3D * np.sqrt(h['kappa_p_MHz'] / iq['kappa_p_MHz']), cjq * h['J_MHz'] / iq['J_MHz'])]
        starts += [(cp0, cj0) for cp0 in (44, 48, 52) for cj0 in (0.4, 0.5, 0.6)]
        sol = None
        for cp0, cj0 in starts:          # several starts: keep the first exact solution (else the best)
            sl = least_squares(r, [*Lq, cp0, cj0], x_scale=[1e-3, 1e-3, 1, 0.05], diff_step=1e-5)
            if sol is None or np.abs(sl.fun).max() < np.abs(sol.fun).max():
                sol = sl
            if np.abs(sol.fun).max() < 1e-3:
                break
        Lf, Lr, cpf, cjf = sol.x
        fit = circuit(base, cpf, cjf, Lf, Lr)
        pf = poles(fit, modes_h)
        jf = point_J([z[0] for z in pf], [z[1] / 2e3 for z in pf])
        rows.append(dict(case=name, trim_um=c['trim_um'], hfss_modes=modes_h, hfss=h,
                         q3d=dict(cp=CP_Q3D, cj=cjq, L=Lq, poles=pq, inv=iq),
                         fit=dict(cp=cpf, cj=cjf, L=(Lf, Lr), poles=pf, inv=jf, resid=float(np.abs(sol.fun).max()))))
        L0 = Lq
        print(f"trim {c['trim_um']:+4.0f}: HFSS kp {h['kappa_p_MHz']:5.1f} J {h['J_MHz']:4.2f} det {h['detuning_MHz']:+6.1f} | "
              f"Q3D-C kp {iq['kappa_p_MHz']:5.1f} J {iq['J_MHz']:4.2f} | fit C_p {cpf:5.2f} C_j {cjf:5.3f} fF "
              f"(kp {jf['kappa_p_MHz']:5.1f} J {jf['J_MHz']:4.2f}, max resid {np.abs(sol.fun).max():.2g} MHz)", flush=True)

    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    db = lambda s: 20 * np.log10(np.abs(s))
    fig, axs = plt.subplots(1, len(rows), figsize=(5.2 * len(rows), 4.8), sharey=True)
    for ax, rw in zip(np.atleast_1d(axs), rows):
        h = rw['hfss']; fc = 0.5 * (h['f_p_GHz'] + h['f_r_GHz'])
        f = np.linspace(fc - a.span / 1e3, fc + a.span / 1e3, 8001); x = (f - fc) * 1e3
        q = circuit(base, rw['q3d']['cp'], rw['q3d']['cj'], *rw['q3d']['L'])
        fit = circuit(base, rw['fit']['cp'], rw['fit']['cj'], *rw['fit']['L'])
        ax.plot(x, db(s21(f, q)), color='C0', lw=1.5,
                label=f"theory, Q3D C: C_p {CP_Q3D:.2f}, C_j {rw['q3d']['cj']:.3f} fF\n"
                      f"  κ_p {rw['q3d']['inv']['kappa_p_MHz']:.1f}, J {rw['q3d']['inv']['J_MHz']:.2f} MHz")
        ax.plot(x, db(s21(f, fit)), color='C3', lw=1.5,
                label=f"fit to HFSS eigen (poles = HFSS modes): C_p {rw['fit']['cp']:.2f}, C_j {rw['fit']['cj']:.3f} fF\n"
                      f"  κ_p {rw['fit']['inv']['kappa_p_MHz']:.1f}, J {rw['fit']['inv']['J_MHz']:.2f} MHz")
        for fm, _ in rw['hfss_modes']:
            ax.axvline((fm - fc) * 1e3, color='0.6', lw=.7, ls=':')
        ax.set_title(f"trim {rw['trim_um']:+.0f} µm: f_p {h['f_p_GHz']:.4f}, f_r {h['f_r_GHz']:.4f} GHz "
                     f"(Δ {h['detuning_MHz']:+.1f} MHz)", fontsize=8.5)
        ax.set_xlabel(f'f − {fc:.4f} GHz (MHz)'); ax.grid(alpha=.3); ax.legend(fontsize=6.6, loc='lower left')
    np.atleast_1d(axs)[0].set_ylabel('|S21| (dB)'); np.atleast_1d(axs)[0].set_ylim(-55, -12)
    fig.suptitle(f'400 µm spacing, hammer gap {a.gap:g} µm, bonds on every CPW — S21: Q3D-capacitor theory vs '
                 f'capacitor model fitted to HFSS eigenmodes\n'
                 f'full channel with C_in (−16 dB baseline); dotted = HFSS eigenfrequencies '
                 f'(poles) — the S21 dips sit off the poles because of the C_in standing wave', fontsize=9.5)
    fig.tight_layout()
    png = os.path.join(a.res, 's21_model_vs_eigen.png'); fig.savefig(png, dpi=140)
    json.dump(dict(gap_um=a.gap, rows=rows), open(os.path.join(a.res, 's21_model_vs_eigen.json'), 'w'), indent=1)
    print('figure:', png)


if __name__ == '__main__':
    main()
