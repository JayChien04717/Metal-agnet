"""Design6: choose C_p, C_j and the two lengths so the S21 looks like Cai 2024 Fig. 3(f)
(omega_r = omega_p, two dips split by ~2J with a peak back to the baseline), using the circuit model fitted to HFSS.

Loaded modes = complex zeros of the total admittance at node n0 (branch + both feed sides, ports = 50 Ohm).
  f = f' + i f''  ->  linewidth kappa/2pi = 2 |f''|.
  kappa_p: filter mode with C_j -> 0.   Alignment: the two hybrid modes have equal linewidths.
  J: at alignment f_pm = f0 - i kappa/4 +- sqrt(J^2 - kappa^2/16)  ->  J = sqrt((split/2)^2 + (kappa_p/4)^2).
Usage: python design_target.py [fit json]   (default results/fit_c6_T3h_eig.json)
"""
import os, sys, json, copy, warnings
warnings.filterwarnings('ignore')
import numpy as np
from scipy.optimize import least_squares
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import fit_theory as F
from theory_abcd import branch_Y, s21_theory

TARGET = dict(f0=6.70, kappa=45.0, J=8.8)          # GHz, MHz, MHz (paper: kappa_p 45, J 8.8; our channel at 6.7 GHz)


def Y_n0(f, P):
    """Total admittance at n0 for complex f (GHz)."""
    Z0 = P['Z0']; Y0 = 1 / Z0; w = 2 * np.pi * f * 1e9; b = w / P['v']
    cs, ca, cb = P['cin']
    Ya = Y0 + 1j * w * ca * 1e-15                            # P1 side: matched line -> Y0, plus electrode a
    Ys = 1j * w * cs * 1e-15 * Ya / (1j * w * cs * 1e-15 + Ya) + 1j * w * cb * 1e-15
    t = np.tan(b * P['l_feed_c'] * 1e-3)
    Yl = Y0 * (Ys + 1j * Y0 * t) / (Y0 + 1j * Ys * t)
    return branch_Y(f, P) + Yl + Y0                         # P2 side: matched line


def root(P, g):
    """Complex zero of Y_n0 near guess g (secant)."""
    a, b = complex(g), complex(g) + 1e-4 + 1e-4j
    fa, fb = Y_n0(a, P), Y_n0(b, P)
    for _ in range(80):
        if fb == fa: break
        c = b - fb * (b - a) / (fb - fa)
        a, fa, b, fb = b, fb, c, Y_n0(c, P)
        if abs(b - a) < 1e-10: break
    return b


def modes(P, lo=6.3, hi=7.1):
    """All loaded modes in [lo, hi]: start the secant from a grid, keep distinct converged roots."""
    out = []
    for g in np.arange(lo, hi, 0.01):
        for gi in (0.002, 0.02):
            r = root(P, g + 1j * gi)
            if lo < r.real < hi and abs(Y_n0(r, P)) < 1e-9 and all(abs(r - o) > 1e-5 for o in out):
                out.append(r)
    return sorted(out, key=lambda z: z.real)


def metrics(P, guess=None):
    """f_p, kappa_p (filter alone), the two hybrid modes, J estimate."""
    Pf = copy.deepcopy(P); Pf['cj'] = (1e-9, P['cj'][1], P['cj'][2])
    mf = [m for m in modes(Pf) if abs(m.imag) > 1e-4]        # the filter is the only lossy mode
    fp = min(mf, key=lambda z: abs(z.real - (guess or z.real)))
    kp = 2e3 * abs(fp.imag)
    hm = sorted([m for m in modes(P, fp.real - 0.15, fp.real + 0.15) if abs(m.imag) > 1e-9], key=lambda z: abs(z.real - fp.real))[:2]
    hm = sorted(hm, key=lambda z: z.real)
    split = 1e3 * abs(hm[1].real - hm[0].real) if len(hm) == 2 else np.nan
    J = np.sqrt((split / 2) ** 2 + (kp / 4) ** 2)
    return dict(f_p=fp.real, kappa_p=kp, hybrid=[(m.real, 2e3 * abs(m.imag)) for m in hm], split=split, J=J)


def main(fit_json):
    fit = json.load(open(fit_json)); x0 = np.array([fit['params'][k] for k in F.NAMES])
    P0, f, s11, s21 = F.load_run(fit['tags'][0])
    if os.environ.get('LFC'): P0['l_feed_c'] = float(os.environ['LFC'])      # C_in -> n0 line length (mm), shape knob
    N = F.NAMES.index
    cur = metrics(F.apply(P0, x0))
    print('current (fit):', {k: (np.round(v, 4) if not isinstance(v, list) else [tuple(np.round(t, 4)) for t in v]) for k, v in cur.items()})

    from scipy.optimize import brentq
    x = x0.copy(); f0 = TARGET['f0']

    seed = [complex(cur['f_p'], cur['kappa_p'] / 2e3)]
    fpole = [f0]                        # target of the loaded filter pole; shifted so the two S21 zeros centre on f0

    def filt(xx):                       # filter-alone loaded mode, followed by continuation from the last one
        P = F.apply(P0, xx); P['cj'] = (1e-9, P['cj'][1], P['cj'][2])
        r = root(P, seed[0])
        if not (abs(Y_n0(r, P)) < 1e-9 and abs(r.imag) > 1e-4 and abs(r - seed[0]) < 0.3):
            r = min([m for m in modes(P, 5.8, 7.6) if abs(m.imag) > 1e-4], key=lambda z: abs(z - seed[0]))
        seed[0] = r
        return r

    def secant(fun, a, b, tol, n=30):   # 1-D secant, fun(a), fun(b) evaluated in order
        fa, fb = fun(a), fun(b)
        for _ in range(n):
            if abs(fb) < tol or fb == fa: break
            a, fa, b = b, fb, b - fb * (b - a) / (fb - fa); fb = fun(b)
        return b

    def set_dLf(xx):                    # filter length so that f_p = f0 (df/dL ~ -1.5 GHz/mm)
        def g(d):
            xx[N('dLf')] = d; return filt(xx).real - fpole[0]
        d0 = xx[N('dLf')]; xx[N('dLf')] = secant(g, d0, d0 + (filt(xx).real - fpole[0]) / 1.5, 1e-5)

    def hyb(xx):                        # the two hybrid modes nearest f0 (fine seeds)
        P = F.apply(P0, xx); out = []
        for g in np.arange(f0 - 0.08, f0 + 0.08, 0.002):
            for gi in (1e-4, 0.01):
                r = root(P, g + 1j * gi)
                if abs(Y_n0(r, P)) < 1e-9 and abs(r.real - f0) < 0.1 and all(abs(r - o) > 1e-5 for o in out): out.append(r)
        return sorted(sorted(out, key=lambda z: abs(z.real - f0))[:2], key=lambda z: z.real)

    def zeros(xx, w=0.15):
        """S21 transmission zeros = the branch shorts n0 (poles of branch_Y) = lossless modes of filter + resonator
        with n0 grounded. Their minimum splitting vs resonator length is 2J (independent of kappa_p)."""
        P = F.apply(P0, xx); ff = np.linspace(f0 - w, f0 + w, 30001)
        X = (1 / branch_Y(ff, P)).imag                     # reactance: zero of X <-> pole of Y; Foster: X goes - -> +
        i = np.where((X[:-1] < 0) & (X[1:] > 0))[0]
        z = [ff[k] - X[k] * (ff[k + 1] - ff[k]) / (X[k + 1] - X[k]) for k in i]
        return sorted(z, key=lambda q: abs(q - f0))[:2]

    def zsplit(xx):
        z = zeros(xx); return abs(z[1] - z[0]) if len(z) == 2 else 1.0

    def align(xx):                      # minimise the splitting over the resonator length
        from scipy.optimize import minimize_scalar
        ds = np.linspace(-0.35, 0.15, 51); v = []
        for d in ds:
            xx[N('dLr')] = d; v.append(zsplit(xx))
        k = int(np.argmin(v)); lo, hi = ds[max(k - 1, 0)], ds[min(k + 1, 50)]
        r = minimize_scalar(lambda d: (xx.__setitem__(N('dLr'), d), zsplit(xx))[1], bounds=(lo, hi), method='bounded',
                            options=dict(xatol=1e-6))
        xx[N('dLr')] = r.x

    for it in range(5):                 # C_p/C_j barely interact; the pole target is re-centred each sweep
        if it:
            z = zeros(x); fpole[0] += f0 - (z[0] + z[1]) / 2
        def gk(s):                      # kappa_p(s_p) - target, with f_p re-centred each time
            x[N('s_p')] = s; set_dLf(x); return 2e3 * abs(filt(x).imag) - TARGET['kappa']
        s0 = x[N('s_p')]; x[N('s_p')] = secant(gk, s0, s0 * 1.05, 0.05); set_dLf(x)
        print(f'  C_p step: s_p {x[N("s_p")]:.4f}, dLf {x[N("dLf")]:.4f}, f_p {filt(x).real:.5f}, kappa_p {2e3 * abs(filt(x).imag):.2f}')

        align(x)                        # dLr at the anticrossing of the two S21 zeros (omega_r = omega_p)

        def gj(s):                      # J = half the minimum splitting of the two S21 zeros
            x[N('s_j')] = s; align(x); return 1e3 * zsplit(x) / 2 - TARGET['J']
        s0 = x[N('s_j')]; x[N('s_j')] = secant(gj, s0, s0 * TARGET['J'] / max(1e3 * zsplit(x) / 2, 1e-3), 0.02)
        print(f'iter {it}:', dict(zip(F.NAMES, np.round(x, 4))), f'zeros {np.round(zeros(x), 5)}')

    class S: pass
    sol = S(); sol.x = x[[N('s_p'), N('s_j'), N('dLf'), N('dLr')]]
    P = F.apply(P0, x); fp = filt(x); h = hyb(x); kp = 2e3 * abs(fp.imag); z = zeros(x)
    m = dict(f_p=fp.real, kappa_p=kp, hybrid=[(q.real, 2e3 * abs(q.imag)) for q in h], s21_zeros=z,
             split=1e3 * zsplit(x), J=1e3 * zsplit(x) / 2)
    print('solution:', dict(zip(['s_p', 's_j', 'dLf', 'dLr'], np.round(sol.x, 4))))
    print('  ->', {k: np.round(np.array(v, dtype=float), 4).tolist() for k, v in m.items()})
    # geometry targets: same HFSS/Q3D ratio as the fitted run
    sp0, sj0 = x0[N('s_p')], x0[N('s_j')]
    geo = dict(C_p_Q3D_fF=P0['cp'][0] * sol.x[0] / sp0, C_j_Q3D_fF=P0['cj'][0] * sol.x[1] / sj0,
               dL_filter_mm=sol.x[2] - x0[N('dLf')], dL_res_mm=sol.x[3] - x0[N('dLr')])
    print('geometry change (vs c6_T3):', {k: round(float(v), 4) for k, v in geo.items()})
    out = os.path.join(HERE, 'results', 'design_target' + (f"_lfc{os.environ['LFC']}" if os.environ.get('LFC') else ''))
    json.dump(dict(target=TARGET, fit=fit_json, y=dict(zip(['s_p', 's_j', 'dLf', 'dLr'], map(float, sol.x))),
                   metrics=dict(f_p=m['f_p'], kappa_p=m['kappa_p'], J=m['J'], split=m['split'], hybrid=m['hybrid'], s21_zeros=m['s21_zeros'], f_pole_target=fpole[0]),
                   geometry=geo), open(out + '.json', 'w'), indent=1)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    ff = np.linspace(TARGET['f0'] - 0.075, TARGET['f0'] + 0.075, 6001)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(ff, np.abs(s21_theory(ff, P)[1]), 'k', label='model, design target')
    ax.set_xlabel('f (GHz)'); ax.set_ylabel('|S21|'); ax.grid(alpha=.3); ax.legend()
    ax.set_title(f"kappa_p {m['kappa_p']:.1f} MHz, J {m['J']:.1f} MHz, split {m['split']:.1f} MHz")
    fig.tight_layout(); fig.savefig(out + '.png', dpi=120); print('saved', out + '.png')
    return P, sol.x, geo


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'results', 'fit_c6_T3h_eig.json'))
