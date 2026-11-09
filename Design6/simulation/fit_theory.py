"""Fit the Design6 ABCD model (theory_abcd.s21_theory) to HFSS S11/S21 of a finished run.

Free parameters: dLf, dLr (mm, effective-length corrections of filter / resonator: bends, T-junction, cap metal),
s_in, s_p, s_j (scale of the Q3D capacitances, penalised towards 1), kv (scale of the phase velocity).
Residual = complex S11 and S21 in linear scale (deep nulls do not dominate), optional feed-line delay is not fitted
(the model already contains the port-to-port line lengths).

Usage: python fit_theory.py <tag> [<tag> ...]     (tags = run folders in results/, e.g. c6_T3 c6_T3g)
Outputs: results/fit_<tags>.json, results/fit_<tags>.png ; prints fitted params and poles (f, kappa) of the model.
"""
import os, sys, json, glob, copy, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.optimize import least_squares, differential_evolution
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results')
sys.path.insert(0, HERE)
from theory_abcd import s21_theory, branch_Y

EIG = []              # HFSS eigenmode frequencies (GHz) of the cut structure (n0 open), set with --eig f1,f2
EIG_SIGMA = 0.010     # GHz
N0_HALF = 0.15        # mm, each open half of the n0 prime line in the eigenmode model


def model_eig(P, fmin=5.9, fmax=7.3, n=20001):
    """Eigenfrequencies of the model with the feed line cut at n0 (as hfss_eigen6.py): zeros of the total
    susceptance at n0 (Foster: B goes - -> + at a zero, + -> - at a pole)."""
    f = np.linspace(fmin, fmax, n); Z0 = P['Z0']; beta = 2 * np.pi * f * 1e9 / P['v']
    B = (branch_Y(f, P) + 2j / Z0 * np.tan(beta * N0_HALF * 1e-3)).imag
    i = np.where((B[:-1] < 0) & (B[1:] > 0))[0]
    return [float(f[k] - B[k] * (f[k + 1] - f[k]) / (B[k + 1] - B[k])) for k in i]

NAMES = ['dLf', 'dLr', 's_in', 's_p', 's_j', 'kv']
BOUNDS = [(-0.6, 0.6), (-0.6, 0.6), (0.7, 1.3), (0.7, 1.3), (0.5, float(os.environ.get('SJ_MAX', 1.6))), (0.93, 1.07)]
PRIOR_SIGMA = dict(s_in=0.15, s_p=0.15, s_j=float(os.environ.get('SJ_SIGMA', 0.4)), kv=0.03)       # soft prior around the Q3D / nominal values


def load_run(tag):
    d = glob.glob(os.path.join(RES, f'*_{tag}'))
    d = [x for x in d if os.path.isdir(x)][-1]
    P = json.load(open(os.path.join(d, f'{tag}_theory_params.json')))
    for k in ('cin', 'cp', 'cj'):
        P[k] = tuple(P[k])
    S = pd.read_csv(os.path.join(d, f'{tag}_S.csv'), index_col=0).map(lambda z: complex(z))
    return P, S.index.values.astype(float), S.S11.values, S.S21.values


def apply(P0, x):
    p = dict(zip(NAMES, x)); P = copy.deepcopy(P0)
    P['L_filt'] = P0['L_filt'] + p['dLf']
    P['L_res'] = P0['L_res'] + p['dLr']
    P['cin'] = tuple(c * p['s_in'] for c in P0['cin'])
    P['cp'] = tuple(c * p['s_p'] for c in P0['cp'])
    P['cj'] = (P0['cj'][0] * p['s_j'], P0['cj'][1], P0['cj'][2])        # only the series C_j is scaled
    P['v'] = P0['v'] * p['kv']
    return P


def residual(x, datasets):
    r = []
    for P0, f, s11, s21 in datasets:
        t11, t21 = s21_theory(f, apply(P0, x))
        r += [(t21 - s21).real, (t21 - s21).imag, 0.3 * (t11 - s11).real, 0.3 * (t11 - s11).imag]
    r = np.concatenate(r) / np.sqrt(sum(len(d[1]) for d in datasets))
    p = dict(zip(NAMES, x))
    if EIG:
        fm = model_eig(apply(datasets[0][0], x))
        e = [min((abs(m - fe) for m in fm), default=1.0) / EIG_SIGMA for fe in EIG]   # nearest model mode
        r = np.concatenate([r, 0.3 * np.array(e)])
    pr = [(p[k] - 1) / s for k, s in PRIOR_SIGMA.items() if k != 'kv'] + [(p['kv'] - 1) / PRIOR_SIGMA['kv']]
    return np.concatenate([r * 30, 0.02 * np.array(pr)])               # data weight vs soft prior


def model_poles(P, fmin=6.3, fmax=7.1):
    """Resonances of the model seen as dips of |S21| and their widths (full width at -3 dB below the local ridge)."""
    f = np.linspace(fmin, fmax, 80001); _, s = s21_theory(f, P); a = 20 * np.log10(np.abs(s))
    from scipy.signal import find_peaks
    i, pr = find_peaks(-a, prominence=3)
    return [(float(f[k]), float(a[k])) for k in i]


def main(tags):
    datasets = [load_run(t) for t in tags]
    P0 = datasets[0][0]
    lo, hi = zip(*BOUNDS)
    de = differential_evolution(lambda x: float(np.sum(residual(x, datasets) ** 2)), BOUNDS, seed=1, maxiter=150, popsize=20,
                                tol=1e-8, polish=False)
    sol = least_squares(residual, de.x, args=(datasets,), bounds=(lo, hi), x_scale=[0.05, 0.05, .05, .05, .1, .01])
    x = sol.x; P = apply(P0, x)
    print('fit params:', {k: round(float(v), 4) for k, v in zip(NAMES, x)})
    rms = lambda xx: [float(np.sqrt(np.mean(np.abs(s21_theory(f, apply(P0, xx))[1] - s21) ** 2))) for _, f, _, s21 in datasets]
    print('rms |dS21| before/after:', rms(np.array([0, 0, 1, 1, 1, 1.])), rms(x))
    print('model dips before:', model_poles(P0)); print('model dips after :', model_poles(P))
    print('model eigen (n0 cut) before:', model_eig(P0), ' after:', model_eig(P), ' HFSS:', EIG)
    out = os.path.join(RES, 'fit_' + '_'.join(tags) + ('_eig' if EIG else ''))
    json.dump(dict(tags=tags, eig_hfss=EIG, params=dict(zip(NAMES, map(float, x))), dips_after=model_poles(P),
                   dips_before=model_poles(P0), eig_before=model_eig(P0), eig_after=model_eig(P), cost=float(sol.cost)),
              open(out + '.json', 'w'), indent=1)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(len(datasets), 1, figsize=(10, 3.6 * len(datasets)), squeeze=False)
    for a, (P_, f, s11, s21), t in zip(ax[:, 0], datasets, tags):
        a.plot(f, 20 * np.log10(np.abs(s21)), label='HFSS S21')
        a.plot(f, 20 * np.log10(np.abs(s21_theory(f, P_)[1])).clip(-90), ':', label='theory (Q3D, before)')
        a.plot(f, 20 * np.log10(np.abs(s21_theory(f, apply(P_, x))[1])).clip(-90), '--', label='theory (fit)')
        a.set_title(t); a.set_ylabel('S21 (dB)'); a.grid(alpha=.3); a.legend()
    ax[-1, 0].set_xlabel('f (GHz)'); fig.tight_layout(); fig.savefig(out + '.png', dpi=120)
    print('saved', out + '.png')


if __name__ == '__main__':
    args = sys.argv[1:]
    if '--eig' in args:
        k = args.index('--eig'); EIG[:] = [float(s) for s in args[k + 1].split(',')]; del args[k:k + 2]
    main(args or ['c6_T3'])
