"""Re-solve the W1 filter / resonator lengths with the Q3D pi networks of the drawn W1 capacitors (theory_abcd model,
raw Q3D, no fit factors) so that omega_r = 6.700 GHz and omega_p = omega_r + 19 MHz (Wallraff Q35 detuning).

omega_p, omega_r, J, kappa_p are extracted from the branch admittance seen at n0 exactly as in wallraff_model.extract
(single channel on a matched line: S21 = 2 / (2 + Z0 Y_branch); omega_r = real zero of Y, poles = hybrid modes).
Then the full-channel S21 (C_in, feed lines) is computed for the solved lengths.

Usage: D6_VARIANT=W1 python solve_w1.py    -> prints DLF_EXTRA / DLR_EXTRA, results/solve_w1.json, results/solve_w1.png
"""
import os, sys, json, copy, warnings
warnings.filterwarnings('ignore')
os.environ.setdefault('D6_VARIANT', 'W1')
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results'); sys.path.insert(0, HERE)
import theory_abcd as T


def croot(fun, z, h=0.02, n=80):
    """Complex secant root near z (MHz) with a small first step h (MHz) and step limiting (narrow modes)."""
    z0, z1 = complex(z), complex(z) + h * (1 + 1j)
    f0, f1 = fun(z0), fun(z1)
    for _ in range(n):
        if f1 == f0: break
        dz = -f1 * (z1 - z0) / (f1 - f0)
        if abs(dz) > 5: dz *= 5 / abs(dz)
        z0, f0, z1 = z1, f1, z1 + dz; f1 = fun(z1)
        if abs(dz) < 1e-8: break
    return z1

TARGET_R, DP = 6.700, 0.019          # GHz


def extract(P, center, span=0.6):
    Z0 = P['Z0']
    Y = lambda z: T.branch_Y(z, P)
    D = lambda z: 2 + Z0 * Y(z)
    f = np.linspace(center - span / 2, center + span / 2, 60001)
    a = np.abs(2 / D(f)); i, _ = find_peaks(-a, prominence=1e-3)
    i = sorted(sorted(i, key=lambda k: a[k])[:2], key=lambda k: f[k])
    if len(i) < 2:
        return None
    l = [croot(lambda z: D(z / 1e3), f[k] * 1e3) / 1e3 for k in i]      # croot works in MHz
    if any(abs(z - f[k]) > 0.1 for z, k in zip(l, i)):
        return None
    B = Y(f).imag; m = (f > l[0].real - 0.03) & (f < l[1].real + 0.03); fm, Bm = f[m], B[m]
    zc = [k for k in np.where(np.sign(Bm[:-1]) != np.sign(Bm[1:]))[0] if abs(Bm[k]) < 1 and abs(Bm[k + 1]) < 1]
    zeros = [fm[k] - Bm[k] * (fm[k + 1] - fm[k]) / (Bm[k + 1] - Bm[k]) for k in zc]
    if not zeros:
        return None
    wr = min(zeros, key=lambda z: abs(z - (l[0].real + l[1].real) / 2))
    s = l[0] + l[1]; wp = s.real - wr
    J2 = (wp + 1j * s.imag) * wr - l[0] * l[1]
    return dict(wp=float(wp), wr=float(wr), J=float(np.sqrt(abs(J2.real))) * 1e3, kp=float(2 * abs(s.imag)) * 1e3,
                J2_imag=float(J2.imag) * 1e6, modes=[(float(z.real), float(2 * abs(z.imag)) * 1e3) for z in l])


def with_dl(P0, dlf, dlr):
    P = copy.deepcopy(P0); P['L_filt'] = P0['L_filt'] + dlf; P['L_res'] = P0['L_res'] + dlr; return P


if __name__ == '__main__':
    P0 = T.params_from_layout()
    print('W1 params (Q3D):', {k: (tuple(round(x, 2) for x in v) if isinstance(v, tuple) else round(v, 4)) for k, v in P0.items()})
    r0 = extract(P0, TARGET_R + DP / 2, span=1.2); print('as drawn:', r0)

    def res(x):
        r = extract(with_dl(P0, *x), TARGET_R + DP / 2, span=1.0)
        return [(r['wr'] - TARGET_R) * 1e3, (r['wp'] - TARGET_R - DP) * 1e3] if r else [500, 500]
    # start: shift each line by its frequency error (dL = -L df / f)
    x = [0.0, 0.0]
    if r0:
        x = [-P0['L_filt'] * (TARGET_R + DP - r0['wp']) / TARGET_R, -P0['L_res'] * (TARGET_R - r0['wr']) / TARGET_R]
    sol = least_squares(res, x, x_scale=[0.01, 0.01], diff_step=1e-4)
    dlf, dlr = map(float, sol.x); P = with_dl(P0, dlf, dlr); r = extract(P, TARGET_R + DP / 2)
    print(f'DLF_EXTRA = {dlf:+.4f} mm (filter {P["L_filt"]:.4f}), DLR_EXTRA = {dlr:+.4f} mm (resonator {P["L_res"]:.4f})')
    print('solved:', r)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    f = np.linspace(6.4, 7.0, 12001)
    s_shunt = 2 / (2 + P['Z0'] * T.branch_Y(f, P)); _, s_full = T.s21_theory(f, P)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(f, 20 * np.log10(np.abs(s_shunt)), label='branch only (matched line)')
    ax.plot(f, 20 * np.log10(np.abs(s_full)), '--', label='full channel (C_in, feed lines)')
    ax.axvline(TARGET_R, color='k', lw=.6, ls=':')
    ax.set_title(f'W1 with Q3D C: wr {r["wr"]:.4f}, wp {r["wp"]:.4f} GHz, J {r["J"]:.1f}, κp {r["kp"]:.1f} MHz')
    ax.set_xlabel('f (GHz)'); ax.set_ylabel('S21 (dB)'); ax.grid(alpha=.3); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(RES, 'solve_w1.png'), dpi=120)
    dips = lambda s: [round(float(f[k]), 4) for k in find_peaks(-20 * np.log10(np.abs(s)), prominence=3)[0]]
    print('dips branch-only', dips(s_shunt), ' full channel', dips(s_full))
    json.dump(dict(DLF_EXTRA=dlf, DLR_EXTRA=dlr, L_filt=P['L_filt'], L_res=P['L_res'], as_drawn=r0, solved=r,
                   cp=list(P0['cp']), cj=list(P0['cj']), dips_full=dips(s_full)), open(os.path.join(RES, 'solve_w1.json'), 'w'), indent=1)
    print('saved', os.path.join(RES, 'solve_w1.png'))
