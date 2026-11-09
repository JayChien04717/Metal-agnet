"""Design6 circuit model of the verification layouts (no C_in: the branch sits on a matched feed line).

  S21 = e^{-j beta l} * 2 / (2 + Z0 Y_branch)          (l = feed line pad to pad; only a phase)
  Y_branch = stub - C_p pi - filter (D_cj to the C_j tap, then shorted) - C_j pi - resonator (open / shorted ends)

VARIANT filter: no resonator lines; the C_j hammer stays, its far plate floating -> the tap sees C_jg1 + (C_j series C_jg2).
Predictions straight from the capacitances (no fitting):
  filter : omega_p, kappa_p = pole of 2/Z0 + Y_filter  (the Purcell <-> feed-line coupling, set by C_p)
  full   : omega_p, omega_r, J, kappa_p from modes5.extract (J = resonator <-> Purcell coupling, set by C_j)

Capacitances: params(source) - 'design5' = Design5 whole-channel Q3D (with C_in, same hugger / hammer geometry),
or a Design6 Q3D json (results/<date>/<HHMM>_q3d6_<variant>/q3d6_caps.json) once that run exists.
Usage: python verify_model.py       -> prints the predictions, results/verify_model.png
"""
import os, sys, json, warnings
warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results')
HFSS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(HFSS, 'Design5', 'simulation'))
import theory5 as T5, modes5 as M5
LEN_JSON = os.path.join(HFSS, 'Design5', 'simulation', 'results', 'chip6_lengths.json')
CH = 3


def params(caps_json=None):
    """W1 parameters with the ch3 (6.7 GHz) lengths. caps_json: Design6 Q3D result overriding the series C."""
    P = T5.params_from_layout()
    c = {c['ch']: c for c in json.load(open(LEN_JSON))['channels']}[CH]
    P['L_filt'], P['L_res'] = c['L_filter'], c['L_res']
    if caps_json:
        q = json.load(open(caps_json))
        P['cp'] = (q['series']['cp'],) + tuple(P['cp'][1:]); P['cj'] = (q['series']['cj'],) + tuple(P['cj'][1:])
        P['cap_source'] = os.path.relpath(caps_json, HERE)
    return P


def branch_Y(f_ghz, P, resonator=True):
    """theory5.branch_Y with the option of no resonator lines (VARIANT filter). f may be complex."""
    if resonator:
        return T5.branch_Y(f_ghz, P)
    Z0, v = P['Z0'], P['v']
    w = 2 * np.pi * np.asarray(f_ghz) * 1e9; beta = w / v; Y0 = 1 / Z0
    t = lambda l: np.tan(beta * l * 1e-3)
    cjs, cjg1, cjg2 = P['cj']
    jw = lambda c: 1j * w * c * 1e-15
    Yres = jw(cjg2)                                               # floating far plate of the hammer
    Yser = jw(cjs) * Yres / (jw(cjs) + Yres)
    Ytap = jw(cjg1) + Yser - 1j * Y0 / t(P['L_filt'] - P['D_cj'])
    td = t(P['D_cj']); Yf = Y0 * (Ytap + 1j * Y0 * td) / (Y0 + 1j * Ytap * td)
    cps, cpg_a, cpg_b = P['cp']
    Y1 = jw(cpg_b) + Yf; Y2 = jw(cps) * Y1 / (jw(cps) + Y1); Y3 = jw(cpg_a) + Y2
    ts = t(P['l_stub'])
    return Y0 * (Y3 + 1j * Y0 * ts) / (Y0 + 1j * Y3 * ts)


def s21(f_ghz, P, resonator=True):
    """|S21| shape on a matched line (feed-line phase left out; the fit takes a free phase / delay)."""
    return 2 / (2 + P['Z0'] * branch_Y(f_ghz, P, resonator))


def s11(f_ghz, P, resonator=True):
    Y = branch_Y(f_ghz, P, resonator); return -P['Z0'] * Y / (2 + P['Z0'] * Y)


def predict(P):
    """Couplings predicted from the capacitances."""
    Z0 = P['Z0']
    D = lambda z: 2 + Z0 * branch_Y(z / 1e3, P, resonator=False)
    f = np.linspace(6.4, 7.0, 60001); i = int(np.argmin(np.abs(D(f * 1e3))))
    zp = M5.croot(D, f[i] * 1e3) / 1e3
    filt = dict(wp=float(zp.real), kp=float(2 * abs(zp.imag)) * 1e3)
    # full channel: centre the extraction on the two S21 dips of this model (they move with the fitted lengths / C)
    from scipy.signal import find_peaks
    a = -20 * np.log10(np.abs(s21(f, P, True)))
    i, _ = find_peaks(a, prominence=3, distance=40); i = sorted(sorted(i, key=lambda k: -a[k])[:2])
    full = None
    if len(i) == 2:
        c, sep = (f[i[0]] + f[i[1]]) / 2, f[i[1]] - f[i[0]]
        full = M5.extract(lambda z: branch_Y(z, P), Z0, c, span=max(0.2, 6 * sep))
    return dict(filter=filt, full=full)


if __name__ == '__main__':
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    P = params(); pr = predict(P)
    print(f'caps ({P["cap_source"]}): C_p {P["cp"][0]:.2f} fF, C_j {P["cj"][0]:.2f} fF')
    print(f'filter only: omega_p {pr["filter"]["wp"]:.4f} GHz, kappa_p {pr["filter"]["kp"]:.1f} MHz')
    fu = pr['full']; print(f'full channel: omega_p {fu["wp"]:.4f}, omega_r {fu["wr"]:.4f} GHz, J {fu["J"]:.1f} MHz, kappa_p {fu["kp"]:.1f} MHz')
    f = np.linspace(6.55, 6.90, 7001); db = lambda s: 20 * np.log10(np.abs(s))
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(f, db(s21(f, P, False)), label=f'filter only: ω_p {pr["filter"]["wp"]:.4f} GHz, κ_p {pr["filter"]["kp"]:.1f} MHz')
    ax.plot(f, db(s21(f, P, True)), label=f'full: J {fu["J"]:.1f} MHz, ω_r {fu["wr"]:.4f} GHz')
    ax.set_xlabel('f (GHz)'); ax.set_ylabel('|S21| (dB)'); ax.grid(alpha=.3); ax.legend(fontsize=8)
    ax.set_title(f'Design6 verification, model with C_p {P["cp"][0]:.1f} fF, C_j {P["cj"][0]:.2f} fF (no C_in)', fontsize=9)
    os.makedirs(RES, exist_ok=True); fig.tight_layout(); fig.savefig(os.path.join(RES, 'verify_model.png'), dpi=110)
