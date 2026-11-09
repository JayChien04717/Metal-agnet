"""Wallraff-group readout circuit model (thesis sec. 2.5, Fig. 2.29, eqs. 2.122-2.137; screenshots + notes in
SQL HFSS/note/wallraff reference/) and a Design6 version with their parameters.

  1. reproduce Table 2.4: S21 of each channel (shunt Z_rc on a matched 50 Ohm line), fit eq. 2.137 -> omega_p,
     omega_r, J, kappa_p; the CPW phase velocity v of their chip is not given -> one global v is calibrated
     (C_AB air bridges set to 0, absorbed in v).
  2. Design6 'W1': same circuit, our CPW (Z0 50.9 Ohm, v 1.2016e8), no qubit (resonator end = open end),
     Wallraff Q35 capacitances, L_pf = 1.2 mm, L_rr = 1.3 mm; L_pfg, L_rrg solved for omega_r = 6.700 GHz,
     omega_p = omega_r + 19 MHz (Q35: 6690 - 6671).

Usage: python wallraff_model.py   -> results/wallraff_model.json, results/wallraff_model.png
"""
import os, json, warnings
warnings.filterwarnings('ignore')
import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.signal import find_peaks
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results')

Q = ['Q36', 'Q35', 'Q34', 'Q33', 'Q32']
COL = dict(Cpf=[48.57, 45.47, 37.08, 45.10, 37.55], Cpff=[-5.72, -5.88, -6.23, -5.90, -6.21],
           Cpfg=[98.80, 93.33, 78.11, 92.56, 79.00], Cpr=[2.50, 2.77, 2.42, 2.5, 2.66], Cprg=[1.95, 2.39, 1.83, 1.95, 2.21],
           Cqr=[12.26, 17.95, 12.26, 17.95, 12.26], Cqrg=[17.61, 27.03, 17.61, 27.04, 17.61],
           Cqg=[105.36, 99.35, 105.36, 99.35, 105.36],
           Lpf=[1200] * 5, Lpfg=[1920.45, 2268.05, 2202.36, 1844.14, 2296.28], Lpfe=[2.625] * 5,
           Lrr=[1303.6, 1300, 1303.6, 1300, 1303.6], Lrrg=[2489.65, 2749.33, 2603.73, 2316.18, 2720.55], Lrre=[2.625] * 5)
SIM = dict(wp=[7250, 6690, 7050, 7440, 6850], wr=[7238, 6671, 7039, 7417, 6839], J=[15.84, 15.04, 15.80, 15.55, 15.41],
           kp=[35.37, 40.98, 52.07, 63.16, 59.91])
MEAS = dict(wp=[7258, 6705, 7056, 7451, 6865], wr=[7227, 6655, 7018, 7394, 6821], J=[20.54, 18.20, 20.59, 16.27, 17.63],
            kp=[55.58, 28.43, 53.61, 47.01, 60.41])
TABLE = {q: {k: v[i] for k, v in COL.items()} for i, q in enumerate(Q)}


def Y_rc(f_mhz, p, Z0=50.0, v=1.2e8, CAB=0.0, qubit=True):
    """Admittance of the resonator + Purcell filter branch seen at the feed-line t-junction (eqs. 2.122-2.134).
    Capacitances in fF, lengths in um. qubit=False: resonator end is a plain open end (Design6)."""
    w = 2 * np.pi * np.asarray(f_mhz, complex) * 1e6; Y0 = 1 / Z0; b = w / v
    jC = lambda c: 1j * w * c * 1e-15
    t = lambda l: np.tan(b * l * 1e-6)
    line = lambda l, YL: Y0 * (YL + 1j * Y0 * t(l)) / (Y0 + 1j * YL * t(l))          # line of length l loaded by YL
    short = lambda l: -1j * Y0 / t(l)
    ser = lambda Ya, Yb: Ya * Yb / (Ya + Yb)
    # resonator, qubit side (Z_rq,1): qubit = C_qr in series with C_qg; + C_qrg + 2 C_AB; line L_rr; + C_AB
    Yq = (ser(jC(p['Cqr']), jC(p['Cqg'])) + jC(p['Cqrg'])) if qubit else 0 * w
    Yrq1 = line(p['Lrr'], Yq + jC(2 * CAB)) + jC(CAB)
    # resonator, ground side (Z_rq,2): L_rre short, 2 C_AB, L_rrg
    Yrq2 = line(p['Lrrg'], short(p['Lrre']) + jC(2 * CAB))
    # hammer coupler (Z_pr)
    Ypr = ser(jC(p['Cpr']), Yrq1 + Yrq2 + jC(CAB + p['Cprg'])) + jC(p['Cprg'] + CAB)
    # Purcell filter, from the shorted end back to the feed line
    Yf = line(p['Lpfg'], short(p['Lpfe']) + jC(CAB)) + jC(2 * CAB) + Ypr
    Yf = line(p['Lpf'], Yf) + jC(p['Cpfg'] + CAB)
    return ser(jC(p['Cpf']), Yf) + jC(p['Cpff'] + CAB)


def s21(f_mhz, p, Z0=50.0, **kw):
    return 2 / (2 + Z0 * Y_rc(f_mhz, p, Z0=Z0, **kw))


def io_model(f, x):
    """eq. 2.137 (complex), f and x in MHz; x = wp, wr, J, kp, A, phi, |B|, arg B."""
    wp, wr, J, kp, A, phi, B, pb = x
    return B * np.exp(1j * pb) * (1 - A * np.exp(1j * phi) * 2j * kp * (wr - f) / (4 * J ** 2 + (2j * (wp - f) + kp) * 2j * (wr - f)))


def fit_io(f, s):
    """Fit eq. 2.137 to a complex S21 trace. Start: two deepest dips d1 < d2, wr = transmission peak between them,
    wp = d1 + d2 - wr (sum of the hybrid modes)."""
    a = np.abs(s); i, _ = find_peaks(-a, prominence=0.05)
    i = sorted(i, key=lambda k: a[k])[:2]
    if len(i) == 2:
        d1, d2 = sorted(f[i]); m = (f > d1) & (f < d2)
        wr0 = f[m][np.argmax(a[m])] if m.any() else (d1 + d2) / 2; wp0 = d1 + d2 - wr0
        J0 = max(np.sqrt(max((d2 - d1) ** 2 - (wp0 - wr0) ** 2, 1)) / 2, 3.0)
    else:
        wr0 = wp0 = f[np.argmin(a)]; J0 = 15.0
    best = None
    for kp0 in (20, 40, 80):
        x0 = [wp0, wr0, J0, kp0, 1, 0, np.abs(s[0]), np.angle(s[0])]
        r = least_squares(lambda x: np.concatenate([(io_model(f, x) - s).real, (io_model(f, x) - s).imag]), x0,
                          bounds=([f[0], f[0], 0, 1, 0, -np.pi, 0, -2 * np.pi], [f[-1], f[-1], 300, 1000, 2, np.pi, 2, 2 * np.pi]))
        best = r if best is None or r.cost < best.cost else best
    return dict(zip(['wp', 'wr', 'J', 'kp'], map(float, best.x[:4])), rms=float(np.sqrt(2 * best.cost / len(f))), x=best.x)


def croot(fun, z, n=60):
    """Complex secant root of fun near z (MHz)."""
    z0, z1 = complex(z), complex(z) + 0.5 + 0.5j
    f0, f1 = fun(z0), fun(z1)
    for _ in range(n):
        if f1 == f0: break
        z2 = z1 - f1 * (z1 - z0) / (f1 - f0); z0, f0, z1, f1 = z1, f1, z2, fun(z2)
        if abs(z1 - z0) < 1e-7: break
    return z1


def extract(p, v, center, span=600, Z0=50.0, **kw):
    """omega_p, omega_r, J, kappa_p of eq. 2.137 directly from the circuit (no curve fit):
    omega_r = real zero of Y_rc (S21 = B there, the transmission peak between the dips);
    the two complex poles l1, l2 of S21 (zeros of 2 + Z0 Y_rc) are the eigenvalues of [[wp - i kp/2, J], [J, wr]]
    -> wp = Re(l1 + l2) - wr, kp = 2 |Im(l1 + l2)|, J^2 = (wp - i kp/2) wr - l1 l2."""
    f = np.linspace(center - span / 2, center + span / 2, 60001)
    D = lambda z: 2 + Z0 * Y_rc(z, p, Z0=Z0, v=v, **kw)
    a = np.abs(2 / D(f)); i, _ = find_peaks(-a, prominence=1e-3)
    i = sorted(i, key=lambda k: a[k])[:2]
    if len(i) < 2:
        return None, f
    l = sorted([croot(D, f[k]) for k in i], key=lambda z: z.real)
    B = Y_rc(f, p, Z0=Z0, v=v, **kw).imag; m = (f > l[0].real - 30) & (f < l[1].real + 30)
    zc = np.where(np.sign(B[m][:-1]) != np.sign(B[m][1:]))[0]
    fm = f[m]; Bm = B[m]
    zeros = [fm[k] - Bm[k] * (fm[k + 1] - fm[k]) / (Bm[k + 1] - Bm[k]) for k in zc if abs(Bm[k]) < 1 and abs(Bm[k + 1]) < 1]
    wr = min(zeros, key=lambda z: abs(z - (l[0].real + l[1].real) / 2)) if zeros else np.nan
    if any(abs(z - f[k]) > 100 for z, k in zip(l, sorted(i, key=lambda k: f[k]))):
        return None, f                                        # secant ran away
    s = l[0] + l[1]; kp = 2 * abs(s.imag); wp = s.real - wr
    J2 = (wp + 1j * s.imag) * wr - l[0] * l[1]               # poles: l1 + l2 = wp + wr + i Im(s), Im(s) = +-kp/2
    return dict(wp=float(wp), wr=float(wr), J=float(np.sqrt(abs(J2.real))), kp=float(kp), J2_imag=float(J2.imag),
                modes=[(float(z.real), float(2 * abs(z.imag))) for z in l]), f


def calibrate(C_AB_free=True):
    """Phase velocity (and optionally one air-bridge C_AB) for all five channels, matched to omega_r^Sim and
    omega_p^Sim of Table 2.4. J and kappa_p are then an independent check."""
    def res(x):
        v, cab = x[0] * 1e8, (x[1] if C_AB_free else 0.0)
        e = []
        for i, q in enumerate(Q):
            r, _ = extract(TABLE[q], v, center=(SIM['wr'][i] + SIM['wp'][i]) / 2, CAB=cab)
            e += [r['wr'] - SIM['wr'][i], r['wp'] - SIM['wp'][i]] if r else [500, 500]
        return e
    x0 = [1.15, 1.0] if C_AB_free else [1.15]
    s = least_squares(res, x0, bounds=([1.0, 0.0], [1.3, 10.0]) if C_AB_free else ([1.0], [1.3]), diff_step=1e-3)
    return s.x[0] * 1e8, (s.x[1] if C_AB_free else 0.0), np.array(res(s.x))


def design_W1(target_r=6700.0, dp=19.0, Z0=50.9, v=1.2016e8, base='Q35', open_end_C=0.0):
    """Design6 version: Wallraff caps of `base`, our CPW, resonator end open (optional lumped C there),
    L_pfg and L_rrg solved for omega_r = target_r and omega_p = target_r + dp (MHz)."""
    p = dict(TABLE[base]); p['Cqrg'] = open_end_C
    # start: bare lambda/4 lengths minus the capacitive loading, then refine on the extracted wp, wr
    lam4 = v / (4 * target_r * 1e6) * 1e6                         # um
    x0 = [lam4 - p['Lpf'] - p['Lpfe'] - Z0 * v * (p['Cpfg'] + p['Cpf']) * 1e-9, lam4 - p['Lrr'] - p['Lrre'] - Z0 * v * open_end_C * 1e-9]
    def res(x):
        p['Lpfg'], p['Lrrg'] = x
        r, _ = extract(p, v, center=target_r + dp / 2, Z0=Z0, qubit=False)
        return [r['wr'] - target_r, r['wp'] - target_r - dp] if r else [300, 300]
    for _ in range(3):                                            # coarse fixed-point on the bare frequencies
        e = res(x0); x0 = [x0[0] * (1 + e[1] / target_r), x0[1] * (1 + e[0] / target_r)]
    sol = least_squares(res, x0, x_scale=[20, 20], diff_step=1e-5)
    p['Lpfg'], p['Lrrg'] = map(float, sol.x)
    r, f = extract(p, v, center=target_r + dp / 2, Z0=Z0, qubit=False)
    return p, r, f


if __name__ == '__main__':
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    out = {}
    fig, axs = plt.subplots(2, 1, figsize=(11, 9))
    for free in (False,):              # v + C_AB together is degenerate (C_AB stays at 0): v only
        v_w, cab, err = calibrate(C_AB_free=free)
        print(f'calibration {"v + C_AB" if free else "v only, C_AB = 0"}: v = {v_w:.4e} m/s, C_AB = {cab:.2f} fF, '
              f'rms(wr, wp error) = {np.sqrt(np.mean(err ** 2)):.1f} MHz')
        out[f'cal_{"v_cab" if free else "v"}'] = dict(v=v_w, C_AB=cab, err=err.tolist())
    print(f'{"":4s} {"wp model/tab":>14s} {"wr model/tab":>14s} {"J model/tab":>13s} {"kp model/tab":>13s}')
    out['table'] = {}
    for i, q in enumerate(Q):
        c = (SIM['wr'][i] + SIM['wp'][i]) / 2
        r, f = extract(TABLE[q], v_w, center=c, CAB=cab)
        out['table'][q] = r
        print(f'{q:4s} {r["wp"]:7.1f}/{SIM["wp"][i]:<6d} {r["wr"]:7.1f}/{SIM["wr"][i]:<6d} {r["J"]:5.2f}/{SIM["J"][i]:<6.2f} '
              f'{r["kp"]:5.2f}/{SIM["kp"][i]:<6.2f}  (Im J^2 {r["J2_imag"]:.1e})')
        fa = np.linspace(c - 150, c + 150, 6001)
        axs[0].plot(fa, np.abs(s21(fa, TABLE[q], v=v_w, CAB=cab)),
                    label=f'{q}: J {r["J"]:.1f} (tab {SIM["J"][i]}), κp {r["kp"]:.1f} (tab {SIM["kp"][i]}) MHz')
    axs[0].set_title(f'Wallraff Table 2.4 channels, circuit model eqs. 2.122-2.136 (v = {v_w:.4e} m/s, C_AB = {cab:.2f} fF)')
    axs[0].set_xlabel('f (MHz)'); axs[0].set_ylabel('|S21|'); axs[0].legend(fontsize=8); axs[0].grid(alpha=.3)
    # Design6 with Wallraff parameters (our CPW, no qubit, no air bridges)
    p, r, f = design_W1()
    out['W1'] = dict(params=p, modes=r, L_filter_total_mm=(p['Lpf'] + p['Lpfg'] + p['Lpfe']) / 1000,
                     L_res_total_mm=(p['Lrr'] + p['Lrrg'] + p['Lrre']) / 1000)
    print(f'W1 (Design6 CPW, Q35 caps, open resonator end): L_pfg {p["Lpfg"]:.1f} um, L_rrg {p["Lrrg"]:.1f} um -> '
          f'filter {out["W1"]["L_filter_total_mm"]:.4f} mm, resonator {out["W1"]["L_res_total_mm"]:.4f} mm')
    print(f'   wp {r["wp"]:.1f}  wr {r["wr"]:.1f}  J {r["J"]:.2f}  kp {r["kp"]:.2f} MHz; hybrid modes (f, width) {r["modes"]}')
    fw = np.linspace(6550, 6850, 6001); s = s21(fw, p, Z0=50.9, v=1.2016e8, qubit=False)
    axs[1].plot(fw, 20 * np.log10(np.abs(s)), label='W1 circuit model')
    x = [r['wp'], r['wr'], r['J'], r['kp'], 1, 0, 1, 0]
    axs[1].plot(fw, 20 * np.log10(np.abs(io_model(fw, x))), ':', label='eq. 2.137 with the extracted wp, wr, J, κp')
    a = 20 * np.log10(np.abs(s)); i, _ = find_peaks(-a, prominence=3); out['W1']['dips_MHz'] = [float(fw[k]) for k in i]
    print('W1 dips (MHz):', [round(float(fw[k]), 1) for k in i])
    axs[1].axvline(6700, color='k', lw=.6, ls='--')
    axs[1].set_title('Design6 W1: Wallraff Q35 capacitances + Design6 CPW, target ω_r = 6.700 GHz, ω_p = ω_r + 19 MHz')
    axs[1].set_xlabel('f (MHz)'); axs[1].set_ylabel('S21 (dB)'); axs[1].legend(); axs[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(RES, 'wallraff_model.png'), dpi=120)
    json.dump(out, open(os.path.join(RES, 'wallraff_model.json'), 'w'), indent=1, default=float)
    print('saved', os.path.join(RES, 'wallraff_model.png'))
