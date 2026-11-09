"""omega_p, omega_r, J, kappa_p of the 2-mode model (Wallraff thesis eq. 2.137) taken directly from a circuit model,
and the direct fit of eq. 2.137 to an S21 trace.

Circuit extraction (no curve fit): single channel on a matched line, S21 = 2 / (2 + Z0 Y_branch)
  omega_r   = real zero of Y_branch (S21 = background there: the peak between the two dips)
  l1, l2    = complex poles of S21 = eigenvalues of [[wp + i kp/2, J], [J, wr]]   (e^{+jwt}: Im > 0)
  wp = Re(l1 + l2) - wr,  kp = 2 |Im(l1 + l2)|,  J^2 = (wp + i Im(l1 + l2)) wr - l1 l2   (Im J^2 ~ 0 = consistency)
Frequencies in GHz, J and kappa in MHz.
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks


def croot(fun, z, h=0.02, n=80):
    """Complex secant root near z (MHz units), small first step, step limited (narrow modes)."""
    z0, z1 = complex(z), complex(z) + h * (1 + 1j)
    f0, f1 = fun(z0), fun(z1)
    for _ in range(n):
        if f1 == f0: break
        dz = -f1 * (z1 - z0) / (f1 - f0)
        if abs(dz) > 5: dz *= 5 / abs(dz)
        z0, f0, z1 = z1, f1, z1 + dz; f1 = fun(z1)
        if abs(dz) < 1e-8: break
    return z1


def extract(Yfun, Z0, center, span=0.6, Yenv=None):
    """Yenv: admittance of the surroundings at the branch node (default matched line, 2/Z0). With the in-channel
    environment (C_in on one side) kappa_p is the one a fit of the channel S21 sees."""
    Ye = Yenv or (lambda z: 2 / Z0 + 0 * np.asarray(z))
    D = lambda z: Z0 * (Ye(z) + Yfun(z))
    f = np.linspace(center - span / 2, center + span / 2, 60001)
    a = np.abs(Z0 * Ye(f) / D(f)); i, _ = find_peaks(-a, prominence=1e-3, distance=200)
    i = sorted(sorted(i, key=lambda k: a[k])[:2], key=lambda k: f[k])
    if len(i) < 2:
        return None
    l = [croot(lambda z: D(z / 1e3), f[k] * 1e3) / 1e3 for k in i]
    if any(abs(z - f[k]) > 0.1 for z, k in zip(l, i)):
        return None
    B = Yfun(f).imag; m = (f > l[0].real - 0.03) & (f < l[1].real + 0.03); fm, Bm = f[m], B[m]
    zc = [k for k in np.where(np.sign(Bm[:-1]) != np.sign(Bm[1:]))[0] if abs(Bm[k]) < 1 and abs(Bm[k + 1]) < 1]
    zeros = [fm[k] - Bm[k] * (fm[k + 1] - fm[k]) / (Bm[k + 1] - Bm[k]) for k in zc]
    if not zeros:
        return None
    wr = min(zeros, key=lambda z: abs(z - (l[0].real + l[1].real) / 2))
    s = l[0] + l[1]; wp = s.real - wr
    J2 = (wp + 1j * s.imag) * wr - l[0] * l[1]
    return dict(wp=float(wp), wr=float(wr), J=float(np.sqrt(abs(J2.real))) * 1e3, kp=float(2 * abs(s.imag)) * 1e3,
                J2_imag_MHz2=float(J2.imag) * 1e6, modes=[(float(z.real), float(2 * abs(z.imag)) * 1e3) for z in l])


# ---------------- direct fit of eq. 2.137 ----------------
NAMES = ['wp', 'wr', 'J', 'kp', 'A', 'phi', 'B', 'pb', 'tau', 'b1']


def io_model(f, x):
    """eq. 2.137 with a background: B e^{i(pb - 2 pi f tau)} (1 + b1 (f - f_c)).  f in GHz, J / kp in MHz, tau in ns."""
    wp, wr, J, kp, A, phi, B, pb, tau, b1 = x
    fc = (wp + wr) / 2; J, kp = J * 1e-3, kp * 1e-3
    core = 1 - A * np.exp(1j * phi) * 2j * kp * (wr - f) / (4 * J ** 2 + (2j * (wp - f) + kp) * 2j * (wr - f))
    return B * np.exp(1j * (pb - 2 * np.pi * f * tau)) * (1 + b1 * (f - fc)) * core


def _dips(f, s21):
    """Two deepest S21 dips (at least 2 MHz apart, robust against split nulls)."""
    a = 20 * np.log10(np.abs(s21)); df = f[1] - f[0]
    i, _ = find_peaks(-a, prominence=3, distance=max(1, int(0.002 / df)))
    return sorted(sorted(i, key=lambda k: a[k])[:2], key=lambda k: f[k])


def direct_fit(f, s21, window=None):
    """Fit eq. 2.137 (+ background) to complex S21, staged:
    1. background B e^{i(pb - 2 pi f tau)} (1 + b1 (f - fc)) on the data outside the resonance region,
    2. resonance (wp, wr, J, kp, A, phi) with the background fixed, several starts,
    3. all parameters together from the best start.
    HFSS and the ABCD model use e^{+jwt}; eq. 2.137 is written for e^{-iwt} -> the fit is done on conj(S21)
    (otherwise the resonance circle turns the wrong way and kappa runs negative). The returned model x gives
    conj(S21): plot np.conj(io_model(f, x))."""
    s21 = np.conj(np.asarray(s21))
    i = _dips(f, s21); a = np.abs(s21)
    if len(i) == 2:
        d1, d2 = f[i[0]], f[i[1]]; m = (f > d1) & (f < d2)
        wr0 = f[m][np.argmax(a[m])] if m.sum() > 2 else (d1 + d2) / 2; wp0 = d1 + d2 - wr0
        J0 = max(np.sqrt(max((d2 - d1) ** 2 - (wp0 - wr0) ** 2, 1e-8)) / 2 * 1e3, 2.0)
    else:
        d1 = d2 = wr0 = wp0 = f[i[0]] if len(i) else f[np.argmin(a)]; J0 = 15.0
    sep = max(d2 - d1, 0.01); fc = (d1 + d2) / 2
    w = window or max(0.15, 5 * sep)
    sel = (f > fc - w) & (f < fc + w); ff, ss = f[sel], s21[sel]
    # 1. background from |f - fc| > 2.5 sep (and > 60 MHz)
    out = np.abs(ff - fc) > max(2.5 * sep, 0.06)
    if out.sum() < 10: out = np.abs(ff - fc) > 0.6 * w
    ph = np.unwrap(np.angle(ss[out])); tau0 = -np.polyfit(ff[out], ph, 1)[0] / (2 * np.pi)
    def bg(xb, x):
        B, pb, tau, b1 = xb
        return B * np.exp(1j * (pb - 2 * np.pi * x * tau)) * (1 + b1 * (x - fc))
    xb0 = [float(np.median(np.abs(ss[out]))), 0.0, tau0, 0.0]
    xb0[1] = float(np.angle(np.mean(ss[out] / bg(xb0, ff[out]))))
    rb = least_squares(lambda xb: np.r_[(bg(xb, ff[out]) - ss[out]).real, (bg(xb, ff[out]) - ss[out]).imag], xb0)
    B, pb, tau, b1 = rb.x
    # io_model uses (1 + b1 (f - (wp + wr)/2)); keep fc fixed by rewriting b1 later in the joint fit
    norm = ss / bg(rb.x, ff)
    # reparametrised by the two S21 zeros z1, z2 (lossless: (wr - f)(wp - f) = J^2 there), which the data pin down
    # sharply:  wp = z1 + z2 - wr,  J^2 = wp wr - z1 z2
    def to_x(y):
        z1, z2, wr, kp, A, phi = y
        wp = z1 + z2 - wr; J2 = wp * wr - z1 * z2
        return [wp, wr, np.sqrt(max(J2, 0.0)) * 1e3, kp, A, phi]
    core = lambda f_, y: io_model(f_, to_x(y) + [1, 0, 0, 0])
    z1, z2 = d1, d2
    # 2. resonance only, background fixed; wr scanned between the zeros (Fano can move the visible peak)
    best = None
    lo2 = [z1 - 0.003, z2 - 0.003, z1, 0.1, 0, -np.pi]; hi2 = [z1 + 0.003, z2 + 0.003, z2, 3000, 3, np.pi]
    for wr_s in np.unique(np.r_[wr0, np.linspace(z1, z2, 7)[1:-1]]):
        for kp0 in (5, 15, 40, 100, 250):
            for phi0 in (-0.5, 0.0, 0.5):
                y0 = [z1, z2, wr_s, kp0, 1, phi0]
                try:
                    r = least_squares(lambda y: np.r_[(core(ff, y) - norm).real, (core(ff, y) - norm).imag], y0,
                                      bounds=(lo2, hi2), x_scale=[1e-4, 1e-4, 1e-3, 1, .1, .1], max_nfev=2000)
                except Exception:
                    continue
                if best is None or r.cost < best.cost:
                    best = r
    # 3. joint polish with the background (io_model background is centred on (wp + wr)/2)
    y = best.x; xr = to_x(y); fcm = (xr[0] + xr[1]) / 2
    k = 1 + b1 * (fcm - fc)
    def full(f_, q):
        xr_ = to_x(q[:6]); B_, pb_, tau_, b1_ = q[6:]
        fcq = (xr_[0] + xr_[1]) / 2
        return io_model(f_, xr_ + [B_, pb_, tau_, b1_]) if True else None
    q0 = list(y) + [B * abs(k), pb + np.angle(k), tau, b1 / k if k != 0 else b1]
    lo = lo2 + [0, -np.inf, tau - 5, -100]; hi = hi2 + [10, np.inf, tau + 5, 100]
    q0 = [min(max(v, l + 1e-12), h - 1e-12) for v, l, h in zip(q0, lo, hi)]
    r = least_squares(lambda q: np.r_[(full(ff, q) - ss).real, (full(ff, q) - ss).imag], q0, bounds=(lo, hi),
                      x_scale=[1e-4, 1e-4, 1e-3, 1, .1, .1, .1, .1, .01, .1], max_nfev=4000)
    x = np.array(to_x(r.x[:6]) + list(r.x[6:]))
    rms = float(np.sqrt(2 * r.cost / len(ff)))
    return dict(zip(NAMES, map(float, x)), z1=float(r.x[0]), z2=float(r.x[1]), rms=rms,
                rms_rel=rms / float(np.median(np.abs(ss))), window=[float(ff[0]), float(ff[-1])]), x
