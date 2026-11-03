"""ZZ of the FTF set (FX1 - TC - FX2) from the Q3D charging energies / couplings.

Hamiltonian (GHz), charge operators n_i of the sum / difference nodes (note/capacitance_method.md):

    H = sum_i H_i + J_1c n_1 n_c + J_2c n_2 n_c + J_12 n_1 n_2
    fluxonium  H = 4 Ec n^2 + EL phi^2 / 2 - EJ cos(phi - phi_ext)     (phi_ext = pi, half flux)
    coupler    H = 4 Ec n^2 - EJc cos(phi)                              (EJc tuned by the coupler flux)

Ec and the J's come from the Q3D run in results/latest (J = 4 e^2 C~^-1[i, j] / h per n_i n_j).
ZZ = E_101 - E_100 - E_001 + E_000 of the dressed states (|FX1, TC, FX2>, TC in its ground state).

Usage:
    python zz_analysis.py [EJ] [EL]          # GHz, default 4 0.6
"""
import json
import os
import sys

import numpy as np
from scipy.linalg import eigh, expm
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
SUMMARY = os.path.join(HERE, 'results', 'latest', 'fxarray_fl15_summary.json')
FIG_DIR = os.path.join(os.path.dirname(HERE), 'fig', 'zz')

N_F, N_C = 8, 6          # kept levels: fluxonium, coupler


def fluxonium(Ec, EJ, EL, phi_ext=np.pi, n_keep=N_F, n_ho=160):
    """Energies (GHz, E0 = 0) and the n operator in the eigenbasis (real gauge at half flux)."""
    a = np.diag(np.sqrt(np.arange(1, n_ho)), 1)
    phi = (8 * Ec / EL) ** 0.25 * (a + a.T) / np.sqrt(2)
    n = 1j * (EL / (8 * Ec)) ** 0.25 * (a.T - a) / np.sqrt(2)
    U = expm(1j * phi)
    cos_phi, sin_phi = ((U + U.conj().T) / 2).real, ((U - U.conj().T) / 2j).real
    H = (np.sqrt(8 * Ec * EL) * (a.T @ a + 0.5 * np.eye(n_ho))
         - EJ * (np.cos(phi_ext) * cos_phi + np.sin(phi_ext) * sin_phi))
    E, V = eigh(H)
    V = V[:, :n_keep]
    n_e = V.conj().T @ n @ V
    # phases i^k make n real when the states alternate in parity (half flux)
    d = np.array([1j ** (k % 2) for k in range(n_keep)])
    n_r = d.conj()[:, None] * n_e * d[None, :]
    if np.abs(n_r.imag).max() < 1e-6:
        n_e = n_r.real
    return E[:n_keep] - E[0], n_e


def transmon(Ec, EJ, n_keep=N_C, nmax=30):
    n = np.arange(-nmax, nmax + 1)
    H = np.diag(4 * Ec * n ** 2.0) - EJ / 2 * (np.eye(len(n), k=1) + np.eye(len(n), k=-1))
    E, V = eigh(H)
    V = V[:, :n_keep]
    return E[:n_keep] - E[0], V.T @ np.diag(n) @ V


def coupler_EJ(Ec, f01):
    """EJ of the coupler for a given f01 (GHz)."""
    return brentq(lambda EJ: np.diff(transmon(Ec, EJ, 2)[0])[0] - f01, 0.2, 400)


def dressed_zz(fx1, cp, fx2, J1c, J2c, J12):
    """ZZ (MHz) and the dressed 0 -> 1 frequencies of the fluxoniums (GHz). J in MHz."""
    (E1, n1), (Ec_, nc), (E2, n2) = fx1, cp, fx2
    nf, nct = len(E1), len(Ec_)
    I1, Ic, I2 = np.eye(nf), np.eye(nct), np.eye(nf)
    k = np.kron
    H0 = k(k(np.diag(E1), Ic), I2) + k(k(I1, np.diag(Ec_)), I2) + k(k(I1, Ic), np.diag(E2))
    H = (H0 + J1c * 1e-3 * k(k(n1, nc), I2) + J2c * 1e-3 * k(k(I1, nc), n2)
         + J12 * 1e-3 * k(k(n1, Ic), n2))
    w, V = eigh(H)

    def state(a, c, b):
        idx = (a * nct + c) * nf + b
        return w[np.argmax(np.abs(V[idx, :]) ** 2)]

    e000, e100, e001, e101 = state(0, 0, 0), state(1, 0, 0), state(0, 0, 1), state(1, 0, 1)
    return (e101 - e100 - e001 + e000) * 1e3, e100 - e000, e001 - e000


def load_design():
    with open(SUMMARY) as f:
        r = json.load(f)['results']
    return dict(Ec1=r['Ec_FX1 [GHz]'], Ec2=r['Ec_FX2 [GHz]'], Ecc=r['Ec_TC [GHz]'],
                J1c=r['J_1c [MHz]'], J2c=r['J_2c [MHz]'], J12=r['J_12 [MHz]'])


class Device:
    def __init__(self, EJ=4.0, EL=0.6, **kw):
        d = dict(load_design(), **kw)
        self.d, self.EJ, self.EL = d, EJ, EL
        self.fx1 = fluxonium(d['Ec1'], EJ, EL)
        self.fx2 = fluxonium(d['Ec2'], EJ, EL)

    def coupler(self, fc):
        EJc = coupler_EJ(self.d['Ecc'], fc)
        return transmon(self.d['Ecc'], EJc), EJc

    def zz(self, fc, J1c=None, J2c=None, J12=None):
        cp, EJc = self.coupler(fc)
        d = self.d
        return dressed_zz(self.fx1, cp, self.fx2, d['J1c'] if J1c is None else J1c,
                          d['J2c'] if J2c is None else J2c, d['J12'] if J12 is None else J12)

    def g_eff(self, fc):
        """Matrix-element couplings in MHz: g_12 = |J12 n01 n01|, g_1c = |J1c n01 n01_c|."""
        cp, _ = self.coupler(fc)
        n1, n2, nc = self.fx1[1][0, 1], self.fx2[1][0, 1], cp[1][0, 1]
        d = self.d
        return abs(d['J12'] * n1 * n2), abs(d['J1c'] * n1 * nc), abs(d['J2c'] * n2 * nc)


def main():
    EJ = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
    EL = float(sys.argv[2]) if len(sys.argv) > 2 else 0.6
    dev = Device(EJ, EL)
    d = dev.d
    print(f'design: Ec FX1/FX2/TC = {d["Ec1"]:.3f}/{d["Ec2"]:.3f}/{d["Ecc"]:.3f} GHz, '
          f'J1c/J2c/J12 = {d["J1c"]:.1f}/{d["J2c"]:.1f}/{d["J12"]:.1f} MHz')
    print(f'fluxonium EJ = {EJ}, EL = {EL} GHz, half flux')
    E, n = dev.fx1
    print('FX1 levels [GHz]:', np.round(E[:6], 3), ' f01 = %.3f GHz, f12 = %.3f GHz' % (E[1], E[2] - E[1]))
    print('FX1 |n01| = %.3f, |n12| = %.3f, |n03| = %.3f' % (abs(n[0, 1]), abs(n[1, 2]), abs(n[0, 3])))
    print()
    print(' fc [GHz]   EJc   |  ZZ [kHz]: full   no J12   no coupler(J12 only) | g12  g1c  g2c [MHz] | f01 FX1 / FX2 [MHz]')
    rows = []
    for fc in (3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0):
        z, f1, f2 = dev.zz(fc)
        z_nodirect = dev.zz(fc, J12=0.0)[0]
        z_direct = dev.zz(fc, J1c=0.0, J2c=0.0)[0]
        g12, g1c, g2c = dev.g_eff(fc)
        EJc = dev.coupler(fc)[1]
        rows.append((fc, EJc, z * 1e3, z_nodirect * 1e3, z_direct * 1e3, g12, g1c, g2c, f1 * 1e3, f2 * 1e3))
        print(f' {fc:6.2f}  {EJc:6.1f}  | {z * 1e3:9.1f} {z_nodirect * 1e3:9.1f} {z_direct * 1e3:9.1f}'
              f'      | {g12:5.2f} {g1c:5.1f} {g2c:5.1f} | {f1 * 1e3:6.1f} / {f2 * 1e3:6.1f}')
    os.makedirs(FIG_DIR, exist_ok=True)
    np.savetxt(os.path.join(FIG_DIR, f'zz_vs_fc_EJ{EJ:g}_EL{EL:g}.csv'), np.array(rows), delimiter=',',
               header='fc_GHz,EJc_GHz,ZZ_kHz,ZZ_noJ12_kHz,ZZ_J12only_kHz,g12_MHz,g1c_MHz,g2c_MHz,f01_FX1_MHz,f01_FX2_MHz')


if __name__ == '__main__':
    main()


# ---------------------------------------------------------------------------------------------
# extra analyses: exchange coupling, readout dispersive shift, ZZ map, tolerance
# ---------------------------------------------------------------------------------------------
def trunc(op, nf):
    E, n = op
    return E[:nf], n[:nf, :nf]


def exchange_g(dev, fc):
    """f-f exchange coupling g (MHz): half the |100> <-> |001> splitting at zero detuning."""
    cp, _ = dev.coupler(fc)
    d = dev.d
    fx2 = (dev.fx1[0], dev.fx1[1])                      # identical FX2 -> resonant
    (E1, n1), (Ec_, nc) = dev.fx1, cp
    nf, nct = len(E1), len(Ec_)
    k, I1, Ic = np.kron, np.eye(nf), np.eye(nct)
    H0 = k(k(np.diag(E1), Ic), I1) + k(k(I1, np.diag(Ec_)), I1) + k(k(I1, Ic), np.diag(E1))
    H = (H0 + d['J1c'] * 1e-3 * k(k(n1, nc), I1) + d['J1c'] * -1 * 1e-3 * k(k(I1, nc), n1)
         + d['J12'] * 1e-3 * k(k(n1, Ic), n1))
    w, V = eigh(H)
    pick = lambda idx: np.argsort(-np.abs(V[idx, :]) ** 2)[:2]
    ia, ib = (1 * nct + 0) * nf + 0, (0 * nct + 0) * nf + 1
    pair = sorted(set(pick(ia)) & set(pick(ib)) or set(pick(ia)))[:2]
    return abs(w[pair[1]] - w[pair[0]]) / 2 * 1e3


def dispersive(fx, g_n, f_r, nres=7, nlev=6):
    """Resonator frequency shift (MHz) when FX is in |0> / |1>, and chi = f_r|1> - f_r|0>."""
    E, n = trunc(fx, nlev)
    a = np.diag(np.sqrt(np.arange(1, nres)), 1)
    k = np.kron
    H = k(np.diag(E), np.eye(nres)) + f_r * k(np.eye(nlev), a.T @ a) + g_n * 1e-3 * k(n, a + a.T)
    w, V = eigh(H)
    st = lambda i, m: w[np.argmax(np.abs(V[i * nres + m, :]) ** 2)]
    fr = [(st(i, 1) - st(i, 0)) for i in range(2)]
    return [(x - f_r) * 1e3 for x in fr], (fr[1] - fr[0]) * 1e3


def zz_map(dev, fc=5.0, xmax=200, ymax=700, nx=41, ny=36, nf=7, nc_=5):
    cp, _ = dev.coupler(fc)
    fx1, fx2, cp = trunc(dev.fx1, nf), trunc(dev.fx2, nf), trunc(cp, nc_)
    xs, ys = np.linspace(0, xmax, nx), np.linspace(0, ymax, ny)
    ratio = dev.d['J2c'] / dev.d['J1c']
    Z = np.zeros((ny, nx))
    for iy, y in enumerate(ys):
        for ix, x in enumerate(xs):
            Z[iy, ix] = dressed_zz(fx1, cp, fx2, -y, y * abs(ratio), -x)[0] * 1e3   # kHz
    return xs, ys, Z


def plot_map(dev, fc=5.0, path=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    xs, ys, Z = zz_map(dev, fc)
    d = dev.d
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    im = ax.pcolormesh(xs, ys, np.abs(Z), norm=LogNorm(1e-3, 1e1), cmap='Blues', shading='auto')

    ax.plot(abs(d['J12']), abs(d['J1c']), 'o', mfc='none', mec='red', ms=7)
    ax.annotate('Design1 (%.1f, %.0f)' % (abs(d['J12']), abs(d['J1c'])), (abs(d['J12']), abs(d['J1c'])),
                (abs(d['J12']) + 20, abs(d['J1c']) + 60), color='red', fontsize=9)
    ax.plot(40, 236, 'x', color='gray')
    ax.annotate('paper device (40, 236)', (40, 236), (60, 280), color='gray', fontsize=8)
    ax.set_xlabel(r'$J_{12}/2\pi$ (MHz)  [coefficient of $n_1 n_2$]')
    ax.set_ylabel(r'$J_{1c}/2\pi$ (MHz)  [coefficient of $n_1 n_c$]')
    ax.set_title(f'|ZZ| of FX-TC-FX, EJ={dev.EJ:g}, EL={dev.EL:g}, f_c={fc:g} GHz', fontsize=9)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(r'|ZZ|/2$\pi$ (kHz)')
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150)
    return Z, xs, ys


def extras(dev):
    d = dev.d
    print('\n--- readout dispersive shift (full numerics, g per n from Q3D, f_r design) ---')
    for name, fx, g_n, f_r in (('FX1-R1', dev.fx1, 161.3, 6.3), ('FX2-R2', dev.fx2, 160.2, 6.5)):
        sh, chi = dispersive(fx, g_n, f_r)
        print(f'  {name}: f_r {f_r} GHz, shift |0> {sh[0]:+.2f} MHz, |1> {sh[1]:+.2f} MHz, chi = {chi:+.2f} MHz')
    print('\n--- exchange coupling |100> <-> |001> (MHz), identical FX1 / FX2 ---')
    for fc in (1.5, 2.0, 3.0, 4.0, 5.0, 6.0):
        print(f'  fc {fc:4.1f} GHz  g_exch = {exchange_g(dev, fc):.3f} MHz   (direct J12 n01^2 = '
              f'{abs(d["J12"] * dev.fx1[1][0, 1] ** 2):.3f} MHz)')
    print('\n--- sensitivity (fc = 5 GHz): ZZ [kHz] ---')
    for dEJ in (-0.1, 0, 0.1):
        for dEL in (-0.1, 0, 0.1):
            dv = Device(dev.EJ * (1 + dEJ), dev.EL * (1 + dEL))
            print(f'  EJ {dev.EJ * (1 + dEJ):.2f}  EL {dev.EL * (1 + dEL):.2f}:  ZZ = {dv.zz(5.0)[0] * 1e3:+.2f},'
                  f'  f01 = {dv.zz(5.0)[1] * 1e3:.0f} MHz')
    print('\n--- ZZ-nulling J12 for the design J1c, fc = 5 GHz ---')
    cp, _ = dev.coupler(5.0)
    f = lambda x: dressed_zz(dev.fx1, cp, dev.fx2, d['J1c'], d['J2c'], -x)[0]
    xs = np.linspace(0, 60, 61)
    zs = [f(x) for x in xs]
    for i in range(len(xs) - 1):
        if zs[i] * zs[i + 1] < 0:
            print(f'  ZZ = 0 at |J12| = {brentq(f, xs[i], xs[i + 1]):.2f} MHz (design {abs(d["J12"]):.1f} MHz)')
    print('  ZZ(J12) [kHz] at 0, 10, 20, 30, 40, 60 MHz:', [round(f(x) * 1e3, 2) for x in (0, 10, 20, 30, 40, 60)])
    print(f'  ZZ(J12 = 0) = {zs[0] * 1e3:.2f} kHz, ZZ(J12 = design) = {f(abs(d["J12"])) * 1e3:.2f} kHz')
    print('\n--- Ec sweep (EJ, EL fixed; n01 and transitions) ---')
    for Ec in (0.8, 1.0, 1.2, 1.32, 1.5, 1.8):
        E, n = fluxonium(Ec, dev.EJ, dev.EL)
        print(f'  Ec {Ec:4.2f}: f01 {E[1]:.3f}  f12 {E[2] - E[1]:.3f}  f03 {E[3]:.3f} GHz,  |n01| {abs(n[0, 1]):.3f}')
    os.makedirs(FIG_DIR, exist_ok=True)
    Z, xs, ys = plot_map(dev, 5.0, os.path.join(FIG_DIR, 'zz_map_J12_J1c.png'))
    print('\nfigure:', os.path.join(FIG_DIR, 'zz_map_J12_J1c.png'))


if __name__ == '__main__' and '--extras' in sys.argv:
    sys.argv.remove('--extras')
    _dev = Device(float(sys.argv[1]) if len(sys.argv) > 1 else 4.0, float(sys.argv[2]) if len(sys.argv) > 2 else 0.6)
    extras(_dev)
