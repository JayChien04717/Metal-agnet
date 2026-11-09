"""Analytic (ABCD) model of the Design3 + Purcell filter chain, to debug the HFSS sweeps.

    C_cap - TL - [tee1 hanger] - TL - [tee2] - TL - [tee3] - TL - C_cap     (ports 50.9 Ohm)

C_cap = Cc of the interdigital caps; each hanger = resonator (open-open line of length
RESONATORS[k]['length'], open at the far end) in series with its tee coupling cap Cc_k.
The extra electrical length delta of the filter (caps' metal, taper) is calibrated so the
filter alone (no tees) peaks at the HFSS value of s21_v3 (6.73 GHz).

Usage: python abcd_model.py
"""
import os
import sys
import warnings

warnings.filterwarnings('ignore')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
import design3_purcell as dp  # noqa: E402
import purcell_layout as pl  # noqa: E402

Z0, V = pl.CPW_Z0, pl.CPW_V
F = np.linspace(6.2, 7.2, 20001) * 1e9
W = 2 * np.pi * F


def tl(l_mm):
    th = W * l_mm * 1e-3 / V
    c, s = np.cos(th), 1j * np.sin(th)
    return np.stack([np.stack([c, Z0 * s], -1), np.stack([s / Z0, c], -1)], -2)


def series_c(C_fF):
    n = len(W)
    M = np.zeros((n, 2, 2), complex)
    M[:, 0, 0] = M[:, 1, 1] = 1
    M[:, 0, 1] = 1 / (1j * W * C_fF * 1e-15)
    return M


def shunt(Y):
    n = len(W)
    M = np.zeros((n, 2, 2), complex)
    M[:, 0, 0] = M[:, 1, 1] = 1
    M[:, 1, 0] = Y
    return M


def hanger(Cc_fF, l_mm):
    """Resonator (open-open line, length l) coupled by series Cc: shunt admittance."""
    th = W * l_mm * 1e-3 / V
    z_in = -1j * Z0 / np.tan(th)                       # open far end
    return shunt(1 / (1 / (1j * W * Cc_fF * 1e-15) + z_in))


def s21(chain):
    M = chain[0]
    for m in chain[1:]:
        M = M @ m
    A, B, C, D = M[:, 0, 0], M[:, 0, 1], M[:, 1, 0], M[:, 1, 1]
    return 2 / (A + B / Z0 + C * Z0 + D)


def db(x):
    return 20 * np.log10(np.abs(x))


def peak(s):
    i = np.argmax(db(s))
    return F[i] / 1e9, db(s)[i]


if __name__ == '__main__':
    d = dp.build_design()
    names = list(dp.RESONATORS)
    cap = pl.CAP_FF
    # calibrate delta on the standalone filter (pin path PIN_LENGTH, HFSS f0 = 6.73 GHz)
    best = None
    for delta in np.linspace(0, 3, 601):
        f0, _ = peak(s21([series_c(cap), tl(pl.PIN_LENGTH + delta), series_c(cap)]))
        if best is None or abs(f0 - 6.73) < best[0]:
            best = (abs(f0 - 6.73), delta)
    delta = best[1]
    print(f'standalone filter: C = {cap:.1f} fF, calibrated extra length delta = {delta:.3f} mm -> f0 {peak(s21([series_c(cap), tl(pl.PIN_LENGTH + delta), series_c(cap)]))[0]:.3f} GHz')

    # node positions along the filter line, from cap_in.b
    mlen = pl.path_length(d, 'filter_in')
    xs = {n: d.components[f'tee_{n}'].options.pos_x for n in names}
    xs = {n: float(str(v).replace('mm', '')) for n, v in xs.items()}
    cl = {n: dp.RESONATORS[n]['coupling_length'] for n in names}
    x1_start = xs[names[0]] - cl[names[0]]            # start of tee1 prime line
    nodes = []
    for n in names:                                    # node at the centre of the tee prime line
        nodes.append(mlen + (xs[n] - x1_start))
    chain = [series_c(cap), tl(delta / 2 + nodes[0])]
    for a, b in zip(nodes[:-1], nodes[1:]):
        chain += [hanger(dp.RESONATORS[names[nodes.index(a)]]['Cc'], dp.RESONATORS[names[nodes.index(a)]]['length']),
                  tl(b - a)]
    last = names[-1]
    chain += [hanger(dp.RESONATORS[last]['Cc'], dp.RESONATORS[last]['length']),
              tl(pl.PIN_LENGTH - nodes[-1] + delta / 2), series_c(cap)]
    # fix: the first hanger of the loop above is tee1; rebuild explicitly to be safe
    chain = [series_c(cap), tl(delta / 2 + nodes[0])]
    for i, n in enumerate(names):
        chain.append(hanger(dp.RESONATORS[n]['Cc'], dp.RESONATORS[n]['length']))
        chain.append(tl((nodes[i + 1] - nodes[i]) if i < len(names) - 1 else pl.PIN_LENGTH - nodes[-1] + delta / 2))
    chain.append(series_c(cap))
    S = s21(chain)
    f0, p0 = peak(S)
    print(f'with 3 hangers: broad peak/max |S21| {f0:.4f} GHz ({p0:.2f} dB); node positions {np.round(nodes, 3)} mm')
    s = db(S)
    from scipy.signal import find_peaks
    idx, pr = find_peaks(-s, prominence=0.3)
    for i, p in zip(idx, pr['prominences']):
        print(f'  dip {F[i] / 1e9:.4f} GHz  {s[i]:.1f} dB  (prominence {p:.1f} dB)')
    sm = np.convolve(s, np.ones(401) / 401, mode='same')
    print(f'smoothed-peak (50 MHz window) at {F[np.argmax(sm[2000:-2000]) + 2000] / 1e9:.3f} GHz')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(F / 1e9, s, label='ABCD model S21')
    ax.plot(F / 1e9, db(s21([series_c(cap), tl(pl.PIN_LENGTH + delta), series_c(cap)])), '--', label='filter alone')
    for r in dp.RESONATORS.values():
        ax.axvline(r['f'], color='gray', ls=':')
    ax.set_xlabel('f (GHz)'); ax.set_ylabel('dB'); ax.set_ylim(-40, 2); ax.legend(); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, 'results', 'abcd_model.png'), dpi=120)
