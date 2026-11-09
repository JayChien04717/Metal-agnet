"""Chip-level circuit model of the Design5 6-channel chip (chip6_layout.py), like Wallraff thesis eq. 2.135:
P1 - feed_in - [C_in pi] - feed line ... shunt Y_branch(ch) at every n0 ... feed line - P2, all lengths from the
drawn layout, each channel = theory5.branch_Y with its own filter / resonator length (same capacitors).

Usage: python chip6_theory.py   -> results/chip6_theory.png, results/chip6_theory.json
"""
import os, sys, json, copy, warnings
warnings.filterwarnings('ignore')
import numpy as np
from scipy.signal import find_peaks
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results')
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, HERE)
import theory5 as T, modes5 as M
import chip6_layout as C6, ftf_layout as f3


def chip_s21(f, P0, chans, segs, s11=False):
    """chans: list of P dicts in feed-line order; segs: feed-line lengths (mm) before ch1, between, after ch6.
    s11=True returns (S11, S21)."""
    Z0, v = P0['Z0'], P0['v']
    w = 2 * np.pi * np.asarray(f) * 1e9; beta = w / v; Y0 = 1 / Z0
    one = np.zeros((2, 2, len(w)), complex); one[0, 0] = one[1, 1] = 1
    mul = lambda A, B: np.einsum('ijn,jkn->ikn', A, B)
    def tl(l):
        M_ = one.copy(); c, s = np.cos(beta * l * 1e-3), np.sin(beta * l * 1e-3)
        M_[0, 0] = M_[1, 1] = c; M_[0, 1] = 1j * Z0 * s; M_[1, 0] = 1j * Y0 * s; return M_
    def shunt(Y):
        M_ = one.copy(); M_[1, 0] = Y; return M_
    def series(Z):
        M_ = one.copy(); M_[0, 1] = Z; return M_
    cis, cga, cgb = P0['cin']
    Mx = mul(tl(P0['l_feed_in']), shunt(1j * w * cga * 1e-15))
    Mx = mul(mul(Mx, series(1 / (1j * w * cis * 1e-15))), shunt(1j * w * cgb * 1e-15))
    for k, P in enumerate(chans):
        Mx = mul(mul(Mx, tl(segs[k])), shunt(T.branch_Y(f, P)))
    Mx = mul(Mx, tl(segs[-1]))
    A, B, Cc, D = Mx[0, 0], Mx[0, 1], Mx[1, 0], Mx[1, 1]
    den = A + B / Z0 + Cc * Z0 + D
    if s11:
        return (A + B / Z0 - Cc * Z0 - D) / den, 2 / den
    return 2 / den


if __name__ == '__main__':
    d = C6.build_chip()
    P0 = T.params_from_layout()                                    # caps, Z0, v, stub, taps of the W1 channel
    L = {c['ch']: c for c in json.load(open(C6.LEN_JSON))['channels']}
    pl = lambda n: f3.path_length(d, n)
    P0['l_feed_in'] = pl('feed_in') + 0.04
    segs = [pl('fl_0') + 0.04 + 0.15] + [pl(f'fl_{k}') + 0.3 for k in range(1, 6)] + [pl('fl_6') + 0.15]
    chans = []
    for ch, _, _ in C6.ORDER:
        P = copy.deepcopy(P0); P['L_filt'], P['L_res'] = L[ch]['L_filter'], L[ch]['L_res']; chans.append(P)
    f = np.linspace(6.1, 7.5, 28001)
    s21 = chip_s21(f, P0, chans, segs)
    a = 20 * np.log10(np.abs(s21)); idx, _ = find_peaks(-a, prominence=3, distance=40)
    dips = [round(float(f[k]), 4) for k in idx]
    print('feed-line segments (mm):', [round(x, 3) for x in segs])
    print('S21 dips (GHz):', dips)
    out = dict(segments_mm=segs, order=[o[0] for o in C6.ORDER], dips=dips, channels=[])
    for (ch, _, _), P in zip(C6.ORDER, chans):
        r = M.extract(lambda z: T.branch_Y(z, P), P['Z0'], L[ch]['f_r'] + 0.0095, span=0.5)
        out['channels'].append(dict(ch=ch, f_r=L[ch]['f_r'], L_filter=P['L_filt'], L_res=P['L_res'], **{k: r[k] for k in ('wp', 'wr', 'J', 'kp')}))
        print(f'ch{ch}: wr {r["wr"]:.4f} wp {r["wp"]:.4f} GHz, J {r["J"]:.1f}, kp (matched) {r["kp"]:.1f} MHz')
    json.dump(out, open(os.path.join(RES, 'chip6_theory.json'), 'w'), indent=1)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(12, 8))
    ax[0].plot(f, np.abs(s21) / np.abs(s21).max(), color='#5b5ea6')
    for c in out['channels']:
        ax[0].text(c['f_r'], 1.03, f'ch{c["ch"]}\n{c["f_r"]:.1f}', ha='center', fontsize=8)
    ax[0].set_ylim(0, 1.12); ax[0].set_ylabel('|S21| (normalised)'); ax[0].grid(alpha=.3)
    ax[0].set_title('Design5 chip: 6 readout + Purcell channels on one feed line (circuit model, Q3D C, drawn lengths)')
    ax[1].plot(f, a, color='#5b5ea6'); ax[1].set_ylabel('S21 (dB)'); ax[1].set_xlabel('f (GHz)'); ax[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(RES, 'chip6_theory.png'), dpi=120)
    print('saved', os.path.join(RES, 'chip6_theory.png'))
