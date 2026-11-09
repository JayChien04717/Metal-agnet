"""kappa_p of every channel INSIDE the chip: the branch of channel k sees at its node n0_k the environment
Y_env = Y_left (toward P1: feed line, other channels, C_in, P1) + Y_right (toward P2: feed line, other channels, P2).
omega_p, omega_r, J, kappa_p from modes5.extract with that Y_env. Compared with
  * the chip without C_in (C_in replaced by a through line),
  * the matched-line value (Wallraff Table 2.4 definition, Y_env = 2/Z0).
Also: C_p needed per channel for kappa_p = KAPPA_TARGET in the chip (kappa_p ~ C_p^2, iterated).

Physics: C_in (~580 Ohm at 6.7 GHz) is almost an open end for the feed line. Looking back toward it from n0 at
electrical distance theta, Y_left ~ j Y0 tan(theta + theta_Cin), so Re Z_env = Z0 / (1 + tan^2) oscillates with
period lambda/2: ~Z0 (kappa ~2x matched) near C_in and at lambda/2, -> 0 (kappa -> 0, voltage node) at lambda/4.

Usage: python chip6_kappa_env.py   -> results/chip6_kappa_env.json / .png
"""
import os, sys, json, copy, warnings
warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results')
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, HERE)
import theory5 as T, modes5 as M
import chip6_layout as C6, ftf_layout as f3

KAPPA_TARGET = 40.0     # MHz, in-chip


def abcd_chain(f, P0, items):
    """items: list of ('tl', l_mm) | ('Y', P) | ('cin', None); returns (A, B, C, D) arrays (f may be complex)."""
    f = np.atleast_1d(np.asarray(f, complex)); Z0, v = P0['Z0'], P0['v']; Y0 = 1 / Z0
    w = 2 * np.pi * f * 1e9; beta = w / v
    A = np.ones_like(w); B = np.zeros_like(w); C = np.zeros_like(w); D = np.ones_like(w)
    def mul(m):
        nonlocal A, B, C, D
        a, b, c, d = m
        A, B, C, D = A * a + B * c, A * b + B * d, C * a + D * c, C * b + D * d
    for kind, x in items:
        if kind == 'tl':
            cs, sn = np.cos(beta * x * 1e-3), np.sin(beta * x * 1e-3); mul((cs, 1j * Z0 * sn, 1j * Y0 * sn, cs))
        elif kind == 'Y':
            mul((1, 0, T.branch_Y(f, x), 1))
        elif kind == 'cin':
            cis, cga, cgb = P0['cin']; jC = lambda c: 1j * w * c * 1e-15
            mul((1, 0, jC(cga), 1)); mul((1, 1 / jC(cis), 0, 1)); mul((1, 0, jC(cgb), 1))
    return A, B, C, D


def chip_items(P0, chans, segs, with_cin=True):
    """Feed-line items from P1 to P2; chans: list of P in feed-line order."""
    items = [('tl', P0['l_feed_in'])] + ([('cin', None)] if with_cin else [('tl', 0.0)])
    for k, P in enumerate(chans):
        items += [('tl', segs[k]), ('node', k)]
    items += [('tl', segs[-1])]
    return items


def env_Y(k, P0, chans, segs, with_cin):
    """Admittance at n0_k without branch k: left part (P1 side, terminated Z0) + right part (P2 side, terminated Z0)."""
    items = chip_items(P0, chans, segs, with_cin)
    i = items.index(('node', k))
    left = [('Y', chans[j]) if it[0] == 'node' else it for it in items[:i] for j in [it[1] if it[0] == 'node' else 0]]
    right = [('Y', chans[j]) if it[0] == 'node' else it for it in items[i + 1:] for j in [it[1] if it[0] == 'node' else 0]]
    Z0 = P0['Z0']
    def Y(f):
        A, B, C, D = abcd_chain(f, P0, left); Zl = (D * Z0 + B) / (C * Z0 + A)     # seen from n0 back to P1
        A, B, C, D = abcd_chain(f, P0, right); Zr = (A * Z0 + B) / (C * Z0 + D)    # seen from n0 to P2
        out = 1 / Zl + 1 / Zr
        return out if np.ndim(f) else out[0]
    return Y


if __name__ == '__main__':
    d = C6.build_chip(); P0 = T.params_from_layout()
    L = {c['ch']: c for c in json.load(open(C6.LEN_JSON))['channels']}
    pl = lambda n: f3.path_length(d, n); nl = len(C6.ORDER)
    P0['l_feed_in'] = pl('feed_in') + 0.04
    segs = [pl('fl_0') + 0.04 + 0.15] + [pl(f'fl_{k}') + 0.3 for k in range(1, nl)] + [pl(f'fl_{nl}') + 0.15]
    chans = []
    for ch, _, _ in C6.ORDER:
        P = copy.deepcopy(P0); P['L_filt'], P['L_res'] = L[ch]['L_filter'], L[ch]['L_res']; chans.append(P)
    dist = np.cumsum([0.0] + segs[:-1])[1:] + 0.0          # feed-line distance C_in -> n0_k (mm)
    rows = []
    for k, (ch, _, _) in enumerate(C6.ORDER):
        P = chans[k]; fc = L[ch]['f_r'] + 0.0095
        Yb = lambda z, P=P: T.branch_Y(z, P)
        r_cin = M.extract(Yb, P['Z0'], fc, span=0.5, Yenv=env_Y(k, P0, chans, segs, True))
        r_no = M.extract(Yb, P['Z0'], fc, span=0.5, Yenv=env_Y(k, P0, chans, segs, False))
        r_m = M.extract(Yb, P['Z0'], fc, span=0.5)
        lam4 = P0['v'] / (4 * L[ch]['f_r'] * 1e9) * 1e3
        # C_p for kappa = target in the chip (with C_in): kappa ~ C_p^2, two fixed-point steps
        Pc = copy.deepcopy(P); cp = P['cp'][0]
        for _ in range(3):
            rr = M.extract(lambda z: T.branch_Y(z, Pc), P['Z0'], fc, span=0.5, Yenv=env_Y(k, P0, [Pc if j == k else chans[j] for j in range(nl)], segs, True))
            if rr is None or rr['kp'] < 0.5: cp = None; break
            cp = cp * np.sqrt(KAPPA_TARGET / rr['kp']); Pc['cp'] = (cp,) + tuple(P['cp'][1:])
        rows.append(dict(ch=ch, f_r=L[ch]['f_r'], dist_from_Cin_mm=float(dist[k]), dist_over_lambda4=float(dist[k] / lam4),
                         kp_with_Cin=r_cin and r_cin['kp'], kp_no_Cin=r_no and r_no['kp'], kp_matched=r_m and r_m['kp'],
                         J_with_Cin=r_cin and r_cin['J'], Cp_for_target_fF=cp))
        print(f'ch{ch} f {L[ch]["f_r"]:.1f}  d(C_in->n0) {dist[k]:.2f} mm = {dist[k] / lam4:.2f} lambda/4   '
              f'kappa_p: with C_in {rows[-1]["kp_with_Cin"]:.1f}, no C_in {rows[-1]["kp_no_Cin"]:.1f}, matched {rows[-1]["kp_matched"]:.1f} MHz; '
              f'C_p for {KAPPA_TARGET:.0f} MHz in chip: {cp if cp is None else round(cp, 1)} fF (now {P["cp"][0]:.1f})')
    json.dump(rows, open(os.path.join(RES, 'chip6_kappa_env.json'), 'w'), indent=1, default=float)
    # kappa vs distance for one channel moved along the line (ch3 branch), to show the lambda/2 period
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    ds = np.linspace(0.2, 12, 80); kk = []
    k3 = [o[0] for o in C6.ORDER].index(3); P = chans[k3]; fc3 = L[3]['f_r'] + 0.0095   # the 6.7 GHz channel
    for dd in ds:
        Y = env_Y(0, P0, [P], [dd, 2.0], True)
        r = M.extract(lambda z: T.branch_Y(z, P), P['Z0'], fc3, span=0.5, Yenv=Y)
        kk.append(r['kp'] if r else np.nan)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(ds, kk, 'o-', ms=3, label='single 6.7 GHz channel at distance d from C_in')
    ax.axhline(rows[k3]['kp_matched'], color='k', ls=':', label='6.7 GHz, matched line (no C_in)')
    lam4 = P0['v'] / (4 * 6.7e9) * 1e3
    for m in (1, 2):
        ax.axvline(m * lam4, color='gray', lw=.6, ls='--')
    ax.text(lam4, ax.get_ylim()[1] * 0.95, ' λ/4', fontsize=8); ax.text(2 * lam4, ax.get_ylim()[1] * 0.95, ' λ/2', fontsize=8)
    for r in rows:
        ax.plot(r['dist_from_Cin_mm'], r['kp_with_Cin'], 's', ms=7, label=f'ch{r["ch"]} {r["f_r"]:.1f} GHz (in chip)')
    ax.set_xlabel('feed-line distance C_in → n0 (mm)'); ax.set_ylabel('κ_p / 2π (MHz)'); ax.grid(alpha=.3); ax.legend(fontsize=7, ncol=2)
    ax.set_title('Purcell-filter κ_p vs distance from C_in (circuit model)')
    fig.tight_layout(); fig.savefig(os.path.join(RES, 'chip6_kappa_env.png'), dpi=120)
    print('saved', os.path.join(RES, 'chip6_kappa_env.png'))
