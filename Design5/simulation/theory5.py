"""Design5 circuit (ABCD / admittance) model (copied from Design6/simulation/theory_abcd.py).

    P1 - TL - [C_in pi] - TL - node n0 -- TL - P2
                           |  n0 shunt branch: TL - [C_p pi] - filter TL (open at C_p, shorted at the far end)
                                                              | tap at D_CJ from C_p: [C_j pi] - resonator TL (open end at d_r, short at L_res)

Capacitances (fF), source selected by env CAP_SOURCE:
  'full' (default when results/q3d5_full_caps.json exists): series C_in, C_p, C_j from the Q3D of the WHOLE channel
         (q3d_full5.py, all capacitors and lines together); the pi shunts (electrode-to-ground) stay from the part runs,
         because in the full run the net-to-ground C also contains the line capacitance that the TL model already has.
         Extra mutuals of the full run (e.g. filter <-> resonator lines) are listed in the json, not used by default.
  'part': everything from the single-component Q3D runs (results/q3d5_part_<name>_cmat_fF.csv).
Lines: Z0 and v of the CPW 10/6 on Si, lossless.
"""
import os, sys, json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results')
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, HERE)
FULL_JSON = os.path.join(RES, 'q3d5_full_caps.json')


def cap_source():
    s = os.environ.get('CAP_SOURCE', '')
    return s if s else ('full' if os.path.exists(FULL_JSON) else 'part')


def q3d_caps(name):
    """(series C, ground C of electrode 1, ground C of electrode 2) in fF; electrode order of the part-run matrix."""
    cm = pd.read_csv(os.path.join(RES, f'q3d5_part_{name}_cmat_fF.csv'), index_col=0)
    sig = [c for c in cm.index if 'ground' not in c.lower()]
    ser = -cm.loc[sig[0], sig[1]]
    g1, g2 = cm.loc[sig[0], sig[0]] - ser, cm.loc[sig[1], sig[1]] - ser
    if cap_source() == 'full':
        ser = json.load(open(FULL_JSON))['series'][name]
    return float(ser), float(g1), float(g2), sig


def branch_Y(f_ghz, P):
    """Admittance of the shunt branch (stub - C_p - filter - C_j - resonator) seen at node n0 (f may be complex)."""
    Z0, v = P['Z0'], P['v']
    w = 2 * np.pi * np.asarray(f_ghz) * 1e9
    beta = w / v; Y0 = 1 / Z0
    t = lambda l: np.tan(beta * l * 1e-3)
    Yopen = lambda l: 1j * Y0 * t(l)
    Yshort = lambda l: -1j * Y0 / t(l)
    cjs, cjg1, cjg2 = P['cj']
    Yres = 1j * w * cjg2 * 1e-15 + Yopen(P['d_r']) + Yshort(P['L_res'] - P['d_r'])
    Yser = (1j * w * cjs * 1e-15 * Yres) / (1j * w * cjs * 1e-15 + Yres)
    Ytap = 1j * w * cjg1 * 1e-15 + Yser + Yshort(P['L_filt'] - P['D_cj'])
    td = t(P['D_cj'])
    Yf = Y0 * (Ytap + 1j * Y0 * td) / (Y0 + 1j * Ytap * td)
    cps, cpg_a, cpg_b = P['cp']                       # (series, ground on the n0 side, ground on the filter side)
    Y1 = 1j * w * cpg_b * 1e-15 + Yf
    Y2 = (1j * w * cps * 1e-15 * Y1) / (1j * w * cps * 1e-15 + Y1)
    Y3 = 1j * w * cpg_a * 1e-15 + Y2
    ts = t(P['l_stub'])
    return Y0 * (Y3 + 1j * Y0 * ts) / (Y0 + 1j * Y3 * ts)


def n0_env_Y(f_ghz, P):
    """Admittance the branch sees at n0 inside the channel: toward P1 (line l_feed_c, C_in pi, line l_feed_in, 50 Ohm
    port) in parallel with toward P2 (matched line). f may be complex. A matched line on both sides would be 2/Z0."""
    Z0, v = P['Z0'], P['v']; Y0 = 1 / Z0
    w = 2 * np.pi * np.asarray(f_ghz) * 1e9; beta = w / v
    line = lambda l, YL: Y0 * (YL + 1j * Y0 * np.tan(beta * l * 1e-3)) / (Y0 + 1j * YL * np.tan(beta * l * 1e-3))
    jC = lambda c: 1j * w * c * 1e-15
    cis, cig_a, cig_b = P['cin']
    Y = line(P['l_feed_in'], Y0 + 0 * w) + jC(cig_a)
    Y = jC(cis) * Y / (jC(cis) + Y) + jC(cig_b)
    return line(P['l_feed_c'], Y) + Y0


def s21_theory(f_ghz, P):
    """Full channel: complex (S11, S21)."""
    Z0, v = P['Z0'], P['v']
    w = 2 * np.pi * np.asarray(f_ghz, float) * 1e9; beta = w / v; Y0 = 1 / Z0
    tl = lambda l: np.array([[np.cos(beta * l * 1e-3), 1j * Z0 * np.sin(beta * l * 1e-3)],
                             [1j * Y0 * np.sin(beta * l * 1e-3), np.cos(beta * l * 1e-3)]])
    mul = lambda A, B: np.einsum('ijn,jkn->ikn', A, B)
    one = np.zeros((2, 2, len(w)), complex); one[0, 0] = 1; one[1, 1] = 1

    def shunt(Y):
        M = one.copy(); M[1, 0] = Y; return M

    def series(Z):
        M = one.copy(); M[0, 1] = Z; return M

    def pi(cser, cga, cgb):
        return mul(mul(shunt(1j * w * cga * 1e-15), series(1 / (1j * w * cser * 1e-15))), shunt(1j * w * cgb * 1e-15))

    cis, cig_a, cig_b = P['cin']
    M = mul(tl(P['l_feed_in']), pi(cis, cig_a, cig_b))
    M = mul(M, tl(P['l_feed_c'])); M = mul(M, shunt(branch_Y(f_ghz, P))); M = mul(M, tl(P['l_out']))
    A, B, C, D = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    den = A + B / Z0 + C * Z0 + D
    return (A + B / Z0 - C * Z0 - D) / den, 2 / den


def params_from_layout():
    import design5_layout as d5, ftf_layout as f3
    d = d5.build_design()
    P = dict(Z0=d5.Z0, v=d5.V, cap_source=cap_source())
    for n in ('cin', 'cp', 'cj'):
        s, g1, g2, sig = q3d_caps(n)
        P[n] = (s, g1, g2)
    # EnclosingCap: frame net on pin a, plate (lead_b) on pin b; part-run net order = [frame, lead_b].
    # T type: pin a toward the line (n0 / P1 side); C type: pin a toward the filter -> swap so the tuple is
    # (series, ground on the line / n0 side, ground on the filter / node side)
    P['cin'] = P['cin'] if d5.CIN_TYPE == 'T' else (P['cin'][0], P['cin'][2], P['cin'][1])
    P['cp'] = P['cp'] if d5.CP_TYPE == 'T' else (P['cp'][0], P['cp'][2], P['cp'][1])
    pl = lambda n: f3.path_length(d, n)
    P['l_feed_in'] = pl('feed_in') + 0.04
    P['l_feed_c'] = pl('feed_c') + 0.08 + 0.15
    P['l_out'] = pl('out_line') + 0.15
    P['l_stub'] = 0.1 + pl('stub_p') + 0.04
    P['L_filt'] = d5.L_FILTER; P['D_cj'] = d5.D_CJ
    P['L_res'] = d5.L_RES_TOTAL; P['d_r'] = d5.VERT + d5.OPEN_DY - d5.D_CJ
    return P
