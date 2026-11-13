"""Lengths of the 6 readout channels of the Design5 chip (chip6_layout.py): same capacitors as the W1 channel
(C_p hugger 165 um plate, C_j hammer 19 um gap; whole-channel Q3D caps), same tap positions (L_pf 1.2 mm, L_rr 1.3 mm),
filter / resonator lengths solved per channel for omega_r = f_k and omega_p = f_k + 19 MHz (Wallraff Q35 detuning),
omega_p / omega_r from the single-channel circuit on a matched line (Wallraff Table 2.4 definition), modes5.extract.

Usage: python solve_chip6.py   -> results/chip6_lengths.json (read by chip6_layout.py)
"""
import os, sys, json, copy, warnings
warnings.filterwarnings('ignore')
import numpy as np
from scipy.optimize import least_squares
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results'); sys.path.insert(0, HERE)
import theory5 as T, modes5 as M

F_CHIP = [6.0, 6.2, 6.4, 6.6, 6.8, 7.0]       # GHz, 200 MHz spacing, 6-7 GHz (Jay 2026-11-12)
DP = 0.019                                    # omega_p - omega_r (GHz)


def solve(P0, fr):
    def res(x):
        P = copy.deepcopy(P0); P['L_filt'], P['L_res'] = x
        r = M.extract(lambda z: T.branch_Y(z, P), P['Z0'], fr + DP / 2, span=0.8)
        return [(r['wr'] - fr) * 1e3, (r['wp'] - fr - DP) * 1e3] if r else [300, 300]
    # start: scale the 6.70 GHz W1 lengths (line part only) by 6.70 / f
    s = 6.70 / fr
    x0 = [P0['L_filt'] * s - 0.15 * (1 - s), P0['L_res'] * s]
    for _ in range(4):
        e = res(x0)
        x0 = [x0[0] * (1 + e[1] * 1e-3 / fr), x0[1] * (1 + e[0] * 1e-3 / fr)]
    sol = least_squares(res, x0, x_scale=[0.01, 0.01], diff_step=1e-4)
    P = copy.deepcopy(P0); P['L_filt'], P['L_res'] = map(float, sol.x)
    rm = M.extract(lambda z: T.branch_Y(z, P), P['Z0'], fr + DP / 2, span=0.8)
    return P, rm


if __name__ == '__main__':
    P0 = T.params_from_layout()
    print('caps', P0['cap_source'], {k: tuple(round(c, 2) for c in P0[k]) for k in ('cin', 'cp', 'cj')})
    out = dict(cap_source=P0['cap_source'], caps={k: list(P0[k]) for k in ('cin', 'cp', 'cj')}, dp_GHz=DP, channels=[])
    for k, fr in enumerate(F_CHIP, 1):
        P, r = solve(P0, fr)
        out['channels'].append(dict(ch=k, f_r=fr, L_filter=P['L_filt'], L_res=P['L_res'],
                                    wp=r['wp'], wr=r['wr'], J=r['J'], kp_matched=r['kp']))
        print(f'ch{k} f_r {fr:.2f}: L_filter {P["L_filt"]:.4f} mm, L_res {P["L_res"]:.4f} mm -> wp {r["wp"]:.4f} wr {r["wr"]:.4f} '
              f'J {r["J"]:.1f} kp {r["kp"]:.1f} MHz')
    json.dump(out, open(os.path.join(RES, 'chip6_lengths.json'), 'w'), indent=1)
    print('saved', os.path.join(RES, 'chip6_lengths.json'))
