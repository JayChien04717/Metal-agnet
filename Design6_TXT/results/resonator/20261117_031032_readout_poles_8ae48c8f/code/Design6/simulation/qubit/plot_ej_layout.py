"""Junction E_J on the layout: GDS zooms of the three SQUIDs with each junction's L_J / E_J / I_c / R_n, plus the
effective E_J(Phi) = E_J,sum sqrt(cos^2(pi Phi) + d^2 sin^2(pi Phi)) and f01(Phi).

Input : results/2026-11-15/ej_table/ej_table.json (written with scqenv, see README section 4: HFSS-retargeted L_J,
        HFSS Ec, R_n from Ambegaokar-Baratoff with Al Delta = 180 ueV), gds/Design6_2026-11-15_sim.gds
Output: fig/Design6_2026-11-15_ej_layout.png
Usage (qmetal): python simulation/plot_ej_layout.py
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import json, os
import numpy as np
import gdstk
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

SIM = _S('Design6/simulation'); DESIGN = _S('Design6')
def result(*path):
    with open(_relocate(os.path.join(SIM, 'results', *path)), encoding='utf-8') as source:
        return json.load(source)


D = result('2026-11-15', 'ej_table', 'ej_table.json')
CZ = result('2026-11-14', 'cz_hfss_v8', 'cz_summary.json')
P = gdstk.read_gds(_S('Design6/gds/Design6_2026-11-15_sim.gds')).top_level()[0].flatten().polygons
# junction sites (um, from the build_design() junction table) and their keys in ej_table.json
SQ = dict(Q1=dict(c=(-450, -205), jj=[(-465, 'Q1a', 'a'), (-435, 'Q1b', 'b')]),
          TC=dict(c=(0, -200), jj=[(-15, 'TC', 'a'), (15, 'TC', 'b')]),
          Q2=dict(c=(450, -205), jj=[(435, 'Q2a', 'a'), (465, 'Q2b', 'b')]))
COL = dict(Q1='#1f77b4', TC='#2ca02c', Q2='#d62728')
METAL, ETCH = '#e6e6e6', '#58606e'


def draw(ax, xl, yl):
    ax.set_facecolor(ETCH)
    for p in P:
        (x0, y0), (x1, y1) = p.bounding_box()
        if p.layer == 1 and x1 > xl[0] and x0 < xl[1] and y1 > yl[0] and y0 < yl[1]:
            ax.add_patch(Polygon(p.points, fc=METAL, ec='none'))
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect('equal')


def main():
    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], hspace=0.28, wspace=0.18)
    J, S = D['junctions'], D['squids']
    params, best = CZ['params'], CZ['best']
    for device, sq in S.items():
        keys = (device + 'a', device + 'b') if device != 'TC' else ('TC',)
        values = [sq['EJsum_GHz'], sq['EC_GHz'], sq['d'], *(J[key]['L_nH'] for key in keys)]
        expected = [params['EJsum_GHz'][device], params['Ec_GHz'][device],
                    params['d'][device], *(params['Lj_nH'][key] for key in keys)]
        np.testing.assert_allclose(values, expected, rtol=1e-8, atol=1e-10,
                                   err_msg=f'{device}: junction table differs from the CZ model')
    cz_flux = {'Q1': best['a_q'], 'TC': CZ['idle']['coupler_flux'] + best['a_c'], 'Q2': 0.0}
    for i, n in enumerate(('Q1', 'TC', 'Q2')):
        ax = fig.add_subplot(gs[0, i]); cx, cy = SQ[n]['c']
        draw(ax, (cx - 42, cx + 42), (cy - 40, cy + 26))
        for x, key, lab in SQ[n]['jj']:
            j = J[key]; y = cy
            ax.plot([x, x], [y - 2.5, y + 2.5], color='#ff7f0e', lw=5, solid_capstyle='butt')
            left = x < cx
            ax.annotate(f"JJ {lab}\nL_J {j['L_nH']:.2f} nH\nE_J {j['EJ_GHz']:.2f} GHz\nI_c {j['Ic_nA']:.1f} nA\n"
                        f"R_n {j['Rn_kohm']:.1f} kΩ",
                        xy=(x, y), xytext=(cx - 40 if left else cx + 40, cy - 37), ha='left' if left else 'right',
                        va='bottom', fontsize=9.5, family='monospace',
                        bbox=dict(boxstyle='round', fc='white', ec=COL[n], alpha=0.93),
                        arrowprops=dict(arrowstyle='->', color='#ff7f0e', lw=1.4))
        s = S[n]; fmax, fmin = s['f01_GHz'][0], s['f01_GHz'][-1]
        ax.set_title(f"{n} SQUID:  E_J,Σ {s['EJsum_GHz']:.2f} GHz,  d = {s['d']:.3f}\n"
                     f"E_J,eff {s['EJsum_GHz']:.2f} → {s['EJeff_GHz'][-1]:.2f} GHz (Φ = 0 → Φ0/2)\n"
                     f"Bare f01 {fmax:.2f} → {fmin:.2f} GHz",
                     fontsize=10.5, color=COL[n])
        ax.set_xlabel('x (µm)')
    # E_J,eff(Phi) and f01(Phi)
    fl = np.array(D['flux'])
    a1 = fig.add_subplot(gs[1, 0]); a2 = fig.add_subplot(gs[1, 1])
    for n in ('Q1', 'TC', 'Q2'):
        s = S[n]; eff = np.array(s['EJeff_GHz']); f = np.array(s['f01_GHz'])
        a1.plot(fl, eff, color=COL[n], lw=2, label=f"{n}: {eff[0]:.2f} → {eff[-1]:.2f} GHz")
        ok = eff / s['EC_GHz'] > 10          # transmon regime; below: dotted (charge regime, not used)
        a2.plot(fl[ok], f[ok], color=COL[n], lw=2, label=n)
        a2.plot(fl[~ok], f[~ok], color=COL[n], lw=1.2, ls=':')
        bias = cz_flux[n]
        phase = np.pi * bias
        ej_cz = params['EJsum_GHz'][n] * np.hypot(np.cos(phase), params['d'][n] * np.sin(phase))
        a1.plot(bias, ej_cz, 'o', color=COL[n])
        a2.plot(bias, np.interp(bias, fl, f), 'o', color=COL[n])
    a1.set_xlabel('Φ / Φ0'); a1.set_ylabel('E_J,eff / h (GHz)'); a1.legend(); a1.grid(alpha=0.3)
    a1.set_title('E_J,eff(Φ) = E_J,Σ √(cos²πΦ + d² sin²πΦ)')
    a2.set_xlabel('Φ / Φ0'); a2.set_ylabel('f01 (GHz)'); a2.legend(); a2.grid(alpha=0.3)
    a2.set_title('Bare f01(Φ), scqubits TunableTransmon (dotted: E_J,eff/E_C < 10)')
    # summary table
    a3 = fig.add_subplot(gs[1, 2]); a3.axis('off')
    rows = [['', 'E_C', 'a : b', 'E_J,Σ', 'E_J,eff\nat Φ0/2', 'f_max', 'f_min']]
    for n in ('Q1', 'Q2', 'TC'):
        s = S[n]; ea, eb = (J[n + 'a'], J[n + 'b']) if n != 'TC' else (J['TC'], J['TC'])
        rows.append([n, f"{s['EC_GHz']:.3f}", f"{ea['EJ_GHz'] / eb['EJ_GHz']:.2f} : 1", f"{s['EJsum_GHz']:.2f}",
                     f"{s['EJeff_GHz'][-1]:.2f}", f"{s['f01_GHz'][0]:.2f}",
                     f"{s['f01_GHz'][-1]:.2f}"])
    t = a3.table(cellText=rows, loc='upper center', cellLoc='center'); t.auto_set_font_size(False); t.set_fontsize(10.5)
    t.scale(1, 2.4)
    a3.text(0.0, 0.49, f"{best['T_ns']:g} ns CZ bias (markers, Φ/Φ0)\n"
            f"Q1 {cz_flux['Q1']:.3f} · TC {cz_flux['TC']:.3f} · Q2 {cz_flux['Q2']:.3f}\n"
            f"Dressed f01 at CZ · Q1 {best['f1_gate_GHz']:.3f} · TC {best['fc_gate_GHz']:.3f} GHz\n"
            'Table values use GHz (E/h and f01).\n'
            'Symmetric TC reaches the charge regime near half flux.\n'
            'R_n uses the low-T Ambegaokar–Baratoff estimate for Al Δ = 180 µeV.',
            fontsize=9.5, va='top', transform=a3.transAxes)
    fig.suptitle(f"Design6 · {best['T_ns']:g} ns CZ Hamiltonian on the simulation GDS", fontsize=14)
    out = _S('Design6/fig/Design6_2026-11-15_ej_layout.png')
    fig.savefig(out, dpi=110, bbox_inches='tight'); print('saved', out)


if __name__ == '__main__':
    main()
