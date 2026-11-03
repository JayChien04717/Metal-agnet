"""Dispersive shift of the readout resonators with scqubits (cross-check of zz_analysis.dispersive).

Fluxonium (EJ, EL, Ec from Q3D, half flux) + resonator, H_int = g n (a + a^+), g per charge number
from the Q3D matrix (g = 2 e beta V_zpf). chi = (f_r|1> - f_r|0>), dressed energies via the scqubits lookup.

Usage:  python scq_dispersive.py [f_r1 f_r2]      (GHz, default 6.3 6.5 = design)
"""
import sys

import numpy as np
import scqubits as scq

scq.settings.PROGRESSBAR_DISABLED = True
EJ, EL = 4.0, 0.6
# Q3D ro_w10 (readout pad 10 um, gap 10 / 11.45 um): Ec [GHz], g per n [GHz]
QUBITS = {'FX1': dict(EC=1.3402, g=0.08095, fr=6.3), 'FX2': dict(EC=1.3418, g=0.08051, fr=6.5)}


def chi(EC, g, fr, EJ=EJ, EL=EL, nlev=10, nosc=8):
    fx = scq.Fluxonium(EJ=EJ, EC=EC, EL=EL, flux=0.5, cutoff=110, truncated_dim=nlev)
    osc = scq.Oscillator(E_osc=fr, truncated_dim=nosc)
    hs = scq.HilbertSpace([fx, osc])
    hs.add_interaction(g_strength=g, op1=(fx.n_operator(), fx), op2=(osc.creation_operator() + osc.annihilation_operator(), osc),
                       add_hc=False)
    hs.generate_lookup()
    E = lambda q, r: hs.energy_by_bare_index((q, r))
    f01 = E(1, 0) - E(0, 0)
    fr0, fr1 = E(0, 1) - E(0, 0), E(1, 1) - E(1, 0)
    return dict(f01=f01, fr0=fr0, fr1=fr1, chi=(fr1 - fr0) * 1e3, shift0=(fr0 - fr) * 1e3, shift1=(fr1 - fr) * 1e3)


if __name__ == '__main__':
    frs = [float(a) for a in sys.argv[1:3]] or [None, None]
    for (name, q), fr_arg in zip(QUBITS.items(), frs):
        fr = fr_arg or q['fr']
        r = chi(q['EC'], q['g'], fr)
        print(f"{name}: Ec {q['EC']} g {q['g'] * 1e3:.1f} MHz f_r {fr} GHz -> f01 {r['f01'] * 1e3:.1f} MHz, "
              f"f_r|0> {r['fr0']:.4f}, f_r|1> {r['fr1']:.4f} GHz, shifts {r['shift0']:+.2f} / {r['shift1']:+.2f} MHz, "
              f"chi = {r['chi']:+.2f} MHz")
    print('\nchi [MHz] vs f_r (FX1 parameters):')
    q = QUBITS['FX1']
    for fr in (5.5, 5.8, 6.0, 6.3, 6.6, 6.8, 7.0, 7.3, 7.6, 8.0):
        print(f'  f_r {fr:.1f} GHz: chi = {chi(q["EC"], q["g"], fr)["chi"]:+.2f}')
