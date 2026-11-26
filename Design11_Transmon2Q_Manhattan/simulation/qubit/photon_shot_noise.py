# -*- coding: utf-8 -*-
r"""Thermal-photon shot-noise dephasing of Q1 / Q2 through their readout resonators, and the resulting T2 estimate.

Environment: scqenv (analysis only).
    & "C:\Users\QEL\miniforge-pypy3\envs\scqenv\python.exe" -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/qubit/photon_shot_noise.py

Model (dispersive regime, resonator with total linewidth kappa and full dispersive shift 2chi; Clerk & Utami,
PRA 75, 042302 (2007); Sears et al., PRB 86, 180504 (2012); Rigetti et al., PRB 86, 100506 (2012)):
    Gamma_phi = (kappa / 2) Re[ sqrt( (1 + i 2chi/kappa)^2 + i 4 (2chi) n_th / kappa ) - 1 ]      (angular units)
    low-n_th limit: Gamma_phi ~= n_th kappa (2chi)^2 / (kappa^2 + (2chi)^2)  -> n_th kappa / 2 at |2chi| = kappa
    n_th = 1 / (exp(h f_r / k_B T_eff) - 1)
    1/T2 = 1/(2 T1) + Gamma_phi (photon shot noise only; flux noise at the sweet spot, TLS and 1/f charge noise not included)
Inputs: |2chi| from the adopted scqubits run (design.json adopted_runs.readout_chi), kappa_r = 3.0 MHz design target
(r15 driven fit 2.93 / 2.97 MHz; r17 not simulated), f_r = 6.35 / 6.15 GHz. T1 values are assumptions (50 / 100 us,
README section 0) plus the Purcell limit; T_eff of the readout mode is unknown and scanned.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general')]

import numpy as np

from run_record import create_run, finish_run

DESIGN = 'Design11_Transmon2Q_Manhattan'
H, KB = 6.62607015e-34, 1.380649e-23
TP = 2 * np.pi
KAPPA_MHZ = 3.0
T_EFF_MK = [30, 40, 50, 60, 70, 80]
N_TH = [1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2]
T1_US = [50, 100]


def gamma_phi(n, kappa_mhz, two_chi_mhz):
    """Shot-noise dephasing rate (1/s) for thermal occupation n."""
    k, x = TP * kappa_mhz * 1e6, TP * abs(two_chi_mhz) * 1e6
    return float(k / 2 * np.real(np.sqrt((1 + 1j * x / k) ** 2 + 4j * x * n / k) - 1))


def n_thermal(f_ghz, t_mk):
    return float(1 / np.expm1(H * f_ghz * 1e9 / (KB * t_mk * 1e-3)))


def main():
    dj = json.loads((ROOT / DESIGN / 'design.json').read_text(encoding='utf-8'))
    chi_run = dj['adopted_runs']['readout_chi']
    two_chi = chi_run['results']['two_chi_MHz']
    t1p_ms = chi_run['results']['T1_purcell_sete_ms']
    f_r = dict(Q1=6.35, Q2=6.15)
    out = dict(inputs=dict(kappa_r_MHz=KAPPA_MHZ, two_chi_MHz=two_chi, f_r_GHz=f_r, readout_chi_run=chi_run['run_id'],
                           T1_assumed_us=T1_US, T1_purcell_sete_ms=t1p_ms), qubits={})
    for q in ('Q1', 'Q2'):
        rows_t, rows_n = [], []
        for t in T_EFF_MK:
            n = n_thermal(f_r[q], t)
            g = gamma_phi(n, KAPPA_MHZ, two_chi[q])
            rows_t.append(dict(T_eff_mK=t, n_th=n, Gamma_phi_per_s=g, T_phi_us=1e6 / g if g > 0 else float('inf'),
                               T2_us={str(t1): 1e6 / (1 / (2 * t1 * 1e-6) + g) for t1 in T1_US}))
        for n in N_TH:
            g = gamma_phi(n, KAPPA_MHZ, two_chi[q])
            rows_n.append(dict(n_th=n, Gamma_phi_per_s=g, T_phi_us=1e6 / g,
                               T2_us={str(t1): 1e6 / (1 / (2 * t1 * 1e-6) + g) for t1 in T1_US}))
        out['qubits'][q] = dict(vs_T_eff=rows_t, vs_n_th=rows_n,
                                low_n_slope_T_phi_us_times_n=1e6 / (TP * KAPPA_MHZ * 1e6 * (two_chi[q] / KAPPA_MHZ) ** 2
                                                                     / (1 + (two_chi[q] / KAPPA_MHZ) ** 2)))
    folder = create_run(DESIGN, 'qubit', 'photon_shot_noise', entrypoint=__file__, inputs=[ROOT / DESIGN / 'design.json'],
                        parameters=out['inputs'])
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    (folder / 'summaries' / 'summary.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    finish_run(folder, status='completed', summary=out)
    print('RUN_FOLDER=' + str(folder))
    for q, v in out['qubits'].items():
        print(q, f"|2chi| {abs(two_chi[q]):.3f} MHz, kappa {KAPPA_MHZ} MHz; T_phi ~ {v['low_n_slope_T_phi_us_times_n']:.3g} us / n_th")
        for r in v['vs_T_eff']:
            print(f"   T_eff {r['T_eff_mK']:3d} mK  n_th {r['n_th']:.2e}  T_phi {r['T_phi_us']:9.1f} us  "
                  f"T2 (T1=50) {r['T2_us']['50']:6.1f} us  T2 (T1=100) {r['T2_us']['100']:6.1f} us")


if __name__ == '__main__':
    main()
