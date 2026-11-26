# -*- coding: utf-8 -*-
r"""Readout option B: g target for |2chi| = kappa_r (scqubits exact transmon + oscillator), n_crit, Purcell T1 and
single-shot assignment estimate; optionally the |2chi| realised by a Q3D core run's g_qr.

Environment: scqenv (quantum analysis).
    & "C:\Users\QEL\miniforge-pypy3\envs\scqenv\python.exe" -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/qubit/readout_optB_chi.py [--core <results/qubit/<q3d_core run>>]

Conventions:
  * g is the transmon 0-1 coupling (repo q3d_core.py: g = (C_qr / 2 sqrt(C_sigma C_r)) sqrt(f_q f_r)); scqubits coupling
    g * n / <0|n|1> so that the 0-1 matrix element equals g.
  * |2chi| = (E11 - E10) - (E01 - E00) of the dressed spectrum (ncut 30, 6 x 6 levels).
  * Purcell T1, two estimates: (a) Lorentzian filter tail gamma = kappa_r(Delta_qp) (g / Delta)^2 with
    kappa_r(Delta) = 4 J^2 kappa_p / (kappa_p^2 + 4 Delta^2); (b) Sete et al. PRA 92, 012325 (2015) band-pass formula
    gamma = kappa_r (g/Delta)^2 (w_q / w_r) (w_r / (2 Q_F Delta)), Q_F = f_p / kappa_p. (b) is the conservative bound.
  * Assignment: homodyne, optimal phase, ring-up included, SNR^2 = 2 eta kappa int |alpha_e - alpha_g|^2 dt,
    eps_sep = erfc(SNR / 2 sqrt 2) / 2, eps_T1 = tau / (4 T1); minimised over tau.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general')]

import numpy as np
import scqubits as scq
from scipy.optimize import brentq
from scipy.special import erfc

from run_record import create_run, finish_run

scq.settings.T1_DEFAULT_WARNING = False
DESIGN = 'Design11_Transmon2Q_Manhattan'
TP = 2 * np.pi
KAPPA_R, KAPPA_P = 3.0, 20.0
# E_J sum (design input, user-confirmed CZ target set) and Q3D E_C of the adopted r8 core (E_C held fixed)
QUBITS = dict(Q1=dict(EJ=14.10, EC=0.2198, fr=6.35), Q2=dict(EJ=11.93, EC=0.2197, fr=6.15))


def spectrum(EJ, EC, fr, g_mhz):
    tm = scq.Transmon(EJ=EJ, EC=EC, ng=0, ncut=30, truncated_dim=6)
    n01 = abs(tm.matrixelement_table('n_operator', evals_count=2)[0, 1])
    osc = scq.Oscillator(E_osc=fr, truncated_dim=6)
    hs = scq.HilbertSpace([tm, osc])
    hs.add_interaction(g=g_mhz * 1e-3 / n01, op1=tm.n_operator, op2=osc.creation_operator, add_hc=True)
    hs.generate_lookup()
    E = lambda q, n: hs.energy_by_bare_index((q, n))
    return abs(((E(1, 1) - E(1, 0)) - (E(0, 1) - E(0, 0))) * 1e3), float(E(1, 0) - E(0, 0))


def assignment(kr, two_chi, nbar, eta, T1_us):
    k, chi = TP * kr * 1e6, TP * two_chi / 2 * 1e6
    t = np.linspace(0, 3e-6, 30001)
    dt = t[1] - t[0]
    den = lambda s: k / 2 + 1j * s * chi
    eps = np.sqrt(nbar) * abs(den(1))
    a = lambda s: -1j * eps / den(s) * (1 - np.exp(-den(s) * t))
    snr2 = 2 * eta * k * np.cumsum(abs(a(1) - a(-1)) ** 2) * dt
    best = None
    for tau in np.arange(20, 3000, 2):
        snr = np.sqrt(snr2[int(tau * 1e-9 / dt)])
        e = 0.5 * erfc(snr / (2 * np.sqrt(2))) + tau / (4 * T1_us * 1e3)
        if best is None or e < best[0]:
            best = (float(e), float(tau), float(snr))
    return dict(F=1 - best[0], tau_ns=best[1], SNR=best[2])


def channel(q, p, g):
    two_chi, fq = spectrum(p['EJ'], p['EC'], p['fr'], g)
    D = (p['fr'] - fq) * 1e3
    J = np.sqrt(KAPPA_R * KAPPA_P) / 2
    kr_q = 4 * J ** 2 * KAPPA_P / (KAPPA_P ** 2 + 4 * D ** 2)
    T1a = 1 / (TP * kr_q * 1e6 * (g / D) ** 2)
    QF = p['fr'] * 1e3 / KAPPA_P
    T1b = 1 / (TP * KAPPA_R * 1e6 * (g / D) ** 2 * (fq / p['fr']) * (p['fr'] * 1e3 / (2 * QF * D)))
    ncrit = (D / (2 * g)) ** 2
    ro = {f'{int(fr * 100)}pct_ncrit_eta{eta}_T1_{T1}us': dict(nbar=fr * ncrit, **assignment(KAPPA_R, two_chi, fr * ncrit, eta, T1))
          for fr, eta, T1 in ((0.1, 0.25, 50), (0.2, 0.5, 100))}
    return dict(g_MHz=float(g), two_chi_MHz=-two_chi, ratio_2chi_over_kr=two_chi / KAPPA_R, f_q_dressed_GHz=fq,
                Delta_MHz=float(D), n_crit=float(ncrit), T1_purcell_lorentz_s=float(T1a), T1_purcell_sete_ms=float(T1b * 1e3),
                readout=ro)


def main():
    global KAPPA_P
    ap = argparse.ArgumentParser()
    ap.add_argument('--core', default='', help='Q3D core run whose g_qr gives the realised |2chi|')
    ap.add_argument('--kappa-p', type=float, default=KAPPA_P)
    a = ap.parse_args()
    KAPPA_P = float(a.kappa_p)
    out = dict(targets=dict(kappa_r_MHz=KAPPA_R, two_chi_MHz=-KAPPA_R, kappa_p_MHz=KAPPA_P), target_g={}, realised={})
    for q, p in QUBITS.items():
        g = brentq(lambda g: spectrum(p['EJ'], p['EC'], p['fr'], g)[0] - KAPPA_R, 60, 300, xtol=1e-3)
        out['target_g'][q] = channel(q, p, g)
    if a.core:
        s = json.loads((ROOT / a.core / 'summaries' / 'summary.json').read_text(encoding='utf-8'))
        for q, p in QUBITS.items():
            pp = dict(p, EC=s['Ec_GHz'][q])
            out['realised'][q] = dict(channel(q, pp, s['readout'][q]['g_qr_MHz']), Ec_GHz=s['Ec_GHz'][q],
                                      C_qr_fF=s['readout'][q]['C_qr_fF'])
    folder = create_run(DESIGN, 'qubit', 'readout_optB_chi', entrypoint=__file__, inputs=[],
                        parameters=dict(core=a.core, qubits=QUBITS, kappa_r=KAPPA_R, kappa_p=KAPPA_P))
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    (folder / 'summaries' / 'summary.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    finish_run(folder, status='completed', summary=out)
    print('RUN_FOLDER=' + str(folder))
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
