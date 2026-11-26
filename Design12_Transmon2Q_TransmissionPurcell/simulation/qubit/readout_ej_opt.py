# -*- coding: utf-8 -*-
r"""Best qubit E_J for the decided readout resonators (|2chi| = kappa_r) and thermal-photon shot-noise T_phi.

User request (2026-11-26): the resonator frequencies are decided (r4), so find the E_J that is best for readout and
compute the photon shot-noise dephasing. E_J is a design input; the user explicitly authorised this retarget.

Model (scqubits exact transmon + oscillator, as Design11 simulation/qubit/readout_optB_chi.py):
  * E_C, C_qr, C_sigma: core Q3D (Design11 r9 core, unchanged in Design12).
  * g(f_q, f_r) = C_qr / (2 sqrt(C_sigma C_r)) sqrt(f_q f_r), C_r = 1 / (8 f_r Z0) (lambda/4, Z0 = 50 Ohm), f_q = bare
    sweet-spot f01 (the convention of the Q3D core run, which reproduces its g_qr).
  * |2chi| = (E11 - E10) - (E01 - E00) of the dressed spectrum.
  * Best E_J: |2chi| = kappa_r (SNR-optimal dispersive readout). The SQUID branch ratio (2:1) is kept; E_J sum scaled.
  * Purcell T1 (shared band-pass filter, Sete / Jeffrey bound): gamma = kappa_r (g/Delta)^2 f_q / (2 Q_F (f_PF - f_q)).
  * Thermal-photon dephasing (Clerk & Utami PRA 75 042302; Rigetti et al. PRB 86 100506), chi = |2chi| / 2:
        Gamma_phi = (kappa / 2) Re[ sqrt((1 + 2 i chi / kappa)^2 + 8 i chi n_th / kappa) - 1 ]
    (-> n_th kappa (2chi)^2 / (kappa^2 + (2chi)^2) for n_th << 1), evaluated over a list of n_th.

Environment: scqenv.
    & "C:\Users\QEL\miniforge-pypy3\envs\scqenv\python.exe" -X utf8 -B Design12_Transmon2Q_TransmissionPurcell/simulation/qubit/readout_ej_opt.py --inputs <json>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code' / 'general'))
import project_paths
project_paths.configure(__file__)

import numpy as np
import scqubits as scq
from scipy.optimize import brentq

from run_record import create_run, finish_run

scq.settings.T1_DEFAULT_WARNING = False
TP = 2 * np.pi
Z0 = 50.0
H, KB = 6.62607015e-34, 1.380649e-23


def g_mhz(c_qr, c_sig, fq, fr):
    c_r = 1.0 / (8 * fr * 1e9 * Z0) * 1e15                     # fF
    return c_qr / (2 * np.sqrt(c_sig * c_r)) * np.sqrt(fq * fr) * 1e3


def transmon(EJ, EC):
    return scq.Transmon(EJ=EJ, EC=EC, ng=0, ncut=30, truncated_dim=6)


def dressed(EJ, EC, fr, c_qr, c_sig):
    tm = transmon(EJ, EC)
    ev = tm.eigenvals(evals_count=3)
    fq_bare, alpha = ev[1] - ev[0], (ev[2] - ev[1]) - (ev[1] - ev[0])
    g = g_mhz(c_qr, c_sig, fq_bare, fr)
    n01 = abs(tm.matrixelement_table('n_operator', evals_count=2)[0, 1])
    osc = scq.Oscillator(E_osc=fr, truncated_dim=6)
    hs = scq.HilbertSpace([tm, osc])
    hs.add_interaction(g=g * 1e-3 / n01, op1=tm.n_operator, op2=osc.creation_operator, add_hc=True)
    hs.generate_lookup()
    E = lambda q, n: hs.energy_by_bare_index((q, n))
    two_chi = ((E(1, 1) - E(1, 0)) - (E(0, 1) - E(0, 0))) * 1e3
    return dict(EJ=EJ, f_q_bare_GHz=float(fq_bare), f_q_dressed_GHz=float(E(1, 0) - E(0, 0)), alpha_MHz=float(alpha * 1e3),
                g_MHz=float(g), two_chi_MHz=float(two_chi), Delta_MHz=float((fr - fq_bare) * 1e3))


def gamma_phi(kappa_mhz, two_chi_mhz, n_th):
    k, chi = TP * kappa_mhz * 1e6, TP * abs(two_chi_mhz) / 2 * 1e6
    return float((k / 2) * np.real(np.sqrt((1 + 2j * chi / k) ** 2 + 8j * chi * n_th / k) - 1))


def n_th_at(f_ghz, t_mk):
    return 1.0 / np.expm1(H * f_ghz * 1e9 / (KB * t_mk * 1e-3))


def analyse(q, p, core, pf, n_list, t_list):
    EC, c_qr, c_sig = core['Ec_GHz'][q], core['readout'][q]['C_qr_fF'], core['C_sigma_fF'][q]
    fr, kr = p['f_r_GHz'], p['kappa_r_MHz']
    f = lambda EJ: abs(dressed(EJ, EC, fr, c_qr, c_sig)['two_chi_MHz']) - kr
    # scan upward from the dispersive side (f_q well below f_r) and take the first |2chi| = kappa crossing; the
    # straddling regime near f_q ~ f_r is excluded (Delta >= 600 MHz)
    grid = [EJ for EJ in np.arange(5.0, 30.0, 0.5)
            if p['f_r_GHz'] - (transmon(EJ, EC).eigenvals(2)[1] - transmon(EJ, EC).eigenvals(2)[0]) >= 0.6]
    vals = [f(EJ) for EJ in grid]
    i = next(k for k in range(len(grid) - 1) if vals[k] < 0 <= vals[k + 1])
    EJ_best = brentq(f, grid[i], grid[i + 1], xtol=1e-4)
    best = dressed(EJ_best, EC, fr, c_qr, c_sig)
    qf = pf['f_GHz'] / (pf['kappa_MHz'] * 1e-3)

    def extras(d):
        D = d['Delta_MHz']
        t1p = 1 / (TP * kr * 1e6 * (d['g_MHz'] / D) ** 2 * d['f_q_bare_GHz'] / (2 * qf * (pf['f_GHz'] - d['f_q_bare_GHz'])))
        tphi = {f'n_th={n:g}': 1e6 / gamma_phi(kr, d['two_chi_MHz'], n) for n in n_list}
        tphi_T = {f'T_eff={t:g}mK (n_th={n_th_at(fr, t):.2e})': 1e6 / gamma_phi(kr, d['two_chi_MHz'], n_th_at(fr, t))
                  for t in t_list}
        return dict(n_crit=float((D / (2 * d['g_MHz'])) ** 2), T1_purcell_us=float(t1p * 1e6),
                    T_phi_shot_us=tphi, T_phi_shot_us_by_temperature=tphi_T)

    best.update(extras(best), ratio_2chi_over_kappa=abs(best['two_chi_MHz']) / kr)
    ra = p['branch_ratio']
    best['EJ_branch_GHz'] = [EJ_best * ra / (1 + ra), EJ_best / (1 + ra)]
    cur = dressed(p['EJ_current'], EC, fr, c_qr, c_sig)
    cur.update(extras(cur), ratio_2chi_over_kappa=abs(cur['two_chi_MHz']) / kr)
    sweep = []
    for EJ in np.linspace(0.85 * EJ_best, 1.15 * EJ_best, 7):
        d = dressed(EJ, EC, fr, c_qr, c_sig)
        e = extras(d)
        sweep.append(dict(EJ=float(EJ), f_q_bare_GHz=d['f_q_bare_GHz'], two_chi_MHz=d['two_chi_MHz'],
                          T1_purcell_us=e['T1_purcell_us'], T_phi_shot_us_n1e3=1e6 / gamma_phi(kr, d['two_chi_MHz'], 1e-3)))
    return dict(inputs=dict(p, EC_GHz=EC, C_qr_fF=c_qr, C_sigma_fF=c_sig), best=best, current_EJ=cur, sweep=sweep)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', required=True, help='JSON file with core, qubits, filter, n_th, T_eff lists')
    ap.add_argument('--tag', default='r4')
    a = ap.parse_args()
    cfg = json.loads(Path(a.inputs).read_text(encoding='utf-8'))
    core = json.loads((ROOT / cfg['core_summary']).read_text(encoding='utf-8'))
    folder = create_run(project_paths.design_dir(__file__).name, 'qubit', f'readout_ej_opt_{a.tag}', entrypoint=__file__,
                        inputs=[Path(a.inputs), ROOT / cfg['core_summary']], parameters=cfg)
    out = dict(method='|2chi| = kappa_r (scqubits dressed); Purcell Sete/Jeffrey; shot-noise Clerk-Utami exact',
               qubits={q: analyse(q, p, core, cfg['filter'], cfg['n_th'], cfg['T_eff_mK']) for q, p in cfg['qubits'].items()})
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    (folder / 'summaries' / 'summary.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    finish_run(folder, status='completed', summary=out)
    print('RUN_FOLDER=' + str(folder))
    for q, r in out['qubits'].items():
        b, c = r['best'], r['current_EJ']
        print(f"{q}: current EJ {c['EJ']:.2f} GHz -> f_q {c['f_q_bare_GHz']:.4f}, 2chi {c['two_chi_MHz']:+.3f} MHz "
              f"(|2chi|/kappa {c['ratio_2chi_over_kappa']:.2f}), T1_P {c['T1_purcell_us']:.0f} us")
        print(f"    best EJ {b['EJ']:.3f} GHz (branches {b['EJ_branch_GHz'][0]:.3f} / {b['EJ_branch_GHz'][1]:.3f}) -> "
              f"f_q bare {b['f_q_bare_GHz']:.4f} / dressed {b['f_q_dressed_GHz']:.4f} GHz, alpha {b['alpha_MHz']:.1f} MHz, "
              f"g {b['g_MHz']:.1f} MHz, 2chi {b['two_chi_MHz']:+.3f} MHz, n_crit {b['n_crit']:.1f}, T1_P {b['T1_purcell_us']:.0f} us")
        for k, v in b['T_phi_shot_us'].items():
            print(f"      T_phi shot ({k}): {v:.1f} us")
        for k, v in b['T_phi_shot_us_by_temperature'].items():
            print(f"      T_phi shot ({k}): {v:.1f} us")


if __name__ == '__main__':
    main()
