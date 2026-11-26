# -*- coding: utf-8 -*-
r"""HFSS driven S21 (P7 -> P8) of the Design12 shared transmission lambda/4 Purcell filter + RR_Q1 / RR_Q2 only.

Fixture: build_design(readout_only=True). Rendered: feed_in (P7 side), IDC_in, PF, PF_tap, PF_short_arm, PF_short,
feed_out (P8 side), RR_Q1 / RR_Q2 with their shorts and etched open ends (no qubits, TC or control lines; the open
ends are extended by (C_claw - C_open) / C' so the bare resonators keep their target f_r). 50 Ohm lumped ports at the
P7 end of feed_in (port 1, input) and at the P8 end of feed_out (port 2, output).

Solver (skill/em_simulation.md): DrivenModal, first-order basis, explicit 10 um length mesh on every rendered metal
object, --passes adaptive passes (user 2026-11-26: 3; repository default 4), one wide Fast sweep with 1 MHz step.
Wire bonds every 350 um only on feed_in / feed_out (no bonds on the filter or resonator CPWs: no air bridges there).

Analysis (no fitted geometry): HFSS uses e^{+j w t}, so the model is fitted to conj(S21) (note/purcell_filter_method.md
section 3). S21 zeros = bare f_r (lossless resonators): the two deepest notches, lower one -> lower-target resonator.
3-mode input-output fit
    S21(f) = A e^{i 2 pi f tau} / (kappa_p / 2 - i (f - f_p) + sum_j J_j^2 / (-i (f - f_rj))) + b
gives f_p, kappa_p, J_j, f_rj; kappa_r,eff_j = J_j^2 kappa_p / ((kappa_p / 2)^2 + (f_rj - f_p)^2);
kappa_in, kappa_out from kappa_p and |A| = sqrt(kappa_in kappa_out). Driven f_r is not converged evidence
(skill/em_simulation.md section 4); use it for J, kappa and lineshape.

Usage (repository root, qmetal):
    python -X utf8 -B Design12_Transmon2Q_TransmissionPurcell/simulation/resonator/hfss_purcell_s21.py --passes 3
    python -X utf8 -B .../hfss_purcell_s21.py --analyse <run folder>      # re-analysis only
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code' / 'general'))
import project_paths
project_paths.configure(__file__)
DESIGN = project_paths.design_dir(__file__)
sys.path.insert(0, str(DESIGN / 'layout'))

import numpy as np

CHANNELS = ('Q1', 'Q2')
READOUT_COMPS = ['feed_in', 'IDC_in', 'PF', 'PF_tap', 'PF_short_arm', 'PF_short', 'feed_out',
                 'RR_Q1', 'RR_Q1_short', 'RR_Q1_open', 'RR_Q2', 'RR_Q2_short', 'RR_Q2_open']


def build_fixture():
    import build as B
    d = B.build_design(variant='sim', readout_only=True, check_keepout_um=None)
    t = d.qgeometry.tables['path']
    t['hfss_wire_bonds'] = False
    bonded = [d.components[n].id for n in ('feed_in', 'feed_out')]
    t.loc[t['component'].isin(bonded) & t['subtract'].astype(bool), 'hfss_wire_bonds'] = True
    ports = [('feed_in', 'start', 50), ('feed_out', 'end', 50)]
    return d, B, READOUT_COMPS, ports, []


# ---- analysis ------------------------------------------------------------------------------------------------------
def model(f, p):
    fp, kp, f1, f2, j1, j2, ar, ai, tau, br, bi = p
    # pole-free form: numerator and denominator multiplied by (f - f1)(f - f2) (zeros of S21 at the bare f_r)
    num = (f - f1) * (f - f2)
    den = (kp / 2 - 1j * (f - fp)) * num + 1j * j1 ** 2 * (f - f2) + 1j * j2 ** 2 * (f - f1)
    return (ar + 1j * ai) * np.exp(2j * np.pi * f * tau) * num / den + (br + 1j * bi)


def fit_s21(f, s, targets):
    from scipy.optimize import least_squares
    from scipy.signal import find_peaks
    db = 20 * np.log10(np.maximum(abs(s), 1e-12))
    idx, _ = find_peaks(-db, prominence=6)
    names = sorted(targets, key=targets.get)          # lower target <-> lower notch
    deep = sorted(sorted(idx, key=lambda i: db[i])[:len(names)], key=lambda i: f[i])
    if len(deep) < len(names):
        raise ValueError(f'only {len(deep)} S21 notches found for {len(names)} resonators')
    zeros = {q: float(f[i]) for q, i in zip(names, deep)}
    i_pk = int(np.argmax(abs(s)))
    kp0 = 0.30
    a0 = abs(s[i_pk]) * kp0 / 2
    ph0 = np.angle(s[i_pk])
    best = None
    for tau0, j0 in [(t, j) for t in np.linspace(-3, 3, 13) for j in (0.01, 0.03)]:   # tau: lead phase, ns
        p0 = [f[i_pk], kp0, zeros['Q1'], zeros['Q2'], j0, j0,
              a0 * np.cos(ph0), a0 * np.sin(ph0), tau0, 0.0, 0.0]
        lo = [f.min(), 0.02, zeros['Q1'] - 0.02, zeros['Q2'] - 0.02, 0, 0, -10, -10, -20, -1, -1]
        hi = [f.max(), 1.5, zeros['Q1'] + 0.02, zeros['Q2'] + 0.02, 0.1, 0.1, 10, 10, 20, 1, 1]

        def res(p):
            r = model(f, p) - s
            return np.r_[r.real, r.imag]
        try:
            sol = least_squares(res, p0, bounds=(lo, hi), x_scale='jac', max_nfev=20000)
        except Exception:
            continue
        if best is None or sol.cost < best.cost:
            best = sol
    p = best.x
    fp, kp, f1, f2, j1, j2 = p[:6]
    amp = abs(p[6] + 1j * p[7])
    disc = max(kp ** 2 - 4 * amp ** 2, 0.0)
    k_hi, k_lo = (kp + np.sqrt(disc)) / 2, (kp - np.sqrt(disc)) / 2
    resid = float(np.linalg.norm(model(f, p) - s) / np.linalg.norm(s))
    out = dict(f_p_GHz=float(fp), kappa_p_MHz=float(kp * 1e3),
               kappa_out_MHz_if_larger=float(k_hi * 1e3), kappa_in_MHz_if_smaller=float(k_lo * 1e3),
               peak_S21_dB=float(db[i_pk]), peak_f_GHz=float(f[i_pk]), relative_residual=resid,
               tau_ns=float(p[8]), channels={})
    for q, fr, j in (('Q1', f1, j1), ('Q2', f2, j2)):
        dlt = fr - fp
        k_eff = j ** 2 * kp / ((kp / 2) ** 2 + dlt ** 2)
        out['channels'][q] = dict(f_r_fit_GHz=float(fr), f_r_zero_GHz=zeros[q], target_f_GHz=targets[q],
                                  df_r_MHz=float((fr - targets[q]) * 1e3), J_MHz=float(j * 1e3),
                                  f_r_minus_f_p_MHz=float(dlt * 1e3), kappa_r_eff_MHz=float(k_eff * 1e3))
    return out, p


def analyse(folder: Path, targets: dict, title: str) -> dict:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd
    S = pd.read_csv(folder / 'raw' / 'S_wide.csv', index_col=0)
    conv = lambda c: S[c].map(lambda v: complex(str(v).replace('i', 'j'))).to_numpy()
    f, s21, s11 = S.index.to_numpy(float), np.conj(conv('S21')), np.conj(conv('S11'))   # e^{+jwt} -> e^{-iwt}
    out, p = fit_s21(f, s21, targets)
    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax[0].plot(f, 20 * np.log10(abs(s21)), 'k-', lw=1.4, label='HFSS |S21|')
    ax[0].plot(f, 20 * np.log10(abs(s11)), color='0.6', lw=1.0, label='HFSS |S11|')
    ax[0].plot(f, 20 * np.log10(abs(model(f, p))), 'r--', lw=1.0, label='3-mode fit')
    for q, c in out['channels'].items():
        ax[0].axvline(c['target_f_GHz'], color='C0', ls=':', lw=0.8)
        ax[0].text(c['target_f_GHz'], 1, f"{q} target", fontsize=8, ha='center')
    ax[0].set_ylabel('dB')
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3)
    ax[1].plot(f, np.unwrap(np.angle(s21)), 'k-', lw=1.2, label='HFSS arg S21')
    ax[1].plot(f, np.unwrap(np.angle(model(f, p))), 'r--', lw=1.0, label='fit')
    ax[1].set_xlabel('Frequency (GHz)')
    ax[1].set_ylabel('rad')
    ax[1].grid(alpha=0.3)
    ax[1].legend(fontsize=8)
    ax[0].set_title(title)
    fig.tight_layout()
    (folder / 'figures').mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / 'figures' / 's21_fit.png', dpi=150)
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--passes', type=int, default=3)
    ap.add_argument('--band', default='5.75,6.75')
    ap.add_argument('--tag', default='r1_3p')
    ap.add_argument('--analyse', default=None, help='existing run folder: re-run the analysis only')
    a = ap.parse_args()
    lo, hi = map(float, a.band.split(','))
    import build as B
    targets = {q: B.RR_CHANNELS[q]['f_GHz'] for q in CHANNELS}

    if a.analyse:                                      # re-analysis: new run folder referencing the parent run
        import shutil
        from run_record import create_run, finish_run
        parent = Path(a.analyse).resolve()
        raw = json.loads((parent / 'raw' / 'driven.json').read_text(encoding='utf-8'))
        folder = create_run(DESIGN.name, 'resonator', f'fit_s21_{a.tag}', entrypoint=__file__, upstream=[str(parent.relative_to(ROOT))],
                            parameters=dict(parent_run=parent.name, method='3-mode input-output fit on conj(S21)'))
        (folder / 'raw').mkdir(parents=True, exist_ok=True)
        for fn in ('S_wide.csv', 'driven.json'):
            shutil.copy2(parent / 'raw' / fn, folder / 'raw' / fn)
        out = analyse(folder, targets, f"{DESIGN.name} {parent.name[-12:]} S21 ({raw['passes']} passes, "
                                       f"dS={raw['last_delta_s']:.3g}), re-analysis")
        summary = dict(parent_run=parent.name, hfss_driven=raw, fit=out)
        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(json.dumps(summary, indent=2, default=float), encoding='utf-8')
        finish_run(folder, status='completed', summary=summary)
        print('RUN_FOLDER=' + str(folder))
        print(json.dumps(out, indent=2))
        return

    from run_record import create_run, finish_run
    import hfss_dirs
    import save_aedt
    import purcell_resonator_pair_sim as S
    folder = create_run(
        DESIGN.name, 'resonator', f'purcell_s21_driven_{a.tag}', entrypoint=__file__,
        inputs=[DESIGN / 'design.json', DESIGN / 'layout' / 'build.py'],
        parameters=dict(fixture='readout_only (PF + RR_Q1 + RR_Q2, no qubits / TC / control lines)',
                        ports='50 Ohm lumped: port 1 feed_in start (P7, input), port 2 feed_out end (P8, output)',
                        max_passes=a.passes, pass_override='user 2026-11-26: 3 passes (repository default 4)',
                        band_GHz=[lo, hi], wide_step_GHz=0.001, sweep='Fast', basis_order=1,
                        mesh_length=S.MESH, mesh_targets='every rendered metal object',
                        wire_bonds=dict(S.WB, on='feed_in, feed_out only'),
                        geometry_revision='2026-11-26_r1', pf=B.PF, idc=B.IDC_OPT, rr=B.RR_CHANNELS))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        hfss_dirs.use_d_drive()
        project = DESIGN / 'simulation' / 'aedt' / 'D12.aedt'
        project.parent.mkdir(parents=True, exist_ok=True)
        save_aedt.use_project(project)
        d, _, comps, ports, opens = build_fixture()
        raw = S.run_driven(folder, d, comps, ports, opens, f0=0.5 * (lo + hi), band=(lo, hi), passes=a.passes,
                           fine=False, wide_step_ghz=0.001, tag='S12')
        out = analyse(folder, targets, f"{DESIGN.name} driven S21 ({raw['passes']} passes, dS={raw['last_delta_s']:.3g})")
        summary = dict(hfss_driven=raw, fit=out)
        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(json.dumps(summary, indent=2, default=float), encoding='utf-8')
        finish_run(folder, status='completed', summary=summary)
        print(json.dumps(summary, indent=2, default=float), flush=True)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
