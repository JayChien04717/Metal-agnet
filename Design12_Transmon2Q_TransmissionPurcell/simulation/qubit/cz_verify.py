# -*- coding: utf-8 -*-
r"""CZ gate and idle ZZ verification for Design12_Transmon2Q_TransmissionPurcell (scqubits).

Copied from Design11_Transmon2Q_Manhattan/simulation/qubit/cz_verify.py (2026-11-26); the CZ model is unchanged.
Design12 uses the Design11 r9 core unchanged, so the input is the Design11 r9 core Q3D summary. Added for Design12:
a T1-limited CZ estimate that includes the Purcell T1 of the shared transmission filter (--purcell), because the
shared filter is what differs from Design11.

Reads a Q3D (or HFSS) summary.json containing Ec_GHz and g_n_MHz, computes:
  1. Static dressed spectrum & idle ZZ vs coupler flux (qubits at sweet spot).
  2. Flux-activated CZ (Q1 + TC cosine-ramp flat-top pulse) across requested gate times.
  3. Saves run.json, summaries/cz_summary.json, figures/zz_vs_flux.png, figures/cz_vs_time.png.

Environment: scqenv (scqubits 4.3).
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code' / 'general'))
import project_paths
project_paths.configure(__file__)
from run_record import create_run, finish_run

import numpy as np
import scqubits as scq
from scipy.linalg import eigh
from scipy.optimize import brentq, minimize

scq.settings.T1_DEFAULT_WARNING = False

DESIGN_DIR = Path(__file__).resolve().parents[2]
E_CHARGE = 1.602176634e-19
H_PLANCK = 6.62607015e-34
PHI0 = H_PLANCK / (2 * E_CHARGE)
NCUT, NQ, NC = 15, 4, 5
DT = 0.05
FC_MARGIN = 0.40  # GHz
DEFAULT_LJ = dict(Q1a=17.39, Q1b=34.79, Q2a=20.55, Q2b=41.11, TC=14.03)


def ej_ghz(L_nH: float) -> float:
    return (PHI0 / (2 * np.pi)) ** 2 / (L_nH * 1e-9) / H_PLANCK / 1e9


def retarget_lj(Ec: dict, lj: dict, targets: dict) -> dict:
    out = dict(lj)
    for n, keys in (('Q1', ('Q1a', 'Q1b')), ('Q2', ('Q2a', 'Q2b')), ('TC', ('TC', 'TC'))):
        f_t = targets[n]
        ejs = brentq(lambda x: scq.Transmon(EJ=x, EC=Ec[n], ng=0, ncut=NCUT).E01() - f_t, 2.0, 120.0)
        if n == 'TC':
            out['TC'] = (PHI0 / (2 * np.pi)) ** 2 / (ejs / 2 * 1e9 * H_PLANCK) * 1e9
        else:
            ra = ej_ghz(lj[keys[0]]) / (ej_ghz(lj[keys[0]]) + ej_ghz(lj[keys[1]]))
            out[keys[0]] = (PHI0 / (2 * np.pi)) ** 2 / (ra * ejs * 1e9 * H_PLANCK) * 1e9
            out[keys[1]] = (PHI0 / (2 * np.pi)) ** 2 / ((1 - ra) * ejs * 1e9 * H_PLANCK) * 1e9
    return out


def load_params(path: str, retarget=None) -> dict:
    r = json.loads(Path(path).read_text(encoding='utf-8'))
    gg = r['g_n_MHz']
    g = {
        'Q1-TC': gg.get('Q1-TC', gg.get('TC-Q1')),
        'Q2-TC': gg.get('Q2-TC', gg.get('TC-Q2')),
        'Q1-Q2': gg.get('Q1-Q2', gg.get('Q2-Q1')),
    }
    ec = r.get('Ec_GHz') or {k: v / 1e3 for k, v in r['Ec_MHz'].items()}
    lj = {k: float(str(v).replace('nH', '')) for k, v in (r.get('Lj_drawn_nH') or DEFAULT_LJ).items()}
    if retarget:
        lj = retarget_lj(ec, lj, retarget)
    ej = {k: ej_ghz(v) for k, v in lj.items()}
    return dict(
        source=r.get('source', 'Q3D capacitance matrix'),
        file=str(path),
        Ec_GHz=dict(ec),
        g_n_GHz={k: v / 1e3 for k, v in g.items()},
        Lj_nH=lj,
        EJsum_GHz=dict(Q1=ej['Q1a'] + ej['Q1b'], Q2=ej['Q2a'] + ej['Q2b'], TC=2 * ej['TC']),
        d=dict(
            Q1=(ej['Q1a'] - ej['Q1b']) / (ej['Q1a'] + ej['Q1b']),
            Q2=(ej['Q2a'] - ej['Q2b']) / (ej['Q2a'] + ej['Q2b']),
            TC=0.0,
        ),
    )


class Model:
    def __init__(self, p: dict):
        self.p = p
        E, J = p['Ec_GHz'], p['EJsum_GHz']
        mk = lambda n, nt: scq.TunableTransmon(
            EJmax=J[n], EC=E[n], d=p['d'][n], flux=0.0, ng=0, ncut=NCUT, truncated_dim=nt
        )
        self.sys = dict(Q1=mk('Q1', NQ), TC=mk('TC', NC), Q2=mk('Q2', NQ))
        self.dims = (NQ, NC, NQ)
        ops = {}
        for n, s in self.sys.items():
            ev, V = s.eigensys(evals_count=s.truncated_dim)
            pr = lambda op: V.conj().T @ op @ V
            cos = pr(s.cos_phi_operator())
            ops[n] = dict(H0=np.diag(ev), n=pr(s.n_operator()), cos=cos, kin=np.diag(ev) + J[n] * cos)
        I = {n: np.eye(self.sys[n].truncated_dim) for n in self.sys}
        k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
        G = p['g_n_GHz']
        self.A = (
            k3(ops['Q1']['kin'], I['TC'], I['Q2'])
            + k3(I['Q1'], ops['TC']['kin'], I['Q2'])
            + k3(I['Q1'], I['TC'], ops['Q2']['kin'])
            + G['Q1-TC'] * k3(ops['Q1']['n'], ops['TC']['n'], I['Q2'])
            + G['Q2-TC'] * k3(I['Q1'], ops['TC']['n'], ops['Q2']['n'])
            + G['Q1-Q2'] * k3(ops['Q1']['n'], I['TC'], ops['Q2']['n'])
        )
        self.B = dict(
            Q1=-k3(ops['Q1']['cos'], I['TC'], I['Q2']),
            TC=-k3(I['Q1'], ops['TC']['cos'], I['Q2']),
            Q2=-k3(I['Q1'], I['TC'], ops['Q2']['cos']),
        )

    def EJ(self, n: str, flux: float) -> float:
        p = self.p
        return p['EJsum_GHz'][n] * np.sqrt(np.cos(np.pi * flux) ** 2 + p['d'][n] ** 2 * np.sin(np.pi * flux) ** 2)

    def H(self, f1=0.0, fc=0.0, f2=0.0):
        return self.A + self.EJ('Q1', f1) * self.B['Q1'] + self.EJ('TC', fc) * self.B['TC'] + self.EJ('Q2', f2) * self.B['Q2']

    def dressed(self, **fl):
        w, v = eigh(self.H(**fl))
        labels, used = {}, set()
        for s in [(i, c, j) for i in range(3) for c in range(2) for j in range(3)]:
            idx = np.ravel_multi_index(s, self.dims)
            for k in np.argsort(-np.abs(v[idx, :]) ** 2):
                if k not in used:
                    labels[s] = k
                    used.add(k)
                    break
        return w, v, labels

    def static(self, **fl):
        w, _, lab = self.dressed(**fl)
        E = lambda s: w[lab[s]] - w[lab[(0, 0, 0)]]
        return dict(
            f1=float(E((1, 0, 0))),
            f2=float(E((0, 0, 1))),
            fc=float(E((0, 1, 0))),
            zz=float(E((1, 0, 1)) - E((1, 0, 0)) - E((0, 0, 1))),
            a1=float(E((2, 0, 0)) - 2 * E((1, 0, 0))),
            a2=float(E((0, 0, 2)) - 2 * E((0, 0, 1))),
        )


def envelope(t, T, ramp):
    t = np.asarray(t, float)
    w = np.ones_like(t)
    r1, r2 = t < ramp, t > T - ramp
    w[r1] = 0.5 * (1 - np.cos(np.pi * t[r1] / ramp))
    w[r2] = 0.5 * (1 - np.cos(np.pi * (T - t[r2]) / ramp))
    return np.clip(w, 0, 1)


COMP = [(0, 0, 0), (0, 0, 1), (1, 0, 0), (1, 0, 1)]


def evolve(m: Model, idle: float, a_q: float, a_c: float, T: float, ramp: float, record=False):
    w0, v0, lab = m.dressed(fc=idle)
    tm = (np.arange(int(round(T / DT))) + 0.5) * DT
    env = envelope(tm, T, ramp)
    psi = np.stack([v0[:, lab[s]] for s in COMP], axis=1).astype(complex)
    rec = []
    if record:
        for k, wk in enumerate(env):
            ev, V = eigh(m.H(f1=a_q * wk, fc=idle + a_c * wk))
            psi = V @ (np.exp(-2j * np.pi * ev * DT)[:, None] * (V.conj().T @ psi))
            if k % 10 == 0:
                rec.append((tm[k], np.abs(v0.conj().T @ psi[:, 3]) ** 2))
    else:
        changes = np.flatnonzero(np.diff(env)) + 1
        starts = np.r_[0, changes]
        counts = np.diff(np.r_[starts, len(env)])
        eigensystems = {}
        for k, count in zip(starts, counts):
            wk = env[k]
            key = round(float(wk), 14)
            if key not in eigensystems:
                eigensystems[key] = eigh(m.H(f1=a_q * wk, fc=idle + a_c * wk))
            ev, V = eigensystems[key]
            psi = V @ (np.exp(-2j * np.pi * ev * (count * DT))[:, None] * (V.conj().T @ psi))
    c = (v0.conj().T @ psi) * np.exp(2j * np.pi * w0 * T)[:, None]
    U = np.array([[c[lab[s], k] for k in range(4)] for s in COMP])
    if record:
        sel = {'11': (1, 0, 1), '02': (0, 0, 2), '20': (2, 0, 0), '1,c,0': (1, 1, 0), '0,c,1': (0, 1, 1)}
        pops = dict(t=[r[0] for r in rec], **{k: [float(r[1][lab[s]]) for r in rec] for k, s in sel.items()})
        return U, pops
    return U


def analyse_U(U):
    ph = np.angle(np.diag(U))
    cond = (ph[3] - ph[2] - ph[1] + ph[0]) % (2 * np.pi)
    leak = 1 - np.mean(np.sum(np.abs(U) ** 2, axis=0))
    Z = np.diag(np.exp(-1j * np.array([ph[0], ph[1], ph[2], ph[1] + ph[2] - ph[0]])))
    Uc = Z @ U
    Ucz = np.diag([1, 1, 1, -1]).astype(complex)
    F = (abs(np.trace(Ucz.conj().T @ Uc)) ** 2 + np.trace(Uc.conj().T @ Uc).real) / 20
    return dict(cond_phase_deg=float(np.degrees(cond)), leakage=float(leak), fidelity=float(F))


def optimise(m: Model, idle: float, T: float, ramp: float, a_q0: float, a_c_max: float):
    unpack = lambda x: (x[0], a_c_max * np.sin(x[1]) ** 2)

    def cost(x):
        aq, ac = unpack(x)
        return 1 - analyse_U(evolve(m, idle, aq, ac, T, ramp))['fidelity']

    us = np.linspace(0.25, 1.45, 10)
    x0 = min([[a_q0, u] for u in us], key=cost)
    simplex = np.array([x0, [x0[0] + 0.004, x0[1]], [x0[0], x0[1] + 0.08]])
    sol = minimize(cost, x0, method='Nelder-Mead',
                   options=dict(xatol=1e-5, fatol=1e-7, maxiter=150, initial_simplex=simplex))
    aq, ac = unpack(sol.x)
    r = analyse_U(evolve(m, idle, aq, ac, T, ramp))
    return dict(T_ns=T, ramp_ns=ramp, a_q=float(aq), a_c=float(ac), **r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--params', required=True)
    ap.add_argument('--tag', default='q3d')
    ap.add_argument('--times', default='55,65,70,75,80,85')
    ap.add_argument('--ramp', type=float, default=8.0)
    ap.add_argument('--retarget', default='', help='Optional f01 targets Q1,Q2,TC (GHz)')
    ap.add_argument('--purcell', default='',
                    help='JSON: {"f_pf": GHz, "kappa_pf_MHz": .., "kappa_r_MHz": .., "f_r": {"Q1": GHz, "Q2": GHz}, '
                         '"source": ".."}; adds Purcell T1 and T1-limited CZ estimates')
    ap.add_argument('--t1-assumed-us', default='50,100')
    a = ap.parse_args()

    tg = dict(zip(('Q1', 'Q2', 'TC'), map(float, a.retarget.split(',')))) if a.retarget else None
    p = load_params(a.params, tg)
    p['retarget_GHz'] = tg
    m = Model(p)

    run_dir = create_run(
        project_paths.design_dir(__file__).name,
        'qubit',
        f'cz_verify_{a.tag}',
        entrypoint=__file__,
        inputs=[Path(a.params)],
        parameters={'retarget_GHz': tg, 'times_ns': a.times, 'ramp_ns': a.ramp,
                    'purcell': json.loads(a.purcell) if a.purcell else None, 't1_assumed_us': a.t1_assumed_us},
    )
    (run_dir / 'summaries').mkdir(parents=True, exist_ok=True)
    (run_dir / 'figures').mkdir(parents=True, exist_ok=True)

    fl = np.linspace(0, 0.40, 81)
    st = [m.static(fc=f) for f in fl]
    zz = np.array([s['zz'] for s in st]) * 1e6
    fcs = np.array([s['fc'] for s in st])
    idle = 0.0
    sgn = [i for i in range(len(fl) - 1) if np.sign(zz[i]) != np.sign(zz[i + 1]) and fcs[i + 1] > st[0]['f1'] + 0.8]
    if sgn:
        i = sgn[0]
        idle = float(brentq(lambda f: m.static(fc=f)['zz'], fl[i], fl[i + 1]))
    si = m.static(fc=idle)
    print(f"idle: coupler_flux={idle:.4f}, f1={si['f1']:.4f}, f2={si['f2']:.4f}, fc={si['fc']:.4f} GHz, "
          f"alpha={si['a1']*1e3:.1f}/{si['a2']*1e3:.1f} MHz, ZZ={si['zz']*1e6:.2f} kHz", flush=True)

    target = si['f2'] - si['a1']
    fq = np.linspace(0, 0.3, 61)
    f1s = np.array([m.static(f1=x, fc=idle)['f1'] for x in fq])
    a_q0 = float(np.interp(target, f1s[::-1], fq[::-1]))
    ok = [f for f, fc_ in zip(fl, fcs) if fc_ >= target + FC_MARGIN and f >= idle]
    acq = max(ok) - idle

    res = []
    for T in [float(x) for x in a.times.split(',')]:
        r = optimise(m, idle, T, a.ramp, a_q0, acq)
        sg = m.static(f1=r['a_q'], fc=idle + r['a_c'])
        r.update(f1_gate_GHz=sg['f1'], fc_gate_GHz=sg['fc'], at_coupler_bound=bool(r['a_c'] > 0.99 * acq))
        res.append(r)
        print(f"[Q1+TC] T={T:5.1f} ns: a_q={r['a_q']:.4f} (f1={sg['f1']:.3f}), a_c={r['a_c']:+.4f} (fc={sg['fc']:.3f}) -> "
              f"phase={r['cond_phase_deg']:.2f} deg, leak={r['leakage']:.2e}, F={r['fidelity']:.6f}", flush=True)

    best = max(res, key=lambda r: r['fidelity'])
    _, pops = evolve(m, idle, best['a_q'], best['a_c'], best['T_ns'], a.ramp, record=True)

    # Purcell T1 through the shared band-pass filter (Sete 2015 / Jeffrey 2014 bound, estimate only):
    #   gamma = kappa_r (g / Delta_qr)^2 * f_q / (2 Q_F Delta_qf),  Q_F = f_PF / kappa_PF
    # and the T1-limited CZ estimate of independent amplitude damping on both qubits during the gate:
    #   eps_T1 = (2 / 5) T (1 / T1_Q1 + 1 / T1_Q2)   ((4 / 5) T / T1 for equal T1; d = 4), 1 / T1 = 1 / T1_int + 1 / T1_P.
    purcell = None
    if a.purcell:
        pc = json.loads(a.purcell)
        rd = json.loads(Path(a.params).read_text(encoding='utf-8')).get('readout', {})
        fq = dict(Q1=si['f1'], Q2=si['f2'])
        qf = pc['f_pf'] / (pc['kappa_pf_MHz'] * 1e-3)
        t1p = {}
        for q in ('Q1', 'Q2'):
            g = rd[q]['g_qr_MHz'] * 1e-3
            d_qr = pc['f_r'][q] - fq[q]
            kr = pc['kappa_r_MHz'][q] if isinstance(pc['kappa_r_MHz'], dict) else pc['kappa_r_MHz']   # scalar or per qubit
            gam = kr * 1e-3 * (g / d_qr) ** 2 * fq[q] / (2 * qf * (pc['f_pf'] - fq[q]))   # GHz
            t1p[q] = 1.0 / (2 * np.pi * gam * 1e9) * 1e6                                             # us
        table = []
        for t1i in [float(x) for x in a.t1_assumed_us.split(',')] + [float('inf')]:
            t1 = {q: 1.0 / (1.0 / t1i + 1.0 / t1p[q]) for q in t1p}
            eps = 0.4 * best['T_ns'] * 1e-3 * (1 / t1['Q1'] + 1 / t1['Q2'])
            table.append(dict(T1_intrinsic_us=t1i, T1_total_us=t1, eps_T1=eps,
                              F_T1_limited=best['fidelity'] - eps))
        purcell = dict(inputs=pc, Q_F=qf, f01_idle_GHz=fq, T1_purcell_us=t1p, cz_T_ns=best['T_ns'],
                       cz_coherent_F=best['fidelity'], t1_limited=table,
                       method='Sete/Jeffrey band-pass Purcell bound; eps_T1 = (2/5) T (1/T1_Q1 + 1/T1_Q2); estimate, '
                              'not a simulation (no T_phi, no |2> / coupler decay, no flux noise)')
        print(f"Purcell T1 (shared filter, Q_F={qf:.1f}): Q1 {t1p['Q1']:.0f} us, Q2 {t1p['Q2']:.0f} us", flush=True)
        for r in table:
            print(f"  T1_int={r['T1_intrinsic_us']:g} us -> T1 {r['T1_total_us']['Q1']:.1f}/{r['T1_total_us']['Q2']:.1f} us, "
                  f"CZ F_T1-limited = {r['F_T1_limited']:.5f}", flush=True)

    summary = dict(
        params=p,
        idle=dict(coupler_flux=idle, **si),
        Q1_flux_for_11_20=a_q0,
        Q1_gate_target_GHz=target,
        cz_q1_plus_coupler=res,
        best=best,
        purcell_t1_estimate=purcell,
    )
    (run_dir / 'summaries' / 'cz_summary.json').write_text(json.dumps(summary, indent=2, default=float), encoding='utf-8')
    finish_run(run_dir, status='completed', summary=dict(idle=summary['idle'], best=best, purcell_t1_estimate=purcell))

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(fl, fcs, label='TC')
    ax[0].plot(fl, [s['f1'] for s in st], label='Q1')
    ax[0].plot(fl, [s['f2'] for s in st], label='Q2')
    ax[0].axvline(idle, ls=':', c='k')
    ax[0].set_xlabel(r'Coupler flux ($\Phi_0$), qubits at sweet spot')
    ax[0].set_ylabel('Dressed frequency (GHz)')
    ax[0].legend()
    ax[1].semilogy(fl, np.abs(zz))
    ax[1].axvline(idle, ls=':', c='k')
    ax[1].set_xlabel(r'Coupler flux ($\Phi_0$)')
    ax[1].set_ylabel(r'$|ZZ|/2\pi$ (kHz)')
    fig.suptitle(f'{project_paths.design_dir(__file__).name} static ZZ ({a.tag})')
    fig.tight_layout()
    fig.savefig(run_dir / 'figures' / 'zz_vs_flux.png', dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].semilogy([r['T_ns'] for r in res], [1 - r['fidelity'] for r in res], 'o-', label='1 - F (Q1 + TC)')
    ax[0].semilogy([r['T_ns'] for r in res], [max(r['leakage'], 1e-9) for r in res], 's--', label='leakage')
    ax[0].set_xlabel('Gate time T (ns)')
    ax[0].legend()
    ax[0].set_title('CZ coherent error')
    for k, v in pops.items():
        if k != 't':
            ax[1].plot(pops['t'], v, label=f'|{k}>')
    ax[1].set_xlabel('t (ns)')
    ax[1].set_ylabel('Population, start |11>')
    ax[1].legend()
    ax[1].set_title(f"best: T={best['T_ns']:.0f} ns, F={best['fidelity']:.5f}")
    fig.tight_layout()
    fig.savefig(run_dir / 'figures' / 'cz_vs_time.png', dpi=140)
    plt.close(fig)
    print('CZ_VERIFY_OK:', run_dir)


if __name__ == '__main__':
    main()
