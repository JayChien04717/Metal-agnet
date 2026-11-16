"""CZ gate of Design6 (three flux-tunable transmons: Q1 - TC - Q2) with scqubits, from the HFSS-derived Hamiltonian.

Hamiltonian (charge basis, scqubits operators):
    H = sum_i [4 Ec_i n_i^2 - EJ_i(Phi_i) cos phi_i]
        + g_1c n_1 n_c + g_2c n_2 n_c + g_12 n_1 n_2            (g = 4 e^2 Cinv_ij / h, "charge units")
    Q1, Q2: asymmetric SQUIDs, EJ(Phi) = EJsum sqrt(cos^2(pi Phi) + d^2 sin^2(pi Phi)), d = (EJa - EJb)/(EJa + EJb)
    TC    : symmetric SQUID, EJ(Phi) = EJmax |cos(pi Phi)|
Inputs (--params): results/<date>/hfss_core_<tag>/fit.json (HFSS eigenmode sweep of the coupler Lj -> fitted Cinv,
    i.e. Ec_i and g_ij), EJ from the drawn junction inductances (EJ = (Phi0/2pi)^2 / Lj). A Q3D summary.json is
    accepted too (comparison).

Steps
 1. static: dressed spectrum and ZZ = E11 - E10 - E01 + E00 vs coupler flux, qubits at their sweet spots (idle).
    Idle coupler flux = zero ZZ if one exists with f_c >= f_Q1 + 0.8 GHz, else the coupler sweet spot.
 2. CZ, two controls with the same flat-top envelope w(t) (cosine ramps t_ramp, total T):
        Q1 flux   Phi_1(t) = a_q  w(t)   (Q1 down to the |11> <-> |20> condition f_Q1 = f_Q2 - alpha_1, no crossing)
        TC flux   Phi_c(t) = Phi_idle + a_c w(t)   (turns the effective Q1-Q2 coupling on)
    Propagation: piecewise-constant (dt = 0.05 ns) in a fixed truncated basis (Q eigenstates x TC eigenstates x Q
    eigenstates at zero flux): H(t) = A + EJ1(t) B1 + EJc(t) Bc + EJ2 B2, B = -cos(phi) (exact operators).
    Gate in the dressed computational basis at idle (rotating frame), virtual Z removed;
    conditional phase, leakage, F_avg = (|Tr(U_CZ^dag U)|^2 + Tr(U^dag U)) / 20.
    (a_q, a_c) optimised (Nelder-Mead on 1 - F) for each T. Also the coupler-only CZ (a_q = 0) for comparison.

Output: results/<date>/cz_<tag>/ (cz_summary.json, zz_vs_flux.png, cz_vs_time.png, cz_populations.png)
Usage:  python cz_scqubits.py --params <fit.json> [--tag hfss] [--times 40,50,60,80]
Python: conda env scqenv (scqubits 4.3).
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse, json, os, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
import scqubits as scq
from scipy.linalg import eigh
from scipy.optimize import brentq, minimize

scq.settings.T1_DEFAULT_WARNING = False
SIM = _S('Design6/simulation')
e, h = 1.602176634e-19, 6.62607015e-34
PHI0 = h / (2 * e)
NCUT, NQ, NC = 15, 4, 5
DT = 0.05          # ns (exact exponential of H on each step)
FC_MARGIN = 0.4    # GHz, the coupler stays this far above the (pulsed) Q1 during a gate (Q1 and TC pulsed together)
DEFAULT_LJ = dict(Q1a=15.89, Q1b=31.77, Q2a=18.78, Q2b=37.56, TC=12.85)   # design6_layout.LJ (nH)


def ej_ghz(L_nH):
    return (PHI0 / 2 / np.pi) ** 2 / (L_nH * 1e-9) / h / 1e9


def retarget_lj(Ec, lj, targets):
    """Junction inductances (keeping the a:b ratio) so that the bare f01 at the sweet spot = targets (GHz), Q1/Q2/TC."""
    out = dict(lj)
    for n, keys in (('Q1', ('Q1a', 'Q1b')), ('Q2', ('Q2a', 'Q2b')), ('TC', ('TC', 'TC'))):
        f_t = targets[n]
        g = lambda ejs: scq.Transmon(EJ=ejs, EC=Ec[n], ng=0, ncut=NCUT).E01() - f_t
        ejs = brentq(g, 2.0, 120.0)
        if n == 'TC':
            out['TC'] = (PHI0 / 2 / np.pi) ** 2 / (ejs / 2 * 1e9 * h) * 1e9
        else:
            ra = ej_ghz(lj[keys[0]]) / (ej_ghz(lj[keys[0]]) + ej_ghz(lj[keys[1]]))
            out[keys[0]] = (PHI0 / 2 / np.pi) ** 2 / (ra * ejs * 1e9 * h) * 1e9
            out[keys[1]] = (PHI0 / 2 / np.pi) ** 2 / ((1 - ra) * ejs * 1e9 * h) * 1e9
    return out


def load_params(path, retarget=None):
    r = json.load(open(path))
    if 'Cinv_per_F' in r:
        src = 'HFSS eigenmode (coupler Lj sweep) -> fitted Cinv'
        g = r['g_n_MHz']; g = {'Q1-TC': g['Q1-TC'], 'Q2-TC': g['Q2-TC'], 'Q1-Q2': g['Q1-Q2']}
    else:
        src = 'Q3D capacitance matrix'
        gg = r['g_n_MHz']; g = {'Q1-TC': gg['Q1-TC'], 'Q2-TC': gg['TC-Q2'], 'Q1-Q2': gg['Q1-Q2']}
    lj = {k: float(str(v).replace('nH', '')) for k, v in (r.get('Lj_drawn_nH') or DEFAULT_LJ).items()}
    if retarget:          # re-choose the junctions (fab parameter) for the target sweet-spot frequencies with these Ec
        lj = retarget_lj(r['Ec_GHz'], lj, retarget)
    ej = {k: ej_ghz(v) for k, v in lj.items()}
    p = dict(source=src, file=path, Ec_GHz=dict(r['Ec_GHz']), g_n_GHz={k: v / 1e3 for k, v in g.items()}, Lj_nH=lj,
             EJsum_GHz=dict(Q1=ej['Q1a'] + ej['Q1b'], Q2=ej['Q2a'] + ej['Q2b'], TC=2 * ej['TC']),
             d=dict(Q1=(ej['Q1a'] - ej['Q1b']) / (ej['Q1a'] + ej['Q1b']), Q2=(ej['Q2a'] - ej['Q2b']) / (ej['Q2a'] + ej['Q2b']),
                    TC=0.0))
    return p


class Model:
    def __init__(self, p):
        self.p = p
        E, J = p['Ec_GHz'], p['EJsum_GHz']
        mk = lambda n, nt: scq.TunableTransmon(EJmax=J[n], EC=E[n], d=p['d'][n], flux=0.0, ng=0, ncut=NCUT, truncated_dim=nt)
        self.sys = dict(Q1=mk('Q1', NQ), TC=mk('TC', NC), Q2=mk('Q2', NQ))
        self.dims = (NQ, NC, NQ)
        ops = {}
        for n, s in self.sys.items():
            ev, V = s.eigensys(evals_count=s.truncated_dim)
            pr = lambda op: V.conj().T @ op @ V
            cos = pr(s.cos_phi_operator())
            ops[n] = dict(H0=np.diag(ev), n=pr(s.n_operator()), cos=cos, kin=np.diag(ev) + J[n] * cos)  # 4Ec n^2
        I = {n: np.eye(self.sys[n].truncated_dim) for n in self.sys}
        k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
        G = p['g_n_GHz']
        self.A = (k3(ops['Q1']['kin'], I['TC'], I['Q2']) + k3(I['Q1'], ops['TC']['kin'], I['Q2'])
                  + k3(I['Q1'], I['TC'], ops['Q2']['kin'])
                  + G['Q1-TC'] * k3(ops['Q1']['n'], ops['TC']['n'], I['Q2'])
                  + G['Q2-TC'] * k3(I['Q1'], ops['TC']['n'], ops['Q2']['n'])
                  + G['Q1-Q2'] * k3(ops['Q1']['n'], I['TC'], ops['Q2']['n']))
        self.B = dict(Q1=-k3(ops['Q1']['cos'], I['TC'], I['Q2']), TC=-k3(I['Q1'], ops['TC']['cos'], I['Q2']),
                      Q2=-k3(I['Q1'], I['TC'], ops['Q2']['cos']))

    def EJ(self, n, flux):
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
                    labels[s] = k; used.add(k); break
        return w, v, labels

    def static(self, **fl):
        w, v, lab = self.dressed(**fl)
        E = lambda s: w[lab[s]] - w[lab[(0, 0, 0)]]
        return dict(f1=E((1, 0, 0)), f2=E((0, 0, 1)), fc=E((0, 1, 0)), zz=E((1, 0, 1)) - E((1, 0, 0)) - E((0, 0, 1)),
                    a1=E((2, 0, 0)) - 2 * E((1, 0, 0)), a2=E((0, 0, 2)) - 2 * E((0, 0, 1)))


def envelope(t, T, ramp):
    t = np.asarray(t, float)
    w = np.ones_like(t)
    r1, r2 = t < ramp, t > T - ramp
    w[r1] = 0.5 * (1 - np.cos(np.pi * t[r1] / ramp))
    w[r2] = 0.5 * (1 - np.cos(np.pi * (T - t[r2]) / ramp))
    return np.clip(w, 0, 1)


COMP = [(0, 0, 0), (0, 0, 1), (1, 0, 0), (1, 0, 1)]          # |Q1 Q2> = 00, 01, 10, 11


def evolve(m, idle, a_q, a_c, T, ramp, record=False):
    w0, v0, lab = m.dressed(fc=idle)
    tm = (np.arange(int(round(T / DT))) + 0.5) * DT
    env = envelope(tm, T, ramp)
    psi = np.stack([v0[:, lab[s]] for s in COMP], axis=1).astype(complex)
    rec = []
    for k, wk in enumerate(env):
        Hk = m.H(f1=a_q * wk, fc=idle + a_c * wk)
        ev, V = eigh(Hk)
        psi = V @ (np.exp(-2j * np.pi * ev * DT)[:, None] * (V.conj().T @ psi))
        if record and k % 10 == 0:
            rec.append((tm[k], np.abs(v0.conj().T @ psi[:, 3]) ** 2))
    c = (v0.conj().T @ psi) * np.exp(2j * np.pi * w0 * T)[:, None]       # idle rotating frame
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


def optimise(m, idle, T, ramp, a_q0, a_c_max, use_q=True):
    """Nelder-Mead on 1 - F. The coupler amplitude is bounded to [0, a_c_max] (a_c = a_c_max sin^2 u), so the coupler
    never comes closer than FC_MARGIN to the qubits; started from the best point of a coarse grid."""
    unpack = lambda x: ((x[0], a_c_max * np.sin(x[1]) ** 2) if use_q else (0.0, a_c_max * np.sin(x[0]) ** 2))

    def cost(x):
        aq, ac = unpack(x)
        return 1 - analyse_U(evolve(m, idle, aq, ac, T, ramp))['fidelity']
    us = np.linspace(0.35, 1.5, 6)
    grid = [[a_q0, u] for u in us] if use_q else [[u] for u in us]
    x0 = min(grid, key=cost)
    simplex = (np.array([x0, [x0[0] + 0.004, x0[1]], [x0[0], x0[1] + 0.08]]) if use_q
               else np.array([x0, [x0[0] + 0.08]]))
    sol = minimize(cost, x0, method='Nelder-Mead',
                   options=dict(xatol=1e-5, fatol=1e-7, maxiter=150, initial_simplex=simplex))
    aq, ac = unpack(sol.x)
    r = analyse_U(evolve(m, idle, aq, ac, T, ramp))
    return dict(T_ns=T, ramp_ns=ramp, a_q=float(aq), a_c=float(ac), **r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--params', required=True)
    ap.add_argument('--tag', default='hfss')
    ap.add_argument('--times', default='30,40,50,60,80')
    ap.add_argument('--ramp', type=float, default=8.0)
    ap.add_argument('--tc-times', default='150', help='coupler-only CZ gate times (ns), empty: skip')
    ap.add_argument('--date', default=time.strftime('%Y-%m-%d'))
    ap.add_argument('--retarget', default='', help='f01 targets Q1,Q2,TC (GHz): re-choose Lj for these with the given Ec')
    a = ap.parse_args()
    out = _relocate(os.path.join(SIM, 'results', a.date, f'cz_{a.tag}'))
    os.makedirs(out, exist_ok=True)
    tg = dict(zip(('Q1', 'Q2', 'TC'), map(float, a.retarget.split(',')))) if a.retarget else None
    p = load_params(a.params, tg)
    p['retarget_GHz'] = tg
    m = Model(p)
    print('parameters:', json.dumps({k: v for k, v in p.items() if k != 'file'}, default=float), flush=True)

    # ---- 1. static ----
    fl = np.linspace(0, 0.40, 81)
    st = [m.static(fc=f) for f in fl]
    zz = np.array([s['zz'] for s in st]) * 1e6
    fcs = np.array([s['fc'] for s in st])
    idle = 0.0
    sgn = [i for i in range(len(fl) - 1) if np.sign(zz[i]) != np.sign(zz[i + 1]) and fcs[i + 1] > st[0]['f1'] + 0.8]
    if sgn:
        i = sgn[0]
        idle = brentq(lambda f: m.static(fc=f)['zz'], fl[i], fl[i + 1])
    si = m.static(fc=idle)
    print(f"idle: coupler flux {idle:.4f}, f1 {si['f1']:.4f}, f2 {si['f2']:.4f}, fc {si['fc']:.4f} GHz, "
          f"alpha {si['a1']*1e3:.1f} / {si['a2']*1e3:.1f} MHz, ZZ {si['zz']*1e6:.2f} kHz", flush=True)
    # Q1 flux for the |11> <-> |20> condition (coupler idle): 2 f1 + alpha1 = f1 + f2  ->  f1 = f2 - alpha1.
    # Q1 moves DOWN by |Delta + alpha1| ~ 0.2 GHz and stays above Q2 (the |11> <-> |02> condition f1 = f2 + alpha2 would
    # make Q1 cross Q2 during the ramps: |01> <-> |10> swap error, F capped at 0.996 in run cz_hfss_v8 v1).
    target = si['f2'] - si['a1']
    fq = np.linspace(0, 0.3, 61)
    f1s = np.array([m.static(f1=x, fc=idle)['f1'] for x in fq])
    a_q0 = float(np.interp(target, f1s[::-1], fq[::-1]))
    print(f'Q1 flux for f1 = f2 - alpha1 ({target:.4f} GHz, |11> <-> |20>): {a_q0:.4f} Phi0', flush=True)

    # ---- 2. CZ ----
    def ac_max(f_min_c):            # largest coupler excursion with f_c >= f_min_c (qubits at idle flux)
        ok = [f for f, fc_ in zip(fl, fcs) if fc_ >= f_min_c and f >= idle]
        return max(ok) - idle
    acq = ac_max(target + FC_MARGIN)
    acc = ac_max(si['f1'] + FC_MARGIN)
    print(f'coupler bound: f_c >= {target + FC_MARGIN:.3f} GHz (Q1+TC gate, a_c <= {acq:.4f}), '
          f'>= {si["f1"] + FC_MARGIN:.3f} GHz (TC-only, a_c <= {acc:.4f})', flush=True)
    res, res_c = [], []
    for T in [float(x) for x in a.times.split(',')]:
        r = optimise(m, idle, T, a.ramp, a_q0, acq)
        sg = m.static(f1=r['a_q'], fc=idle + r['a_c'])
        r.update(f1_gate_GHz=sg['f1'], fc_gate_GHz=sg['fc'])
        res.append(r)
        print(f"[Q1+TC] T {T:5.1f} ns: a_q {r['a_q']:.4f} (f1 {sg['f1']:.3f}), a_c {r['a_c']:+.4f} (fc {sg['fc']:.3f}) -> "
              f"phase {r['cond_phase_deg']:.2f} deg, leakage {r['leakage']:.2e}, F {r['fidelity']:.6f}", flush=True)
    for T in [float(x) for x in a.tc_times.split(',') if x]:
        r = optimise(m, idle, T, a.ramp, 0.0, acc, use_q=False)
        sg = m.static(fc=idle + r['a_c'])
        r.update(fc_gate_GHz=sg['fc'], zz_gate_MHz=sg['zz'] * 1e3)
        res_c.append(r)
        print(f"[TC only] T {T:5.1f} ns: a_c {r['a_c']:+.4f} (fc {sg['fc']:.3f}, ZZ {sg['zz']*1e3:.2f} MHz) -> "
              f"phase {r['cond_phase_deg']:.2f}, leakage {r['leakage']:.2e}, F {r['fidelity']:.6f}", flush=True)
    best = max(res, key=lambda r: r['fidelity'])
    _, pops = evolve(m, idle, best['a_q'], best['a_c'], best['T_ns'], a.ramp, record=True)
    summary = dict(params=p, idle=dict(coupler_flux=idle, **si), Q1_flux_for_11_20=a_q0, Q1_gate_target_GHz=target, cz_q1_plus_coupler=res,
                   cz_coupler_only=res_c, best=best,
                   zz_vs_flux=dict(flux=fl.tolist(), fc_GHz=fcs.tolist(), zz_kHz=zz.tolist(),
                                   f1_GHz=[s['f1'] for s in st], f2_GHz=[s['f2'] for s in st]))
    json.dump(summary, open(_relocate(os.path.join(out, 'cz_summary.json')), 'w'), indent=1, default=float)

    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(fl, fcs, label='TC'); ax[0].plot(fl, [s['f1'] for s in st], label='Q1')
    ax[0].plot(fl, [s['f2'] for s in st], label='Q2'); ax[0].axvline(idle, ls=':', c='k')
    ax[0].set_xlabel('coupler flux ($\\Phi_0$), qubits at sweet spot'); ax[0].set_ylabel('dressed frequency (GHz)'); ax[0].legend()
    ax[1].semilogy(fl, np.abs(zz)); ax[1].axvline(idle, ls=':', c='k')
    ax[1].set_xlabel('coupler flux ($\\Phi_0$)'); ax[1].set_ylabel('|ZZ| / 2$\\pi$ (kHz)')
    fig.suptitle(f"Design6 static ZZ ({p['source']})"); fig.tight_layout()
    fig.savefig(_relocate(os.path.join(out, 'zz_vs_flux.png')), dpi=140); plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].semilogy([r['T_ns'] for r in res], [1 - r['fidelity'] for r in res], 'o-', label='1 - F (Q1 + TC pulse)')
    ax[0].semilogy([r['T_ns'] for r in res], [max(r['leakage'], 1e-9) for r in res], 's--', label='leakage (Q1 + TC)')
    ax[0].semilogy([r['T_ns'] for r in res_c], [1 - r['fidelity'] for r in res_c], '^-', label='1 - F (TC only)')
    ax[0].set_xlabel('gate time T (ns)'); ax[0].legend(); ax[0].set_title('CZ coherent error')
    for k, v in pops.items():
        if k != 't':
            ax[1].plot(pops['t'], v, label=f'|{k}>')
    ax[1].set_xlabel('t (ns)'); ax[1].set_ylabel('population, start |11>'); ax[1].legend()
    ax[1].set_title(f"best: T = {best['T_ns']:.0f} ns, F = {best['fidelity']:.5f}")
    fig.tight_layout(); fig.savefig(_relocate(os.path.join(out, 'cz_vs_time.png')), dpi=140); plt.close(fig)
    print('saved', out)


if __name__ == '__main__':
    main()
