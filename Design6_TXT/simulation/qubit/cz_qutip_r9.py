"""Independent QuTiP check of the Design6_TXT r9 Q1 + TC CZ (scqubits Hamiltonian, no EM rerun), as a new run.

Hamiltonian: identical inputs to cz_r9.py 'layout' set -- adopted r9 core Q3D (run 20261117_161936_q3d_core_r9_3b00f1f5,
pass 20, Ec and charge couplings g_n) and the drawn L_J in layout/design6_layout.py (read from the source text).
Pulse: the parent run's best sampled Q1 + TC point (run 20261117_181834_cz_r9_20f32abc), not re-optimised unless stated.

Checks
 1. coherent: QuTiP sesolve (adaptive ODE, continuous cosine-ramp envelope) vs the parent piecewise-constant propagator;
 2. truncation: same pulse with larger local level counts and charge cutoff (pulse re-optimised at the largest basis);
 3. pulse sensitivity: gate time and flux-amplitude detuning (parent propagator, fast);
 4. open system: QuTiP mesolve on a reduced 3 x 4 x 3 basis, Lindblad T1 (a_i) and pure dephasing (a_i^dag a_i) on all three transmons;
    the decoherence-added error (vs the same-basis run without collapse operators) is added to the 4 x 5 x 4 coherent error,
    constant rates (flux-noise increase away from the sweet spot NOT modelled). T1/T2 are ASSUMED scan values.
Average gate fidelity F = (sum_j Tr[U_cz P_j^dag U_cz^dag E(P_j)] + d^2) / (d^2 (d + 1)), d = 4, two-qubit Paulis P_j
embedded on the dressed idle computational states, idle rotating frame, single-qubit virtual Z taken from the coherent U.

Python: scqenv (scqubits 4.1, qutip 5.0).
"""
from argparse import ArgumentParser
from pathlib import Path
import itertools
import json
import re
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(Path(__file__).parent), str(ROOT / 'code/general')]
import cz_scqubits as C            # noqa: E402
from run_record import create_run, finish_run  # noqa: E402
import qutip as qt                 # noqa: E402

DESIGN = Path(__file__).resolve().parents[2]
Q3D = DESIGN / 'results/qubit/20261117_161936_q3d_core_r9_3b00f1f5'
PARENT = DESIGN / 'results/qubit/20261117_181834_cz_r9_20f32abc'
LAYOUT = DESIGN / 'layout/design6_layout.py'
TWO_PI = 2 * np.pi


def layout_lj():
    m = re.search(r'^LJ = dict\((.*?)\)', LAYOUT.read_text(encoding='utf-8'), re.M)
    return {k: float(v) for k, v in re.findall(r"(\w+)='([\d.]+)nH'", m.group(1))}


def params_file(folder):
    f = json.loads((Q3D / 'summary_out.json').read_text())['final']
    p = dict(source_run=Q3D.name, pass_=20, Ec_GHz=f['Ec_GHz'], C_sigma_fF=f['C_sigma_fF'], g_n_MHz=f['g_n_MHz'],
             Lj_drawn_nH=layout_lj())
    path = folder / 'params_r9_q3d.json'
    path.write_text(json.dumps(p, indent=1))
    return path


def model(pf, nq, nc, ncut):
    C.NQ, C.NC, C.NCUT = nq, nc, ncut
    return C.Model(C.load_params(str(pf)))


def idle_frame(m):
    w0, v0, lab = m.dressed(fc=0.0)
    comp = [lab[s] for s in C.COMP]
    return w0, v0, comp


def h_qobjevo(m, a_q, a_c, T, ramp):
    """H(t) in rad/ns: A + EJ1(t) B1 + EJc(t) Bc + EJ2 B2, envelope as in cz_scqubits (exact cosine ramps)."""
    dims = [list(m.dims), list(m.dims)]
    Q = lambda M: qt.Qobj(TWO_PI * M, dims=dims).to('CSR')    # sparse: keeps the mesolve Liouvillian small

    def w(t):                                                  # scalar copy of C.envelope
        if t < ramp:
            return 0.5 * (1 - np.cos(np.pi * t / ramp))
        if t > T - ramp:
            return 0.5 * (1 - np.cos(np.pi * max(T - t, 0.0) / ramp))
        return 1.0
    H0 = Q(m.A + m.EJ('Q2', 0.0) * m.B['Q2'])
    return qt.QobjEvo([H0, [Q(m.B['Q1']), lambda t: m.EJ('Q1', a_q * w(t))],
                       [Q(m.B['TC']), lambda t: m.EJ('TC', a_c * w(t))]])


OPEN_LEVELS = (3, 4, 15)   # (nq, nc, ncut) for the mesolve check
ODE = dict(atol=1e-11, rtol=1e-9, nsteps=10 ** 7, max_step=0.02)


def coherent_qutip(m, a_q, a_c, T, ramp):
    w0, v0, comp = idle_frame(m)
    H = h_qobjevo(m, a_q, a_c, T, ramp)
    U = np.zeros((4, 4), complex)
    for k, ci in enumerate(comp):
        psi0 = qt.Qobj(v0[:, ci], dims=[list(m.dims), [1] * 3])
        psi = qt.sesolve(H, psi0, [0, T], options=ODE).states[-1].full().ravel()
        c = (v0.conj().T @ psi) * np.exp(1j * TWO_PI * w0 * T)
        U[:, k] = c[comp]
    return U


def z_correction(U):
    ph = np.angle(np.diag(U))
    return np.diag(np.exp(-1j * np.array([ph[0], ph[1], ph[2], ph[1] + ph[2] - ph[0]])))


PAULI = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0])]
P2 = [np.kron(a, b) for a, b in itertools.product(PAULI, PAULI)]
UCZ = np.diag([1, 1, 1, -1]).astype(complex)


def collapse_ops(m, T1_us, T2_us):
    """Per-transmon T1 (sqrt(k) ladder) and pure dephasing (number op) in the local eigenbasis, rates in 1/ns."""
    if T1_us is None:
        return []
    ops = []
    Tphi = 1.0 / (1.0 / T2_us - 1.0 / (2 * T1_us))
    for i, n in enumerate(m.dims):
        a = np.diag(np.sqrt(np.arange(1, n)), 1)
        def emb(op, i=i):
            mats = [op if k == i else np.eye(d) for k, d in enumerate(m.dims)]
            return qt.Qobj(np.kron(np.kron(mats[0], mats[1]), mats[2]), dims=[list(m.dims), list(m.dims)]).to('CSR')
        ops.append(np.sqrt(1.0 / (T1_us * 1e3)) * emb(a))
        ops.append(np.sqrt(2.0 / (Tphi * 1e3)) * emb(a.T @ a))
    return ops


def open_fidelity(m, a_q, a_c, T, ramp, T1_us, T2_us, Zc):
    w0, v0, comp = idle_frame(m)
    H = h_qobjevo(m, a_q, a_c, T, ramp)
    cops = collapse_ops(m, T1_us, T2_us)
    Piso = v0[:, comp]                             # 80 x 4 isometry onto dressed comp states
    D = np.exp(1j * TWO_PI * w0 * T)
    # E on the 16 matrix units |i><j| (linear), then assemble Paulis
    E = {}
    for i in range(4):
        for j in range(4):
            rho0 = qt.Qobj(np.outer(Piso[:, i], Piso[:, j].conj()), dims=[list(m.dims), list(m.dims)])
            r = qt.mesolve(H, rho0, [0, T], c_ops=cops, options=ODE).states[-1].full()
            rd = (D[:, None] * (v0.conj().T @ r @ v0) * D.conj()[None, :])[np.ix_(comp, comp)]
            E[i, j] = Zc @ rd @ Zc.conj().T
    tot = 0.0
    for P in P2:
        out = sum(P[i, j] * E[i, j] for i in range(4) for j in range(4))
        tot += np.trace(UCZ @ P.conj().T @ UCZ.conj().T @ out).real
    return (tot + 16) / (16 * 5)


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--t1-t2-us', default='200:200,100:100,50:50,30:30,20:20',
                    help='assumed T1:T2 pairs (us) for all three transmons')
    ap.add_argument('--skip-open', action='store_true')
    a = ap.parse_args()
    parent = json.loads((PARENT / 'cz_summary.json').read_text())['layout']
    best = max(parent['q1_plus_tc']['cz'], key=lambda r: r['fidelity'])
    pairs = [tuple(None if x == 'inf' else float(x) for x in s.split(':')) for s in a.t1_t2_us.split(',')]
    folder = create_run(DESIGN.name, 'qubit', 'cz_qutip_r9', entrypoint=__file__,
                        inputs=[Q3D / 'summary_out.json', LAYOUT, PARENT / 'cz_summary.json'],
                        upstream=[str(Q3D), str(PARENT)],
                        parameters=dict(EM_rerun=False, parent_pulse=best, t1_t2_us_assumed=a.t1_t2_us,
                                        decoherence_model='constant Lindblad T1 + pure dephasing on Q1, TC, Q2; '
                                                          'no flux-noise increase off sweet spot; assumed values',
                                        ode=ODE, qutip=qt.__version__))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        pf = params_file(folder)
        T, ramp, aq, ac = best['T_ns'], best['ramp_ns'], best['a_q'], best['a_c']
        out = dict(parent_run=PARENT.name, hamiltonian_run=Q3D.name, pulse=best)

        # 1. coherent: parent propagator vs QuTiP ODE (4 x 5 x 4, ncut 15)
        m = model(pf, 4, 5, 15)
        U_pc = C.evolve(m, 0.0, aq, ac, T, ramp)
        U_qt = coherent_qutip(m, aq, ac, T, ramp)
        out['coherent'] = dict(piecewise=C.analyse_U(U_pc), qutip_sesolve=C.analyse_U(U_qt),
                               max_abs_U_diff=float(np.max(np.abs(U_pc - U_qt))))
        print('[coherent]', out['coherent'], flush=True)
        Zc = z_correction(U_qt)

        # 2. truncation (parent propagator, same pulse; re-optimised on the largest basis)
        trunc = []
        for nq, nc, ncut in ((4, 5, 15), (5, 6, 15), (6, 7, 20)):
            mm = model(pf, nq, nc, ncut)
            st = mm.static(fc=0.0)
            r = C.analyse_U(C.evolve(mm, 0.0, aq, ac, T, ramp))
            trunc.append(dict(levels=[nq, nc, nq], ncut=ncut, idle_zz_kHz=st['zz'] * 1e6, f1=st['f1'], f2=st['f2'],
                              fc=st['fc'], same_pulse=r))
            print(f'[trunc] {nq}x{nc}x{nq} ncut {ncut}: ZZ {st["zz"]*1e6:+.3f} kHz, {r}', flush=True)
        mm = model(pf, 6, 7, 20)
        ac_max = max(r['a_c'] for r in parent['q1_plus_tc']['cz'])   # parent coupler bound (f_c >= target + 0.4 GHz)
        reopt = C.optimise(mm, 0.0, T, ramp, aq, ac_max)
        trunc.append(dict(levels=[6, 7, 6], ncut=20, reoptimised=reopt))
        print('[trunc reopt 6x7x6]', reopt, flush=True)
        out['truncation'] = trunc

        # 3. sensitivity (4 x 5 x 4)
        m = model(pf, 4, 5, 15)
        sens = dict(gate_time=[], a_q=[], a_c=[])
        for dT in (-2, -1, -0.5, 0, 0.5, 1, 2):
            sens['gate_time'].append(dict(dT_ns=dT, **C.analyse_U(C.evolve(m, 0.0, aq, ac, T + dT, ramp))))
        for rel in (-1e-2, -3e-3, -1e-3, 0, 1e-3, 3e-3, 1e-2):
            sens['a_q'].append(dict(rel=rel, **C.analyse_U(C.evolve(m, 0.0, aq * (1 + rel), ac, T, ramp))))
            sens['a_c'].append(dict(rel=rel, **C.analyse_U(C.evolve(m, 0.0, aq, ac * (1 + rel), T, ramp))))
        out['sensitivity'] = sens
        for k, v in sens.items():
            print(f'[sens {k}]', [(list(x.values())[0], round(x['fidelity'], 6)) for x in v], flush=True)

        # 4. open system on a reduced 3 x 4 x 3 basis (Liouvillian 1296^2); the decoherence-added error is taken
        #    relative to the same-basis mesolve without collapse operators and added to the 4 x 5 x 4 coherent error.
        mo = model(pf, *OPEN_LEVELS)
        Zo = z_correction(coherent_qutip(mo, aq, ac, T, ramp))
        F_coh_full = out['coherent']['qutip_sesolve']['fidelity']
        opn = []
        if not a.skip_open:
            F_ref = None
            for T1, T2 in [(None, None)] + [x for x in pairs if x[0] is not None]:
                t0 = time.time()
                F = float(open_fidelity(mo, aq, ac, T, ramp, T1, T2, Zo))
                F_ref = F if T1 is None else F_ref
                added = F_ref - F
                opn.append(dict(T1_us=T1, T2_us=T2, levels=list(mo.dims), fidelity_reduced_basis=F,
                                decoherence_error=added, fidelity_estimate=F_coh_full - added,
                                minutes=(time.time() - t0) / 60))
                print(f'[open] T1 {T1} T2 {T2} us: F(3x4x3) {F:.6f}, added {added:.2e}, '
                      f'F est {F_coh_full - added:.6f} ({(time.time()-t0)/60:.1f} min)', flush=True)
        out['open_system'] = opn

        (folder / 'cz_qutip_summary.json').write_text(json.dumps(out, indent=1, default=float))
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(14, 4))
        ax[0].semilogy([x['dT_ns'] for x in sens['gate_time']], [1 - x['fidelity'] for x in sens['gate_time']], 'o-')
        ax[0].set_xlabel('gate-time offset (ns)'); ax[0].set_ylabel('1 - F (coherent)'); ax[0].grid(alpha=.3)
        for k in ('a_q', 'a_c'):
            ax[1].semilogy([x['rel'] * 100 for x in sens[k]], [1 - x['fidelity'] for x in sens[k]], 'o-', label=k)
        ax[1].set_xlabel('flux amplitude error (%)'); ax[1].legend(); ax[1].grid(alpha=.3)
        oo = [x for x in opn if x['T1_us']]
        if oo:
            ax[2].loglog([x['T1_us'] for x in oo], [1 - x['fidelity_estimate'] for x in oo], 'o-', label='coherent + mesolve decoherence')
            ax[2].axhline(1 - out['coherent']['qutip_sesolve']['fidelity'], ls=':', c='k', label='coherent')
            ax[2].set_xlabel('T1 = T2 (us, assumed)'); ax[2].set_ylabel('1 - F'); ax[2].legend(); ax[2].grid(alpha=.3)
        fig.suptitle(f'Design6 r9 Q1+TC CZ, T = {T:.0f} ns (QuTiP check)'); fig.tight_layout()
        fig.savefig(folder / 'cz_qutip_r9.png', dpi=130); plt.close(fig)
        finish_run(folder, status='completed',
                   summary=dict(coherent=out['coherent'], open_system=opn,
                                truncation_reopt=trunc[-1]['reoptimised']))
    except Exception:
        import traceback
        (folder / 'error.log').write_text(traceback.format_exc())
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()}); raise


if __name__ == '__main__':
    main()
