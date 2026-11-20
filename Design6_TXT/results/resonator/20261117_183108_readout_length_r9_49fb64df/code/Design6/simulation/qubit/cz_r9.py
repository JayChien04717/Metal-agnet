"""Design6 CZ on the r9 core Q3D Hamiltonian (scqubits, no EM rerun), written to a new run folder.

Hamiltonian source: run 20261117_161936_q3d_core_r9_3b00f1f5 (pass 20, Delta 0.275 %, not converged to 0.1 %; the user
accepts ~3 % error on Ec / g for design decisions). Ec and g_n (charge units) from that Q3D matrix.
Junctions: the fab inputs in design6_layout.LJ (read from the source text, not retargeted) for the 'layout' set; the
'strad_*' sets re-solve the junctions for hypothetical straddling targets (proposals only, EJ is a user input).

Schemes (functions reused from cz_scqubits / cz_coupler_flux, same propagation and fidelity definition):
  * Q1 + TC pulse (the rev-A design scheme): Q1 to the |11>-|20> condition, coupler towards the qubits.
  * coupler-only: coupler pulsed from its ZZ-off (idle) bias, qubits at their sweet spots.
"""
from argparse import ArgumentParser
from pathlib import Path
import json
import re
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(Path(__file__).parent), str(ROOT / 'code/general')]
import cz_scqubits as C            # noqa: E402
import cz_coupler_flux as CF       # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

Q3D = ROOT / 'Design6/results/qubit/20261117_161936_q3d_core_r9_3b00f1f5'
LAYOUT = ROOT / 'Design6/layout/design6_layout.py'
SETS = {
    'layout': dict(retarget=None, g_qc=None, g12=None, note='r9 Q3D Ec/g, layout LJ as drawn'),
    'strad_native_g': dict(retarget=dict(Q1=4.55, Q2=4.40, TC=6.80), g_qc=None, g12=None,
                           note='straddling targets, r9 couplings unchanged'),
    'strad_g61_g12_4': dict(retarget=dict(Q1=4.55, Q2=4.40, TC=6.80), g_qc=61.0, g12=4.0,
                            note='straddling proposal of the 2026-11-17 exploration, on r9 Ec'),
}


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


def build(path, s):
    p = C.load_params(str(path), s['retarget'])
    p['source'] = 'Q3D r9 core (pass 20)'
    if s['g_qc'] is not None:
        for k in ('Q1-TC', 'Q2-TC'):
            p['g_n_GHz'][k] = s['g_qc'] / 1e3
    if s['g12'] is not None:
        p['g_n_GHz']['Q1-Q2'] = s['g12'] / 1e3
    return C.Model(p), p


def q1tc(m, times, ramp, fl):
    """Q1 + TC scheme exactly as cz_scqubits.main (idle at the ZZ zero, else the coupler sweet spot)."""
    st = [m.static(fc=f) for f in fl]
    zz = np.array([s['zz'] for s in st]); fcs = np.array([s['fc'] for s in st])
    idle = 0.0
    sgn = [i for i in range(len(fl) - 1) if np.sign(zz[i]) != np.sign(zz[i + 1]) and fcs[i + 1] > st[0]['f1'] + 0.8]
    if sgn:
        from scipy.optimize import brentq
        idle = brentq(lambda f: m.static(fc=f)['zz'], fl[sgn[0]], fl[sgn[0] + 1])
    si = m.static(fc=idle)
    target = si['f2'] - si['a1']
    fq = np.linspace(0, 0.3, 61)
    f1s = np.array([m.static(f1=x, fc=idle)['f1'] for x in fq])
    a_q0 = float(np.interp(target, f1s[::-1], fq[::-1]))
    acq = max(f for f, c in zip(fl, fcs) if c >= target + C.FC_MARGIN and f >= idle) - idle
    res = []
    for T in times:
        r = C.optimise(m, idle, T, ramp, a_q0, acq)
        sg = m.static(f1=r['a_q'], fc=idle + r['a_c'])
        r.update(f1_gate_GHz=sg['f1'], fc_gate_GHz=sg['fc'])
        res.append(r)
        print(f"  [Q1+TC] T {T:4.0f} ns: f1 {sg['f1']:.3f} fc {sg['fc']:.3f} phase {r['cond_phase_deg']:7.2f} "
              f"leak {r['leakage']:.1e} F {r['fidelity']:.5f}", flush=True)
    return dict(idle=dict(coupler_flux=idle, **si), Q1_gate_target_GHz=target, a_q0=a_q0, cz=res)


def tc_only(m, times, ramp, fl, idle_margin=0.8, gate_margin=0.3):
    idle, how, st = CF.idle_point(m, idle_margin, fl)
    si = m.static(fc=idle); fcs = np.array([x['fc'] for x in st])
    ac_max = max(f for f, c in zip(fl, fcs) if c >= si['f1'] + gate_margin and f >= idle) - idle
    son = m.static(fc=idle + ac_max)
    print(f"  [TC only] idle ({how}) fc {si['fc']:.3f} ZZ {si['zz']*1e6:+.2f} kHz | bound fc {son['fc']:.3f} "
          f"ZZ_on {son['zz']*1e3:+.3f} MHz", flush=True)
    res = []
    for T in times:
        if T <= 2 * ramp:
            continue
        r, _ = CF.cz_at(m, idle, T, ramp, ac_max)
        r['fc_on_GHz'] = m.static(fc=idle + r['a_c'])['fc']
        r.pop('pop_from_11_final', None)
        res.append(r)
        print(f"  [TC only] T {T:4.0f} ns: fc_on {r['fc_on_GHz']:.3f}{' (bound)' if r['at_bound'] else ''} "
              f"phase {r['cond_phase_deg']:7.2f} leak {r['leakage']:.1e} F {r['fidelity']:.5f}", flush=True)
    return dict(idle=dict(coupler_flux=idle, method=how, **si), on_bound=dict(coupler_flux=idle + ac_max, **son),
                zz_vs_flux=dict(fc_GHz=fcs.tolist(), zz_kHz=[x['zz'] * 1e6 for x in st]), cz=res)


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--sets', default=','.join(SETS))
    ap.add_argument('--q1tc-times', default='40,50,60,80')
    ap.add_argument('--q1tc-ramp', type=float, default=8.0)
    ap.add_argument('--tc-times', default='100,130,160,200')
    ap.add_argument('--tc-times-layout', default='200,400,700')
    ap.add_argument('--tc-ramp', type=float, default=20.0)
    a = ap.parse_args()
    names = a.sets.split(',')
    folder = create_run('Design6', 'qubit', 'cz_r9', entrypoint=__file__,
                        inputs=[Q3D / 'summary_out.json', LAYOUT], upstream=[str(Q3D)],
                        parameters=dict(EM_rerun=False, sets={n: SETS[n] for n in names}, q1tc_times_ns=a.q1tc_times,
                                        q1tc_ramp_ns=a.q1tc_ramp, tc_times_ns=a.tc_times,
                                        tc_times_layout_ns=a.tc_times_layout, tc_ramp_ns=a.tc_ramp,
                                        levels=dict(Q=C.NQ, TC=C.NC, ncut=C.NCUT), dt_ns=C.DT, dissipation='none',
                                        q1tc_coupler_margin_GHz=C.FC_MARGIN,
                                        accuracy='Q3D not converged (0.275 %); ~3 % on Ec/g accepted by the user'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        pf = params_file(folder)
        fl = np.linspace(0, 0.42, 85)
        out = {}
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
        for n in names:
            s = SETS[n]; m, p = build(pf, s); t0 = time.time()
            print(f"[{n}] Lj {p['Lj_nH']}", flush=True)
            r1 = q1tc(m, [float(x) for x in a.q1tc_times.split(',')], a.q1tc_ramp, fl)
            print(f"  idle f1 {r1['idle']['f1']:.4f} f2 {r1['idle']['f2']:.4f} fc {r1['idle']['fc']:.4f} GHz "
                  f"alpha {r1['idle']['a1']*1e3:.0f}/{r1['idle']['a2']*1e3:.0f} MHz ZZ {r1['idle']['zz']*1e6:+.2f} kHz",
                  flush=True)
            tt = a.tc_times_layout if n == 'layout' else a.tc_times
            r2 = tc_only(m, [float(x) for x in tt.split(',')], a.tc_ramp, fl)
            out[n] = dict(set=s, Lj_nH=p['Lj_nH'], EJsum_GHz=p['EJsum_GHz'], d=p['d'], Ec_GHz=p['Ec_GHz'],
                          g_n_MHz={k: v * 1e3 for k, v in p['g_n_GHz'].items()}, q1_plus_tc=r1, coupler_only=r2,
                          minutes=(time.time() - t0) / 60)
            ln, = ax[0].semilogy(r2['zz_vs_flux']['fc_GHz'], np.abs(r2['zz_vs_flux']['zz_kHz']), label=n)
            ax[1].semilogy([r['T_ns'] for r in r1['cz']], [1 - r['fidelity'] for r in r1['cz']], 'o-',
                           color=ln.get_color(), label=f'{n}: Q1+TC')
            ax[1].semilogy([r['T_ns'] for r in r2['cz']], [1 - r['fidelity'] for r in r2['cz']], 's--',
                           color=ln.get_color(), label=f'{n}: TC only')
        ax[0].set_xlabel('coupler f (GHz), qubits at sweet spots'); ax[0].set_ylabel('|ZZ| (kHz)')
        ax[0].invert_xaxis(); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3); ax[0].set_title('static ZZ, r9 Q3D')
        ax[1].set_xlabel('gate time (ns)'); ax[1].set_ylabel('1 - F (coherent)'); ax[1].legend(fontsize=7)
        ax[1].grid(alpha=.3); ax[1].set_title('CZ on r9 Q3D Hamiltonian')
        fig.tight_layout(); fig.savefig(folder / 'cz_r9.png', dpi=130); plt.close(fig)
        (folder / 'cz_summary.json').write_text(json.dumps(out, indent=1, default=float))
        best = lambda rs: max(rs, key=lambda r: r['fidelity']) if rs else None
        finish_run(folder, status='completed',
                   summary={n: dict(idle_q1tc=v['q1_plus_tc']['idle'], best_q1tc=best(v['q1_plus_tc']['cz']),
                                    idle_tc=v['coupler_only']['idle'], on_bound_tc=v['coupler_only']['on_bound'],
                                    best_tc=best(v['coupler_only']['cz'])) for n, v in out.items()})
    except Exception:
        import traceback
        (folder / 'error.log').write_text(traceback.format_exc())
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()}); raise


if __name__ == '__main__':
    main()
