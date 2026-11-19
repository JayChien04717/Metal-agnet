"""Coupler-flux-only CZ for Design6 parameter exploration (scqubits, no EM rerun).

Gate: the qubits stay at their sweet spots; only the coupler flux is pulsed from the ZZ-off (idle) bias towards the
qubits and back (flat top, cosine ramps). The conditional phase comes from the ZZ accumulated at the on point.

Hamiltonian and propagation: cz_scqubits.Model / evolve / analyse_U (charge couplings g n_i n_j, 4 x 5 x 4 levels,
exact piecewise exponentials dt = 0.05 ns, dressed computational basis at idle, virtual Z removed, coherent only).

Parameter sets are HYPOTHETICAL design targets built on the HFSS v8 core fit (revision-A geometry, needs_revalidation):
  * Ec from that fit (assumed unchanged when g is scaled; a larger C_qc lowers Ec by a few %),
  * qubit and coupler sweet-spot f01 targets -> junctions re-solved with scqubits (EJ is a user design input;
    nothing is written back to the layout or the junction table),
  * g_qc = fit value x g_scale, g_12 set explicitly (MHz, charge units).
Idle coupler bias: the ZZ zero closest to the coupler sweet spot if one exists with f_c >= f_Q1 + idle_margin,
otherwise the bias of minimum |ZZ| in that range. Coupler bound during the gate: f_c >= f_Q1 + gate_margin.
"""
from argparse import ArgumentParser
from pathlib import Path
import json
import sys
import time
import numpy as np
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(Path(__file__).parent), str(ROOT / 'code/general')]
import cz_scqubits as C            # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

FIT = ROOT / 'Design6/results/qubit/2026-11-14/hfss_core_v8/fit.json'
SETS = {
    'baseline': dict(f1=4.80, f2=4.40, fc=6.30, g_scale=1.0, g12=None,
                     note='current Hamiltonian targets (Q1+TC CZ design), coupler-only operation'),
    'strad_g61_g12_4_fc68': dict(f1=4.55, f2=4.40, fc=6.80, g_scale=1.3, g12=4.0,
                                 note='straddling candidate from the static scan'),
    'strad_g61_g12_4_fc65': dict(f1=4.55, f2=4.40, fc=6.50, g_scale=1.3, g12=4.0,
                                 note='same, coupler f_max 6.5 GHz (more readout margin)'),
}


def build(s):
    p = C.load_params(str(FIT), dict(Q1=s['f1'], Q2=s['f2'], TC=s['fc']))
    p['g_n_GHz'] = {k: v * (s['g_scale'] if 'TC' in k else 1.0) for k, v in p['g_n_GHz'].items()}
    if s['g12'] is not None:
        p['g_n_GHz']['Q1-Q2'] = s['g12'] / 1e3
    return C.Model(p), p


def idle_point(m, idle_margin, fl):
    st = [m.static(fc=f) for f in fl]
    zz = np.array([x['zz'] for x in st]); fcs = np.array([x['fc'] for x in st]); f1 = st[0]['f1']
    ok = fcs >= f1 + idle_margin
    zeros = [i for i in range(len(fl) - 1) if ok[i + 1] and np.sign(zz[i]) != np.sign(zz[i + 1])]
    if zeros:
        i = zeros[0]
        return brentq(lambda f: m.static(fc=f)['zz'], fl[i], fl[i + 1]), 'zz_zero', st
    i = int(np.argmin(np.where(ok, np.abs(zz), np.inf)))
    return float(fl[i]), 'min_abs_zz', st


def cz_at(m, idle, T, ramp, ac_max):
    """Coupler amplitude for a 180 deg conditional phase (1-D bounded search on 1 - F)."""
    cost = lambda a: 1 - C.analyse_U(C.evolve(m, idle, 0.0, a, T, ramp))['fidelity']
    grid = np.linspace(0.15, 1.0, 8) * ac_max
    a0 = min(grid, key=cost)
    lo, hi = max(0.0, a0 - 0.15 * ac_max), min(ac_max, a0 + 0.15 * ac_max)
    sol = minimize_scalar(cost, bounds=(lo, hi), method='bounded', options=dict(xatol=1e-5))
    U, pops = C.evolve(m, idle, 0.0, sol.x, T, ramp, record=True)
    r = C.analyse_U(U)
    final = {k: float(v[-1]) for k, v in pops.items() if k != 't'}
    peak = {k: float(np.max(v)) for k, v in pops.items() if k != 't'}
    return dict(T_ns=T, ramp_ns=ramp, a_c=float(sol.x), at_bound=bool(sol.x > 0.999 * ac_max), **r,
                pop_from_11_final=final, pop_from_11_peak=peak), pops


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--sets', default=','.join(SETS))
    ap.add_argument('--times', default='60,80,100,130,160')
    ap.add_argument('--baseline-times', default='200,400,700')
    ap.add_argument('--ramps', default='10,20')
    ap.add_argument('--idle-margin', type=float, default=0.8)
    ap.add_argument('--gate-margin', type=float, default=0.3)
    a = ap.parse_args()
    names = a.sets.split(',')
    folder = create_run('Design6', 'qubit', 'cz_coupler_flux', entrypoint=__file__, inputs=[FIT],
                        upstream=[str(FIT.parent)],
                        parameters=dict(EM_rerun=False, sets={n: SETS[n] for n in names}, times_ns=a.times,
                                        baseline_times_ns=a.baseline_times, ramps_ns=a.ramps,
                                        idle_margin_GHz=a.idle_margin, gate_margin_GHz=a.gate_margin,
                                        levels=dict(Q=C.NQ, TC=C.NC, ncut=C.NCUT), dt_ns=C.DT, dissipation='none',
                                        status_of_inputs='hypothetical targets on the revision-A HFSS v8 fit '
                                                         '(needs_revalidation); junctions re-solved, not adopted'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fl = np.linspace(0, 0.42, 85)
        out, figs = {}, {}
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
        for name in names:
            s = SETS[name]; m, p = build(s); t0 = time.time()
            idle, how, st = idle_point(m, a.idle_margin, fl)
            si = m.static(fc=idle)
            fcs = np.array([x['fc'] for x in st]); f1 = si['f1']
            ok = [f for f, c in zip(fl, fcs) if c >= f1 + a.gate_margin and f >= idle]
            ac_max = max(ok) - idle
            son = m.static(fc=idle + ac_max)
            print(f"[{name}] idle flux {idle:.4f} ({how}): f1 {si['f1']:.4f} f2 {si['f2']:.4f} fc {si['fc']:.4f} GHz, "
                  f"ZZ {si['zz']*1e6:+.2f} kHz | bound fc {son['fc']:.3f} GHz ZZ_on {son['zz']*1e3:+.3f} MHz", flush=True)
            times = a.baseline_times if name == 'baseline' else a.times
            res, best_pops = [], None
            for ramp in map(float, a.ramps.split(',')):
                for T in map(float, times.split(',')):
                    if T <= 2 * ramp:
                        continue
                    r, pops = cz_at(m, idle, T, ramp, ac_max)
                    r['fc_on_GHz'] = m.static(fc=idle + r['a_c'])['fc']
                    res.append(r)
                    if best_pops is None or r['fidelity'] >= max(x['fidelity'] for x in res):
                        best_pops = (r, pops)
                    print(f"  ramp {ramp:4.0f} T {T:5.0f} ns: fc_on {r['fc_on_GHz']:.3f}{' (bound)' if r['at_bound'] else ''} "
                          f"phase {r['cond_phase_deg']:7.2f} leak {r['leakage']:.1e} F {r['fidelity']:.5f} "
                          f"| peak |20> {r['pop_from_11_peak']['20']:.3f} |02> {r['pop_from_11_peak']['02']:.3f}", flush=True)
            out[name] = dict(set=s, Lj_nH=p['Lj_nH'], EJsum_GHz=p['EJsum_GHz'], d=p['d'], Ec_GHz=p['Ec_GHz'],
                             g_n_MHz={k: v * 1e3 for k, v in p['g_n_GHz'].items()},
                             idle=dict(coupler_flux=idle, method=how, **si),
                             on_bound=dict(coupler_flux=idle + ac_max, **son), cz=res,
                             best=best_pops[0] if best_pops else None,
                             zz_vs_flux=dict(flux=fl.tolist(), fc_GHz=fcs.tolist(), zz_kHz=[x['zz'] * 1e6 for x in st]),
                             minutes=(time.time() - t0) / 60)
            ln, = ax[0].semilogy(fcs, np.abs([x['zz'] * 1e6 for x in st]), label=name)
            ax[0].axvline(si['fc'], color=ln.get_color(), ls=':')
            for ramp in map(float, a.ramps.split(',')):
                rr = [r for r in res if r['ramp_ns'] == ramp]
                ax[1].semilogy([r['T_ns'] for r in rr], [1 - r['fidelity'] for r in rr], 'o-' if ramp == 10 else 's--',
                               color=ln.get_color(), label=f'{name}, ramp {ramp:.0f} ns')
            if name != 'baseline' and best_pops:
                r, pops = best_pops
                for k in ('11', '20', '02', '1,c,0', '0,c,1'):
                    ax[2].plot(pops['t'], pops[k], label=f'{name[:14]} |{k}>' if k != '11' else f'{name[:14]} |11>')
        ax[0].set_xlabel('coupler f (GHz), qubits at sweet spots'); ax[0].set_ylabel('|ZZ| (kHz)'); ax[0].invert_xaxis()
        ax[0].set_title('static ZZ (dotted: idle bias)'); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)
        ax[1].set_xlabel('gate time T (ns)'); ax[1].set_ylabel('1 - F (coherent)'); ax[1].set_title('coupler-only CZ')
        ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
        ax[2].set_xlabel('t (ns)'); ax[2].set_ylabel('population from |11>'); ax[2].set_title('best gate per set')
        ax[2].legend(fontsize=6); ax[2].grid(alpha=.3)
        fig.tight_layout(); fig.savefig(folder / 'cz_coupler_flux.png', dpi=130); plt.close(fig)
        (folder / 'cz_summary.json').write_text(json.dumps(out, indent=1, default=float))
        finish_run(folder, status='completed', summary={n: dict(idle=v['idle'], on_bound=v['on_bound'], best=v['best'])
                                                         for n, v in out.items()})
    except Exception:
        import traceback
        (folder / 'error.log').write_text(traceback.format_exc())
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()}); raise


if __name__ == '__main__':
    main()
