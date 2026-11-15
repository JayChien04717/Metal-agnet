"""HFSS eigenmode of the Design6 core (Q1, TC, Q2 with claws), coupler junction inductance swept.

Model: build_design(channels=False, control_lines=False): Q1 / Q2 Xmon (asymmetric SQUID: two lumped JJs Lj_Q*a/b,
the linear model uses their parallel inductance),
TC claw coupler (two lumped SQUID JJs, Lj_TC each -> L_c = Lj_TC / 2), readout claws ending in their 40 um CPW stub
(open). Box = components + buffer (qiskit box_plus_buffer). Eigenmode, 3 modes from 3 GHz, first order,
6 passes max (skill.md section 12), 10 um length mesh on every qubit / coupler / claw electrode.

Each Lj_TC of the sweep is a variation of ONE HFSS design (project variable Lj_TC), solved in turn.

Linear circuit fit (the HFSS-derived Hamiltonian, used by simulation/cz):
    C phi'' + L^-1 phi = 0  ->  omega^2 = eig(Cinv @ diag(1/L_Q1, 1/L_c, 1/L_Q2))
6 free parameters = the symmetric Cinv; fitted to all (3 modes x N Lj_TC) HFSS frequencies.
Ec_i = e^2 Cinv_ii / 2h,  g_ij(n) = 4 e^2 Cinv_ij / h  (H_int = g n_i n_j).

Usage:  python hfss_core.py [--lj 9,11.3,14,17,21] [--tag base]     (nH, each SQUID JJ)
        python hfss_core.py --fit-only --tag base
Output: simulation/results/<date>/hfss_core_<tag>/{modes.json, convergence_<Lj>.csv, fit.json}
"""
import argparse, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np

SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DESIGN = os.path.dirname(SIM); REPO = os.path.dirname(DESIGN)
sys.path[:0] = [os.path.join(DESIGN, 'pattern design'), os.path.join(REPO, 'Design5', 'simulation')]
e, h = 1.602176634e-19, 6.62607015e-34
PHI0 = h / (2 * e)


def ej_ghz(L_nH):
    return (PHI0 / 2 / np.pi) ** 2 / (L_nH * 1e-9) / h / 1e9


def run(args, out):
    import design6_layout as L
    from qiskit_metal import Dict
    from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim
    import save_aedt
    d = L.build_design(channels=False, control_lines=False,
                       lj=dict(Q1a='Lj_Q1a', Q1b='Lj_Q1b', Q2a='Lj_Q2a', Q2b='Lj_Q2b', TC='Lj_TC'),
                       **{**json.loads(args.overrides), 'coupler': dict(json.loads(args.overrides).get('coupler', {}), flux_as_ground='True', squid_contact_overlap='0um', squid_shift_toward_stem='0um')})
    comps = ['Q1', 'Q2', 'TC']
    lj0 = {k: float(v.replace('nH', '')) for k, v in L.LJ.items()}
    lj_list = [float(x) for x in args.lj.split(',')]
    save_aedt.use_project(os.path.join(SIM, 'aedt', 'Design6.aedt'))

    metal = []
    id2n = {c.id: n for n, c in d.components.items()}
    for kind in ('poly', 'path'):
        for r in d.qgeometry.tables[kind].itertuples():
            if not bool(r.subtract) and id2n[int(r.component)] in comps:
                metal.append(f'{r.name}_{id2n[int(r.component)]}')

    class Eig(EigenmodeSim):
        def _analyze(self):
            self.renderer.modeler.mesh_length('qubit_mesh', metal, MaxLength=args.mesh)
            print(f'mesh {args.mesh} on {len(metal)} objects', flush=True)
            super()._analyze()

    eig = Eig(d, 'hfss')
    vars0 = Dict(Lj_Q1a=f"{lj0['Q1a']}nH", Lj_Q1b=f"{lj0['Q1b']}nH", Lj_Q2a=f"{lj0['Q2a']}nH",
                 Lj_Q2b=f"{lj0['Q2b']}nH", Lj_TC=f'{lj_list[0]}nH')
    for k, v in dict(min_freq_ghz=3.0, n_modes=3, max_delta_f=0.1, max_passes=args.maxp, min_passes=2,
                     min_converged=1, pct_refinement=30, basis_order=1, vars=vars0).items():
        eig.setup[k] = v
    t0 = time.time()
    dname = save_aedt.design_name(f'core_{args.tag}')
    eig.run(name=dname, components=comps, open_terminations=[('Q1', 'readout'), ('Q2', 'readout')])
    pinfo = eig.renderer.pinfo
    setup = pinfo.setup
    rows = []
    for k, lj in enumerate(lj_list):
        if k > 0:
            pinfo.design.set_variable('Lj_TC', f'{lj}nH')
            setup.analyze()
        f, im = setup.get_solutions().eigenmodes()
        conv = setup.get_convergence()[0]
        conv.to_csv(os.path.join(out, f'convergence_LjTC{lj:g}nH.csv'))
        row = dict(Lj_TC_nH=lj, L_c_nH=lj / 2, f_GHz=sorted(float(x) for x in f), passes=int(len(conv)),
                   last_delta_f_pct=float(conv.iloc[-1, -1]), minutes=round((time.time() - t0) / 60, 1))
        rows.append(row)
        print(f"Lj_TC {lj:5.2f} nH: modes {', '.join(f'{x:.4f}' for x in row['f_GHz'])} GHz, "
              f"{row['passes']} passes, dF {row['last_delta_f_pct']:.3f} %, {row['minutes']} min", flush=True)
        json.dump(dict(design=dname + '_hfss', Lj_nH=lj0, mesh=args.mesh, max_passes=args.maxp,
                       overrides=json.loads(args.overrides), sweep=rows),
                  open(os.path.join(out, 'modes.json'), 'w'), indent=1)
    try:
        save_aedt.snapshot(out, f'core_{args.tag}')
    except Exception as ex:
        print('snapshot failed:', repr(ex)[:120])
    return rows


def model_freqs(cinv6, Lq1, Lc, Lq2):
    a, b, c, d_, e_, f_ = cinv6               # Cinv entries (1/fF): [[a, d, e], [d, b, f], [e, f, c]]
    Ci = np.array([[a, d_, e_], [d_, b, f_], [e_, f_, c]]) * 1e15
    Linv = np.diag([1 / Lq1, 1 / Lc, 1 / Lq2]) * 1e9
    w2 = np.linalg.eigvals(Ci @ Linv)
    return np.sort(np.sqrt(np.abs(w2.real))) / 2 / np.pi / 1e9


def q3d_prior(path):
    """Cinv (1/fF) of Q1, TC, Q2 from a q3d_core summary (other nets grounded), in the fit parameter order."""
    C3 = json.load(open(path))['C_maxwell_fF']
    n = ['Q1', 'TC', 'Q2']
    Ci = np.linalg.inv(np.array([[C3[a][b] for b in n] for a in n]))
    return [Ci[0, 0], Ci[1, 1], Ci[2, 2], Ci[0, 1], Ci[0, 2], Ci[1, 2]]


def fit(out, cinv_guess=None, prior=None, sig_diag=0.10, sig_off=0.30, prior_mhz=5.0):
    """Least squares on the HFSS mode frequencies. prior (Q3D Cinv, q3d_prior): Gaussian prior with relative widths
    sig_diag / sig_off, weighted like a prior_mhz frequency error per parameter. With Q1 - Q2 detuned by ~380 MHz the
    direct Q1-Q2 term only enters at second order, so without the prior it is not determined (first fit: g12 32 MHz)."""
    from scipy.optimize import least_squares
    m = json.load(open(os.path.join(out, 'modes.json')))
    par = lambda a, b: a * b / (a + b)
    Lq1, Lq2 = par(m['Lj_nH']['Q1a'], m['Lj_nH']['Q1b']), par(m['Lj_nH']['Q2a'], m['Lj_nH']['Q2b'])
    data = [(r['L_c_nH'], np.array(r['f_GHz'])) for r in m['sweep']]
    if cinv_guess is None:
        cinv_guess = list(prior) if prior is not None else [1 / 97, 1 / 97, 1 / 97, 3.6 / 97 ** 2, 0.15 / 97 ** 2, 3.6 / 97 ** 2]
    pr = np.asarray(prior) if prior is not None else None
    sig = np.array([sig_diag] * 3 + [sig_off] * 3)

    def res(x):
        return np.concatenate([(model_freqs(x, Lq1, Lc, Lq2) - f) * 1e3 for Lc, f in data])   # MHz

    def res_all(x):
        r = res(x)
        if pr is None:
            return r
        return np.concatenate([r, prior_mhz * (np.asarray(x) - pr) / (sig * np.abs(pr))])

    sol = least_squares(res_all, cinv_guess, x_scale=np.abs(cinv_guess))
    a, b, c, d_, e_, f_ = sol.x
    Ci = np.array([[a, d_, e_], [d_, b, f_], [e_, f_, c]]) * 1e15   # 1/F
    names = ['Q1', 'TC', 'Q2']
    C = np.linalg.inv(Ci) * 1e15
    r = dict(names=names, Cinv_per_F=Ci.tolist(), C_fF=C.tolist(),
             Ec_GHz={n: float(e ** 2 * Ci[i, i] / 2 / h / 1e9) for i, n in enumerate(names)},
             g_n_MHz={'Q1-TC': float(4 * e ** 2 * Ci[0, 1] / h / 1e6), 'Q2-TC': float(4 * e ** 2 * Ci[2, 1] / h / 1e6),
                      'Q1-Q2': float(4 * e ** 2 * Ci[0, 2] / h / 1e6)},
             EJ_GHz={'Q1': ej_ghz(Lq1), 'Q2': ej_ghz(Lq2), 'TC_max (2 JJs at nominal Lj_TC)': None},
             Lj_drawn_nH=m['Lj_nH'], L_qubit_parallel_nH=dict(Q1=Lq1, Q2=Lq2), rms_MHz=float(np.sqrt(np.mean(res(sol.x) ** 2))),
             residual_MHz=res(sol.x).reshape(-1, 3).tolist(), sweep=m['sweep'],
             prior=None if pr is None else dict(Cinv_q3d_per_fF=pr.tolist(), sig_diag=sig_diag, sig_off=sig_off,
                                                weight_MHz=prior_mhz,
                                                pull_pct=(100 * (sol.x - pr) / np.abs(pr)).tolist()))
    r['bare_linear_GHz'] = {n: float(np.sqrt(Ci[i, i] / (L * 1e-9)) / 2 / np.pi / 1e9)
                            for i, (n, L) in enumerate(zip(names, (Lq1, m['Lj_nH']['TC'] / 2, Lq2)))}
    json.dump(r, open(os.path.join(out, 'fit.json'), 'w'), indent=1)
    print('fit rms %.2f MHz' % r['rms_MHz'])
    print('Ec (GHz):', {k: round(v, 4) for k, v in r['Ec_GHz'].items()})
    print('g_n (MHz):', {k: round(v, 2) for k, v in r['g_n_MHz'].items()})
    print('bare linear f (GHz) at nominal Lj:', {k: round(v, 4) for k, v in r['bare_linear_GHz'].items()})
    return r


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--lj', default='9,11.3,14,17,21')
    ap.add_argument('--tag', default='base')
    ap.add_argument('--date', default=time.strftime('%Y-%m-%d'))
    ap.add_argument('--mesh', default='10um')
    ap.add_argument('--maxp', type=int, default=6)
    ap.add_argument('--overrides', default='{}')
    ap.add_argument('--fit-only', action='store_true')
    ap.add_argument('--prior', default=os.path.join(SIM, 'results', '2026-11-14', 'q3d_v8', 'summary.json'),
                    help="q3d_core summary.json used as prior ('' = none)")
    a = ap.parse_args()
    out = os.path.join(SIM, 'results', a.date, f'hfss_core_{a.tag}')
    os.makedirs(out, exist_ok=True)
    if not a.fit_only:
        run(a, out)
    fit(out, prior=q3d_prior(a.prior) if a.prior else None)
