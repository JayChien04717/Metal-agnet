"""Eigenmode of each Design5 chip channel (chip6_layout.py, 6-7 GHz), one channel at a time.

Per channel k only its own components (*_ch<k>) and the two feed-line segments next to its T junction are rendered;
the far ends of those two segments are 50 ohm lumped ports (matched line, Wallraff Table 2.4 kappa_p definition).
Resonator open end listed in open_terminations; filter / resonator far ends shorted.
Solver rule (Jay 2026-11-12): 6 passes, 10 um length mesh on every line and capacitor; qiskit auto wire bonds (350 um).
2 modes from f_r - 0.15 GHz -> point_J inversion: f_p, f_r, kappa_p, J.

Usage: python eigenmode/eig_chip6.py [--ch 1,2,3,4,5,6] [--maxp 6] [--out results/<date>/eig_chip6]
Output: <out>/ch<k>/eigen_modes.json, <out>/summary.json
"""
import argparse, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np

SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DESIGN = os.path.dirname(SIM)
sys.path[:0] = [SIM, os.path.join(DESIGN, 'pattern design')]
MESH = '10um'
WB_OPTIONS = dict(wb_threshold='350um', wb_offset='0um', wb_size=2)


def point_J(freqs, imags):
    """2-mode inversion of [[f_p + i k_p/2, J], [J, f_r]] from the two complex eigenfrequencies."""
    from scipy.optimize import least_squares
    z = sorted(np.array(freqs) + 1j * np.array(imags), key=lambda c: c.real)
    s_, p_ = z[0] + z[1], z[0] * z[1]

    def r(x):
        fp, kp, fr, J = x; a = fp + 0.5j * kp
        e1, e2 = a + fr - s_, a * fr - J ** 2 - p_
        return 1e3 * np.array([e1.real, e1.imag, e2.real, e2.imag])
    wide = int(np.argmax([c.imag for c in z]))
    fp, kp, fr, J = least_squares(r, [z[wide].real, 2 * (z[0].imag + z[1].imag), z[1 - wide].real, 0.002]).x
    return dict(f_p_GHz=fp, kappa_p_MHz=1e3 * kp, f_r_GHz=fr, detuning_MHz=1e3 * (fr - fp), J_MHz=1e3 * abs(J))


def channel_subset(d, C6, ch):
    i = [o[0] for o in C6.ORDER].index(ch)
    comps = [n for n, c in d.components.items() if n.endswith(f'_ch{ch}') and type(c).__name__ != 'ResEnd']
    comps += [f'fl_{i}', f'fl_{i + 1}']
    ports = [(f'fl_{i}', 'start', 50), (f'fl_{i + 1}', 'end', 50)]
    return comps, ports


def run_ch(d, C6, ch, target, args, out):
    from qiskit_metal import Dict
    from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim
    import save_aedt
    comps, ports = channel_subset(d, C6, ch)
    t = d.qgeometry.tables['path']
    t['hfss_wire_bonds'] = False
    id2n = {c.id: n for n, c in d.components.items()}
    sel = t['component'].map(lambda c: id2n[int(c)] in comps) & t['subtract'].astype(bool)
    t.loc[sel, 'hfss_wire_bonds'] = True
    os.makedirs(out, exist_ok=True)
    json.dump(dict(ch=ch, target=target, components=comps, ports=ports, mesh=MESH, wire_bonds=WB_OPTIONS),
              open(os.path.join(out, 'layout.json'), 'w'), indent=1)

    class Eig(EigenmodeSim):
        def _analyze(self):
            objs = [f'{r.name}_{id2n[int(r.component)]}' for kind in ('path', 'poly')
                    for r in d.qgeometry.tables[kind].itertuples()
                    if id2n[int(r.component)] in comps and not bool(r.subtract)]
            self.renderer.modeler.mesh_length('mesh10', objs, MaxLength=MESH)
            print(f'ch{ch}: mesh {MESH} on {len(objs)} objects', flush=True)
            super()._analyze()

    eig = Eig(d, 'hfss')
    eig.renderer._options.update(WB_OPTIONS)
    minf = target['f_r'] - 0.15
    for k, v in dict(min_freq_ghz=minf, n_modes=2, max_delta_f=0.1, max_passes=args.maxp, min_passes=3,
                     min_converged=1, pct_refinement=20, basis_order=1, vars=Dict()).items():
        eig.setup[k] = v
    t0 = time.time()
    name = save_aedt.design_name(f'eig_ch{ch}')
    eig.run(name=name, components=comps, open_terminations=[(f'res_a_ch{ch}', 'start')], port_list=ports)
    setup = eig.renderer.pinfo.setup
    conv = setup.get_convergence()[0]; conv.to_csv(os.path.join(out, 'eigen_convergence.csv'))
    f, im = setup.get_solutions().eigenmodes()
    res = dict(frequencies_GHz=[float(x) for x in f], imag_GHz=[float(x) for x in im], passes=int(len(conv)),
               elements=int(conv.iloc[-1, 0]), last_delta_f_pct=float(conv.iloc[-1, -1]),
               minutes=round((time.time() - t0) / 60, 1), aedt_design=name + '_hfss')
    res.update(point_J(res['frequencies_GHz'], res['imag_GHz']))
    json.dump(res, open(os.path.join(out, 'eigen_modes.json'), 'w'), indent=1)
    print(f"ch{ch}: modes {f[0]:.4f}, {f[1]:.4f} GHz -> f_p {res['f_p_GHz']:.4f} f_r {res['f_r_GHz']:.4f} GHz, "
          f"kappa_p {res['kappa_p_MHz']:.1f} J {res['J_MHz']:.1f} MHz; {res['passes']} passes {res['minutes']} min", flush=True)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ch', default='1,2,3,4,5,6')
    ap.add_argument('--maxp', type=int, default=6, help='eigenmode passes: 6 (Jay 2026-11-12)')
    ap.add_argument('--out', default=os.path.join(SIM, 'results', time.strftime('%Y-%m-%d'), 'eig_chip6'))
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    import chip6_layout as C6
    d = C6.build_chip()
    L = {c['ch']: c for c in json.load(open(C6.LEN_JSON))['channels']}
    chans = [int(x) for x in args.ch.split(',')]
    for ch in chans:
        comps, ports = channel_subset(d, C6, ch)
        print(f'ch{ch} f_r {L[ch]["f_r"]}: {len(comps)} components, ports {ports}', flush=True)
    if args.dry_run:
        return
    import hfss_dirs, save_aedt
    hfss_dirs.use_d_drive(); save_aedt.use_project(save_aedt.project_file(SIM))
    summ = {}
    sp = os.path.join(args.out, 'summary.json')
    if os.path.exists(sp):
        summ = json.load(open(sp))
    for ch in chans:
        try:
            r = run_ch(d, C6, ch, L[ch], args, os.path.join(args.out, f'ch{ch}'))
            summ[f'ch{ch}'] = dict(target_f_r=L[ch]['f_r'], target_f_p=L[ch]['wp'], model_J=L[ch]['J'],
                                   model_kp=L[ch]['kp_matched'], **{k: r[k] for k in
                                   ('f_p_GHz', 'f_r_GHz', 'kappa_p_MHz', 'J_MHz', 'detuning_MHz', 'passes', 'minutes')})
            json.dump(summ, open(sp, 'w'), indent=1)
        except Exception as e:
            print(f'!!! ch{ch} failed: {e!r}', flush=True)
    try:
        save_aedt.snapshot(args.out, 'eig_chip6')
    except Exception as e:
        print('not snapshotted', repr(e)[:120])
    print('summary:', sp)


if __name__ == '__main__':
    main()
