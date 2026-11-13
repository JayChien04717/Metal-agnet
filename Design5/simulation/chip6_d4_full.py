"""Whole Design5 chip (chip6_d4.py: 6 Design4/Cai channels, 6.0-7.0 GHz) simulated in ONE model, all 6 channels together.

  q3d : capacitance matrix of the full chip, 8 passes (Jay 2026-11-12). Nets: input, feed, filter1..6, res1..6
        (shorted ends rendered open). Series C per channel: C_p = -C[feed, filter_k], C_j = -C[filter_k, res_k];
        C_in = -C[input, feed]; cross-channel mutuals kept.
  eig : eigenmode of the full chip, 12 modes from 5.8 GHz, 6 passes, 10 um length mesh on every line and capacitor,
        50 ohm lumped ports on the two feed-line ends (feed_in start, last fl end).

Ground: qiskit auto wire bonds every 350 um (P.WB_OPTIONS) + per channel the hammer-stub and T-junction straps
(purcell_6500.ground_straps 'all'), same as Design4.

Usage: python chip6_d4_full.py q3d|eig [--maxp N] [--out DIR]
Output: results/<date>/chip6_d4_full/{q3d,eig}/...
"""
import argparse, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
SIM = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(SIM)
sys.path[:0] = [SIM, os.path.join(ROOT, 'pattern design'), os.path.join(os.path.dirname(ROOT), 'Design4', 'pattern design')]
import chip6_d4 as C6
P = C6.P
CHS = [o[0] for o in C6.ORDER]
SIDE = dict(frame='a', lead_a='a', head_a='a', plate='b', lead_b='b', head_b='b')
MESH, MESH_LINE = '10um', '30um'


class ChannelView:
    """Lets purcell_6500's strap functions (hard-coded names cj, jf, jr, n0) work on channel k of the chip."""
    def __init__(self, d, ch):
        self.components = {n: d.components[f'{n}_ch{ch}'] for n in ('cj', 'jf', 'jr', 'n0')}
        self.qgeometry, self.parse_value = d.qgeometry, d.parse_value


def chip_straps(d):
    out = []
    for ch in CHS:
        out += [dict(s, where=f"ch{ch}.{s.get('where', '')}") for s in P.ground_straps(ChannelView(d, ch), 'all')]
    return out


def draw_straps(renderer, straps):
    from pyEPR.ansys import parse_units
    return [renderer.modeler.draw_wirebond(pos=parse_units(np.asarray(s['pos'])), ori=np.asarray(s['ori']),
                                           width=parse_units(s['span']), height=parse_units(s['span']), z=0,
                                           wire_diameter=f"{s['diameter']}mm", NumSides=6, name='g_wb',
                                           material='pec', solve_inside=False) for s in straps]


def caps_of(comp):
    return comp == 'cin' or comp.startswith(('cp_ch', 'cj_ch'))


def nets(d):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    id2n = {c.id: n for n, c in d.components.items()}
    node = lambda comp, pin: (comp, pin) if caps_of(comp) else comp
    for _, g in d.net_info.groupby('net_id'):
        ms = [node(id2n[int(r.component_id)], r.pin_name) for r in g.itertuples()]
        for m in ms[1:]:
            parent[find(m)] = find(ms[0])
    label = {find('feed_in'): 'input', find('fl_0'): 'feed'}
    for ch in CHS:
        label[find(f'filter_tail_ch{ch}')] = f'filter{ch}'; label[find(f'res_tail_ch{ch}')] = f'res{ch}'
    objmap = {}
    for kind in ('path', 'poly'):
        for r in d.qgeometry.tables[kind].itertuples():
            if bool(r.subtract):
                continue
            comp = id2n[int(r.component)]
            k = (comp, SIDE.get(r.name, 'a')) if caps_of(comp) else comp
            objmap[f'{r.name}_{comp}'] = label.get(find(k), '?')
    return objmap


def comps_of(d):
    return [n for n, c in d.components.items() if type(c).__name__ != 'ResEnd']


def mesh_all(renderer, d, comps):
    id2n = {c.id: n for n, c in d.components.items()}
    # 24 GB RAM: 10 um everywhere runs out of memory at pass 2 (2026-11-13) -> 10 um on capacitors / hammer region,
    # MESH_LINE on the long CPW lines (Jay 2026-11-13: "coarser line mesh")
    fine = lambda n: n == 'cin' or n.startswith(('cp_ch', 'cj_ch', 'jf_ch', 'jr_ch', 'stub_f_ch', 'stub_r_ch', 'n0_ch', 'stub_p_ch'))
    objs = [(f'{r.name}_{id2n[int(r.component)]}', fine(id2n[int(r.component)])) for kind in ('path', 'poly')
            for r in d.qgeometry.tables[kind].itertuples()
            if id2n[int(r.component)] in comps and not bool(r.subtract) and not id2n[int(r.component)].startswith('LP')]
    caps = [o for o, f in objs if f]; lines = [o for o, f in objs if not f]
    renderer.modeler.mesh_length('mesh_cap', caps, MaxLength=MESH)
    renderer.modeler.mesh_length('mesh_line', lines, MaxLength=MESH_LINE)
    print(f'mesh {MESH} on {len(caps)} capacitor objects, {MESH_LINE} on {len(lines)} line objects', flush=True)


def run_q3d(args, out):
    from qiskit_metal.analyses.quantization import LOManalysis
    import save_aedt
    # no bonds / straps in Q3D: on the full chip the merged 'g_wb' net was rejected as not fully connected (2026-11-12)
    d = C6.build_chip(); P.set_wire_bonds(d, False)
    objmap, comps, straps = nets(d), comps_of(d), []
    lom = LOManalysis(d, 'q3d')
    lom.sim.renderer._options.update(P.WB_OPTIONS)
    orig = lom.sim._analyze

    def _analyze():
        print('straps drawn:', len(draw_straps(lom.sim.renderer, straps)), flush=True)
        orig()
    lom.sim._analyze = _analyze
    for k, v in dict(max_passes=args.maxp or 8, min_passes=2, percent_error=0.3, freq_ghz=6.5,
                     solution_order='Medium', auto_increase_solution_order=False).items():
        lom.sim.setup[k] = v
    t0 = time.time()
    opens = [(f'{t}_tail_ch{ch}', 'end') for ch in CHS for t in ('filter', 'res')]
    lom.sim.run(name=save_aedt.design_name('chip6_d4_q3d'), components=comps, open_terminations=opens)
    cm = lom.sim.capacitance_matrix
    os.makedirs(out, exist_ok=True); cm.to_csv(os.path.join(out, 'q3d_matrix_fF.csv'))
    isgnd = lambda n: 'ground' in n.lower() or n.startswith('g_wb')
    lab = {n: objmap.get(n, '?') for n in cm.index if not isgnd(n)}
    inv = {v: k for k, v in lab.items()}
    C = lambda a, b: float(-cm.loc[inv[a], inv[b]])
    per = {f'ch{ch}': dict(cp=C('feed', f'filter{ch}'), cj=C(f'filter{ch}', f'res{ch}')) for ch in CHS}
    names = [v for v in inv if v != '?']
    cross = {f'{a}-{b}': C(a, b) for i, a in enumerate(names) for b in names[i + 1:]
             if a[-1].isdigit() and b[-1].isdigit() and a[-1] != b[-1] and C(a, b) > 0.01}
    res = dict(minutes=round((time.time() - t0) / 60, 1), cin=C('input', 'feed'), channels=per,
               cross_channel_fF=cross, net_names=lab, straps=len(straps))
    json.dump(res, open(os.path.join(out, 'q3d_caps.json'), 'w'), indent=1)
    print(f"C_in {res['cin']:.2f} fF; " + '; '.join(f"{k}: C_p {v['cp']:.2f} C_j {v['cj']:.3f}" for k, v in per.items()))
    print('largest cross-channel C (fF):', sorted(cross.items(), key=lambda kv: -kv[1])[:5], flush=True)


def run_eig(args, out):
    from qiskit_metal import Dict
    from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim
    import save_aedt
    d = C6.build_chip(); P.set_wire_bonds(d)
    comps, straps = comps_of(d), chip_straps(d)

    class Eig(EigenmodeSim):
        def _analyze(self):
            print('straps drawn:', len(draw_straps(self.renderer, straps)), flush=True)
            mesh_all(self.renderer, d, comps)
            super()._analyze()
    eig = Eig(d, 'hfss')
    eig.renderer._options.update(P.WB_OPTIONS)
    for k, v in dict(min_freq_ghz=5.8, n_modes=12, max_delta_f=0.1, max_passes=args.maxp or 6, min_passes=3,
                     min_converged=1, pct_refinement=20, basis_order=1, vars=Dict()).items():
        eig.setup[k] = v
    t0 = time.time()
    nl = len(C6.ORDER)
    eig.run(name=save_aedt.design_name('chip6_d4_eig'), components=comps, open_terminations=[],
            port_list=[('feed_in', 'start', 50), (f'fl_{nl}', 'end', 50)])
    setup = eig.renderer.pinfo.setup
    os.makedirs(out, exist_ok=True)
    conv = setup.get_convergence()[0]; conv.to_csv(os.path.join(out, 'eigen_convergence.csv'))
    f, im = setup.get_solutions().eigenmodes()
    res = dict(frequencies_GHz=[float(x) for x in f], kappa_MHz=[2e3 * float(x) for x in im], passes=int(len(conv)),
               elements=int(conv.iloc[-1, 0]), last_delta_f_pct=float(conv.iloc[-1, -1]),
               minutes=round((time.time() - t0) / 60, 1), straps=len(straps))
    json.dump(res, open(os.path.join(out, 'eigen_modes.json'), 'w'), indent=1)


def point_J(freqs, imags):
    from scipy.optimize import least_squares
    z = sorted(np.array(freqs) + 1j * np.array(imags), key=lambda c: c.real)
    s_, p_ = z[0] + z[1], z[0] * z[1]

    def r(x):
        fp, kp, fr, J = x; a = fp + 0.5j * kp
        e1, e2 = a + fr - s_, a * fr - J ** 2 - p_
        return 1e3 * np.array([e1.real, e1.imag, e2.real, e2.imag])
    wide = int(np.argmax([c.imag for c in z]))
    fp, kp, fr, J = least_squares(r, [z[wide].real, 2 * (z[0].imag + z[1].imag), z[1 - wide].real, 0.002]).x
    return dict(f_p_GHz=fp, kappa_p_MHz=1e3 * kp, f_r_GHz=fr, J_MHz=1e3 * abs(J))


def pair_modes(out):
    m = json.load(open(os.path.join(out, 'eigen_modes.json')))
    z = sorted(zip(m['frequencies_GHz'], [k / 2e3 for k in m['kappa_MHz']]))
    rows = {}
    for i, f in enumerate(C6.F_CHIP):
        pair = sorted(z, key=lambda t: abs(t[0] - f))[:2]
        rows[f'ch{i + 1}'] = dict(target_GHz=f, modes_GHz=[p[0] for p in sorted(pair)], **point_J(*zip(*pair)))
    json.dump(rows, open(os.path.join(out, 'channels.json'), 'w'), indent=1)
    for k, r in rows.items():
        print(f"{k} target {r['target_GHz']:.2f}: modes {r['modes_GHz'][0]:.4f} {r['modes_GHz'][1]:.4f} -> f_p {r['f_p_GHz']:.4f} "
              f"f_r {r['f_r_GHz']:.4f} GHz, kappa_p {r['kappa_p_MHz']:.1f} J {r['J_MHz']:.1f} MHz", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('what', choices=('q3d', 'eig', 'pair', 'dry'))
    ap.add_argument('--maxp', type=int)
    ap.add_argument('--out', default=os.path.join(SIM, 'results', time.strftime('%Y-%m-%d'), 'chip6_d4_full'))
    args = ap.parse_args()
    if args.what == 'dry':
        d = C6.build_chip(); st = chip_straps(d)
        print(len(st), 'straps;', len(comps_of(d)), 'components;', 'nets:', sorted(set(nets(d).values()))); return
    if args.what == 'pair':
        pair_modes(os.path.join(args.out, 'eig')); return
    import hfss_dirs, save_aedt
    hfss_dirs.use_d_drive(); save_aedt.use_project(save_aedt.project_file(SIM))
    if args.what == 'q3d':
        run_q3d(args, os.path.join(args.out, 'q3d'))
    else:
        run_eig(args, os.path.join(args.out, 'eig'))
        pair_modes(os.path.join(args.out, 'eig'))


if __name__ == '__main__':
    main()
