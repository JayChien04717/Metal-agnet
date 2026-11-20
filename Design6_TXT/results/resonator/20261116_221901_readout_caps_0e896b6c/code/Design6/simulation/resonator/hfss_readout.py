"""HFSS eigenmode + pyEPR of one Design6 readout channel with its qubit.

Model: build_design() (full chip geometry), rendered: the qubit (Xmon, SQUID = two lumped JJs Lj_<Q>a / Lj_<Q>b), its Purcell channel
(n0 tee, C_p, filter, hammer C_j, resonator, open run into the claw) and the feed line from the launch pad to the
other channel's tee; the two feed ends are 50 Ohm lumped ports (so filter / resonator modes carry kappa).
Shorted ends (filter / resonator ResEnd) are shorted by the renderer. Box = rendered parts + buffer.
PEC ground straps on every straight CPW piece, <= 350 um apart (Design4 'straps cpw'; no qiskit auto bonds).
Eigenmode: 3 modes from 4 GHz (qubit, two hybrid filter / resonator modes), first order, 6 passes,
10 um length mesh on every CPW, capacitor and qubit electrode (skill.md section 12).

Results:
  * complex eigenfrequencies -> f_q, the two readout modes, their kappa; 2-mode inversion -> f_p, kappa_p, f_r, J
  * pyEPR (junction Lj_<Q>): anharmonicity alpha_q, cross-Kerr chi_qr / chi_qp (full dispersive shift 2chi),
    Lamb shifts, Purcell-limited T1 of the qubit (kappa of the qubit mode).

Usage:  python hfss_readout.py Q1 [--trim 0] [--tag v4]
Output: simulation/results/<date>/hfss_readout_<Q>_<tag>/{modes.json, epr.json, convergence.csv}
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np

SIM = _S('Design6/simulation')
DESIGN = _S('Design6'); REPO = _S('.')
sys.path[:0] = [_I('Design6/layout'), _I('Design5/simulation'),
                _I('Design4/simulation/resonator')]


def channel_components(d, q):
    s = '_' + q
    names = [n for n in d.components if n.endswith(s) and not n.startswith(('xy', 'cl_', 'fl_'))]
    other = 'Q2' if q == 'Q1' else 'Q1'
    feed = ['feed_in', 'feed_mid'] if q == 'Q1' else ['feed_out', 'feed_mid']
    return [q] + names + feed, other


class _Alias:
    """Design view in which the channel's tee / hammer components carry the plain Design4 names (n0, jf, jr, cj), so
    purcell_6500.ground_straps can be reused. jf must be the LEFT tap (hammer_ground_straps assumes x_jf < x_jr)."""
    def __init__(self, d, mapping):
        self._d = d
        self.components = {**{n: d.components[n] for n in d.components}, **{k: d.components[v] for k, v in mapping.items()}}

    def __getattr__(self, name):
        return getattr(self._d, name)


def strap_list(d, comps, q):
    """PEC ground straps: (1) the Design4 hammer-stub straps + every T-junction arm (n0, jf, jr) of this channel --
    without them the etched stub + hammer pocket cut the ground and HFSS shows a large non-capacitive filter-resonator
    coupling (first Design6 run: J 53.5 MHz; Design4: 68 -> 20 MHz); (2) every straight CPW piece, <= 350 um apart."""
    import purcell_6500 as P
    s = '_' + q
    left, right = sorted(('jf' + s, 'jr' + s), key=lambda n: d.parse_value(d.components[n].options.pos_x))
    alias = _Alias(d, {'n0': 'n0' + s, 'jf': left, 'jr': right, 'cj': 'cj' + s})
    ch = P.ground_straps(alias, 'all')
    ids = {d.components[n].id for n in comps}
    names = {c.id: n for n, c in d.components.items()}
    st = P.cpw_ground_straps(d, pitch_um=350)
    return P._drop_crossing(ch + [x for x in st if x['where'].split('.')[0] in {names[i] for i in ids}])


def run(a, out):
    import design6_layout as L
    import purcell_6500 as P
    from qiskit_metal import Dict
    from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim
    from pyEPR.ansys import parse_units
    import save_aedt
    q = a.qubit
    lj = {k: f'Lj_{k}' for k in ('Q1a', 'Q1b', 'Q2a', 'Q2b', 'TC')}
    lengths = getattr(a, 'lengths', None) or L.lengths_from_params()
    if a.trim:
        lengths[q] = dict(lengths[q], L_res=lengths[q]['L_res'] + a.trim / 1000)
    d = L.build_design(lj=lj, lengths=lengths, readout_hammer=getattr(a, 'readout_hammer', None),
                       coupler=dict(flux_as_ground='True', squid_contact_overlap='0um',
                                                             squid_shift_toward_stem='0um'))
    P.set_wire_bonds(d, False)        # straps only (auto bonds + straps intersect in HFSS: 'Parts g_wb_* intersect')
    comps, other = channel_components(d, q)
    ports = ([('feed_in', 'start', 50), ('feed_mid', 'end', 50)] if q == 'Q1'
             else [('feed_out', 'start', 50), ('feed_mid', 'start', 50)])
    straps = strap_list(d, [c for c in comps if c != q], q)
    save_aedt.use_project(_S('Design6/simulation/aedt/Design6.aedt'))
    id2n = {c.id: n for n, c in d.components.items()}
    metal = [f'{r.name}_{id2n[int(r.component)]}' for kind in ('poly', 'path')
             for r in d.qgeometry.tables[kind].itertuples()
             if not bool(r.subtract) and id2n[int(r.component)] in comps]

    class Eig(EigenmodeSim):
        def _analyze(self):
            for s in straps:
                self.renderer.modeler.draw_wirebond(
                    pos=parse_units(np.asarray(s['pos'])), ori=np.asarray(s['ori']), width=parse_units(s['span']),
                    height=parse_units(s['span']), z=0, wire_diameter=f"{s['diameter']}mm", NumSides=6,
                    name='g_wb', material='pec', solve_inside=False)
            self.renderer.modeler.mesh_length('metal_mesh', metal, MaxLength=a.mesh)
            print(f'{len(straps)} straps, mesh {a.mesh} on {len(metal)} objects', flush=True)
            super()._analyze()

    eig = Eig(d, 'hfss')
    lj0 = {k: v for k, v in L.LJ.items()}
    vars0 = Dict(**{f'Lj_{k}': v for k, v in lj0.items()}, Cj='0fF')
    for k, v in dict(min_freq_ghz=4.0, n_modes=3, max_delta_f=0.1, max_passes=a.maxp, min_passes=2, min_converged=1,
                     pct_refinement=30, basis_order=1, vars=vars0).items():
        eig.setup[k] = v
    t0 = time.time()
    dname = getattr(a, 'design_name', None) or save_aedt.design_name(f'ro{q}_{a.tag}')
    eig.run(name=dname, components=comps, open_terminations=[], port_list=ports)
    actual_design = eig.renderer.pinfo.design.name
    if actual_design != dname + '_hfss':
        raise RuntimeError(f'Unexpected HFSS design: {actual_design}')
    setup = eig.renderer.pinfo.setup
    conv = setup.get_convergence()[0]
    conv.to_csv(_relocate(os.path.join(out, 'convergence.csv')))
    if getattr(a, 'quick', False) and eig.convergence_f is not None:
        eig.convergence_f.to_csv(_relocate(os.path.join(out, 'convergence_f.csv')))
    f, im = setup.get_solutions().eigenmodes()
    modes = dict(design=dname + '_hfss', qubit=q, trim_um=a.trim, lengths=lengths[q], straps=len(straps),
                 f_GHz=[float(x) for x in f], kappa_MHz=[2e3 * float(x) for x in im],
                 passes=int(len(conv)), last_delta_f_pct=float(conv.iloc[-1, -1]), minutes=round((time.time() - t0) / 60, 1))
    json.dump(modes, open(_relocate(os.path.join(out, 'modes.json')), 'w'), indent=1)
    print('modes:', ', '.join(f'{x:.5f} GHz (kappa {k:.3f} MHz)' for x, k in zip(modes['f_GHz'], modes['kappa_MHz'])),
          f"| {modes['passes']} passes, {modes['minutes']} min", flush=True)
    if getattr(a, 'quick', False):
        return modes

    # ---- pyEPR ----
    try:
        from pyEPR import DistributedAnalysis, QuantumAnalysis
        pinfo = eig.renderer.pinfo
        cid = d.components[q].id
        for k, ab in ((1, 'a'), (2, 'b')):        # asymmetric SQUID: two junctions (xmon.py rect_jj1 / rect_jj2)
            pinfo.junctions[f'j{q}{ab}'] = dict(Lj_variable=f'Lj_{q}{ab}', rect=f'JJ_rect_Lj_{q}_rect_jj{k}',
                                                line=f'JJ_Lj_{q}_rect_jj{k}_', Cj_variable='Cj')
        pinfo.validate_junction_info()
        pinfo.dissipative['dielectrics_bulk'] = ['main']
        epra = DistributedAnalysis(pinfo)
        epra.do_EPR_analysis()
        qa = QuantumAnalysis(epra.data_filename)
        r = qa.analyze_all_variations(cos_trunc=8, fock_trunc=7, print_result=False)
        var = list(r.keys())[0]
        rr = r[var]
        chi = np.array(rr['chi_ND'])                     # MHz, numerical diag. cross-Kerr matrix
        fND = np.array(rr['f_ND']) / 1e3 if np.max(rr['f_ND']) > 100 else np.array(rr['f_ND'])
        epr = dict(variation=var, f_ND_GHz=fND.tolist(), f_0_GHz=(np.array(rr['f_0']) / 1e3).tolist(),
                   chi_ND_MHz=chi.tolist(), epr_data=epra.data_filename)
        json.dump(epr, open(_relocate(os.path.join(out, 'epr.json')), 'w'), indent=1, default=float)
        print('EPR f_ND (GHz):', np.round(fND, 5), '\nchi_ND (MHz):\n', np.round(chi, 4), flush=True)
    except Exception as ex:
        import traceback; traceback.print_exc()
        print('EPR failed:', ascii(repr(ex))[:300], flush=True)
    try:
        save_aedt.snapshot(out, f'ro{q}_{a.tag}')
    except Exception as ex:
        print('snapshot failed:', repr(ex)[:120])
    return modes


def summarize(out):
    """Mode assignment and 2-mode inversion of the readout pair (Design4 point_J)."""
    from hfss_hammer_straps import point_J
    m = json.load(open(_relocate(os.path.join(out, 'modes.json'))))
    f, k = np.array(m['f_GHz']), np.array(m['kappa_MHz'])
    iq = int(np.argmin(f))                       # qubit: lowest mode (4.6 - 4.8 GHz)
    pair = [i for i in range(len(f)) if i != iq]
    res = dict(f_q_GHz=float(f[iq]), kappa_q_MHz=float(k[iq]), T1_purcell_ms=float(1 / (2 * np.pi * k[iq] * 1e6) * 1e3),
               readout_modes_GHz=f[pair].tolist(), readout_kappa_MHz=k[pair].tolist())
    if len(pair) == 2:
        res.update(point_J(f[pair], k[pair] / 2e3))
    p = _relocate(os.path.join(out, 'epr.json'))
    if os.path.exists(p):
        e = json.load(open(p))
        chi = np.array(e['chi_ND_MHz'])
        res['alpha_q_MHz'] = float(chi[iq, iq])
        res['chi_q_modes_MHz'] = [float(chi[iq, j]) for j in pair]   # pyEPR chi = full shift 2chi (Kerr convention)
        res['f_ND_GHz'] = e['f_ND_GHz']
    json.dump(res, open(_relocate(os.path.join(out, 'summary.json')), 'w'), indent=1)
    print(json.dumps(res, indent=1))
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('qubit', choices=('Q1', 'Q2'))
    ap.add_argument('--trim', type=float, default=0.0, help='resonator length change, um')
    ap.add_argument('--tag', default='v4')
    ap.add_argument('--date', default=time.strftime('%Y-%m-%d'))
    ap.add_argument('--mesh', default='20um')   # coarse on purpose (Jay 2026-11-14: frequencies need not be precise)
    ap.add_argument('--maxp', type=int, default=4)
    ap.add_argument('--summary-only', action='store_true')
    a = ap.parse_args()
    out = _relocate(os.path.join(SIM, 'results', a.date, f'hfss_readout_{a.qubit}_{a.tag}'))
    os.makedirs(out, exist_ok=True)
    if not a.summary_only:
        run(a, out)
    summarize(out)
