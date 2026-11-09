"""Q3D of the WHOLE Design5 channel in one run (all capacitors, lines, launch pads together), instead of each
capacitor alone: the series capacitances then include the influence of the neighbouring metal, and the extra
mutual capacitances between the nets (e.g. filter <-> resonator along the parallel lines) become visible.

Nets (identified from the layout connectivity, the shorted ends are rendered OPEN so that no net merges with ground):
  P1_feed   : launch pad P1, feed_in, C_in electrode on pin a
  node      : C_in pin b, feed_c, n0, out_line, launch pad P2, stub_p, C_p electrode toward n0
  filter    : C_p electrode toward the filter, feed_jf, jf, stub_f, C_j head_a, filt
  resonator : C_j head_b, stub_r, jr, res_a, res
Series C:  C_in = -C[P1_feed, node],  C_p = -C[node, filter],  C_j = -C[filter, resonator]
Line capacitance check: C_line = L / (Z0 v) per net; excess = C_net,ground - C_line (electrodes + ends).

Usage: python q3d_full5.py [tag]     (env MAXP, PCT)
Output: results/q3d5_full_cmat_fF.csv (+ time-stamped copy), results/q3d5_full_caps.json (series, mutuals, nets)
"""
import os, sys, json, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES = os.path.join(HERE, 'results')   # simulation/ (this script: simulation/<kind>/)
import sys as _sys; _sys.path[:0] = [os.path.join(HERE, _k) for _k in ('driven', 'eigenmode', 'q3d')] + [HERE]
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
import design5_layout as d5, ftf_layout as f3

CAPS = ('cin', 'cp', 'cj')
SIDE = dict(frame='a', lead_a='a', head_a='a', plate='b', lead_b='b', head_b='b')
OPEN_ENDS = [('filt', 'start'), ('res', 'start'), ('res_a', 'start')]      # pins on the (pin-only) ResEnd components


def nets(d):
    """Union-find over the pin connections -> {object name as rendered: net label}, rendered component list."""
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b): parent[find(a)] = find(b)
    id2n = {c.id: n for n, c in d.components.items()}
    node = lambda comp, pin: (comp, pin) if comp in CAPS else comp
    ni = d.net_info
    for nid, g in ni.groupby('net_id'):
        ms = [node(id2n[int(r.component_id)], r.pin_name) for r in g.itertuples()]
        for m in ms[1:]: union(ms[0], m)
    for n in d.components: find(n) if n not in CAPS else (find((n, 'a')), find((n, 'b')))
    label = {}
    for key, name in (('feed_in', 'P1_feed'), ('out_line', 'node'), ('filt', 'filter'), ('res', 'resonator')):
        label[find(key)] = name
    comps = [n for n, c in d.components.items() if type(c).__name__ != 'ResEnd']
    objmap = {}
    for kind in ('path', 'poly'):
        t = d.qgeometry.tables[kind]
        for r in t.itertuples():
            if bool(r.subtract): continue
            comp = id2n[int(r.component)]
            k = (comp, SIDE[r.name]) if comp in CAPS else comp
            objmap[f'{r.name}_{comp}'] = label.get(find(k), '?')
    return objmap, comps


def line_lengths(d):
    pl = lambda n: f3.path_length(d, n)
    return dict(P1_feed=pl('feed_in'), node=pl('feed_c') + 0.3 + 0.1 + pl('out_line') + pl('stub_p'),
                filter=pl('feed_jf') + 0.3 + pl('stub_f') + pl('filt'), resonator=pl('stub_r') + 0.3 + pl('res_a') + pl('res'))


if __name__ == '__main__':
    tag = sys.argv[1] if len(sys.argv) > 1 else 'q3d5_full'
    d = d5.build_design()
    objmap, comps = nets(d)
    print('rendered:', comps); print('nets:', {k: sorted(o for o, v in objmap.items() if v == k) for k in set(objmap.values())})
    from qiskit_metal.analyses.quantization import LOManalysis
    lom = LOManalysis(d, 'q3d')
    for k, v in dict(max_passes=int(os.environ.get('MAXP', 15)), min_passes=2, percent_error=float(os.environ.get('PCT', 0.3)),
                     freq_ghz=6.7, solution_order='Medium', auto_increase_solution_order=False).items():
        lom.sim.setup[k] = v
    t0 = time.time()
    lom.sim.run(name=tag, components=comps, open_terminations=OPEN_ENDS)
    cm = lom.sim.capacitance_matrix
    import runstore
    run_dir = runstore.new_run_dir(tag, RES)                 # kept copy: results/<date>/<HHMM>_<tag>/
    cm.to_csv(os.path.join(RES, 'q3d5_full_cmat_fF.csv')); cm.to_csv(os.path.join(run_dir, 'q3d5_full_cmat_fF.csv'))
    # matrix nets -> labels
    lab = {}
    for n in cm.index:
        if 'ground' in n.lower(): continue
        lab[n] = objmap.get(n, '?')
    print('Q3D nets -> labels:', lab)
    inv = {v: k for k, v in lab.items()}
    C = lambda a, b: float(-cm.loc[inv[a], inv[b]])
    order = ['P1_feed', 'node', 'filter', 'resonator']
    mut = {f'{a}-{b}': C(a, b) for i, a in enumerate(order) for b in order[i + 1:] if a in inv and b in inv}
    Lline = line_lengths(d); Cp_um = 1 / (d5.Z0 * d5.V) * 1e15 * 1e-3       # fF per mm
    gn = [n for n in cm.index if 'ground' in n.lower()][0]
    gnd = {k: float(-cm.loc[inv[k], gn]) for k in order if k in inv}          # C to ground of each net
    out = dict(tag=tag, date=time.strftime('%Y-%m-%d %H:%M'), minutes=round((time.time() - t0) / 60, 1),
               layout=dict(variant=d5.VARIANT, CP_TYPE=d5.CP_TYPE, Cp_plate_um=d5.ENC_P['plate_width'], Cj_gap_um=d5.CJ_HG,
                           L_FILTER=d5.L_FILTER, L_RES_TOTAL=d5.L_RES_TOTAL, D_CJ=d5.D_CJ),
               series=dict(cin=mut.get('P1_feed-node'), cp=mut.get('node-filter'), cj=mut.get('filter-resonator')),
               mutual_fF=mut, ground_fF=gnd, line_mm=Lline,
               ground_excess_fF={k: gnd[k] - Cp_um * Lline[k] for k in gnd}, net_names=lab)
    json.dump(out, open(os.path.join(RES, 'q3d5_full_caps.json'), 'w'), indent=1)
    json.dump(out, open(os.path.join(run_dir, 'q3d5_full_caps.json'), 'w'), indent=1)
    print('SERIES (fF):', out['series']); print('ALL MUTUALS (fF):', mut); print('GROUND (fF):', gnd)
    print('GROUND EXCESS over line C (fF):', out['ground_excess_fF'])
    try: lom.sim.close()
    except Exception: pass
