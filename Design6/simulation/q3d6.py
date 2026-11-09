"""Q3D of a Design6 verification layout, whole channel in one run (all metal together).

Nets (from the pin connectivity; the shorted ends are rendered OPEN so no net merges with ground):
  node      : launch pads, feed_in, n0, feed_out, stub_p, C_p electrode toward n0
  filter    : C_p electrode toward the filter, feed_jf, jf, stub_f, C_j head_a, filt
  resonator : C_j head_b (+ stub_r, jr, res_a, res in VARIANT full; only the floating plate in VARIANT filter)
Series C: C_p = -C[node, filter], C_j = -C[filter, resonator]; ground C of every net = -C[net, ground].

Usage: python q3d6.py [full|filter]        (env MAXP, PCT)
Output: results/<date>/<HHMM>_q3d6_<variant>/q3d6_caps.json + cmat; latest copy results/q3d6_<variant>_caps.json
"""
import os, sys, json, time, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results')
HFSS = os.path.dirname(os.path.dirname(HERE))
os.environ.setdefault('WIREBONDS', '0')                      # no HFSS wire bonds in Q3D
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, os.path.join(HFSS, 'Design5', 'simulation'))
import verify_layout as V
CH = V.CH
CAPS = (f'cp_ch{CH}', f'cj_ch{CH}')
SIDE = dict(frame='a', lead_a='a', head_a='a', plate='b', lead_b='b', head_b='b')


def nets(d):
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b): parent[find(a)] = find(b)
    id2n = {c.id: n for n, c in d.components.items()}
    node = lambda comp, pin: (comp, pin) if comp in CAPS else comp
    for _, g in d.net_info.groupby('net_id'):
        ms = [node(id2n[int(r.component_id)], r.pin_name) for r in g.itertuples()]
        for m in ms[1:]: union(ms[0], m)
    for n in d.components:
        if n in CAPS: find((n, 'a')); find((n, 'b'))
        else: find(n)
    label = {find('feed_in'): 'node', find(f'filt_ch{CH}'): 'filter', find((f'cj_ch{CH}', 'b')): 'resonator'}
    objmap = {}
    for kind in ('path', 'poly'):
        for r in d.qgeometry.tables[kind].itertuples():
            if bool(r.subtract): continue
            comp = id2n[int(r.component)]
            k = (comp, SIDE.get(r.name, 'a')) if comp in CAPS else comp
            objmap[f'{r.name}_{comp}'] = label.get(find(k), '?')
    comps = [n for n, c in d.components.items() if type(c).__name__ != 'ResEnd']
    return objmap, comps


if __name__ == '__main__':
    variant = sys.argv[1] if len(sys.argv) > 1 else 'full'
    tag = f'q3d6_{variant}'
    d = V.build(variant)
    objmap, comps = nets(d)
    opens = [(f'filt_ch{CH}', 'start')] + ([(f'res_ch{CH}', 'start'), (f'res_a_ch{CH}', 'start')] if variant == 'full' else [])
    print('nets:', {k: sorted(o for o, v in objmap.items() if v == k) for k in set(objmap.values())})
    from qiskit_metal.analyses.quantization import LOManalysis
    lom = LOManalysis(d, 'q3d')
    for k, v in dict(max_passes=int(os.environ.get('MAXP', 15)), min_passes=2, percent_error=float(os.environ.get('PCT', 0.3)),
                     freq_ghz=6.7, solution_order='Medium', auto_increase_solution_order=False).items():
        lom.sim.setup[k] = v
    t0 = time.time()
    lom.sim.run(name=tag, components=comps, open_terminations=opens)
    cm = lom.sim.capacitance_matrix
    import runstore
    run_dir = runstore.new_run_dir(tag, RES)
    cm.to_csv(os.path.join(run_dir, 'q3d6_cmat_fF.csv'))
    lab = {n: objmap.get(n, '?') for n in cm.index if 'ground' not in n.lower()}
    print('Q3D nets -> labels:', lab)
    inv = {v: k for k, v in lab.items()}
    gn = [n for n in cm.index if 'ground' in n.lower()][0]
    C = lambda a, b: float(-cm.loc[inv[a], inv[b]])
    order = [k for k in ('node', 'filter', 'resonator') if k in inv]
    out = dict(tag=tag, variant=variant, date=time.strftime('%Y-%m-%d %H:%M'), minutes=round((time.time() - t0) / 60, 1),
               layout=V.info(), series=dict(cp=C('node', 'filter'), cj=C('filter', 'resonator')),
               mutual_fF={f'{a}-{b}': C(a, b) for i, a in enumerate(order) for b in order[i + 1:]},
               ground_fF={k: float(-cm.loc[inv[k], gn]) for k in order}, net_names=lab)
    json.dump(out, open(os.path.join(run_dir, 'q3d6_caps.json'), 'w'), indent=1)
    json.dump(dict(out, run_dir=os.path.relpath(run_dir, HERE)), open(os.path.join(RES, f'q3d6_{variant}_caps.json'), 'w'), indent=1)
    print('SERIES (fF):', out['series']); print('MUTUAL:', out['mutual_fF']); print('GROUND:', out['ground_fF'])
    print('run folder', run_dir)
    try:                                               # keep the AEDT project with the run
        import save_aedt; save_aedt.save_active(run_dir, tag, close=True)
    except Exception as e:
        print('AEDT project not saved:', repr(e)[:150])
    try: lom.sim.close()
    except Exception: pass
