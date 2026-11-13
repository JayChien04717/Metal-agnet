"""Q3D sweep of the paper-hammer gap at a given filter-resonator spacing (purcell_6500 layout, whole channel).

Nets from the pin connectivity (shorted ends rendered OPEN, so no net merges with ground):
  input     : launch pad in, feed_in, C_in frame
  feed      : C_in plate, feed_mid, n0, feed_out, pad out, stub_p, C_p frame
  filter    : C_p plate, feed_jf, jf, filter_tail, stub_f, hammer head_a
  resonator : hammer head_b, stub_r, jr, res_open_run, readout_open, res_tail
Series C: C_in = -C[input, feed], C_p = -C[feed, filter], C_j = -C[filter, resonator].

Usage: python q3d/q3d_hammer_gap.py [--spacing 0.4] [--gaps 9,20,30,45,60,80]
Output: results/<date>/q3d<spacing>_gap_sweep/gap<g>um/{q3d_caps.json, q3d_matrix_fF.csv, q3d_convergence.csv},
        summary.json + cj_vs_gap.png
"""
import argparse, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np

SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DESIGN = os.path.dirname(SIM); HFSS = os.path.dirname(DESIGN)
sys.path[:0] = [os.path.join(DESIGN, 'pattern design'), os.path.join(HFSS, 'Design5', 'simulation')]
import purcell_6500 as P

LENGTHS = P.LENGTHS       # Cai 2024 design, 6.08 GHz
CAPS = ('cin', 'cp', 'cj')
SIDE = dict(frame='a', lead_a='a', head_a='a', plate='b', lead_b='b', head_b='b')


def nets(d):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    def union(a, b):
        parent[find(a)] = find(b)

    id2n = {c.id: n for n, c in d.components.items()}
    node = lambda comp, pin: (comp, pin) if comp in CAPS else comp
    for _, g in d.net_info.groupby('net_id'):
        ms = [node(id2n[int(r.component_id)], r.pin_name) for r in g.itertuples()]
        for m in ms[1:]:
            union(ms[0], m)
    for n in d.components:
        if n in CAPS:
            find((n, 'a')); find((n, 'b'))
        else:
            find(n)
    label = {find('feed_in'): 'input', find('n0'): 'feed', find('filter_tail'): 'filter', find('res_tail'): 'resonator'}
    objmap = {}
    for kind in ('path', 'poly'):
        for r in d.qgeometry.tables[kind].itertuples():
            if bool(r.subtract):
                continue
            comp = id2n[int(r.component)]
            k = (comp, SIDE.get(r.name, 'a')) if comp in CAPS else comp
            objmap[f'{r.name}_{comp}'] = label.get(find(k), '?')
    comps = [n for n, c in d.components.items() if type(c).__name__ != 'ResEnd']
    return objmap, comps


def run_gap(gap, spacing, out, args):
    from qiskit_metal.analyses.quantization import LOManalysis
    import save_aedt
    d = P.build(LENGTHS, paper_hammer=True, paper_hammer_gap_um=gap, spacing_mm=spacing)
    P.set_wire_bonds(d)                       # same auto bonds as HFSS (g_wb nets count as ground)
    objmap, comps = nets(d)
    lom = LOManalysis(d, 'q3d')
    lom.sim.renderer._options.update(P.WB_OPTIONS)
    for k, v in dict(max_passes=args.maxp, min_passes=2, percent_error=args.pct, freq_ghz=6.5,
                     solution_order='Medium', auto_increase_solution_order=False).items():
        lom.sim.setup[k] = v
    t0 = time.time()
    lom.sim.run(name=save_aedt.design_name(f'Q{spacing * 1e3:.0f}g{gap:g}'), components=comps,
                open_terminations=[('filter_tail', 'end'), ('res_tail', 'end')])
    cm = lom.sim.capacitance_matrix
    os.makedirs(out, exist_ok=True)
    cm.to_csv(os.path.join(out, 'q3d_matrix_fF.csv'))
    try:
        conv = lom.sim.renderer.pinfo.setup.get_convergence()
        (conv[0] if isinstance(conv, tuple) else conv).to_csv(os.path.join(out, 'q3d_convergence.csv'))
    except Exception as e:
        print('convergence not saved:', repr(e)[:120])
    isgnd = lambda n: 'ground' in n.lower() or n.startswith('g_wb')
    lab = {n: objmap.get(n, '?') for n in cm.index if not isgnd(n)}
    inv = {v: k for k, v in lab.items()}
    missing = {'input', 'feed', 'filter', 'resonator'} - set(inv)
    if missing:
        raise RuntimeError(f'nets not identified: {missing}; Q3D nets {lab}')
    gn = [n for n in cm.index if isgnd(n)][0]
    C = lambda a, b: float(-cm.loc[inv[a], inv[b]])
    order = ['input', 'feed', 'filter', 'resonator']
    res = dict(gap_um=gap, spacing_um=spacing * 1e3, minutes=round((time.time() - t0) / 60, 1),
               series_fF=dict(cin=C('input', 'feed'), cp=C('feed', 'filter'), cj=C('filter', 'resonator')),
               all_mutual_fF={f'{a}-{b}': C(a, b) for i, a in enumerate(order) for b in order[i + 1:]},
               net_ground_fF={k: float(-cm.loc[inv[k], gn]) for k in order}, net_names=lab)
    json.dump(res, open(os.path.join(out, 'q3d_caps.json'), 'w'), indent=1)
    print(f'gap {gap:g} um: C_j {res["series_fF"]["cj"]:.3f} fF, C_p {res["series_fF"]["cp"]:.2f} fF, '
          f'{res["minutes"]} min', flush=True)
    try:
        lom.sim.close()
    except Exception:
        pass
    return res


def summarize(root, gaps, spacing, targets):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    rows = []
    for g in gaps:
        p = os.path.join(root, f'gap{g:g}um', 'q3d_caps.json')
        if os.path.exists(p):
            rows.append(json.load(open(p)))
    if not rows:
        return
    g = np.array([r['gap_um'] for r in rows]); cj = np.array([r['series_fF']['cj'] for r in rows])
    picks = {}
    for name, c in targets.items():   # interpolate log C_j vs log gap
        if cj.min() <= c <= cj.max():
            o = np.argsort(cj)
            picks[name] = dict(cj_fF=c, gap_um=float(np.exp(np.interp(np.log(c), np.log(cj[o]), np.log(g[o])))))
    s = dict(spacing_um=spacing * 1e3, gap_um=g.tolist(), cj_fF=cj.tolist(),
             cp_fF=[r['series_fF']['cp'] for r in rows], targets=picks)
    json.dump(s, open(os.path.join(root, 'summary.json'), 'w'), indent=1)
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.loglog(g, cj, 'o-', label='Q3D C_j')
    for name, c in targets.items():
        ax.axhline(c, ls='--', lw=1, color='C3' if 'HFSS' in name else 'C2', label=f'{name}: {c:.2f} fF'
                   + (f' → gap {picks[name]["gap_um"]:.0f} µm' if name in picks else ''))
    ax.set_xlabel('hammer gap (µm)'); ax.set_ylabel('C_j (fF)'); ax.grid(alpha=.3, which='both')
    ax.set_title(f'paper hammer, {spacing * 1e3:.0f} µm spacing', fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(root, 'cj_vs_gap.png'), dpi=130)
    print('C_j vs gap:', ', '.join(f'{a:g} um {b:.3f} fF' for a, b in zip(g, cj)))
    for k, v in picks.items():
        print(f'target {k}: C_j {v["cj_fF"]:.2f} fF -> gap {v["gap_um"]:.1f} um')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--spacing', type=float, default=0.4, help='filter-resonator centreline spacing (mm)')
    ap.add_argument('--gaps', default='9,20,30,45,60,80')
    ap.add_argument('--maxp', type=int, default=8, help='Q3D passes: 8 (Jay 2026-11-12)')
    ap.add_argument('--pct', type=float, default=0.3)
    ap.add_argument('--target-j', type=float, default=4.9, help='J/2pi target (MHz)')
    ap.add_argument('--j-par', type=float, default=0.0, help='parasitic J without hammer (MHz), subtracted')
    ap.add_argument('--summary-only', action='store_true')
    ap.add_argument('--out')
    args = ap.parse_args()
    gaps = [float(x) for x in args.gaps.split(',')]
    root = args.out or os.path.join(SIM, 'results', time.strftime('%Y-%m-%d'), f'q3d{args.spacing * 1e3:.0f}_gap_sweep')
    jc = args.target_j - args.j_par
    targets = {'model 6.73 MHz/fF': jc / 6.73, 'HFSS 8.5 MHz/fF': jc / 8.5}
    if not args.summary_only:
        import hfss_dirs, save_aedt
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(SIM))
        for g in gaps:
            out = os.path.join(root, f'gap{g:g}um')
            if os.path.exists(os.path.join(out, 'q3d_caps.json')):
                print(f'gap {g:g} um: done before, skipped', flush=True); continue
            try:
                run_gap(g, args.spacing, out, args)
            except Exception as e:
                print(f'!!! gap {g:g} um failed: {e!r}', flush=True)
    summarize(root, gaps, args.spacing, targets)


if __name__ == '__main__':
    main()
