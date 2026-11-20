"""Full Q3D of the Design6_TXT r9 readout-only geometry (the same components as readout_s21.py), one run.

Purpose: a static reference for the readout circuit model. HFSS run 20261117_071647_readout_s21_clawdl_cf91f657
(+0.1 mm L_res, no claw) is the S21 reference; this Q3D solves the identical geometry (same build(), same components),
so every net-to-net and net-to-ground capacitance of the model can be checked against one solution that includes all
neighbouring metal (Cp frame, tees, hammer, parallel lines, cross-channel).

Nets (shorted ends and both feed ends rendered OPEN so that no net merges with ground):
  feed      : feed_in, feed_mid, feed_out, n0_Q*, stub_p_Q*, cp_Q* side a (frame / lead_a)
  filter_Q* : cp_Q* side b, feed_jf, jf, stub_f, cj side a, filter_tail
  res_Q*    : cj side b, stub_r, jr, res_open, res_tail
Ground: the ground plane only (no straps; Q3D treats every strap net as ground, and the HFSS straps only join ground).
Solver: 20 passes / 0.1 % (the Design6_TXT readout_caps setting).
"""
from argparse import ArgumentParser
from pathlib import Path
import json
import sys
import time
import traceback

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'Design6_TXT/simulation/resonator')]
import readout_s21 as R                       # noqa: E402  (sets sys.path)
from run_record import create_run, finish_run  # noqa: E402

CAPS = ('cp_', 'cj_')


def side(obj):
    return 'a' if obj.startswith(('frame', 'lead_a', 'head_a')) else 'b'


def net_map(d, comps):
    """Union-find over pin connections, capacitors split by electrode side -> {rendered object: net label}."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    id2n = {c.id: n for n, c in d.components.items()}
    key = lambda comp, pin: (comp, pin) if comp.startswith(CAPS) else comp
    for _, g in d.net_info.groupby('net_id'):
        ms = [key(id2n[int(r.component_id)], r.pin_name) for r in g.itertuples() if id2n[int(r.component_id)] in comps]
        for m in ms[1:]:
            parent[find(ms[0])] = find(m)
    label = {find('feed_in'): 'feed'}
    for q in ('Q1', 'Q2'):
        label[find(f'filter_tail_{q}')] = f'filter_{q}'
        label[find(f'res_tail_{q}')] = f'res_{q}'
    objmap = {}
    for kind in ('path', 'poly'):
        for r in d.qgeometry.tables[kind].itertuples():
            comp = id2n[int(r.component)]
            if bool(r.subtract) or comp not in comps:
                continue
            k = (comp, side(r.name)) if comp.startswith(CAPS) else comp
            objmap[f'{r.name}_{comp}'] = label.get(find(k), '?')
    assert '?' not in objmap.values(), {o: v for o, v in objmap.items() if v == '?'}
    return objmap


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--res-extra-mm', type=float, default=0.1)
    ap.add_argument('--max-passes', type=int, default=20)
    ap.add_argument('--percent-error', type=float, default=0.1)
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    d, hammer, lengths = R.build(args.res_extra_mm)
    comps = R.components(d)
    opens = [('feed_in', 'start'), ('feed_out', 'start')] + [(f'{n}_{q}', 'end') for q in ('Q1', 'Q2')
                                                             for n in ('res_open', 'filter_tail', 'res_tail')]
    objmap = net_map(d, comps)
    nets = {v: sorted(o for o, w in objmap.items() if w == v) for v in sorted(set(objmap.values()))}
    folder = create_run('Design6_TXT', 'resonator', 'readout_q3d_full', entrypoint=__file__,
                        inputs=[ROOT / 'Design6_TXT/design.json', Path(R.L._length_file())],
                        upstream=[str(ROOT / 'Design6_TXT/results/resonator/20261117_071647_readout_s21_clawdl_cf91f657')],
                        parameters=dict(geometry='r9 readout-only, identical to readout_s21.build(res_extra_mm)',
                                        res_extra_mm=args.res_extra_mm, lengths_mm=lengths, readout_hammer_um=hammer,
                                        components=comps, open_terminations=opens, nets=nets, ground_straps='none',
                                        Q3D=dict(max_passes=args.max_passes, min_passes=2, percent_error=args.percent_error,
                                                 solution_order='Medium'), dry_run=args.dry_run))
    print('RUN_FOLDER=' + str(folder), flush=True)
    print(json.dumps(nets, indent=1), flush=True)
    try:
        import hammer_deembed as H
        H.write(folder / 'raw/geometry.json', H.geometry_dump(d))
        if args.dry_run:
            finish_run(folder, status='completed', summary=dict(dry_run=True, nets=nets)); return
        import hfss_dirs, save_aedt
        from qiskit_metal.analyses.quantization import LOManalysis
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(ROOT / 'Design6_TXT/simulation'))
        lom = LOManalysis(d, 'q3d')
        lom.sim.setup.update(max_passes=args.max_passes, min_passes=2, percent_error=args.percent_error, freq_ghz=7.0,
                             solution_order='Medium', auto_increase_solution_order=False)
        t0 = time.time()
        lom.sim.run(name='D6Q3DFULL_' + folder.name[-8:], components=comps, open_terminations=opens)
        cm = lom.sim.capacitance_matrix
        if str(lom.sim.units).lower() not in ('ff', 'femtofarad', 'femtofarads'):
            raise RuntimeError(f'Unexpected units {lom.sim.units}')
        cm.to_csv(folder / 'raw/matrix_fF.csv')
        sig = [n for n in cm.index if 'ground' not in n.lower()]
        lab = {n: objmap.get(n, '?') for n in sig}
        assert sorted(lab.values()) == sorted(nets), f'net labels {lab}'
        inv = {v: k for k, v in lab.items()}
        order = ['feed', 'filter_Q1', 'res_Q1', 'filter_Q2', 'res_Q2']
        M = cm.loc[[inv[k] for k in order], [inv[k] for k in order]].to_numpy(float)
        passes = []
        for k, m in lom.sim.capacitance_all_passes.items():
            f = m if isinstance(m, pd.DataFrame) else pd.DataFrame(m, index=cm.index, columns=cm.columns)
            f.to_csv(folder / f'raw/pass_{k}.csv')
            passes.append(f.loc[[inv[x] for x in order], [inv[x] for x in order]].to_numpy(float))
        rel = None
        if len(passes) > 1:
            rel = float(np.max(np.abs(passes[-1] - passes[-2]) / np.maximum(np.abs(passes[-1]), 1e-30)))
        out = dict(units='fF', order=order, maxwell_fF=M.tolist(),
                   mutual_fF={f'{a}-{b}': float(-M[i, j]) for i, a in enumerate(order) for j, b in enumerate(order) if j > i},
                   ground_fF={a: float(M[i].sum()) for i, a in enumerate(order)},
                   net_object_names=lab, solver_converged=bool(lom.sim.is_converged), passes_saved=len(passes),
                   last_pass_max_relative_change=rel, minutes=(time.time() - t0) / 60,
                   aedt_design=lom.sim.renderer.pinfo.design.name)
        lom.sim.renderer.pinfo.project._project.Save()
        (folder / 'summary_out.json').write_text(json.dumps(out, indent=1))
        finish_run(folder, status='completed', summary=out)
        print(json.dumps(out, indent=1), flush=True)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
