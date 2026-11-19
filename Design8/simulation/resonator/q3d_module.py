"""Design8 Q3D: hammer C_j, C_p loads and the static filter-resonator parasitic of the PurcellReadout module.

Four local electrostatic solves (Q3D, 8 passes max, skill.md section 12), one AEDT design each, in
Design8/simulation/aedt/Design8.aedt:
  cj      module coupon (draw_mode='coupon'): two T regions + stubs + hammer, straight main lines +-coupon_lead,
          etched open ends; nets filter / resonator
  ref     the same two straight main lines without T stubs and hammer (reference for the deembedding)
  cp      EnclosingCap C_p alone, open pins (frame = feed side, plate = filter side)
  module  the full PurcellReadout (short_ends='open', open pins): total static filter-resonator mutual C
Derived (summary.json):
  C_j            = mutual(cj) - mutual(ref)                 hammer-only coupling capacitance
  head_excess    = ground(cj) - ground(ref) - C' * s / 2     per side; the TL model carries the stub to the hammer
                   centre (s / 2) as a line, this is the rest (heads, pocket, T discontinuity)
  parasitic_C    = mutual(module) - C_j                      everything that is not the hammer
  cp_mutual, cp_plate_ground                                 filter open-end load = their sum
No ground bonds (finite local coupons, common ground; skill.md section 13 allows Q3D without bonds).

    python -I -B Design8/simulation/resonator/q3d_module.py [--cases cj,ref,cp,module] [--dry-run]
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'Design8' / 'layout')]
import project_paths  # noqa: E402
project_paths.configure(__file__)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from qiskit_metal import designs, draw, Dict  # noqa: E402
from qiskit_metal.qlibrary.core import QComponent  # noqa: E402

import build as B  # noqa: E402
import purcell_readout_model as M  # noqa: E402
from purcell_readout import PurcellReadout  # noqa: E402
from enclosing_cap import EnclosingCap  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

MAX_PASSES, PCT_ERROR = 8, 0.3


class TwoLines(QComponent):
    """Reference: two straight CPWs (x = +-spacing/2, y = +-half_length) with etched open ends; nets 'line_f', 'line_r'."""
    default_options = Dict(spacing='400um', half_length='400um', trace_width='10um', trace_gap='6um')

    def make(self):
        p = self.p
        for nm, x in (('line_f', -p.spacing / 2), ('line_r', p.spacing / 2)):
            ln = draw.LineString([(x, -p.half_length), (x, p.half_length)])
            self.add_qgeometry('poly', {nm: ln.buffer(p.trace_width / 2, cap_style=2)}, layer=1)
            cut = draw.LineString([(x, -p.half_length - p.trace_gap), (x, p.half_length + p.trace_gap)])
            self.add_qgeometry('poly', {nm + '_cut': cut.buffer(p.trace_width / 2 + p.trace_gap, cap_style=2)},
                               layer=1, subtract=True)


def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, default=float) + '\n', encoding='utf-8')


def new_design(size_x, size_y, cy=0.0):
    d = designs.DesignPlanar()
    d.overwrite_enabled = True
    d.chips.main.material = B.GEOMETRY['substrate']['material']
    d.chips.main.size.update(size_x=f'{size_x}mm', size_y=f'{size_y}mm', size_z='-500um',
                             center_x='0mm', center_y=f'{cy}mm')
    return d


def case_design(case):
    g = B.GEOMETRY
    opts = B.module_options(B.model_lengths('open'))
    s = g['module']['line_spacing_um'] / 1000
    if case == 'cj':
        d = new_design(2.0, 2.0)
        PurcellReadout(d, 'cj', options=dict(opts, draw_mode='coupon', coupon_lead='400um'))
        return d, ['cj'], [], ('filter', 'resonator')
    if case == 'ref':
        d = new_design(2.0, 2.0)
        TwoLines(d, 'ref', options=dict(spacing=opts['line_spacing'], half_length='400um',
                                        trace_width=opts['trace_width'], trace_gap=opts['trace_gap']))
        return d, ['ref'], [], ('line_f', 'line_r')
    if case == 'cp':
        d = new_design(1.5, 1.0)
        enc = {k: f'{v}um' for k, v in g['cp_cap'].items() if k != 'source'}
        EnclosingCap(d, 'cp', options=dict(enc, trace_width=opts['trace_width'], trace_gap=opts['trace_gap']))
        return d, ['cp'], [('cp', 'a'), ('cp', 'b')], ('frame|lead_a', 'plate|lead_b')
    if case == 'module':
        d = new_design(4.6, 3.6, -0.05)
        PurcellReadout(d, 'purcell_resonator_pair', options=dict(opts, short_ends='open'))
        return d, ['purcell_resonator_pair'], [('purcell_resonator_pair', 'cp'), ('purcell_resonator_pair', 'qubit')], ('filter', 'resonator')
    raise ValueError(case)


def geometry_dump(d):
    out = {}
    for kind in ('poly', 'path'):
        out[kind] = [dict(name=r.name, subtract=bool(r.subtract), wkt=r.geometry.wkt,
                          width=float(getattr(r, 'width', 0) or 0)) for r in d.qgeometry.tables[kind].itertuples()]
    return out


def solve(case, folder, run_tag):
    from qiskit_metal.analyses.quantization import LOManalysis
    d, comps, opens, keys = case_design(case)
    out = folder / 'raw' / case
    out.mkdir(parents=True, exist_ok=False)
    write(out / 'geometry.json', geometry_dump(d))
    lom = LOManalysis(d, 'q3d')
    lom.sim.renderer.options['x_buffer_width_mm'] = 0
    lom.sim.renderer.options['y_buffer_width_mm'] = 0
    lom.sim.setup.update(max_passes=MAX_PASSES, min_passes=2, percent_error=PCT_ERROR, freq_ghz=6.5,
                         solution_order='Medium', auto_increase_solution_order=False)
    t0 = time.time()
    lom.sim.run(name=f'D8Q_{case}_{run_tag}', components=comps, open_terminations=opens, box_plus_buffer=False)
    cm = lom.sim.capacitance_matrix
    cm.to_csv(out / 'matrix.csv')
    if str(lom.sim.units).lower() not in ('ff', 'femtofarad', 'femtofarads'):
        raise RuntimeError(f'units {lom.sim.units}')
    names = [n for n in cm.index if 'ground' not in n.lower()]
    signals = []
    for key in keys:
        hit = [n for n in names if any(k in n.lower() for k in key.split('|'))]
        if len(hit) != 1:
            raise RuntimeError(f'{case}: net for {key!r} not unique in {names}')
        signals.append(hit[0])
    if len(names) != 2:
        raise RuntimeError(f'{case}: expected 2 signal nets, got {names}')
    sub = cm.loc[signals, signals].to_numpy(float)
    passes = []
    for k, m in lom.sim.capacitance_all_passes.items():
        frame = m if isinstance(m, pd.DataFrame) else pd.DataFrame(m, index=cm.index, columns=cm.columns)
        frame.to_csv(out / f'pass_{k}.csv')
        passes.append(frame.loc[signals, signals].to_numpy(float))
    delta = None
    if len(passes) > 1:
        delta = float(np.max(np.abs(passes[-1] - passes[-2]) / np.maximum(np.abs(passes[-1]), 1e-30)))
    keys = tuple(k.split('|')[0] for k in keys)
    res = dict(case=case, signals=dict(zip(keys, signals)), mutual_fF=float(-sub[0, 1]),
               ground_fF=dict(zip(keys, sub.sum(axis=1).tolist())), maxwell_fF=sub.tolist(),
               max_passes=MAX_PASSES, percent_error=PCT_ERROR, passes_saved=len(passes),
               last_pass_max_relative_change=delta, solver_converged=bool(lom.sim.is_converged),
               minutes=(time.time() - t0) / 60, aedt_design=lom.sim.renderer.pinfo.design.name)
    write(out / 'caps.json', res)
    lom.sim.renderer.pinfo.project._project.Save()
    print(json.dumps(res), flush=True)
    return res


def derive(r):
    g = B.GEOMETRY
    cl = M.line(g['cpw']['model_preset'])['c_line']
    s_mm = g['module']['line_spacing_um'] / 1000
    out = {}
    if 'cj' in r and 'ref' in r:
        cj = r['cj']['mutual_fF'] - r['ref']['mutual_fF']
        ex = [r['cj']['ground_fF'][a] - r['ref']['ground_fF'][b] - cl * s_mm / 2
              for a, b in (('filter', 'line_f'), ('resonator', 'line_r'))]
        out.update(cj=cj, head_excess=float(np.mean(ex)), head_excess_sides=ex,
                   cj_fit_estimate=M.hammer_cj(g['hammer']['head_gap_um'], g['hammer']['head_height_um']),
                   ref_mutual_fF=r['ref']['mutual_fF'],
                   ref_c_line_fF_per_mm=r['ref']['ground_fF']['line_f'] / 0.8)
    if 'cp' in r:
        out.update(cp_mutual=r['cp']['mutual_fF'], cp_plate_ground=r['cp']['ground_fF']['plate'],
                   cp_frame_ground=r['cp']['ground_fF']['frame'])
    if 'module' in r:
        out.update(module_mutual=r['module']['mutual_fF'])
        if 'cj' in out:
            out['parasitic_C'] = r['module']['mutual_fF'] - out['cj']
            out['parasitic_fraction'] = out['parasitic_C'] / out['cj']
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cases', default='cj,ref,cp,module')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    cases = args.cases.split(',')
    params = dict(cases=cases, max_passes=MAX_PASSES, percent_error=PCT_ERROR, geometry_revision=B.GEOMETRY['revision'],
                  module_options=B.module_options(B.model_lengths('open')), cp_cap=B.GEOMETRY['cp_cap'],
                  substrate=B.GEOMETRY['substrate'], ground_bonds='none', short_ends='rendered open (module case)')
    folder = create_run('Design8', 'resonator', 'q3d_module', entrypoint=__file__,
                        inputs=[ROOT / 'Design8/parameters/geometry.json'], parameters=params)
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        for c in cases:
            d, *_ = case_design(c)
            for t in d.qgeometry.tables.values():
                if 'geometry' in t:
                    assert all(gm.is_valid and not gm.is_empty for gm in t.geometry), c
        if args.dry_run:
            finish_run(folder, status='completed', summary=dict(dry_run=True, geometry_checks='passed', solver_run=False))
            return
        import save_aedt
        import hfss_dirs
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(ROOT / 'Design8' / 'simulation'))
        rec = json.loads((folder / 'run.json').read_text(encoding='utf-8'))
        rec['environment']['solver_version'] = save_aedt._desktop().get_version()
        rec['geometry_revision'] = B.GEOMETRY['revision']
        write(folder / 'run.json', rec)
        results = {}
        for c in cases:
            results[c] = solve(c, folder, folder.name[-8:])
        summary = dict(results=results, derived=derive(results),
                       validation='Q3D 8-pass budget; per-case solver_converged and last-pass change recorded; '
                                  'pass count is not convergence evidence')
        print(json.dumps(summary['derived'], indent=1), flush=True)
        finish_run(folder, status='completed', summary=summary)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
