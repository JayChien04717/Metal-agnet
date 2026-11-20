"""Design4 Q3D: hammer C_j, C_p loads and the static filter-resonator parasitic for one stub variant.

Four local electrostatic solves (Q3D, 8 passes max, skill.md section 12), one AEDT design each, in
Design4/simulation/aedt/Design4.aedt:
  cj      module coupon (draw_mode='coupon'): two T regions + stubs + hammer, straight main lines +-400 um,
          etched open ends; nets filter / resonator
  ref     the same two straight main lines without stubs and hammer (deembedding reference)
  cp      Cai T-type EnclosingCap alone, open pins (frame = feed side, plate = filter side)
  module  the full PurcellReadout (short_ends='open', open pins): total static filter-resonator mutual C
Derived (summary.json):
  C_j          = mutual(cj) - mutual(ref)
  head_excess  = ground(cj) - ground(ref) - C' s / 2   per side (the model carries tap -> hammer centre as a line)
  parasitic_C  = mutual(module) - C_j
No ground bonds (finite local coupons, common ground; skill.md section 13).

    python -I -B Design4/simulation/resonator/q3d_module.py --variant ref|half [--cases cj,ref,cp,module] [--dry-run]
"""
from pathlib import Path
import argparse
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator'), str(ROOT / 'Design4' / 'layout')]
import project_paths  # noqa: E402
project_paths.configure(__file__)
from qiskit_metal import designs  # noqa: E402

import build as B  # noqa: E402
import purcell_readout_model as M  # noqa: E402
import purcell_resonator_pair_sim as S  # noqa: E402
from purcell_readout import PurcellReadout  # noqa: E402
from enclosing_cap import EnclosingCap  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

MAX_PASSES, PCT_ERROR = 8, 0.3


def new_design(size_x, size_y, cy=0.0):
    d = designs.DesignPlanar()
    d.overwrite_enabled = True
    d.chips.main.material = B.GEOMETRY['substrate']['material']
    d.chips.main.size.update(size_x=f'{size_x}mm', size_y=f'{size_y}mm', size_z='-500um',
                             center_x='0mm', center_y=f'{cy}mm')
    return d


def case_design(case, variant):
    opts = B.module_options(variant, B.model_lengths(variant))
    if case == 'cj':
        d = new_design(2.0, 2.0)
        PurcellReadout(d, 'cj', options=dict(opts, draw_mode='coupon', coupon_lead='400um'))
        return d, ['cj'], [], ('filter', 'resonator')
    if case == 'ref':
        d = new_design(2.0, 2.0)
        S.TwoLines(d, 'ref', options=dict(spacing=opts['line_spacing'], half_length='400um',
                                          trace_width=opts['trace_width'], trace_gap=opts['trace_gap']))
        return d, ['ref'], [], ('line_f', 'line_r')
    if case == 'cp':
        d = new_design(1.5, 1.0)
        EnclosingCap(d, 'cp', options=dict(B.cp_options(), trace_width=opts['trace_width'], trace_gap=opts['trace_gap']))
        return d, ['cp'], [('cp', 'a'), ('cp', 'b')], ('frame|lead_a', 'plate|lead_b')
    if case == 'module':
        d = new_design(4.6, 3.6, -0.05)
        PurcellReadout(d, 'purcell_resonator_pair', options=dict(opts, short_ends='open'))
        return d, ['purcell_resonator_pair'], [('purcell_resonator_pair', 'cp'), ('purcell_resonator_pair', 'qubit')], \
            ('filter', 'resonator')
    raise ValueError(case)


def derive(r, variant):
    cl = M.line(B.GEOMETRY['cpw']['model_preset'])['c_line']
    s_mm = B.spacing_um(variant) / 1000
    out = {}
    if 'cj' in r and 'ref' in r:
        out.update(S.hammer_loads(r['cj'], r['ref'], cl, s_mm))
        h = B.GEOMETRY['hammer']
        out['cj_fit_estimate'] = M.hammer_cj(h['head_gap_um'], h['head_height_um'])
    if 'cp' in r:
        out.update(cp_mutual=r['cp']['mutual_fF'], cp_plate_ground=r['cp']['ground_fF']['plate'],
                   cp_frame_ground=r['cp']['ground_fF']['frame'])
    if 'module' in r:
        out['module_mutual'] = r['module']['mutual_fF']
        if 'cj' in out:
            out['parasitic_C'] = r['module']['mutual_fF'] - out['cj']
            out['parasitic_fraction'] = out['parasitic_C'] / out['cj']
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--variant', required=True, choices=sorted(B.GEOMETRY['variants']))
    ap.add_argument('--cases', default='cj,ref,cp,module')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    cases = args.cases.split(',')
    v = args.variant
    params = dict(variant=v, line_spacing_um=B.spacing_um(v), stub_metal_um=B.stub_metal_um(v), cases=cases,
                  max_passes=MAX_PASSES, percent_error=PCT_ERROR, geometry_revision=B.GEOMETRY['revision'],
                  module_options=B.module_options(v, B.model_lengths(v)), cp_cap=B.GEOMETRY['cp_cap'],
                  substrate=B.GEOMETRY['substrate'], ground_bonds='none', short_ends='rendered open (module case)')
    folder = create_run('Design4', 'resonator', f'q3d_{v}', entrypoint=__file__,
                        inputs=[ROOT / 'Design4/parameters/geometry.json'], parameters=params)
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        for c in cases:
            d, *_ = case_design(c, v)
            for t in d.qgeometry.tables.values():
                if 'geometry' in t:
                    assert all(gm.is_valid and not gm.is_empty for gm in t.geometry), c
        if args.dry_run:
            finish_run(folder, status='completed', summary=dict(dry_run=True, geometry_checks='passed', solver_run=False))
            return
        import save_aedt
        import hfss_dirs
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(ROOT / 'Design4' / 'simulation'))
        rec = json.loads((folder / 'run.json').read_text(encoding='utf-8'))
        rec['environment']['solver_version'] = save_aedt._desktop().get_version()
        S.write(folder / 'run.json', rec)
        results = {}
        for c in cases:
            d, comps, opens, keys = case_design(c, v)
            results[c] = S.q3d_case(c, d, comps, opens, keys, folder, f'{v}_{folder.name[-8:]}',
                                    max_passes=MAX_PASSES, pct_error=PCT_ERROR, freq_ghz=B.GEOMETRY['f_target_GHz'])
        summary = dict(variant=v, results=results, derived=derive(results, v),
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
