"""Compare mirrored XY termination offsets in a local r9 Q3D fixture."""
from argparse import ArgumentParser
import json
from pathlib import Path
import sys
import traceback

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(Path(__file__).parent), str(ROOT / 'code/general')]
import q3d_core_r9 as core
from run_record import create_run, finish_run

L = core.L
BASE = ROOT / 'Design6_TXT/results/qubit/20261117_161936_q3d_core_r9_3b00f1f5'
COMPS = core.COMPS + ['xy_start_Q1', 'xy_start_Q2']
OPENS = core.OPENS[:2]
CORNER_Y_MM = -0.6
CUT_X_MM = 1.7
SETUP = dict(core.SETUP, max_passes=3, percent_error=3.0)


def fixture(sep_um):
    cp = dict(flux_as_ground='True', squid_contact_overlap='0um', squid_shift_toward_stem='0um')
    d = L.build_design(channels=False, wafer_ports=True, control_lines=False, coupler=cp)
    d.chips.main.size.update(size_x='4mm', size_y='3mm')
    dw, dg = (d.parse_value(L.XMON[x]) for x in ('cross_width', 'cross_gap'))
    end_y = L.YQ - dw / 2 - dg - sep_um / 1000
    if end_y - CORNER_Y_MM <= d.parse_value(L.FILLET):
        raise ValueError('XY terminal segment must clear the r9 bend radius')
    pins = {}
    for q, side in (('Q1', -1), ('Q2', 1)):
        end_x = side * (L.QX + L.XY_END['dx'])
        pin_size = dict(width=L.CPW_CTL['trace_width'], gap=L.CPW_CTL['trace_gap'])
        start = L.OpenToGround(d, 'xy_start_' + q, options=dict(
            L.pos(side * CUT_X_MM, CORNER_Y_MM), orientation='180' if side < 0 else '0', **pin_size))
        end = L.OpenToGround(d, 'xy_end_' + q, options=dict(
            L.pos(end_x, end_y), orientation='90', **pin_size))
        route = L.route(d, 'xy_' + q, (start.name, 'open'), (end.name, 'open'),
                        [(end_x, CORNER_Y_MM)], L.CPW_CTL)
        assert np.allclose(route.pins['end'].middle, end.pins['open'].middle)
        assert np.dot(route.pins['end'].normal, end.pins['open'].normal) < -.99
        pins[q] = dict(start_mm=start.pins['open'].middle.tolist(), end_mm=end.pins['open'].middle.tolist())
    return d, pins


def metrics(matrix):
    r = core.QC.analyse(matrix, verbose=False)
    return dict(charge_lines=r['charge_lines'], xy_crosstalk=r['xy_crosstalk'], nodes=r['nodes'])


def compare_runs(sources):
    cases = {}
    for source in sources:
        for raw in sorted((source / 'raw').glob('*/cmat_fF.csv')):
            gap = int(raw.parent.name)
            if gap in cases:
                raise ValueError(f'Duplicate fixture for gap {gap} um')
            cases[gap] = raw
    inputs = [BASE / 'summary_out.json']
    for matrix in cases.values():
        inputs += [matrix, matrix.with_name('convergence.csv')]
        previous = matrix.with_name('pass_2_fF.csv')
        if previous.exists():
            inputs.append(previous)
    folder = create_run('Design6_TXT', 'qubit', 'xy_gap_compare_r9', entrypoint=__file__,
                        inputs=inputs, upstream=[str(BASE), *map(str, sources)],
                        parameters=dict(parent_geometry=core.REVISION, EM_rerun=False, gap_um=sorted(cases)))
    try:
        variants = {}
        for gap, matrix in sorted(cases.items()):
            matrix_data = pd.read_csv(matrix, index_col=0)
            np.testing.assert_allclose(matrix_data.to_numpy(), matrix_data.to_numpy().T, atol=1e-5)
            current = metrics(matrix_data)
            conv = pd.read_csv(matrix.with_name('convergence.csv'), index_col=0)
            previous = matrix.with_name('pass_2_fF.csv')
            prior = metrics(pd.read_csv(previous, index_col=0)) if previous.exists() else None
            drift = {}
            for q, charge in current['charge_lines'].items():
                key = 'C_line-qubit [fF]'
                drift[q] = None if prior is None else 100 * (charge[key] / prior['charge_lines'][q][key] - 1)
            variants[gap] = dict(final=current, last_delta_pct=float(conv.iloc[-1, -1]),
                                 pass_count=int(conv.index[-1]), pass2_to_last_C_change_pct=drift,
                                 source=str(matrix.relative_to(ROOT)))
        baseline = json.loads((BASE / 'summary_out.json').read_text(encoding='utf-8'))['final']
        combined_t1 = {q: 1 / (.01 + 1 / r['T1 limit [us]']) for q, r in baseline['charge_lines'].items()}
        summary = dict(variants=variants, accepted_r9=baseline, combined_with_intrinsic_100us_us=combined_t1)
        finish_run(folder, status='completed', summary=summary)
        print('RUN_FOLDER=' + str(folder), flush=True)
        for gap, report in variants.items():
            print(gap, report['last_delta_pct'], report['pass2_to_last_C_change_pct'],
                  report['final']['charge_lines'])
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--sep-um', type=int, nargs='+', default=[250, 350, 425])
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--compare', nargs='+', type=Path, help='Analyse existing fixture matrices without COM')
    args = parser.parse_args()
    if args.compare:
        compare_runs([path.resolve() for path in args.compare])
        return
    fixtures = {gap: fixture(gap) for gap in args.sep_um}
    geometry = {gap: pins for gap, (_, pins) in fixtures.items()}
    folder = create_run('Design6_TXT', 'qubit', 'q3d_xy_gap_r9', entrypoint=__file__,
                        inputs=[BASE / 'summary_out.json', ROOT / 'Design6_TXT/design.json'],
                        upstream=[str(BASE)], parameters=dict(
                            parent_geometry=core.REVISION, scope='local candidate fixture',
                            sep_um=args.sep_um, q3d=SETUP, chip_mm=[4, 3, 0.5],
                            sep_definition='vertical offset from the horizontal arm etch at the existing XY x',
                            geometry=geometry, components=COMPS, opens=OPENS,
                            ground_straps='none, matching core Q3D', dry_run=args.dry_run))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        if args.dry_run:
            finish_run(folder, status='completed', summary=geometry)
            print(json.dumps(geometry, indent=2))
            return
        from qiskit_metal.analyses.quantization import LOManalysis
        import hfss_dirs
        import save_aedt

        temp = hfss_dirs.use_d_drive()
        project = save_aedt.project_file(ROOT / 'Design6_TXT/simulation')
        save_aedt.use_project(project)
        manifest = folder / 'run.json'
        record = json.loads(manifest.read_text(encoding='utf-8'))
        record['environment'].update(solver_version=save_aedt._desktop().get_version(), solver_temp=temp['after'])
        record['parameters']['aedt_project'] = str(project)
        manifest.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')

        variants = {}
        for gap, (d, _) in fixtures.items():
            name = f'D6XY_{folder.name[-8:]}_{gap}'
            raw = folder / 'raw' / str(gap)
            raw.mkdir()
            lom = LOManalysis(d, 'q3d')
            lom.sim.setup.update(SETUP)
            lom.sim.run(name=name, components=COMPS, open_terminations=OPENS, box_plus_buffer=False)
            cmat = lom.sim.capacitance_matrix
            assert lom.sim.units.lower() == 'ff'
            cmat.to_csv(raw / 'cmat_fF.csv')
            conv, report = lom.sim.renderer.pinfo.setup.get_convergence()
            (raw / 'convergence.txt').write_text(report, encoding='utf-8')
            if conv is not None:
                conv.to_csv(raw / 'convergence.csv')
            by_pass = {}
            for k, array in lom.sim.capacitance_all_passes.items():
                matrix = pd.DataFrame(np.asarray(array) * 1e15, index=cmat.index, columns=cmat.columns)
                matrix.to_csv(raw / f'pass_{k}_fF.csv')
                try:
                    by_pass[k] = metrics(matrix)
                except ZeroDivisionError:
                    by_pass[k] = dict(error='Zero mutual capacitance in this coarse pass')
            variants[gap] = dict(final=metrics(cmat), by_pass=by_pass, solver_converged=lom.sim.is_converged,
                                 aedt_design=lom.sim.renderer.pinfo.design.name,
                                 delta_pct=None if conv is None else float(conv.iloc[-1, -1]))
            lom.sim.renderer.pinfo.project._project.Save()
            print(f'Gap {gap} um: {variants[gap]["final"]["charge_lines"]}', flush=True)
        finish_run(folder, status='completed', summary=dict(variants=variants))
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
