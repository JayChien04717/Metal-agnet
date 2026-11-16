"""Compare local Design5 hammer couplings with fixed 400 µm CPW centres."""

import json
from pathlib import Path
import sys
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / path) for path in ('code/general', 'Design5/layout')]
from chip6_d4 import P, add_hammer, d5
from run_record import create_run, finish_run


def dimensions(hammer):
    half_span = 500 * P.MIN_CPW_SPACING
    return dict(**hammer, cpw_centres_um=2 * half_span,
                centre_to_outer_head_um=half_span - hammer['gap_um'] / 2 - hammer['width_um'],
                centre_to_inner_head_um=half_span - hammer['gap_um'] / 2)


def solve(hammer, destination, solver_name):
    from qiskit_metal.analyses.quantization import LOManalysis

    destination.mkdir()
    design = d5.designs.DesignPlanar()
    design.chips.main.material = d5.SI_MATERIAL
    design.chips.main.size.update(size_x='2mm', size_y='2mm', size_z='-500um',
                                  center_x='0mm', center_y='0mm')
    add_hammer(design, '', -.2, 0, True, hammer)
    geometry = {kind: [dict(name=row.name, component=int(row.component),
                          subtract=bool(row.subtract), wkt=row.geometry.wkt,
                          width_mm=getattr(row, 'width', None))
                      for row in design.qgeometry.tables[kind].itertuples()]
                for kind in ('poly', 'path')}
    (destination / 'geometry.json').write_text(json.dumps(geometry, indent=2), encoding='utf-8')

    sim = LOManalysis(design, 'q3d').sim
    sim.setup.update(max_passes=8, min_passes=2, percent_error=.3, freq_ghz=6.5,
                     solution_order='Medium', auto_increase_solution_order=False)
    opened = [(tee, end) for tee in ('jf', 'jr') for end in ('prime_start', 'prime_end')]
    sim.run(name=solver_name, open_terminations=opened, box_plus_buffer=False)
    matrix = sim.capacitance_matrix
    matrix.to_csv(destination / 'matrix_fF.csv')
    if str(sim.units).lower() not in ('ff', 'femtofarad', 'femtofarads'):
        raise ValueError(f'Unexpected Q3D units: {sim.units}')
    signals = [name for name in matrix.index if 'ground' not in name.lower()]
    if len(signals) != 2:
        raise RuntimeError(f'Expected two connected CPW electrodes, got {signals}')
    cm = matrix.loc[signals, signals].to_numpy(float)
    assert np.all(np.linalg.eigvalsh(cm) > 0) and cm[0, 1] < 0 and cm[1, 0] < 0
    assert np.allclose(cm, cm.T, atol=.01)
    passes = []
    for number, pass_matrix in sorted(sim.capacitance_all_passes.items(), key=lambda pair: int(pair[0])):
        pass_matrix.to_csv(destination / f'pass_{number}.csv')
        passes.append(abs(float(pass_matrix.loc[signals[0], signals[1]])))
    delta = abs(passes[-1] - passes[-2]) / passes[-1] if len(passes) >= 2 else None
    out = dict(geometry=dimensions(hammer), mutual_fF=float(-(cm[0, 1] + cm[1, 0]) / 2),
               shunt_fF=cm.sum(axis=1).tolist(), signals=signals, solver_converged=bool(sim.is_converged),
               last_mutual_relative_change=delta, adaptive_pass_count=len(passes))
    (destination / 'caps.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    print(json.dumps(out), flush=True)
    return out


def main():
    import hfss_dirs
    import save_aedt

    base = P.PAPER_HAMMER.copy()
    folder = create_run('Design5', 'resonator', 'hammer_quick_q3d', entrypoint=__file__,
                        inputs=[ROOT / 'Design5/design.json', ROOT / 'Design5/parameters/Design5_2026-11-13.json'],
                        parameters=dict(baseline=dimensions(base), domain_mm=[2, 2, .5], eps_r=11.5,
                                        cpw=P.CPW, prime_length_um=300, max_passes=8,
                                        percent_error=.3, frequency_ghz=6.5, ground_bridges=False,
                                        tolerance_relative=.02, user_reference='attached hammer illustration'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(ROOT / 'Design5/simulation'))
        record_path = folder / 'run.json'
        record = json.loads(record_path.read_text(encoding='utf-8'))
        record['environment']['solver_version'] = save_aedt._desktop().get_version()
        record['geometry_revision'] = 'Design5 local hammer fixed 400 um CPW spacing'
        record_path.write_text(json.dumps(record, indent=2), encoding='utf-8')
        rows = []

        def evaluate(hammer):
            label = f'case_{len(rows)}'
            row = solve(hammer, folder / 'raw' / label, f'D5CJ_{folder.name[-8:]}_{label}')
            row['raw_path'] = f'raw/{label}/caps.json'
            rows.append(row)
            return row

        reference = evaluate(base)
        target = reference['mutual_fF']
        # A wider head reduces the narrow arm, and a larger gap reduces the
        # CPW-to-inner-tip reach. Adjust head height to recover the baseline Cj.
        candidate = dict(base, width_um=50, gap_um=50)
        heights = [65.]
        for _ in range(4):
            height = heights[-1]
            row = evaluate(dict(candidate, height_um=height))
            error = row['mutual_fF'] / target - 1
            if abs(error) <= .02:
                break
            if len(heights) == 1:
                new_height = height * target / row['mutual_fF']
            else:
                prior = rows[-2]
                slope = (row['mutual_fF'] - prior['mutual_fF']) / (height - heights[-2])
                new_height = height + (target - row['mutual_fF']) / slope if slope > 0 else height * target / row['mutual_fF']
            new_height = round(float(np.clip(new_height, 30, 100)), 1)
            if new_height in heights:
                break
            heights.append(new_height)
        best = min(rows[1:], key=lambda row: abs(row['mutual_fF'] - target))
        error = best['mutual_fF'] / target - 1
        summary = dict(baseline=reference, selected=best, relative_delta=error, cases=rows,
                       meets_2_percent=abs(error) <= .02,
                       scope='Local electrostatic Cj comparison, with identical 300 um CPW sections and open ends.')
        save_aedt._desktop().get_active_project()._project.Save()
        finish_run(folder, status='completed', validation='local_check', summary=summary)
        print(f'Baseline {target:.4f} fF, best {best["mutual_fF"]:.4f} fF ({error:+.2%})', flush=True)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
