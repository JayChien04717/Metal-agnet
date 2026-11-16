"""Compare straight 190 µm hammer CPWs with the original 400 µm coupon."""

import json
from pathlib import Path
import sys
import traceback

import numpy as np
import pandas as pd
from shapely import get_parts, union_all
from shapely.geometry import box

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / path) for path in ('code/general', 'Design5/layout')]
from chip6_d4 import P, add_hammer, d5
from run_record import create_run, finish_run


def dimensions(hammer):
    half_span = hammer.get('cpw_centres_um', 1000 * P.MIN_CPW_SPACING) / 2
    return dict(hammer, cpw_centres_um=2 * half_span,
                centre_to_outer_head_um=half_span - hammer['gap_um'] / 2 - hammer['width_um'],
                centre_to_inner_head_um=half_span - hammer['gap_um'] / 2)


def verify_coupon(design, opened):
    metal, etches = [], []
    for kind in ('poly', 'path'):
        for row in design.qgeometry.tables[kind].itertuples():
            shape = row.geometry if kind == 'poly' else row.geometry.buffer(row.width / 2, cap_style=2)
            (etches if row.subtract else metal).append(shape)
    for component, name in opened:
        pin = design.components[component].pins[name]
        x, y = pin.middle
        half = pin.width / 2 + pin.gap
        outer_y = y + pin.normal[1] * pin.gap
        etches.append(box(x - half, y, x + half, outer_y))
    conductors = get_parts(union_all(metal, grid_size=.000001))
    ground = box(-1, -1, 1, 1).difference(union_all(etches, grid_size=.000001))
    if len(conductors) != 2 or any(ground.distance(net) < .005 for net in conductors):
        raise ValueError('The Q3D coupon must have two isolated conductors with open ends')


def coupon(hammer, prime_half_mm):
    design = d5.designs.DesignPlanar()
    design.chips.main.material = d5.SI_MATERIAL
    design.chips.main.size.update(size_x='2mm', size_y='2mm', size_z='-500um',
                                  center_x='0mm', center_y='0mm')
    left_x = -dimensions(hammer)['cpw_centres_um'] / 2000
    add_hammer(design, '', left_x, 0, True, hammer, prime_half_mm=prime_half_mm)
    geometry = {kind: [dict(name=row.name, component=int(row.component),
                          subtract=bool(row.subtract), wkt=row.geometry.wkt,
                          width_mm=getattr(row, 'width', None))
                      for row in design.qgeometry.tables[kind].itertuples()]
                for kind in ('poly', 'path')}
    opened = [(tee, end) for tee in ('jf', 'jr') for end in ('prime_start', 'prime_end')]
    verify_coupon(design, opened)
    return design, geometry, opened


def solve(hammer, destination, solver_name, prime_half_mm):
    from qiskit_metal.analyses.quantization import LOManalysis

    destination.mkdir()
    design, geometry, opened = coupon(hammer, prime_half_mm)
    (destination / 'geometry.json').write_text(json.dumps(geometry, indent=2), encoding='utf-8')
    sim = LOManalysis(design, 'q3d').sim
    sim.setup.update(max_passes=8, min_passes=2, percent_error=.3, freq_ghz=6.5,
                     solution_order='Medium', auto_increase_solution_order=False)
    # With components=None the renderer can silently solve the active AEDT design.
    rendered_name, setup_name = sim.run_sim(name=solver_name, components=list(design.components),
                                            open_terminations=opened, box_plus_buffer=False)
    active_name = sim.renderer.pinfo.design.name
    if rendered_name != f'{solver_name}_q3d' or active_name != rendered_name:
        raise RuntimeError(f'Q3D rendered {rendered_name!r} and analyzed {active_name!r}; expected {solver_name}_q3d')
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
        frame = (pass_matrix if isinstance(pass_matrix, pd.DataFrame)
                 else pd.DataFrame(pass_matrix, index=matrix.index, columns=matrix.columns))
        frame.to_csv(destination / f'pass_{number}.csv')
        passes.append(abs(float(frame.loc[signals[0], signals[1]])))
    delta = abs(passes[-1] - passes[-2]) / passes[-1] if len(passes) >= 2 else None
    out = dict(geometry=dimensions(hammer), mutual_fF=float(-(cm[0, 1] + cm[1, 0]) / 2),
               shunt_fF=cm.sum(axis=1).tolist(), signals=signals, solver_converged=bool(sim.is_converged),
               last_mutual_relative_change=delta, adaptive_pass_count=len(passes),
               aedt_design=rendered_name, aedt_setup=setup_name)
    (destination / 'caps.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    print(json.dumps(out), flush=True)
    return out


def main():
    import hfss_dirs
    import save_aedt

    base = {key: P.PAPER_HAMMER[key] for key in ('width_um', 'height_um', 'gap_um', 'ground_um')}
    candidate = dict(base, cpw_centres_um=190, width_um=25, gap_um=40)
    reference_run = ROOT / 'Design5/results/resonator/20261116_160806_hammer_190_q3d_070bbb59'
    reference_files = reference_run / 'raw/case_0'
    ref_manifest = json.loads((reference_run / 'run.json').read_text(encoding='utf-8'))
    prime_half_mm = ref_manifest['parameters']['projected_prime_length_um'] / 2000
    _, current, _ = coupon(base, prime_half_mm)
    baseline_geometry = json.loads((reference_files / 'geometry.json').read_text(encoding='utf-8'))
    if current != baseline_geometry:
        raise ValueError('The straight 400 µm baseline does not match its saved geometry')
    folder = create_run('Design5', 'resonator', 'hammer_parallel_190_q3d', entrypoint=__file__,
                        inputs=[ROOT / 'Design5/design.json', ROOT / 'Design5/parameters/Design5_2026-11-13.json',
                                reference_run / 'run.json', reference_files / 'geometry.json', reference_files / 'caps.json'],
                        upstream=[ref_manifest['run_id']],
                        parameters=dict(baseline=dimensions(base), domain_mm=[2, 2, .5], eps_r=11.5,
                                        cpw=P.CPW, projected_prime_length_um=2000 * prime_half_mm, max_passes=8,
                                        percent_error=.3, frequency_ghz=6.5, ground_bridges=False,
                                        tolerance_relative=.02, user_reference='attached hammer illustration'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(ROOT / 'Design5/simulation'))
        record_path = folder / 'run.json'
        record = json.loads(record_path.read_text(encoding='utf-8'))
        record['environment']['solver_version'] = save_aedt._desktop().get_version()
        record['geometry_revision'] = 'Design5 straight parallel 190 um hammer'
        record_path.write_text(json.dumps(record, indent=2), encoding='utf-8')
        rows = []

        def evaluate(hammer):
            label = f'case_{len(rows)}'
            row = solve(hammer, folder / 'raw' / label, f'D5CJ_{folder.name[-8:]}_{label}', prime_half_mm)
            row['raw_path'] = f'raw/{label}/caps.json'
            rows.append(row)
            return row

        reference = json.loads((reference_files / 'caps.json').read_text(encoding='utf-8'))
        reference['raw_path'] = str((reference_files / 'caps.json').relative_to(ROOT))
        rows.append(reference)
        target = reference['mutual_fF']
        # Retest the previous head dimensions on the straight geometry first.
        heights = [54.]
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
            new_height = round(float(np.clip(new_height, 10, 90)), 1)
            if new_height in heights:
                break
            heights.append(new_height)
        best = min(rows[1:], key=lambda row: abs(row['mutual_fF'] - target))
        error = best['mutual_fF'] / target - 1
        summary = dict(baseline=reference, selected=best, relative_delta=error, cases=rows,
                       meets_2_percent=abs(error) <= .02,
                       scope='Local electrostatic Cj comparison with equal straight CPW lengths and open ends.')
        save_aedt._desktop().get_active_project()._project.Save()
        finish_run(folder, status='completed', validation='local_check', summary=summary)
        print(f'Baseline {target:.4f} fF, best {best["mutual_fF"]:.4f} fF ({error:+.2%})', flush=True)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
