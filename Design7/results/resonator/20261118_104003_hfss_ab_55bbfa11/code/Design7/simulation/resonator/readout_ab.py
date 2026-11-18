"""Run the Design7 parallel/opposing A/B Q3D or driven HFSS study."""

from argparse import ArgumentParser
import json
from pathlib import Path
import sys
import traceback

import numpy as np
import pandas as pd

HOME = Path(__file__).resolve().parents[2]
ROOT = HOME.parent
sys.path[:0] = [str(HOME / 'layout'), str(ROOT / 'Design4/simulation/resonator'),
                str(ROOT / 'code/general')]
import design7_layout as L
import hfss_dirs
import save_aedt
from run_record import create_run, finish_run

P = L.P
SETTINGS = HOME / 'parameters/solver.json'
PORTS = [('feed_in', 'start'), ('feed_out', 'start')]
Q3D_OPENS = PORTS + [('filter_tail', 'end'), ('res_tail', 'end')]
NETS = tuple(L.NET_COMPONENTS)


def numpy_to_json(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f'Cannot serialize {type(value).__name__}')


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False, default=numpy_to_json) + '\n',
                    encoding='utf-8')


def metals(design):
    """Map renderer object names to electrical nodes."""
    groups = {name: net for net, names in L.NET_COMPONENTS.items() for name in names}
    objects = {}
    for kind in ('path', 'poly'):
        for row in design.qgeometry.tables[kind].itertuples():
            if row.subtract:
                continue
            component = design._components[row.component].name
            if component == 'cp':
                net = 'feed' if row.name.startswith(('frame', 'lead_a')) else 'filter'
            elif component == 'cj':
                net = 'filter' if row.name.startswith(('head_a', 'lead_a')) else 'resonator'
            else:
                net = groups[component]
            objects[f'{row.name}_{component}'] = net
    return objects


def check_port_planes(design):
    for (component, pin), terminal in zip(PORTS, ('input', 'output'), strict=True):
        route = design.components[component].pins[pin]
        reference = design.components[terminal].pins['end']
        if (not np.allclose(route.middle, reference.middle, atol=1e-7)
                or not np.allclose(route.normal, -reference.normal, atol=1e-7)):
            raise ValueError(f'{component}.{pin} is not the outward-facing {terminal} reference plane')


def ground_straps(design, variant):
    cpw = P.cpw_ground_straps(design, pitch_um=350)
    if variant == 'parallel':
        return P._drop_crossing(P.ground_straps(design, 'all') + cpw)

    from shapely.geometry import Point

    etch = P._etched_region(design)
    # Feed tee arm positions follow the existing Design4 policy. Opposing taps
    # need explicit straps on the lateral branch as well as both vertical arms.
    fixed = [s for s in P.junction_ground_straps(design) if s['where'].startswith('n0.')]
    span, diameter = .060, .015
    for name in ('jf', 'jr'):
        comp = design.components[name]
        origin = np.array([comp.p.pos_x, comp.p.pos_y])
        for pin, data in comp.pins.items():
            direction = data.normal
            pos = origin + direction * (.068 if pin == 'second_end' else .038)
            ortho = np.array([-direction[1], direction[0]])
            for foot in (pos - ortho * span / 2, pos + ortho * span / 2):
                if Point(*foot).buffer(diameter / 2).intersects(etch):
                    raise ValueError(f'Ground strap foot intersects {name}.{pin} etch')
            fixed.append(dict(pos=pos.tolist(), ori=ortho.tolist(), span=span,
                              diameter=diameter, where=f'{name}.{pin}'))
    return P._drop_crossing(fixed + cpw)


def check_design(sim, rendered, requested, suffix):
    expected = f'{requested}_{suffix}'
    active = sim.renderer.pinfo.design.name
    if active != expected or rendered != expected:
        raise RuntimeError(f'Expected {expected}, rendered {rendered}, active {active}')


def q3d(design, objects, folder, name, setup):
    from qiskit_metal.analyses.quantization import LOManalysis

    sim = LOManalysis(design, 'q3d').sim
    sim.setup.update(max_passes=setup['max_passes'], min_passes=setup['minimum_passes'],
                     percent_error=setup['percent_error_target'], freq_ghz=setup['adaptive_GHz'],
                     solution_order=setup['solution_order'], auto_increase_solution_order=False)
    rendered, setup_name = sim.run_sim(name=name, components=list(design.components),
                                       open_terminations=Q3D_OPENS, box_plus_buffer=False)
    check_design(sim, rendered, name, 'q3d')
    matrix = sim.capacitance_matrix
    matrix.to_csv(folder / 'matrix_fF.csv')
    convergence, convergence_text = sim.renderer.pinfo.setup.get_convergence()
    if convergence is not None:
        convergence.to_csv(folder / 'convergence.csv')
    (folder / 'convergence.txt').write_text(convergence_text, encoding='utf-8')
    if str(sim.units).lower() not in ('ff', 'femtofarad', 'femtofarads'):
        raise ValueError(f'Unexpected Q3D units: {sim.units}')
    signals = {obj: objects.get(obj) for obj in matrix.index if 'ground' not in obj.lower()}
    if (not any('ground' in obj.lower() for obj in matrix.index) or None in signals.values()
            or sorted(signals.values()) != sorted(NETS)):
        raise ValueError(f'Expected ground and three isolated signal nets, got {list(matrix.index)}')
    names = {v: k for k, v in signals.items()}
    order = [names[net] for net in NETS]
    sub = matrix.loc[order, order].to_numpy(float)
    cm = (sub + sub.T) / 2
    off_diagonal = np.triu(cm, 1)
    if (not np.allclose(sub, sub.T, atol=.01) or np.any(off_diagonal > 0)
            or np.any(np.linalg.eigvalsh(cm) <= 0) or np.any(cm.sum(axis=1) <= 0)):
        raise ValueError(f'Invalid signal Maxwell matrix: {sub}')
    passes = []
    for number, pass_matrix in sorted(sim.capacitance_all_passes.items(), key=lambda item: int(item[0])):
        frame = (pass_matrix if isinstance(pass_matrix, pd.DataFrame)
                 else pd.DataFrame(pass_matrix, index=matrix.index, columns=matrix.columns))
        frame.to_csv(folder / f'pass_{number}.csv')
        a, b = names['filter'], names['resonator']
        passes.append(abs(float((frame.loc[a, b] + frame.loc[b, a]) / 2)))
    last_rel = abs(passes[-2] / passes[-1] - 1) if len(passes) > 1 and passes[-1] else None
    return dict(aedt_design=rendered, aedt_setup=setup_name, units='fF', order=NETS,
                net_objects=names, maxwell_fF=cm.tolist(), ground_fF=dict(zip(NETS, cm.sum(axis=1))),
                mutual_fF={f'{a}-{b}': float(-cm[i, j]) for i, a in enumerate(NETS)
                           for j, b in enumerate(NETS) if j > i},
                solver_converged=sim.is_converged, saved_passes=len(passes),
                last_filter_res_mutual_relative_change=last_rel)


def hfss(design, objects, folder, name, setup, straps):
    from pyEPR.ansys import parse_units
    from qiskit_metal import Dict
    from qiskit_metal.analyses.simulation import ScatteringImpedanceSim

    class Driven(ScatteringImpedanceSim):
        def _analyze(self):
            for strap in straps:
                span = parse_units(strap['span'])
                self.renderer.modeler.draw_wirebond(
                    pos=parse_units(np.asarray(strap['pos'])), ori=np.asarray(strap['ori']),
                    width=span, height=span, z=0, wire_diameter=f"{strap['diameter']}mm",
                    NumSides=6, name='g_wb', material='pec', solve_inside=False)
            self.renderer.modeler.mesh_length('signal_mesh', list(objects),
                                              MaxLength=f"{setup['mesh_seed_um']}um")
            print(f'{name}: {len(straps)} ground straps and {len(objects)} seeded metal objects', flush=True)
            super()._analyze()

    sim = Driven(design, 'hfss')
    sim.setup.update(freq_ghz=setup['adaptive_GHz'], max_delta_s=setup['max_delta_s_target'],
                     max_passes=setup['max_passes'], min_passes=setup['minimum_passes'],
                     min_converged=1, pct_refinement=20, basis_order=1, vars=Dict())
    lo, hi = setup['wide_band_GHz']
    sim.setup.sweep_setup.update(name='Wide', start_ghz=lo, stop_ghz=hi,
                                 count=round((hi - lo) * 1000) + 1, type='Fast', save_fields=False)
    ports = [(*pin, setup['port_impedance_ohm']) for pin in PORTS]
    rendered, setup_name = sim.run_sim(name=name, components=list(design.components),
                                       port_list=ports, open_terminations=[], box_plus_buffer=False)
    check_design(sim, rendered, name, 'hfss')
    _, _, wide = sim.renderer.get_params(['S11', 'S21', 'S12', 'S22'])
    wide.to_csv(folder / 'S_wide.csv')
    conv = sim.renderer.pinfo.setup.get_convergence()[0]
    conv.to_csv(folder / 'convergence.csv')
    # Save the adaptive result and wide spectrum before requesting more samples.
    sim.renderer.pinfo.project._project.Save()
    lo, hi = setup['fine_band_GHz']
    sim.renderer.add_sweep(setup_name=setup_name, name='Fine', start_ghz=lo, stop_ghz=hi,
                           count=round((hi - lo) / setup['fine_step_GHz']) + 1,
                           type='Fast', save_fields=False)
    sim.renderer.analyze_sweep('Fine', setup_name)
    _, _, fine = sim.renderer.get_params(['S11', 'S21', 'S12', 'S22'])
    fine.to_csv(folder / 'S_fine.csv')
    after = sim.renderer.pinfo.setup.get_convergence()[0]
    after.to_csv(folder / 'convergence_after_fine.csv')
    return dict(aedt_design=rendered, aedt_setup=setup_name, straps=len(straps), wide_GHz=setup['wide_band_GHz'],
                fine_GHz=[lo, hi], wide_points=len(wide), fine_points=len(fine),
                adaptive_history_unchanged_by_fine=bool(conv.equals(after)))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('solver', choices=('q3d', 'hfss'))
    parser.add_argument('--variant', choices=L.VARIANTS, action='append')
    parser.add_argument('--controls', action='store_true', help='Add Q3D cases with only the hammer heads removed')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if args.controls and args.solver != 'q3d':
        parser.error('--controls applies only to Q3D')
    solver_key = 'q3d' if args.solver == 'q3d' else 'hfss_driven'
    setup = json.loads(SETTINGS.read_text(encoding='utf-8'))[solver_key]
    geometry = json.loads(L.SETTINGS.read_text(encoding='utf-8'))
    variants = args.variant or ('opposing',)
    states = (True, False) if args.controls else (True,)
    case_defs = [(variant, heads) for heads in states for variant in variants]
    folder = create_run(
        'Design7', 'resonator', f'{args.solver}_ab' + ('_dry' if args.dry_run else ''),
        entrypoint=__file__, inputs=[HOME / 'design.json', L.SETTINGS, SETTINGS],
        parameters=dict(geometry_revision=geometry['geometry_revision'], solver=args.solver,
                        setup=setup, cases=case_defs, open_terminations=Q3D_OPENS if args.solver == 'q3d' else [],
                        ports=PORTS if args.solver == 'hfss' else [],
                        domain='fixed 10 by 10 mm substrate, 0.5 mm silicon; renderer default vertical vacuum',
                        dry_run=args.dry_run))
    print(f'RUN_FOLDER={folder}', flush=True)
    manifest_path = folder / 'run.json'
    record = json.loads(manifest_path.read_text(encoding='utf-8'))
    record['geometry_revision'] = geometry['geometry_revision']
    results, errors = {}, {}
    try:
        if not args.dry_run:
            record['environment']['solver_version'] = save_aedt._desktop().get_version()
            record['environment']['aedt_temp'] = hfss_dirs.use_d_drive()
            save_aedt.use_project(save_aedt.project_file(HOME / 'simulation'))
        write(manifest_path, record)
        for variant, heads in case_defs:
            case = f'{variant}_{"full" if heads else "no_heads"}'
            dest = folder / 'raw' / case
            dest.mkdir()
            name = f'D7_{folder.name[-8:]}_{case}'
            print(f'CASE={case}', flush=True)
            try:
                run_case(args, variant, heads, dest, name, setup, results, case)
            except Exception:
                errors[case] = traceback.format_exc()
                (dest / 'error.log').write_text(errors[case], encoding='utf-8')
                print(errors[case], flush=True)
    except Exception:
        errors['initialization'] = traceback.format_exc()
        print(errors['initialization'], flush=True)
    summary = dict(dry_run=args.dry_run, cases=results, errors=errors)
    finish_run(folder, status='failed' if errors else 'completed', summary=summary)
    if errors:
        raise SystemExit(1)


def run_case(args, variant, heads, dest, name, setup, results, case):
    from hammer_deembed import geometry_dump

    design, measured = L.build(variant)
    layout_checks = L.checks(design, measured)
    check_port_planes(design)
    if not heads:
        table = design.qgeometry.tables['poly']
        rows = table.index[(table.component == design.components['cj'].id) &
                           table.name.isin(['head_a', 'head_b'])]
        if len(rows) != 2:
            raise ValueError('Expected exactly two hammer heads')
        table.drop(rows, inplace=True)
    objects = metals(design)
    straps = ground_straps(design, variant) if args.solver == 'hfss' else []
    write(dest / 'geometry.json', geometry_dump(design))
    write(dest / 'boundary.json', dict(straps=straps, metal_nets=objects, layout=layout_checks,
                                       substrate=design.chips.main.material))
    if args.dry_run:
        out = dict(layout_checks=layout_checks, straps=len(straps), metal_objects=len(objects))
    elif args.solver == 'q3d':
        out = q3d(design, objects, dest, name, setup)
    else:
        out = hfss(design, objects, dest, name, setup, straps)
    if not args.dry_run:
        save_aedt._desktop().get_active_project()._project.Save()
    write(dest / 'result.json', out)
    results[case] = out
    print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main()
