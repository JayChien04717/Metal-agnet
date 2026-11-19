"""Resonator-only HFSS eigenmode cases for the Design7 frequency audit."""

from argparse import ArgumentParser
from pathlib import Path
import traceback

import numpy as np
from shapely.geometry import Polygon, box

from readout_ab import HOME, L, P, apply_domain, check_design, create_run, finish_run, save_aedt, hfss_dirs, write
from readout_circuit import CASE, read
from readout_exchange import hammer_admittances, pole

CASES = ('isolated', 'grounded_hammer', 'plain_stub', 'no_tee', 'open_gap', 'short_flare')
RESONATOR = ('jr', 'res_open_run', 'readout_open', 'res_tail', 'res_tail_short')
SETUP = dict(min_freq_ghz=5.5, n_modes=2, max_delta_f=.1, max_passes=6, min_passes=2,
             min_converged=1, pct_refinement=30, basis_order=1)
MESH = '10um'
OPEN_GAP_UM = 45
FLARE_LENGTH_MM = .075
FLARE_WIDTH_MM = .028


def remove_filter_electrode(design):
    cj = design.components['cj'].id
    for kind, name in (('poly', 'head_a'), ('path', 'lead_a')):
        table = design.qgeometry.tables[kind]
        table.drop(table.index[(table.component == cj) & table.name.eq(name)], inplace=True)


def ground_counterelectrode(design):
    cj = design.components['cj']
    table = design.qgeometry.tables['poly']
    pocket, = table.loc[(table.component == cj.id) & table.name.eq('pocket'), 'geometry']
    pin = cj.pins['a'].middle
    from shapely.geometry import LineString

    cj.add_qgeometry('path', {'counterelectrode_short': LineString([pin, (pin[0], pocket.bounds[3])])},
                    width=cj.p.trace_width, layer=cj.p.layer)


def corner_in_place_of_tee(design):
    tee = design.components['jr']
    table = design.qgeometry.tables['poly']
    for row in table.loc[table.component == tee.id].itertuples():
        ceiling = tee.p.pos_y + tee.p.trace_width / 2 + (tee.p.trace_gap if row.subtract else 0)
        table.at[row.Index, 'geometry'] = row.geometry.intersection(box(-5, -5, 5, ceiling))


def build_case(settings, domain, parent_straps, case):
    design, lengths = L.build('opposing', settings)
    L.checks(design, lengths)
    apply_domain(design, domain)
    components = list(RESONATOR)
    if case in ('isolated', 'grounded_hammer'):
        components += ['stub_r', 'cj']
        if case == 'isolated':
            remove_filter_electrode(design)
        else:
            ground_counterelectrode(design)
    elif case == 'plain_stub':
        cj = design.components['cj'].p
        design.delete_component('stub_r')
        P.OpenToGround(design, 'stub_open', options=dict(
            **P.position(cj.pos_x, cj.pos_y), orientation='90',
            width=P.CPW['trace_width'], gap=P.CPW['trace_gap'], termination_gap=P.CPW['trace_gap']))
        P.straight(design, 'stub_r', ('jr', 'prime_end'), ('stub_open', 'open'))
        components += ['stub_r', 'stub_open']
    else:
        corner_in_place_of_tee(design)
        if case == 'open_gap':
            term = design.components['readout_open']
            term.options.termination_gap = f'{OPEN_GAP_UM}um'
            term.rebuild()
        if case == 'short_flare':
            term = design.components['res_tail_short']
            x, y = term.pins['end'].middle
            narrow, wide = term.p.trace_width / 2, FLARE_WIDTH_MM / 2
            flare = Polygon([(x - narrow, y - FLARE_LENGTH_MM), (x + narrow, y - FLARE_LENGTH_MM),
                             (x + wide, y), (x - wide, y)])
            term.add_qgeometry('poly', {'short_flare': flare}, layer=1)

    # Keep the driven fixture's resonator-side bridges. A removed stub also
    # removes its two ground-return bridges. Open and short controls share all bridges.
    straps = [s for s in parent_straps if s['where'].split('.')[0] in components]
    if 'stub_r' not in components:
        straps = [s for s in straps if s['where'] != 'jr.prime_end']
    ids = {design.components[name].id: name for name in components}
    metal = []
    for kind in ('path', 'poly'):
        for row in design.qgeometry.tables[kind].itertuples():
            if row.component not in ids:
                continue
            if row.geometry.is_empty or not row.geometry.is_valid:
                raise ValueError(f'Invalid {case} geometry: {row.name}')
            if not row.subtract:
                metal.append(f'{row.name}_{ids[row.component]}')
    return design, components, straps, metal, lengths['resonator']


def circuit_reference(case, params, length):
    if case in ('no_tee', 'open_gap', 'short_flare'):
        return params['velocity_m_s'] / (4 * length * 1e-3) / 1e9
    cj = params['cj_fF']
    if case == 'isolated':
        cj = [0, *cj[1:]]
    if case == 'plain_stub':
        cj = [0, 0, 0]
    p = dict(params, cj_fF=cj)
    return pole(lambda f: hammer_admittances(f, p)[1], 6.9).real


def solve(design, components, straps, metal, name, out, options):
    from pyEPR.ansys import parse_units
    from qiskit_metal import Dict
    from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim

    class ResonatorEigenmode(EigenmodeSim):
        def _analyze(self):
            for s in straps:
                span = parse_units(s['span'])
                self.renderer.modeler.draw_wirebond(
                    pos=parse_units(np.asarray(s['pos'])), ori=np.asarray(s['ori']), width=span,
                    height=span, z=0, wire_diameter=f"{s['diameter']}mm", NumSides=6,
                    name='g_wb', material='pec', solve_inside=False)
            self.renderer.modeler.mesh_length('signal_mesh', metal, MaxLength=MESH)
            self.sim_setup_name = self.renderer.initialize_eigenmode(**self.setup)
            self.renderer.analyze_setup(self.sim_setup_name)
            self.renderer.get_f_convergence([], save_csv=False)
            self.renderer.pinfo.design._reporter.ExportToFile('Freq. vs. pass', str(out / 'convergence_f.csv'))

    sim = ResonatorEigenmode(design, 'hfss')
    sim.setup.update(**options, vars=Dict())
    rendered, setup_name = sim.run_sim(name=name, components=components, open_terminations=[],
                                       box_plus_buffer=False)
    check_design(sim, rendered, name, 'hfss')
    setup = sim.renderer.pinfo.setup
    convergence, raw = setup.get_convergence()
    convergence.to_csv(out / 'convergence.csv')
    (out / 'convergence.txt').write_text(raw, encoding='utf-8')
    frequencies, kappa_over_2pi = setup.get_solutions().eigenmodes()
    if not any(6.5 < f < 7.5 for f in frequencies):
        raise ValueError(f'No resonator mode near 6.9 GHz: {frequencies}')
    candidates = [i for i, f in enumerate(frequencies) if 6.5 < f < 7.5]
    if len(candidates) != 1:
        raise ValueError(f'Ambiguous mode selection: {frequencies}')
    index, = candidates
    deltas = convergence.iloc[:, -1].dropna()
    sim.renderer.pinfo.project._project.Save()
    return dict(aedt_design=rendered, aedt_setup=setup_name, modes_GHz=list(map(float, frequencies)),
                mode_index_1based=index + 1, resonator_GHz=float(frequencies[index]),
                eigenmode_linewidth_GHz=None if kappa_over_2pi is None else list(map(float, kappa_over_2pi)),
                adaptive_passes=len(convergence), last_delta_f_percent=float(deltas.iloc[-1]) if len(deltas) else None)


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('geometry', type=Path)
    ap.add_argument('driven_run', type=Path)
    ap.add_argument('prediction_run', type=Path)
    ap.add_argument('--case', choices=CASES, action='append')
    ap.add_argument('--min-freq', type=float, default=SETUP['min_freq_ghz'])
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    if not 0 < args.min_freq < 6.5:
        ap.error('The eigenmode search floor must be positive and below the mode-selection band')
    options = dict(SETUP, min_freq_ghz=args.min_freq)
    settings = args.geometry.resolve()
    revision = read(settings)['geometry_revision']
    parent = args.driven_run.resolve()
    prediction = args.prediction_run.resolve()
    if {read(run / 'run.json')['geometry_revision'] for run in (parent, prediction)} != {revision}:
        raise ValueError('The model and driven reference must use the same geometry revision')
    boundary = read(parent / CASE / 'boundary.json')
    cases = args.case or CASES
    folder = create_run(
        'Design7', 'resonator', 'resonator_eigen_audit', entrypoint=__file__,
        inputs=[settings, parent / 'run.json', parent / CASE / 'boundary.json', prediction / 'run.json'],
        upstream=[str(parent), str(prediction)], parameters=dict(
            geometry_revision=revision, cases=cases, eigenmode=options, mesh=MESH, domain=boundary['domain'],
            open_gap_um=OPEN_GAP_UM, short_flare_width_mm=FLARE_WIDTH_MM,
            short_flare_length_mm=FLARE_LENGTH_MM, dry_run=args.dry_run,
            limitation='End controls measure sensitivity; shared offset includes propagation, bends and both ends.'))
    print(f'RUN_FOLDER={folder}', flush=True)
    record = read(folder / 'run.json')
    record['geometry_revision'] = revision
    if not args.dry_run:
        record['environment']['solver_version'] = save_aedt._desktop().get_version()
        record['environment']['aedt_temp'] = hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(HOME / 'simulation'))
    write(folder / 'run.json', record)
    params = read(prediction / 'run.json')['parameters']['circuit']
    results, errors = {}, {}
    from hammer_deembed import geometry_dump

    for case in cases:
        out = folder / 'raw' / case
        out.mkdir()
        print(f'CASE={case}', flush=True)
        try:
            design, comps, straps, metal, length = build_case(settings, boundary['domain'], boundary['straps'], case)
            geometry = geometry_dump(design)
            selected_ids = {design.components[n].id for n in comps}
            geometry['options'] = {k: v for k, v in geometry['options'].items() if k in comps}
            for table in ('path', 'poly'):
                geometry[table] = [row for row in geometry[table] if row['component'] in selected_ids]
            write(out / 'geometry.json', geometry)
            write(out / 'boundary.json', dict(components=comps, straps=straps, metal=metal,
                                              domain=boundary['domain'], ports=[]))
            result = dict(circuit_GHz=circuit_reference(case, params, length), length_mm=length)
            if not args.dry_run:
                result.update(solve(design, comps, straps, metal, f'D7_{folder.name[-8:]}_{case}', out, options))
                result['hfss_minus_circuit_MHz'] = 1000 * (result['resonator_GHz'] - result['circuit_GHz'])
            write(out / 'result.json', result)
            results[case] = result
            print(result, flush=True)
        except Exception:
            errors[case] = traceback.format_exc()
            (out / 'error.log').write_text(errors[case], encoding='utf-8')
            print(errors[case], flush=True)
    finish_run(folder, status='failed' if errors else 'completed', summary=dict(cases=results, errors=errors))
    if errors:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
