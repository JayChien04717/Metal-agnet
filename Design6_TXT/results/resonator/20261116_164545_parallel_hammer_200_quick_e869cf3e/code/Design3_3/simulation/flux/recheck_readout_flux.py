"""Recheck the final Design3_3 readout and flux fixtures in AEDT."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse
import os
from pathlib import Path
import traceback

import pandas as pd
from pyEPR.ansys import CalcObject
from shapely.geometry import LineString, box

from q3d import HERE, base
from layout import build_design

RESULTS = _P('Design3_3/results/general/legacy')
PROJECT = Path.home() / 'AppData/Local/Temp/D33Verify/D33Verify.aedt'
READOUTS = ('FX1', 'FX2')


def run_mutual(case):
    is_tc = case.startswith('tc')
    device = 'TC' if is_tc else case.upper()
    design = build_design(for_process=True, fx_options={'flux_as_ground': 'False'},
                          tmon_options={'flux_as_ground': 'False'})
    component = design.components[device]
    flux_pin = 'flux_bias_line' if is_tc else 'flux_line'
    junction = 'rect_jj'
    if case == 'fx1':
        pin = component.pins[flux_pin]
        x, y = pin.middle
        half = pin.width / 2
        table = design.qgeometry.tables['poly']
        i, = table.index[(table.component == component.id) & (table.name == 'flux_line')]
        leg = table.at[i, 'geometry']
        table.at[i, 'geometry'] = leg.difference(box(leg.bounds[0]-.001, y-half-.001, x, y+half+.001))
        assert table.at[i, 'geometry'].geom_type == 'Polygon'
    if is_tc:
        junction = 'rect_jj1' if case == 'tc_left' else 'rect_jj2'
        table = design.qgeometry.tables['junction']
        other, = table.loc[table.component.eq(component.id) & table.name.ne(junction)].itertuples()
        # Add 0.25 um contact to each end of the shorted junction.
        a, b = tuple(other.geometry.coords)
        length = other.geometry.length
        dx, dy = ((b[i]-a[i])*.00025/length for i in range(2))
        bridge = LineString([(a[0]-dx, a[1]-dy), (b[0]+dx, b[1]+dy)]).buffer(other.width/2, cap_style=2)
        component.add_qgeometry('poly', {'mutual_loop_closure': bridge})
        design.qgeometry.tables['junction'] = table.drop(index=other.Index)
    tag = 'm_' + case
    pins = [(device, pin) for pin in ('readout_line', 'charge_line') if pin in component.pins]
    base.mutual.RESULTS = _S('Design3_3/results/general/legacy')
    base.mutual.run_design(tag, design, device, flux_pin, junction, pins)
    desktop = base.save_aedt._desktop()._desktop
    messages = desktop.GetMessages(PROJECT.stem, f'{tag}_hfss', 0) or ()
    (_relocate(RESULTS / f'{tag}_messages.log')).write_text('\n'.join(messages), encoding='utf-8')


def run_kappa():
    solver = base.kappa
    solver.RESULTS = _S('Design3_3/results/general/legacy')
    solver.build_design = build_design
    solver.RESONATORS = {q: base.RESONATORS[q] for q in READOUTS}

    class NamedLOM(solver.LOManalysis):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            renderer = self.sim.renderer
            def analyze(name):
                s = renderer.pinfo.get_setup(name)
                s.parent._setup_module.EditSetup(name, [
                    f'NAME:{name}', 'AdaptiveFreq:=', '5GHz', 'SaveFields:=', False, 'Enabled:=', True,
                    ['NAME:Cap', 'MaxPass:=', 14, 'MinPass:=', 4, 'MinConvPass:=', 2,
                     'PerError:=', .2, 'PerRefine:=', 25, 'AutoIncreaseSolutionOrder:=', False,
                     'SolutionOrder:=', 'Medium', 'Solver Type:=', 'Iterative']])
                s.analyze(name)
            renderer.analyze_setup = analyze
            original = self.sim.run
            self.sim.run = lambda **kw: original(name='readout_kappa', **kw)

    solver.LOManalysis = NamedLOM
    solver.run('readout_kappa')
    project = base.save_aedt._desktop().get_active_project()
    setup = project.get_design('readout_kappa_q3d').get_setup('Setup')
    conv, _ = setup.get_convergence()
    conv.to_csv(_P('Design3_3/results/general/legacy/readout_kappa_convergence.csv'))
    for index in conv.index[-3:]:
        matrix, *_ = setup.get_matrix(pass_number=int(index), solution_kind='AdaptivePass')
        matrix.to_csv(_relocate(RESULTS / f'readout_kappa_pass{index}_cmat_fF.csv'))


def run_eigen(device):
    design = build_design(for_process=True)
    sim = base.EigenmodeSim(design, 'hfss')
    sim.setup.update(dict(base.SETUP, min_freq_ghz=5.75, n_modes=2, max_passes=15,
                          min_converged=2, pct_refinement=21))
    tee = f'tee_{device}'
    res = f'res_{device}'
    sim.run(name=f'readout_{device}', components=[res, tee, device],
            open_terminations=[(tee, 'prime_start'), (tee, 'prime_end')])
    sim.convergence_f.dropna(how='all').to_csv(_relocate(RESULTS / f'readout_{device}_convergence.csv'))
    sim.convergence_t.to_csv(_relocate(RESULTS / f'readout_{device}_mesh.csv'))
    base.save()
    setup = sim.renderer.pinfo.setup
    sol = setup.get_solutions()
    path = design.qgeometry.tables['path']
    trace = path.loc[(path.component == design.components[res].id) & (path.name == 'trace'), 'geometry'].item()
    xmin, ymin, xmax, ymax = trace.bounds
    cx = design.components[device].p.pos_x
    regions = {
        'resonator': (xmin-.04, xmax+.04, max(.35,ymin), ymax+.04),
        'qubit': (cx-.13, cx+.13, -1.03, -.62),
    }
    for key, (x0,x1,y0,y1) in regions.items():
        sim.renderer.modeler.draw_box_corner([x0*1e-3,y0*1e-3,-.12e-3],
            [(x1-x0)*1e-3,(y1-y0)*1e-3,.22e-3], name='E_'+key, nonmodel=True, wireframe=True)
    frequencies, _ = sol.eigenmodes()
    rows = []
    for index, frequency in enumerate(frequencies, 1):
        sol.set_mode(index)
        row = dict(mode=index, frequency_GHz=frequency)
        for key, (x0,x1,y0,y1) in regions.items():
            setup.parent._fields_calc.CalcStack('Clear')
            E = CalcObject([],setup).getQty('E')
            integral = E.dot(E.conj()).real().integrate_vol('E_'+key).evaluate()
            row[key+'_mean_E2'] = float(integral)/((x1-x0)*(y1-y0)*.22*1e-9)
        rows.append(row)
    pd.DataFrame(rows).to_csv(_relocate(RESULTS / f'readout_{device}_fields.csv'), index=False)
    print(rows, flush=True)


def main():
    global RESULTS, PROJECT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cases', nargs='+', choices=['fx1','fx2','tc_left','tc_right','kappa',*READOUTS])
    parser.add_argument("--results-dir", type=Path, default=RESULTS)
    parser.add_argument("--project", type=Path, default=PROJECT)
    args = parser.parse_args()
    RESULTS, PROJECT = args.results_dir.resolve(), args.project.resolve()
    RESULTS.mkdir(parents=True, exist_ok=True)
    os.chdir(RESULTS)
    base.save_aedt.use_project(PROJECT)
    errors = []
    for case in args.cases:
        print('Starting', case, flush=True)
        try:
            if case in READOUTS:
                run_eigen(case)
            elif case == 'kappa':
                run_kappa()
            else:
                run_mutual(case)
            base.save()
        except Exception:
            traceback.print_exc()
            errors.append(case)
    if errors:
        raise SystemExit('Failed cases: ' + ', '.join(errors))


if __name__ == '__main__':
    main()
