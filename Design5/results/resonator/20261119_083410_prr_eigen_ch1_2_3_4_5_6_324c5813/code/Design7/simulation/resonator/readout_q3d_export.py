"""Recover completed passes from the compact Q3D job after a clean stop."""

from argparse import ArgumentParser
from pathlib import Path
from shutil import copy2
import traceback

from readout_ab import HOME, NETS, create_run, finish_run, save_aedt, write
from readout_circuit import CASE, read


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('aedt_design')
    parser.add_argument('--pass-number', type=int, required=True)
    a = parser.parse_args()
    source = a.source.resolve()
    raw = source / CASE
    boundary = read(raw / 'boundary.json')
    revision = read(source / 'run.json')['geometry_revision']
    folder = create_run(
        'Design7', 'resonator', 'q3d_pass_export', entrypoint=__file__,
        inputs=[source / 'run.json', raw / 'geometry.json', raw / 'boundary.json', raw / 'error.log'],
        upstream=[str(source)], parameters=dict(geometry_revision=revision, EM_rerun=False,
                                               selected_pass=a.pass_number, aedt_design=a.aedt_design))
    print(f'RUN_FOLDER={folder}', flush=True)
    case = folder / CASE
    case.mkdir()
    for filename in ('geometry.json', 'boundary.json'):
        copy2(raw / filename, case / filename)

    try:
        desktop = save_aedt._desktop()
        project, = [p for p in desktop.get_projects() if p.name == HOME.name]
        active_design = project.get_design(a.aedt_design)
        setup = active_design.get_setup('Setup')
        record = read(folder / 'run.json')
        record['geometry_revision'] = revision
        record['environment']['solver_version'] = desktop.get_version()
        write(folder / 'run.json', record)
        conv, conv_text = setup.get_convergence()
        conv.to_csv(case / 'convergence.csv')
        (case / 'convergence.txt').write_text(conv_text, encoding='utf-8')

        for k in (a.pass_number - 1, a.pass_number):
            cm, units, _, _ = setup.get_matrix(pass_number=k, solution_kind='AdaptivePass')
            if str(units).lower() not in ('ff', 'femtofarad', 'femtofarads'):
                raise ValueError(f'Unexpected Q3D units: {units}')
            (cm * 1e-15).to_csv(case / f'pass_{k}.csv')
        cm.to_csv(case / 'matrix_fF.csv')
        names = {boundary['metal_nets'].get(obj): obj for obj in cm.index if 'ground' not in obj.lower()}
        if set(names) != set(NETS) or len(cm.index) != len(NETS) + 1:
            raise ValueError(f'Expected isolated feed, filter, resonator and ground: {list(cm.index)}')
        order = [names[net] for net in NETS]
        matrix = cm.loc[order, order].to_numpy(float)
        result = dict(units='fF', order=NETS, net_objects=names, selected_pass=a.pass_number,
                      aedt_design=a.aedt_design, maxwell_fF=matrix.tolist(),
                      ground_fF=dict(zip(NETS, matrix.sum(axis=1))),
                      mutual_fF={f'{a}-{b}': float(-matrix[i, j]) for i, a in enumerate(NETS)
                                 for j, b in enumerate(NETS) if j > i})
        project._project.Save()
        write(case / 'result.json', result)
    except Exception:
        finish_run(folder, status='failed', summary=dict(error=traceback.format_exc()))
        raise
    finish_run(folder, status='completed', summary=result)
    print(result, flush=True)


if __name__ == '__main__':
    main()
