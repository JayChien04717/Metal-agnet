"""Drive both readout channels with their claws and Xmons on an r9 length candidate."""
from argparse import ArgumentParser, Namespace
import json
from pathlib import Path
import traceback

import readout_s21 as R
from run_record import create_run, finish_run

ROOT = Path(__file__).resolve().parents[3]
DESIGN = ROOT / 'Design6'
LAYOUT = DESIGN / 'layout/design6_layout.py'
CHANNELS = ('Q1', 'Q2')
LENGTH_KEYS = ('L_filter', 'L_res')


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def check_trim_geometry(original, updated):
    tails = {f'{kind}_tail_{q}' for kind in ('filter', 'res') for q in CHANNELS}
    names = list(original['options'])
    if names != list(updated['options']):
        raise RuntimeError('Trim changed the list or order of Qiskit Metal components')
    ids = {names.index(name) + 1 for name in tails}
    for kind in ('poly', 'path'):
        if ([r for r in original[kind] if r['component'] not in ids]
                != [r for r in updated[kind] if r['component'] not in ids]):
            raise RuntimeError(f'Trim changed geometry outside the four tails ({kind})')


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('length_run', type=Path, help='Length model run or completed loaded HFSS run')
    parser.add_argument('--trim-um', type=float, nargs=4, metavar=('Q1_FILTER', 'Q1_RES', 'Q2_FILTER', 'Q2_RES'),
                        help='Signed offsets from a completed loaded HFSS run, in micrometres')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    path = args.length_run.resolve()
    source = read_json(path / 'run.json')
    if source['status'] != 'completed':
        parser.error('The parent run must be completed')
    wafer = read_json(DESIGN / 'design.json')['wafer_candidate']
    hammer = wafer['readout_hammer_um']
    junction_inductance = {k: v for k, v in R.L.LJ.items() if k != 'TC'}
    solver = Namespace(freq=6.9, wide=[6.55, 7.25], max_passes=3, fine_step_ghz=.0001,
                       open_terminations=[(q, 'flux_bias_line') for q in CHANNELS])
    trim = args.trim_um
    inputs = [path / 'run.json', path / 'summary.json', DESIGN / 'design.json', LAYOUT]
    if trim is None:
        lengths = read_json(path / 'summary.json')['lengths_fab_mm']
        prior_lengths = deltas = None
    else:
        parent = source['parameters']
        if (not parent.get('candidate_only') or parent.get('junction_inductance') != junction_inductance
                or parent.get('readout_hammer_um') != hammer or parent.get('parent_geometry') != wafer['revision']):
            parser.error('A trim requires the previous loaded candidate with the same junctions, hammer and wafer revision')
        prior_lengths = parent['lengths_mm']
        deltas = {q: dict(zip(LENGTH_KEYS, trim[2 * i:2 * i + 2])) for i, q in enumerate(CHANNELS)}
        lengths = {q: {key: prior_lengths[q][key] + deltas[q][key] / 1000 for key in LENGTH_KEYS}
                   for q in CHANNELS}
        inputs.append(path / 'raw/geometry.json')
    d = R.L.build_design(wafer_ports=True, readout_hammer=hammer, lengths=lengths,
                         control_lines=False, coupler=dict(flux_as_ground='True'))
    comps = R.components(d) + ['Q1', 'Q2']
    junctions = d.qgeometry.tables['junction']
    for q in CHANNELS:
        assert len(junctions[junctions['component'] == d.components[q].id]) == 2
    tag = 'readout_loaded_r9_trim' if trim is not None else 'readout_loaded_r9'
    folder = create_run('Design6', 'resonator', tag, entrypoint=__file__, inputs=inputs,
                        upstream=[str(path)], parameters=dict(
                            parent_geometry=wafer['revision'], candidate_only=True, lengths_mm=lengths,
                            previous_lengths_mm=prior_lengths, length_offsets_um=deltas,
                            components=comps, readout_hammer_um=hammer, added_claw_dl_mm=0,
                            junction_inductance=junction_inductance,
                            junction_capacitance_fF=0, omitted='TC and control routing',
                            HFSS=vars(solver), mesh_um=10, ports='50 ohm at the feed ends', dry_run=args.dry_run))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        import hammer_deembed as H
        geometry = H.geometry_dump(d)
        H.write(folder / 'raw/geometry.json', geometry)
        if trim is not None:
            check_trim_geometry(read_json(path / 'raw/geometry.json'), geometry)
            if comps != source['parameters']['components']:
                raise RuntimeError('Trim changed the HFSS component selection')
        print('LOADED_LENGTHS_MM=' + json.dumps(lengths), flush=True)
        if args.dry_run:
            finish_run(folder, status='completed', summary=dict(
                components=comps, open_terminations=solver.open_terminations, lengths_mm=lengths))
            return
        info = R.run_hfss(folder, d, comps, solver)
        R.plot(folder, title='Design6 r9 length candidate, HFSS driven with both claws and Xmons')
        finish_run(folder, status='completed', summary=info)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
