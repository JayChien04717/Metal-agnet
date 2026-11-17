"""Drive both readout channels with their claws and Xmons on an r9 length candidate."""
from argparse import ArgumentParser, Namespace
import json
from pathlib import Path
import sys
import traceback

import readout_s21 as R
from run_record import create_run, finish_run

ROOT = Path(__file__).resolve().parents[3]
DESIGN = ROOT / 'Design6'
LAYOUT = DESIGN / 'layout/design6_layout.py'


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('length_run', type=Path)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    path = args.length_run.resolve()
    lengths = json.loads((path / 'summary.json').read_text(encoding='utf-8'))['lengths_fab_mm']
    wafer = json.loads((DESIGN / 'design.json').read_text(encoding='utf-8'))['wafer_candidate']
    hammer = wafer['readout_hammer_um']
    d = R.L.build_design(wafer_ports=True, readout_hammer=hammer, lengths=lengths,
                         control_lines=False, coupler=dict(flux_as_ground='True'))
    comps = R.components(d) + ['Q1', 'Q2']
    opens = [(q, 'flux_bias_line') for q in ('Q1', 'Q2')]
    junctions = d.qgeometry.tables['junction']
    for q in ('Q1', 'Q2'):
        assert len(junctions[junctions['component'] == d.components[q].id]) == 2
    solver = Namespace(freq=6.9, wide=[6.55, 7.25], max_passes=3, fine_step_ghz=.0001, open_terminations=opens)
    folder = create_run('Design6', 'resonator', 'readout_loaded_r9', entrypoint=__file__,
                        inputs=[path / 'summary.json', DESIGN / 'design.json', LAYOUT],
                        upstream=[str(path)], parameters=dict(
                            parent_geometry=wafer['revision'], candidate_only=True, lengths_mm=lengths,
                            components=comps, readout_hammer_um=hammer, added_claw_dl_mm=0,
                            junction_inductance={k: v for k, v in R.L.LJ.items() if k != 'TC'},
                            junction_capacitance_fF=0, omitted='TC and control routing',
                            HFSS=vars(solver), mesh_um=10, ports='50 ohm at the feed ends', dry_run=args.dry_run))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        import hammer_deembed as H
        H.write(folder / 'raw/geometry.json', H.geometry_dump(d))
        if args.dry_run:
            finish_run(folder, status='completed', summary=dict(components=comps, open_terminations=opens))
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
