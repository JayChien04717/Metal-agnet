"""Solve the loaded Design6 r10 fixture with a simulation-only series input capacitor."""
from argparse import ArgumentParser, Namespace
import json
import math
import traceback

from shapely.geometry import LineString
from shapely.ops import substring

import readout_s21 as R
from readout_loaded_r9 import DESIGN, LAYOUT, read_json
from run_record import create_run, finish_run


def add_capacitor(design, capacitance_fF):
    """Split only the signal strip, retaining the baseline ground cut and both port pins."""
    feed = design.components['feed_out']  # Route from the P7 pad toward Q2.
    paths = design.qgeometry.tables['path']
    trace = paths[(paths.component == feed.id) & (paths.name == 'trace')].iloc[0]
    line = trace.geometry
    p7, corner, right, lower, tee = list(line.coords)
    gap_mm = 0.005
    x_mm = right[0] - trace.fillet - 0.03  # Straight, 30 µm past the bend tangent.
    center = math.dist(p7, corner) + x_mm - corner[0]
    start, stop = center - gap_mm / 2, center + gap_mm / 2
    before, after = substring(line, 0, start), substring(line, stop, line.length)
    sheet = LineString([before.coords[-1], after.coords[0]])
    paths.at[trace.name, 'geometry'] = before
    feed.add_qgeometry('path', {'trace_after_cin': after}, width=trace.width, fillet=trace.fillet)
    feed.add_qgeometry('junction', {'Cin': sheet}, width=trace.width,
                       hfss_capacitance=f'{capacitance_fF:g}fF', hfss_inductance=0, hfss_resistance=0)
    jj = design.qgeometry.tables['junction']
    cap = jj.loc[(jj.component == feed.id) & (jj.name == 'Cin')].iloc[0]
    assert cap.hfss_capacitance == f'{capacitance_fF:g}fF' and cap.hfss_inductance == cap.hfss_resistance == 0
    tee_center_x = design.parse_value(design.components['n0_Q2'].options.pos_x)
    arm_mm = line.length - center + abs(tee[0] - tee_center_x) - 2 * trace.fillet * (2 - math.pi / 2)
    return dict(Cin_fF=capacitance_fF, sheet_gap_um=gap_mm * 1000, strip_width_um=trace.width * 1000,
                sheet_wkt=sheet.wkt, arm_from_Q2_tee_mm=arm_mm, port_input='P7', port_output='P8')


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    import hammer_deembed as H

    wafer = read_json(DESIGN / 'design.json')['wafer_candidate']
    baseline = (DESIGN / wafer['readout_lengths_source']).parent
    source = read_json(baseline / 'run.json')
    settings = source['parameters']
    if source['status'] != 'completed' or settings['lengths_mm'] != wafer['readout_lengths_mm']:
        raise RuntimeError('The loaded baseline must match the current wafer readout lengths')
    if settings['junction_inductance'] != {k: v for k, v in R.L.LJ.items() if k != 'TC'}:
        raise RuntimeError('The loaded baseline junction inductances differ from the layout inputs')
    d = R.L.build_design(wafer_ports=True, readout_hammer=wafer['readout_hammer_um'],
                         lengths=wafer['readout_lengths_mm'], control_lines=False,
                         coupler=dict(flux_as_ground='True'))
    components = R.components(d) + ['Q1', 'Q2']
    geometry = json.loads(json.dumps(H.geometry_dump(d), default=str))
    if geometry != read_json(baseline / 'raw/geometry.json') or components != settings['components']:
        raise RuntimeError('The rebuilt fixture differs from the loaded HFSS baseline')
    cap = add_capacitor(d, 73.)
    solver = Namespace(**settings['HFSS'])
    solver.max_passes = 3
    solver.fine_band = [6.72, 7.08]
    solver.sparams = ['S11', 'S21', 'S12', 'S22']
    folder = create_run('Design6', 'resonator', 'readout_cin_r10' + ('_dry' if args.dry_run else ''),
                        entrypoint=__file__, upstream=[str(baseline)],
                        inputs=[DESIGN / 'design.json', LAYOUT, baseline / 'run.json',
                                baseline / 'raw/geometry.json', baseline / 'raw/ground_straps.json'],
                        parameters=dict(parent_geometry=wafer['revision'], candidate_only=True, capacitor=cap,
                                        fixture=settings, HFSS=vars(solver), dry_run=args.dry_run,
                                        ports={'1': 'P8', '2': 'P7'}, ground_straps='reused from baseline'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        H.write(folder / 'raw/geometry.json', H.geometry_dump(d))
        H.write(folder / 'raw/capacitor.json', cap)
        print(json.dumps(cap), flush=True)
        if args.dry_run:
            finish_run(folder, status='completed', summary=dict(capacitor=cap, components=components))
            return
        solver.ground_straps = read_json(baseline / 'raw/ground_straps.json')
        info = R.run_hfss(folder, d, components, solver)
        R.plot(folder, title='Design6 r10 loaded fixture with 73 fF series input capacitor')
        finish_run(folder, status='completed', summary=dict(capacitor=cap, **info))
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
