"""Solve the loaded Design6_TXT r10 fixture with a simulation-only series input capacitor."""
from argparse import ArgumentParser, Namespace
import json
import math
from pathlib import Path
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
    parser.add_argument('--retune-from', type=Path, help='Completed, untrimmed 73 fF HFSS run')
    parser.add_argument('--q2-filter-add-um', type=float, default=0., help='Signed Q2 filter length offset in µm')
    args = parser.parse_args()
    if not math.isfinite(args.q2_filter_add_um) or bool(args.q2_filter_add_um) != bool(args.retune_from):
        parser.error('Supply a finite nonzero Q2 filter offset and its untrimmed parent together')
    import hammer_deembed as H

    wafer = read_json(DESIGN / 'design.json')['wafer_candidate']
    baseline = (DESIGN / wafer['readout_lengths_source']).parent
    source = read_json(baseline / 'run.json')
    settings = source['parameters']
    if source['status'] != 'completed' or settings['lengths_mm'] != wafer['readout_lengths_mm']:
        raise RuntimeError('The loaded baseline must match the current wafer readout lengths')
    if settings['junction_inductance'] != {k: v for k, v in R.L.LJ.items() if k != 'TC'}:
        raise RuntimeError('The loaded baseline junction inductances differ from the layout inputs')
    lengths = {q: dict(values) for q, values in wafer['readout_lengths_mm'].items()}
    lengths['Q2']['L_filter'] += args.q2_filter_add_um / 1000
    d = R.L.build_design(wafer_ports=True, readout_hammer=wafer['readout_hammer_um'],
                         lengths=lengths, control_lines=False,
                         coupler=dict(flux_as_ground='True'))
    components = R.components(d) + ['Q1', 'Q2']
    geometry = json.loads(json.dumps(H.geometry_dump(d), default=str))
    if components != settings['components'] or (not args.retune_from and geometry != read_json(baseline / 'raw/geometry.json')):
        raise RuntimeError('The rebuilt fixture differs from the loaded HFSS baseline')
    cap = add_capacitor(d, 73.)
    upstream = [str(baseline)]
    inputs = [DESIGN / 'design.json', LAYOUT, baseline / 'run.json',
              baseline / 'raw/geometry.json', baseline / 'raw/ground_straps.json']
    straps_source = baseline
    tail_change = None
    if args.retune_from:
        parent = args.retune_from.resolve()
        previous = read_json(parent / 'run.json')
        if (previous['status'] != 'completed' or previous['upstream'] != [str(baseline)]
                or previous['parameters']['dry_run'] or previous['parameters'].get('q2_filter_add_um', 0)):
            raise RuntimeError('The retune parent must be a completed untrimmed capacitor solve for this baseline')
        if cap != read_json(parent / 'raw/capacitor.json'):
            raise RuntimeError('The filter retune changed the capacitor')
        original = read_json(parent / 'raw/geometry.json')
        current = json.loads(json.dumps(H.geometry_dump(d), default=str))
        tail = 'filter_tail_Q2'
        tail_id = d.components[tail].id
        if original['options'].keys() != current['options'].keys():
            raise RuntimeError('The filter retune changed the component names')
        for kind in ('poly', 'path'):
            old = [r for r in original[kind] if r['component'] != tail_id]
            new = [r for r in current[kind] if r['component'] != tail_id]
            if old != new:
                raise RuntimeError(f'The filter retune changed {kind} geometry outside {tail}')
        old_mm = d.parse_value(original['options'][tail]['_actual_length'])
        new_mm = d.parse_value(current['options'][tail]['_actual_length'])
        tail_change = 1000 * (new_mm - old_mm)
        if not math.isclose(tail_change, args.q2_filter_add_um, abs_tol=.01, rel_tol=0):
            raise RuntimeError(f'Expected the Q2 filter tail to change by {args.q2_filter_add_um:g} µm, got {tail_change:g}')
        upstream.append(str(parent))
        inputs.extend(parent / name for name in ('run.json', 'raw/geometry.json', 'raw/capacitor.json', 'raw/ground_straps.json'))
        straps_source = parent
    solver = Namespace(**settings['HFSS'])
    solver.max_passes = 3
    solver.fine_band = [6.72, 7.08]
    solver.sparams = ['S11', 'S21', 'S12', 'S22']
    tag = 'cin_q2trim' if args.retune_from else 'readout_cin_r10'
    folder = create_run('Design6_TXT', 'resonator', tag + ('_dry' if args.dry_run else ''),
                        entrypoint=__file__, upstream=upstream, inputs=inputs,
                        parameters=dict(parent_geometry=wafer['revision'], candidate_only=True, capacitor=cap,
                                        fixture=dict(settings, lengths_mm=lengths), HFSS=vars(solver), dry_run=args.dry_run,
                                        q2_filter_add_um=args.q2_filter_add_um, actual_tail_change_um=tail_change,
                                        ports={'1': 'P8', '2': 'P7'}, ground_straps_source=str(straps_source)))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        H.write(folder / 'raw/geometry.json', H.geometry_dump(d))
        H.write(folder / 'raw/capacitor.json', cap)
        print(json.dumps(cap), flush=True)
        if tail_change is not None:
            print(f'Q2_FILTER_TAIL_CHANGE_UM={tail_change:.6f}', flush=True)
        if args.dry_run:
            finish_run(folder, status='completed', summary=dict(capacitor=cap, components=components))
            return
        solver.ground_straps = read_json(straps_source / 'raw/ground_straps.json')
        info = R.run_hfss(folder, d, components, solver)
        R.plot(folder, title='Design6_TXT r10 loaded fixture with 73 fF series input capacitor')
        finish_run(folder, status='completed', summary=dict(capacitor=cap, **info))
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
