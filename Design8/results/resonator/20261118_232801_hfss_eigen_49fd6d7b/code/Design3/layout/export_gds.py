"""Export and verify the Design3 simulation GDS and fabrication mask."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse
from datetime import date
from pathlib import Path
from tempfile import TemporaryDirectory

import gdstk
import matplotlib
import numpy as np
from shapely import get_parts, union_all
from shapely.affinity import scale
from shapely.geometry import Point, Polygon

from layout import build_design, baseline
from round_mask import round_mask

HERE = _P('Design3')


def make_mask(fab):
    design = build_design(rounded=False, jj_metal=False, include_tc_landing_pads=True,
                          include_test_junctions=False, fab_contacts=fab)
    assert not any(name.startswith('testjj_') for name in design.components)
    table = design.qgeometry.tables['poly']
    assert not table.name.str.contains('array_end|jj_stub').any()
    electrodes = {}
    for device, names in (('FX1', ('pad_top', 'pad_bot')), ('FX2', ('pad_top', 'pad_bot')),
                          ('TC', ('tmon', 'squid_ground'))):
        for name in names:
            shape = table.loc[table.component.eq(design.components[device].id)
                              & table.name.eq(name), 'geometry'].item()
            electrodes[device, name] = scale(shape, 1000, 1000, origin=(0, 0))
    junction_points = [
        scale(row.geometry.interpolate(.5, normalized=True), 1000, 1000, origin=(0, 0))
        for row in design.qgeometry.tables['junction'].itertuples()
        if row.name.startswith('rect_jj')
    ]
    return design, electrodes, junction_points


def export_simulation(filename, **layout_options):
    design = build_design(for_process=True, **layout_options)
    names = {component.id: name for name, component in design.components.items()}
    junctions = list(design.qgeometry.tables['junction'].itertuples())
    sites = {(names[row.component], row.name) for row in junctions}
    expected = {(device, name) for device in ('FX1', 'FX2')
                for name in ('rect_jj', 'array_top', 'array_bot')}
    expected.update({('TC', 'rect_jj1'), ('TC', 'rect_jj2')})
    assert sites == expected and len(junctions) == 8
    assert not any(name.startswith('testjj_') for name in design.components)

    baseline.export_gds(design, str(filename))
    lib = gdstk.read_gds(str(filename))
    assert lib.unit == 1e-6 and lib.precision == 1e-9
    top, = lib.top_level()
    metal = [p for p in top.get_polygons() if p.layer == 1]
    markers = []
    for row in junctions:
        coordinates = np.asarray(row.geometry.coords) * 1000
        marker, = gdstk.FlexPath(coordinates, row.width * 1000, layer=50).to_polygons()
        assert not any(gdstk.inside([coordinates.mean(axis=0)], metal)), (names[row.component], row.name)
        markers.append(marker)
    top.add(*markers)
    lib.write_gds(str(filename))
    top, = gdstk.read_gds(str(filename)).top_level()
    polygons = top.get_polygons()
    assert {p.layer for p in polygons} == {1, 50}
    assert sum((p.layer, p.datatype) == (50, 0) for p in polygons) == 8
    assert top.bounding_box() == ((-5000., -5000.), (5000., 5000.))
    print('Verified 4 JJ and 4 array annotations on simulation layer 50/0.')
    return design


def read_metal(filename):
    lib = gdstk.read_gds(str(filename))
    assert lib.unit == 1e-6 and lib.precision == 1e-9
    top, = lib.cells
    assert top.name == 'TOP' and not top.references and not top.labels
    polygons = top.get_polygons()
    assert {(p.layer, p.datatype) for p in polygons} == {(1, 0)}
    metal = union_all([Polygon(p.points).buffer(0) for p in polygons], grid_size=.001)
    assert metal.is_valid and metal.bounds == (-5000, -5000, 5000, 5000)
    return metal


def validate_masks(normal, fab, electrodes, fab_electrodes, junction_points):
    parts = [list(get_parts(metal)) for metal in (normal, fab)]
    assert len(parts[0]) == len(parts[1]), 'Fab cuts changed conductor count'
    removed = []
    for key, shape in electrodes.items():
        removed.extend(get_parts(shape.difference(fab_electrodes[key])))
    notches = [p for p in removed if p.area > 1]
    assert len(notches) == 8 and all(abs(p.area - 8) < .01 for p in notches)
    for notch in notches:
        point = notch.representative_point()
        assert normal.covers(point) and not fab.covers(point), 'Contact recess differs from requested mask'
    for metal, conductors in zip((normal, fab), parts):
        assert all(not metal.covers(point) for point in junction_points), 'Junction gap contains metal'
        ground = max(range(len(conductors)), key=lambda i: conductors[i].area)
        indices = {}
        for key, shape in fab_electrodes.items():
            probe = shape.representative_point()
            selected, = [i for i, part in enumerate(conductors) if part.covers(probe)]
            indices[key] = selected
        assert indices.pop(('TC', 'squid_ground')) == ground
        assert len(set(indices.values())) == 5 and ground not in indices.values()
    assert sum(len(p.interiors) for p in parts[1]) > sum(len(p.interiors) for p in parts[0])
    print(f'Both masks have {len(parts[0])} connected metal regions. Verified eight fab-only recesses.')


def check_flux_routes(design, metal):
    traces = design.qgeometry.tables['path']
    for name in ('fl_FX1', 'fl_FX2', 'fl_TC'):
        rows = traces.loc[traces.component.eq(design.components[name].id)]
        trace, = rows.loc[~rows['subtract']].itertuples()
        cut, = rows.loc[rows['subtract']].itertuples()
        coordinates = np.asarray(trace.geometry.coords) * 1000
        vectors = np.diff(coordinates, axis=0)
        lengths = np.linalg.norm(vectors, axis=1)
        index = int(np.argmax(lengths))
        midpoint = (coordinates[index] + coordinates[index + 1]) / 2
        dx, dy = vectors[index] / lengths[index]
        gap_center = 250 * (trace.width + cut.width)
        assert metal.covers(Point(midpoint)), (name, 'missing signal')
        for side in (-1, 1):
            probe = midpoint + side * gap_center * np.array([-dy, dx])
            assert not metal.covers(Point(probe)), (name, 'missing CPW slot')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--date', required=True, type=date.fromisoformat)
    args = parser.parse_args()
    stem = f'Design3_{args.date.isoformat()}'
    gds_dir = _P('Design3/gds')
    fig_dir = _P('Design3/fig')
    gds_dir.mkdir(exist_ok=True)
    fig_dir.mkdir(exist_ok=True)
    sim_name = f'{stem}_sim.gds'
    fab_name = f'{stem}_fab.gds'
    normal_name = 'unnotched_reference.gds'
    with TemporaryDirectory(dir=HERE) as tmp:
        staging = Path(tmp)
        normal, normal_electrodes, junctions = make_mask(False)
        fab, fab_electrodes, _ = make_mask(True)
        for design, name, is_fab in ((normal, normal_name, False), (fab, fab_name, True)):
            raw = staging / f'raw_{name}'
            baseline.export_gds(design, str(raw), fab=is_fab)
            round_mask(raw, staging / name)
        fab_metal = read_metal(staging / fab_name)
        validate_masks(read_metal(staging / normal_name), fab_metal,
                       normal_electrodes, fab_electrodes, junctions)
        check_flux_routes(fab, fab_metal)
        simulation = export_simulation(staging / sim_name)
        top, = gdstk.read_gds(str(staging / sim_name)).top_level()
        sim_metal = union_all([Polygon(p.points).buffer(0) for p in top.get_polygons()
                               if p.layer == 1], grid_size=.001)
        check_flux_routes(simulation, sim_metal)
        print('Verified signal and both CPW slots on all three flux routes in both GDS files.')
        matplotlib.rcParams['patch.antialiased'] = False
        figures = {
            f'{stem}_layout.png': dict(ports=baseline.pcb_ports),
            f'{stem}_detail.png': dict(window=(-.65,.65,-1.35,-.45), dpi=140),
        }
        for filename, options in figures.items():
            baseline.plot_gds(str(staging / fab_name), str(staging / filename), **options)
        for name in (sim_name, fab_name):
            (staging / name).replace(gds_dir / name)
        for name in figures:
            (staging / name).replace(fig_dir / name)


if __name__ == '__main__':
    main()
