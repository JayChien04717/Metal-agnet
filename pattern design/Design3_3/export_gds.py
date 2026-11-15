"""Export and verify the two final Design3_3 metal masks."""
import argparse
from datetime import date
from pathlib import Path
from tempfile import TemporaryDirectory

import gdstk
import matplotlib
from shapely import get_parts, union_all
from shapely.affinity import scale
from shapely.geometry import Polygon

from layout import build_design, baseline
from round_mask import round_mask

HERE = Path(__file__).resolve().parent


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--date', required=True, type=date.fromisoformat)
    args = parser.parse_args()
    stem = f'Design3_3_{args.date.isoformat()}'
    gds_dir = HERE / 'gds'
    fig_dir = HERE / 'fig'
    gds_dir.mkdir(exist_ok=True)
    fig_dir.mkdir(exist_ok=True)
    normal_name = f'{stem}_nocheese.gds'
    fab_name = f'{stem}.gds'
    with TemporaryDirectory(dir=HERE) as tmp:
        staging = Path(tmp)
        normal, normal_electrodes, junctions = make_mask(False)
        fab, fab_electrodes, _ = make_mask(True)
        for design, name, is_fab in ((normal, normal_name, False), (fab, fab_name, True)):
            raw = staging / f'raw_{name}'
            baseline.export_gds(design, str(raw), fab=is_fab)
            round_mask(raw, staging / name)
        validate_masks(read_metal(staging / normal_name), read_metal(staging / fab_name),
                       normal_electrodes, fab_electrodes, junctions)
        matplotlib.rcParams['patch.antialiased'] = False
        figures = {
            f'{stem}_layout.png': dict(ports=baseline.pcb_ports),
            f'{stem}_detail.png': dict(window=(-.65,.65,-1.35,-.45), dpi=140),
        }
        for filename, options in figures.items():
            baseline.plot_gds(str(staging / normal_name), str(staging / filename), **options)
        for name in (normal_name, fab_name):
            (staging / name).replace(gds_dir / name)
        for name in figures:
            (staging / name).replace(fig_dir / name)


if __name__ == '__main__':
    main()
