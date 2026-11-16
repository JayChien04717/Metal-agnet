"""Build and check the 10 mm wafer-placement GDS candidates; preserve prior exports."""

import argparse
from datetime import date
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import gdstk

from gds_tools import die_edge_gaps, export_gds, plot_gds, write_wafer_gds
from project_paths import ROOT

sys.path.insert(0, str(ROOT / 'pattern design'))
from pcbport import pcb_ports


# Exact inputs for designs whose launch geometry already matches the PCB anchors.
REFERENCE_PAIRS = {
    'Design1': ('ftf_layout.gds', 'ftf_layout_fab.gds'),
    'Design3': ('FTF00_metal_only.gds', 'FTF00_metal_only_fab.gds'),
    'Design3_2': ('FTF00_metal_only.gds', 'FTF00_metal_only_fab.gds'),
    'Design3_3': ('Design3_3_2026-11-15_sim.gds', 'Design3_3_2026-11-15_fab.gds'),
}
DESIGNS = ('Design1', 'Design2', 'Design3', 'Design3_2', 'Design3_3', 'Design4', 'Design5', 'Design6')
SIDE_PORTS = {
    'IN': dict(pos_x='-4.5mm', pos_y='0.0mm', orientation='0'),
    'OUT': dict(pos_x='4.5mm', pos_y='0.0mm', orientation='180'),
}


def build_pair(name, staging):
    design = ROOT / name
    if name in REFERENCE_PAIRS:
        return [design / 'gds' / file for file in REFERENCE_PAIRS[name]]

    sim, fab = (staging / file for file in ('reference.gds', 'fab.gds'))
    sys.path.insert(0, str(design / 'layout'))
    if name == 'Design2':
        from layout3_full import export

        source = design / 'gds' / 'previous' / 'layout3 ftf.gds'
        export(source, sim, standard_launches=True)
        export(source, fab, fab=True, standard_launches=True)
    elif name == 'Design4':
        from purcell_6500 import LENGTHS, build

        for file, is_fab in ((sim, False), (fab, True)):
            export_gds(build(LENGTHS, wafer=True), str(file), fab=is_fab)
    elif name == 'Design5':
        from chip6_d4_fab import build_fab

        sim = design / 'gds' / 'Design5_2026-11-13_sim.gds'
        export_gds(build_fab(include_test_junctions=False), str(fab), fab=True)
    elif name == 'Design6':
        from design6_layout import export_masks

        export_masks(sim, fab, test_junctions=False, wafer_ports=True, port_options=pcb_ports)
    return sim, fab


def port_probes(port, width, gap):
    """A metal centre and both etched slots, on the pad, taper and neck."""
    x, y = (float(port[f'pos_{axis}'][:-2]) * 1000 for axis in ('x', 'y'))
    angle = int(port['orientation']) // 90
    dx, dy = ((1, 0), (0, 1), (-1, 0), (0, -1))[angle]
    taper_offset = ((width + gap) / 2 + 300) / 2  # midpoint of the slot halfway down the taper
    for radial, offset in ((-150, 210), (100, taper_offset), (265, (width + gap) / 2)):
        cx, cy = x + radial * dx, y + radial * dy
        yield [(cx, cy), (cx - offset * dy, cy + offset * dx),
               (cx + offset * dy, cy - offset * dx)]


def verify_output(file, name, fab, *, ports_override=None):
    lib = gdstk.read_gds(str(file))
    top, = lib.top_level()
    polygons = top.get_polygons()
    assert top.bounding_box() == ((-5000, -5000), (5000, 5000))
    assert {p.layer for p in polygons} <= ({1} if fab else {1, 50})
    metal = [p for p in polygons if p.layer == 1]
    assert all(p.datatype == 0 for p in metal)
    assert sum(p.area() for p in die_edge_gaps(metal)) < .01, (file, 'incomplete metal perimeter')
    ports = ports_override or (SIDE_PORTS if name in ('Design4', 'Design5') else pcb_ports)
    for label, port in ports.items():
        width, gap = ((15, 10) if name == 'Design2' else (16, 9)
                      if name in ('Design4', 'Design5') or (name == 'Design6' and label in ('P7', 'P8'))
                      else (10, 6))
        for points in port_probes(port, width, gap):
            assert gdstk.inside(points, metal) == (True, False, False), (file, label, points)
    if name == 'Design6':
        # The spare joins the ground at the inward end of the launch pocket.
        assert all(gdstk.inside([(-1100, -4210), (-1113, -4210), (-1087, -4210)], metal))
    if name == 'Design2' or (name == 'Design5' and fab):
        from shapely import get_parts
        from shapely.affinity import translate
        from xmark import label_geometry

        label = label_geometry('FTFC00' if name == 'Design2' else 'PFR00', 100, 10, 400, 30)
        left, _, right, _ = label.bounds
        lower_right = translate(label, 4650 - left - right, -4650)
        points = [part.representative_point().coords[0] for part in get_parts(lower_right)]
        assert not any(gdstk.inside(points, metal)), (file, 'blocked die ID')
        if name == 'Design5':
            upper_left = translate(label, -4650, 4650)
            points = [part.representative_point().coords[0] for part in get_parts(upper_left)]
            assert all(gdstk.inside(points, metal)), (file, 'old die ID remains')
    print(f'{file.name}: 10 × 10 mm outline; {len(ports)} launches checked', flush=True)


def export_design(name, tag, *, replace=False):
    home = ROOT / name
    stem = f'{name}_{tag}_wafer'
    file_types = ('reference', 'fab') if name in ('Design1', 'Design2', 'Design3', 'Design3_2') else ('sim', 'fab')
    output = [home / 'gds' / f'{stem}_{kind}.gds' for kind in file_types]
    figure = home / 'fig' / f'{stem}_gds.png'
    figures = {figure: dict(window=(-5.15, 5.15, -5.15, 5.15),
                            ports=SIDE_PORTS if name in ('Design4', 'Design5') else pcb_ports)}
    if name in ('Design2', 'Design5'):
        figures[home / 'fig' / f'{stem}_die_id.png'] = dict(window=(3.85, 5.055, -5.055, -4.40))
    if name == 'Design2':
        figures[home / 'fig' / f'{stem}_p5.png'] = dict(window=(.35, 5.05, -2.20, -1.05),
                                                      figsize=(10, 3.2))
    if not replace and any(path.exists() for path in [*output, *figures]):
        raise FileExistsError(f'{stem} already exists; choose an unused date or use --replace')
    with TemporaryDirectory(dir=home / 'gds') as temp:
        staging = Path(temp)
        staged = [staging / path.name for path in output]
        for source, target, fab in zip(build_pair(name, staging), staged, (False, True)):
            write_wafer_gds(source, target)
            verify_output(target, name, fab)
        for target, options in figures.items():
            plot_gds(str(staged[1]), str(staging / target.name), **options)
        for source, target in zip(staged, output):
            source.replace(target)
        for target in figures:
            (staging / target.name).replace(target)


def main():
    import geopandas as gpd

    # The existing GDS renderer still uses the older GeoDataFrame API.
    if not hasattr(gpd.GeoDataFrame, 'append'):
        gpd.GeoDataFrame.append = gpd.GeoDataFrame._append
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--date', required=True, type=date.fromisoformat)
    parser.add_argument('--design', choices=DESIGNS, nargs='+', default=DESIGNS)
    parser.add_argument('--replace', action='store_true', help='replace this date\'s wafer exports after validation')
    args = parser.parse_args()
    for name in args.design:
        export_design(name, args.date.isoformat(), replace=args.replace)


if __name__ == '__main__':
    main()
