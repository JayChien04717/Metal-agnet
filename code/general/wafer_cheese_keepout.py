"""Check wafer ground perforations and restore the local launch grid."""

import json
import os
from pathlib import Path
from shutil import which
import subprocess
from tempfile import TemporaryDirectory


# Sampled from the eight launch ports in the Design3_3 r5 fabrication mask.
# Each row is (inward offset from the PCB anchor, clearance from the centreline) in µm.
LAUNCH_CLEARANCE_UM = ((50, 397.5), (100, 297.5), (150, 247.5),
                       (200, 197.5), (250, 97.5), (300, 97.5))


def _read_klayout_voids(source, destination, grid_input=None, grid_output=None):
    # Executed by KLayout's embedded Python. A merged region avoids GDS fracture seams.
    import pya

    layout = pya.Layout()
    layout.read(source)
    if layout.dbu != .001:
        raise ValueError('Expected a 1 nm wafer GDS grid')
    metal = pya.Region(layout.top_cell().begin_shapes_rec(layout.layer(1, 0))).merged()
    if grid_input:
        ground = pya.Region(max(metal.each(), key=lambda polygon: polygon.area()))
        openings = pya.Region()
        for x, y in json.loads(Path(grid_input).read_text(encoding='utf-8')):
            openings.insert(pya.Box(x - 2500, y - 2500, x + 2500, y + 2500))
        centres = [polygon.bbox().center() for polygon in openings.inside(ground).each()]
        Path(grid_output).write_text(json.dumps([(p.x, p.y) for p in centres]), encoding='utf-8')
    die = pya.Region(pya.Box(-5000000, -5000000, 5000000, 5000000))
    with Path(destination).open('w', encoding='utf-8') as file:
        for polygon in (die - metal).merged().each():
            bounds = polygon.bbox()
            is_perforation = (bounds.width() <= 5010 and bounds.height() <= 5010
                              and polygon.area() > 100000)
            hull = [(p.x, p.y) for p in polygon.each_point_hull()]
            rings = [[(p.x, p.y) for p in polygon.each_point_hole(i)] for i in range(polygon.holes())]
            file.write(json.dumps([is_perforation, hull, rings]) + '\n')


def _launch_frame(port):
    angle = int(port['orientation']) % 360
    axis = (angle // 90) % 2
    inward = 1 if angle < 180 else -1
    anchor = [round(float(port[f'pos_{coord}'][:-2]) * 1e6) for coord in ('x', 'y')]
    return anchor, axis, inward


def _launch_grid(ports):
    """Sample the 50 µm lattice within 500 µm of each launch centreline."""
    centres = []
    for port in ports.values():
        anchor, axis, inward = _launch_frame(port)
        for depth, clearance in LAUNCH_CLEARANCE_UM:
            for offset in range(-500, 501, 50):
                if abs(offset) - 2.5 >= clearance:
                    point = anchor.copy()
                    point[axis] += inward * depth * 1000
                    point[1 - axis] += offset * 1000
                    centres.append(point)
    return centres


def _launch_fills(holes, ports):
    import numpy as np

    if not holes or not ports:
        return set()
    bounds = np.array([hole.bounds for hole in holes])
    low, high = bounds[:, :2], bounds[:, 2:]
    centres = (low + high) / 2
    distances, widths = np.array(LAUNCH_CLEARANCE_UM).T * 1000
    fills = set()
    for port in ports.values():
        anchor, axis, inward = _launch_frame(port)
        depth = inward * (centres[:, axis] - anchor[axis])
        across = 1 - axis
        gap = np.maximum(low[:, across] - anchor[across], anchor[across] - high[:, across])
        clearance = np.interp(depth, distances, widths)
        fills.update(np.flatnonzero((depth >= 0) & (depth <= distances[-1]) & (gap < clearance - .001)))
    return fills


def _perforation_fills(rows, launch_ports, available_centres, rounded_launch_holes):
    import numpy as np
    from shapely import STRtree, linestrings, set_precision
    from shapely.geometry import Polygon, box

    holes, etches = [], []
    for row in rows:
        is_perforation, hull, rings = json.loads(row)
        if is_perforation:
            holes.append(Polygon(hull, rings))
        else:
            for ring in [hull, *rings]:
                points = np.array(ring)
                etches.extend(np.stack((points, np.roll(points, -1, axis=0)), axis=1))

    allowed = box(-4500000, -4500000, 4500000, 4500000)
    edge = {i for i, hole in enumerate(holes) if not allowed.covers(hole)}
    # A metric query against boundary segments gives the edge-to-edge distance.
    close = set()
    boundaries = STRtree(linestrings(etches)) if etches else None
    if holes and boundaries is not None:
        close.update(boundaries.query(holes, predicate='dwithin', distance=49999.999)[0])
    launches = _launch_fills(holes, launch_ports)
    filled = edge | close | launches
    openings = [box(x - 2500, y - 2500, x + 2500, y + 2500) for x, y in available_centres]
    if rounded_launch_holes:
        openings = [set_precision(opening.buffer(-2250).buffer(2250, quad_segs=8), 1) for opening in openings]
    if openings:
        blocked = {i for i, opening in enumerate(openings) if not allowed.covers(opening)}
        if boundaries is not None:
            blocked.update(boundaries.query(openings, predicate='dwithin', distance=49999.999)[0])
        retained = [hole for i, hole in enumerate(holes) if i not in filled]
        if retained:
            blocked.update(STRtree(retained).query(openings, predicate='dwithin', distance=44999.999)[0])
        openings = [opening for i, opening in enumerate(openings) if i not in blocked]
    return [holes[i] for i in sorted(filled)], len(edge), len(close), len(launches), openings


def check_wafer_cheese_keepout(filename, *, repair=False, launch_ports=None,
                               restore_launch_grid=False, rounded_launch_holes=False):
    """Check wafer perforations and optionally restore eligible holes at launches."""
    import gdstk
    import numpy as np
    from shapely import STRtree
    from shapely.geometry import box

    path = Path(filename).resolve()
    executable = os.environ.get('KLAYOUT_EXE') or which('klayout_app') or which('klayout')
    if not executable and os.name == 'nt':
        executable = Path.home() / 'AppData/Roaming/KLayout/klayout_app.exe'
    if not executable or not Path(executable).is_file():
        raise FileNotFoundError('KLayout executable not found. Set KLAYOUT_EXE to its path.')

    with TemporaryDirectory(dir=path.parent) as folder:
        voids = Path(folder) / 'voids.jsonl'
        macro = Path(__file__).resolve()
        command = [str(executable), '-b', '-r', str(macro),
                   '-rd', f'source_gds={path}', '-rd', f'void_file={voids}']
        grid_input, grid_output = (Path(folder) / name for name in ('grid.json', 'available.json'))
        if restore_launch_grid and launch_ports:
            grid_input.write_text(json.dumps(_launch_grid(launch_ports)), encoding='utf-8')
            command += ['-rd', f'grid_input={grid_input}', '-rd', f'grid_output={grid_output}']
        result = subprocess.run(command, capture_output=True, text=True, errors='replace')
        if result.returncode or not voids.is_file() or (grid_input.is_file() and not grid_output.is_file()):
            raise RuntimeError(f'KLayout wafer keep-out scan failed: {result.stderr or result.stdout}')
        available = json.loads(grid_output.read_text(encoding='utf-8')) if grid_output.is_file() else []
        with voids.open(encoding='utf-8') as rows:
            fills, edge, close, launches, openings = _perforation_fills(
                rows, launch_ports, available, rounded_launch_holes)

        if repair and (fills or openings):
            library = gdstk.read_gds(str(path))
            top, = library.top_level()
            if openings:
                cuts = [gdstk.Polygon(np.asarray(opening.exterior.coords[:-1]) * .001) for opening in openings]
                lookup = STRtree([box(*low, *high) for low, high in (cut.bounding_box() for cut in cuts)])
                flattened = gdstk.Library(unit=library.unit, precision=library.precision)
                new_top = flattened.new_cell(top.name)
                for polygon in top.get_polygons():
                    low, high = polygon.bounding_box()
                    hits = lookup.query(box(*low, *high))
                    if len(hits) and (polygon.layer, polygon.datatype) == (1, 0):
                        new_top.add(*gdstk.boolean(polygon, [cuts[i] for i in hits], 'not', precision=.001, layer=1))
                    else:
                        new_top.add(polygon)
                library, top = flattened, new_top
            if any(hole.interiors for hole in fills):
                raise ValueError('Keep-out repair encountered a nested perforation')
            for hole in fills:
                top.add(gdstk.Polygon(np.asarray(hole.exterior.coords[:-1]) * .001, layer=1))
            staged = Path(folder) / path.name
            library.write_gds(str(staged), max_points=199)
            staged.replace(path)

    return {'holes': len(fills), 'edge_holes': edge, 'etch_holes': close,
            'launch_holes': launches, 'restorable_launch_holes': len(openings)}


if __name__ == '__main__':
    _read_klayout_voids(source_gds, void_file, globals().get('grid_input'), globals().get('grid_output'))
