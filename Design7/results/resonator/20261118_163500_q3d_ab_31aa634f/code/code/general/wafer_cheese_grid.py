"""Replace ground perforations on the finished fabrication geometry."""

import json
import os
from pathlib import Path
from shutil import which
import subprocess
from tempfile import TemporaryDirectory


HOLE_NM = 5000
PITCH_NM = 20000
PHASE_NM = 10000
WAFER_EDGE_UM = 547.5  # ±4452.5 µm, the Design3_3 reference's outer hole edges.
ETCH_CLEARANCE_NM = 50000
TILE_NM = 500000


def _encode(polygon):
    return ([(p.x, p.y) for p in polygon.each_point_hull()],
            [[(p.x, p.y) for p in polygon.each_point_hole(i)] for i in range(polygon.holes())])


def _decode(pya, data):
    hull, holes = data
    polygon = pya.Polygon([pya.Point(x, y) for x, y in hull])
    for ring in holes:
        polygon.insert_hole([pya.Point(x, y) for x, y in ring])
    return polygon


def _inspect(pya, source, record, spec):
    layout = pya.Layout()
    layout.read(source)
    if layout.dbu != .001:
        raise ValueError('Expected GDS on a 1 nm database grid')
    top = layout.top_cell()
    layers = {}
    for index in layout.layer_indices():
        info = layout.get_info(index)
        region = pya.Region(top.begin_shapes_rec(index)).merged()
        if not region.is_empty():
            layers[f'{info.layer}/{info.datatype}'] = list(region.each())
    shapes = layers[spec]
    main = max(shapes, key=lambda shape: shape.area())
    shapes.remove(main)
    hull, rings = _encode(main)
    filtered = []
    for ring in rings:
        bounds = pya.Polygon([pya.Point(x, y) for x, y in ring]).bbox()
        if bounds.width() > HOLE_NM + 10 or bounds.height() > HOLE_NM + 10:
            filtered.append(ring)
    data = {'ground': (hull, filtered), 'spec': spec, 'top': top.name,
            'old_holes': len(rings) - len(filtered),
            'layers': {key: [_encode(polygon) for polygon in parts] for key, parts in layers.items()}}
    Path(record).write_text(json.dumps(data), encoding='utf-8')


def _grid(ground, edge_um, launches):
    import numpy as np
    from shapely import STRtree, box, contains_xy, linestrings, prepare
    from wafer_cheese_keepout import LAUNCH_CLEARANCE_UM, _launch_frame

    margin = round(edge_um * 1000) + HOLE_NM // 2
    low_x, low_y, high_x, high_y = ground.bounds

    def axis(low, high):
        lo = int(np.ceil((low + margin - PHASE_NM) / PITCH_NM))
        hi = int(np.floor((high - margin - PHASE_NM) / PITCH_NM))
        return np.arange(lo, hi + 1, dtype=np.int64) * PITCH_NM + PHASE_NM

    grid_x, grid_y = np.meshgrid(axis(low_x, high_x), axis(low_y, high_y))
    x, y = grid_x.ravel(), grid_y.ravel()
    prepare(ground)
    allowed = contains_xy(ground, x, y)
    offset, width = np.array(LAUNCH_CLEARANCE_UM).T * 1000
    coords = (x, y)
    for options in launches.values():
        anchor, along, inward = _launch_frame(options)
        depth = inward * (coords[along] - anchor[along])
        across = 1 - along
        gap = abs(coords[across] - anchor[across]) - HOLE_NM // 2
        allowed &= ~((depth >= 0) & (depth <= offset[-1]) & (gap < np.interp(depth, offset, width)))
    x, y = x[allowed], y[allowed]

    segments = []
    for ring in [ground.exterior, *ground.interiors]:
        points = np.asarray(ring.coords)
        segments.extend(np.stack((points[:-1], points[1:]), axis=1))
    boundary = STRtree(linestrings(segments))
    centres = []
    half = HOLE_NM // 2
    for start in range(0, len(x), 20000):
        bx, by = x[start:start + 20000], y[start:start + 20000]
        windows = box(bx - half, by - half, bx + half, by + half)
        valid = np.ones(len(bx), dtype=bool)
        valid[boundary.query(windows, predicate='dwithin', distance=ETCH_CLEARANCE_NM - .001)[0]] = False
        centres.extend(np.stack((bx[valid], by[valid]), axis=1).tolist())
    return centres


def _write(pya, record, locations, destination, rounded, verify):
    from collections import defaultdict

    source = json.loads(Path(record).read_text(encoding='utf-8'))
    centres = json.loads(Path(locations).read_text(encoding='utf-8'))
    ground = pya.Region(_decode(pya, source['ground']))
    layout = pya.Layout()
    layout.dbu = .001
    top = layout.create_cell(source['top'])

    def layer(key):
        return layout.layer(*map(int, key.split('/')))

    for key, parts in source['layers'].items():
        for polygon in parts:
            top.shapes(layer(key)).insert(_decode(pya, polygon))

    holes_by_tile = defaultdict(list)
    for x, y in centres:
        holes_by_tile[x // TILE_NM, y // TILE_NM].append((x, y))
    h = HOLE_NM // 2
    template = pya.Polygon(pya.Box(-h, -h, h, h))
    if rounded:
        template = template.round_corners(2250, 2250, 32)

    ground_layer = layer(source['spec'])
    bounds = ground.bbox()
    for i in range(bounds.left // TILE_NM, (bounds.right - 1) // TILE_NM + 1):
        for j in range(bounds.bottom // TILE_NM, (bounds.top - 1) // TILE_NM + 1):
            tile = pya.Region(pya.Box(i * TILE_NM, j * TILE_NM, (i + 1) * TILE_NM, (j + 1) * TILE_NM))
            base = ground & tile
            cuts = pya.Region()
            for x, y in holes_by_tile[i, j]:
                cuts.insert(template.transformed(pya.Trans(x, y)))
            if not (cuts - base).is_empty():
                raise ValueError(f'Cheese opening crosses non-ground metal in tile {(i, j)}')
            if not base.is_empty():
                cell = layout.create_cell(f'G_{i}_{j}')
                cell.shapes(ground_layer).insert(base - cuts)
                top.insert(pya.CellInstArray(cell.cell_index(), pya.Trans()))

    options = pya.SaveLayoutOptions()
    options.gds2_max_vertex_count = 8190
    layout.write(destination, options)

    written = pya.Layout()
    written.read(destination)
    target_layer = written.layer(*map(int, source['spec'].split('/')))
    rounding_area = 0
    for i, j in holes_by_tile:
        cell = written.cell(f'G_{i}_{j}')
        tile = pya.Region(pya.Box(i * TILE_NM, j * TILE_NM, (i + 1) * TILE_NM, (j + 1) * TILE_NM))
        expected = ground & tile
        if cell is None:
            if not expected.is_empty() or holes_by_tile[i, j]:
                raise ValueError(f'Missing ground tile {(i, j)}')
            continue
        cuts = pya.Region()
        for x, y in holes_by_tile[i, j]:
            cuts.insert(template.transformed(pya.Trans(x, y)))
        actual = pya.Region(cell.begin_shapes_rec(target_layer))
        delta = actual ^ (expected - cuts)
        if not delta.is_empty():
            # GDS polygon fracture can round a split point on a sloped edge.
            # The openings themselves must read back exactly.
            area = delta.area()
            edge_band = expected.sized(3) - expected.sized(-3)
            if (delta - edge_band).area() or not delta.interacting(cuts).is_empty():
                raise ValueError(f'GDS geometry readback differs in tile {(i, j)} by {area} nm²')
            rounding_area += area
    Path(verify).write_text(json.dumps({'count': len(centres), 'rounding_area_nm2': rounding_area}), encoding='utf-8')


def replace_wafer_cheese(filename, *, launches=None, rounded=False, ground_spec=(1, 0), edge_um=WAFER_EDGE_UM):
    """Rebuild a 5 × 5 µm / 20 µm ground lattice and return verified counts."""
    from shapely.geometry import Polygon

    path = Path(filename).resolve()
    executable = os.environ.get('KLAYOUT_EXE') or which('klayout_app') or which('klayout')
    if not executable and os.name == 'nt':
        executable = Path.home() / 'AppData/Roaming/KLayout/klayout_app.exe'
    if not executable or not Path(executable).is_file():
        raise FileNotFoundError('KLayout executable not found. Set KLAYOUT_EXE to its path.')
    macro = Path(__file__).resolve()

    def klayout(**settings):
        command = [str(executable), '-b', '-r', str(macro)]
        for key, value in settings.items():
            command.extend(('-rd', f'{key}={value}'))
        result = subprocess.run(command, capture_output=True, text=True, errors='replace')
        if result.returncode:
            raise RuntimeError(f'KLayout {settings["mode"]} failed: {result.stderr or result.stdout}')

    with TemporaryDirectory(dir=path.parent) as tmp:
        folder = Path(tmp)
        sheet, grid, staged, verify = (folder / item for item in ('sheet.json', 'grid.json', path.name, 'verified.txt'))
        klayout(mode='inspect', source=path, record=sheet, spec='/'.join(map(str, ground_spec)))
        original = json.loads(sheet.read_text(encoding='utf-8'))
        centres = _grid(Polygon(*original['ground']), edge_um, launches or {})
        grid.write_text(json.dumps(centres), encoding='utf-8')
        klayout(mode='write', record=sheet, locations=grid, destination=staged, rounded=int(rounded), verify=verify)
        report = json.loads(verify.read_text(encoding='utf-8'))
        if report['count'] != len(centres):
            raise ValueError('Incorrect cheese count on GDS readback')
        staged.replace(path)
    return {'old_holes': original['old_holes'], 'new_holes': len(centres),
            'hole_um': 5, 'pitch_um': 20, 'edge_um': edge_um,
            'gds_rounding_area_um2': report['rounding_area_nm2'] / 1e6,
            'max_hole_edge_um': max(max(abs(x), abs(y)) for x, y in centres) / 1000 + 2.5 if centres else None}


if __name__ == '__main__':
    import pya

    if mode == 'inspect':
        _inspect(pya, source, record, spec)
    else:
        _write(pya, record, locations, destination, bool(int(rounded)), verify)
