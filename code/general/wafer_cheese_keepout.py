"""Fill ground perforations that cross wafer fabrication keep-outs."""

import json
import os
from pathlib import Path
from shutil import which
import subprocess
from tempfile import TemporaryDirectory


def _read_klayout_voids(source, destination):
    # Executed by KLayout's embedded Python. A merged region avoids GDS fracture seams.
    import pya

    layout = pya.Layout()
    layout.read(source)
    if layout.dbu != .001:
        raise ValueError('Expected a 1 nm wafer GDS grid')
    metal = pya.Region(layout.top_cell().begin_shapes_rec(layout.layer(1, 0))).merged()
    die = pya.Region(pya.Box(-5000000, -5000000, 5000000, 5000000))
    with Path(destination).open('w', encoding='utf-8') as file:
        for polygon in (die - metal).merged().each():
            bounds = polygon.bbox()
            is_perforation = (bounds.width() <= 5010 and bounds.height() <= 5010
                              and polygon.area() > 100000)
            hull = [(p.x, p.y) for p in polygon.each_point_hull()]
            rings = [[(p.x, p.y) for p in polygon.each_point_hole(i)] for i in range(polygon.holes())]
            file.write(json.dumps([is_perforation, hull, rings]) + '\n')


def _perforation_fills(rows):
    import numpy as np
    from shapely import STRtree, linestrings
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
    if holes and etches:
        boundaries = STRtree(linestrings(etches))
        close.update(boundaries.query(holes, predicate='dwithin', distance=49999.999)[0])
    return [holes[i] for i in sorted(edge | close)], len(edge), len(close)


def check_wafer_cheese_keepout(filename, *, repair=False):
    """Measure keep-outs on 10 mm layer 1/0 metal and optionally fill offending holes."""
    import gdstk
    import numpy as np

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
        result = subprocess.run(command, capture_output=True, text=True, errors='replace')
        if result.returncode or not voids.is_file():
            raise RuntimeError(f'KLayout wafer keep-out scan failed: {result.stderr or result.stdout}')
        with voids.open(encoding='utf-8') as rows:
            fills, edge, close = _perforation_fills(rows)

        if repair and fills:
            library = gdstk.read_gds(str(path))
            top, = library.top_level()
            if any(hole.interiors for hole in fills):
                raise ValueError('Keep-out repair encountered a nested perforation')
            for hole in fills:
                top.add(gdstk.Polygon(np.asarray(hole.exterior.coords[:-1]) * .001, layer=1))
            staged = Path(folder) / path.name
            library.write_gds(str(staged), max_points=199)
            staged.replace(path)

    return {'holes': len(fills), 'edge_holes': edge, 'etch_holes': close}


if __name__ == '__main__':
    _read_klayout_voids(source_gds, void_file)
