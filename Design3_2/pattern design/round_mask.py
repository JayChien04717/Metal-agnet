"""Round the final positive-metal mask after resolving polygon seams."""
from pathlib import Path
from math import floor, ceil
import argparse

import gdstk
from shapely import get_parts, union_all, intersects_xy
from shapely.geometry import Polygon
from shapely.ops import clip_by_rect


import numpy as np
from shapely.geometry import Polygon

def rounded_ring(ring, radius, preserve=None):
    points = np.asarray(ring.coords)[:-1]
    incoming = np.roll(points, 1, axis=0) - points
    outgoing = np.roll(points, -1, axis=0) - points
    before = np.linalg.norm(incoming, axis=1)
    after = np.linalg.norm(outgoing, axis=1)
    incoming /= before[:, None]
    outgoing /= after[:, None]
    right_angle = np.abs(np.sum(incoming * outgoing, axis=1)) < 1e-8
    if preserve is not None:
        right_angle &= ~intersects_xy(preserve, points[:, 0], points[:, 1])
    # Reserve space for both fillets on short edges, including 2 um fab notches.
    short_edge = np.minimum(before, after)
    setback = np.where(short_edge >= 2 * radius - 1e-9, radius, .45 * short_edge)
    angles = np.linspace(0, np.pi/2, 9)
    result = []
    for point, u, v, r, is_corner in zip(points, incoming, outgoing, setback, right_angle):
        if is_corner and r >= .005:
            result.extend(point + r * ((1 - np.sin(angles[:, None])) * u
                                       + (1 - np.cos(angles[:, None])) * v))
        else:
            result.append(point)
    return result

def rounded_region(region, radius):
    region = region.simplify(0, preserve_topology=True)
    return Polygon(rounded_ring(region.exterior, radius),
                   [rounded_ring(ring, radius) for ring in region.interiors])

def round_mask(source, output, radius_um=5.0):
    """Apply convex and concave fillets to each connected metal region."""
    if radius_um <= 0:
        raise ValueError('Radius must be positive.')
    source_lib = gdstk.read_gds(str(source))
    if source_lib.unit != 1e-6:
        raise ValueError('Input GDS coordinates must be in micrometres.')
    top, = source_lib.top_level()
    polygons = top.get_polygons()
    if {p.layer for p in polygons} != {1}:
        raise ValueError('Expected a positive-metal mask on layer 1.')
    metal = union_all([Polygon(p.points).buffer(0) for p in polygons], grid_size=.001)
    original = list(get_parts(metal))
    rounded = []
    for conductor in original:
        region = rounded_region(conductor, radius_um)
        if not region.is_valid:
            raise ValueError('Rounded conductor boundary intersects itself.')
        if region.geom_type != 'Polygon' or region.is_empty:
            raise ValueError('Rounding broke or removed a metal region.')
        if len(region.interiors) != len(conductor.interiors):
            raise ValueError('Rounding changed the number of holes in a conductor.')
        rounded.append(region)
    if len(get_parts(union_all(rounded, grid_size=.001))) != len(original):
        raise ValueError('Rounding changed metal connectivity.')
    lib = gdstk.Library(unit=1e-6, precision=1e-9)
    cell = lib.new_cell('TOP')
    for region in rounded:
        xmin, ymin, xmax, ymax = region.bounds
        # Tile before GDS fracturing so the perforated ground is never passed
        # to the fracture algorithm as a single polygon with thousands of holes.
        for x in range(floor(xmin / 500), ceil(xmax / 500)):
            for y in range(floor(ymin / 500), ceil(ymax / 500)):
                clipped = clip_by_rect(region, x * 500, y * 500, (x + 1) * 500, (y + 1) * 500)
                for part in get_parts(clipped):
                    if part.geom_type != 'Polygon' or part.is_empty:
                        continue
                    outer = gdstk.Polygon(part.exterior.coords, layer=1)
                    holes = [gdstk.Polygon(ring.coords) for ring in part.interiors]
                    cell.add(*gdstk.boolean([outer], holes, 'not', precision=.001, layer=1))
    lib.write_gds(str(output), max_points=199)
    print(f'{output}: {len(original)} connected metal regions, '
          f'{sum(len(p.interiors) for p in rounded)} holes, R={radius_um:g} µm', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--radius-um', type=float, default=5.0)
    args = parser.parse_args()
    round_mask(args.source, args.output, args.radius_um)
