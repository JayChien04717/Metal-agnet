"""Replace wafer corner marks with the DesignTestJunction fabrication coupons."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import gdstk
from shapely import STRtree, get_parts
from shapely.affinity import scale, translate
from shapely.geometry import Point, Polygon, box

from export_wafer_masks import verify_output
from gds_tools import plot_gds
from project_paths import ROOT

sys.path.insert(0, str(ROOT / 'DesignTestJunction' / 'layout'))
from build import build_design, GEOMETRY, read_region
from xmark import CORNERS


DESIGNS = ('Design2', 'Design3_3', 'Design6')
GRID = .001  # µm
AREA_TOL = .2  # µm², allows for vertex rounding on the 1 nm grid


def polygons(shape):
    """Convert a region in µm to positive metal, preserving holes."""
    result = []
    for part in get_parts(shape):
        if part.geom_type != 'Polygon' or not part.area:
            continue
        shell = gdstk.Polygon(part.exterior.coords[:-1], layer=1)
        holes = [gdstk.Polygon(ring.coords[:-1]) for ring in part.interiors]
        result.extend(gdstk.boolean(shell, holes, 'not', layer=1, precision=GRID) if holes else [shell])
    return result


def region_in(window, metal, index):
    fragments = [metal[i] for i in index.query(window)]
    clip = gdstk.rectangle(window.bounds[:2], window.bounds[2:])
    return read_region(gdstk.boolean(fragments, clip, 'and', precision=GRID))


def coupon_shapes(variant, x, y):
    design = build_design(variant, fab_contacts=True)
    table = design.qgeometry.tables['poly']
    placed = {row.name: translate(scale(row.geometry, 1000, 1000, origin=(0, 0)), x, y)
              for row in table.itertuples()}
    return placed['pocket'], placed['letter'], (placed['pad_l'], placed['pad_r'])


def add_corner_coupons(source, output, placements):
    lib = gdstk.read_gds(str(source))
    if (lib.unit, lib.precision) != (1e-6, 1e-9):
        raise ValueError('Expected parent wafer GDS with µm units and 1 nm precision')
    top, = lib.top_level()
    metal = top.get_polygons()
    assert {(p.layer, p.datatype) for p in metal} == {(1, 0)}
    index = STRtree([box(*low, *high) for low, high in (p.bounding_box() for p in metal)])

    cuts, additions, checks = [], [], []
    for corner, sx, sy in CORNERS:
        x, y = placements[corner]['centre_um']
        variant = placements[corner]['variant']
        assert x * sx > 0 and y * sy > 0, corner
        pocket, letter, pads = coupon_shapes(variant, x, y)
        etch = pocket.union(letter)
        assert box(-4500, -4500, 4500, 4500).covers(etch), corner
        keepout = etch.buffer(50)
        window = etch.envelope.buffer(65, join_style=2)
        existing = region_in(window, metal, index)
        voids = get_parts(window.difference(existing))
        holes = [hole for hole in voids if hole.intersects(keepout)]
        for hole in holes:
            x0, y0, x1, y1 = hole.bounds
            assert hole.area <= 30 and max(x1 - x0, y1 - y0) <= 6, (corner, 'existing feature in coupon keep-out')
            additions.extend(polygons(hole.difference(etch)))
        cuts.extend(polygons(etch))
        for pad in pads:
            additions.extend(polygons(pad))

        # Fill only the old cross window. The die ID sits outside this window.
        cross = box(4450, 4450, 4850, 4850)
        cross = scale(cross, sx, sy, origin=(0, 0))
        assert not cross.intersects(keepout)
        additions.extend(polygons(cross.difference(region_in(cross, metal, index))))
        checks.append((corner, variant, x, y, pocket, letter, pads, cross, window))

    # Keep the original positive metal and cut only polygons that meet an etch.
    cut_shapes = [Polygon(cut.points) for cut in cuts]
    hits = index.query(cut_shapes)
    by_polygon = {}
    for cut_i, metal_i in zip(*hits):
        by_polygon.setdefault(metal_i, []).append(cuts[cut_i])
    result = gdstk.Library(unit=lib.unit, precision=lib.precision)
    cell = result.new_cell('TOP')
    for i, poly in enumerate(metal):
        if i in by_polygon:
            cell.add(*gdstk.boolean(poly, by_polygon[i], 'not', precision=GRID, layer=1))
        else:
            cell.add(poly)
    cell.add(*additions)
    result.write_gds(str(output), max_points=199)
    return checks


def check_corners(filename, checks):
    lib = gdstk.read_gds(str(filename))
    top, = lib.top_level()
    metal = top.get_polygons()
    index = STRtree([box(*low, *high) for low, high in (p.bounding_box() for p in metal)])
    views = {}
    for corner, variant, x, y, pocket, letter, pads, cross, window in checks:
        local = region_in(window, metal, index)
        etch = pocket.union(letter)
        pads_in_mask = local.intersection(pocket)
        assert pads_in_mask.symmetric_difference(pads[0].union(pads[1])).area < AREA_TOL, corner
        assert len(get_parts(pads_in_mask)) == 2 and not local.covers(Point(x, y)), corner
        assert local.intersection(letter).area < AREA_TOL, (corner, 'letter etch')
        assert cross.difference(region_in(cross, metal, index)).area < AREA_TOL, (corner, 'old cross')
        keepout = etch.buffer(50).difference(etch)
        assert keepout.difference(local).area < AREA_TOL, (corner, 'holes within keep-out')
        envelope = GEOMETRY['variants'][variant].get('array_envelope_um')
        if envelope:
            array = translate(box(*envelope), x, y)
            assert not local.intersects(array) and abs(array.distance(pocket.boundary) - 55) < GRID, corner
        views[corner] = (variant, local, window.bounds)
    return views


def plot_corners(views, filename):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as Patch

    fig, axes = plt.subplots(2, 2, figsize=(6.5, 8.2), layout='constrained')
    for corner, ax in zip(('UL', 'UR', 'LL', 'LR'), axes.flat):
        variant, metal, (x0, y0, x1, y1) = views[corner]
        for part in sorted(get_parts(metal), key=lambda part: part.area, reverse=True):
            ax.add_patch(Patch(part.exterior.coords, facecolor='#939da8', edgecolor='none'))
            for hole in part.interiors:
                ax.add_patch(Patch(hole.coords, facecolor='white', edgecolor='none'))
        ax.set(xlim=(x0, x1), ylim=(y0, y1), aspect='equal',
               title=f'{corner}  {variant.capitalize()}')
        ax.set_xticks([])
        ax.set_yticks([])
    fig.savefig(filename, dpi=180)
    plt.close(fig)


def export_design(name):
    home = ROOT / name
    candidate = json.loads((home / 'design.json').read_text(encoding='utf-8'))['wafer_candidate']
    source = home / candidate['parent_fabrication_gds']
    if sha256(source.read_bytes()).hexdigest() != candidate['parent_gds_sha256']:
        raise ValueError(f'{source.name} differs from the recorded parent GDS')
    targets = [home / candidate[key] for key in ('fabrication_gds', 'gds_preview', 'corner_preview')]
    assert source not in targets
    with TemporaryDirectory(dir=home / 'gds') as temp:
        mask, preview, corners = (Path(temp) / path.name for path in targets)
        checks = add_corner_coupons(source, mask, candidate['corner_test_junctions'])
        views = check_corners(mask, checks)
        verify_output(mask, name, True)
        plot_gds(str(mask), str(preview), window=(-5.15, 5.15, -5.15, 5.15))
        plot_corners(views, corners)
        for staged, target in zip((mask, preview, corners), targets):
            staged.replace(target)
    print(f'{name}: ' + ', '.join(f'{corner}={variant}' for corner, (variant, *_) in views.items()),
          flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--design', choices=DESIGNS, nargs='+', default=DESIGNS)
    for name in parser.parse_args().design:
        export_design(name)
