"""Export wafer fabrication masks with the current side launch coordinates."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import subprocess
from tempfile import TemporaryDirectory

import gdstk
from shapely import STRtree, get_parts
from shapely.affinity import scale, translate
from shapely.geometry import Point, Polygon, box

from export_wafer_masks import verify_output
from gds_tools import export_gds, write_wafer_gds
from project_paths import ROOT

sys.path.insert(0, str(ROOT / 'DesignTestJunction' / 'layout'))
from build import build_design, GEOMETRY, read_region
from wafer_pcbport import wafer_pcb_ports
from xmark import CORNERS


DESIGNS = ('Design1', 'Design2', 'Design3', 'Design3_2', 'Design3_3', 'Design5', 'Design6')
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


def coupon_shapes(variant, x, y, label_um):
    label_options = {key: f'{value}um' for key, value in label_um.items()}
    design = build_design(variant, fab_contacts=True, label_options=label_options)
    table = design.qgeometry.tables['poly']
    placed = {row.name: translate(scale(row.geometry, 1000, 1000, origin=(0, 0)), x, y)
              for row in table.itertuples()}
    return placed['pocket'], placed['letter'], (placed['pad_l'], placed['pad_r'])


def add_corner_coupons(source, output, placements, label_um):
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
        pocket, letter, pads = coupon_shapes(variant, x, y, label_um)
        etch = pocket.union(letter)
        assert box(-4500, -4500, 4500, 4500).covers(pocket), corner
        assert box(-4700, -4700, 4700, 4700).covers(letter), corner
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

        # Replace the entire old mark and a 5 µm margin with one metal polygon.
        # Keeping the original cross polygon would leave its outline visible in CAD.
        cross = box(4450, 4450, 4850, 4850)
        cross = scale(cross, sx, sy, origin=(0, 0))
        marker_void = cross.difference(region_in(cross, metal, index))
        had_cross = marker_void.area > cross.area / 2
        if had_cross:
            repair = cross.buffer(5, join_style=2)
            assert not repair.intersects(keepout), corner
            margin = repair.difference(cross)
            assert margin.difference(region_in(repair, metal, index)).area < AREA_TOL, corner
            cuts.extend(polygons(repair))
            additions.extend(polygons(repair))
        checks.append((corner, variant, x, y, pocket, letter, pads, cross, had_cross, window))

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
    for corner, variant, x, y, pocket, letter, pads, cross, had_cross, window in checks:
        local = region_in(window, metal, index)
        etch = pocket.union(letter)
        pads_in_mask = local.intersection(pocket.buffer(-GRID))
        assert pads_in_mask.symmetric_difference(pads[0].union(pads[1])).area < AREA_TOL, corner
        assert len(get_parts(pads_in_mask)) == 2 and not local.covers(Point(x, y)), corner
        assert local.intersection(letter).area < AREA_TOL, (corner, 'letter etch')
        if had_cross:
            # An etched silhouette would fail the area check. A filled silhouette
            # made of touching polygons would fail the single-polygon check.
            assert cross.difference(region_in(cross, metal, index)).area < AREA_TOL, (corner, 'old cross')
            interiors = [Polygon(metal[i].points).buffer(0).intersection(cross)
                         for i in index.query(cross)]
            covering = [part for part in interiors if part.area > AREA_TOL]
            assert len(covering) == 1 and covering[0].covers(cross), (corner, 'old marker seam')
        keepout = etch.buffer(50).difference(etch)
        assert keepout.difference(local).area < AREA_TOL, (corner, 'holes within keep-out')
        envelope = GEOMETRY['variants'][variant].get('array_envelope_um')
        if envelope:
            array = translate(box(*envelope), x, y)
            assert not local.intersects(array) and abs(array.distance(pocket.boundary) - 55) < GRID, corner


def build_launch_mask(name, staging):
    """Rebuild the fab mask from the Design's layout with the ±2 mm anchors."""
    sys.path.insert(0, str(ROOT / name / 'layout'))
    raw = staging / 'launches_raw.gds'
    rounded = staging / 'launches_rounded.gds'
    mask = staging / 'launches_wafer.gds'
    if name == 'Design2':
        from layout3_full import export

        source = ROOT / name / 'gds' / 'previous' / 'layout3 ftf.gds'
        export(source, raw, fab=True, standard_launches=True, port_options=wafer_pcb_ports)
    else:
        options = dict(jj_metal=False, fab_contacts=True)
        if name == 'Design6':
            from design6_layout import build_design as build_chip

            options['wafer_ports'] = True
        else:
            options.update(include_tc_landing_pads=True, include_test_junctions=False)
            if name == 'Design3_3':
                from layout import build_design as build_chip

                options['rounded'] = False
            else:
                from ftf_layout import build_design as build_chip
        export_gds(build_chip(**options), str(raw), fab=True)
    if name in ('Design3_2', 'Design3_3', 'Design6'):
        sys.path.insert(0, str(ROOT / 'Design3_2' / 'layout'))
        from round_mask import round_mask

        round_mask(raw, rounded)
        raw = rounded
    write_wafer_gds(raw, mask)
    return mask


def export_design(name):
    home = ROOT / name
    candidate = json.loads((home / 'design.json').read_text(encoding='utf-8'))['wafer_candidate']
    if name == 'Design5':
        from shutil import copyfile

        source = home / candidate['fabrication_gds']
        target = home / candidate['delivery_fabrication_gds']
        verify_output(source, name, True, expect_die_id=bool(candidate.get('die_id')))
        copyfile(source, target)
        assert sha256(source.read_bytes()).digest() == sha256(target.read_bytes()).digest()
        return
    target = home / candidate['fabrication_gds']
    target.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=home / 'gds') as temp:
        mask = Path(temp) / target.name
        launch_mask = build_launch_mask(name, Path(temp))
        placements = candidate.get('corner_test_junctions')
        if placements:
            checks = add_corner_coupons(launch_mask, mask, placements, candidate['corner_label_um'])
            check_corners(mask, checks)
        else:
            launch_mask.replace(mask)
        verify_output(mask, name, True, ports_override=wafer_pcb_ports)
        mask.replace(target)
    if placements:
        print(f'{name}: ' + ', '.join(f'{corner}={item["variant"]}' for corner, item in placements.items()),
              flush=True)


if __name__ == '__main__':
    import geopandas as gpd

    if not hasattr(gpd.GeoDataFrame, 'append'):
        gpd.GeoDataFrame.append = gpd.GeoDataFrame._append
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--design', choices=DESIGNS, nargs='+', default=DESIGNS)
    selected = parser.parse_args().design
    if len(selected) == 1:
        export_design(selected[0])
    else:
        # Several designs use modules with the same name in their own layout directories.
        for name in selected:
            subprocess.run([sys.executable, '-X', 'utf8', '-B', __file__, '--design', name], check=True)
