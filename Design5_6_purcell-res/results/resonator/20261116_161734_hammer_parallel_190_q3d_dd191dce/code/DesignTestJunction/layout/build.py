"""Export and verify the two standalone test-junction layouts."""

from json import loads
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import gdstk
from qiskit_metal import designs, MetalGUI
from shapely import get_parts, union_all
from shapely.affinity import rotate, scale
from shapely.geometry import Point, Polygon, box


DESIGN = Path(__file__).resolve().parents[1]
ROOT = DESIGN.parent
sys.path[:0] = [str(ROOT / 'pattern design'), str(ROOT / 'code' / 'general'),
                str(DESIGN / 'layout' / 'components')]

from coupon_mark import CouponMark
from gds_tools import export_gds
from test_junction import HookedTestJunction

GEOMETRY = loads((DESIGN / 'parameters' / 'geometry.json').read_text(encoding='utf-8'))
CHEESED_GROUND = (1, 100)


def add_test_junction(design, variant, name='TestJJ', *, fab_contacts=False, **placement):
    """Place an aligned or offset coupon into any planar design."""
    options = (GEOMETRY['common_options'] | GEOMETRY['variants'][variant]['component'] |
               {'orientation': '90'})
    return HookedTestJunction(design, name, options=options | placement |
                              {'fab_contacts': str(fab_contacts)})


def build_design(variant, *, fab_contacts=False, label_options=None):
    design = designs.DesignPlanar()
    design.chips.main.size.update(GEOMETRY['chip'])
    add_test_junction(design, variant, fab_contacts=fab_contacts)
    CouponMark(design, 'CouponMark', options=GEOMETRY['label'] | (label_options or {}) |
               {'letter': GEOMETRY['variants'][variant]['label']})
    return design


def read_region(polys):
    return union_all([Polygon(poly.points).buffer(0) for poly in polys], grid_size=.001)


def clean_perforations(filename):
    """Fill clipped holes and any hole extending into the 500 µm edge keep-out."""
    lib = gdstk.read_gds(str(filename))
    top, = lib.top_level()
    ground_cell, = [cell for cell in lib.cells if any(
        (poly.layer, poly.datatype) == CHEESED_GROUND for poly in cell.polygons)]
    ground_polys = [poly for poly in ground_cell.polygons
                    if (poly.layer, poly.datatype) == CHEESED_GROUND]
    (x0, y0), (x1, y1) = top.bounding_box()
    die_inner = box(x0 + 500, y0 + 500, x1 - 500, y1 - 500)
    ground = read_region(ground_polys)
    fills = []
    for void in ground.interiors:
        hole = Polygon(void)
        # The pocket and the letter are larger intentional openings.
        if hole.area < 26 and (hole.area < 25 - .001 or not die_inner.covers(hole)):
            fills.append(gdstk.Polygon(list(void.coords)[:-1],
                                       layer=CHEESED_GROUND[0], datatype=CHEESED_GROUND[1]))
    ground_cell.add(*fills)
    lib.write_gds(str(filename), max_points=0)


def check_masks(reference_file, fab_file, design, variant):
    """Read back the mask in µm and verify pattern transfer and connectivity."""
    p = design.components['TestJJ'].p
    table = design.qgeometry.tables['poly'].set_index('name')
    assert len(table) == 4 and table.at['pocket', 'subtract'] and table.at['letter', 'subtract']
    assert design.qgeometry.tables['junction'].empty
    assert all(g.is_valid and g.area > 0 for g in table.geometry)

    hx, hy = (500 * design.parse_value(GEOMETRY['chip'][dim]) for dim in ('size_x', 'size_y'))
    px, py = 500 * p.pocket_height, 500 * p.pocket_width
    pocket = box(-px, -py, px, py)
    letter = scale(table.at['letter', 'geometry'], 1000, 1000, origin=(0, 0))
    envelope_um = GEOMETRY['variants'][variant].get('array_envelope_um')
    pocket_metal, pad_pairs, grounds = [], [], []
    half_gap = 500 * p.landing_gap
    body_edge = 500 * p.pad_body_gap

    for filename in (reference_file, fab_file):
        lib = gdstk.read_gds(str(filename))
        assert (lib.unit, lib.precision) == (1e-6, 1e-9)
        top, = lib.top_level()
        polys = top.get_polygons()
        assert {poly.layer for poly in polys} == {1}
        assert top.bounding_box() == ((-hx, -hy), (hx, hy))
        metal = read_region(polys)
        assert metal.is_valid
        parts = sorted(get_parts(metal), key=lambda part: part.area)
        assert len(parts) == 3, 'Expected two floating pads and one connected ground'
        ground = parts[2]
        assert ground.distance(pocket) < .001
        assert all(ground.distance(pad) >= 49.999 for pad in parts[:2])
        assert min(hx - px, hy - py) >= 500
        assert parts[0].distance(parts[1]) >= 2 * half_gap - .001
        assert not metal.covers(Point(0, 0)), 'Junction site must stay open'
        assert letter.intersection(metal).area < .001, 'Etched letter missing or blocked'
        assert ground.covers(Point(letter.bounds[2] + 5, letter.centroid.y))
        if envelope_um:
            envelope = box(*envelope_um)
            assert envelope.distance(ground) >= 55 - .001
            assert not envelope.intersects(metal), 'Reserved array envelope crosses base metal'

        pads = sorted((rotate(part, -90, origin=(0, 0)) for part in parts[:2]),
                      key=lambda part: part.centroid.x)
        for pad, body_window in zip(pads, [box(-py, -px, -body_edge, px),
                                          box(body_edge, -px, py, px)]):
            x0, y0, x1, y1 = pad.intersection(body_window).bounds
            assert abs((x1 - x0) - 1000 * p.pad_length) < .001
            assert abs((y1 - y0) - 1000 * p.pad_width) < .001
        pocket_metal.append(metal.intersection(box(-px + .001, -py + .001, px - .001, py - .001)))
        pad_pairs.append(pads)
        grounds.append(ground)

    holes = grounds[0].difference(grounds[1])
    assert not holes.is_empty, 'Fabrication ground must be perforated'
    assert all(abs(hole.area - 25) < .001 for hole in get_parts(holes))
    assert holes.distance(pocket.union(letter)) >= 50
    assert box(-hx + 500, -hy + 500, hx - 500, hy - 500).covers(holes)
    assert not pocket_metal[1].difference(pocket_metal[0]).area
    half_hook, depth, shift = 500 * p.hook_width, 1000 * p.hook_depth, 1000 * p.hook_offset
    expected_cuts = (box(-half_gap - depth, -half_hook, -half_gap, half_hook),
                     box(half_gap, shift - half_hook, half_gap + depth, shift + half_hook))
    for ref, fab, cut in zip(*pad_pairs, expected_cuts):
        assert ref.difference(fab).symmetric_difference(cut).area < .001
    print(f'{reference_file.stem}: verified pad dimensions, two isolated pads, open junction, '
          f'top recess x offset {-shift:g} um, etched {GEOMETRY["variants"][variant]["label"]}, '
          'and intact fab holes with keep-outs.')


def screenshots(design, folder, stem):
    gui = MetalGUI(design)
    gui.main_window.resize(1280, 850)
    try:
        gui.rebuild()
        gui.zoom_on_components(['TestJJ', 'CouponMark'])
        gui.qApp.processEvents()
        gui.screenshot(name=str(folder / f'{stem}_layout'), display=False)
        gui.canvas.zoom_to_rectangle((-0.026, -0.023, 0.026, 0.023))
        gui.qApp.processEvents()
        gui.screenshot(name=str(folder / f'{stem}_hook'), display=False)
    finally:
        gui.main_window.close()


def main(only_variant=None):
    revision = GEOMETRY['geometry_revision']
    with TemporaryDirectory(dir=DESIGN) as tmp:
        staged = Path(tmp)
        for variant in (only_variant,) if only_variant else GEOMETRY['variants']:
            stem = f'{DESIGN.name}_{revision}_{variant}'
            reference = staged / f'{stem}_reference.gds'
            fab = staged / f'{stem}_fab.gds'
            export_gds(build_design(variant), str(reference))
            design = build_design(variant, fab_contacts=True)
            export_gds(design, str(fab), fab=True)
            clean_perforations(fab)
            check_masks(reference, fab, design, variant)
            screenshots(design, staged, stem)
        for artifact in staged.iterdir():
            target = 'gds' if artifact.suffix == '.gds' else 'fig'
            artifact.replace(DESIGN / target / artifact.name)


if __name__ == '__main__':
    main()
