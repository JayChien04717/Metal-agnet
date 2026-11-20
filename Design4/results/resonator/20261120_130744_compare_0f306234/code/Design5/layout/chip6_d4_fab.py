"""Fabrication GDS of the Design5 six-channel Purcell chip (chip6_d4.py), 10 x 10 mm.

Lab fab rules (skill/gds_export.md, skill/chip_drc.md C10/C13), same as Design1/3:
  - 5 x 5 um square cheese holes on a 20 um pitch in the ground (layer 1),
  - 50 um to etched features and a 547.5 um unperforated wafer margin,
  - optional test-junction pockets in the four corners (800 x 300 um, pads 225/295 x 129 um, gap 211 um),
  - no die markings by default, with optional legacy marks.
GDS written in um (chip -5000 .. 5000).

Usage: python chip6_d4_fab.py --test-junctions no [--die-id PFR00]
  -> gds/Design5_<date>_fab.gds, fig/fab/Design5_<date>_{chip,testjj,holes}.png;
  fab settings: section `fab` of parameters/Design5_<date>.json
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import json, os, sys, warnings
from pathlib import Path
from tempfile import TemporaryDirectory
warnings.filterwarnings('ignore')
HERE = _S('Design5/layout'); ROOT = _S('Design5')
sys.path.insert(0, HERE)
import chip6_d4 as C
import design5_layout as d5
from gds_tools import export_gds, write_wafer_gds
from xmark import XMark, corner_mark_centres

f3 = d5.f3                                   # Design3 ftf_layout: TEST_JJ, marks, fab options
FAB_GDS = _relocate(os.path.join(ROOT, 'gds', C.TAG + '_fab.gds'))
FIG_DIR = _S('Design5/fig/fab')


def wafer_hammer():
    config = json.loads((Path(ROOT) / 'design.json').read_text(encoding='utf-8'))
    return config['wafer_candidate'].get('hammer_geometry_um')


def build_fab(die_id='', *, include_test_junctions=False, include_crosses=False,
              hammer=None, lengths=None):
    if hammer is None and not include_crosses:
        hammer = wafer_hammer()
    d = C.build_chip(lengths, hammer=hammer)
    if include_test_junctions:
        raise ValueError('include_test_junctions: the old corner-pocket TestJunction was retired on 2026-11-20; corner F/T coupons (HookedTestJunction, pattern design/junction_test_coupon.py) are added by code/general/export_wafer_test_junctions.py')
    if include_crosses:
        f3.add_corner_marks(d, die_id=die_id, inset=f3.MARK_INSET,
                            label_corner='LR', label_left='True', **f3.MARK)
    elif die_id:
        x, y = corner_mark_centres(10, 10, f3.MARK_INSET)['LR']
        XMark(d, 'FabDieId', options=dict(pos_x=f'{x}mm', pos_y=f'{y}mm', label=die_id,
                                        label_left='True', include_cross='False', **f3.MARK))
    return d


def export(d, filename):
    export_gds(d, filename, fab=True)


def export_wafer(output, *, die_id='', hammer_summary=None):
    """Write a validated positive-metal wafer GDS without corner structures."""
    import gdstk
    from export_wafer_masks import SIDE_PORTS, verify_output
    from wafer_cheese_grid import replace_wafer_cheese
    from shapely import STRtree, get_parts, union_all
    from shapely.geometry import Polygon, box

    output = Path(output).resolve()
    lengths_file = Path(ROOT) / 'parameters' / 'Design5_2026-11-13.json'
    lengths = json.loads(lengths_file.read_text(encoding='utf-8'))['channels']
    hammer = wafer_hammer()
    if hammer_summary:
        comparison = json.loads(Path(hammer_summary).read_text(encoding='utf-8'))
        if not comparison['meets_2_percent']:
            raise ValueError('The selected hammer missed the 2% capacitance tolerance')
        hammer = comparison['selected']['geometry']

    design = build_fab(die_id, hammer=hammer, lengths=lengths)
    assert design.qgeometry.tables['junction'].empty
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as temp:
        unnormalized = Path(temp) / 'raw_fab.gds'
        mask = Path(temp) / output.name
        export_gds(design, str(unnormalized), fab=True, cheese=False)
        write_wafer_gds(unnormalized, mask)
        replace_wafer_cheese(mask, launches=SIDE_PORTS)
        verify_output(mask, 'Design5', True, expect_die_id=bool(die_id))

        library = gdstk.read_gds(str(mask))
        top, = library.top_level()
        polygons = top.get_polygons()
        index = STRtree([box(*low, *high) for low, high in (p.bounding_box() for p in polygons)])

        def metal_in(window):
            clip = gdstk.rectangle(window.bounds[:2], window.bounds[2:])
            pieces = gdstk.boolean([polygons[i] for i in index.query(window)], clip, 'and', precision=.001)
            return union_all([Polygon(p.points).buffer(0) for p in pieces], grid_size=.001)

        for x, y in corner_mark_centres(10, 10, f3.MARK_INSET).values():
            former_window = box(x * 1000 - 200, y * 1000 - 200, x * 1000 + 200, y * 1000 + 200)
            voids = former_window.difference(metal_in(former_window))
            assert all(hole.area <= 30 for hole in get_parts(voids)), 'Former cross window is still etched'
        for sx in (-1, 1):
            for sy in (-1, 1):
                x = sx * f3.TEST_JJ['x'] * 1000
                y = sy * f3.TEST_JJ['y'] * 1000
                half_w, half_h = (500 * size for size in f3.TEST_JJ['pocket'])
                pocket = box(x - half_w, y - half_h, x + half_w, y + half_h)
                voids = pocket.difference(metal_in(pocket))
                assert all(hole.area <= 30 for hole in get_parts(voids)), 'Test-junction pocket appears in fab'

        if hammer:
            for ch, x, up in C.ORDER:
                spacing = C.hammer_spacing_mm(hammer) * 1000
                cx = x * 1000 + spacing / 2
                cy = C.P.TAP_Y * (1000 if up else -1000)
                half_gap = hammer['gap_um'] / 2
                half_height = hammer['height_um'] / 2
                gap = box(cx - half_gap, cy - half_height, cx + half_gap, cy + half_height)
                assert metal_in(gap).area < .001, (ch, 'Hammer gap is blocked')
                for sign in (-1, 1):
                    px = cx + sign * (half_gap + hammer['width_um'] / 2)
                    pad_probe = box(px - 1, cy - half_height + 1, px + 1, cy + half_height - 1)
                    assert pad_probe.difference(metal_in(pad_probe)).area < .001, (ch, 'Hammer head size mismatch')
                    trace_x = cx + sign * spacing / 2
                    for lo, hi in ((-350, -65), (65, 350)):
                        window = box(trace_x - 16, cy + lo, trace_x + 16, cy + hi)
                        expected = box(trace_x - 8, cy + lo, trace_x + 8, cy + hi)
                        assert metal_in(window).symmetric_difference(expected).area < .001, (
                            ch, 'CPW does not follow the straight axis')
        mask.replace(output)
    print(f'Fab GDS: {output}', flush=True)


def plots(filename, *, include_test_junctions=False):
    import gdspy, matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    top = [c for c in gdspy.GdsLibrary(infile=filename).top_level() if c.name == 'TOP'][0]
    spec = top.get_polygons(by_spec=True)
    os.makedirs(FIG_DIR, exist_ok=True)
    for tag, win, size in (('chip', (-5100, 5100, -5100, 5100), (9, 9)),
                           ('testjj', (-4450, -3450, 3950, 4550), (8, 5)),
                           ('holes', (600, 1100, 300, 800), (8, 8))):
        if tag == 'testjj' and not include_test_junctions:
            continue
        fig, ax = plt.subplots(figsize=size)
        for (layer, dt), ps in sorted(spec.items()):
            if layer == 1:
                ax.add_collection(PolyCollection(ps, facecolor='#9aa3ad' if dt == 100 else '#4b5560', edgecolor='none'))
        ax.set_xlim(win[:2]); ax.set_ylim(win[2:]); ax.set_aspect('equal')
        ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)'); ax.set_title(f'{C.TAG} ({tag})')
        fig.tight_layout(); fig.savefig(_relocate(os.path.join(FIG_DIR, f'{C.TAG}_{tag}.png')), dpi=100); plt.close(fig)
    bb = top.get_bounding_box()
    n = {k: len(v) for k, v in spec.items()}
    print('bbox (um):', bb.tolist(), '\npolygons per layer/datatype:', n)


if __name__ == '__main__':
    import geopandas as gpd
    if not hasattr(gpd.GeoDataFrame, 'append'):
        gpd.GeoDataFrame.append = gpd.GeoDataFrame._append
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--die-id', default='')
    parser.add_argument('--output', type=Path, help='Write a wafer fab mask only')
    parser.add_argument('--hammer-summary', type=Path, help='Use the selected local Q3D hammer')
    parser.add_argument('--test-junctions', choices=('yes', 'no'), default='no',
                        help='retired: only "no" (Design5 has no corner test junctions)')
    args = parser.parse_args()
    include_tests = args.test_junctions == 'yes'
    if args.output:
        if include_tests:
            parser.error('--output wafer export uses no corner test junctions')
        export_wafer(args.output, die_id=args.die_id, hammer_summary=args.hammer_summary)
        sys.exit(0)
    if args.hammer_summary:
        parser.error('--hammer-summary requires --output')
    d = build_fab(args.die_id, include_test_junctions=include_tests)
    export(d, FAB_GDS)
    plots(FAB_GDS, include_test_junctions=include_tests)
    print('gds', FAB_GDS)
