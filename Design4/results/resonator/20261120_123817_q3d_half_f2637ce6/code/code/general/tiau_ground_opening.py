"""Open a rectangular window in the ground metal at every chip-level TiAu marker site of a wafer assembly.

The chip markers (Mask2 layer 2/0) are fixed by the existing TiAu photomask. So that each marker sits on
bare substrate, the assembly metal (1/0) is removed inside a rectangle around every marker site:
the marker bounding box (84 x 68 um) grown by ``--margin`` on each side.

The assembly must place every die on the Mask2 grid without rotation (``assemble_wafer_mask2_grid.py``).
Then every die of one design sees its markers at the same local positions, so the rectangles are cut once
from each flat source cell (D2, D3, ...). The script checks this, and checks that no die-ID polygon and no
die outside the grid is touched.

Usage (repository root, qmetal):
    python -X utf8 -B code/general/tiau_ground_opening.py TiAu_batch/summery_mask2grid.gds \
        --mask TiAu_batch/inputs/Fluxonium_2Q_bandpass_filter_Mask2.gds --out-stem summery_mask2grid_open
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import gdstk
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from tiau_marker_overlay import CHIP_LAYER, METAL, group_sites  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('assembly', type=Path, help='assembly GDS with <stem>.json next to it')
    parser.add_argument('--mask', required=True, type=Path, help='TiAu marker photomask (Mask2)')
    parser.add_argument('--margin', type=float, default=50.0, help='clearance around the marker bbox (um)')
    parser.add_argument('--out-stem', required=True)
    args = parser.parse_args()

    gds_in = ROOT / args.assembly
    manifest = json.loads(gds_in.with_suffix('.json').read_text(encoding='utf-8'))
    out_gds = gds_in.with_name(f'{args.out_stem}.gds')

    mask_top, = gdstk.read_gds(str(ROOT / args.mask)).top_level()
    sites = group_sites(mask_top.get_polygons(layer=CHIP_LAYER[0], datatype=CHIP_LAYER[1]), gap=40.0)

    lib = gdstk.read_gds(str(gds_in))
    top, = lib.top_level()
    cells = {c.name: c for c in lib.cells}
    dies = {d['cell']: d for d in manifest['dies']}

    # Local rectangle per design, from the marker sites of every die of that design.
    half = np.asarray(manifest['die_slot_size_um']) / 2
    local = {}
    for ref in top.references:
        die = dies[ref.cell.name]
        assert ref.rotation == 0 and not ref.x_reflection and ref.magnification == 1, ref.cell.name
        origin = np.asarray(ref.origin)
        boxes = sorted(tuple(np.round(np.concatenate(s['bbox']) - np.tile(origin, 2), 3)) for s in sites
                       if np.all(np.abs(s['centre'] - origin) <= half))
        if len(boxes) != 4:
            raise SystemExit(f"{die['die_id']}: {len(boxes)} marker sites, expected 4")
        if local.setdefault(die['design'], boxes) != boxes:
            raise SystemExit(f"{die['die_id']}: marker positions differ from other {die['design']} dies")

    m = args.margin
    openings, removed = {}, {}
    for design, boxes in local.items():
        rects = [gdstk.rectangle((x0 - m, y0 - m), (x1 + m, y1 + m)) for x0, y0, x1, y1 in boxes]
        cell = cells[design]
        assert not cell.references, design
        # Only metal polygons whose bounding box reaches an opening are cut; the rest stay untouched.
        def hits(p):
            (a, b), (c, d) = p.bounding_box()
            return any(a < x1 + m and c > x0 - m and b < y1 + m and d > y0 - m for x0, y0, x1, y1 in boxes)
        touched = [p for p in cell.polygons if (p.layer, p.datatype) == METAL and hits(p)]
        before = sum(p.area() for p in touched)
        cut = gdstk.boolean(touched, rects, 'not', precision=1e-3, layer=METAL[0], datatype=METAL[1])
        removed[design] = round(before - sum(p.area() for p in cut), 1)
        cell.remove(*touched)
        cell.add(*cut)
        openings[design] = [[round(float(v), 3) for v in (*r.bounding_box()[0], *r.bounding_box()[1])]
                            for r in rects]
        # The die-ID polygons live in the WAFER_Dxx cells; they must stay clear of the openings.
        for name, die in dies.items():
            if die['design'] == design and cells[name].polygons:
                if gdstk.boolean(cells[name].polygons, rects, 'and', precision=1e-3):
                    raise SystemExit(f'{name}: die ID overlaps a ground opening')

    lib.write_gds(str(out_gds), max_points=8190)

    # Check on the written file: no metal is left in any opening of any die (source cell + die-ID cell).
    out_cells = {c.name: c for c in gdstk.read_gds(str(out_gds)).cells}
    worst = 0.0
    for name, die in dies.items():
        polys = out_cells[name].get_polygons(layer=METAL[0], datatype=METAL[1])
        rects = [gdstk.rectangle((x0, y0), (x1, y1)) for x0, y0, x1, y1 in openings[die['design']]]
        near = [p for p in polys if any(p.bounding_box()[0][0] < r.bounding_box()[1][0]
                                        and p.bounding_box()[1][0] > r.bounding_box()[0][0]
                                        and p.bounding_box()[0][1] < r.bounding_box()[1][1]
                                        and p.bounding_box()[1][1] > r.bounding_box()[0][1] for r in rects)]
        if near:
            worst = max(worst, sum(p.area() for p in gdstk.boolean(near, rects, 'and', precision=1e-3)))
    if worst > 1e-3:
        raise SystemExit(f'metal left in an opening: {worst} um^2')

    manifest = dict(manifest)
    manifest['description'] = manifest['description'] + ' Ground opened at the chip TiAu marker sites.'
    manifest['tiau_ground_opening'] = {
        'parent_gds': gds_in.name, 'parent_sha256': sha256(gds_in),
        'mask': args.mask.as_posix(), 'mask_sha256': sha256(ROOT / args.mask),
        'marker_layer': CHIP_LAYER, 'metal_layer': METAL,
        'rule': f'metal 1/0 removed inside the marker bounding box grown by {m:g} um per side',
        'margin_um': m, 'opening_size_um': [round(b[2] - b[0], 3) for b in [openings[next(iter(openings))][0]]]
        + [round(b[3] - b[1], 3) for b in [openings[next(iter(openings))][0]]],
        'openings_local_um_by_design': openings,
        'metal_removed_per_die_um2': removed,
        'sites_opened': 4 * len(top.references),
        'check_metal_left_in_openings_um2': round(worst, 6),
        'output_gds': out_gds.name, 'output_sha256': sha256(out_gds),
    }
    out_gds.with_suffix('.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(manifest['tiau_ground_opening'], indent=2))


if __name__ == '__main__':
    main()
