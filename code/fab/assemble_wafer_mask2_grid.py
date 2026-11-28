"""Re-place the wafer dies on the old single-chip grid so that the fixed TiAu marker mask (Mask2) lines up.

The old wafer (``Fluxonium_2Q_bandpass_filter_*`` laser-cut file and Mask2) has 61 chips of 10 x 10 mm on a
10.2 mm pitch: an odd 9 x 9 grid with a chip centred at (0, 0), rows of 3, 7, 7, 9, 9, 9, 7, 7, 3 chips and
200 um streets. This script uses exactly those 61 centres for the same six die sources as ``assemble_wafer``.

Differences from ``assemble_wafer`` (20261119 summery):
- pitch 10200 um in x and y (was 10380 x 10410) and a chip, not a street crossing, at the wafer centre;
- 61 dies instead of 52. Counts are set in COUNTS. No two copies of one design touch side by side.
- Chip4 and PCM (Allspice, up to 10280 x 10310 um) are wider than the 10 mm chips. They never touch each other,
  even diagonally,
  or on the four outermost corner sites. The script checks that no two actual die outlines overlap.
- Dicing crosses: device-only wafers (no Chip4/PCM) get the same 250 x 250 um, 5 um-stroke metal plus as the NCKU /
  standard wafer (assemble_wafer.alignment_cross) at every die-corner street intersection; the streets stay metal-free.

Usage (repository root, qmetal):
    python -X utf8 -B code/general/assemble_wafer_mask2_grid.py 20261119 --out TiAu_batch \
        --input D2=... --input D3=... --input D5=... --input D6=... [--source-folder 20261116]
"""

import argparse
import json
import random
import sys
from collections import Counter
from math import hypot
from os.path import relpath
from pathlib import Path
from tempfile import TemporaryDirectory

import gdstk
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from assemble_wafer import (COLORS, CROSS_SPAN, CROSS_STROKE, GDS_MAX_POINTS, INPUTS, METAL, WAFER_GUIDE,  # noqa: E402
                            WAFER_RADIUS, add_die_ids, alignment_cross, check_die_ids, load_die, polygon_signature, wafer_outline)

PITCH = 10200.0
CHIP = 10000.0
CROSS_CELL = 'ALIGNMENT_CROSS'   # 250 x 250 um, 5 um-stroke metal plus (same as the NCKU / standard wafer)
GLOBAL_CROSS_CELL = 'GLOBAL_MARK_CROSS'
ROW_HALF_WIDTH = {3: 3, 2: 3, 1: 3, 0: 3, -1: 3, -2: 3, -3: 3}   # inner 7 x 7 (outer clamp strips omitted)
CORNER_SLOTS = {(i, j) for i in (-3, 3) for j in (-3, 3)}   # 7 x 7 corners, left empty -> 45 dies
ROW_HALF_WIDTH_61 = {4: 1, 3: 3, 2: 3, 1: 4, 0: 4, -1: 4, -2: 3, -3: 3, -4: 1}   # full 61-slot TiAu mask grid
COUNTS = {'D2': 13, 'D3': 13, 'D5': 13, 'D6': 12, 'Chip4': 5, 'PCM': 5}
COUNTS_DEVICES_45 = {'D2': 12, 'D3': 11, 'D5': 11, 'D6': 11}
ALLSPICE = ('Chip4', 'PCM')
SEED = 20261119


def die_corners(dies):
    """Street intersections at the corners of the populated dies (dicing-cross sites), sorted."""
    return sorted({(round(d['centre_um'][0] + sx * PITCH / 2, 3), round(d['centre_um'][1] + sy * PITCH / 2, 3))
                   for d in dies for sx in (-1, 1) for sy in (-1, 1)})


def sites(row_half_width=None):
    rhw = ROW_HALF_WIDTH if row_half_width is None else row_half_width
    skip = CORNER_SLOTS if rhw is ROW_HALF_WIDTH else set()   # inner grid: four corner slots stay empty
    return [(i, j) for j in sorted(rhw, reverse=True)
            for i in range(-rhw[j], rhw[j] + 1) if (i, j) not in skip]


def region(i, j):
    """Centre 3 x 3 block, then a pinwheel of four quadrants."""
    if abs(i) <= 1 and abs(j) <= 1:
        return 'A'
    if i <= 0 < j:
        return 'B'
    if i > 0 <= j:
        return 'C'
    if i < 0 and j <= 0:
        return 'D'
    return 'E'


def assign(target_counts=None, row_half_width=None):
    """Backtracking assignment with exact global counts, balanced regions and the adjacency rules."""
    target_counts = dict(COUNTS if target_counts is None else target_counts)
    grid = sites(row_half_width)
    if sum(target_counts.values()) != len(grid):
        raise ValueError(f'Target counts sum to {sum(target_counts.values())}, but grid has {len(grid)} sites')
    regions = {s: region(*s) for s in grid}
    size = Counter(regions.values())
    allspice = [n for n in ALLSPICE if target_counts.get(n, 0) > 0]
    devices = [n for n in target_counts if n not in ALLSPICE and target_counts[n] > 0]
    caps = {r: {**{n: 1 for n in allspice},
                **{n: -(-(size[r] - len(allspice)) // len(devices)) for n in devices}} for r in size}
    corners = {(i, j) for i in (-3, 3) for j in (-3, 3)}
    rng = random.Random(SEED)
    placed, total, per_region = {}, Counter(), {r: Counter() for r in size}

    # Step 1: optional Chip4 and PCM per region. Two Allspice dies overlap even diagonally (their
    # corners reach past 5100 um), so no two may be among each other's eight neighbours.
    if allspice:
        def near(a, b):
            return max(abs(a[0] - b[0]), abs(a[1] - b[1])) <= 1
        by_region = {r: [s for s in grid if regions[s] == r and s not in corners] for r in sorted(size)}
        for _ in range(200_000):
            chosen = {r: rng.sample(by_region[r], len(allspice)) for r in by_region}
            flat = [s for pair in chosen.values() for s in pair]
            if all(not near(a, b) for k, a in enumerate(flat) for b in flat[k + 1:]):
                break
        else:
            raise RuntimeError('No Chip4/PCM placement found')
        for r, picked in chosen.items():
            for s, name in zip(picked, allspice):
                placed[s] = name
                total[name] += 1
                per_region[r][name] += 1

    def allowed(s, name):
        r = regions[s]
        if total[name] >= target_counts[name] or per_region[r][name] >= caps[r][name]:
            return False
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            if placed.get((s[0] + di, s[1] + dj)) == name:
                return False
        return True

    # Step 2: the device designs, row by row, with backtracking.
    order = [s for s in grid if s not in placed]

    def solve(k):
        if k == len(order):
            return True
        s = order[k]
        names = list(devices)
        rng.shuffle(names)
        names.sort(key=lambda n: total[n] - target_counts[n])
        for name in names:
            if allowed(s, name):
                placed[s] = name
                total[name] += 1
                per_region[regions[s]][name] += 1
                if solve(k + 1):
                    return True
                del placed[s]
                total[name] -= 1
                per_region[regions[s]][name] -= 1
        return False

    if not solve(0):
        raise RuntimeError('No assignment satisfies the counts and adjacency rules')
    dies = [{'design': placed[s], 'region': regions[s], 'grid_index': list(s),
             'centre_um': [s[0] * PITCH, s[1] * PITCH]} for s in grid]
    dies.sort(key=lambda d: (-d['centre_um'][1], d['centre_um'][0]))
    for index, die in enumerate(dies):
        die['die_id'] = f'D{index:02d}'
        die['cell'] = f'WAFER_{die["die_id"]}'
    counts = {r: dict(per_region[r]) for r in sorted(size)}
    counts['Total'] = dict(total)
    return dies, counts


def outline(die, sources, cells):
    low, high = np.asarray(cells[die['design']].bounding_box())
    return low + die['centre_um'], high + die['centre_um']


def make_cross(span_um, stroke_um):
    """Build a centred metal plus polygon with total span ``span_um`` and bar width ``stroke_um``."""
    arm, half_w = span_um / 2, stroke_um / 2
    horizontal = gdstk.rectangle((-arm, -half_w), (arm, half_w))
    vertical = gdstk.rectangle((-half_w, -arm), (half_w, arm))
    cross, = gdstk.boolean(horizontal, vertical, 'or', precision=.001, layer=METAL[0], datatype=METAL[1])
    return cross


def verify(path, cells, dies, sources, crosses=(), global_marks=()):
    result = gdstk.read_gds(str(path))
    top, = result.top_level()
    assert (result.unit, result.precision) == (1e-6, 1e-9)
    assert result.layers_and_datatypes() == {METAL, WAFER_GUIDE}
    expected_refs = Counter((d['cell'], tuple(d['centre_um'])) for d in dies)
    expected_refs.update((CROSS_CELL, tuple(c)) for c in crosses)
    expected_refs.update((GLOBAL_CROSS_CELL, tuple(c)) for c in global_marks)
    assert Counter((r.cell_name, tuple(r.origin)) for r in top.references) == expected_refs
    assert all(r.rotation == 0 and r.magnification == 1 and not r.x_reflection for r in top.references)
    by_name = {c.name: c for c in result.cells}
    for name in sources:
        assert polygon_signature(by_name[name]) == polygon_signature(cells[name]), name
    boxes = [outline(d, sources, cells) for d in dies]
    for a in range(len(dies)):
        lo, hi = boxes[a]
        assert hypot(max(abs(lo[0]), abs(hi[0])), max(abs(lo[1]), abs(hi[1]))) < WAFER_RADIUS, dies[a]['die_id']
        for b in range(a + 1, len(dies)):
            lo2, hi2 = boxes[b]
            overlap = (lo[0] < hi2[0] and lo2[0] < hi[0] and lo[1] < hi2[1] and lo2[1] < hi[1])
            assert not overlap, (dies[a]['die_id'], dies[b]['die_id'])
    gaps = []
    for a in range(len(dies)):
        for b in range(a + 1, len(dies)):
            (lo, hi), (lo2, hi2) = boxes[a], boxes[b]
            if abs(dies[a]['centre_um'][1] - dies[b]['centre_um'][1]) < 1 and \
                    abs(dies[a]['centre_um'][0] - dies[b]['centre_um'][0]) - PITCH < 1:
                gaps.append(max(lo2[0] - hi[0], lo[0] - hi2[0]))
            if abs(dies[a]['centre_um'][0] - dies[b]['centre_um'][0]) < 1 and \
                    abs(dies[a]['centre_um'][1] - dies[b]['centre_um'][1]) - PITCH < 1:
                gaps.append(max(lo2[1] - hi[1], lo[1] - hi2[1]))
    radius = max(hypot(max(abs(lo[0]), abs(hi[0])), max(abs(lo[1]), abs(hi[1]))) for lo, hi in boxes)
    return {'minimum_neighbour_gap_um': round(min(gaps), 3), 'maximum_die_corner_radius_um': round(radius, 3),
            'minimum_radial_margin_um': round(WAFER_RADIUS - radius, 3)}


def plot(path, dies, counts, cells, active_names=None, crosses=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Rectangle

    fig, (ax, legend) = plt.subplots(1, 2, figsize=(12.4, 8), gridspec_kw={'width_ratios': [2.35, 1]},
                                    layout='constrained')
    ax.add_patch(Circle((0, 0), WAFER_RADIUS / 1000, fill=False, ec='#26394b', lw=1.4))
    for die in dies:
        (lo, hi) = (v / 1000 for v in outline(die, None, cells))
        ax.add_patch(Rectangle(lo, *(hi - lo), fc=COLORS[die['design']], ec='#243445', lw=.3))
        x, y = (v / 1000 for v in die['centre_um'])
        ax.text(x, y, f"{die['design']}\n{die['die_id']}", ha='center', va='center', fontsize=7.5, color='#152436')
    for (i, j) in [(-1.5, -1.5), (-1.5, 1.5), (1.5, -1.5), (1.5, 1.5)]:
        ax.add_patch(Rectangle((-15.3, -15.3), 30.6, 30.6, fill=False, ls='--', ec='#526274', lw=.8))
    ax.set(xlim=(-53, 53), ylim=(-53, 53), aspect='equal', xlabel='x (mm)', ylabel='y (mm)')
    ax.set_title(f'4-inch wafer on the NTHU TiAu grid  |  {len(dies)} dies  |  10 mm @ 10.2 mm', loc='left',
                 weight='bold')
    legend.set_axis_off()
    legend.text(0, .97, 'Dies per region', weight='bold', fontsize=13)
    names = list(active_names or COUNTS)
    rows = [*sorted(r for r in counts if r != 'Total'), 'Total']
    table = legend.table(cellText=[[counts[r].get(n, 0) for n in names] for r in rows], colLabels=names,
                         rowLabels=rows, colColours=[COLORS[n] for n in names], cellLoc='center',
                         bbox=[.07, .62, .90, .285])
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    quad_size = (len(dies) - 9) // 4
    legend.text(0, .52, f'A: centre 3 x 3 chips (dashed)\nB-E: pinwheel quadrants of {quad_size} chips\n\n'
                        f'{len(dies)} chip centres on the NTHU TiAu mask\n'
                        'grid (10.2 mm pitch, die at (0, 0), outer\n'
                        'clamp strips and 7 x 7 corners omitted, 45 dies).\n\n'
                        'Metal 1/0  ·  Wafer guide 99/0\n'
                        + ('Dicing crosses: 250 µm metal plus (5 µm)\nat every die corner' if crosses
                           else 'No street alignment crosses'), fontsize=9, va='top')
    fig.savefig(path, dpi=185, facecolor='white')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('folder', type=Path, help='wafer folder (default source folder)')
    parser.add_argument('--out', type=Path, required=True, help='output folder')
    parser.add_argument('--stem', default='summery_mask2grid')
    parser.add_argument('--source-folder', type=Path)
    parser.add_argument('--input', action='append', default=[], metavar='KEY=PATH')
    parser.add_argument('--count', action='append', default=[], metavar='KEY=N',
                        help='override target die count per design (e.g. D2=13)')
    parser.add_argument('--devices-only', action='store_true',
                        help='populate inner 7 x 7 minus its 4 corners (45 dies) with D2, D3, D5, D6 only (omit Chip4 and PCM)')
    parser.add_argument('--grid-61', action='store_true',
                        help='use all 61 slots (outer clamp slots and 7 x 7 corners included) instead of the 45-die inner grid')
    parser.add_argument('--global-cross', action='append', default=[], metavar='X,Y',
                        help='wafer-level global mark cross centre in um (repeatable, e.g. 0,40800)')
    parser.add_argument('--global-cross-span', type=float, default=1200.0,
                        help='total span of each wafer global mark cross in um (default 1200)')
    parser.add_argument('--global-cross-stroke', type=float, default=4.0,
                        help='stroke width of each wafer global mark cross in um (default 4)')
    parser.add_argument('--global-cross-window-margin', type=float, default=100.0,
                        help='bare-window margin around each wafer global mark cross in etch in um (default 100)')
    parser.add_argument('--replace', action='store_true')
    args = parser.parse_args()

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    outputs = [out / f'{args.stem}.{s}' for s in ('gds', 'png', 'json')] + [out / f'{args.stem}_die_id_check.png']
    if not args.replace and any(p.exists() for p in outputs):
        raise SystemExit('Output exists; pass --replace')
    rhw = ROW_HALF_WIDTH_61 if args.grid_61 else ROW_HALF_WIDTH
    if args.count:
        target_counts = {k: int(v) for item in args.count for k, v in [item.split('=', 1)]}
    elif args.devices_only or not args.grid_61:
        target_counts = dict(COUNTS_DEVICES_45 if args.devices_only else
                             {'D2': 9, 'D3': 9, 'D5': 9, 'D6': 8, 'Chip4': 5, 'PCM': 5})
    else:
        target_counts = dict(COUNTS)
    active_names = [k for k, v in target_counts.items() if v > 0]
    inputs = dict(item.split('=', 1) for item in args.input)
    cells, sources = {}, {}
    for name in active_names:
        filename = INPUTS[name]
        path = Path(inputs[name]) if name in inputs else (args.source_folder or args.folder) / filename
        cells[name], sources[name] = load_die(path, name)
        sources[name]['path_from_assembly'] = Path(relpath(path, out)).as_posix()

    dies, counts = assign(target_counts=target_counts, row_half_width=rhw)
    wafer_die_ids = add_die_ids(cells, sources, dies)
    marks = [] if any(n in active_names for n in ALLSPICE) else die_corners(dies)
    global_centres = [tuple(float(v) for v in item.split(',', 1)) for item in args.global_cross]
    lib = gdstk.Library('WAFER', unit=1e-6, precision=1e-9)
    top = lib.new_cell('WAFER_4_INCH')
    lib.add(*cells.values())
    top.add(*(gdstk.Reference(cells[d['cell']], origin=d['centre_um']) for d in dies))
    if marks:
        cross = lib.new_cell(CROSS_CELL)
        cross.add(alignment_cross())
        top.add(*(gdstk.Reference(cross, origin=c) for c in marks))
    if global_centres:
        gcross = lib.new_cell(GLOBAL_CROSS_CELL)
        gcross.add(make_cross(args.global_cross_span, args.global_cross_stroke))
        top.add(*(gdstk.Reference(gcross, origin=c) for c in global_centres))
    top.add(*wafer_outline())
    size = tuple(max(info['size_um'][axis] for info in sources.values()) for axis in (0, 1))
    gwin = args.global_cross_span + 2 * args.global_cross_window_margin
    manifest = {
        'description': f'4-inch wafer on the NTHU TiAu mask grid ({len(dies)} dies, 10 mm @ 10.2 mm pitch).',
        'wafer_diameter_um': WAFER_RADIUS * 2, 'wafer_guide_layer': WAFER_GUIDE,
        'die_slot_size_um': size, 'chip_size_um': [CHIP, CHIP], 'pitch_um': [PITCH, PITCH],
        'street_width_um': PITCH - CHIP, 'grid_origin': 'chip centred at (0, 0)',
        'alignment_cross': {
            'cell': CROSS_CELL, 'layer': METAL, 'span_um': CROSS_SPAN, 'stroke_um': CROSS_STROKE,
            'polarity': 'metal in positive and etch', 'sites': 'every die-corner street intersection',
            'style': 'same plus as the NCKU / standard wafer (assemble_wafer.alignment_cross)',
            'centres_um': [list(c) for c in marks],
        } if marks else None,
        'global_marks': {
            'cell': GLOBAL_CROSS_CELL, 'layer': METAL,
            'span_um': args.global_cross_span, 'stroke_um': args.global_cross_stroke,
            'window_margin_um': args.global_cross_window_margin, 'window_size_um': [gwin, gwin],
            'polarity': 'metal cross inside exposed bare window (high contrast in positive and etch)',
            'centres_um': [list(c) for c in global_centres],
        } if global_centres else None,
        'street_frame': None,
        'assignment_seed': SEED, 'target_counts': target_counts,
        'sources': sources, 'wafer_die_ids': wafer_die_ids, 'counts': counts, 'dies': dies,
    }
    with TemporaryDirectory(dir=out) as temp:
        gds, png, record, id_png = [Path(temp) / p.name for p in outputs]
        lib.write_gds(str(gds), max_points=GDS_MAX_POINTS)
        manifest['checks'] = verify(gds, cells, dies, sources, crosses=marks, global_marks=global_centres)
        check_die_ids(id_png, cells, dies, wafer_die_ids, sources)
        plot(png, dies, counts, cells, active_names=active_names, crosses=bool(marks))
        record.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        for source, target in zip((gds, png, record, id_png), outputs):
            source.replace(target)
    print(json.dumps({'counts': counts['Total'], **manifest['checks']}, indent=2))


if __name__ == '__main__':
    main()
