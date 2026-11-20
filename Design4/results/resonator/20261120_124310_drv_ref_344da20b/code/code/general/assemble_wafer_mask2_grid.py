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
- No 250 um alignment crosses in the streets: the old wafer had none, and the TiAu wafer markers take that role.

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
from assemble_wafer import (COLORS, GDS_MAX_POINTS, INPUTS, METAL, WAFER_GUIDE, WAFER_RADIUS,  # noqa: E402
                            add_die_ids, check_die_ids, load_die, polygon_signature, wafer_outline)

PITCH = 10200.0
CHIP = 10000.0
ROW_HALF_WIDTH = {4: 1, 3: 3, 2: 3, 1: 4, 0: 4, -1: 4, -2: 3, -3: 3, -4: 1}   # old wafer, rows north to south
COUNTS = {'D2': 13, 'D3': 13, 'D5': 13, 'D6': 12, 'Chip4': 5, 'PCM': 5}
ALLSPICE = ('Chip4', 'PCM')
SEED = 20261119


def sites():
    return [(i, j) for j in sorted(ROW_HALF_WIDTH, reverse=True)
            for i in range(-ROW_HALF_WIDTH[j], ROW_HALF_WIDTH[j] + 1)]


def region(i, j):
    """Centre 3 x 3 block, then a pinwheel of four 13-site quadrants."""
    if abs(i) <= 1 and abs(j) <= 1:
        return 'A'
    if i <= 0 < j:
        return 'B'
    if i > 0 <= j:
        return 'C'
    if i < 0 and j <= 0:
        return 'D'
    return 'E'


def assign():
    """Backtracking assignment with exact global counts, balanced regions and the adjacency rules."""
    grid = sites()
    regions = {s: region(*s) for s in grid}
    size = Counter(regions.values())
    # Per region: one Chip4 and one PCM, then the four device designs within one of each other.
    devices = [n for n in COUNTS if n not in ALLSPICE]
    caps = {r: {**{n: 1 for n in ALLSPICE},
                **{n: -(-(size[r] - 2) // len(devices)) for n in devices}} for r in size}
    corners = {(i, j) for i in (-3, 3) for j in (-3, 3)}
    rng = random.Random(SEED)
    placed, total, per_region = {}, Counter(), {r: Counter() for r in size}

    # Step 1: one Chip4 and one PCM per region. Two Allspice dies overlap even diagonally (their
    # corners reach past 5100 um), so no two may be among each other's eight neighbours.
    def near(a, b):
        return max(abs(a[0] - b[0]), abs(a[1] - b[1])) <= 1
    by_region = {r: [s for s in grid if regions[s] == r and s not in corners] for r in sorted(size)}
    for _ in range(200_000):
        chosen = {r: rng.sample(by_region[r], 2) for r in by_region}
        flat = [s for pair in chosen.values() for s in pair]
        if all(not near(a, b) for k, a in enumerate(flat) for b in flat[k + 1:]):
            break
    else:
        raise RuntimeError('No Chip4/PCM placement found')
    for r, (chip4, pcm) in chosen.items():
        for s, name in ((chip4, 'Chip4'), (pcm, 'PCM')):
            placed[s] = name
            total[name] += 1
            per_region[r][name] += 1

    def allowed(s, name):
        r = regions[s]
        if total[name] >= COUNTS[name] or per_region[r][name] >= caps[r][name]:
            return False
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            if placed.get((s[0] + di, s[1] + dj)) == name:
                return False
        return True

    # Step 2: the four device designs, row by row, with backtracking.
    order = [s for s in grid if s not in placed]

    def solve(k):
        if k == len(order):
            return True
        s = order[k]
        names = list(devices)
        rng.shuffle(names)
        names.sort(key=lambda n: total[n] - COUNTS[n])
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


def verify(path, cells, dies, sources):
    result = gdstk.read_gds(str(path))
    top, = result.top_level()
    assert (result.unit, result.precision) == (1e-6, 1e-9)
    assert result.layers_and_datatypes() == {METAL, WAFER_GUIDE}
    assert Counter((r.cell_name, tuple(r.origin)) for r in top.references) == \
        Counter((d['cell'], tuple(d['centre_um'])) for d in dies)
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


def plot(path, dies, counts, cells):
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
    ax.set_title(f'4-inch wafer on the Mask2 grid  |  {len(dies)} dies  |  10 mm @ 10.2 mm', loc='left',
                 weight='bold')
    legend.set_axis_off()
    legend.text(0, .97, 'Dies per region', weight='bold', fontsize=13)
    names = list(COUNTS)
    rows = [*sorted(r for r in counts if r != 'Total'), 'Total']
    table = legend.table(cellText=[[counts[r].get(n, 0) for n in names] for r in rows], colLabels=names,
                         rowLabels=rows, colColours=[COLORS[n] for n in names], cellLoc='center',
                         bbox=[.07, .62, .90, .285])
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    legend.text(0, .52, 'A: centre 3 x 3 chips (dashed)\nB-E: pinwheel quadrants of 13 chips\n\n'
                        'Same 61 chip centres as the laser-cut file\n'
                        'and Mask2 (TiAu markers line up with\nevery 10 mm chip).\n\n'
                        'Metal 1/0  ·  Wafer guide 99/0\nNo street alignment crosses', fontsize=9, va='top')
    fig.savefig(path, dpi=185, facecolor='white')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('folder', type=Path, help='wafer folder (default source folder)')
    parser.add_argument('--out', type=Path, required=True, help='output folder')
    parser.add_argument('--stem', default='summery_mask2grid')
    parser.add_argument('--source-folder', type=Path)
    parser.add_argument('--input', action='append', default=[], metavar='KEY=PATH')
    parser.add_argument('--replace', action='store_true')
    args = parser.parse_args()

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    outputs = [out / f'{args.stem}.{s}' for s in ('gds', 'png', 'json')] + [out / f'{args.stem}_die_id_check.png']
    if not args.replace and any(p.exists() for p in outputs):
        raise SystemExit('Output exists; pass --replace')
    inputs = dict(item.split('=', 1) for item in args.input)
    cells, sources = {}, {}
    for name, filename in INPUTS.items():
        path = Path(inputs[name]) if name in inputs else (args.source_folder or args.folder) / filename
        cells[name], sources[name] = load_die(path, name)
        sources[name]['path_from_assembly'] = Path(relpath(path, out)).as_posix()

    dies, counts = assign()
    wafer_die_ids = add_die_ids(cells, sources, dies)
    lib = gdstk.Library('WAFER', unit=1e-6, precision=1e-9)
    top = lib.new_cell('WAFER_4_INCH')
    lib.add(*cells.values())
    top.add(*(gdstk.Reference(cells[d['cell']], origin=d['centre_um']) for d in dies))
    top.add(*wafer_outline())
    size = tuple(max(info['size_um'][axis] for info in sources.values()) for axis in (0, 1))
    manifest = {
        'description': 'Same dies as the 20261119 summery, re-placed on the old Mask2 / laser-cut chip grid.',
        'wafer_diameter_um': WAFER_RADIUS * 2, 'wafer_guide_layer': WAFER_GUIDE,
        'die_slot_size_um': size, 'chip_size_um': [CHIP, CHIP], 'pitch_um': [PITCH, PITCH],
        'street_width_um': PITCH - CHIP, 'grid_origin': 'chip centred at (0, 0)',
        'alignment_cross': None, 'assignment_seed': SEED, 'target_counts': COUNTS,
        'sources': sources, 'wafer_die_ids': wafer_die_ids, 'counts': counts, 'dies': dies,
    }
    with TemporaryDirectory(dir=out) as temp:
        gds, png, record, id_png = [Path(temp) / p.name for p in outputs]
        lib.write_gds(str(gds), max_points=GDS_MAX_POINTS)
        manifest['checks'] = verify(gds, cells, dies, sources)
        check_die_ids(id_png, cells, dies, wafer_die_ids, sources)
        plot(png, dies, counts, cells)
        record.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        for source, target in zip((gds, png, record, id_png), outputs):
            source.replace(target)
    print(json.dumps({'counts': counts['Total'], **manifest['checks']}, indent=2))


if __name__ == '__main__':
    main()
