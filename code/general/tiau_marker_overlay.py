"""Extract the fixed TiAu alignment markers from the existing photomask and overlay a wafer assembly.

The TiAu marker photomask (``Fluxonium_2Q_bandpass_filter_Mask2.gds``) is already fabricated, so every
marker coordinate is copied unchanged. The script:

1. writes the wafer-level marks (Mask2 layer 3/0) and the chip-level marks (Mask2 layer 2/0) to separate
   GDS files plus one combined file, all on the original layers and coordinates;
2. overlays them on a wafer assembly (``<folder>/summery.gds``, metal 1/0) with both drawn about the
   same wafer centre (0, 0);
3. classifies every marker site by where it lands (die, street, outside the dies) and how much metal
   of the assembly covers it, and renders overview and close-up figures.

Usage (repository root, qmetal):
    python -X utf8 -B code/general/tiau_marker_overlay.py 20261119 --mask <path to Mask2 GDS>
"""

import argparse
import hashlib
import json
import shutil
from collections import defaultdict
from pathlib import Path

import gdstk
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CHIP_LAYER = (2, 0)     # Mask2 chip markers
WAFER_LAYER = (3, 0)    # Mask2 wafer markers (frame_M2)
METAL = (1, 0)
STREET = 100.0
SITE_MARGIN = 30.0      # window around a marker site used for the local metal fraction (um)
CHEESE_MAX = 5.0        # widest opening treated as a cheese hole (um); 5 um holes in the fab GDS


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def flatten_layer(cell, layer):
    return cell.get_polygons(layer=layer[0], datatype=layer[1])


def group_sites(polygons, gap):
    """Cluster polygons whose bounding boxes lie within ``gap`` of each other (one cluster = one site)."""
    boxes = [np.asarray(p.bounding_box()) for p in polygons]
    parent = list(range(len(polygons)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    order = np.argsort([b[0][0] for b in boxes])
    for a_index, a in enumerate(order):
        for b in order[a_index + 1:]:
            if boxes[b][0][0] - boxes[a][1][0] > gap:
                break
            ba, bb = boxes[a], boxes[b]
            if (bb[0][0] - ba[1][0] <= gap and ba[0][0] - bb[1][0] <= gap
                    and bb[0][1] - ba[1][1] <= gap and ba[0][1] - bb[1][1] <= gap):
                parent[find(a)] = find(b)
    groups = defaultdict(list)
    for i in range(len(polygons)):
        groups[find(i)].append(polygons[i])
    sites = []
    for polys in groups.values():
        lo = np.min([p.bounding_box()[0] for p in polys], axis=0)
        hi = np.max([p.bounding_box()[1] for p in polys], axis=0)
        sites.append({'polygons': polys, 'bbox': (lo, hi), 'centre': (lo + hi) / 2})
    sites.sort(key=lambda s: (-round(s['centre'][1], 3), round(s['centre'][0], 3)))
    return sites


def write_gds(path, cell_name, layers):
    lib = gdstk.Library(unit=1e-6, precision=1e-9)
    cell = lib.new_cell(cell_name)
    for polygons in layers:
        cell.add(*[p.copy() for p in polygons])
    lib.write_gds(str(path), max_points=8190)  # one polygon per wafer ring, as in summery.gds


def area(polygons):
    return sum(p.area() for p in polygons)


class MetalIndex:
    """Bounding-box prefilter for the flattened assembly metal."""

    def __init__(self, polygons):
        self.polygons = polygons
        boxes = np.array([p.bounding_box() for p in polygons])
        self.lo, self.hi = boxes[:, 0, :], boxes[:, 1, :]

    def window(self, lo, hi):
        mask = ((self.lo[:, 0] <= hi[0]) & (self.hi[:, 0] >= lo[0])
                & (self.lo[:, 1] <= hi[1]) & (self.hi[:, 1] >= lo[1]))
        return [self.polygons[i] for i in np.flatnonzero(mask)]


def metal_in(index, polygons, lo, hi):
    candidates = index.window(lo, hi)
    if not candidates:
        return 0.0, []
    clipped = gdstk.boolean(polygons, candidates, 'and', precision=1e-3)
    return area(clipped), candidates


def locate(centre, dies, size):
    """Return (die, local xy) if the point lies inside a die slot, otherwise (None, None)."""
    half = np.asarray(size) / 2
    for die in dies:
        local = np.asarray(centre) - die['centre_um']
        if np.all(np.abs(local) <= half):
            return die, local
    return None, None


def classify(site, index, dies, size, die_extent):
    lo, hi = site['bbox']
    marker_area = area(site['polygons'])
    covered, _ = metal_in(index, site['polygons'], lo, hi)
    win_lo, win_hi = lo - SITE_MARGIN, hi + SITE_MARGIN
    window = [gdstk.rectangle(win_lo, win_hi)]
    window_metal, _ = metal_in(index, window, win_lo, win_hi)
    window_fraction = window_metal / window[0].area()
    die, local = locate(site['centre'], dies, size)
    record = {
        'centre_um': [round(float(v), 3) for v in site['centre']],
        'bbox_um': [round(float(v), 3) for v in (*lo, *hi)],
        'marker_area_um2': round(marker_area, 1),
        'metal_over_marker_fraction': round(covered / marker_area, 4),
        'metal_in_window_fraction': round(window_fraction, 4),
        'window_margin_um': SITE_MARGIN,
    }
    if die is None:
        radius = float(np.hypot(*site['centre']))
        record.update(location='outside_dies' if radius < 50_800 - 1_000 else 'wafer_edge',
                      die_id=None, design=None, local_um=None)
    else:
        record.update(die_id=die['die_id'], design=die['design'],
                      local_um=[round(float(v), 3) for v in local])
        # Source dies are 10 x 10 mm inside the larger slot: the rim between them is unpatterned.
        ext = die_extent[die['design']]
        inside = (ext[0][0] <= local[0] <= ext[1][0]) and (ext[0][1] <= local[1] <= ext[1][1])
        record['location'] = 'die' if inside else 'slot_rim'
    # Openings narrower than CHEESE_MAX (the 5 um cheese holes) vanish under an erosion by half that
    # width; what survives is a real opening (CPW gap, trace clearance, die edge, unpatterned area).
    openings = gdstk.boolean(window, index.window(win_lo, win_hi), 'not', precision=1e-3)
    real = gdstk.offset(openings, -CHEESE_MAX / 2 - 0.1, precision=1e-3) if openings else []
    record['opening_area_in_window_um2'] = round(area(openings), 1)
    record['opening_area_after_erosion_um2'] = round(area(real), 1)
    frac, wfrac = record['metal_over_marker_fraction'], record['metal_in_window_fraction']
    if wfrac == 0:
        record['status'] = 'clear'                 # nothing patterned within the window
    elif not real:
        record['status'] = 'ground_plane'          # solid ground, at most cheese holes
    else:
        record['status'] = 'conflict'              # CPW gap, trace, die edge or other feature nearby
    return record


def plot_tiau_summery(path, chip_sites, wafer_sites):
    """Wafer view of the TiAu marker layer alone, with close-ups of one wafer site and one chip site."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Rectangle, Polygon as MplPolygon

    fig = plt.figure(figsize=(14, 9), layout='constrained')
    grid = fig.add_gridspec(2, 2, width_ratios=[1.6, 1])
    ax = fig.add_subplot(grid[:, 0])
    ax.add_patch(Circle((0, 0), 50.8, fill=False, ec='#26394b', lw=1.4))
    for i in range(-4, 5):
        for j in range(-4, 5):
            ax.add_patch(Rectangle((i * 10.2 - 5, j * 10.2 - 5), 10, 10, fill=False, ec='#9a9a9a',
                                   lw=.5, ls=':'))
    for site in chip_sites:
        ax.scatter(*(site['centre'] / 1000), s=14, marker='P', color='#c98a00', edgecolors='none')
    for site in wafer_sites:
        ax.scatter(*(site['centre'] / 1000), s=110, marker='X', color='#8a3b12', edgecolors='none')
    ax.set(xlim=(-53, 53), ylim=(-53, 53), aspect='equal', xlabel='x (mm)', ylabel='y (mm)')
    ax.set_title(f'TiAu marker summery  |  {len(wafer_sites)} wafer sites (3/0)  |  '
                 f'{len(chip_sites)} chip sites (2/0)', loc='left', weight='bold')
    ax.legend(handles=[plt.Line2D([], [], ls='', marker='X', ms=10, color='#8a3b12', label='wafer marker site'),
                       plt.Line2D([], [], ls='', marker='P', ms=7, color='#c98a00', label='chip marker site'),
                       plt.Line2D([], [], ls=':', color='#9a9a9a', label='Mask2 die grid (10 mm @ 10.2 mm)')],
              loc='lower left', fontsize=9)

    def closeup(axis, site, pad, title):
        lo, hi = site['bbox']
        c = site['centre']
        r = max(hi - lo) / 2 + pad
        for p in site['polygons']:
            axis.add_patch(MplPolygon(p.points, fc='#e3a21a', ec='#7a4f00', lw=.5))
        axis.set(xlim=(c[0] - r, c[0] + r), ylim=(c[1] - r, c[1] + r), aspect='equal',
                 xlabel='x (µm)', ylabel='y (µm)')
        axis.set_title(title, fontsize=10, loc='left')

    upper = [s for s in wafer_sites if s['centre'][1] > 0]
    big = max(upper, key=lambda s: area(s['polygons']))
    near = [s for s in wafer_sites if abs(s['centre'][0] - big['centre'][0]) < 1 and s is not big
            and abs(s['centre'][1] - big['centre'][1]) < 3000]
    merged = {'polygons': big['polygons'] + sum((s['polygons'] for s in near), [])}
    lo = np.min([p.bounding_box()[0] for p in merged['polygons']], axis=0)
    hi = np.max([p.bounding_box()[1] for p in merged['polygons']], axis=0)
    merged.update(bbox=(lo, hi), centre=(lo + hi) / 2)
    closeup(fig.add_subplot(grid[0, 1]), merged, 300,
            f"Wafer marker site at ({merged['centre'][0]:.0f}, {merged['centre'][1]:.0f}) µm")
    chip = min(chip_sites, key=lambda s: np.hypot(*(s['centre'] - (2529.25, 4102.0))))
    closeup(fig.add_subplot(grid[1, 1]), chip, 15,
            f"Chip marker site at ({chip['centre'][0]:.0f}, {chip['centre'][1]:.0f}) µm")
    fig.savefig(path, dpi=170, facecolor='white')
    plt.close(fig)


def plot_overview(path, dies, size, chip_sites, wafer_sites, records, title='summery'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Rectangle
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    from assemble_wafer import COLORS

    colours = {'clear': '#2a9d4b', 'ground_plane': '#3769c4', 'conflict': '#d1342f'}
    fig, ax = plt.subplots(figsize=(11, 11), layout='constrained')
    ax.add_patch(Circle((0, 0), 50.8, fill=False, ec='#26394b', lw=1.4))
    w, h = (v / 1000 for v in size)
    for die in dies:
        x, y = (v / 1000 for v in die['centre_um'])
        ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, fc=COLORS[die['design']], ec='#243445',
                               lw=.3, alpha=.55))
        ax.text(x, y, f"{die['design']}\n{die['die_id']}", ha='center', va='center', fontsize=7,
                color='#152436')
    # Old die outlines from the Mask2 grid (10 mm dies on a 10.2 mm pitch).
    for i in range(-4, 5):
        for j in range(-4, 5):
            ax.add_patch(Rectangle((i * 10.2 - 5, j * 10.2 - 5), 10, 10, fill=False, ec='#7a7a7a',
                                   lw=.5, ls=':'))
    for kind, sites, marker, size_pt in (('chip', chip_sites, 'P', 34), ('wafer', wafer_sites, 'X', 120)):
        for site, rec in zip(sites, records[kind]):
            x, y = site['centre'] / 1000
            ax.scatter(x, y, s=size_pt, marker=marker, color=colours[rec['status']],
                       edgecolors='black', linewidths=.4, zorder=5)
    handles = [plt.Line2D([], [], ls='', marker='P', ms=8, mfc=c, mec='black', label=k.replace('_', ' '))
               for k, c in colours.items()]
    handles += [plt.Line2D([], [], ls='', marker='X', ms=11, mfc='white', mec='black',
                           label='wafer marker site (Mask2 3/0)'),
                plt.Line2D([], [], ls='', marker='P', ms=8, mfc='white', mec='black',
                           label='chip marker site (Mask2 2/0)'),
                plt.Line2D([], [], ls=':', color='#7a7a7a', label='Mask2 die grid (10 mm @ 10.2 mm)')]
    ax.legend(handles=handles, loc='lower left', fontsize=9, framealpha=.95)
    ax.set(xlim=(-53, 53), ylim=(-53, 53), aspect='equal', xlabel='x (mm)', ylabel='y (mm)')
    ax.set_title(f'TiAu markers (Mask2, fixed) on {title} dies', loc='left', weight='bold')
    fig.savefig(path, dpi=170, facecolor='white')
    plt.close(fig)


def plot_closeups(path, sites, records, index, title, ncols=6, max_sites=36, half=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon

    chosen = list(zip(sites, records))[:max_sites]
    nrows = max(1, -(-len(chosen) // ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.6 * ncols, 2.75 * nrows), layout='constrained',
                             squeeze=False)
    for ax in axes.flat:
        ax.set_axis_off()
    for ax, (site, rec) in zip(axes.flat, chosen):
        lo, hi = site['bbox']
        c = site['centre']
        r = half if half is not None else max(hi - lo) / 2 + SITE_MARGIN + 40
        win_lo, win_hi = c - r, c + r
        ax.add_patch(plt.Rectangle(win_lo, 2 * r, 2 * r, fc='#f2f2f2', ec='none'))
        for p in index.window(win_lo, win_hi):
            ax.add_patch(MplPolygon(p.points, fc='#9aa9b8', ec='none'))
        for p in site['polygons']:
            ax.add_patch(MplPolygon(p.points, fc='#e3a21a', ec='#7a4f00', lw=.4, alpha=.9))
        ax.set(xlim=(win_lo[0], win_hi[0]), ylim=(win_lo[1], win_hi[1]), aspect='equal')
        ax.set_axis_on()
        ax.tick_params(labelsize=5)
        where = rec['die_id'] + ' ' + rec['design'] if rec['die_id'] else rec['location']
        ax.set_title(f"{where}\n({c[0]:.0f}, {c[1]:.0f}) {rec['status'].replace('_', ' ')}",
                     fontsize=6.5, color={'clear': '#1d6b33', 'ground_plane': '#24498a'}
                     .get(rec['status'], '#a0221e'))
    fig.suptitle(title + '   (grey: summery metal 1/0, gold: TiAu marker)', fontsize=10, x=0.01, ha='left')
    fig.savefig(path, dpi=200, facecolor='white')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('folder', help='wafer folder with summery.gds and summery.json, e.g. 20261119')
    parser.add_argument('--mask', required=True, type=Path, help='existing TiAu marker photomask GDS (Mask2)')
    parser.add_argument('--out', default='TiAugds', help='output subfolder')
    parser.add_argument('--stem', default='summery',
                        help='assembly stem: reads <assembly-dir>/<stem>.gds and <stem>.json')
    parser.add_argument('--assembly-dir', type=Path, help='folder of the assembly (default: the wafer folder)')
    parser.add_argument('--replace', action='store_true')
    args = parser.parse_args()

    folder = ROOT / args.folder
    out = folder / args.out
    stem = args.stem
    prefix = '' if stem == 'summery' else f'{stem}_'   # keep the original names for summery
    assembly_dir = ROOT / args.assembly_dir if args.assembly_dir else folder
    overlay_gds = out / f'{stem}_with_TiAu_marker.gds'
    if overlay_gds.exists() and not args.replace:
        raise SystemExit(f'{overlay_gds} exists; pass --replace or choose --out')
    (out / 'inputs').mkdir(parents=True, exist_ok=True)

    mask_copy = out / 'inputs' / 'Fluxonium_2Q_bandpass_filter_Mask2.gds'
    if not (mask_copy.exists() and mask_copy.resolve() == args.mask.resolve()):
        shutil.copyfile(args.mask, mask_copy)
    mask_lib = gdstk.read_gds(str(mask_copy))
    mask_top, = mask_lib.top_level()
    chip_polys = flatten_layer(mask_top, CHIP_LAYER)
    wafer_polys = flatten_layer(mask_top, WAFER_LAYER)
    other = {(p.layer, p.datatype) for p in mask_top.get_polygons()} - {CHIP_LAYER, WAFER_LAYER}
    if other:
        raise SystemExit(f'Mask2 has unexpected layers {sorted(other)}')
    # Every per-chip marker cell must be identical, otherwise the chip markers differ by die type.
    sigs = {c.name: sorted((p.layer, p.datatype, tuple(np.round(p.points, 3).ravel())) for p in c.polygons)
            for c in mask_lib.cells if c.name.startswith('TOP_MAIN_2')}
    if len({json.dumps(v) for v in sigs.values()}) != 1:
        raise SystemExit('TOP_MAIN_2* marker cells differ')

    chip_sites = group_sites(chip_polys, gap=40.0)
    wafer_sites = group_sites(wafer_polys, gap=600.0)

    wafer_gds = out / 'TiAu_marker_wafer.gds'
    chip_gds = out / 'TiAu_marker_chip.gds'
    both_gds = out / 'TiAu_marker.gds'
    write_gds(wafer_gds, 'TIAU_MARKER_WAFER', [wafer_polys])
    write_gds(chip_gds, 'TIAU_MARKER_CHIP', [chip_polys])
    write_gds(both_gds, 'TIAU_MARKER', [chip_polys, wafer_polys])

    # Wafer-level TiAu summery: the markers plus the same 4-inch wafer guide (99/0) as summery.gds.
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    from assemble_wafer import wafer_outline
    tiau_summery_gds = out / 'summery_TiAu_marker.gds'
    write_gds(tiau_summery_gds, 'WAFER_4_INCH_TIAU', [chip_polys, wafer_polys, wafer_outline()])
    plot_tiau_summery(out / 'summery_TiAu_marker.png', chip_sites, wafer_sites)

    # Round-trip check: written polygons equal the flattened Mask2 polygons on each layer.
    def signature(polys):
        return sorted(tuple(np.round(p.points, 3).ravel()) for p in polys)
    rt = gdstk.read_gds(str(both_gds)).top_level()[0]
    for layer, ref in ((CHIP_LAYER, chip_polys), (WAFER_LAYER, wafer_polys)):
        if signature(flatten_layer(rt, layer)) != signature(ref):
            raise SystemExit(f'Round-trip mismatch on {layer}')

    summery_gds = assembly_dir / f'{stem}.gds'
    manifest = json.loads((assembly_dir / f'{stem}.json').read_text(encoding='utf-8'))
    dies = manifest['dies']
    for die in dies:
        die['centre_um'] = np.asarray(die['centre_um'], float)
    size = manifest['die_slot_size_um']
    lib = gdstk.read_gds(str(summery_gds))
    top, = lib.top_level()
    cells = {c.name: c for c in lib.cells}
    die_extent = {name: np.asarray(cells[name].bounding_box()) for name in manifest['sources']}
    metal = flatten_layer(top, METAL)
    index = MetalIndex(metal)

    records = {'chip': [classify(s, index, dies, size, die_extent) for s in chip_sites],
               'wafer': [classify(s, index, dies, size, die_extent) for s in wafer_sites]}

    # Overlay GDS: summery as it is, plus the markers on their Mask2 layers in a separate cell.
    marker_cell = gdstk.Cell('TIAU_MARKER')
    marker_cell.add(*[p.copy() for p in chip_polys], *[p.copy() for p in wafer_polys])
    lib.add(marker_cell)
    top.add(gdstk.Reference(marker_cell))
    lib.write_gds(str(overlay_gds), max_points=8190)

    plot_overview(out / f'{stem}_with_TiAu_marker.png', dies, size, chip_sites, wafer_sites, records, title=stem)
    in_die = [(s, r) for s, r in zip(chip_sites, records['chip']) if r['status'] == 'conflict']
    plot_closeups(out / f'{prefix}chip_marker_conflicts.png', [s for s, _ in in_die], [r for _, r in in_die],
                  index, f'Mask2 chip marker sites in conflict with {stem} patterns', ncols=6, max_sites=200)
    plot_closeups(out / f'{prefix}wafer_marker_sites.png', wafer_sites, records['wafer'], index,
                  'Mask2 wafer marker sites', ncols=4, max_sites=8)

    def tally(recs):
        out_counts = defaultdict(int)
        for r in recs:
            out_counts[f"{r['location']}/{r['status']}"] += 1
        return dict(sorted(out_counts.items()))

    report = {
        'description': 'Fixed TiAu marker mask (Mask2) extracted unchanged and overlaid on the wafer assembly.',
        'assumption': f'Mask2 and {stem}.gds share the wafer centre (0, 0) and orientation; no offset applied.',
        'mask': {'source_path': str(args.mask), 'snapshot': mask_copy.relative_to(folder).as_posix(),
                 'sha256': sha256(mask_copy),
                 'chip_marker_layer': CHIP_LAYER, 'wafer_marker_layer': WAFER_LAYER,
                 'chip_marker_polygons': len(chip_polys), 'wafer_marker_polygons': len(wafer_polys),
                 'chip_marker_sites': len(chip_sites), 'wafer_marker_sites': len(wafer_sites),
                 'mask2_die_pitch_um': 10200.0, 'mask2_die_size_um': 10000.0},
        'assembly': {'gds': summery_gds.relative_to(folder).as_posix(), 'sha256': sha256(summery_gds), 'metal_layer': METAL,
                     'die_slot_size_um': size, 'pitch_um': manifest['pitch_um']},
        'outputs': {p.relative_to(folder).as_posix(): sha256(p)
                    for p in (wafer_gds, chip_gds, both_gds, overlay_gds, tiau_summery_gds)},
        'classification': {
            'clear': 'no assembly metal within the window around the site',
            'ground_plane': 'site inside ground; only openings no wider than cheese_max_um (cheese holes) in the window',
            'conflict': 'a larger opening (CPW gap, trace clearance, die edge, unpatterned area) within the window',
            'cheese_max_um': CHEESE_MAX,
            'window_margin_um': SITE_MARGIN},
        'summary': {'chip': tally(records['chip']), 'wafer': tally(records['wafer'])},
        'sites': records,
    }
    (out / f'{prefix}tiau_marker_overlay.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report['summary'], indent=2))


if __name__ == '__main__':
    main()
