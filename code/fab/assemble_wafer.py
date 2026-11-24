"""Assemble the six fabrication layouts on one circular 4-inch wafer."""

import argparse
from collections import Counter
import hashlib
import json
from math import hypot
from os.path import relpath
from pathlib import Path
from tempfile import TemporaryDirectory

import gdstk
import numpy as np


INPUTS = {
    'D2': 'Design2_2026-11-16_wafer_r8_fab.gds',
    'D3': 'Design3_3_2026-11-16_wafer_r7_fab.gds',   # main FTF (former Design3_3)
    'D5': 'Design5_2026-11-16_wafer_r7_fab.gds',
    'D6': 'Design6_2026-11-16_wafer_r9_fab.gds',
    'Chip4': 'AllspiceR3_v5p4_Chip4.oas',
    'PCM': 'AllspiceR3_v5p4_PCMdie_M1.oas',
}
COLORS = dict(zip(INPUTS, ('#80bfd0', '#92b4eb', '#ecc586', '#aacbaf', '#cfb1d7', '#e4aaa6')))

# Rows are north to south. Each region has one PCM and one Chip4. The four
# other designs share the remaining sites. Dashes mark the circular edge.
ROWS = (
    '-     -     D2    D3    D6    D3    -     -',
    '-     D5    PCM   Chip4 D2    D5    D2    -',
    'D6    D3    D2    D5    D3    D6    Chip4 PCM',
    'D2    D5    D6    PCM   D2    D5    D3    D6',
    'D6    PCM   D3    Chip4 D6    D3    D2    D5',
    'D5    D2    D5    D2    D3    D2    D6    D3  ',
    '-     D6    D3    Chip4 D6    D5    PCM   -',
    '-     -     D2    D5    D3    Chip4 -     -',
)
DEVICE_KEYS = ('D2', 'D3', 'D5', 'D6')
# Balanced 52-die devices-only grid (13 each of D2, D3, D5, D6; no adjacent duplicates;
# 4 each in centre region A and 2-3 each in outer regions B-E).
ROWS_DEVICES_ONLY = (
    '-     -     D2    D3    D5    D6    -     -',
    '-     D6    D5    D6    D2    D3    D5    -',
    'D2    D3    D2    D3    D5    D6    D2    D3',
    'D5    D6    D5    D6    D2    D3    D5    D6',
    'D2    D3    D2    D3    D5    D6    D2    D3',
    'D5    D6    D5    D6    D2    D3    D5    D6',
    '-     D3    D2    D3    D5    D6    D2    -',
    '-     -     D5    D6    D2    D3    -     -',
)

WAFER_RADIUS = 50_800.0  # 4 inches in diameter
EDGE_MARGIN = 1_000.0
STREET = 100.0
CROSS_SPAN = 250.0
CROSS_STROKE = 5.0
SOURCE_CROSS_COUNTS = {'Chip4': 5, 'PCM': 9}
METAL = (1, 0)
WAFER_GUIDE = (99, 0)
GDS_MAX_POINTS = 8190  # accommodates the closed polygon in one XY record
DIE_ID_SIZE = 150.0   # wafer die ID height (um); 100 um before 2026-11-18
DIE_ID_INSET = 50.0


def alignment_cross():
    """The same thin plus as the original Allspice alignment markers."""
    arm, half_width = CROSS_SPAN / 2, CROSS_STROKE / 2
    horizontal = gdstk.rectangle((-arm, -half_width), (arm, half_width))
    vertical = gdstk.rectangle((-half_width, -arm), (half_width, arm))
    cross, = gdstk.boolean(horizontal, vertical, 'or', precision=.001, layer=METAL[0])
    return cross


def cross_centre(polygon, vertices):
    """Identify only the exact thin plus, independent of vertex ordering."""
    if len(polygon.points) != len(vertices):
        return None
    low, high = np.asarray(polygon.bounding_box())
    if not np.allclose(high - low, CROSS_SPAN, rtol=0, atol=.001):
        return None
    centre = (low + high) / 2
    normalized = {tuple(row) for row in np.rint((polygon.points - centre) * 1000)}
    return centre if normalized == vertices else None


def load_die(path, name):
    contents = path.read_bytes()
    oasis = contents.startswith(b'%SEMI-OASIS\r\n')
    read = gdstk.read_oas if oasis else gdstk.read_gds
    lib = read(str(path))
    if (lib.unit, lib.precision) != (1e-6, 1e-9):
        raise ValueError(f'{path.name}: expected coordinates in µm on a 1 nm grid')
    source, = lib.top_level()
    if source.get_labels():
        raise ValueError(f'{path.name}: cannot map GDS labels to metal')
    polygons = source.get_polygons()
    expected_spec = (2, 0) if name == 'Chip4' else METAL
    if {(p.layer, p.datatype) for p in polygons} != {expected_spec}:
        raise ValueError(f'{path.name}: unexpected source layers for the M1 mask')

    low, high = np.asarray(source.bounding_box())
    centre = (low + high) / 2
    cell = gdstk.Cell(name)
    cross_vertices = {tuple(row) for row in np.rint(alignment_cross().points * 1000)}
    removed_crosses = []
    for polygon in polygons:
        if name in SOURCE_CROSS_COUNTS:
            old_centre = cross_centre(polygon, cross_vertices)
            if old_centre is not None:
                removed_crosses.append(np.round(old_centre, 3).tolist())
                continue
        area_before = polygon.area()
        polygon.translate(*(-centre))
        polygon.layer, polygon.datatype = METAL
        if len(polygon.points) > GDS_MAX_POINTS:
            parts = polygon.fracture(max_points=GDS_MAX_POINTS, precision=.001)
            if abs(sum(p.area() for p in parts) - area_before) > 1:
                raise ValueError(f'{path.name}: polygon fracture changed source metal area')
            cell.add(*parts)
        else:
            cell.add(polygon)
    if len(removed_crosses) != SOURCE_CROSS_COUNTS.get(name, 0):
        raise ValueError(f'{path.name}: unexpected number of Allspice alignment marks')
    return cell, {
        'file': path.name,
        'format': 'OASIS' if oasis else 'GDSII',
        'sha256': hashlib.sha256(contents).hexdigest(),
        'source_metal_layer': expected_spec,
        'output_metal_layer': METAL,
        'source_centre_um': np.round(centre, 3).tolist(),
        'size_um': np.round(high - low, 3).tolist(),
        'removed_alignment_cross_centres_um': sorted(removed_crosses),
    }


def place_dies(size, *, devices_only=False):
    width, height = size
    pitch_x, pitch_y = width + STREET, height + STREET
    rows_grid = ROWS_DEVICES_ONLY if devices_only else ROWS
    active_keys = DEVICE_KEYS if devices_only else tuple(INPUTS)
    occupied = {}
    for j, row in enumerate(rows_grid):
        for i, name in enumerate(row.split()):
            col, row_index = i - 4, 3 - j
            x, y = (col + .5) * pitch_x, (row_index + .5) * pitch_y
            corner_radius = hypot(abs(x) + width / 2, abs(y) + height / 2)
            fits = corner_radius <= WAFER_RADIUS - EDGE_MARGIN
            if fits != (name != '-'):
                raise ValueError(f'Grid position {(col, row_index)} has unexpected wafer coverage')
            if not fits:
                continue
            if -2 <= col <= 1 and -2 <= row_index <= 1:
                region = 'A'
            else:
                region = ('B' if col < 0 else 'C') if row_index >= 0 else ('D' if col < 0 else 'E')
            occupied[col, row_index] = {'design': name, 'region': region, 'centre_um': [x, y]}

    counts = {r: Counter(d['design'] for d in occupied.values() if d['region'] == r) for r in 'ABCDE'}
    counts['Total'] = Counter(d['design'] for d in occupied.values())
    for region, count in counts.items():
        if devices_only:
            device_counts = [count[name] for name in active_keys]
            if set(count) != set(active_keys) or max(device_counts) - min(device_counts) > 1:
                raise ValueError(f'Unbalanced assignment in {region}: {count}')
        else:
            pcm_chip4_count = 5 if region == 'Total' else 1
            device_counts = [count[name] for name in INPUTS if name not in ('Chip4', 'PCM')]
            if (set(count) != set(INPUTS)
                    or any(count[name] != pcm_chip4_count for name in ('Chip4', 'PCM'))
                    or max(device_counts) - min(device_counts) > 1):
                raise ValueError(f'Unbalanced assignment in {region}: {count}')
    corners = Counter()
    for (i, j), die in occupied.items():
        for di, dj in ((1, 0), (0, 1)):
            if occupied.get((i + di, j + dj), {}).get('design') == die['design']:
                raise ValueError(f'Adjacent copies of {die["design"]} at {(i, j)}')
        corners.update((i + di, j + dj) for di in (0, 1) for dj in (0, 1))
    marks = sorted((i * pitch_x, j * pitch_y) for (i, j), count in corners.items() if count >= 2)
    dies = sorted(occupied.values(), key=lambda die: (-die['centre_um'][1], die['centre_um'][0]))
    for index, die in enumerate(dies):
        die['die_id'] = f'D{index:02d}'
        die['cell'] = f'WAFER_{die["die_id"]}'
    return dies, counts, marks


def die_id_polygons(text, size):
    """Place polygon text 50 µm inside the actual die's upper-right bounds."""
    polygons = gdstk.text(text, DIE_ID_SIZE, (0, 0), layer=METAL[0])
    right = max(p.bounding_box()[1][0] for p in polygons)
    top = max(p.bounding_box()[1][1] for p in polygons)
    offset = np.round(np.asarray(size) / 2 - DIE_ID_INSET - (right, top), 3)
    return [gdstk.Polygon(np.round(polygon.points + offset, 3), layer=METAL[0])
            for polygon in polygons]


def add_die_ids(cells, sources, dies):
    """Keep shared die geometry and put a unique polygon ID in each placement."""
    windows, styles, glyphs = {}, {}, {}
    for die in dies:
        glyphs[die['die_id']] = die_id_polygons(die['die_id'], sources[die['design']]['size_um'])
    for name in sources:
        boxes = [p.bounding_box() for die in dies if die['design'] == name
                 for p in glyphs[die['die_id']]]
        low = np.floor((np.min([b[0] for b in boxes], axis=0) - 5) * 1000) / 1000
        high = np.ceil((np.max([b[1] for b in boxes], axis=0) + 5) * 1000) / 1000
        window = gdstk.rectangle(low, high, layer=METAL[0])
        windows[name] = [*low.tolist(), *high.tolist()]
        cell = cells[name]
        overlaps = []
        for polygon in cell.polygons:
            a, b = polygon.bounding_box()
            if a[0] < high[0] and b[0] > low[0] and a[1] < high[1] and b[1] > low[1]:
                overlaps.append(polygon)
        metal = gdstk.boolean(overlaps, window, 'and', precision=.001, layer=METAL[0])
        area = sum(p.area() for p in metal)
        if abs(area - window.area()) < .01:
            styles[name] = 'etched_in_ground'
        elif area < .01:
            styles[name] = 'metal_on_blank'
        else:
            raise ValueError(f'{name}: upper-right die-ID window contains existing features')
        for polygon in overlaps:
            cell.remove(polygon)
            cell.add(*gdstk.boolean(polygon, window, 'not', precision=.001, layer=METAL[0]))

    for die in dies:
        name = die['design']
        text = glyphs[die['die_id']]
        wrapper = gdstk.Cell(die['cell'])
        wrapper.add(gdstk.Reference(cells[name]))
        if styles[name] == 'etched_in_ground':
            x0, y0, x1, y1 = windows[name]
            patch = gdstk.boolean(gdstk.rectangle((x0, y0), (x1, y1)), text, 'not',
                                  precision=.001, layer=METAL[0])
        else:
            patch = text
        wrapper.add(*patch)
        cells[wrapper.name] = wrapper
        die['die_id_bounds_um'] = [
            min(p.bounding_box()[0][0] for p in text), min(p.bounding_box()[0][1] for p in text),
            max(p.bounding_box()[1][0] for p in text), max(p.bounding_box()[1][1] for p in text),
        ]
    return {'order': 'left_to_right_then_top_to_bottom', 'first': dies[0]['die_id'],
            'last': dies[-1]['die_id'], 'font_size_um': DIE_ID_SIZE,
            'top_right_inset_um': DIE_ID_INSET, 'layer': METAL,
            'windows_um': windows, 'styles': styles}


def wafer_outline():
    # An even polygonal ring includes all four cardinal extrema exactly.
    angles = np.arange(2048) * (2 * np.pi / 2048)
    unit_circle = np.column_stack((np.cos(angles), np.sin(angles)))
    outer = gdstk.Polygon(WAFER_RADIUS * unit_circle)
    inner = gdstk.Polygon((WAFER_RADIUS - 2) * unit_circle)
    return gdstk.boolean(outer, inner, 'not', precision=.001, layer=WAFER_GUIDE[0])


def polygon_signature(cell):
    return Counter((p.layer, p.datatype, np.rint(p.points * 1000).astype('int64').tobytes())
                   for p in cell.polygons)


def verify_gds(path, cells, dies, marks, size, sources, die_id_settings):
    result = gdstk.read_gds(str(path))
    top, = result.top_level()
    assert (result.unit, result.precision) == (1e-6, 1e-9)
    assert result.layers_and_datatypes() == {METAL, WAFER_GUIDE}
    assert top.bounding_box() == ((-WAFER_RADIUS, -WAFER_RADIUS), (WAFER_RADIUS, WAFER_RADIUS))
    expected = Counter((d.get('cell', d['design']), tuple(d['centre_um'])) for d in dies)
    expected.update(('ALIGNMENT_CROSS', origin) for origin in marks)
    assert Counter((r.cell_name, tuple(r.origin)) for r in top.references) == expected
    assert all(r.rotation == 0 and r.magnification == 1 and not r.x_reflection for r in top.references)
    for cell in result.cells:
        if cell.name in cells:
            assert polygon_signature(cell) == polygon_signature(cells[cell.name]), cell.name
            low, high = np.asarray(cell.bounding_box())
            assert np.all(low >= -np.asarray(size) / 2 - .001)
            assert np.all(high <= np.asarray(size) / 2 + .001)
    by_name = {cell.name: cell for cell in result.cells}
    for die in dies:
        cell = by_name[die['cell']]
        assert len(cell.references) == 1 and cell.references[0].cell_name == die['design']
        name = die['design']
        label = die_id_polygons(die['die_id'], sources[name]['size_um'])
        if die_id_settings['styles'][name] == 'etched_in_ground':
            x0, y0, x1, y1 = die_id_settings['windows_um'][name]
            expected_patch = gdstk.boolean(gdstk.rectangle((x0, y0), (x1, y1)), label, 'not',
                                          precision=.001, layer=METAL[0])
        else:
            expected_patch = label
        delta = gdstk.boolean(cell.polygons, expected_patch, 'xor', precision=.001)
        assert sum(p.area() for p in delta) < .01, die['die_id']
    cross = next(cell for cell in result.cells if cell.name == 'ALIGNMENT_CROSS')
    reach, half_stroke = CROSS_SPAN / 2, CROSS_STROKE / 2
    assert cross.bounding_box() == ((-reach, -reach), (reach, reach))
    assert cross.polygons[0].area() == 2 * CROSS_SPAN * CROSS_STROKE - CROSS_STROKE ** 2
    # The arms can extend past the junction square along a street. Test each
    # bar against each die slot to ensure neither arm enters a die.
    half_width, half_height = np.asarray(size) / 2
    for mx, my in marks:
        for die in dies:
            x, y = die['centre_um']
            dx, dy = abs(mx - x), abs(my - y)
            assert dx >= half_width + reach or dy >= half_height + half_stroke
            assert dx >= half_width + half_stroke or dy >= half_height + reach


def plot_wafer(path, size, dies, counts, marks):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Polygon, Rectangle

    width, height = (value / 1000 for value in size)
    radius = WAFER_RADIUS / 1000
    fig, (ax, legend) = plt.subplots(1, 2, figsize=(12.4, 8),
                                    gridspec_kw={'width_ratios': [2.35, 1]}, layout='constrained')
    ax.add_patch(Circle((0, 0), radius, fill=False, edgecolor='#26394b', linewidth=1.4))
    a_x, a_y = 2 * (width + STREET / 1000), 2 * (height + STREET / 1000)
    ax.add_patch(Rectangle((-a_x, -a_y), 2 * a_x, 2 * a_y, fill=False, ls='--', ec='#526274'))
    for start, end in (((0, a_y), (0, radius)), ((0, -a_y), (0, -radius)),
                       ((a_x, 0), (radius, 0)), ((-a_x, 0), (-radius, 0))):
        ax.plot((start[0], end[0]), (start[1], end[1]), '--', color='#526274', lw=.8)
    for die in dies:
        x, y = (v / 1000 for v in die['centre_um'])
        name = die['design']
        ax.add_patch(Rectangle((x - width / 2, y - height / 2), width, height,
                               facecolor=COLORS[name], edgecolor='#243445', linewidth=.3))
        ax.text(x, y, f'{name}\n{die.get("die_id", "")}', ha='center', va='center', fontsize=7.5,
                color='#152436')
    for label, origin in zip('ABCDE', ((0, 0), (-33.5, 32.5), (33.5, 32.5),
                                       (-33.5, -32.5), (33.5, -32.5))):
        ax.annotate(label, origin, ha='center', va='center', weight='bold', fontsize=10,
                    bbox={'boxstyle': 'circle,pad=.26', 'fc': 'white', 'ec': '#526274', 'lw': .6})
    ax.set(xlim=(-53, 53), ylim=(-53, 53), aspect='equal', xlabel='x (mm)', ylabel='y (mm)')
    ax.set_title(f'4-inch wafer  |  Ø101.6 mm  |  {len(dies)} dies', loc='left', weight='bold')
    legend.set_axis_off()
    legend.text(0, .97, 'Dies per region', weight='bold', fontsize=13)
    names = [n for n in INPUTS if counts['Total'].get(n, 0) > 0]
    table = legend.table(cellText=[[counts[r][name] for name in names] for r in (*'ABCDE', 'Total')],
                         colLabels=names, rowLabels=[*'ABCDE', 'Total'],
                         colColours=[COLORS[name] for name in names], cellLoc='center',
                         bbox=[.07, .62, .90, .285])
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    for cell in table.get_celld().values():
        cell.set_edgecolor('#c8cfd5')
    legend.text(0, .55, 'Shared street intersection', weight='bold', fontsize=12)
    inset = legend.inset_axes([.12, .09, .75, .43])
    half = STREET / 2
    corner_size = 175
    for x in (-half - corner_size, half):
        for y in (-half - corner_size, half):
            inset.add_patch(Rectangle((x, y), corner_size, corner_size, color='#d0d8de', ec='none'))
    inset.add_patch(Polygon(alignment_cross().points, fc='#26394b', ec='none'))
    inset.annotate('', xy=(-half, 171), xytext=(half, 171),
                   arrowprops={'arrowstyle': '<->', 'color': '#26394b'})
    inset.text(0, 198, f'{STREET:g} µm', ha='center', fontsize=9)
    inset.set(xlim=(-210, 210), ylim=(-210, 240), aspect='equal')
    inset.set_axis_off()
    legend.text(0, .065, f'{len(marks)} metal crosses, {CROSS_SPAN:g} × {CROSS_SPAN:g} µm, '
                         f'{CROSS_STROKE:g} µm strokes\n'
                         'Metal and crosses 1/0   ·   Wafer guide 99/0', fontsize=9)
    fig.savefig(path, dpi=185, facecolor='white')
    plt.close(fig)


def check_die_ids(path, cells, dies, wafer_die_ids, sources):
    """Assert unique, ordered die IDs and draw a contact sheet: each die's die ID (upper right) next to the
    design-ID corner (upper left) of its design, cropped from the assembled cells."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    ids = [die['die_id'] for die in dies]
    if ids != [f'D{i:02d}' for i in range(len(dies))]:
        raise ValueError(f'Die IDs are not unique and sequential: {ids}')

    def crop(polygons, window):
        x0, y0, x1, y1 = window
        near = [p for p in polygons
                if p.bounding_box()[0][0] < x1 and p.bounding_box()[1][0] > x0
                and p.bounding_box()[0][1] < y1 and p.bounding_box()[1][1] > y0]
        if not near:
            return []
        return [p.points for p in gdstk.boolean(near, gdstk.rectangle((x0, y0), (x1, y1)), 'and', precision=.001)]

    design_corner = {}
    for name, info in sources.items():
        w, h = info['size_um']
        window = (-w / 2, h / 2 - 260, -w / 2 + 900, h / 2)
        design_corner[name] = (window, crop(cells[name].get_polygons(), window))
    cols = 6
    rows = -(-len(dies) // cols)
    fig, axes = plt.subplots(rows, cols * 2, figsize=(cols * 4.2, rows * 1.1), squeeze=False)
    for ax in axes.flat:
        ax.set_axis_off()
    for k, die in enumerate(dies):
        name = die['design']
        x0, y0, x1, y1 = wafer_die_ids['windows_um'][name]
        id_window = (x0 - 40, y0 - 40, x1 + 40, y1 + 40)
        shapes = crop(cells[die['cell']].polygons + cells[name].get_polygons(), id_window)
        panels = ((design_corner[name][0], design_corner[name][1], f'{name} design ID'),
                  (id_window, shapes, die['die_id']))
        for offset, (window, polys, title) in enumerate(panels):
            ax = axes[k // cols][2 * (k % cols) + offset]
            ax.set_axis_on(); ax.set_xticks([]); ax.set_yticks([])
            ax.add_collection(PolyCollection(polys, facecolor='#59636e', edgecolor='none'))
            ax.set_xlim(window[0], window[2]); ax.set_ylim(window[1], window[3]); ax.set_aspect('equal')
            ax.set_title(title, fontsize=7)
    fig.suptitle(f'{Path(path).name}: {len(dies)} dies, die ID height {DIE_ID_SIZE:g} um (dark = metal)',
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor='white')
    plt.close(fig)


def assemble(folder, *, replace=False, source_folder=None, design6_gds=None, output_stem='summery',
             inputs=None, devices_only=False):
    """inputs: {INPUTS key: GDS path} overrides (e.g. new Design masks); design6_gds is the older D6-only form."""
    folder.mkdir(parents=True, exist_ok=True)
    inputs = dict(inputs or {})
    if design6_gds:
        inputs.setdefault('D6', design6_gds)
    active_inputs = {k: INPUTS[k] for k in DEVICE_KEYS} if devices_only else INPUTS
    unknown = set(inputs) - set(active_inputs)
    if unknown:
        raise ValueError(f'Unknown input keys {sorted(unknown)}; use {sorted(active_inputs)}')
    outputs = [folder / f'{output_stem}.{suffix}' for suffix in ('gds', 'png', 'json')]
    outputs.append(folder / f'{output_stem}_die_id_check.png')
    if not replace and any(path.exists() for path in outputs):
        raise FileExistsError('Wafer output exists. Use --replace to refresh it.')
    cells, sources = {}, {}
    for name, filename in active_inputs.items():
        path = Path(inputs[name]) if name in inputs else (source_folder or folder) / filename
        cells[name], sources[name] = load_die(path, name)
        if path.parent.resolve() != folder.resolve():
            sources[name]['path_from_assembly'] = Path(relpath(path, folder)).as_posix()
    size = tuple(max(info['size_um'][axis] for info in sources.values()) for axis in (0, 1))
    dies, counts, marks = place_dies(size, devices_only=devices_only)
    wafer_die_ids = add_die_ids(cells, sources, dies)

    lib = gdstk.Library('WAFER', unit=1e-6, precision=1e-9)
    top = lib.new_cell('WAFER_4_INCH')
    lib.add(*cells.values())
    cross = lib.new_cell('ALIGNMENT_CROSS')
    cross.add(alignment_cross())
    top.add(*(gdstk.Reference(cells[d['cell']], origin=d['centre_um']) for d in dies))
    top.add(*(gdstk.Reference(cross, origin=origin) for origin in marks))
    top.add(*wafer_outline())

    width, height = size
    max_radius = max(hypot(abs(d['centre_um'][0]) + width / 2, abs(d['centre_um'][1]) + height / 2)
                     for d in dies)
    manifest = {
        'wafer_diameter_um': WAFER_RADIUS * 2,
        'wafer_guide_layer': WAFER_GUIDE,
        'die_slot_size_um': size,
        'pitch_um': (width + STREET, height + STREET),
        'street_width_um': STREET,
        'minimum_radial_margin_to_slot_corner_um': round(WAFER_RADIUS - max_radius, 3),
        'alignment_cross': {'layer': METAL, 'span_um': CROSS_SPAN, 'stroke_um': CROSS_STROKE,
                            'centres_um': marks},
        'sources': sources,
        'wafer_die_ids': wafer_die_ids,
        'counts': counts,
        'dies': dies,
    }
    with TemporaryDirectory(dir=folder) as temp:
        gds, png, record, id_png = [Path(temp) / path.name for path in outputs]
        lib.write_gds(str(gds), max_points=GDS_MAX_POINTS)
        verify_gds(gds, cells, dies, marks, size, sources, wafer_die_ids)
        check_die_ids(id_png, cells, dies, wafer_die_ids, sources)
        plot_wafer(png, size, dies, counts, marks)
        record.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        for source, target in zip((gds, png, record, id_png), outputs):
            source.replace(target)
    print(f'Wrote {outputs[0]} with {len(dies)} dies and {len(marks)} crosses.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path, help='output folder and default source folder')
    parser.add_argument('--source-folder', type=Path, help='read the standard inputs from a separate folder')
    parser.add_argument('--design6-gds', type=Path, help='use this Design6_TXT fabrication mask in the assembly')
    parser.add_argument('--replace', action='store_true')
    parser.add_argument('--output-stem', default='summery', help='stem of the output GDS, map and manifest')
    parser.add_argument('--devices-only', action='store_true',
                        help='populate all 52 slots with D2, D3, D5, D6 only (13 each, no Chip4/PCM)')
    parser.add_argument('--input', action='append', default=[], metavar='KEY=PATH',
                        help='replace one source mask, KEY in ' + ', '.join(INPUTS) + ' (repeatable)')
    args = parser.parse_args()
    assemble(args.folder, replace=args.replace, source_folder=args.source_folder, design6_gds=args.design6_gds,
             output_stem=args.output_stem, inputs=dict(item.split('=', 1) for item in args.input),
             devices_only=args.devices_only)
