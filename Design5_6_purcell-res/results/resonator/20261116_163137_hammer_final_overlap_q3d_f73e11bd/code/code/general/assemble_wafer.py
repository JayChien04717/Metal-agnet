"""Assemble the six 20261116 fabrication layouts on one circular 4-inch wafer."""

import argparse
from collections import Counter
import hashlib
import json
from math import hypot
from pathlib import Path
from tempfile import TemporaryDirectory

import gdstk
import numpy as np


INPUTS = {
    'D2': 'Design2_2026-11-16_wafer_r4_fab.gds',
    'D3_3': 'Design3_3_2026-11-16_wafer_r4_fab.gds',
    'D5': 'Design5_2026-11-16_wafer_r4_fab.gds',
    'D6': 'Design6_2026-11-16_wafer_r4_fab.gds',
    'Chip4': 'AllspiceR3_v5p4_Chip4.oas',
    'PCM': 'AllspiceR3_v5p4_PCMdie_M1.oas',
}
COLORS = dict(zip(INPUTS, ('#80bfd0', '#92b4eb', '#ecc586', '#aacbaf', '#cfb1d7', '#e4aaa6')))

# Rows are north to south. Regions A–E and the six designs are balanced even
# across diagonal neighbors. A dash is an empty position at the circular edge.
ROWS = (
    '-     -     D2    D3_3  D6    Chip4 -     -',
    '-     D5    PCM   Chip4 PCM   D5    D2    -',
    'D6    D3_3  D2    D5    D3_3  D6    Chip4 PCM',
    'D2    D5    D6    PCM   D2    D5    D3_3  D6',
    'Chip4 PCM   D3_3  Chip4 D6    PCM   D2    PCM',
    'D5    D2    D5    D2    D3_3  Chip4 D6    D3_3',
    '-     D6    D3_3  Chip4 D6    D5    PCM   -',
    '-     -     D2    D5    D3_3  Chip4 -     -',
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


def place_dies(size):
    width, height = size
    pitch_x, pitch_y = width + STREET, height + STREET
    occupied = {}
    for j, row in enumerate(ROWS):
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
        if set(count) != set(INPUTS) or max(count.values()) - min(count.values()) > 1:
            raise ValueError(f'Unbalanced assignment in {region}: {count}')
    corners = Counter()
    for (i, j), die in occupied.items():
        for di, dj in ((1, 0), (0, 1), (1, 1), (1, -1)):
            if occupied.get((i + di, j + dj), {}).get('design') == die['design']:
                raise ValueError(f'Adjacent copies of {die["design"]} at {(i, j)}')
        corners.update((i + di, j + dj) for di in (0, 1) for dj in (0, 1))
    marks = sorted((i * pitch_x, j * pitch_y) for (i, j), count in corners.items() if count >= 2)
    return list(occupied.values()), counts, marks


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


def verify_gds(path, cells, dies, marks, size):
    result = gdstk.read_gds(str(path))
    top, = result.top_level()
    assert (result.unit, result.precision) == (1e-6, 1e-9)
    assert result.layers_and_datatypes() == {METAL, WAFER_GUIDE}
    assert top.bounding_box() == ((-WAFER_RADIUS, -WAFER_RADIUS), (WAFER_RADIUS, WAFER_RADIUS))
    expected = Counter((d['design'], tuple(d['centre_um'])) for d in dies)
    expected.update(('ALIGNMENT_CROSS', origin) for origin in marks)
    assert Counter((r.cell_name, tuple(r.origin)) for r in top.references) == expected
    assert all(r.rotation == 0 and r.magnification == 1 and not r.x_reflection for r in top.references)
    for cell in result.cells:
        if cell.name in cells:
            assert polygon_signature(cell) == polygon_signature(cells[cell.name]), cell.name
            low, high = np.asarray(cell.bounding_box())
            assert np.all(low >= -np.asarray(size) / 2 - .001)
            assert np.all(high <= np.asarray(size) / 2 + .001)
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
        ax.text(x, y, name, ha='center', va='center', fontsize=8, color='#152436')
    for label, origin in zip('ABCDE', ((0, 0), (-33.5, 32.5), (33.5, 32.5),
                                       (-33.5, -32.5), (33.5, -32.5))):
        ax.annotate(label, origin, ha='center', va='center', weight='bold', fontsize=10,
                    bbox={'boxstyle': 'circle,pad=.26', 'fc': 'white', 'ec': '#526274', 'lw': .6})
    ax.set(xlim=(-53, 53), ylim=(-53, 53), aspect='equal', xlabel='x (mm)', ylabel='y (mm)')
    ax.set_title(f'4-inch wafer  |  Ø101.6 mm  |  {len(dies)} dies', loc='left', weight='bold')
    legend.set_axis_off()
    legend.text(0, .97, 'Dies per region', weight='bold', fontsize=13)
    names = list(INPUTS)
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


def assemble(folder, *, replace=False):
    outputs = [folder / f'summery.{suffix}' for suffix in ('gds', 'png', 'json')]
    if not replace and any(path.exists() for path in outputs):
        raise FileExistsError('Wafer output exists. Use --replace to refresh it.')
    cells, sources = {}, {}
    for name, filename in INPUTS.items():
        cells[name], sources[name] = load_die(folder / filename, name)
    size = tuple(max(info['size_um'][axis] for info in sources.values()) for axis in (0, 1))
    dies, counts, marks = place_dies(size)

    lib = gdstk.Library('WAFER', unit=1e-6, precision=1e-9)
    top = lib.new_cell('WAFER_4_INCH')
    lib.add(*cells.values())
    cross = lib.new_cell('ALIGNMENT_CROSS')
    cross.add(alignment_cross())
    top.add(*(gdstk.Reference(cells[d['design']], origin=d['centre_um']) for d in dies))
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
        'counts': counts,
        'dies': dies,
    }
    with TemporaryDirectory(dir=folder) as temp:
        gds, png, record = [Path(temp) / path.name for path in outputs]
        lib.write_gds(str(gds), max_points=GDS_MAX_POINTS)
        verify_gds(gds, cells, dies, marks, size)
        plot_wafer(png, size, dies, counts, marks)
        record.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        for source, target in zip((gds, png, record), outputs):
            source.replace(target)
    print(f'Wrote {outputs[0]} with {len(dies)} dies and {len(marks)} crosses.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path, help='folder containing the six source masks')
    parser.add_argument('--replace', action='store_true')
    args = parser.parse_args()
    assemble(args.folder, replace=args.replace)
