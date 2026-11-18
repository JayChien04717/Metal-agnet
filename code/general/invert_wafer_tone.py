"""Write and verify wafer etch artwork as GDS and DXF."""

from argparse import ArgumentParser
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
from shutil import which
import subprocess
from tempfile import TemporaryDirectory


PROCESS = (1, 0)
GUIDE = (99, 0)
CROSS_CELL = 'ALIGNMENT_CROSS'
GDS_EDGE_TOLERANCE_NM = 3


def _invert(source, manifest, output, dxf_output, *, scope):
    import pya

    info = json.loads(Path(manifest).read_text(encoding='utf-8'))
    expose_wafer = scope == 'wafer'
    expose_crosses = scope == 'dies-and-crosses'
    expose_streets = scope == 'dies-and-streets'
    expose_outside = expose_wafer or expose_streets

    def read(filename):
        layout = pya.Layout()
        layout.read(filename)
        if layout.dbu != .001:
            raise ValueError('Expected a 1 nm GDS database unit')
        return layout

    def nm(coords):
        return tuple(round(x * 1000) for x in coords)

    def placements(layout, name=None):
        result = Counter()
        cell = layout.cell(name) if name else layout.top_cell()
        for instance in cell.each_inst():
            trans = instance.trans
            if trans.angle != 0 or trans.is_mirror() or instance.is_regular_array():
                raise ValueError('Unexpected rotated, reflected or arrayed wafer placement')
            result[instance.cell.name, (trans.disp.x, trans.disp.y)] += 1
        return result

    def region(layout, name, layer):
        return pya.Region(layout.cell(name).shapes(layout.layer(*layer))).merged()

    def die_field(width, height):
        x, y = nm((width / 2, height / 2))
        return pya.Region(pya.Box(-x, -y, x, y))

    positive = read(source)
    top = positive.top_cell()
    top_name = top.name
    die_placements = Counter((die.get('cell', die['design']), nm(die['centre_um'])) for die in info['dies'])
    cross_points = info['alignment_cross']['centres_um']
    expected = die_placements + Counter((CROSS_CELL, nm(point)) for point in cross_points)
    if placements(positive) != expected:
        raise ValueError('Source GDS placements differ from the assembly JSON')
    process = positive.layer(*PROCESS)
    if not top.shapes(process).is_empty():
        raise ValueError('Unexpected top-level process geometry outside the referenced cells')

    fields = {name: die_field(*die['size_um']) for name, die in info['sources'].items()}
    label_settings = info.get('wafer_die_ids')
    label_windows, label_fields, label_children = {}, {}, {}
    if label_settings:
        for name, bounds in label_settings['windows_um'].items():
            label_windows[name] = pya.Region(pya.Box(*nm(bounds)))
            if not (label_windows[name] - fields[name]).is_empty():
                raise ValueError(f'{name}: die ID extends beyond the die')
        ordered = sorted(info['dies'], key=lambda die: (-die['centre_um'][1], die['centre_um'][0]))
        if [die['die_id'] for die in ordered] != [f'D{i:02d}' for i in range(len(ordered))]:
            raise ValueError('Wafer die IDs do not follow the selected row order')
        for die in ordered:
            name = die['cell']
            label_fields[name] = label_windows[die['design']]
            label_children[name] = Counter({(die['design'], (0, 0)): 1})
            if placements(positive, name) != label_children[name]:
                raise ValueError(f'{name}: unexpected source-die references')
    guide = region(positive, top_name, GUIDE)
    outline, = guide.each_merged()
    wafer = pya.Region(pya.Polygon(list(outline.each_point_hull())))
    die_fields = pya.Region()
    for die in info['dies']:
        die_fields += fields[die['design']].moved(*nm(die['centre_um']))
    total_area = sum(fields[die['design']].area() for die in info['dies'])
    if not (die_fields - wafer).is_empty() or die_fields.merged().area() != total_area:
        raise ValueError('Inversion boundaries must be disjoint and inside the wafer')
    cross = region(positive, CROSS_CELL, PROCESS)
    if cross.is_empty() or any(positive.cell(CROSS_CELL).each_inst()):
        raise ValueError('Expected a flattened alignment cross on the process layer')
    fiducials = pya.Region()
    for centre in cross_points:
        fiducials += cross.moved(*nm(centre))
    if not (fiducials - wafer).is_empty() or not (fiducials & die_fields).is_empty():
        raise ValueError('Alignment crosses must lie on the wafer outside the dies')

    street_field = pya.Region()
    street_x, street_y = [], []
    if expose_streets:
        street_width, = nm((info['street_width_um'],))
        pitch_x, pitch_y = nm(info['pitch_um'])
        slot_width, slot_height = nm(info['die_slot_size_um'])
        if pitch_x - slot_width != street_width or pitch_y - slot_height != street_width:
            raise ValueError('Die pitch and slot dimensions disagree with the street width')
        if street_width != 100_000:
            raise ValueError('This export requires the selected 100 µm street width')
        street_x = sorted({nm(die['centre_um'])[0] + sign * pitch_x // 2
                           for die in info['dies'] for sign in (-1, 1)})
        street_y = sorted({nm(die['centre_um'])[1] + sign * pitch_y // 2
                           for die in info['dies'] for sign in (-1, 1)})
        half = street_width // 2
        extent = wafer.bbox()
        for x in street_x:
            street_field.insert(pya.Box(x - half, extent.bottom, x + half, extent.top))
        for y in street_y:
            street_field.insert(pya.Box(extent.left, y - half, extent.right, y + half))
        street_field = street_field.merged() & wafer
        if not (street_field & die_fields).is_empty():
            raise ValueError('A dicing street intersects a chip')
        if not (fiducials - street_field).is_empty():
            raise ValueError('A metal cross extends beyond the dicing streets')

    outside_dies = pya.Region()
    if expose_wafer:
        outside_dies = wafer - (die_fields | fiducials)
    elif expose_streets:
        outside_dies = street_field - fiducials
    if expose_outside:
        top.shapes(process).insert(outside_dies)
    else:
        top.shapes(positive.layer(*GUIDE)).clear()
        guide = pya.Region()

    for name, field in fields.items():
        cell = positive.cell(name)
        if any(cell.each_inst()):
            raise ValueError(f'{name}: expected a flattened die cell')
        metal = region(positive, name, PROCESS)
        # The assembly dimensions include the vendor marks removed before placement.
        if metal.is_empty() or not (metal - field).is_empty():
            raise ValueError(f'{name}: expected source metal inside the die inversion boundary')
        etch = field - metal - label_windows.get(name, pya.Region())
        cell.shapes(process).clear()
        cell.shapes(process).insert(etch)
    for name, window in label_fields.items():
        cell = positive.cell(name)
        metal = region(positive, name, PROCESS)
        if not (metal - window).is_empty():
            raise ValueError(f'{name}: local label patch extends outside its window')
        cell.shapes(process).clear()
        cell.shapes(process).insert(window - metal)
    if not expose_crosses:
        positive.delete_cell(positive.cell(CROSS_CELL).cell_index())
    active_placements = expected if expose_crosses else die_placements

    options = pya.SaveLayoutOptions()
    options.gds2_max_vertex_count = 8190
    positive.write(output, options)
    del positive

    positive, etched = read(source), read(output)
    if placements(etched) != active_placements:
        raise ValueError('Etch GDS placements differ from the expected assembly')
    rounding_area = {}

    def verify_complement(name, complement, artwork):
        delta = artwork ^ complement
        if not delta.is_empty():
            # Fracturing a polygon for GDS can round a cut on a sloped edge.
            tolerance = GDS_EDGE_TOLERANCE_NM
            window = pya.Region(delta.bbox().enlarged(2 * tolerance, 2 * tolerance))
            local = complement & window
            if not (delta - local.sized(tolerance)).is_empty() or not (delta & local.sized(-tolerance)).is_empty():
                raise ValueError(f'{name}: complement differs beyond {tolerance} nm of its boundary')
            rounding_area[name] = delta.area()

    wafer_etch = region(etched, top_name, PROCESS)
    if not expose_outside and not wafer_etch.is_empty():
        raise ValueError('Found top-level exposure outside the allowed references')
    verify_complement(top_name, outside_dies, wafer_etch)
    for name, field in fields.items():
        die_etch = region(etched, name, PROCESS)
        if not expose_outside and not (die_etch - field).is_empty():
            raise ValueError(f'{name}: exposure extends beyond the die boundary')
        complement = field - region(positive, name, PROCESS) - label_windows.get(name, pya.Region())
        verify_complement(name, complement, die_etch)
    for name, window in label_fields.items():
        if placements(etched, name) != label_children[name]:
            raise ValueError(f'{name}: source-die reference changed in the etch GDS')
        verify_complement(name, window - region(positive, name, PROCESS), region(etched, name, PROCESS))
    process_cells = [*fields, *label_fields]
    if expose_crosses:
        if not (cross ^ region(etched, CROSS_CELL, PROCESS)).is_empty():
            raise ValueError('Alignment cross contours changed in the etch GDS')
        process_cells.append(CROSS_CELL)
    if not (guide ^ region(etched, top_name, GUIDE)).is_empty():
        raise ValueError('Wafer guide geometry changed')
    if expose_streets:
        if not (wafer_etch & fiducials).is_empty():
            raise ValueError('Etch exposure removes a retained metal cross')
        retained = (street_field - wafer_etch).merged()
        if not (retained ^ fiducials).is_empty() or retained.count() != len(cross_points):
            raise ValueError('Unexpected retained metal in the dicing streets')
    del positive

    # The DXF writer stores closed POLYLINE contours in micrometres and reuses die blocks.
    etched.write(dxf_output)
    drawing = read(dxf_output)
    if placements(drawing) != active_placements:
        raise ValueError('DXF die placement changed')
    for name, expected_children in label_children.items():
        if placements(drawing, name) != expected_children:
            raise ValueError(f'{name}: DXF source-die reference changed')
    drawing_top = drawing.top_cell().name
    for gds_name, dxf_name, layer in [(top_name, drawing_top, PROCESS), (top_name, drawing_top, GUIDE),
                                      *((name, name, PROCESS) for name in process_cells)]:
        if not (region(etched, gds_name, layer) ^ region(drawing, dxf_name, layer)).is_empty():
            raise ValueError(f'{dxf_name}: DXF geometry differs from the verified GDS on layer {layer}')

    return {'source_gds': Path(source).name,
            'source_gds_sha256': hashlib.sha256(Path(source).read_bytes()).hexdigest(),
            'source_manifest_sha256': hashlib.sha256(Path(manifest).read_bytes()).hexdigest(),
            'scope': scope, 'unit': 'µm', 'database_grid_nm': 1,
            'process_layer': PROCESS, 'guide_layer': GUIDE if expose_outside else None,
            'street_width_um': info['street_width_um'] if expose_streets else None,
            'street_centrelines_x_um': [x / 1000 for x in street_x],
            'street_centrelines_y_um': [y / 1000 for y in street_y],
            'dies': len(info['dies']), 'die_types': len(info['sources']),
            'wafer_die_ids': {'count': len(label_fields), 'first': info['dies'][0]['die_id'],
                             'last': info['dies'][-1]['die_id']} if label_fields else None,
            'retained_metal_crosses': len(cross_points) if expose_outside else 0,
            'retained_cross_span_um': info['alignment_cross']['span_um'] if expose_outside else None,
            'retained_cross_stroke_um': info['alignment_cross']['stroke_um'] if expose_outside else None,
            'dxf_matches_gds': True, 'em_rerun': False,
            'exposed_crosses': len(cross_points) if expose_crosses else 0,
            'rounding_area_nm2': rounding_area}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path, help='Folder with summery.gds and summery.json')
    parser.add_argument('--source-stem', default='summery', help='stem of the source wafer GDS and JSON')
    scope_options = parser.add_mutually_exclusive_group()
    scope_options.add_argument('--dies-only', dest='scope', action='store_const', const='dies-only',
                               help='Expose only inside the die rectangles')
    scope_options.add_argument('--dies-and-crosses', dest='scope', action='store_const', const='dies-and-crosses',
                               help='Expose inside the die rectangles and on the alignment cross strokes')
    scope_options.add_argument('--dies-and-streets', dest='scope', action='store_const', const='dies-and-streets',
                               help='Etch chip openings and 100 µm streets to the wafer edge; retain metal crosses')
    scope_options.add_argument('--wafer', dest='scope', action='store_const', const='wafer',
                               help='Reproduce the whole-wafer complement, including the blank perimeter')
    parser.set_defaults(scope='dies-and-streets')
    parser.add_argument('--replace', action='store_true', help='Replace an earlier etch export after verification')
    args = parser.parse_args()
    folder = args.folder.resolve()
    source, manifest = (folder / f'{args.source_stem}.{suffix}' for suffix in ('gds', 'json'))
    stem = args.source_stem + '_etch'
    if args.scope != 'wafer':
        stem += '_' + args.scope.replace('-', '_')
    targets = [folder / f'{stem}.{extension}' for extension in ('gds', 'dxf')]
    verification = folder / f'{stem}_verification.json' if args.scope == 'dies-and-streets' else None
    for target in [*targets, *([verification] if verification else [])]:
        if target.exists() and not args.replace:
            raise FileExistsError(target)

    klayout = os.environ.get('KLAYOUT_EXE') or which('klayout_app') or which('klayout')
    if not klayout and os.name == 'nt':
        klayout = Path.home() / 'AppData/Roaming/KLayout/klayout_app.exe'
    if not klayout or not Path(klayout).is_file():
        raise FileNotFoundError('Set KLAYOUT_EXE to the KLayout executable path')

    # Windows may briefly hold an open directory handle after the DXF readback.
    with TemporaryDirectory(dir=folder, ignore_cleanup_errors=True) as tmp:
        staged_gds, staged_dxf = (Path(tmp) / target.name for target in targets)
        proof = Path(tmp) / 'verification.json'
        command = [str(klayout), '-b', '-r', str(Path(__file__).resolve())]
        for key, value in dict(positive_gds=source, assembly_json=manifest, etch_gds=staged_gds,
                               etch_dxf=staged_dxf, proof_file=proof, exposure_scope=args.scope).items():
            command.extend(('-rd', f'{key}={value}'))
        result = subprocess.run(command, capture_output=True, text=True, errors='replace')
        if result.returncode or not all(path.is_file() for path in (proof, staged_gds, staged_dxf)):
            raise RuntimeError(f'KLayout inversion failed: {result.stderr or result.stdout}')
        counts = json.loads(proof.read_text(encoding='utf-8'))
        for staged, target in zip((staged_gds, staged_dxf), targets):
            if args.replace:
                staged.replace(target)
            else:
                staged.rename(target)
        if verification:
            proof.replace(verification)
    print(f'Verified {stem} GDS and DXF in {folder}: {counts["die_types"]} die types, '
          f'{counts["dies"]} die placements and {counts["exposed_crosses"]} exposed alignment crosses.')
    if counts['rounding_area_nm2']:
        print(f'Rounding within {GDS_EDGE_TOLERANCE_NM} nm of the complement boundary '
              f'(nm^2 per cell): {counts["rounding_area_nm2"]}')


if __name__ == '__main__':
    if 'positive_gds' in globals():
        proof = _invert(positive_gds, assembly_json, etch_gds, etch_dxf, scope=exposure_scope)
        Path(proof_file).write_text(json.dumps(proof), encoding='utf-8')
    else:
        main()
