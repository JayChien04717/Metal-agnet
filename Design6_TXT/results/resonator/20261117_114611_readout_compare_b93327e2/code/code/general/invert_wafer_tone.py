"""Write and verify wafer etch artwork as GDS and DXF."""

from argparse import ArgumentParser
from collections import Counter
import json
import os
from pathlib import Path
from shutil import which
import subprocess
from tempfile import TemporaryDirectory


PROCESS = (1, 0)
GUIDE = (99, 0)
GDS_EDGE_TOLERANCE_NM = 3


def _invert(source, manifest, output, dxf_output, *, dies_only):
    import pya

    info = json.loads(Path(manifest).read_text(encoding='utf-8'))

    def read(filename):
        layout = pya.Layout()
        layout.read(filename)
        if layout.dbu != .001:
            raise ValueError('Expected a 1 nm GDS database unit')
        return layout

    def nm(coords):
        return tuple(round(x * 1000) for x in coords)

    def placements(layout):
        result = Counter()
        for instance in layout.top_cell().each_inst():
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
    expected = Counter((die['design'], nm(die['centre_um'])) for die in info['dies'])
    expected.update(('ALIGNMENT_CROSS', nm(point)) for point in info['alignment_cross']['centres_um'])
    if placements(positive) != expected:
        raise ValueError('Source GDS placements differ from the assembly JSON')
    process = positive.layer(*PROCESS)
    if not top.shapes(process).is_empty():
        raise ValueError('Unexpected top-level process geometry outside the referenced cells')

    fields = {name: die_field(*die['size_um']) for name, die in info['sources'].items()}
    guide = region(positive, top_name, GUIDE)
    outline, = guide.each_merged()
    wafer = pya.Region(pya.Polygon(list(outline.each_point_hull())))
    die_fields = pya.Region()
    for die in info['dies']:
        die_fields += fields[die['design']].moved(*nm(die['centre_um']))
    total_area = sum(fields[die['design']].area() for die in info['dies'])
    if not (die_fields - wafer).is_empty() or die_fields.merged().area() != total_area:
        raise ValueError('Inversion boundaries must be disjoint and inside the wafer')
    outside_dies = pya.Region()
    if dies_only:
        top.shapes(positive.layer(*GUIDE)).clear()
        guide = pya.Region()
    else:
        cross = region(positive, 'ALIGNMENT_CROSS', PROCESS)
        fiducials = pya.Region()
        for centre in info['alignment_cross']['centres_um']:
            fiducials += cross.moved(*nm(centre))
        if not (fiducials - wafer).is_empty() or not (fiducials & die_fields).is_empty():
            raise ValueError('Alignment crosses must lie on the wafer outside the dies')
        outside_dies = wafer - (die_fields | fiducials)
        top.shapes(process).insert(outside_dies)

    for name, field in fields.items():
        cell = positive.cell(name)
        if any(cell.each_inst()):
            raise ValueError(f'{name}: expected a flattened die cell')
        metal = region(positive, name, PROCESS)
        # The assembly dimensions include the vendor marks removed before placement.
        if metal.is_empty() or not (metal - field).is_empty():
            raise ValueError(f'{name}: expected source metal inside the die inversion boundary')
        etch = field - metal
        cell.shapes(process).clear()
        cell.shapes(process).insert(etch)
    # The cross locations lie outside the dies and are omitted from exposure.
    positive.delete_cell(positive.cell('ALIGNMENT_CROSS').cell_index())
    active_placements = Counter({site: count for site, count in expected.items() if site[0] in fields})

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
    if dies_only and not wafer_etch.is_empty():
        raise ValueError('The die-only GDS contains exposure outside the die references')
    verify_complement(top_name, outside_dies, wafer_etch)
    for name, field in fields.items():
        die_etch = region(etched, name, PROCESS)
        if dies_only and not (die_etch - field).is_empty():
            raise ValueError(f'{name}: exposure extends beyond the die boundary')
        verify_complement(name, field - region(positive, name, PROCESS), die_etch)
    if not (guide ^ region(etched, top_name, GUIDE)).is_empty():
        raise ValueError('Wafer guide geometry changed')
    del positive

    # The DXF writer stores closed POLYLINE contours in micrometres and reuses die blocks.
    etched.write(dxf_output)
    drawing = read(dxf_output)
    if placements(drawing) != active_placements:
        raise ValueError('DXF die placement changed')
    drawing_top = drawing.top_cell().name
    for gds_name, dxf_name, layer in [(top_name, drawing_top, PROCESS), (top_name, drawing_top, GUIDE),
                                      *((name, name, PROCESS) for name in fields)]:
        if not (region(etched, gds_name, layer) ^ region(drawing, dxf_name, layer)).is_empty():
            raise ValueError(f'{dxf_name}: DXF geometry differs from the verified GDS on layer {layer}')

    return {'dies': len(info['dies']), 'die_types': len(info['sources']),
            'rounding_area_nm2': rounding_area}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path, help='Folder with summery.gds and summery.json')
    parser.add_argument('--dies-only', action='store_true', help='Expose only inside the die rectangles')
    parser.add_argument('--replace', action='store_true', help='Replace an earlier etch export after verification')
    args = parser.parse_args()
    folder = args.folder.resolve()
    source, manifest = (folder / name for name in ('summery.gds', 'summery.json'))
    stem = 'summery_etch_dies_only' if args.dies_only else 'summery_etch'
    targets = [folder / f'{stem}.{extension}' for extension in ('gds', 'dxf')]
    for target in targets:
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
                               etch_dxf=staged_dxf, proof_file=proof, die_scope=int(args.dies_only)).items():
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
    print(f'Verified {stem} GDS and DXF in {folder}: {counts["die_types"]} die types, '
          f'{counts["dies"]} placements.')
    if counts['rounding_area_nm2']:
        print(f'Rounding within {GDS_EDGE_TOLERANCE_NM} nm of the complement boundary '
              f'(nm² per cell): {counts["rounding_area_nm2"]}')


if __name__ == '__main__':
    if 'positive_gds' in globals():
        proof = _invert(positive_gds, assembly_json, etch_gds, etch_dxf, dies_only=die_scope == '1')
        Path(proof_file).write_text(json.dumps(proof), encoding='utf-8')
    else:
        main()
