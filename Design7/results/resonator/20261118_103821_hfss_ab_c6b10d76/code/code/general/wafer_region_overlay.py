"""Replace bounded positive-metal regions in a parent wafer GDS using KLayout."""
import json
import os
from pathlib import Path
from shutil import which
import subprocess
from tempfile import TemporaryDirectory


def overlay_wafer_regions(parent, updated, boxes_nm):
    """Apply selected regions of an exported GDS to a parent on the same 1 nm grid."""
    updated = Path(updated).resolve()
    exe = os.environ.get('KLAYOUT_EXE') or which('klayout_app') or which('klayout')
    if not exe and os.name == 'nt':
        exe = Path.home() / 'AppData/Roaming/KLayout/klayout_app.exe'
    if not exe or not Path(exe).is_file():
        raise FileNotFoundError('KLayout executable not found. Set KLAYOUT_EXE to its path.')
    with TemporaryDirectory(dir=updated.parent) as temp:
        folder = Path(temp)
        spec, staged = folder / 'boxes.json', folder / updated.name
        spec.write_text(json.dumps(list(boxes_nm)), encoding='utf-8')
        command = [str(exe), '-b', '-r', str(Path(__file__).resolve()),
                   '-rd', f'parent_file={Path(parent).resolve()}', '-rd', f'updated_file={updated}',
                   '-rd', f'boxes_file={spec}', '-rd', f'output_file={staged}']
        result = subprocess.run(command, capture_output=True, text=True, errors='replace')
        if result.returncode or not staged.is_file():
            raise RuntimeError(f'KLayout regional overlay failed: {result.stderr or result.stdout}')
        staged.replace(updated)


def _overlay(pya, parent, updated, boxes, output):
    original, replacement = pya.Layout(), pya.Layout()
    original.read(parent)
    replacement.read(updated)
    if original.dbu != .001 or replacement.dbu != .001:
        raise ValueError('Expected 1 nm database units in both GDS files')
    top, new_top = original.top_cell(), replacement.top_cell()
    if top.bbox() != new_top.bbox():
        raise ValueError('Parent and replacement die bounds must agree')
    metal, new_metal = original.find_layer(1, 0), replacement.find_layer(1, 0)
    region = pya.Region()
    for box in boxes:
        region.insert(pya.Box(*box))
    region.merge()
    new_shapes = pya.Region(new_top.begin_shapes_rec(new_metal)) & region

    # Keep unaffected fabricated ground tiles as separate cells. The small
    # simulation geometry can be flattened without touching junction annotations.
    instances = list(top.each_inst())
    tiles = [original.cell(inst.cell_index) for inst in instances]
    simple = (len({cell.cell_index() for cell in tiles}) == len(tiles)
              and all(inst.trans == pya.Trans() and not inst.is_regular_array() for inst in instances)
              and all(not list(cell.each_inst()) for cell in tiles))
    if not simple:
        top.flatten(True)
        tiles = []
    for cell in [top, *tiles]:
        if not cell.bbox().overlaps(region.bbox()):
            continue
        shapes = cell.shapes(metal)
        data = pya.Region(shapes)
        if not (data & region).is_empty():
            remaining = data - region
            shapes.clear()
            shapes.insert(remaining)
    top.shapes(metal).insert(new_shapes)

    options = pya.SaveLayoutOptions()
    options.gds2_max_vertex_count = 8190
    original.write(output, options)


if __name__ == '__main__':
    import pya

    _overlay(pya, parent_file, updated_file, json.loads(Path(boxes_file).read_text(encoding='utf-8')), output_file)
