"""Write and verify the die-local etch tone of an assembled wafer GDS."""

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


def _invert(source, manifest, output):
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
        return pya.Region(layout.cell(name).begin_shapes_rec(layout.layer(*layer))).merged()

    def die_field(width, height):
        x, y = nm((width / 2, height / 2))
        return pya.Region(pya.Box(-x, -y, x, y))

    positive = read(source)
    top = positive.top_cell()
    expected = Counter((die['design'], nm(die['centre_um'])) for die in info['dies'])
    expected.update(('ALIGNMENT_CROSS', nm(point)) for point in info['alignment_cross']['centres_um'])
    if placements(positive) != expected:
        raise ValueError('Source GDS placements differ from the assembly JSON')
    process = positive.layer(*PROCESS)
    if not top.shapes(process).is_empty():
        raise ValueError('Unexpected top-level process geometry outside the referenced cells')

    for name, die in info['sources'].items():
        cell = positive.cell(name)
        if any(cell.each_inst()):
            raise ValueError(f'{name}: expected a flattened die cell')
        metal = region(positive, name, PROCESS)
        field = die_field(*die['size_um'])
        if metal.bbox() != field.bbox() or not (metal - field).is_empty():
            raise ValueError(f'{name}: source metal disagrees with the die inversion boundary')
        etch = field - metal
        cell.shapes(process).clear()
        cell.shapes(process).insert(etch)

    options = pya.SaveLayoutOptions()
    options.gds2_max_vertex_count = 8190
    positive.write(output, options)
    del positive

    positive, etched = read(source), read(output)
    if placements(etched) != expected or not etched.top_cell().shapes(etched.layer(*PROCESS)).is_empty():
        raise ValueError('Etch GDS placements or street geometry differ from the expected assembly')
    for name, die in info['sources'].items():
        metal = region(positive, name, PROCESS)
        etch = region(etched, name, PROCESS)
        if not (etch ^ (die_field(*die['size_um']) - metal)).is_empty():
            raise ValueError(f'{name}: the written etch shapes are not the exact complement of the metal')
    for name, layer in (('ALIGNMENT_CROSS', PROCESS), (top.name, GUIDE)):
        if not (region(positive, name, layer) ^ region(etched, name, layer)).is_empty():
            raise ValueError(f'{name}: guide geometry changed')

    return {'dies': len(info['dies']), 'die_types': len(info['sources']),
            'crosses': len(info['alignment_cross']['centres_um'])}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path, help='Folder with summery.gds and summery.json')
    args = parser.parse_args()
    folder = args.folder.resolve()
    source, manifest, target = (folder / name for name in ('summery.gds', 'summery.json', 'summery_etch.gds'))
    if target.exists():
        raise FileExistsError(target)

    klayout = os.environ.get('KLAYOUT_EXE') or which('klayout_app') or which('klayout')
    if not klayout and os.name == 'nt':
        klayout = Path.home() / 'AppData/Roaming/KLayout/klayout_app.exe'
    if not klayout or not Path(klayout).is_file():
        raise FileNotFoundError('Set KLAYOUT_EXE to the KLayout executable path')

    with TemporaryDirectory(dir=folder) as tmp:
        staged, proof = Path(tmp) / target.name, Path(tmp) / 'verification.json'
        command = [str(klayout), '-b', '-r', str(Path(__file__).resolve())]
        for key, value in dict(positive_gds=source, assembly_json=manifest, etch_gds=staged, proof_file=proof).items():
            command.extend(('-rd', f'{key}={value}'))
        result = subprocess.run(command, capture_output=True, text=True, errors='replace')
        if result.returncode or not proof.is_file() or not staged.is_file():
            raise RuntimeError(f'KLayout inversion failed: {result.stderr or result.stdout}')
        counts = json.loads(proof.read_text(encoding='utf-8'))
        staged.rename(target)
    print(f'{target}: verified {counts["die_types"]} inverted die types, '
          f'{counts["dies"]} placements and {counts["crosses"]} alignment crosses.')


if __name__ == '__main__':
    if 'positive_gds' in globals():
        Path(proof_file).write_text(json.dumps(_invert(positive_gds, assembly_json, etch_gds)), encoding='utf-8')
    else:
        main()
