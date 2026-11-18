"""Export the standalone junction/array mask draft, without chip ground."""
import json
import hashlib
import argparse
from pathlib import Path
import sys

import gdspy
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.affinity import scale
from qiskit_metal import designs

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [*[str(ROOT / 'pattern design' / g) for g in ('', 'fluxonium', 'transmon', 'resonator', 'marker')], str(ROOT / 'code/general')]
from fluxonium_junction import FluxoniumJunction
from gds_tools import export_component_polygons


def main():
    folder = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--num-cells', type=int, default=43, help='Cell count per arm')
    parser.add_argument('--junction-gap', default='0.075um', help='Shared arrow clearance, with units')
    parser.add_argument('--junction-width', default='0.25um', help='Right-hand junction bar y width, with units')
    parser.add_argument('--output', type=Path, default=folder / 'fluxonium_junction.gds')
    args = parser.parse_args()
    output = args.output.resolve()
    record_path = output.with_name(output.stem + '_record.json')
    if output.name == 'fluxonium_junction.gds':
        record_path = output.with_name('export_record.json')
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite {output}; choose a new output filename')
    if record_path.exists():
        raise FileExistsError(f'Refusing to overwrite {record_path}')
    design = designs.DesignPlanar()
    component = FluxoniumJunction(design, 'fluxonium_junction', options=dict(
        num_cells=args.num_cells, junction_gap=args.junction_gap,
        junction_width=args.junction_width))
    assert component.status == 'good', component.status
    expected = design.qgeometry.tables['poly'].copy()
    export_component_polygons(design, output, ['fluxonium_junction'])
    lib = gdspy.GdsLibrary(infile=str(output))
    assert abs(lib.unit-1e-6) < 1e-15
    assert abs(lib.precision-1e-9) < 1e-18
    polygons = lib.cells['TOP'].get_polygons(by_spec=True)
    assert set(polygons) == {(i, 0) for i in (1, 2, 3, 4)}, list(polygons)
    for layer in (1, 2, 3, 4):
        source = unary_union(list(expected[expected.layer == layer].geometry))
        source = scale(source, 1000, 1000, origin=(0, 0))
        exported = unary_union([Polygon(points) for points in polygons[layer, 0]])
        assert source.symmetric_difference(exported).area < 1e-7, layer
    cells = [Polygon(points) for points in polygons[3, 0]
             if abs(Polygon(points).area-0.998*2) < 1e-7]
    assert len(cells) == 2*args.num_cells, len(cells)
    for cell in cells:
        x0, y0, x1, y1 = cell.bounds
        assert abs(x1-x0-0.998) < 1e-9 and abs(y1-y0-2) < 1e-9
    gap_um = component.p.junction_gap / design.parse_value('1um')
    for metal, recess in ((0, 0), (1, 0), (1, 9), (2, 11), (2, 9), (2, 10)):
        measured = Polygon(polygons[1, 0][metal]).distance(Polygon(polygons[2, 0][recess]))
        assert abs(measured-gap_um) < 1e-9, (metal, recess, measured, gap_um)
    record = dict(status='draft; not fabrication or EM validated',
                  component='FluxoniumJunction', options=dict(component.options),
                  source_sha256=component.reference()['source_sha256'],
                  component_sha256=hashlib.sha256((ROOT / 'pattern design/fluxonium/fluxonium_junction.py').read_bytes()).hexdigest(),
                  gds_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                  unit_m=lib.unit, precision_m=lib.precision,
                  bounds_um=lib.cells['TOP'].get_bounding_box().tolist(),
                  polygons_per_layer={str(k): len(v) for k, v in polygons.items()},
                  checks=f'four layers; geometry round trip; {len(cells)} fixed-size cells; six shared gaps',
                  simulations_rerun=False)
    record_path.write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
    print(output)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
