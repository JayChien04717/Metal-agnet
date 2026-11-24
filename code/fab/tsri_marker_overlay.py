# -*- coding: utf-8 -*-
r"""TSRI A12 e-beam marks on a 4-inch wafer assembly: references, site checks, and TSRI e-beam GDS.

Shared delivery tool per ``skill/fabrication_output.md`` §2 (TSRI A12 marks).
Uses ``pattern design/marker/tsri_a12_mark.py`` (`a12_groups()`):
- 8 chip-mark groups per die on layer 2/0 (4 fine-cross groups + 4 coarse crosses within +-4.8 mm)
- Optional wafer-level complete A12 fields on layer 3/0 (`--wafer-field X,Y`)
- Bare-rectangle opening check (group bbox + 50 um margin per side) on every die
- Writes `--deliver-gds` with layer 1/0 base metal vertex-identical to the positive wafer GDS

Usage:
    python -I -X utf8 -B code/general/tsri_marker_overlay.py --write-references tmp/
    python -I -X utf8 -B code/general/tsri_marker_overlay.py <delivery_dir> \
        --assembly-dir <tmp_dir> --stem <positive_stem> --out <report_dir> \
        --deliver-gds "<delivery_dir>/TSRI ebeam gds/<date>_<source>_ebeam_TSRI.gds"
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
for p in (ROOT / 'code' / 'general', ROOT / 'pattern design' / 'marker'):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import gdstk
import shapely
from shapely.geometry import Polygon, box

from tsri_a12_mark import a12_groups

CHIP_LAYER = 2
WAFER_LAYER = 3
BASE_LAYER = 1
MARGIN_UM = 50.0
AREA_TOL_UM2 = 0.1


def _to_gdstk_polys(shapely_polys, layer, datatype=0):
    out = []
    for g in shapely_polys:
        for p in shapely.get_parts(g):
            if p.is_empty or p.area <= 0:
                continue
            out.append(gdstk.Polygon(list(p.exterior.coords)[:-1], layer=layer, datatype=datatype))
    return out


def _a12_squares_um(margin_um=MARGIN_UM):
    out = {}
    for name, gs in a12_groups().items():
        x0, y0, x1, y1 = shapely.union_all(gs).bounds
        out[name] = box(x0 - margin_um, y0 - margin_um, x1 + margin_um, y1 + margin_um)
    return out


def write_references(out_dir, wafer_fields=()):
    """Write ``A12_chip_single_die.gds`` (2/0) and optional ``A12_wafer_fields.gds`` (3/0)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    chip_polys = [g for gs in a12_groups().values() for g in gs]

    chip_lib = gdstk.Library(unit=1e-6, precision=1e-9)
    chip_top = chip_lib.new_cell('A12_CHIP')
    chip_top.add(*_to_gdstk_polys(chip_polys, layer=CHIP_LAYER))
    chip_path = out_dir / 'A12_chip_single_die.gds'
    chip_lib.write_gds(str(chip_path))

    written = {'chip_reference': str(chip_path)}
    if wafer_fields:
        full_polys = [g for gs in a12_groups(include_corners=True).values() for g in gs]
        w_lib = gdstk.Library(unit=1e-6, precision=1e-9)
        w_top = w_lib.new_cell('A12_WAFER_FIELDS')
        field_cell = w_lib.new_cell('A12_FIELD')
        field_cell.add(*_to_gdstk_polys(full_polys, layer=WAFER_LAYER))
        for x_mm, y_mm in wafer_fields:
            w_top.add(gdstk.Reference(field_cell, origin=(x_mm * 1000.0, y_mm * 1000.0)))
        w_path = out_dir / 'A12_wafer_fields.gds'
        w_lib.write_gds(str(w_path))
        written['wafer_reference'] = str(w_path)
    return written


def _classify_die_cells(gds_path, dies, squares):
    """Check metal area inside every A12 opening per unique design cell."""
    lib = gdstk.read_gds(str(gds_path))
    by_name = {c.name: c for c in lib.cells}
    design_status = {}
    for design_name in sorted({d['design'] for d in dies}):
        cell = by_name[design_name]
        metal_polys = [Polygon(p.points).buffer(0) for p in cell.get_polygons() if p.layer == BASE_LAYER]
        metal_u = shapely.union_all(metal_polys) if metal_polys else Polygon()
        per_group = {}
        for gname, sq in squares.items():
            area = float(metal_u.intersection(sq).area)
            per_group[gname] = 'clear' if area <= AREA_TOL_UM2 else f'conflict({area:.1f}um2)'
        design_status[design_name] = per_group

    total_counts = Counter()
    by_design = {}
    for d in dies:
        dn = d['design']
        by_design.setdefault(dn, Counter())
        for gname, st in design_status[dn].items():
            tag = 'clear' if st == 'clear' else 'conflict'
            total_counts[f'die/{tag}'] += 1
            by_design[dn][tag] += 1
    return {
        'chip': dict(total_counts),
        'wafer': {},
        'chip_by_design': {k: dict(v) for k, v in by_design.items()},
    }


def _polygon_signature(cell, layer):
    return Counter(
        tuple((round(x * 1000), round(y * 1000)) for x, y in p.points)
        for p in cell.get_polygons() if p.layer == layer
    )


def run_overlay(assembly_dir, stem, out_dir, deliver_gds, classify_stem=None,
                wafer_fields=(), replace=False):
    assembly_dir = Path(assembly_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pos_gds = assembly_dir / f'{stem}.gds'
    asm_json = assembly_dir / f'{stem}.json'
    asm = json.loads(asm_json.read_text(encoding='utf-8'))
    dies = asm['dies']
    squares = _a12_squares_um()

    classify_gds = assembly_dir / f'{classify_stem or stem}.gds'
    summary_before = _classify_die_cells(classify_gds, dies, squares)
    summary_after = _classify_die_cells(pos_gds, dies, squares)

    result = {
        'summary_before_opening': summary_before,
        'summary_after_opening': summary_after,
    }

    if deliver_gds is not None:
        deliver_path = Path(deliver_gds)
        if deliver_path.exists() and not replace:
            raise FileExistsError(f'{deliver_path} exists; pass --replace to overwrite')
        deliver_path.parent.mkdir(parents=True, exist_ok=True)

        lib = gdstk.read_gds(str(pos_gds))
        top, = lib.top_level()
        pos_sig = _polygon_signature(top, BASE_LAYER)

        chip_cell = lib.new_cell('TSRI_A12_CHIP')
        chip_polys = [g for gs in a12_groups().values() for g in gs]
        chip_cell.add(*_to_gdstk_polys(chip_polys, layer=CHIP_LAYER))
        for d in dies:
            top.add(gdstk.Reference(chip_cell, origin=tuple(d['centre_um'])))

        if wafer_fields:
            wf_cell = lib.new_cell('TSRI_A12_WAFER')
            full_polys = [g for gs in a12_groups(include_corners=True).values() for g in gs]
            wf_cell.add(*_to_gdstk_polys(full_polys, layer=WAFER_LAYER))
            for x_mm, y_mm in wafer_fields:
                top.add(gdstk.Reference(wf_cell, origin=(x_mm * 1000.0, y_mm * 1000.0)))

        lib.write_gds(str(deliver_path), max_points=8190)

        # Read back and verify layer 1/0 base metal is vertex-identical to positive GDS
        chk_lib = gdstk.read_gds(str(deliver_path))
        chk_top, = chk_lib.top_level()
        eb_sig = _polygon_signature(chk_top, BASE_LAYER)
        matches = (eb_sig == pos_sig)
        if not matches:
            raise RuntimeError('E-beam base metal (layer 1/0) does not match positive GDS')
        result['ebeam_base_metal_matches_positive'] = matches
        result['ebeam_sha256'] = sha256(deliver_path.read_bytes()).hexdigest()

    rep_path = out_dir / f'{stem}_tsri_a12.json'
    rep_path.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    return result


def _parse_fields(items):
    out = []
    for s in items:
        x_s, y_s = s.split(',', 1)
        out.append((float(x_s), float(y_s)))
    return out


def main():
    ap = argparse.ArgumentParser(description='TSRI A12 e-beam marks on a wafer assembly.')
    ap.add_argument('folder', nargs='?', type=Path, help='delivery folder (optional)')
    ap.add_argument('--write-references', type=Path, help='only write A12_chip_single_die.gds (2/0) to DIR')
    ap.add_argument('--wafer-field', action='append', default=[], metavar='X,Y',
                    help='wafer A12 field centre in mm (repeatable; default: none)')
    ap.add_argument('--assembly-dir', type=Path, help='folder of <stem>.gds/.json')
    ap.add_argument('--stem', help='opened (positive) assembly stem')
    ap.add_argument('--classify-stem', help='assembly before openings (default: --stem)')
    ap.add_argument('--out', type=Path, help='folder for report JSON (tmp/)')
    ap.add_argument('--deliver-gds', type=Path, help='TSRI e-beam GDS output path')
    ap.add_argument('--replace', action='store_true')
    args = ap.parse_args()

    fields = _parse_fields(args.wafer_field)
    if args.write_references:
        res = write_references(args.write_references, wafer_fields=fields)
        print(json.dumps(res, indent=2))
        return

    if not args.assembly_dir or not args.stem or not args.out:
        ap.error('--assembly-dir, --stem, and --out are required unless --write-references is used')
    res = run_overlay(
        args.assembly_dir, args.stem, args.out, args.deliver_gds,
        classify_stem=args.classify_stem, wafer_fields=fields, replace=args.replace,
    )
    print(json.dumps(res, indent=2))


if __name__ == '__main__':
    main()
