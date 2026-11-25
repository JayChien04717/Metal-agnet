"""Design2_FTF_Andersen r12 reference: r11 reference with the fluxonia drawn by FluxoniumDelftDolan.

The full wafer pipeline (code/general/export_wafer_test_junctions.py) cannot run at present
(it imports DesignTestJunction/layout/build.py, removed in commit 6c55641). This script takes
the confirmed r11 reference GDS and applies layout3_full.use_fluxonium_delft(), the same swap
recover() now applies: inside the FluxoniumDelftDolan pockets / line gaps the r11 metal is replaced by
the component metal (the TC coupler crescents are kept). Everything else is r11 unchanged.

    python -X utf8 -B Design2_FTF_Andersen/layout/export_fluxonium_delft_r12.py

Writes gds/Design2_2026-11-22_wafer_r12_reference.gds (refuses to overwrite) and prints the
metal difference to r11. No fab GDS (Stage 3).
"""
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import configure as _configure
_configure(__file__)

import sys
from pathlib import Path

import gdstk
import numpy as np
from shapely import get_parts
from shapely.geometry import Polygon
from shapely.ops import unary_union

sys.path.insert(0, str(Path(__file__).resolve().parent))
from layout3_full import use_fluxonium_delft  # noqa: E402

HOME = Path(__file__).resolve().parents[1]
PARENT = HOME / 'gds' / 'Design2_2026-11-18_wafer_r11_reference.gds'
OUTPUT = HOME / 'gds' / 'Design2_2026-11-22_wafer_r12_reference.gds'


def read_metal(path):
    lib = gdstk.read_gds(str(path))
    top, = lib.top_level()
    other = sorted({p.layer for p in top.get_polygons()} - {1})
    if other:
        raise ValueError(f'{path.name}: unexpected layers {other}')
    return unary_union([Polygon(p.points).buffer(0) for p in top.get_polygons()])


def write_metal(metal, path):
    lib = gdstk.Library(unit=1e-6, precision=1e-9)
    top = lib.new_cell('TOP')
    for part in get_parts(metal):
        shell = gdstk.Polygon(np.asarray(part.exterior.coords)[:-1], layer=1)
        holes = [gdstk.Polygon(np.asarray(r.coords)[:-1]) for r in part.interiors]
        for poly in (gdstk.boolean(shell, holes, 'not', layer=1, precision=.001) if holes else [shell]):
            top.add(*poly.fracture(max_points=199, precision=.001))
    lib.write_gds(str(path), max_points=199)


def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    parent = read_metal(PARENT)
    metal = use_fluxonium_delft(parent, fab=False)
    write_metal(metal, OUTPUT)
    written = read_metal(OUTPUT)
    removed, added = parent.difference(written), written.difference(parent)
    print(f'{OUTPUT.name}: metal removed vs r11 {removed.area:.3f} um^2, added {added.area:.3f} um^2')
    for part in sorted(get_parts(removed), key=lambda g: -g.area):
        if part.area > 1e-2:
            print('  removed', round(part.area, 3), [round(v, 2) for v in part.bounds])
    for part in sorted(get_parts(added), key=lambda g: -g.area):
        if part.area > 1e-2:
            print('  added', round(part.area, 3), [round(v, 2) for v in part.bounds])


if __name__ == '__main__':
    main()
