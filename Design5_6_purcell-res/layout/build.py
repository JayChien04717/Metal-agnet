# -*- coding: utf-8 -*-
r"""Single layout entry point for Design5_6_purcell-res, rebuilt from the saved simulation GDS (2026-11-25).

The original layout scripts (design5_layout.py, chip6_*.py) were removed with Design4_Single_purcell-res. The layout
geometry comes from the revision's simulation GDS (GdsLayout component); masks are rendered from it with the current
repository GDS settings (skill/gds_export.md):
    purcell_resonator_pair_r1 -> gds/Design5_2026-11-19_purcell_resonator_pair_r1_sim.gds
  variant='sim': render with fab=False (no cheese), for Q3D / HFSS and figures;
  variant='fab': same geometry, exported with fab=True and no renderer cheese, merged to positive metal on 1/0, then
                 perforated by code/fab/wafer_cheese_grid.replace_wafer_cheese (5 x 5 um holes, 40 um pitch, lattice
                 symmetric about the die centre, hole centres within +-4500 um, 50 um etch keep-out, launch clearance profile;
                 skill/gds_export.md section 3), the same flow as layout/export_tsri.py.
The delivered PR01 mask (20 um cheese, die ID 'PR01') is not regenerated here. Geometry only: no pins or part names.

    python -X utf8 -B Design5_6_purcell-res/layout/build.py [--variant sim|fab] [--out <gds>] [--check]

--check (or --out) exports the design and compares its metal with the source sim GDS (cheese holes filled) and,
for the fab variant, reports the hole count, size and pitch.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'layout'), str(ROOT / 'code' / 'fab'),
                str(HERE / 'components')]
import project_paths
project_paths.configure(__file__)

import numpy as np
from qiskit_metal import designs
from shapely.geometry import Polygon
from gds_layout import GdsLayout, read_gds_metal

DESIGN_DIR = HERE.parent
CHEESE_HOLE_UM2 = 100.0           # holes below this area are cheese (5 x 5 um = 25 um^2)
REVISIONS = {
    'purcell_resonator_pair_r1': dict(gds=DESIGN_DIR / 'gds' / 'Design5_2026-11-19_purcell_resonator_pair_r1_sim.gds'),
}


def build_design(variant='sim', overrides=None, revision='purcell_resonator_pair_r1'):
    """Rebuild Design5_6 from the revision's sim GDS. The variant only changes how the GDS is exported."""
    if variant not in ('sim', 'fab'):
        raise ValueError(f'variant {variant!r}')
    gds = Path((overrides or {}).get('gds_file', REVISIONS[revision]['gds']))
    d = designs.DesignPlanar()
    d.overwrite_enabled = True
    d.chips.main.material = 'silicon_eps11p5'      # Design5_6 used Si eps_r 11.5; run solvers with SI_EPS=11.5 (save_aedt)
    d.chips.main.size.update(size_x='10mm', size_y='10mm', size_z='-500um', center_x='0mm', center_y='0mm')
    GdsLayout(d, 'chip', options=dict(gds_file=str(gds), source_format='sim'))
    return d


# Design5 is a 2-port chip (launches centred left / right at +-4.5, 0 mm), as in layout/export_tsri.py.
TWO_PORT = {
    'P_L': {'pos_x': '-4.5mm', 'pos_y': '0.0mm', 'orientation': '0'},
    'P_R': {'pos_x': '4.5mm', 'pos_y': '0.0mm', 'orientation': '180'},
}


def export_fab(design, out_gds):
    """Fab GDS: un-cheesed fab render -> merged positive metal 1/0 -> wafer_cheese_grid 40 um lattice."""
    import tempfile
    import gdstk
    from gds_tools import export_gds
    from wafer_cheese_grid import replace_wafer_cheese
    with tempfile.TemporaryDirectory() as tmp:
        raw = Path(tmp) / 'raw_fab.gds'
        export_gds(design, str(raw), fab=True, cheese=False)
        _, metal = read_gds_metal(raw, source_format='all')            # mm
    lib = gdstk.Library(unit=1e-6, precision=1e-9)
    top = lib.new_cell('TOP')
    for poly in getattr(metal, 'geoms', [metal]):
        ext = gdstk.Polygon([(x * 1e3, y * 1e3) for x, y in list(poly.exterior.coords)[:-1]], layer=1, datatype=0)
        holes = [gdstk.Polygon([(x * 1e3, y * 1e3) for x, y in list(h.coords)[:-1]]) for h in poly.interiors]
        top.add(*(gdstk.boolean([ext], holes, 'not', precision=0.001, layer=1, datatype=0) if holes else [ext]))
    Path(out_gds).parent.mkdir(parents=True, exist_ok=True)
    lib.write_gds(str(out_gds), max_points=199)
    return replace_wafer_cheese(out_gds, launches=TWO_PORT)


def compare(revision, out_gds):
    """Metal of the export (all layer-1 datatypes, cheese holes filled) vs the source sim GDS, plus hole statistics."""
    _, src = read_gds_metal(REVISIONS[revision]['gds'], source_format='sim')
    _, raw = read_gds_metal(out_gds, source_format='all')
    _, out = read_gds_metal(out_gds, source_format='all', fill_holes_um2=CHEESE_HOLE_UM2)
    x = src.symmetric_difference(out)
    holes = [Polygon(h) for g in getattr(raw, 'geoms', [raw]) for h in g.interiors if Polygon(h).area * 1e6 < CHEESE_HOLE_UM2]
    res = dict(xor_area_um2=float(x.area * 1e6), metal_area_mm2=[float(src.area), float(out.area)],
               ok=bool(x.area * 1e6 <= 1.0), n_cheese_holes=len(holes))
    if holes:
        c = np.array([[h.centroid.x, h.centroid.y] for h in holes]) * 1e3          # um
        res.update(hole_area_um2=float(np.median([h.area for h in holes]) * 1e6),
                   pitch_um=[float(np.median(np.diff(np.unique(np.round(c[:, i], 2))))) for i in (0, 1)],
                   hole_extent_um=[float(c[:, 0].min()), float(c[:, 0].max()), float(c[:, 1].min()), float(c[:, 1].max())],
                   phase_mod_pitch_um=sorted({round(float(v), 1) for v in np.round(c[:, 0], 1) % 40.0}))
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--revision', default='purcell_resonator_pair_r1', choices=sorted(REVISIONS))
    ap.add_argument('--variant', default='sim', choices=('sim', 'fab'))
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--out', default='', help='GDS path for the export (default: system temp folder)')
    a = ap.parse_args()
    d = build_design(variant=a.variant, revision=a.revision)
    t = d.qgeometry.tables['poly']
    sub = t['subtract'].astype(bool)
    print(f'{a.revision} {a.variant}: {int((~sub).sum())} metal islands, {int(sub.sum())} etch pieces '
          f'from {REVISIONS[a.revision]["gds"].name}')
    if a.check or a.out:
        import tempfile
        from gds_tools import export_gds
        out = Path(a.out) if a.out else Path(tempfile.gettempdir()) / f'Design5_6_{a.revision}_{a.variant}_rebuilt.gds'
        if a.variant == 'fab':
            cheese = export_fab(d, out)
            print('CHEESE', {k: cheese[k] for k in ('old_holes', 'new_holes', 'pitch_um', 'hole_um', 'edge_um') if k in cheese})
        else:
            export_gds(d, str(out), fab=False, cheese=False)
        print('CHECK', compare(a.revision, out), '->', out)
