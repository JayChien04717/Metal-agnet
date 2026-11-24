# -*- coding: utf-8 -*-
r"""Export Design3 TSRI single-die fab GDS and review figure.

Uses ``build.py`` with ``placement_tsri.json`` and ``routes_tsri.json``:
- TSRI A12 bare-rectangle ground openings (via ``TSRIA12Mark``, no mark metal in 1/0)
- 4 ``HookedTestJunction`` coupons (2x F, 2x T) in open ground regions
- Upper-left die ID ``FTF00`` (clear of ``coarse_UL``)
- 40 um wafer cheese lattice via ``code/general/wafer_cheese_grid.py``
- Keep-out clearance >= 100 um verified in QGeometry before export

Usage (qmetal env):
    python -X utf8 -B Design3_FTF_NTHU/layout/export_tsri.py
"""
from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'code' / 'general' / 'project_paths.py').is_file())
sys.path.insert(0, str(_ROOT / 'code' / 'general'))
import project_paths
project_paths.configure(__file__)

DESIGN_DIR = project_paths.design_dir(__file__)
LAYOUT_DIR = Path(__file__).resolve().parent
if str(LAYOUT_DIR) not in sys.path:
    sys.path.insert(0, str(LAYOUT_DIR))

import geopandas as gpd
if not hasattr(gpd.GeoDataFrame, 'append'):
    gpd.GeoDataFrame.append = gpd.GeoDataFrame._append

import gdstk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shapely
from shapely.affinity import translate
from shapely.geometry import Polygon, box

import build as d3_build
from gds_tools import export_gds, write_wafer_gds
from tsri_a12_mark import a12_groups
from wafer_cheese_grid import replace_wafer_cheese
from wafer_pcbport import wafer_pcb_ports
from xmark import LABEL_SIZE_UM, label_geometry

OUT_GDS = DESIGN_DIR / 'gds' / 'Design3_2026-11-24_tsri_r12_fab.gds'
OUT_PNG = DESIGN_DIR / 'fig' / 'Design3_2026-11-24_tsri_r12_review.png'


def _etch_upper_left_die_id(gds_path, text='FTF00', inset_um=50.0, height_um=LABEL_SIZE_UM):
    """Etch the standalone die ID in the upper-left margin (layer 1/0)."""
    lib = gdstk.read_gds(str(gds_path))
    top, = lib.top_level()
    glyph = label_geometry(text, height_um, 10, 400, 30)
    x0, y0, x1, y1 = glyph.bounds
    placed = translate(glyph, -5000.0 + inset_um - x0, 5000.0 - inset_um - y1)
    cut_polys = [gdstk.Polygon(list(p.exterior.coords)[:-1], layer=1, datatype=0)
                 for p in shapely.get_parts(placed) if p.area > 0]
    metal = [p for p in top.get_polygons() if p.layer == 1]
    updated = gdstk.boolean(metal, cut_polys, 'not', precision=0.001, layer=1, datatype=0)
    out_lib = gdstk.Library(unit=1e-6, precision=1e-9)
    out_top = out_lib.new_cell('TOP')
    out_top.add(*updated)
    out_lib.write_gds(str(gds_path), max_points=199)
    return [round(v, 2) for v in placed.bounds]


def _render_review_png(gds_path, png_path, keepout_report, title):
    """Render whole-chip metal + 8 TSRI A12 square zoom-ins with blue A12 mark overlay."""
    lib = gdstk.read_gds(str(gds_path))
    top, = lib.top_level()
    metal = [Polygon(p.points).buffer(0) for p in top.get_polygons() if p.layer == 1]
    metal_u = shapely.union_all(metal)
    groups = a12_groups()
    squares = {}
    for name, gs in groups.items():
        x0, y0, x1, y1 = shapely.union_all(gs).bounds
        squares[name] = box(x0 - 50.0, y0 - 50.0, x1 + 50.0, y1 + 50.0)

    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as MplPath
    import numpy as np

    def _poly_to_path(poly):
        verts, codes = [], []
        for ring in (poly.exterior, *poly.interiors):
            coords = np.asarray(ring.coords, dtype=float) / 1000.0
            if len(coords) < 3:
                continue
            verts.append(coords)
            c = np.full(len(coords), MplPath.LINETO, dtype=np.uint8)
            c[0] = MplPath.MOVETO
            c[-1] = MplPath.CLOSEPOLY
            codes.append(c)
        if not verts:
            return None
        return MplPath(np.vstack(verts), np.concatenate(codes))

    def _draw_geom(ax, geom, facecolor='#8b96a3', edgecolor='none', lw=0.0, zorder=1):
        for poly in shapely.get_parts(geom):
            if poly.is_empty or poly.area <= 0:
                continue
            p = _poly_to_path(poly)
            if p is not None:
                ax.add_patch(PathPatch(p, facecolor=facecolor, edgecolor=edgecolor,
                                       linewidth=lw, zorder=zorder))

    fig = plt.figure(figsize=(18, 9), dpi=140)
    gs = fig.add_gridspec(2, 6, width_ratios=[2.2, 2.2, 1, 1, 1, 1], wspace=0.28, hspace=0.28)
    ax_full = fig.add_subplot(gs[:, 0:2])
    _draw_geom(ax_full, metal_u, facecolor='#8b96a3')
    for name, sq in squares.items():
        sx, sy = sq.exterior.xy
        ax_full.plot([x / 1000 for x in sx], [y / 1000 for y in sy], color='#c62828', lw=0.9, zorder=3)
        for m in groups[name]:
            mx, my = m.exterior.xy
            ax_full.fill([x / 1000 for x in mx], [y / 1000 for y in my],
                         facecolor='#1565c0', edgecolor='#1565c0', lw=0.6, zorder=4)
    ax_full.set_xlim(-5.05, 5.05)
    ax_full.set_ylim(-5.05, 5.05)
    ax_full.set_aspect('equal')
    ax_full.set_xlabel('x (mm)')
    ax_full.set_ylabel('y (mm)')
    ax_full.set_title(title, fontsize=11)

    order = ['coarse_UL', 'fine_UL', 'fine_UR', 'coarse_UR',
             'coarse_LL', 'fine_LL', 'fine_LR', 'coarse_LR']
    for idx, name in enumerate(order):
        r, c = divmod(idx, 4)
        ax = fig.add_subplot(gs[r, c + 2])
        sq = squares[name]
        x0, y0, x1, y1 = sq.bounds
        pad = 220.0
        win = box(x0 - pad, y0 - pad, x1 + pad, y1 + pad)
        _draw_geom(ax, metal_u.intersection(win), facecolor='#8b96a3')
        sx, sy = sq.exterior.xy
        ax.plot([x / 1000 for x in sx], [y / 1000 for y in sy], color='#c62828', lw=1.0, zorder=3)
        for m in groups[name]:
            mx, my = m.exterior.xy
            ax.fill([x / 1000 for x in mx], [y / 1000 for y in my],
                    facecolor='#1565c0', edgecolor='#1565c0', lw=0.8, zorder=4)
        cl = keepout_report.get(f'TSRI_A12:{name}', {}).get('distance_um')
        ax.set_xlim((x0 - pad) / 1000, (x1 + pad) / 1000)
        ax.set_ylim((y0 - pad) / 1000, (y1 + pad) / 1000)
        ax.set_aspect('equal')
        ax.set_title(f'{name}: {cl} um' if cl is not None else name, fontsize=9)
        ax.tick_params(labelsize=7)

    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(png_path), bbox_inches='tight')
    plt.close(fig)


def main():
    design = d3_build.build_design(
        variant='fab',
        placement_path=LAYOUT_DIR / 'placement_tsri.json',
        routes_path=LAYOUT_DIR / 'routes_tsri.json',
        check_keepout_um=100.0,
    )
    OUT_GDS.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=OUT_GDS.parent) as tmp:
        raw = Path(tmp) / 'raw.gds'
        export_gds(design, str(raw), fab=False)
        write_wafer_gds(raw, OUT_GDS)
    die_id_bounds = _etch_upper_left_die_id(OUT_GDS, text='FTF00')
    cheese = replace_wafer_cheese(OUT_GDS, launches=wafer_pcb_ports)
    _render_review_png(
        OUT_GDS, OUT_PNG, design.metadata['keepout_check']['clearances'],
        'Design3 (FTF00) TSRI r12: A12 openings + 4 test junctions + 40 um cheese',
    )
    summary = {
        'gds': str(OUT_GDS.relative_to(_ROOT)),
        'sha256': sha256(OUT_GDS.read_bytes()).hexdigest(),
        'review_png': str(OUT_PNG.relative_to(_ROOT)),
        'die_id_bounds_um': die_id_bounds,
        'route_lengths': design.metadata['route_lengths'],
        'keepout_check': design.metadata['keepout_check'],
        'cheese': {k: cheese[k] for k in ('new_holes', 'pitch_um', 'hole_um', 'edge_um')},
    }
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
