# -*- coding: utf-8 -*-
r"""Export Design6 (TXT02) TSRI single-die fab GDS and review figure.

Starts from the adopted Design6_TXT wafer mask, removes the old corner test
junctions (which overlapped TSRI A12 coarse crosses), verifies >= 100 um
clearance to BOTH TSRI A12 (8 squares) and NTHU TiAu (4 squares), cuts the
8 TSRI A12 bare-rectangle ground openings, places 4 vertical transmon-style
``HookedTestJunction`` coupons in open ground regions, rebuilds the 40 um wafer
cheese lattice, and renders a review PNG.

Usage (qmetal env):
    python -X utf8 -B Design6_TXT/layout/export_tsri.py
"""
from hashlib import sha256
import json
from pathlib import Path
import sys

_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'code' / 'general' / 'project_paths.py').is_file())
sys.path.insert(0, str(_ROOT / 'code' / 'general'))
import project_paths
project_paths.configure(__file__)

DESIGN_DIR = project_paths.design_dir(__file__)

import geopandas as gpd
if not hasattr(gpd.GeoDataFrame, 'append'):
    gpd.GeoDataFrame.append = gpd.GeoDataFrame._append

import gdstk
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath
import matplotlib.pyplot as plt
import numpy as np
import shapely
from shapely.affinity import scale
from shapely.geometry import Polygon, box

from qiskit_metal import Dict, designs
from junction_test_coupon import HookedTestJunction
from nthu_tiau_marker import chip_marker_sites, site_opening
from tsri_a12_mark import a12_groups
from wafer_cheese_grid import replace_wafer_cheese
from wafer_pcbport import wafer_pcb_ports

SRC_GDS = DESIGN_DIR / 'gds' / 'Design6_2026-11-23_wafer_purcell_resonator_pair_r1_txt02_fab.gds'
OUT_GDS = DESIGN_DIR / 'gds' / 'Design6_2026-11-24_tsri_r2_txt02_fab.gds'
OUT_PNG = DESIGN_DIR / 'fig' / 'Design6_2026-11-24_tsri_r2_txt02_review.png'

# 4 vertical transmon (T) HookedTestJunction coupons in open ground regions (TXT is transmon-TC-transmon).
COUPONS = {
    'TJ_UL': dict(pos_x='-4.05mm', pos_y='0.00mm', orientation='90', landing_gap='5um', hook_offset='2um'),
    'TJ_UR': dict(pos_x='4.05mm', pos_y='0.00mm', orientation='90', landing_gap='5um', hook_offset='2um'),
    'TJ_LL': dict(pos_x='-2.10mm', pos_y='-3.10mm', orientation='90', landing_gap='5um', hook_offset='2um'),
    'TJ_LR': dict(pos_x='2.10mm', pos_y='-3.10mm', orientation='90', landing_gap='5um', hook_offset='2um'),
}
COUPON_COMMON = dict(
    pocket_width='820um', pocket_height='230um',
    pad_length='300um', pad_width='100um', pad_corner_radius='10um',
    pad_body_gap='100um', stem_width='10um',
    hook_width='2um', hook_depth='4um', fab_contacts='True',
)


def _coupon_geoms_um():
    d = designs.DesignPlanar()
    d.overwrite_enabled = True
    pockets, pads = {}, []
    for name, spec in COUPONS.items():
        c = HookedTestJunction(d, name, options=Dict(**COUPON_COMMON, **spec))
        pockets[name] = scale(c.keepouts()['pocket'], 1000.0, 1000.0, origin=(0, 0))
    for row in d.qgeometry.tables['poly'].itertuples():
        if not row.subtract:
            pads.append(scale(row.geometry, 1000.0, 1000.0, origin=(0, 0)))
    return shapely.union_all(list(pockets.values())), shapely.union_all(pads), pockets


def _site_squares_um(margin_um=50.0):
    tsri = {}
    for name, gs in a12_groups().items():
        x0, y0, x1, y1 = shapely.union_all(gs).bounds
        tsri[name] = box(x0 - margin_um, y0 - margin_um, x1 + margin_um, y1 + margin_um)
    nthu = {name: site_opening(gs, margin_um) for name, gs in chip_marker_sites().items()}
    return tsri, nthu


def _to_gdstk(geom, layer=1, datatype=0):
    out = []
    for poly in shapely.get_parts(geom):
        if poly.is_empty or poly.area <= 0:
            continue
        ext = gdstk.Polygon(list(poly.exterior.coords)[:-1], layer=layer, datatype=datatype)
        if poly.interiors:
            holes = [gdstk.Polygon(list(h.coords)[:-1], layer=layer, datatype=datatype)
                     for h in poly.interiors]
            out.extend(gdstk.boolean([ext], holes, 'not', precision=0.001, layer=layer, datatype=datatype))
        else:
            out.append(ext)
    return out


def _render_review_png(gds_path, png_path, clearance_report, title):
    lib = gdstk.read_gds(str(gds_path))
    top, = lib.top_level()
    metal_u = shapely.union_all([Polygon(p.points).buffer(0) for p in top.get_polygons() if p.layer == 1])
    groups = a12_groups()
    tsri_sq, _ = _site_squares_um()

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
        return MplPath(np.vstack(verts), np.concatenate(codes)) if verts else None

    def _draw_geom(ax, geom, facecolor='#8b96a3'):
        for poly in shapely.get_parts(geom):
            if poly.is_empty or poly.area <= 0:
                continue
            p = _poly_to_path(poly)
            if p is not None:
                ax.add_patch(PathPatch(p, facecolor=facecolor, edgecolor='none', linewidth=0, zorder=1))

    fig = plt.figure(figsize=(18, 9), dpi=140)
    gs = fig.add_gridspec(2, 6, width_ratios=[2.2, 2.2, 1, 1, 1, 1], wspace=0.28, hspace=0.28)
    ax_full = fig.add_subplot(gs[:, 0:2])
    _draw_geom(ax_full, metal_u)
    for name, sq in tsri_sq.items():
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
        sq = tsri_sq[name]
        x0, y0, x1, y1 = sq.bounds
        pad = 220.0
        win = box(x0 - pad, y0 - pad, x1 + pad, y1 + pad)
        _draw_geom(ax, metal_u.intersection(win))
        sx, sy = sq.exterior.xy
        ax.plot([x / 1000 for x in sx], [y / 1000 for y in sy], color='#c62828', lw=1.0, zorder=3)
        for m in groups[name]:
            mx, my = m.exterior.xy
            ax.fill([x / 1000 for x in mx], [y / 1000 for y in my],
                    facecolor='#1565c0', edgecolor='#1565c0', lw=0.8, zorder=4)
        cl = clearance_report.get(f'TSRI:{name}', {}).get('distance_um')
        ax.set_xlim((x0 - pad) / 1000, (x1 + pad) / 1000)
        ax.set_ylim((y0 - pad) / 1000, (y1 + pad) / 1000)
        ax.set_aspect('equal')
        ax.set_title(f'{name}: {cl} um' if cl is not None else name, fontsize=9)
        ax.tick_params(labelsize=7)

    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(png_path), bbox_inches='tight')
    plt.close(fig)


def main():
    lib = gdstk.read_gds(str(SRC_GDS))
    top, = lib.top_level()
    raw_metal = shapely.union_all([Polygon(p.points).buffer(0) for p in top.get_polygons() if p.layer == 1])
    clean_parts = []
    for p in shapely.get_parts(raw_metal):
        clean_parts.append(Polygon(p.exterior, [r for r in p.interiors if Polygon(r).area >= 30.0]).buffer(0))
    metal = shapely.union_all(clean_parts)

    # Fill the 4 old corner test-junction boxes (at +-4.20, +-4.00 mm + F/T letters,
    # staying below y=4.75 mm so upper-left TXT02 die ID at y=4.83..4.95 mm is untouched).
    old_corners = shapely.union_all([
        box(sx * 3.55e3, sy * 3.45e3, sx * 4.75e3, sy * 4.65e3)
        for sx in (-1, 1) for sy in (-1, 1)
    ])
    metal = metal.union(old_corners)

    tsri_sq, nthu_sq = _site_squares_um()
    coupon_cut, coupon_pads, coupon_pockets = _coupon_geoms_um()

    active_etch = box(-5000, -5000, 5000, 5000).difference(metal)
    all_circuitry = active_etch.union(coupon_cut)
    clearances = {}
    violations = []
    for prefix, sqs in (('TSRI', tsri_sq), ('NTHU', nthu_sq)):
        for name, sq in sqs.items():
            d = round(float(sq.distance(all_circuitry)), 2)
            clearances[f'{prefix}:{name}'] = {'distance_um': d}
            if d < 100.0:
                violations.append((f'{prefix}:{name}', d))
    for name, pk in coupon_pockets.items():
        d = round(float(pk.distance(active_etch)), 2)
        clearances[f'{name}:to_circuitry'] = {'distance_um': d}
        if d < 100.0:
            violations.append((f'{name}:to_circuitry', d))
    if violations:
        raise ValueError(f'Design6 TSRI clearance below 100 um: {violations}')

    metal = metal.difference(shapely.union_all(list(tsri_sq.values())))
    metal = metal.difference(coupon_cut).union(coupon_pads)

    OUT_GDS.parent.mkdir(parents=True, exist_ok=True)
    out_lib = gdstk.Library(unit=1e-6, precision=1e-9)
    out_top = out_lib.new_cell('TOP')
    out_top.add(*_to_gdstk(metal, layer=1, datatype=0))
    out_lib.write_gds(str(OUT_GDS), max_points=199)

    cheese = replace_wafer_cheese(OUT_GDS, launches=wafer_pcb_ports)
    _render_review_png(
        OUT_GDS, OUT_PNG, clearances,
        'Design6 (TXT02) TSRI r2: A12 openings + 4 test junctions + 40 um cheese',
    )
    print(json.dumps({
        'gds': str(OUT_GDS.relative_to(_ROOT)),
        'sha256': sha256(OUT_GDS.read_bytes()).hexdigest(),
        'review_png': str(OUT_PNG.relative_to(_ROOT)),
        'clearances': clearances,
        'cheese': {k: cheese[k] for k in ('new_holes', 'pitch_um', 'hole_um', 'edge_um')},
    }, indent=2))


if __name__ == '__main__':
    main()
