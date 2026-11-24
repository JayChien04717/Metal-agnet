# -*- coding: utf-8 -*-
r"""Export Design2 TSRI single-die fab GDS and review figure.

Rebuilds Design2 from ``layout3_full.recover`` with:
- TSRI feedline clearance: ``readout_in`` / ``readout_out`` routed at ``|x| = 2750 um``
  down to ``y = 2500 um`` (257.5 um clear of ``fine_UL`` / ``fine_UR``, vs 7.5 um before)
- 8 TSRI A12 bare-rectangle ground openings (``TSRIA12Mark``, no mark metal in 1/0)
- 4 ``HookedTestJunction`` coupons (2x F, 2x T) in open ground regions
- Upper-left die ID ``FTFC00`` (clear of ``coarse_UL``)
- 40 um wafer cheese lattice via ``code/general/wafer_cheese_grid.py``

Usage (qmetal env):
    python -X utf8 -B Design2_FTF_Andersen/layout/export_tsri.py
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
from shapely.affinity import scale, translate
from shapely.geometry import Polygon, box

from qiskit_metal import Dict, designs
from junction_test_coupon import HookedTestJunction
import layout3_full as L2
from tsri_a12_mark import a12_groups
from wafer_cheese_grid import replace_wafer_cheese
from wafer_pcbport import wafer_pcb_ports
from xmark import LABEL_SIZE_UM, label_geometry

OUT_GDS = DESIGN_DIR / 'gds' / 'Design2_2026-11-24_tsri_r12_fab.gds'
OUT_PNG = DESIGN_DIR / 'fig' / 'Design2_2026-11-24_tsri_r12_review.png'

# 4 vertical HookedTestJunction coupons (2x F, 2x T) in open ground regions.
COUPONS = {
    'TJ_UL': dict(pos_x='-2.60mm', pos_y='0.40mm', orientation='90', landing_gap='5um', hook_offset='2um'),
    'TJ_UR': dict(pos_x='2.60mm', pos_y='0.40mm', orientation='90', landing_gap='30um', hook_offset='0um'),
    'TJ_LL': dict(pos_x='-2.10mm', pos_y='-3.10mm', orientation='90', landing_gap='30um', hook_offset='0um'),
    'TJ_LR': dict(pos_x='2.10mm', pos_y='-3.10mm', orientation='90', landing_gap='5um', hook_offset='2um'),
}
COUPON_COMMON = dict(
    pocket_width='820um', pocket_height='230um',
    pad_length='300um', pad_width='100um', pad_corner_radius='10um',
    pad_body_gap='100um', stem_width='10um',
    hook_width='2um', hook_depth='4um', fab_contacts='True',
)


def _coupon_geoms_um():
    """Return (pocket_union_um, pad_union_um, {name: pocket_um}) from HookedTestJunction."""
    d = designs.DesignPlanar()
    d.overwrite_enabled = True
    pockets = {}
    pads = []
    for name, spec in COUPONS.items():
        c = HookedTestJunction(d, name, options=Dict(**COUPON_COMMON, **spec))
        pockets[name] = scale(c.keepouts()['pocket'], 1000.0, 1000.0, origin=(0, 0))
    for row in d.qgeometry.tables['poly'].itertuples():
        if not row.subtract:
            pads.append(scale(row.geometry, 1000.0, 1000.0, origin=(0, 0)))
    return shapely.union_all(list(pockets.values())), shapely.union_all(pads), pockets


def _a12_squares_um(margin_um=50.0):
    out = {}
    for name, gs in a12_groups().items():
        x0, y0, x1, y1 = shapely.union_all(gs).bounds
        out[name] = box(x0 - margin_um, y0 - margin_um, x1 + margin_um, y1 + margin_um)
    return out


def _measure_clearances(base_metal_um, squares, coupon_pockets):
    """Measure clearance [um] from every A12 square and coupon pocket to active etch."""
    die = box(-5000, -5000, 5000, 5000)
    # Remove old small cheese holes (< 5.2 um wide) before measuring active circuitry etch.
    raw_etch = die.difference(base_metal_um)
    active_etch = [e for e in shapely.get_parts(raw_etch)
                   if e.area > 0.2 and not e.buffer(-2.6).is_empty]
    etch_u = shapely.union_all(active_etch)

    report = {}
    violations = []
    for name, sq in squares.items():
        inside = etch_u.intersection(sq).area
        dist = round(float(sq.distance(etch_u)), 2)
        report[f'TSRI_A12:{name}'] = {'inside_um2': round(float(inside), 3), 'distance_um': dist}
        if inside > 0.1 or dist < 100.0:
            violations.append((f'TSRI_A12:{name}', dist))
    sq_union = shapely.union_all(list(squares.values()))
    for name, pk in coupon_pockets.items():
        d_etch = round(float(pk.distance(etch_u)), 2)
        d_sq = round(float(pk.distance(sq_union)), 2)
        dist = min(d_etch, d_sq)
        report[f'{name}:pocket'] = {'to_etch_um': d_etch, 'to_a12_um': d_sq, 'distance_um': dist}
        if dist < 100.0:
            violations.append((f'{name}:pocket', dist))
    if violations:
        raise ValueError(f'Design2 TSRI clearance below 100 um: {violations}')
    return report


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
    metal = [Polygon(p.points).buffer(0) for p in top.get_polygons() if p.layer == 1]
    metal_u = shapely.union_all(metal)
    groups = a12_groups()
    squares = _a12_squares_um()

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

    def _draw_geom(ax, geom, facecolor='#8b96a3'):
        for poly in shapely.get_parts(geom):
            if poly.is_empty or poly.area <= 0:
                continue
            p = _poly_to_path(poly)
            if p is not None:
                ax.add_patch(PathPatch(p, facecolor=facecolor, edgecolor='none',
                                       linewidth=0, zorder=1))

    fig = plt.figure(figsize=(18, 9), dpi=140)
    gs = fig.add_gridspec(2, 6, width_ratios=[2.2, 2.2, 1, 1, 1, 1], wspace=0.28, hspace=0.28)
    ax_full = fig.add_subplot(gs[:, 0:2])
    _draw_geom(ax_full, metal_u)
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
        _draw_geom(ax, metal_u.intersection(win))
        sx, sy = sq.exterior.xy
        ax.plot([x / 1000 for x in sx], [y / 1000 for y in sy], color='#c62828', lw=1.0, zorder=3)
        for m in groups[name]:
            mx, my = m.exterior.xy
            ax.fill([x / 1000 for x in mx], [y / 1000 for y in my],
                    facecolor='#1565c0', edgecolor='#1565c0', lw=0.8, zorder=4)
        cl = clearance_report.get(f'TSRI_A12:{name}', {}).get('distance_um')
        ax.set_xlim((x0 - pad) / 1000, (x1 + pad) / 1000)
        ax.set_ylim((y0 - pad) / 1000, (y1 + pad) / 1000)
        ax.set_aspect('equal')
        ax.set_title(f'{name}: {cl} um' if cl is not None else name, fontsize=9)
        ax.tick_params(labelsize=7)

    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(png_path), bbox_inches='tight')
    plt.close(fig)


def _shift_feed_corners_inward(metal):
    """Shift the two outer feedline bends from |x|=3000 um inward to |x|=2750 um.

    Fills the old bends inside |x| in [2600, 3050], y in [2450, 3320] back to
    solid ground and draws matching 15/10 um CPW bends (bend_radius=75 um) from
    (+-2600, 3270) -> (+-2750, 3270) -> (+-2750, 2500) -> (+-2600, 2500), giving
    257.5 um clearance to fine_UL / fine_UR (at |x|=3025 um).
    """
    for sx in (-1, 1):
        x0, x1 = sorted((sx * 2600.0, sx * 3050.0))
        patch = box(x0, 2450.0, x1, 3320.0)
        metal = metal.union(patch)
        pts = [(sx * 2595.0, 3270.0), (sx * 2750.0, 3270.0),
               (sx * 2750.0, 2500.0), (sx * 2595.0, 2500.0)]
        trace_poly, = gdstk.FlexPath(pts, 15.0, bend_radius=75.0, tolerance=0.005).to_polygons()
        cut_poly, = gdstk.FlexPath(pts, 35.0, bend_radius=75.0, tolerance=0.005).to_polygons()
        metal = metal.difference(Polygon(cut_poly.points).buffer(0)).union(Polygon(trace_poly.points).buffer(0))
    return metal


def main():
    source_gds = _ROOT / 'Design2_FTF_Andersen' / 'gds' / 'previous' / 'layout3 ftf.gds'
    if not source_gds.is_file():
        source_gds = next(_ROOT.rglob('layout3 ftf.gds'))
    ground, nets = L2.recover(
        str(source_gds), fab=True, standard_launches=True,
        port_options=wafer_pcb_ports, tc_flux_offset_um=3.5, squid_ground_tail_um=7,
    )
    # Strip legacy 20 um cheese holes before boolean ops (replace_wafer_cheese rebuilds 40 um lattice).
    clean_ground = Polygon(ground.exterior, [r for r in ground.interiors if Polygon(r).area >= 30.0]).buffer(0)
    metal = shapely.union_all([clean_ground, *nets.values()])
    metal = _shift_feed_corners_inward(metal)

    squares = _a12_squares_um()
    coupon_cut, coupon_pads, coupon_pockets = _coupon_geoms_um()
    clearances = _measure_clearances(metal, squares, coupon_pockets)

    # Cut the 8 TSRI A12 bare squares, add 4 test-junction coupons, and etch upper-left die ID.
    metal = metal.difference(shapely.union_all(list(squares.values())))
    metal = metal.difference(coupon_cut).union(coupon_pads)
    glyph = label_geometry('FTFC00', LABEL_SIZE_UM, 10, 400, 30)
    x0, y0, x1, y1 = glyph.bounds
    die_id = translate(glyph, -5000.0 + 50.0 - x0, 5000.0 - 50.0 - y1)
    metal = metal.difference(die_id)

    OUT_GDS.parent.mkdir(parents=True, exist_ok=True)
    lib = gdstk.Library(unit=1e-6, precision=1e-9)
    top = lib.new_cell('TOP')
    top.add(*_to_gdstk(metal, layer=1, datatype=0))
    lib.write_gds(str(OUT_GDS), max_points=199)

    cheese = replace_wafer_cheese(OUT_GDS, launches=wafer_pcb_ports)
    _render_review_png(
        OUT_GDS, OUT_PNG, clearances,
        'Design2 (FTFC00) TSRI r12: shifted feedline + A12 openings + 4 test junctions + 40 um cheese',
    )
    summary = {
        'gds': str(OUT_GDS.relative_to(_ROOT)),
        'sha256': sha256(OUT_GDS.read_bytes()).hexdigest(),
        'review_png': str(OUT_PNG.relative_to(_ROOT)),
        'die_id_bounds_um': [round(v, 2) for v in die_id.bounds],
        'clearances': clearances,
        'cheese': {k: cheese[k] for k in ('new_holes', 'pitch_um', 'hole_um', 'edge_um')},
    }
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
