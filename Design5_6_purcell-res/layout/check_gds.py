# -*- coding: utf-8 -*-
r"""GDS check of the Design5_6 layout rebuilt from GDS (layout/build.py) against the saved and delivered masks.

    python -X utf8 -B Design5_6_purcell-res/layout/check_gds.py [--revision purcell_resonator_pair_r1]

Checks (results/general/<run>/summaries/gds_check.json, figures/*.png; skill/gds_export.md section 3):
  1. re-export of build_design() (temporary GDS, deleted afterwards): units um, 10 x 10 mm outline;
  2. metal XOR fab re-export (40 um cheese filled) vs the source sim GDS of build.py;
  3. metal XOR re-export vs every delivered fab mask of the revision (layer 1/0, cheese holes filled first:
     interior holes below --hole-um2 are removed), with the differing regions listed and plotted;
  4. metal narrower than 4 um in the re-export (chip_drc C12a). The GDS carries no component names, so the
     resonator exemption cannot be applied automatically: regions are listed for review.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'layout'), str(HERE)]
import project_paths
project_paths.configure(__file__)

import gdstk
from shapely.geometry import Polygon
from shapely.ops import unary_union

import build as D56
import check_metal_width as CMW
from gds_tools import export_gds
from run_record import create_run, finish_run

DESIGN = 'Design5_6_purcell-res'
FAB = {
    'purcell_resonator_pair_r1': [
        ROOT / '20261119' / 'chips' / 'Design5_2026-11-19_purcell_resonator_pair_r1_pr01_fab.gds',   # delivered 20261119
        ROOT / DESIGN / 'gds' / 'Design5_2026-11-23_purcell_resonator_pair_r1_pr01_fab.gds',        # 40 um cheese re-mask
    ],
    'purcell_resonator_pair_r2': [],
}


def metal_um(path, datatypes):
    lib = gdstk.read_gds(str(path))
    s = lib.unit / 1e-6
    top = lib.top_level()[0]
    polys = [p for p in top.get_polygons() if p.layer == 1 and p.datatype in datatypes]
    return unary_union([Polygon(p.points * s).buffer(0) for p in polys]), lib, top


def fill_small_holes(geom, max_area):
    out = []
    for g in getattr(geom, 'geoms', [geom]):
        out.append(Polygon(g.exterior, [h for h in g.interiors if Polygon(h).area > max_area]))
    return unary_union(out)


def regions(x, min_area=1.0, limit=40):
    parts = sorted([g for g in getattr(x, 'geoms', [x]) if g.area >= min_area], key=lambda g: -g.area)
    return [dict(area_um2=round(g.area, 2), bounds_um=[round(v, 1) for v in g.bounds]) for g in parts[:limit]], len(parts)


def plot_xor(a, b, x, png, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 9))
    for g in getattr(a, 'geoms', [a]):
        ax.fill(*g.exterior.xy, color='#c8ced6', lw=0)
    for g in getattr(x, 'geoms', [x]):
        if g.area >= 1.0:
            ax.fill(*g.exterior.xy, color='red', lw=0)
            x0, y0, x1, y1 = g.bounds
            ax.add_patch(matplotlib.patches.Rectangle((x0 - 60, y0 - 60), x1 - x0 + 120, y1 - y0 + 120, fill=False, ec='red', lw=1))
    ax.set_xlim(-5000, 5000); ax.set_ylim(-5000, 5000); ax.set_aspect('equal')
    ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)'); ax.set_title(title, fontsize=10)
    fig.tight_layout(); fig.savefig(png, dpi=110); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--revision', default='purcell_resonator_pair_r1', choices=sorted(D56.REVISIONS))
    ap.add_argument('--hole-um2', type=float, default=100.0, help='fab cheese holes below this area are filled (5 x 5 um = 25)')
    a = ap.parse_args()
    src = D56.REVISIONS[a.revision]['gds']
    variant = 'fab'
    fabs = [f for f in FAB[a.revision] if f.exists()]
    folder = create_run(DESIGN, 'general', f'gds_check_{a.revision}', entrypoint=__file__,
                        inputs=[HERE / 'build.py', HERE / 'components' / 'gds_layout.py', src] + fabs,
                        parameters=dict(revision=a.revision, variant=variant, hole_um2=a.hole_um2))
    figs = folder / 'figures'; figs.mkdir(parents=True, exist_ok=True)
    out = dict(revision=a.revision, variant=variant, source_gds=str(src.relative_to(ROOT)))
    with tempfile.TemporaryDirectory() as tmp:
        rebuilt = Path(tmp) / 'rebuilt.gds'
        export_gds(D56.build_design(variant=variant, revision=a.revision), str(rebuilt), fab=True, cheese=True)
        m_reb, lib, top = metal_um(rebuilt, (0, 10, 11, 100))
        (x0, y0), (x1, y1) = top.bounding_box()
        out['units'] = dict(unit_m=lib.unit, ok=abs(lib.unit - 1e-6) < 1e-12)
        out['outline_um'] = dict(bbox=[round(v, 3) for v in (x0, y0, x1, y1)], ok=all(abs(abs(v) - 5000) < 1e-3 for v in (x0, y0, x1, y1)))
        m_src, _, _ = metal_um(src, (0, 10, 11))
        x = m_src.symmetric_difference(fill_small_holes(m_reb, a.hole_um2))
        reg, n = regions(x)
        out['xor_vs_source'] = dict(xor_area_um2=round(x.area, 3), n_regions_ge_1um2=n, metal_area_mm2=round(m_reb.area / 1e6, 4),
                                        ok=x.area <= 1.0, regions=reg)
        out['xor_vs_fab'] = []
        for f in fabs:
            m_fab, _, _ = metal_um(f, (0,))
            m_fab_f = fill_small_holes(m_fab, a.hole_um2)
            m_reb_f = fill_small_holes(m_reb, a.hole_um2)     # compare cheese-free metal on both sides
            x = m_fab_f.symmetric_difference(m_reb_f)
            reg, n = regions(x)
            extra, missing = m_fab_f.difference(m_reb_f), m_reb_f.difference(m_fab_f)
            png = figs / f'xor_fab{len(out["xor_vs_fab"]) + 1}.png'   # short name: MAX_PATH
            plot_xor(m_reb_f, m_fab_f, x, png, f'{a.revision} rebuilt vs {f.name}\nred: metal differs (cheese holes < {a.hole_um2:g} um2 filled)')
            out['xor_vs_fab'].append(dict(fab_gds=str(f.relative_to(ROOT)), xor_area_um2=round(x.area, 2),
                                          metal_only_in_fab_um2=round(extra.area, 2), metal_only_in_rebuilt_um2=round(missing.area, 2),
                                          n_regions_ge_1um2=n, regions=reg, figure=str(png.relative_to(folder))))
        thin = CMW.thin_metal(m_reb, 4.0)
        out['metal_width_4um'] = dict(n_regions=len(thin), note='no component names in the GDS: resonator exemption not applied; review',
                                      regions=[dict(width_um=round(w, 3), bounds_um=[round(v, 1) for v in g.bounds]) for g, w in sorted(thin, key=lambda t: t[1])][:40])
    out['passed_geometry'] = out['units']['ok'] and out['outline_um']['ok'] and out['xor_vs_source']['ok']
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    (folder / 'summaries' / 'gds_check.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    finish_run(folder, status='completed', summary=dict(passed_geometry=out['passed_geometry']))
    print('RUN_FOLDER=' + str(folder))
    print(json.dumps(out, indent=1)[:7000])


if __name__ == '__main__':
    main()
