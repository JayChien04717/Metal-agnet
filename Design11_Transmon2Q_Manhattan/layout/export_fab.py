# -*- coding: utf-8 -*-
r"""Stage 3 single-die export for Design11 (NCKU): matching *_sim.gds and *_fab.gds of the current revision.

NCKU delivery (2026-11-27): single 10 mm die; e-beam site NCKU (no alignment markers); junction sites left
EMPTY (no drawn e-beam junction layers; NCKU draws them); no test-junction coupons; die ID 'IPRT2Q00'.

Usage (repository root, qmetal):
    python -X utf8 -B Design11_Transmon2Q_Manhattan/layout/export_fab.py
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'code' / 'general'))
import project_paths
project_paths.configure(__file__)
DESIGN = project_paths.design_dir(__file__)
LAYOUT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(LAYOUT_DIR))
sys.path.insert(0, str(ROOT / 'Design3_FTF_NTHU' / 'layout' / 'components'))

import gdstk
import numpy as np
import shapely
from shapely.affinity import translate

import build as B
import check_metal_width as cmw
from gds_tools import export_gds, write_wafer_gds
from round_mask import round_mask
from run_record import create_run, finish_run
from wafer_cheese_grid import replace_wafer_cheese
from wafer_pcbport import wafer_pcb_ports
from xmark import LABEL_SIZE_UM, label_geometry

DIE_ID = 'IPRT2Q00'
DATE = '2026-11-27'
STEM = f'Design11_{DATE}_{B.REVISION}'
SIM_GDS = DESIGN / 'gds' / f'{STEM}_sim.gds'
FAB_GDS = DESIGN / 'gds' / f'{STEM}_fab.gds'


def layer1(path):
    lib = gdstk.read_gds(str(path))
    assert lib.unit == 1e-6 and lib.precision == 1e-9, 'expected um units, 1 nm precision'
    top, = lib.top_level()
    return top, [p for p in top.get_polygons() if p.layer == 1]


def covered(polys, pts):
    return [bool(v) for v in gdstk.inside([tuple(map(float, q)) for q in pts], polys)]


def annotate_junctions(jt, path):
    lib = gdstk.read_gds(str(path))
    top, = lib.top_level()
    for row in jt.itertuples():
        g = row.geometry.buffer(row.width / 2.0, cap_style=2)
        top.add(gdstk.Polygon([(x * 1e3, y * 1e3) for x, y in g.exterior.coords[:-1]], layer=50, datatype=0))
    lib.write_gds(str(path), max_points=199)


def etch_die_id(path, text=DIE_ID, inset_um=50.0, height_um=LABEL_SIZE_UM):
    lib = gdstk.read_gds(str(path))
    top, = lib.top_level()
    glyph = label_geometry(text, height_um, 10, 400, 30)
    x0, y0, x1, y1 = glyph.bounds
    placed = translate(glyph, -5000.0 + inset_um - x0, 5000.0 - inset_um - y1)
    cuts = [gdstk.Polygon(list(p.exterior.coords)[:-1], layer=1, datatype=0)
            for p in shapely.get_parts(placed) if p.area > 0]
    metal = [p for p in top.get_polygons() if p.layer == 1]
    out = gdstk.Library(unit=1e-6, precision=1e-9)
    cell = out.new_cell('TOP')
    cell.add(*gdstk.boolean(metal, cuts, 'not', precision=0.001, layer=1, datatype=0))
    out.write_gds(str(path), max_points=199)
    return placed


def route_probes(design, polys):
    t = design.qgeometry.tables['path']
    names = {c.id: n for n, c in design.components.items()}
    probes, labels = [], []
    for comp_id in t.component.unique():
        rows = t[t.component == comp_id]
        sig, cut = rows[~rows['subtract']], rows[rows['subtract']]
        for s_row in sig.itertuples():
            match = cut[cut.name == s_row.name + '_cut']
            c_rows = match if not match.empty else cut
            if c_rows.empty:
                continue
            c_row = c_rows.iloc[0]
            pts = np.asarray(s_row.geometry.coords) * 1e3
            seg = np.diff(pts, axis=0)
            L = np.linalg.norm(seg, axis=1)
            i = int(np.argmax(L))
            if L[i] < 30:
                continue
            mid = (pts[i] + pts[i + 1]) / 2
            u = seg[i] / L[i]
            nrm = np.array([-u[1], u[0]])
            gap_c = 250 * (s_row.width + c_row.width)
            probes += [mid, mid + gap_c * nrm, mid - gap_c * nrm]
            labels.append(f'{names[comp_id]}:{s_row.name}')
    inside = covered(polys, probes)
    bad = [lab for k, lab in enumerate(labels)
           if not (inside[3 * k] and not inside[3 * k + 1] and not inside[3 * k + 2])]
    return len(labels), bad


def main():
    folder = create_run(DESIGN.name, 'general', f'gds_export_{B.REVISION}_ncku', entrypoint=__file__,
                        inputs=[DESIGN / 'design.json', LAYOUT_DIR / 'build.py'],
                        parameters=dict(die_id=DIE_ID, site='NCKU (no markers)', junction_layers='none (sites empty)',
                                        test_junctions='2 x ManhattanTestJunction at (+-1.50, -3.10) mm', revision=B.REVISION, round_radius_um=5.0,
                                        cheese='wafer_cheese_grid 5 um / 40 um, launch clearance profile'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    checks = {}
    try:
        SIM_GDS.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(dir=SIM_GDS.parent) as tmp:
            tmp = Path(tmp)
            d_sim = B.build_design(variant='sim', check_keepout_um=100.0)
            jt = d_sim.qgeometry.tables['junction'].copy()
            export_gds(d_sim, str(tmp / 'sim_raw.gds'), fab=True, cheese=False)
            write_wafer_gds(tmp / 'sim_raw.gds', SIM_GDS)
            annotate_junctions(jt, SIM_GDS)

            d_fab = B.build_design(variant='fab', check_keepout_um=100.0)
            checks['keepouts'] = 'build_design check_keepouts (min 100 um, NTHU + TSRI sites) passed'
            export_gds(d_fab, str(tmp / 'fab_raw.gds'), fab=True, cheese=False)
            round_mask(tmp / 'fab_raw.gds', tmp / 'fab_round.gds', radius_um=5.0)
            stage = tmp / FAB_GDS.name
            write_wafer_gds(tmp / 'fab_round.gds', stage)

            glyph = label_geometry(DIE_ID, LABEL_SIZE_UM, 10, 400, 30)
            gx0, gy0, gx1, gy1 = glyph.bounds
            w = (-5000.0 + 50 - 30, 5000.0 - 50 - (gy1 - gy0) - 30,
                 -5000.0 + 50 + (gx1 - gx0) + 30, 5000.0 - 50 + 30)
            _, polys = layer1(stage)
            win = gdstk.rectangle(w[:2], w[2:])
            metal_in = gdstk.boolean(polys, [win], 'and', precision=0.001)
            etched = win.area() - sum(q.area() for q in metal_in)
            die = etch_die_id(stage)
            checks['die_id'] = dict(text=DIE_ID, bounds_um=[round(v, 1) for v in die.bounds],
                                    other_etch_in_window_um2=round(etched, 3))
            assert etched < 1.0, 'die ID window already contains etch'

            for mw in (4.0, 2.0):
                thin = cmw.check(stage, min_width=mw, exempt=cmw.design_exemptions(d_fab, min_width=mw))
                checks[f'metal_width_{mw:g}um'] = dict(
                    violations=len(thin), checked_on='rounded mask with die ID, before cheese',
                    regions=[dict(width_um=round(wd, 2), bounds=[round(v) for v in r.bounds]) for r, wd in thin[:20]])
            cheese = replace_wafer_cheese(stage, launches=wafer_pcb_ports, rounded=True)
            checks['cheese'] = {k: cheese[k] for k in ('new_holes', 'pitch_um', 'hole_um', 'edge_um') if k in cheese}
            checks['cheese']['islands'] = len(cheese.get('islands', []))
            stage.replace(FAB_GDS)

        for tag, path, design in (('sim', SIM_GDS, d_sim), ('fab', FAB_GDS, d_fab)):
            top, polys = layer1(path)
            bb = top.bounding_box()
            assert bb == ((-5000.0, -5000.0), (5000.0, 5000.0)), (tag, bb)
            n, bad = route_probes(design, polys)
            checks[f'{tag}_outline'] = 'um units, 10 x 10 mm (-5000 ... 5000 um)'
            checks[f'{tag}_etch_survival'] = dict(routes_probed=n, failed=bad)
            assert not bad, (tag, bad)
        top, fab_polys = layer1(FAB_GDS)
        names = {c.id: n for n, c in d_sim.components.items()}
        sites = [(f'{names[r.component]}:{r.name}', r.geometry.interpolate(.5, normalized=True))
                 for r in jt.itertuples()]
        hit = covered(fab_polys, [(pt.x * 1e3, pt.y * 1e3) for _, pt in sites])
        filled = [lab for (lab, _), h in zip(sites, hit) if h]
        checks['junction_sites_empty'] = dict(sites=len(sites), metal_at_site=filled)
        assert not filled, filled
        checks['fab_layers'] = sorted({(q.layer, q.datatype) for q in top.get_polygons()})
        summary = dict(sim_gds=str(SIM_GDS.relative_to(ROOT)), sim_sha256=sha256(SIM_GDS.read_bytes()).hexdigest(),
                       fab_gds=str(FAB_GDS.relative_to(ROOT)), fab_sha256=sha256(FAB_GDS.read_bytes()).hexdigest(),
                       checks=checks)
        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(json.dumps(summary, indent=2, default=str), encoding='utf-8')
        finish_run(folder, status='completed', summary=summary)
        print(json.dumps(summary, indent=2, default=str))
    except Exception as exc:
        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(json.dumps(dict(checks=checks, error=repr(exc)), indent=2,
                                                                      default=str), encoding='utf-8')
        finish_run(folder, status='failed', summary=dict(checks=checks, error=repr(exc)))
        raise


if __name__ == '__main__':
    main()
