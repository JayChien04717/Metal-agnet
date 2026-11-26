# -*- coding: utf-8 -*-
r"""GDS checks for a Design11 export (skill/gds_export.md section 3, skill/chip_drc.md C12a).

    python -X utf8 -B Design11_Transmon2Q_Manhattan/layout/check_gds.py [--gds <file>]

Checks (results in results/general/<run>/summaries/gds_check.json, violation plots in figures/):
  1. units um and the 10 x 10 mm die outline (-5000 ... +5000 um);
  2. metal width: Design11 local override 2 um (README), with the general 4 um rule reported for reference;
     resonator parts exempt (check_metal_width.design_exemptions); metal = layer 1 datatypes 0, 10, 11 merged;
  3. marker keep-outs (check_keepouts.check_gds_keepouts, >= 100 um, NTHU / TSRI);
  4. etch survival probe: a representative point of every etched (subtract) qgeometry region of
     build_design() must be etched in the GDS, and of every metal region must be metal.
No GDS is written.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'layout'), str(HERE)]
import project_paths
project_paths.configure(__file__)

import gdstk
from shapely.affinity import scale
from shapely.ops import unary_union

import build as D11
import check_metal_width as CMW
from check_keepouts import check_gds_keepouts
from run_record import create_run, finish_run

DESIGN = 'Design11_Transmon2Q_Manhattan'
DEFAULT_GDS = ROOT / DESIGN / 'gds' / 'Design11_2026-11-25_r9_sim.gds'
MIN_W_OVERRIDE, MIN_W_GENERAL = 2.0, 4.0


def merged_metal(path):
    return unary_union([CMW.gds_metal(path, 1, dt) for dt in (0, 10, 11)])


def probes(design, metal):
    unit = design.parse_value('1um')
    id2n = {c.id: n for n, c in design.components.items()}
    bad, n = [], 0
    for kind in ('poly', 'path'):
        for r in design.qgeometry.tables[kind].itertuples():
            g = r.geometry
            if kind == 'path':
                # qgeometry paths are unfilleted (the renderer adds the fillet), so probe the midpoint of the longest
                # straight segment (furthest from any fillet), offset sideways into the gap (0.75 x half the cut width)
                c = list(g.coords)
                segs = [(c[i], c[i + 1]) for i in range(len(c) - 1)]
                (xa, ya), (xb, yb) = max(segs, key=lambda s: (s[1][0] - s[0][0]) ** 2 + (s[1][1] - s[0][1]) ** 2)
                L = ((xb - xa) ** 2 + (yb - ya) ** 2) ** 0.5
                off = 0.75 * r.width / 2 if bool(r.subtract) else 0.0
                from shapely.geometry import Point
                pt = Point((xa + xb) / 2 - (yb - ya) / L * off, (ya + yb) / 2 + (xb - xa) / L * off)
                pt = scale(pt, 1 / unit, 1 / unit, origin=(0, 0))
                n += 1
                if bool(r.subtract) == metal.contains(pt):
                    bad.append(dict(component=id2n[int(r.component)], name=r.name, subtract=bool(r.subtract),
                                    x_um=round(pt.x, 2), y_um=round(pt.y, 2)))
                continue
            g = scale(g, 1 / unit, 1 / unit, origin=(0, 0))
            if g.is_empty or g.area < 1e-6:
                continue
            if bool(r.subtract):
                # an etch is fully covered by other (non-subtract) metal by design, e.g. CPW cut under its trace:
                # probe the part of the etch not covered by the design's own metal
                g = g.difference(design_metal_um)
                if g.is_empty or g.area < 1e-3:
                    continue
            pt = g.representative_point()
            n += 1
            in_metal = metal.contains(pt)
            if bool(r.subtract) == in_metal:
                bad.append(dict(component=id2n[int(r.component)], name=r.name, subtract=bool(r.subtract),
                                x_um=round(pt.x, 2), y_um=round(pt.y, 2)))
    return n, bad


def main():
    global design_metal_um
    ap = argparse.ArgumentParser()
    ap.add_argument('--gds', default=str(DEFAULT_GDS))
    a = ap.parse_args()
    gds = Path(a.gds)
    folder = create_run(DESIGN, 'general', 'gds_check', entrypoint=__file__, inputs=[gds, HERE / 'build.py'],
                        parameters=dict(gds=str(gds), min_width_override_um=MIN_W_OVERRIDE, min_width_general_um=MIN_W_GENERAL))
    out = dict(gds=gds.name)
    lib = gdstk.read_gds(str(gds))
    (x0, y0), (x1, y1) = lib.top_level()[0].bounding_box()
    out['units'] = dict(unit_m=lib.unit, ok=abs(lib.unit - 1e-6) < 1e-12)
    out['outline_um'] = dict(bbox=[x0, y0, x1, y1], ok=all(abs(abs(v) - 5000) < 1e-3 for v in (x0, y0, x1, y1)))

    design = D11.build_design(variant='sim', channels=True, check_keepout_um=None)
    unit = design.parse_value('1um')
    design_metal_um = unary_union(
        [scale(r.geometry if k == 'poly' else r.geometry.buffer(r.width / 2, cap_style=2), 1 / unit, 1 / unit, origin=(0, 0))
         for k in ('poly', 'path') for r in design.qgeometry.tables[k].itertuples() if not bool(r.subtract)])
    metal = merged_metal(gds)
    exempt = CMW.design_exemptions(design, MIN_W_GENERAL)
    (folder / 'figures').mkdir(parents=True, exist_ok=True)
    for label, mw in (('override_2um', MIN_W_OVERRIDE), ('general_4um', MIN_W_GENERAL)):
        found = [(g, w) for g, w in CMW.thin_metal(metal, mw) if exempt is None or not g.within(exempt)]
        out[f'metal_width_{label}'] = dict(
            min_width_um=mw, n_violations=len(found), ok=not found,
            violations=[dict(width_um=round(w, 3), bounds_um=[round(v, 1) for v in g.bounds]) for g, w in sorted(found, key=lambda t: t[1])][:40])
    try:
        out['keepouts'] = check_gds_keepouts(str(gds), min_um=100.0, sites=('NTHU', 'TSRI'), raise_on_violation=False)
    except Exception as e:  # no site markers in a sim GDS
        out['keepouts'] = dict(passed=None, note=f'not evaluated: {e}')
    n, bad = probes(design, metal)
    out['etch_survival'] = dict(n_probes=n, n_failed=len(bad), ok=not bad, failed=bad[:40])
    out['passed'] = all(out[k]['ok'] for k in ('units', 'outline_um', 'metal_width_override_2um', 'etch_survival')) \
        and out['keepouts'].get('passed') is not False
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    (folder / 'summaries' / 'gds_check.json').write_text(json.dumps(out, indent=2, default=str), encoding='utf-8')
    finish_run(folder, status='completed', summary=dict(passed=out['passed']))
    print('RUN_FOLDER=' + str(folder))
    print(json.dumps({k: (v if k != 'keepouts' else {kk: v.get(kk) for kk in ('passed', 'min_clearance_um', 'note') if kk in v})
                      for k, v in out.items()}, indent=1, default=str)[:6000])


if __name__ == '__main__':
    design_metal_um = None
    main()
