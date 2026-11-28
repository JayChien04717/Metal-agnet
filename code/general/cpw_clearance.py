# -*- coding: utf-8 -*-
"""CPW centerline clearance / crossing checker (skill/airbridge.md §1, §3, §4a).

All coordinates in mm. A "net" is a list of polylines (several sub-routes split at
bridges belong to the same net) plus a priority (0 = readout / Purcell, >= 1 = control).

Checks
------
- clearance_violations: true minimum centerline distance between distinct nets
  (< d_min), excluding a disc of radius d_min / sin(theta) around each legitimate
  (>= min_cross_deg) crossing.
- crossing_report: every transverse crossing with its angle; oblique ones flagged.
- bridge_spacing_violations: consecutive crossings along one control net closer than
  min_gap (arc length).
"""
from __future__ import annotations

import math

from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import unary_union


def _segments(poly):
    for k in range(len(poly) - 1):
        a, b = poly[k], poly[k + 1]
        if math.hypot(b[0] - a[0], b[1] - a[1]) > 1e-6:
            yield a, b


def crossings(net_a, net_b, min_angle_parallel=15.0):
    """Return [(x, y, angle_deg)] of transverse crossings between two nets (lists of polylines)."""
    out = []
    for pa in net_a:
        for a0, a1 in _segments(pa):
            sa = LineString([a0, a1])
            ang_a = math.atan2(a1[1] - a0[1], a1[0] - a0[0])
            for pb in net_b:
                for b0, b1 in _segments(pb):
                    sb = LineString([b0, b1])
                    if not sa.intersects(sb):
                        continue
                    inter = sa.intersection(sb)
                    if inter.geom_type != 'Point':
                        continue
                    ang_b = math.atan2(b1[1] - b0[1], b1[0] - b0[0])
                    d = abs((math.degrees(ang_a - ang_b) + 90.0) % 180.0 - 90.0)
                    if d < min_angle_parallel:
                        continue
                    out.append((float(inter.x), float(inter.y), d))
    # de-duplicate (crossing exactly at a vertex is reported by two segments)
    ded = []
    for c in out:
        if not any(math.hypot(c[0] - e[0], c[1] - e[1]) < 1e-4 for e in ded):
            ded.append(c)
    return ded


def _bridge_pairs(nets, geoms, bridges, pin_tol=0.045, under_tol=0.006):
    """Map each crossover bridge centre to (overpass net, underpass net): the overpass net has a sub-route
    end at the bridge pins, the underpass net passes (continuous) through the centre."""
    out = []
    for bx, by in bridges:
        c = Point(bx, by)
        over = [n for n in nets if any(math.hypot(p[i][0] - bx, p[i][1] - by) < pin_tol
                                       for p in nets[n]['polys'] for i in (0, -1))]
        under = [n for n in nets if n not in over and geoms[n].distance(c) < under_tol]
        for a in over:
            for b in under:
                out.append((bx, by, a, b))
    return out


def clearance_violations(nets, d_min=0.18, min_cross_deg=85.0, skip_pairs=None, connected=None, bridges=None):
    """nets: {name: dict(polys=[...], priority=int)}. Returns list of
    (dist_um, name_a, name_b, x, y). Pairs where both are priority 0 are skipped
    (resonator-to-own-Purcell coupling is by design); pairs in skip_pairs too.
    connected: {name: set(names)} nets that share a pin (touching is expected near that pin).
    bridges: [(x, y)] crossover-bridge centres; the overpass is drawn split there, so the crossing is
    recovered from the bridge (perpendicular by DRC) and exempted like any legitimate crossing."""
    skip_pairs = skip_pairs or set()
    connected = connected or {}
    names = list(nets)
    geoms = {n: MultiLineString([LineString(p) for p in nets[n]['polys'] if len(p) >= 2]) for n in names}
    bpairs = _bridge_pairs(nets, geoms, bridges or [])
    bad = []
    for i, na in enumerate(names):
        ga = geoms[na]
        for nb in names[i + 1:]:
            if nets[na]['priority'] == 0 and nets[nb]['priority'] == 0:
                continue
            if frozenset((na, nb)) in skip_pairs:
                continue
            gb = geoms[nb]
            if ga.distance(gb) >= d_min:
                continue
            # exclusion discs around legitimate crossings and shared pins
            excl = []
            for x, y, ang in crossings(nets[na]['polys'], nets[nb]['polys']):
                if ang >= min_cross_deg:
                    excl.append(Point(x, y).buffer(d_min / max(math.sin(math.radians(ang)), 0.5) + 0.01))
            for bx, by, o, u in bpairs:
                if {o, u} == {na, nb}:
                    excl.append(Point(bx, by).buffer(d_min + 0.01))
            if nb in connected.get(na, set()) or na in connected.get(nb, set()):
                for p in nets[na]['polys']:
                    for e in (p[0], p[-1]):
                        if gb.distance(Point(e)) < 0.02:
                            excl.append(Point(e).buffer(d_min + 0.02))
            ex = unary_union(excl) if excl else None
            ra = ga.difference(ex) if ex is not None else ga
            rb = gb.difference(ex) if ex is not None else gb
            if ra.is_empty or rb.is_empty:
                continue
            d = ra.distance(rb)
            if d < d_min - 1e-6:
                # locate: nearest point on ra to rb
                from shapely.ops import nearest_points
                p1, _ = nearest_points(ra, rb)
                bad.append((round(d * 1000.0, 1), na, nb, round(p1.x, 3), round(p1.y, 3)))
    bad.sort()
    return bad


def _arc_position(poly_list, x, y):
    """Arc length of point (x, y) along the concatenated net polylines (closest)."""
    s_acc = 0.0
    best = (1e9, 0.0)
    for poly in poly_list:
        ls = LineString(poly)
        d = ls.distance(Point(x, y))
        if d < best[0]:
            best = (d, s_acc + ls.project(Point(x, y)))
        s_acc += ls.length
    return best[1]


def crossing_report(nets, min_cross_deg=85.0):
    """Return (all_crossings, oblique) with entries (name_a, name_b, x, y, angle)."""
    names = list(nets)
    allc, obl = [], []
    for i, na in enumerate(names):
        for nb in names[i + 1:]:
            if nets[na]['priority'] == 0 and nets[nb]['priority'] == 0:
                continue
            for x, y, ang in crossings(nets[na]['polys'], nets[nb]['polys']):
                rec = (na, nb, round(x, 3), round(y, 3), round(ang, 1))
                allc.append(rec)
                if ang < min_cross_deg:
                    obl.append(rec)
    return allc, obl


def bridge_spacing_violations(nets, all_crossings, min_gap=0.50):
    """Consecutive crossings along each control net closer than min_gap (mm arc length)."""
    per = {}
    for na, nb, x, y, _ in all_crossings:
        for n in (na, nb):
            if nets[n]['priority'] >= 1:
                per.setdefault(n, []).append((_arc_position(nets[n]['polys'], x, y), x, y))
    bad = []
    for n, lst in per.items():
        lst.sort()
        for k in range(len(lst) - 1):
            gap = lst[k + 1][0] - lst[k][0]
            if gap < min_gap:
                bad.append((n, round(gap * 1000.0, 0), lst[k][1], lst[k][2]))
    return bad
