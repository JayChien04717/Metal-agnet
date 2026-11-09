"""Plot a qiskit-metal design straight from its qgeometry tables (no GDS export): metal = dark, etched gap = light."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from shapely.geometry import MultiPolygon, Polygon


def _polys(g):
    if isinstance(g, Polygon):
        return [g]
    if isinstance(g, MultiPolygon):
        return list(g.geoms)
    return [p for p in getattr(g, 'geoms', []) if isinstance(p, Polygon)]


def _fillet(line, R, n=12):
    """Round the interior vertices of a polyline with arcs of radius R (limited to half the adjacent segments)."""
    from shapely.geometry import LineString
    P = np.array(line.coords)
    if len(P) < 3:
        return line
    out = [P[0]]
    for i in range(1, len(P) - 1):
        a, b, c = P[i - 1], P[i], P[i + 1]
        u, v = (a - b), (c - b)
        lu, lv = np.hypot(*u), np.hypot(*v)
        u, v = u / lu, v / lv
        cosang = np.clip(np.dot(u, v), -1, 1)
        th = np.arccos(cosang)
        if th > np.pi - 1e-6 or th < 1e-6:
            out.append(b); continue
        t = min(R / np.tan(th / 2), lu / 2, lv / 2)
        r = t * np.tan(th / 2)
        p1, p2 = b + u * t, b + v * t
        bis = (u + v); bis /= np.hypot(*bis)
        cen = b + bis * (r / np.sin(th / 2))
        a1 = np.arctan2(*(p1 - cen)[::-1]); a2 = np.arctan2(*(p2 - cen)[::-1])
        d = (a2 - a1 + np.pi) % (2 * np.pi) - np.pi
        for s in np.linspace(0, 1, n):
            ang = a1 + d * s
            out.append(cen + r * np.array([np.cos(ang), np.sin(ang)]))
    out.append(P[-1])
    return LineString(out)


def plot_design(design, png, window, comps=None, dpi=130, labels=True):
    metal, gap = [], []
    ids = {c.id: n for n, c in design.components.items()}
    for tname in ('path', 'poly'):
        t = design.qgeometry.tables[tname]
        for _, r in t.iterrows():
            if comps is not None and ids.get(r['component']) not in comps:
                continue
            g = r.geometry
            if tname == "path":
                if r.get("fillet", 0):
                    g = _fillet(g, r["fillet"])
                g = g.buffer(r['width'] / 2, cap_style=2)
            (gap if r['subtract'] else metal).extend(_polys(g))
    fig, ax = plt.subplots(figsize=(9, 9 * (window[3] - window[2]) / (window[1] - window[0])))
    ax.set_facecolor('#c9ced6')                       # ground plane
    ax.add_collection(PolyCollection([np.array(p.exterior.coords) for p in gap], facecolor='white', edgecolor='none'))
    ax.add_collection(PolyCollection([np.array(p.exterior.coords) for p in metal], facecolor='#33415c', edgecolor='none'))
    if labels:
        for n, c in design.components.items():
            try:
                pts = [p['middle'] for p in c.pins.values()]
            except Exception:
                continue
    ax.set_xlim(window[0], window[1]); ax.set_ylim(window[2], window[3]); ax.set_aspect('equal')
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    fig.tight_layout(); fig.savefig(png, dpi=dpi); plt.close(fig)
