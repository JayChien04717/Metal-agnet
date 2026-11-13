"""Minimum centre-line spacing between CPW lines of a qiskit-metal design (Design4 rule, 2026-11-12).

Rule: every two CPW runs that lie side by side (parallel, with overlapping projections, not on the
same straight line) must be >= 400 um apart, centre to centre. Not counted:
  * capacitor electrode leads (cin / cp / cj lead_*), they are the capacitor, not a line;
  * collinear pieces of one line (e.g. feed_in | C_in | feed_mid, feed_jf | jf | filter_tail);
  * perpendicular segments (they cannot run side by side; tees / corners);
  * the pairs in ALLOW (the hammer stubs stub_f / stub_r form the coupler C_j on purpose).

Usage:  from cpw_spacing_check import report; report(design, min_um=400)
        python cpw_spacing_check.py         (checks purcell_6500.build() default layout)
"""
import numpy as np
from shapely.geometry import LineString

ALLOW = {frozenset(('stub_f', 'stub_r'))}
SKIP_NAMES = ('lead_a', 'lead_b', 'lead_a_sub', 'lead_b_sub')


def _segments(design):
    names = {c.id: n for n, c in design.components.items()}
    segs = []
    for row in design.qgeometry.tables['path'].itertuples():
        if row.subtract or row.helper or row.name in SKIP_NAMES:
            continue
        pts = np.asarray(row.geometry.coords)
        for k in range(len(pts) - 1):
            a, b = pts[k], pts[k + 1]
            if np.linalg.norm(b - a) < 1e-6:
                continue
            segs.append(dict(a=a, b=b, comp=names[row.component], where=f'{names[row.component]}.{row.name}[{k}]'))
    return segs


def violations(design, min_um=400):
    segs, lim = _segments(design), min_um / 1000
    out = []
    for i, s in enumerate(segs):
        u = (s['b'] - s['a']) / np.linalg.norm(s['b'] - s['a'])
        n = np.array([-u[1], u[0]])
        for t in segs[i + 1:]:
            if frozenset((s['comp'], t['comp'])) in ALLOW:
                continue
            v = (t['b'] - t['a']) / np.linalg.norm(t['b'] - t['a'])
            if abs(np.cross(u, v)) > 0.05:                    # not parallel
                continue
            offset = abs(np.dot(t['a'] - s['a'], n))           # centre-to-centre distance of the two lines
            if offset < 1e-3 or offset >= lim:                 # collinear (same line) or far enough
                continue
            p0, p1 = sorted([np.dot(s['a'], u), np.dot(s['b'], u)])
            q0, q1 = sorted([np.dot(t['a'], u), np.dot(t['b'], u)])
            overlap = min(p1, q1) - max(p0, q0)
            if overlap <= 0:                                   # not side by side
                continue
            out.append((round(1000 * offset, 1), round(1000 * overlap), s['where'], t['where']))
    return sorted(out)


def report(design, min_um=400):
    v = violations(design, min_um)
    if not v:
        print(f'CPW spacing OK: every side-by-side pair of CPW runs >= {min_um} um')
    for d, ov, a, b in v:
        print(f'  {d:6.1f} um (side by side over {ov} um)  {a}  <->  {b}')
    return v


if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    import purcell_6500 as P
    for wide in (False, True):
        print('wide' if wide else 'old folded', 'layout:')
        report(P.build(P.LENGTHS, paper_hammer=True,
                       paper_hammer_gap_um=70, spacing_mm=0.4, wide=wide))
