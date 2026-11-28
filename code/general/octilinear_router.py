# -*- coding: utf-8 -*-
"""Clearance-aware octilinear grid (A*) CPW router (skill/airbridge.md §4a).

Coordinates are in mm. The chip is rasterized on a square grid of pitch `h`.
Each cell carries a bit mask:

- bits 0..3  : an already-routed CPW of orientation class o (0: 0 deg, 1: 45 deg,
               2: 90 deg, 3: 135 deg) passes within the clearance radius.
- bit 4      : HARD obstacle (component pocket / pad / resonator meander hull / die margin).
- bit 5      : SOFT zone (inside the preferred clearance but outside the hard radius).

A move in direction d (0..7, multiples of 45 deg) may enter a cell holding
orientation bits only when that cell's bits are exactly the class perpendicular
to d, and no turn is allowed while inside such a buffer. Therefore a new route
can only cross an existing CPW straight and at 90 deg, never runs parallel to it
closer than the clearance radius, and never crosses near its corners (where two
orientation classes overlap).

Every obstacle / line is stored as an item with an `owner`; when routing net N
the items owned by N are ignored (its own reserved stubs, its own fixed pieces).
"""
from __future__ import annotations

import heapq
import math

import numpy as np
import shapely
from shapely.geometry import LineString

DX = (1, 1, 0, -1, -1, -1, 0, 1)
DY = (0, 1, 1, 1, 0, -1, -1, -1)
HARD = 1 << 4
SOFT = 1 << 5
ORI_MASK = 0x0F
SQRT2 = math.sqrt(2.0)


def dir_index(vec, tol_deg=3.0):
    """Octilinear direction index (0..7) of a vector, or None if not octilinear."""
    ang = math.degrees(math.atan2(vec[1], vec[0]))
    k = round(ang / 45.0)
    if abs(ang - 45.0 * k) > tol_deg:
        return None
    return int(k) % 8


class RouteFailed(RuntimeError):
    pass


class OctiGrid:
    def __init__(self, xmin, xmax, ymin, ymax, h=0.05):
        self.h = float(h)
        self.x0, self.y0 = float(xmin), float(ymin)
        self.nx = int(math.ceil((xmax - xmin) / h)) + 1
        self.ny = int(math.ceil((ymax - ymin) / h)) + 1
        self.items = []  # (owner, flat_idx ndarray, bits ndarray|int)
        xs = self.x0 + self.h * np.arange(self.nx)
        ys = self.y0 + self.h * np.arange(self.ny)
        self._X, self._Y = np.meshgrid(xs, ys)  # [iy, ix]

    # ------------------------------------------------------------------ geometry
    def xy(self, ix, iy):
        return (self.x0 + ix * self.h, self.y0 + iy * self.h)

    def cell_ahead(self, x, y, d):
        """Grid cell at or ahead of (x, y) along direction d, so a bend placed there never lands behind
        the continuous start point after exact re-intersection."""
        def _c(v, v0, comp, axis_aligned):
            f = (v - v0) / self.h
            if comp > 0:
                return int(math.ceil(f + (0.5 if axis_aligned else 0.0) - 1e-9))
            if comp < 0:
                return int(math.floor(f - (0.5 if axis_aligned else 0.0) + 1e-9))
            return int(round(f))
        ax = (d % 2 == 0)
        return (_c(x, self.x0, DX[d], ax), _c(y, self.y0, DY[d], ax))

    def cell(self, x, y):
        return (int(round((x - self.x0) / self.h)), int(round((y - self.y0) / self.h)))

    def _window(self, bounds, pad):
        x_lo, y_lo, x_hi, y_hi = bounds
        ix0 = max(0, int(math.floor((x_lo - pad - self.x0) / self.h)))
        ix1 = min(self.nx - 1, int(math.ceil((x_hi + pad - self.x0) / self.h)))
        iy0 = max(0, int(math.floor((y_lo - pad - self.y0) / self.h)))
        iy1 = min(self.ny - 1, int(math.ceil((y_hi + pad - self.y0) / self.h)))
        return ix0, ix1, iy0, iy1

    def add_hard(self, geom, owner=None):
        """Mark every cell whose centre lies inside shapely `geom` as HARD."""
        if geom.is_empty:
            return
        ix0, ix1, iy0, iy1 = self._window(geom.bounds, self.h)
        if ix1 < ix0 or iy1 < iy0:
            return
        X = self._X[iy0:iy1 + 1, ix0:ix1 + 1]
        Y = self._Y[iy0:iy1 + 1, ix0:ix1 + 1]
        inside = shapely.contains_xy(geom, X, Y) | shapely.intersects_xy(geom.boundary, X, Y)
        iy, ix = np.nonzero(inside)
        flat = (iy + iy0) * self.nx + (ix + ix0)
        self.items.append((owner, flat.astype(np.int64), HARD))

    def add_hard_rect_outside(self, x_lo, x_hi, y_lo, y_hi, owner=None):
        """Mark everything outside the rectangle as HARD (die margin)."""
        out = (self._X < x_lo) | (self._X > x_hi) | (self._Y < y_lo) | (self._Y > y_hi)
        iy, ix = np.nonzero(out)
        self.items.append((owner, (iy * self.nx + ix).astype(np.int64), HARD))

    def add_line(self, poly, owner, radius, soft_radius=None):
        """Oriented clearance buffer around a polyline (non-octilinear segments become HARD)."""
        soft_radius = soft_radius or radius
        for k in range(len(poly) - 1):
            a, b = np.asarray(poly[k], float), np.asarray(poly[k + 1], float)
            L = float(np.linalg.norm(b - a))
            if L < 1e-6:
                continue
            d = dir_index(b - a)
            bounds = (min(a[0], b[0]), min(a[1], b[1]), max(a[0], b[0]), max(a[1], b[1]))
            ix0, ix1, iy0, iy1 = self._window(bounds, soft_radius + self.h)
            X = self._X[iy0:iy1 + 1, ix0:ix1 + 1]
            Y = self._Y[iy0:iy1 + 1, ix0:ix1 + 1]
            u = (b - a) / L
            t = np.clip((X - a[0]) * u[0] + (Y - a[1]) * u[1], 0.0, L)
            dist = np.hypot(X - (a[0] + t * u[0]), Y - (a[1] + t * u[1]))
            bit = HARD if d is None else (1 << (d % 4))
            iy, ix = np.nonzero(dist < radius)
            self.items.append((owner, ((iy + iy0) * self.nx + (ix + ix0)).astype(np.int64), bit))
            if soft_radius > radius:
                iy, ix = np.nonzero((dist >= radius) & (dist < soft_radius))
                self.items.append((owner, ((iy + iy0) * self.nx + (ix + ix0)).astype(np.int64), SOFT))

    def remove_owner(self, owner):
        self.items = [it for it in self.items if it[0] != owner]

    def mask_for(self, ignore_owners=()):
        ign = set(ignore_owners)
        m = np.zeros(self.nx * self.ny, np.uint8)
        for owner, flat, bits in self.items:
            if owner in ign:
                continue
            if isinstance(bits, int):
                m[flat] |= np.uint8(bits)
            else:
                np.bitwise_or.at(m, flat, bits)
        return m

    # ------------------------------------------------------------------ routing
    def _march_free(self, mask, p, d, max_mm=2.0, lookahead=3):
        """Step from continuous point p along direction d until the snapped cell (and the next two) are free.

        HARD cells are only tolerated while still inside the pin's own component zone; after that the stub
        may only cross other lines perpendicularly (otherwise RouteFailed -> caller rips up the blocker)."""
        step = self.h * (SQRT2 if d % 2 else 1.0)
        ux, uy = DX[d] / (SQRT2 if d % 2 else 1.0), DY[d] / (SQRT2 if d % 2 else 1.0)
        perp_bit = 1 << ((d + 2) % 4)
        left_hard = False
        n = 0
        while n * step <= max_mm:
            q = (p[0] + ux * n * step, p[1] + uy * n * step)
            ix, iy = self.cell(*q)
            if not (0 <= ix < self.nx and 0 <= iy < self.ny):
                break
            mm = mask[iy * self.nx + ix]
            if mm & HARD:
                if left_hard:
                    raise RouteFailed(f'stub from {p} dir {d} runs into a hard obstacle at {q}')
                n += 1
                continue
            left_hard = True
            ob = mm & ORI_MASK
            if ob and ob != perp_bit:
                raise RouteFailed(f'stub from {p} dir {d} would run non-perpendicular through a line at {q}')
            ok = True
            for j in range(lookahead):  # this cell and the next (lookahead - 1) along d must be free
                jx, jy = ix + j * DX[d], iy + j * DY[d]
                if not (0 <= jx < self.nx and 0 <= jy < self.ny):
                    ok = False
                    break
                mj = mask[jy * self.nx + jx]
                if (mj & HARD) or (mj & ORI_MASK):
                    ok = False
                    break
            if ok:
                return q, (ix, iy)
            n += 1
        raise RouteFailed(f'no free cell within {max_mm} mm from {p} along dir {d}')

    def route(self, owner, p_start, d_start, p_end, d_end, ignore_owners=(), *, start_is_anchor=False,
              turn45=0.12, turn90=0.60, cross_pen=1.2, soft_pen=0.6,
              min_run45=2, min_run90=3, max_expand=4_000_000):
        """Route from pin p_start (leaving along d_start) to pin p_end (arriving along d_end).

        Returns the exact octilinear polyline [p_start, ..., p_end] (mm)."""
        mask_np = self.mask_for(set(ignore_owners) | {owner})
        mask = bytes(mask_np)
        nx, ny = self.nx, self.ny

        # Fixed straight stubs out of the start pin and into the end pin
        if start_is_anchor:
            # route continues from a fixed anchor (not a component pin): may turn right away
            q_s = tuple(p_start)
        else:
            try:
                q_s, _c = self._march_free(mask, p_start, d_start)
            except RouteFailed:
                # crowded pin: accept the first free cell and let A* turn away earlier
                try:
                    q_s, _c = self._march_free(mask, p_start, d_start, lookahead=1)
                except RouteFailed as exc:
                    raise RouteFailed(f'{owner}: start {exc}')
        sx, sy = self.cell_ahead(q_s[0], q_s[1], d_start)
        m0 = mask[sy * nx + sx] if (0 <= sx < nx and 0 <= sy < ny) else HARD
        if (m0 & HARD) or ((m0 & ORI_MASK) and (m0 & ORI_MASK) != (1 << ((d_start + 2) % 4))):
            raise RouteFailed(f'{owner}: start cell at {self.xy(sx, sy)} is blocked')
        d_back = (d_end + 4) % 8
        try:
            q_e, (ex, ey) = self._march_free(mask, p_end, d_back)
        except RouteFailed as exc:
            raise RouteFailed(f'{owner}: end {exc}')
        start_c = sy * nx + sx
        goal_c = ey * nx + ex
        h = self.h

        def heur(c):
            dx = abs(c % nx - ex)
            dy = abs(c // nx - ey)
            return h * (max(dx, dy) + (SQRT2 - 1.0) * min(dx, dy))

        RUNCAP = 3
        def key(c, d, r):
            return (c * 8 + d) * 4 + r

        k0 = key(start_c, d_start, RUNCAP)  # start cell is snapped ahead of the start point -> may turn here
        gbest = {k0: 0.0}
        parent = {k0: None}
        heap = [(heur(start_c), 0.0, k0)]
        n_exp = 0
        found = None
        while heap:
            f, g, k = heapq.heappop(heap)
            if g > gbest.get(k, 1e18) + 1e-12:
                continue
            r = k % 4
            d = (k // 4) % 8
            c = k // 32
            if c == goal_c and d == d_end:
                found = k
                break
            n_exp += 1
            if n_exp > max_expand:
                break
            cx, cy = c % nx, c // nx
            m_here = mask[c]
            in_buf = m_here & ORI_MASK
            for dd, tc, need in ((0, 0.0, 0), (1, turn45, min_run45), (-1, turn45, min_run45),
                                 (2, turn90, min_run90), (-2, turn90, min_run90)):
                if dd != 0 and (in_buf or r < need):
                    continue
                nd = (d + dd) % 8
                nxi, nyi = cx + DX[nd], cy + DY[nd]
                if nxi < 0 or nyi < 0 or nxi >= nx or nyi >= ny:
                    continue
                nc = nyi * nx + nxi
                mm = mask[nc]
                if mm & HARD:
                    continue
                ob = mm & ORI_MASK
                step = h * (SQRT2 if nd % 2 else 1.0)
                cost = step + tc
                if ob:
                    if ob != (1 << ((nd + 2) % 4)):
                        continue
                    if not in_buf:
                        cost += cross_pen
                if mm & SOFT:
                    cost += soft_pen * step
                nr = 0 if dd != 0 else min(RUNCAP, r + 1)
                nk = key(nc, nd, nr)
                ng = g + cost
                if ng < gbest.get(nk, 1e18) - 1e-12:
                    gbest[nk] = ng
                    parent[nk] = k
                    heapq.heappush(heap, (ng + heur(nc), ng, nk))
        if found is None:
            raise RouteFailed(f'{owner}: no path (expanded {n_exp})')

        # Reconstruct cell/direction chain
        chain = []
        k = found
        while k is not None:
            chain.append(((k // 32), (k // 4) % 8))
            k = parent[k]
        chain.reverse()
        # Lines: (point, dir) for each straight run
        runs = []
        for c, d in chain:
            if not runs or runs[-1][1] != d:
                runs.append((self.xy(c % nx, c // nx), d))
        # runs[0] direction == d_start; force its line through p_start, last through p_end
        runs[0] = (tuple(p_start), runs[0][1])
        if len(runs) == 1:
            # single straight run: check alignment
            d0 = runs[0][1]
            ux, uy = DX[d0], DY[d0]
            off = abs((p_end[0] - p_start[0]) * uy - (p_end[1] - p_start[1]) * ux) / math.hypot(ux, uy)
            if off > 1e-6:
                raise RouteFailed(f'{owner}: straight path misaligned by {off*1000:.1f} um')
            return [tuple(p_start), tuple(p_end)]
        runs[-1] = (tuple(p_end), runs[-1][1])
        verts = [tuple(p_start)]
        for k in range(1, len(runs)):
            (pa, da), (pb, db) = runs[k - 1], runs[k]
            verts.append(_intersect(pa, da, pb, db))
        verts.append(tuple(p_end))
        return verts


def _intersect(pa, da, pb, db):
    ax, ay = DX[da], DY[da]
    bx, by = DX[db], DY[db]
    den = ax * by - ay * bx
    if abs(den) < 1e-12:
        raise RouteFailed('parallel consecutive runs')
    t = ((pb[0] - pa[0]) * by - (pb[1] - pa[1]) * bx) / den
    return (pa[0] + t * ax, pa[1] + t * ay)


def min_segment_mm(poly):
    return min(math.hypot(poly[k + 1][0] - poly[k][0], poly[k + 1][1] - poly[k][1]) for k in range(len(poly) - 1))


def polyline_ok(poly, min_seg=0.03):
    """Sanity: positive-length segments, all octilinear, no turn > 90 deg."""
    for k in range(len(poly) - 1):
        v = (poly[k + 1][0] - poly[k][0], poly[k + 1][1] - poly[k][1])
        if math.hypot(*v) < min_seg or dir_index(v, 0.5) is None:
            return False
    for k in range(1, len(poly) - 1):
        d1 = dir_index((poly[k][0] - poly[k - 1][0], poly[k][1] - poly[k - 1][1]), 0.5)
        d2 = dir_index((poly[k + 1][0] - poly[k][0], poly[k + 1][1] - poly[k][1]), 0.5)
        if min((d1 - d2) % 8, (d2 - d1) % 8) > 2:
            return False
    return True
