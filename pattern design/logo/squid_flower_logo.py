# -*- coding: utf-8 -*-
"""Parametric 6-petal SQUID-flower logo QComponent.

Geometry (all dimensions normalized to overall tip-to-tip height ``size``):
- 6 identical hollow diamond petals at 60 deg intervals (90, 30, -30, -90,
  -150, 150 deg), radial span [0.211, 0.500], tangential half-width 0.135.
- Central rounded-square SQUID loop (outer half-size 0.188, corner radius
  0.055), clipping the inner tips of the 4 side petals cleanly.
- Left and right square Josephson-junction [X] symbols straddling the SQUID
  sidewalls at (+-0.194, 0), outer size 0.122 x 0.122, masking the underlying
  SQUID wall inside each box.
- Two touching 4-pointed concave stars (sparkles) inside the SQUID loop at
  (+-0.052, 0).

By default (``tone='etched'``) the strokes are subtracted from layer 1 ground
metal (matching dark-background / light-line logos). Set ``tone='metal'`` to
place positive metal strokes instead.
"""
from __future__ import annotations

import math
import numpy as np
from shapely.geometry import LineString, Polygon, box
from shapely.ops import unary_union
from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent


def _rounded_square_outline(half_size: float, radius: float, resolution: int = 16) -> Polygon:
    inner = box(-(half_size - radius), -(half_size - radius),
                half_size - radius, half_size - radius)
    return inner.buffer(radius, resolution=resolution)


def _stroke_ring(outer: Polygon, stroke: float) -> Polygon:
    inner = outer.buffer(-stroke, join_style=2)
    return outer.difference(inner)


def _diamond_petal(r_in: float, r_out: float, half_w: float,
                   angle_deg: float, stroke: float) -> Polygon:
    r_mid = 0.5 * (r_in + r_out)
    pts_local = [
        (r_in, 0.0),
        (r_mid, half_w),
        (r_out, 0.0),
        (r_mid, -half_w),
    ]
    ang = math.radians(angle_deg)
    ca, sa = math.cos(ang), math.sin(ang)
    pts = [(x * ca - y * sa, x * sa + y * ca) for x, y in pts_local]
    outer = Polygon(pts)
    # Use bevel/mitre inset so sharp diamond tips keep clean uniform wall thickness
    inner = outer.buffer(-stroke, join_style=2, mitre_limit=4.0)
    ring = outer.difference(inner)
    # Slight corner rounding on outer/inner tips keeps fab strokes smooth
    return ring


def _jj_box(cx: float, cy: float, outer_size: float, stroke: float) -> tuple[Polygon, Polygon]:
    """Return (jj_strokes, jj_mask_box) for a square [X] junction symbol."""
    hs = 0.5 * outer_size
    r = 0.10 * outer_size
    mask_box = box(cx - (hs - r), cy - (hs - r), cx + (hs - r), cy + (hs - r)).buffer(r, resolution=8)
    frame = _stroke_ring(mask_box, stroke)
    # Diagonal X inside the box (slightly thinner than frame for crisp negative triangles)
    x_sw = 0.78 * stroke
    d_in = hs - 0.5 * stroke
    diag1 = LineString([(cx - d_in, cy - d_in), (cx + d_in, cy + d_in)]).buffer(
        0.5 * x_sw, cap_style=2)
    diag2 = LineString([(cx - d_in, cy + d_in), (cx + d_in, cy - d_in)]).buffer(
        0.5 * x_sw, cap_style=2)
    strokes = unary_union([frame, diag1, diag2]).intersection(mask_box)
    return strokes, mask_box


def _four_point_sparkle(cx: float, cy: float, rx: float, ry: float,
                        stroke: float, n_arc: int = 16) -> Polygon:
    """Concave 4-pointed star outline built from 4 quadratic Bezier curves."""
    tips = [
        (cx + rx, cy),
        (cx, cy + ry),
        (cx - rx, cy),
        (cx, cy - ry),
    ]
    # Control point close to center gives slender, crisp 4-point sparkle tips
    ctrl_frac = 0.06
    ctrls = [
        (cx + ctrl_frac * rx, cy + ctrl_frac * ry),
        (cx - ctrl_frac * rx, cy + ctrl_frac * ry),
        (cx - ctrl_frac * rx, cy - ctrl_frac * ry),
        (cx + ctrl_frac * rx, cy - ctrl_frac * ry),
    ]
    ts = np.linspace(0.0, 1.0, n_arc, endpoint=False)
    pts = []
    for i in range(4):
        p0 = np.array(tips[i])
        pc = np.array(ctrls[i])
        p1 = np.array(tips[(i + 1) % 4])
        for t in ts:
            pt = (1 - t) ** 2 * p0 + 2 * (1 - t) * t * pc + t ** 2 * p1
            pts.append((float(pt[0]), float(pt[1])))
    pts.append(pts[0])
    # Buffer the closed curve symmetrically so the star is a clean stroke
    return LineString(pts).buffer(0.5 * stroke, cap_style=1, join_style=1)


class SquidFlowerLogo(QComponent):
    """6-petal SQUID-flower decorative logo for chip ground planes."""

    default_options = Dict(
        pos_x='0mm',
        pos_y='0mm',
        orientation='0',
        size='800um',          # overall tip-to-tip height
        stroke_width='20um',   # main stroke width (petals, SQUID, JJ boxes)
        sparkle_stroke='10um', # inner 4-pointed stars stroke width
        tone='etched',         # 'etched' (subtract=True) or 'metal' (subtract=False)
        layer='1',
        chip='main',
    )

    component_metadata = Dict(
        short_name='Logo',
        _qgeometry_table_poly='True',
    )

    def make(self):
        p = self.p
        H = float(p.size)
        sw = float(p.stroke_width)
        ssw = float(p.sparkle_stroke)

        # 1. Central rounded-square SQUID loop
        sq_half = 0.188 * H
        sq_rad = 0.055 * H
        sq_outer = _rounded_square_outline(sq_half, sq_rad, resolution=16)
        sq_ring = _stroke_ring(sq_outer, sw)

        # 2. Six diamond petals (top/bottom at +-90 deg; 4 side petals at +-25 / +-155 deg)
        r_in = 0.211 * H
        r_out = 0.500 * H
        petal_hw = 0.132 * H
        petals = [
            _diamond_petal(r_in, r_out, petal_hw, deg, sw)
            for deg in (90.0, 25.0, -25.0, -90.0, -155.0, 155.0)
        ]
        petals_union = unary_union(petals).difference(sq_outer.buffer(-0.5 * sw))

        # 3. Left & right square JJ [X] boxes (mask SQUID wall behind them)
        jj_cx = 0.194 * H
        jj_size = 0.122 * H
        jj_l, mask_l = _jj_box(-jj_cx, 0.0, jj_size, sw)
        jj_r, mask_r = _jj_box(+jj_cx, 0.0, jj_size, sw)
        jj_masks = unary_union([mask_l, mask_r])

        # Combine petals + SQUID (masked under JJ boxes) + JJ boxes
        body = unary_union([
            petals_union.difference(jj_masks),
            sq_ring.difference(jj_masks),
            jj_l,
            jj_r,
        ])

        # 4. Two touching 4-pointed sparkles inside the SQUID loop
        sp_cx = 0.048 * H
        sp_rx = 0.052 * H
        sp_ry = 0.084 * H
        sp_l = _four_point_sparkle(-sp_cx, 0.0, sp_rx, sp_ry, ssw)
        sp_r = _four_point_sparkle(+sp_cx, 0.0, sp_rx, sp_ry, ssw)

        logo_geom = unary_union([body, sp_l, sp_r]).buffer(0)

        # Apply orientation and position
        logo_geom = draw.rotate(logo_geom, p.orientation, origin=(0, 0))
        logo_geom = draw.translate(logo_geom, p.pos_x, p.pos_y)

        subtract = (str(p.tone).strip().lower() != 'metal')
        self.add_qgeometry(
            'poly',
            {'logo': logo_geom},
            subtract=subtract,
            layer=p.layer,
            chip=p.chip,
        )
