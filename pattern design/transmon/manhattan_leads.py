# -*- coding: utf-8 -*-
r"""Shared base-metal landing lead geometry for Manhattan-style transmon junctions.

Matches the tapered + 45-degree-bent rounded-tip base-metal leads used with
``TransmonJunctionManhattan`` (see ``pattern design/ebeam_junction/transmon_junction_manhattan.py``):

Sketch (local frame centred at (0, 0), orientation 0, not to scale)::

         \             /   y = +y_root (width = lead_base_width, default 13 um)
          \           /    linear taper down to y = +elbow_y (10 um)
           \         /
            |       |      width = lead_tip_width (default 3 um)
             \       \
              \       \    45-deg diagonal arm along u = +arm_u (6.5 um)
               \       )   rounded semicircular tip (radius = lead_tip_width / 2)
                `-----'
                  (0, 0)   <-- Manhattan JJ centre (pads land across the two 45-deg arms)
        ,-----.
       (       \
        \       \          lower lead = exact 180-deg rotation of upper lead
         \       \
          |       |
         /         \
        /           \      y = -y_root (width = lead_base_width)

All dimensions passed to ``manhattan_lead_pair()`` are in the design's units
(normally mm) unless noted; ``unit`` is ``design.parse_value('1um')``.
"""
from math import cos, pi, sin, sqrt

from qiskit_metal import Dict, draw
from shapely.geometry import LineString, MultiPolygon, Polygon


DEFAULT_MANHATTAN_LEAD_OPTIONS = Dict(
    lead_base_width='13um',
    lead_tip_width='3um',
    lead_elbow_y='10um',
    lead_arm_u='6.5um',
    lead_tip_v='-2.755um',
    lead_elbow_radius='1.5um',
)


def _clean(geom, min_area=1e-10):
    """Drop zero-area slivers left by boolean ops (Ansys cannot draw them)."""
    parts = [g for g in getattr(geom, 'geoms', [geom]) if g.area > min_area]
    return parts[0] if len(parts) == 1 else MultiPolygon(parts)


def _upper_lead(y_root, base_w, tip_w, elbow_y, arm_u, tip_v, elbow_r, resolution=16):
    """Build the upper Manhattan base-metal lead in local coordinates centred at (0, 0).

    Starts at ``y = +y_root`` with width ``base_w``, tapers linearly to
    ``y = +elbow_y`` with width ``tip_w``, bends smoothly 45 degrees toward
    ``(+x, -y)`` along ``u = (x + y)/sqrt(2) = arm_u``, and ends in a rounded
    semicircular tip of radius ``tip_w / 2`` at ``v = (-x + y)/sqrt(2) = tip_v``.
    """
    if y_root <= elbow_y:
        raise ValueError(f'Manhattan lead y_root ({y_root}) must exceed lead_elbow_y ({elbow_y})')
    if tip_w <= 0 or base_w < tip_w:
        raise ValueError(f'Manhattan lead requires base_w ({base_w}) >= tip_w ({tip_w}) > 0')

    s2 = sqrt(2.0)
    # Corner where the vertical centreline (x = 0) meets the 45-deg arm centreline (x + y = arm_u * sqrt(2))
    y_corner = arm_u * s2
    # Tip centre in (x, y) from (u, v) = (arm_u, tip_v)
    x_tip = (arm_u - tip_v) / s2
    y_tip = (arm_u + tip_v) / s2

    # Smooth circular arc of radius elbow_r on the centreline at (0, y_corner), turning 45 deg
    # from (0, -1) to (+1/sqrt(2), -1/sqrt(2)). Tangent distance along each leg = elbow_r * tan(pi/8).
    tan_225 = sqrt(2.0) - 1.0
    t_dist = min(elbow_r * tan_225, max(0.0, y_root - y_corner) * 0.5)
    r_eff = t_dist / tan_225 if tan_225 > 0 else 0.0

    # Tapered trunk from y_root down to y_taper_end (top of the elbow arc)
    y_taper_end = max(elbow_y, y_corner + t_dist)
    taper = Polygon([
        (-base_w / 2.0, y_root),
        (+base_w / 2.0, y_root),
        (+tip_w / 2.0, y_taper_end),
        (-tip_w / 2.0, y_taper_end),
    ])

    # Centreline from y_taper_end through the filleted 45-deg elbow to (x_tip, y_tip)
    pts = [(0.0, y_taper_end)]
    if r_eff > 1e-9:
        cx_arc = r_eff
        cy_arc = y_corner + t_dist
        # Angle from pi (point (0, cy_arc)) to 5*pi/4 (tangent to 45-deg arm)
        n_arc = max(6, resolution // 2)
        for i in range(n_arc + 1):
            ang = pi + (pi / 4.0) * (i / n_arc)
            pts.append((cx_arc + r_eff * cos(ang), cy_arc + r_eff * sin(ang)))
    else:
        pts.append((0.0, y_corner))
    pts.append((x_tip, y_tip))

    arm = LineString(pts).buffer(tip_w / 2.0, cap_style=1, join_style=1, resolution=resolution)
    # Trim any tiny semicircular overhang above y_taper_end by unioning with taper
    return _clean(draw.union(taper, arm))


def manhattan_lead_pair(y_root, p, x_offset=0.0):
    """Return ``(upper_lead, lower_lead, jj_line)`` for one Manhattan junction site.

    The site is centred at ``(x_offset, 0)`` in the local frame:
    - ``upper_lead`` attaches to the upper electrode at ``y = +y_root``.
    - ``lower_lead`` is the 180-deg rotation of ``upper_lead``, attaching at ``y = -y_root``.
    - ``jj_line`` is the HFSS lumped-junction ``LineString`` bridging the inner edges of
      the two facing 45-deg arms along the +45-deg axis (passing through ``(x_offset, 0)``).
    """
    base_w = p.lead_base_width
    tip_w = p.lead_tip_width
    elbow_y = p.lead_elbow_y
    arm_u = p.lead_arm_u
    tip_v = p.lead_tip_v
    elbow_r = p.lead_elbow_radius

    upper = _upper_lead(y_root, base_w, tip_w, elbow_y, arm_u, tip_v, elbow_r)
    lower = draw.rotate(upper, 180.0, origin=(0.0, 0.0))

    # Inner edge of upper arm sits at u = arm_u - tip_w / 2; bridge along v = 0 (x = y)
    u_inner = arm_u - tip_w / 2.0
    if u_inner <= 0:
        raise ValueError(f'lead_arm_u ({arm_u}) must exceed lead_tip_width/2 ({tip_w/2})')
    c = u_inner / sqrt(2.0)
    jj_line = LineString([(c, c), (-c, -c)])

    if x_offset != 0.0:
        upper, lower, jj_line = draw.translate([upper, lower, jj_line], x_offset, 0.0)
    return upper, lower, jj_line
