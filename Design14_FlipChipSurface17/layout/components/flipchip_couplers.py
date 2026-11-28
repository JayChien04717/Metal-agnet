# -*- coding: utf-8 -*-
"""Flip-chip C-tier (interposer) vertical coupling QComponents for Design14.

- VerticalReadoutCoupler: circular capacitive pad directly beneath a Q-tier qubit island,
  coupling across the ~7 um flip-chip vacuum gap and exiting radially at 45 deg.
- VerticalFluxLoop: shorted T-bar / micro-loop directly beneath a Q-tier FloatingTunableCoupler
  SQUID hub for vertical inductive flux bias across the ~7 um gap.
- IndiumBumpPad: 20 x 20 um UBM / Indium micro-bump pad (Layer 30) stitching Q-tier and C-tier
  ground planes together to suppress parallel-plate chip modes.
"""
from __future__ import annotations

import math
import numpy as np
from shapely.geometry import Point, Polygon, LineString

from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent


class VerticalReadoutCoupler(QComponent):
    """Circular vertical-coupling readout pad on the C-tier directly beneath a Q-tier qubit."""

    default_options = Dict(
        pos_x='0mm',
        pos_y='0mm',
        pad_diameter='100um',
        pocket_gap='30um',
        ro_corner='NE',
        lead_length='45um',
        trace_width='10um',
        trace_gap='6um',
        layer='1',
    )

    def make(self):
        p = self.p
        r_pad = 0.5 * p.pad_diameter
        r_pock = r_pad + p.pocket_gap
        tw = p.trace_width
        tg = p.trace_gap
        lead = p.lead_length

        corner_ang = {'NE': 45.0, 'NW': 135.0, 'SW': -135.0, 'SE': -45.0}[str(self.options.ro_corner)]
        ang = math.radians(corner_ang)
        ca, sa = math.cos(ang), math.sin(ang)

        pad = Point(0.0, 0.0).buffer(r_pad, resolution=32)
        pocket = Point(0.0, 0.0).buffer(r_pock, resolution=32)

        # Short radial CPW launch from pad edge to r_pock + lead
        p_in = (r_pad * 0.85 * ca, r_pad * 0.85 * sa)
        p_out = ((r_pock + lead) * ca, (r_pock + lead) * sa)
        stem_metal = LineString([p_in, p_out]).buffer(0.5 * tw, cap_style=2)
        stem_cut = LineString([p_in, p_out]).buffer(0.5 * tw + tg, cap_style=2)

        pad_full = pad.union(stem_metal)
        pocket_full = pocket.union(stem_cut)

        geom = draw.translate([pad_full, pocket_full], p.pos_x, p.pos_y)
        self.add_qgeometry('poly', {'ro_vpad': geom[0]}, layer=p.layer, subtract=False)
        self.add_qgeometry('poly', {'ro_vpocket': geom[1]}, layer=p.layer, subtract=True)

        pin_end = np.array([p.pos_x + p_out[0], p.pos_y + p_out[1]], float)
        pin_in = np.array([p.pos_x + (r_pock + lead - 0.01) * ca, p.pos_y + (r_pock + lead - 0.01) * sa], float)
        self.add_pin('readout', points=np.array([pin_in, pin_end]), width=tw, input_as_norm=True)


class VerticalFluxLoop(QComponent):
    """Shorted T-shaped flux bias termination on C-tier directly beneath a Q-tier coupler SQUID."""

    default_options = Dict(
        pos_x='0mm',
        pos_y='0mm',
        orientation='0',
        bar_length='80um',
        bar_width='6um',
        bar_gap='6um',
        lead_length='55um',
        trace_width='10um',
        trace_gap='6um',
        layer='1',
    )

    def make(self):
        p = self.p
        bl = p.bar_length
        bw = p.bar_width
        bg = p.bar_gap
        ll = p.lead_length
        tw = p.trace_width
        tg = p.trace_gap

        # Local frame: T-bar runs along x at y=0 (underneath SQUID), CPW feed exits along -y!
        t_cut = draw.rectangle(bl + 2.0 * bg, bw + 2.0 * bg, 0.0, 0.0)
        stem_cut = draw.rectangle(tw + 2.0 * tg, ll + bg, 0.0, -0.5 * (ll + bg))
        cut = t_cut.union(stem_cut)

        # Center conductor shorts into ground at the two ends of the T-bar (inductive loop!)
        t_metal = draw.rectangle(bl + 2.0 * bg + 0.004, bw, 0.0, 0.0)
        stem_metal = draw.rectangle(tw, ll, 0.0, -0.5 * ll)
        metal = t_metal.union(stem_metal)

        geoms = [metal, cut]
        geoms = draw.rotate(geoms, p.orientation, origin=(0, 0))
        geoms = draw.translate(geoms, p.pos_x, p.pos_y)
        self.add_qgeometry('poly', {'flux_vmetal': geoms[0]}, layer=p.layer, subtract=False)
        self.add_qgeometry('poly', {'flux_vcut': geoms[1]}, layer=p.layer, subtract=True)

        ori_rad = math.radians(float(p.orientation))
        nx, ny = math.sin(ori_rad), -math.cos(ori_rad)
        pin_end = np.array([p.pos_x + ll * nx, p.pos_y + ll * ny], float)
        pin_in = np.array([p.pos_x + (ll - 0.01) * nx, p.pos_y + (ll - 0.01) * ny], float)
        self.add_pin('flux', points=np.array([pin_in, pin_end]), width=tw, input_as_norm=True)
