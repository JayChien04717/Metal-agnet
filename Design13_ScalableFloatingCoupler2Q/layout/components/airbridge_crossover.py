# -*- coding: utf-8 -*-
r"""Planar superconducting CPW airbridge crossover and ground-stitching bridge.

Used in 2D planar multi-qubit lattices (e.g. rotated d=3 Surface-17) to route
interior XY, coupler flux, and readout lines across orthogonal CPWs without
flip-chip bump bonding and without any base-metal (layer 1) shorts.

Layers:
- Layer 1 (base metal): continuous underpass CPW (along x in local frame) and
  two broken overpass CPW landing pads (along y in local frame), separated by
  a >= 8 um ground strip from the underpass trench.
- Layer 20 (airbridge span): suspended Al bridge metal spanning between the
  two overpass landing pads (or across a CPW for ground stitching).
- Layer 21 (airbridge piers): via / anchor feet on the landing pads.

ASCII sketch of AirbridgeCrossover (orientation 0, centred at (0, 0))::

                       over_top (+y)
                           | |
                        +--+ +--+   <- landing pad + pier (layer 21)
                        |  ###  |
            ========    |  ###  |    ========
   under_a ------------===#####===------------ under_b   (layer 1 continuous)
            ========    |  ###  |    ========
                        |  ###  |   <- airbridge span (layer 20)
                        +--+ +--+
                           | |
                       over_bot (-y)
"""
from __future__ import annotations

import math

from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent
from shapely.geometry import LineString, box
from shapely.ops import unary_union

from metal_width import min_metal_width, require_strip


class AirbridgeCrossover(QComponent):
    """Orthogonal CPW-over-CPW airbridge crossover (zero layer-1 shorts)."""

    default_options = Dict(
        pos_x='0mm',
        pos_y='0mm',
        orientation='0',         # 0: underpass along x, overpass along y
        trace_width='10um',
        trace_gap='6um',
        under_length='90um',     # total length of underpass CPW segment (layer 1)
        over_length='100um',     # total span between over_bot and over_top pins
        ground_strip='10um',     # continuous ground strip between underpass trench and overpass pad
        pad_width='18um',        # overpass landing pad width
        pad_length='18um',       # overpass landing pad length
        bridge_width='12um',     # suspended airbridge span width (layer 20)
        pier_size='12um',        # anchor pier square size (layer 21)
        layer='1',
        bridge_layer='20',
        pier_layer='21',
        chip='main',
    )

    component_metadata = Dict(
        short_name='ABCross',
        _qgeometry_table_poly='True',
    )

    def make(self):
        p = self.p
        unit = self.design.parse_value('1um')
        w = float(p.trace_width)
        g = float(p.trace_gap)
        u_len = float(p.under_length)
        o_len = float(p.over_length)
        g_strip = float(p.ground_strip)
        pw = float(p.pad_width)
        pl = float(p.pad_length)
        bw = float(p.bridge_width)
        pier = float(p.pier_size)

        require_strip(self.name, 'ground strip between underpass and overpass', g_strip, min_metal_width(unit))

        # 1. Underpass CPW (continuous along x, layer 1)
        u_half = 0.5 * u_len
        u_trench_half = 0.5 * w + g
        u_metal = box(-u_half, -0.5 * w, u_half, 0.5 * w)
        u_etch = box(-u_half, -u_trench_half, u_half, u_trench_half)

        # 2. Overpass landing pads + CPW stubs (along y, layer 1, broken across underpass)
        o_half = 0.5 * o_len
        y_etch_in = u_trench_half + g_strip
        y_pad_in = y_etch_in + g
        y_pad_out = y_pad_in + pl
        if y_pad_out >= o_half - 2.0 * unit:
            o_half = y_pad_out + 12.0 * unit

        top_pad = unary_union([
            box(-0.5 * pw, y_pad_in, 0.5 * pw, y_pad_out),
            box(-0.5 * w, y_pad_out, 0.5 * w, o_half),
        ])
        bot_pad = unary_union([
            box(-0.5 * pw, -y_pad_out, 0.5 * pw, -y_pad_in),
            box(-0.5 * w, -o_half, 0.5 * w, -y_pad_out),
        ])
        top_etch = unary_union([
            box(-0.5 * pw - g, y_etch_in, 0.5 * pw + g, y_pad_out + g),
            box(-0.5 * w - g, y_pad_out, 0.5 * w + g, o_half),
        ])
        bot_etch = unary_union([
            box(-0.5 * pw - g, -y_pad_out - g, 0.5 * pw + g, -y_etch_in),
            box(-0.5 * w - g, -o_half, 0.5 * w + g, -y_pad_out),
        ])

        # 3. Suspended airbridge metal (layer 20) & anchor piers (layer 21)
        y_pier_c = y_pad_in + 0.5 * pier
        bridge_span = box(-0.5 * bw, -(y_pier_c + 0.5 * pier), 0.5 * bw, +(y_pier_c + 0.5 * pier))
        pier_top = box(-0.5 * pier, y_pier_c - 0.5 * pier, 0.5 * pier, y_pier_c + 0.5 * pier)
        pier_bot = box(-0.5 * pier, -y_pier_c - 0.5 * pier, 0.5 * pier, -y_pier_c + 0.5 * pier)

        ang = float(p.orientation)
        px, py = float(p.pos_x), float(p.pos_y)

        def xform(geom):
            return draw.translate(draw.rotate(geom, ang, origin=(0, 0)), px, py)

        self.add_qgeometry(
            'poly', {'u_etch': xform(u_etch), 'top_etch': xform(top_etch), 'bot_etch': xform(bot_etch)},
            subtract=True, layer=p.layer, chip=p.chip,
        )
        self.add_qgeometry(
            'poly', {'under': xform(u_metal), 'over_top': xform(top_pad), 'over_bot': xform(bot_pad)},
            subtract=False, layer=p.layer, chip=p.chip,
        )
        self.add_qgeometry(
            'poly', {'bridge': xform(bridge_span)},
            subtract=False, layer=p.bridge_layer, chip=p.chip,
        )
        self.add_qgeometry(
            'poly', {'pier_top': xform(pier_top), 'pier_bot': xform(pier_bot)},
            subtract=False, layer=p.pier_layer, chip=p.chip,
        )

        # 4 outward pins: under_a (-x), under_b (+x), over_bot (-y), over_top (+y)
        step = 10.0 * unit
        pins = [
            ('under_a', [(-u_half + step, 0.0), (-u_half, 0.0)]),
            ('under_b', [(+u_half - step, 0.0), (+u_half, 0.0)]),
            ('over_bot', [(0.0, -o_half + step), (0.0, -o_half)]),
            ('over_top', [(0.0, +o_half - step), (0.0, +o_half)]),
        ]
        for pname, pts in pins:
            p0, p1 = xform(LineString(pts)).coords
            self.add_pin(pname, points=[p0, p1], width=p.trace_width, gap=p.trace_gap, input_as_norm=True)


class AirbridgeGround(QComponent):
    """Ground-stitching airbridge across a CPW line to suppress slotline modes."""

    default_options = Dict(
        pos_x='0mm',
        pos_y='0mm',
        orientation='0',         # 0: bridges across a horizontal CPW (span along y)
        span_length='48um',      # total bridge length across CPW trenches
        bridge_width='12um',
        pier_size='12um',
        bridge_layer='20',
        pier_layer='21',
        chip='main',
    )

    component_metadata = Dict(
        short_name='ABGnd',
        _qgeometry_table_poly='True',
    )

    def make(self):
        p = self.p
        s_half = 0.5 * float(p.span_length)
        bw = float(p.bridge_width)
        pier = float(p.pier_size)

        bridge = box(-0.5 * bw, -s_half, 0.5 * bw, +s_half)
        pier_top = box(-0.5 * pier, s_half - pier, 0.5 * pier, s_half)
        pier_bot = box(-0.5 * pier, -s_half, 0.5 * pier, -s_half + pier)

        ang = float(p.orientation)
        px, py = float(p.pos_x), float(p.pos_y)

        def xform(geom):
            return draw.translate(draw.rotate(geom, ang, origin=(0, 0)), px, py)

        self.add_qgeometry('poly', {'bridge': xform(bridge)}, subtract=False, layer=p.bridge_layer, chip=p.chip)
        self.add_qgeometry(
            'poly', {'pier_top': xform(pier_top), 'pier_bot': xform(pier_bot)},
            subtract=False, layer=p.pier_layer, chip=p.chip,
        )


class AirbridgeSpan(QComponent):
    """Long-span CPW airbridge (layer 20) jumping cleanly over a coupler pocket without layer-1 contact."""

    default_options = Dict(
        pos_x='0mm',
        pos_y='0mm',
        orientation='0',         # 0: spans vertically (along y) from bot (-y) to top (+y)
        trace_width='10um',
        trace_gap='6um',
        clear_span='34um',       # layer-1 open gap along overpass axis (centerline)
        total_length='74um',     # pin-to-pin span (over_bot to over_top)
        pad_width='10um',        # flush with 10 um CPW trace -> zero impedance step!
        pad_length='10um',
        bridge_width='8um',
        pier_size='8um',
        under_angle='0',         # underpass angle relative to local +x (0 = 90 deg crossing)
        under_keepout='15um',    # half-width of underpass trench (11 um) + unbroken ground strip (4 um)
        end_gap='2.5um',         # isolation slot between ground strip and overpass pad -> 35.0 um pier-to-pier arch!
        layer='1',
        bridge_layer='20',
        pier_layer='21',
        chip='main',
    )

    component_metadata = Dict(
        short_name='ABSpan',
        _qgeometry_table_poly='True',
    )

    def make(self):
        p = self.p
        unit = self.design.parse_value('1um')
        w = float(p.trace_width)
        g = float(p.trace_gap)
        c_half = 0.5 * float(p.clear_span)
        o_half = 0.5 * float(p.total_length)
        pw = float(p.pad_width)
        pl = float(p.pad_length)
        bw = float(p.bridge_width)
        pier = float(p.pier_size)
        u_ang = math.radians(float(p.under_angle))
        u_ko = float(p.under_keepout)
        eg = float(p.end_gap) if 'end_gap' in p else g

        # Ensure overpass etch NEVER intrudes into underpass trench + ground strip (even at 45 deg!)
        cos_u = max(abs(math.cos(u_ang)), 0.35)
        etch_hw = 0.5 * pw + g
        min_c_half = (u_ko + etch_hw * abs(math.sin(u_ang))) / cos_u
        y_etch_in = max(c_half, min_c_half)
        y_pad_in = y_etch_in + eg
        y_pad_out = y_pad_in + pl
        if y_pad_out >= o_half - 4.0 * unit:
            o_half = y_pad_out + 8.0 * unit

        top_pad = unary_union([
            box(-0.5 * pw, y_pad_in, 0.5 * pw, y_pad_out),
            box(-0.5 * w, y_pad_out, 0.5 * w, o_half),
        ])
        bot_pad = unary_union([
            box(-0.5 * pw, -y_pad_out, 0.5 * pw, -y_pad_in),
            box(-0.5 * w, -o_half, 0.5 * w, -y_pad_out),
        ])
        top_etch = unary_union([
            box(-0.5 * pw - g, y_etch_in, 0.5 * pw + g, y_pad_out + g),
            box(-0.5 * w - g, y_pad_out, 0.5 * w + g, o_half),
        ])
        bot_etch = unary_union([
            box(-0.5 * pw - g, -y_pad_out - g, 0.5 * pw + g, -y_etch_in),
            box(-0.5 * w - g, -o_half, 0.5 * w + g, -y_pad_out),
        ])
        # Clip overpass etch against the underpass protection corridor (parallel to underpass!)
        if u_ko > 0:
            u_corr = draw.rotate(box(-0.20, -u_ko, 0.20, u_ko), math.degrees(u_ang), origin=(0, 0))
            top_etch = top_etch.difference(u_corr)
            bot_etch = bot_etch.difference(u_corr)

        y_pier_c = y_pad_in + 0.5 * pier
        bridge_span = box(-0.5 * bw, -(y_pier_c + 0.5 * pier), 0.5 * bw, +(y_pier_c + 0.5 * pier))
        pier_top = box(-0.5 * pier, y_pier_c - 0.5 * pier, 0.5 * pier, y_pier_c + 0.5 * pier)
        pier_bot = box(-0.5 * pier, -y_pier_c - 0.5 * pier, 0.5 * pier, -y_pier_c + 0.5 * pier)

        ang = float(p.orientation)
        px, py = float(p.pos_x), float(p.pos_y)

        def xform(geom):
            return draw.translate(draw.rotate(geom, ang, origin=(0, 0)), px, py)

        self.add_qgeometry('poly', {'top_etch': xform(top_etch), 'bot_etch': xform(bot_etch)},
                           subtract=True, layer=p.layer, chip=p.chip)
        self.add_qgeometry('poly', {'over_top': xform(top_pad), 'over_bot': xform(bot_pad)},
                           subtract=False, layer=p.layer, chip=p.chip)
        self.add_qgeometry('poly', {'bridge': xform(bridge_span)},
                           subtract=False, layer=p.bridge_layer, chip=p.chip)
        self.add_qgeometry('poly', {'pier_top': xform(pier_top), 'pier_bot': xform(pier_bot)},
                           subtract=False, layer=p.pier_layer, chip=p.chip)

        step = 10.0 * unit
        for pname, pts in (('over_bot', [(0.0, -o_half + step), (0.0, -o_half)]),
                           ('over_top', [(0.0, +o_half - step), (0.0, +o_half)])):
            p0, p1 = xform(LineString(pts)).coords
            self.add_pin(pname, points=[p0, p1], width=p.trace_width, gap=p.trace_gap, input_as_norm=True)
