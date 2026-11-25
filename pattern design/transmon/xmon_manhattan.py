# -*- coding: utf-8 -*-
r"""Flux-tunable (or fixed) Xmon qubit with Manhattan-style junction leads.

Subclasses ``XmonDolan`` and replaces the Dolan SQUID (or single junction) at the
south arm tip with tapered + 45-degree-bent Manhattan base-metal leads compatible
with ``TransmonJunctionManhattan``.

Sketch (orientation = 0, south arm ends in Manhattan SQUID at -y)::

                      +---+
                      |   |   north arm (cross_length, cross_width)
         +------------+   +------------+
         |   west arm       east arm   |   <- readout claw wraps chosen arm tip
         +------------+   +------------+
                      |   |   south arm
                   ___|___|___   island crossbar
                    \ \   \ \    2 translated Manhattan lead pairs (or 1 if squid.make='False')
                    JJ1   JJ2
                   \ \   \ \
                   |_________|   ground crossbar -> chip ground
                    ______|      L-shaped flux bias line

Pins: readout claw pins from ``connection_pads`` + ``flux_bias_line`` (when ``squid.make='True'``).
Call ``manhattan_sites()`` after construction to get placed site centres for
``TransmonJunctionManhattan`` / ``add_ebeam_junctions()``.
"""
from qiskit_metal import Dict, draw
from shapely.geometry import Point

from manhattan_leads import DEFAULT_MANHATTAN_LEAD_OPTIONS, _clean, manhattan_lead_pair
from metal_width import min_metal_width, require_strip
from xmon_dolan import XmonDolan, _true


class XmonManhattan(XmonDolan):
    """Flux-tunable or fixed-frequency Xmon qubit with Manhattan junction leads."""

    default_options = Dict(
        XmonDolan.default_options,
        **DEFAULT_MANHATTAN_LEAD_OPTIONS,
        jj_width='3um',
        squid=Dict(
            XmonDolan.default_options.squid,
            jj_gap='48um',
            squid_spacing='30um',
            squid_bar_height='6um',
            contact_overlap='3um',
            shift_toward_arm='0um',
            ground_tail='',
        ),
    )
    component_metadata = Dict(XmonDolan.component_metadata, short_name='XmonManhattan')

    def make_pocket(self):
        p = self.p
        unit = self.design.parse_value('1um')
        min_w = min_metal_width(unit)
        cl, cw, cg = p.cross_length, p.cross_width, p.cross_gap
        cross_line = draw.shapely.ops.unary_union([
            draw.LineString([(0, cl), (0, -cl)]),
            draw.LineString([(cl, 0), (-cl, 0)]),
        ])
        cross = cross_line.buffer(cw / 2.0, cap_style=2)
        cross_etch = cross.buffer(cg, cap_style=3, join_style=2)

        sq = p.squid
        is_squid = _true(sq.make)
        jg, ov = sq.jj_gap, sq.contact_overlap
        bar_h = sq.squid_bar_height
        sp = sq.squid_spacing
        yc = cl + jg / 2.0 - sq.shift_toward_arm

        if is_squid:
            require_strip(self.name, 'Manhattan SQUID crossbar height', bar_h, min_w)
            if sp < p.lead_base_width:
                raise ValueError(f'{self.name}: squid_spacing ({sp}) must be >= lead_base_width ({p.lead_base_width})')
            y_root = jg / 2.0 - bar_h / 2.0
        else:
            y_root = jg / 2.0

        if y_root <= p.lead_elbow_y:
            raise ValueError(
                f'{self.name}: squid.jj_gap ({jg}) too small for Manhattan leads '
                f'(y_root={y_root} <= lead_elbow_y={p.lead_elbow_y})'
            )

        if is_squid:
            x_jj = sp / 2.0
            bar_w = max(sp + p.lead_base_width, cw)
            u1, l1, jj1 = manhattan_lead_pair(y_root, p, x_offset=-x_jj)
            u2, l2, jj2 = manhattan_lead_pair(y_root, p, x_offset=+x_jj)
            u1, l1, jj1, u2, l2, jj2 = draw.scale([u1, l1, jj1, u2, l2, jj2], 1, -1, origin=(0, 0))
            u1, l1, jj1, u2, l2, jj2 = draw.translate([u1, l1, jj1, u2, l2, jj2], 0.0, yc)

            isl_y0 = min(cl - ov, yc - y_root - bar_h)
            isl_y1 = yc - y_root + 0.5 * unit
            island = _clean(draw.union(draw.box(-bar_w / 2.0, isl_y0, bar_w / 2.0, isl_y1), u1, u2))
            gnd_y0 = yc + y_root - 0.5 * unit
            gnd_y1 = max(cl + jg + ov, yc + y_root + bar_h)
            ground = _clean(draw.union(draw.box(-bar_w / 2.0, gnd_y0, bar_w / 2.0, gnd_y1), l1, l2))

            ground_edge = cl + jg
            self._squid_etch_edge = ground_edge
            win_half = max(cw / 2.0 + cg, bar_w / 2.0 + cg)
            window = draw.box(-win_half, 0.0, win_half, ground_edge)
            cross_etch = _clean(draw.union(cross_etch, window))

            jj_metal = _true(p.jj_metal)
            if jj_metal:
                cross = _clean(draw.union(cross, island))
            cross, cross_etch, ground, jj1, jj2 = self._place([cross, cross_etch, ground, jj1, jj2])
            self.add_qgeometry('poly', dict(cross=cross), chip=p.chip)
            if jj_metal:
                self.add_qgeometry('poly', dict(squid_ground=ground), chip=p.chip)
            self.add_qgeometry('poly', dict(cross_etch=cross_etch), subtract=True, chip=p.chip)
            self.add_qgeometry('junction', dict(rect_jj1=jj1), width=p.jj_width, chip=p.chip, hfss_inductance=p.Lj_a)
            self.add_qgeometry('junction', dict(rect_jj2=jj2), width=p.jj_width, chip=p.chip, hfss_inductance=p.Lj_b)
        else:
            u1, l1, jj1 = manhattan_lead_pair(y_root + ov, p, x_offset=0.0)
            u1, l1, jj1 = draw.scale([u1, l1, jj1], 1, -1, origin=(0, 0))
            u1, l1, jj1 = draw.translate([u1, l1, jj1], 0.0, yc)
            ground_edge = cl + jg
            self._squid_etch_edge = ground_edge
            window = draw.box(-(cw / 2.0 + cg), 0.0, cw / 2.0 + cg, ground_edge)
            cross_etch = _clean(draw.union(cross_etch, window))

            jj_metal = _true(p.jj_metal)
            if jj_metal:
                cross = _clean(draw.union(cross, u1))
            cross, cross_etch, l1, jj1 = self._place([cross, cross_etch, l1, jj1])
            self.add_qgeometry('poly', dict(cross=cross), chip=p.chip)
            if jj_metal:
                self.add_qgeometry('poly', dict(jj_ground=l1), chip=p.chip)
            self.add_qgeometry('poly', dict(cross_etch=cross_etch), subtract=True, chip=p.chip)
            self.add_qgeometry('junction', dict(rect_jj=jj1), width=p.jj_width, chip=p.chip, hfss_inductance=p.Lj_a)

    def manhattan_sites(self):
        """Return placed Manhattan JJ site metadata (in design units, mm) for e-beam placement."""
        p = self.p
        sq = p.squid
        yc = p.cross_length + sq.jj_gap / 2.0 - sq.shift_toward_arm
        ori = float(p.orientation)
        if _true(sq.make):
            x_jj = sq.squid_spacing / 2.0
            c, p1, p2 = self._place([Point(0.0, yc), Point(-x_jj, yc), Point(+x_jj, yc)])
            return dict(
                squid=True,
                squid_spacing=sq.squid_spacing,
                center=(c.x, c.y),
                orientation=ori,
                sites=[
                    dict(name='jj1', center=(p1.x, p1.y), orientation=ori),
                    dict(name='jj2', center=(p2.x, p2.y), orientation=ori),
                ],
            )
        (c,) = self._place([Point(0.0, yc)])
        return dict(
            squid=False,
            center=(c.x, c.y),
            orientation=ori,
            sites=[dict(name='jj1', center=(c.x, c.y), orientation=ori)],
        )

    def contact_hooks(self):
        raise NotImplementedError(
            f'{self.name}: Manhattan Xmon uses manhattan_sites() instead of Dolan contact_hooks()'
        )
