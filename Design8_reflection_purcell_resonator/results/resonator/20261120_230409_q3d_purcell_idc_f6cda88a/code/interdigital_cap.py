# -*- coding: utf-8 -*-
r"""Tapered interdigital finger capacitor (C_in / Purcell filter input coupler)
for ``design6_reflection_purcell_resonator``.

Two facing combs along x with hexagonal sloped tapers transitioning smoothly from
the CPW center trace (`trace_width`) up to the comb height `H`, surrounded by a
conformal hexagonal etched ground pocket (`ground_gap`).

At ``orientation = 0``:
    - Pin ``a`` sits at ``-x`` (left lead, connected to the Purcell filter bus).
    - Pin ``b`` sits at ``+x`` (right lead, connected to the input launch pad).

Sketch (orientation = 0, N = 6 shown)::

              ________________________________________
             /      ____                    ____      \     <- hexagonal ground pocket
            /      /|  |====================|  |\      \       (ground_gap around comb,
           /      / |  |  ==================|  | \      \       tapered to trace_gap)
    a ====+======|  |bus====================|bus  |======+==== b
           \      \ |  |  ==================|  | /      /
            \      \|__|====================|__|/      /
             \________________________________________/
           |<-->|<->|  |<---- overlap ----->|  |<->|<-->|
           lead taper bus                      bus taper lead

Local coordinates (centered at (0, 0), design units mm via ``self.p``):
    H       = N * w + (N - 1) * g                 comb height (N fingers, width w, gap g)
    F       = overlap + end_gap                   length of one finger from its bus inner face
    S       = F + end_gap = overlap + 2 * end_gap inner spacing between the two bus bars
    x_bus   = S / 2 + bus_width                   outer x edge of each bus bar
    x_lead  = x_bus + taper_length                start of straight CPW lead
    x_pin   = x_lead + lead_length                pin plane (outward normal along -x / +x)
"""

import numpy as np
from shapely.geometry import MultiPolygon
from shapely.ops import unary_union
from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent


def _clean(geom, min_area=1e-10):
    """Drop zero-area slivers left by boolean unions (Ansys cannot import zero-area polygons)."""
    parts = [g for g in getattr(geom, 'geoms', [geom]) if g.area > min_area]
    return parts[0] if len(parts) == 1 else MultiPolygon(parts)


def alley_cap_fF(finger_count, overlap_um, finger_width_um=6.0, h_um=500.0, eps_r=11.9):
    """Analytic interdigital capacitance [fF] from the Alley / Bahl closed-form estimate.

    Parameters:
        finger_count: total number of interleaved fingers N (both combs combined, N >= 3)
        overlap_um: finger overlap length l (um)
        finger_width_um: finger width W (um, assuming equal finger width and gap)
        h_um: substrate thickness (um)
        eps_r: relative permittivity of substrate

    Returns:
        Estimated series capacitance C in fF (first estimate; replace with Q3D for final tuning).
    """
    n = int(finger_count)
    a1 = 4.409e-6 * np.tanh(0.55 * (h_um / finger_width_um) ** 0.45)
    a2 = 9.92e-6 * np.tanh(0.52 * (h_um / finger_width_um) ** 0.5)
    return float((eps_r + 1.0) * overlap_um * ((n - 3) * a1 + a2) * 1e3)


def overlap_for_cap_fF(target_fF, finger_count=10, finger_width_um=6.0, h_um=500.0, eps_r=11.9):
    """Return finger overlap [um] required for ``target_fF`` under the Alley / Bahl formula."""
    c_per_um = alley_cap_fF(finger_count, 1.0, finger_width_um=finger_width_um, h_um=h_um, eps_r=eps_r)
    return float(target_fF / c_per_um)


idc_capacitance = alley_cap_fF


class InterdigitalCap(QComponent):
    """Tapered interdigital finger capacitor between two CPW pins (`a` at -x, `b` at +x for orientation 0).

    Default Options:
        * finger_count: '10' -- total number of fingers across both combs
        * finger_width: '6um' -- width of each finger
        * finger_gap: '6um' -- edge-to-edge gap between adjacent fingers
        * overlap: '200um' -- finger-to-finger overlap length along x
        * end_gap: '6um' -- gap between each finger tip and the opposite bus bar
        * bus_width: '20um' -- x-width of each vertical bus bar
        * taper_length: '100um' -- hexagonal linear taper length from comb height H to trace_width
        * lead_length: '45um' -- straight CPW lead length between taper end and pin
        * trace_width: '10um' -- CPW center-conductor width at pins a and b
        * trace_gap: '6um' -- CPW etch gap at pins a and b
        * ground_gap: '24um' -- ground pocket clearance above/below the comb and along the taper
    """

    component_metadata = Dict(
        short_name='idc',
        _qgeometry_table_poly='True',
        _qgeometry_table_path='True',
    )

    default_options = Dict(
        pos_x='0um',
        pos_y='0um',
        orientation='0',
        chip='main',
        layer='1',
        finger_count='10',
        finger_width='6um',
        finger_gap='6um',
        overlap='200um',
        end_gap='6um',
        bus_width='20um',
        taper_length='100um',
        lead_length='45um',
        trace_width='10um',
        trace_gap='6um',
        ground_gap='24um',
    )

    TOOLTIP = 'Hexagon-tapered interdigital finger capacitor with two CPW pins (a, b).'

    @staticmethod
    def comb_geometry(n_fingers, w, g, overlap, end_gap, bus_w):
        """Return (H, S, left_comb, right_comb) in local coordinates centered at (0, 0)."""
        h_comb = n_fingers * w + (n_fingers - 1) * g
        f_len = overlap + end_gap
        s_inner = f_len + end_gap
        left_parts = [draw.rectangle(bus_w, h_comb, -s_inner / 2.0 - bus_w / 2.0, 0.0)]
        right_parts = [draw.rectangle(bus_w, h_comb, +s_inner / 2.0 + bus_w / 2.0, 0.0)]
        for i in range(n_fingers):
            y_i = -h_comb / 2.0 + w / 2.0 + i * (w + g)
            if i % 2 == 0:
                left_parts.append(draw.rectangle(f_len, w, -s_inner / 2.0 + f_len / 2.0, y_i))
            else:
                right_parts.append(draw.rectangle(f_len, w, +s_inner / 2.0 - f_len / 2.0, y_i))
        return h_comb, s_inner, unary_union(left_parts), unary_union(right_parts)

    def center_to_pin_mm(self):
        """Return distance from component center (0, 0) to either pin (`a` or `b`) in design units (mm)."""
        p = self.p
        s_inner = p.overlap + 2.0 * p.end_gap
        return s_inner / 2.0 + p.bus_width + p.taper_length + p.lead_length

    def make(self):
        p = self.p
        n_fingers = int(p.finger_count)
        if n_fingers < 2:
            raise ValueError(f'{self.name}: finger_count must be >= 2, got {n_fingers}')

        h_comb, s_inner, left_comb, right_comb = self.comb_geometry(
            n_fingers, p.finger_width, p.finger_gap, p.overlap, p.end_gap, p.bus_width
        )

        x_bus = s_inner / 2.0 + p.bus_width
        t_len = p.taper_length
        x_lead = x_bus + t_len
        x_pin = x_lead + p.lead_length
        tw_half = p.trace_width / 2.0
        tg = p.trace_gap
        gg = p.ground_gap

        if t_len > 0:
            def _trapezoid(sgn, hw_inner, hw_outer):
                return draw.Polygon([
                    (sgn * x_bus, -hw_inner),
                    (sgn * x_lead, -hw_outer),
                    (sgn * x_lead, +hw_outer),
                    (sgn * x_bus, +hw_inner),
                ])

            taper_l = _trapezoid(-1.0, h_comb / 2.0, tw_half)
            taper_r = _trapezoid(+1.0, h_comb / 2.0, tw_half)
            left_comb = _clean(unary_union([left_comb, taper_l]))
            right_comb = _clean(unary_union([right_comb, taper_r]))

            center_pocket = draw.rectangle(2.0 * x_bus, h_comb + 2.0 * gg, 0.0, 0.0)
            etch_l = _trapezoid(-1.0, h_comb / 2.0 + gg, tw_half + tg)
            etch_r = _trapezoid(+1.0, h_comb / 2.0 + gg, tw_half + tg)
            pocket = _clean(unary_union([center_pocket, etch_l, etch_r]))
        else:
            pocket = draw.rectangle(2.0 * (x_bus + gg), h_comb + 2.0 * gg, 0.0, 0.0)

        lead_a = draw.LineString([(-x_lead + 1e-6, 0.0), (-x_pin, 0.0)])
        lead_b = draw.LineString([(+x_lead - 1e-6, 0.0), (+x_pin, 0.0)])

        items = [left_comb, right_comb, pocket, lead_a, lead_b]
        items = draw.rotate(items, p.orientation, origin=(0.0, 0.0))
        items = draw.translate(items, p.pos_x, p.pos_y)
        left_comb, right_comb, pocket, lead_a, lead_b = items

        self.add_qgeometry('poly', {'comb_a': left_comb, 'comb_b': right_comb},
                           chip=p.chip, layer=p.layer)
        self.add_qgeometry('poly', {'pocket': pocket},
                           chip=p.chip, layer=p.layer, subtract=True)

        for nm, ln in (('lead_a', lead_a), ('lead_b', lead_b)):
            self.add_qgeometry('path', {nm: ln}, width=p.trace_width,
                               chip=p.chip, layer=p.layer)
            self.add_qgeometry('path', {f'{nm}_sub': ln},
                               width=p.trace_width + 2.0 * p.trace_gap,
                               chip=p.chip, layer=p.layer, subtract=True)

        # Pins: points = [inner, outer] so outward normal points from inner -> outer
        for pin_name, sgn in (('a', -1.0), ('b', +1.0)):
            pts = draw.LineString([
                (sgn * (x_pin - p.trace_width), 0.0),
                (sgn * x_pin, 0.0),
            ])
            pts = draw.rotate(pts, p.orientation, origin=(0.0, 0.0))
            pts = draw.translate(pts, p.pos_x, p.pos_y)
            self.add_pin(pin_name, points=np.asarray(pts.coords),
                         width=p.trace_width, gap=p.trace_gap, input_as_norm=True)
