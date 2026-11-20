"""Interdigital finger coupling capacitor for the transmission Purcell filter (Design4).

Two combs facing each other along x; the left comb is connected to CPW pin 'a', the right
comb to pin 'b'. Fingers alternate left / right from bottom to top.

    pin a --[bus]=|=|=|=|=[bus]-- pin b          (top view, orientation 0)
                  fingers interleave

Geometry (local coordinates, centre at 0, 0):
    H  = N w + (N - 1) g                          comb height (N fingers, width w, gap g)
    F  = overlap + end_gap                        length of one finger, measured from its bus
    S  = F + end_gap = overlap + 2 end_gap        inner distance between the two buses
    each bus is bus_width wide and H tall, the CPW lead (trace_width) leaves the bus centre
    for lead_length; the ground pocket is ground_gap around the combs, the leads carry their
    own CPW gap cut.
"""
import numpy as np
from shapely.ops import unary_union
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import QComponent


class InterdigitalCap(QComponent):
    """Interdigital capacitor between two CPW pins (a: left, b: right at orientation 0).

    Options:
        * finger_count: '10' -- total number of fingers (both combs together)
        * finger_width: '6um' -- width of each finger
        * finger_gap: '6um' -- gap between neighbouring fingers
        * overlap: '200um' -- length over which neighbouring fingers face each other
        * end_gap: '6um' -- gap between a finger tip and the opposite bus bar
        * bus_width: '20um' -- width of the bus bars the fingers hang on
        * taper_length: '0um' -- hexagon-style taper (sloped edges) from the CPW width up to the
            comb height between each bus and its lead; 0 = abrupt step. The etch follows it.
        * lead_length: '40um' -- CPW centre trace between the taper and its pin
        * trace_width: '10um' / trace_gap: '6um' -- CPW of the leads (and of the pins)
        * ground_gap: '6um' -- etch around the combs
    """
    component_metadata = Dict(short_name='idc', _qgeometry_table_poly='True',
                              _qgeometry_table_path='True')
    default_options = Dict(
        finger_count='10', finger_width='6um', finger_gap='6um',
        overlap='200um', end_gap='6um', bus_width='20um', taper_length='0um', lead_length='40um',
        trace_width='10um', trace_gap='6um', ground_gap='6um', layer='1')
    TOOLTIP = 'Interdigital finger capacitor with two CPW pins.'

    @staticmethod
    def geometry(N, w, g, overlap, eg, bw):
        """Return (H, S, left comb, right comb) in local coordinates (floats in design units)."""
        H = N * w + (N - 1) * g
        F = overlap + eg
        S = F + eg
        left = [draw.rectangle(bw, H, -S / 2 - bw / 2, 0)]
        right = [draw.rectangle(bw, H, S / 2 + bw / 2, 0)]
        for i in range(N):
            y = -H / 2 + w / 2 + i * (w + g)
            if i % 2 == 0:
                left.append(draw.rectangle(F, w, -S / 2 + F / 2, y))
            else:
                right.append(draw.rectangle(F, w, S / 2 - F / 2, y))
        return H, S, unary_union(left), unary_union(right)

    def make(self):
        p = self.p
        N = int(p.finger_count)
        H, S, left, right = self.geometry(N, p.finger_width, p.finger_gap, p.overlap,
                                          p.end_gap, p.bus_width)
        x_bus = S / 2 + p.bus_width                    # outer edge of a bus
        T = p.taper_length
        x_lead = x_bus + T                             # start of the straight CPW lead
        x_pin = x_lead + p.lead_length
        pocket = draw.rectangle(2 * (x_bus + p.ground_gap), H + 2 * p.ground_gap, 0, 0)
        tw, tg = p.trace_width / 2, p.trace_gap
        tap_l = tap_r = None
        if T > 0:
            def trap(sgn, hw0, hw1):
                return draw.Polygon([(sgn * x_bus, -hw0), (sgn * x_lead, -hw1),
                                     (sgn * x_lead, hw1), (sgn * x_bus, hw0)])
            tap_l, tap_r = trap(-1, H / 2, tw), trap(1, H / 2, tw)
            etch = [trap(s_, H / 2 + p.ground_gap, tw + tg) for s_ in (-1, 1)]
            pocket = unary_union([pocket] + etch)
            left, right = unary_union([left, tap_l]), unary_union([right, tap_r])
        lead_l = draw.LineString([(-x_lead + 1e-6, 0), (-x_pin, 0)])
        lead_r = draw.LineString([(x_lead - 1e-6, 0), (x_pin, 0)])

        items = [left, right, pocket, lead_l, lead_r]
        items = draw.rotate(items, p.orientation, origin=(0, 0))
        items = draw.translate(items, p.pos_x, p.pos_y)
        left, right, pocket, lead_l, lead_r = items

        self.add_qgeometry('poly', {'comb_a': left, 'comb_b': right}, layer=p.layer)
        self.add_qgeometry('poly', {'pocket': pocket}, layer=p.layer, subtract=True)
        for nm, ln in (('lead_a', lead_l), ('lead_b', lead_r)):
            self.add_qgeometry('path', {nm: ln}, width=p.trace_width, layer=p.layer)
            self.add_qgeometry('path', {nm + '_sub': ln}, layer=p.layer, subtract=True,
                               width=p.trace_width + 2 * p.trace_gap)

        # pins: points = [inner, outer]; the normal points outward (inner -> outer)
        for nm, sgn in (('a', -1), ('b', 1)):
            pts = np.array([(sgn * (x_pin - p.trace_width), 0), (sgn * x_pin, 0)])
            pts = np.array(draw.translate(draw.rotate(draw.LineString(pts), p.orientation,
                                                      origin=(0, 0)), p.pos_x, p.pos_y).coords)
            self.add_pin(nm, points=pts, width=p.trace_width, gap=p.trace_gap, input_as_norm=True)
