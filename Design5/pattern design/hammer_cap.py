"""Hammer capacitor C_j (Cai et al., Chin. Phys. B 33, 068501 (2024), Fig. 3(e) inset): the CPW centre trace of the
bridge ends in a rectangular hammer head on each side; the two heads are PARALLEL plates that face each other
across a narrow gap (plates perpendicular to the bridge, i.e. parallel to the two long lines it connects).

    pin a ---handle===|head_a|  gap  |head_b|===handle--- pin b        (top view, orientation 0)

head_height is across the line (about three trace widths), head_width along it. The etched pocket follows the
outline at ground_gap, the handles carry the normal CPW gap; head corners are rounded (corner_radius).
Local frame at orientation 0: pin a at -x, pin b at +x, centre at (0, 0).
"""
import numpy as np
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import QComponent


class HammerCap(QComponent):
    """Options:
        head_width: thickness of each head plate along the line     head_height: plate length across the line
        head_gap: gap between the two parallel plates               lead_length: handle (CPW trace) from head to pin
        trace_width / trace_gap: CPW of the handles and the pins    ground_gap: etch around the heads
        corner_radius: rounding of the head corners
    """
    component_metadata = Dict(short_name='ham', _qgeometry_table_poly='True', _qgeometry_table_path='True')
    default_options = Dict(head_width='25um', head_height='65um', head_gap='15um', lead_length='10um',
                           trace_width='10um', trace_gap='6um', ground_gap='16um', corner_radius='2um', layer='1')

    def make(self):
        p = self.p
        hw, hh, g = p.head_width, p.head_height, p.head_gap
        tw, tg, gg, r = p.trace_width, p.trace_gap, p.ground_gap, p.corner_radius
        left = draw.rectangle(hw, hh, -g / 2 - hw / 2, 0)
        right = draw.rectangle(hw, hh, g / 2 + hw / 2, 0)
        if r > 0:                                        # no sharp metal corners (skill/chip_drc.md C1)
            left = left.buffer(-r, join_style=1).buffer(r, join_style=1)
            right = right.buffer(-r, join_style=1).buffer(r, join_style=1)
        x_head = g / 2 + hw                              # outer edge of a head
        x_pin = x_head + p.lead_length
        pocket = draw.rectangle(2 * (x_head + gg), hh + 2 * gg, 0, 0)
        if r > 0:
            pocket = pocket.buffer(-r, join_style=1).buffer(r, join_style=1)
        lead_a = draw.LineString([(-x_head + 1e-6, 0), (-x_pin, 0)])
        lead_b = draw.LineString([(x_head - 1e-6, 0), (x_pin, 0)])

        items = [left, right, pocket, lead_a, lead_b]
        items = draw.rotate(items, p.orientation, origin=(0, 0))
        items = draw.translate(items, p.pos_x, p.pos_y)
        left, right, pocket, lead_a, lead_b = items

        self.add_qgeometry('poly', {'head_a': left, 'head_b': right}, layer=p.layer)
        self.add_qgeometry('poly', {'pocket': pocket}, layer=p.layer, subtract=True)
        for nm, ln in (('lead_a', lead_a), ('lead_b', lead_b)):
            self.add_qgeometry('path', {nm: ln}, width=tw, layer=p.layer)
            self.add_qgeometry('path', {nm + '_sub': ln}, layer=p.layer, subtract=True, width=tw + 2 * tg)
        for nm, sgn in (('a', -1), ('b', 1)):
            pts = np.array([(sgn * (x_pin - tw), 0), (sgn * x_pin, 0)])
            pts = np.array(draw.translate(draw.rotate(draw.LineString(pts), p.orientation, origin=(0, 0)),
                                          p.pos_x, p.pos_y).coords)
            self.add_pin(nm, points=pts, width=tw, gap=tg, input_as_norm=True)
