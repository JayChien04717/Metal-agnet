"""T-type enclosing capacitor (Cai et al., Chin. Phys. B 33, 068501 (2024), Fig. 3(c)).

A rectangular plate (central electrode, connected to pin b through the opening) sits inside a U-shaped
frame (surrounding electrode, connected to pin a); the opening of the frame points toward pin b
(the filter side). Local coordinates at orientation 0: pin a at -x, pin b at +x.

    a --[ U frame |  [plate]  )---- b          (frame legs reach past the plate; plate lead leaves
                                                through the opening)
"""
import numpy as np
from shapely.ops import unary_union
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import QComponent


class EnclosingCap(QComponent):
    """Options (all lengths in um unless noted):
        plate_width / plate_length: size of the central plate (width across y, length along x)
        gap: spacing between plate and frame      wall: frame back-wall thickness (pin a side)
        wall_side: thickness of the frame side walls (the frame is a thick solid block with a slot, as in the paper)
        leg_ext: how far the frame legs reach beyond the plate end (the opening)
        lead_length: CPW centre trace between the structure and each pin
        trace_width / trace_gap: CPW of the leads;  ground_gap: etch around the frame
    """
    component_metadata = Dict(short_name='enc', _qgeometry_table_poly='True', _qgeometry_table_path='True')
    default_options = Dict(plate_width='193um', plate_length='41um', gap='8um', wall='52um', wall_side='57um', leg_ext='72um',
                           lead_length='40um', trace_width='10um', trace_gap='6um', ground_gap='12um', corner_radius='8um', layer='1')

    def make(self):
        p = self.p
        pw, pl_, g, t = p.plate_width, p.plate_length, p.gap, p.wall
        ts = p.wall_side
        ext, tw, tg, gg = p.leg_ext, p.trace_width, p.trace_gap, p.ground_gap
        plate = draw.rectangle(pl_, pw, 0, 0)
        x0, x1 = -pl_ / 2 - g - t, pl_ / 2 + ext            # outer frame extent along x
        yo = pw / 2 + g + ts                                    # outer half height
        yi = pw / 2 + g                                        # inner (cavity) half height
        outer = draw.rectangle(x1 - x0, 2 * yo, (x0 + x1) / 2, 0)
        # cavity = chamber around the plate + a narrow channel (just the CPW stem and its gaps) to the opening
        xc1 = pl_ / 2 + g
        chamber = draw.rectangle(xc1 - (x0 + t), 2 * yi, ((x0 + t) + xc1) / 2, 0)
        e = 0.001                                            # 1 um overlap (design units are mm)
        channel = draw.rectangle(x1 + e - (xc1 - e), 2 * (tw / 2 + tg), ((xc1 - e) + x1 + e) / 2, 0)
        cavity = unary_union([chamber, channel])
        frame = outer.difference(cavity)
        # no sharp metal corners (skill/chip_drc.md C1): round convex and concave corners of plate and frame
        r = p.corner_radius
        plate = plate.buffer(-r, join_style=1).buffer(r, join_style=1)
        frame = frame.buffer(-r, join_style=1).buffer(r, join_style=1).buffer(r, join_style=1).buffer(-r, join_style=1)
        x_pin = x1 + p.lead_length
        lead_a = draw.LineString([(x0 + 1e-6, 0), (-x_pin, 0)])
        lead_b = draw.LineString([(pl_ / 2 - 1e-6, 0), (x_pin, 0)])
        pocket = outer.buffer(gg, cap_style=2, join_style=1)      # rounded pocket corners
        items = [plate, frame, pocket, lead_a, lead_b]
        items = draw.rotate(items, p.orientation, origin=(0, 0))
        items = draw.translate(items, p.pos_x, p.pos_y)
        plate, frame, pocket, lead_a, lead_b = items
        self.add_qgeometry('poly', {'plate': plate, 'frame': frame}, layer=p.layer)
        self.add_qgeometry('poly', {'pocket': pocket}, layer=p.layer, subtract=True)
        for nm, ln in (('lead_a', lead_a), ('lead_b', lead_b)):
            self.add_qgeometry('path', {nm: ln}, width=tw, layer=p.layer)
            self.add_qgeometry('path', {nm + '_sub': ln}, layer=p.layer, subtract=True, width=tw + 2 * tg)
        for nm, sgn in (('a', -1), ('b', 1)):
            pts = np.array([(sgn * (x_pin - tw), 0), (sgn * x_pin, 0)])
            pts = np.array(draw.translate(draw.rotate(draw.LineString(pts), p.orientation, origin=(0, 0)),
                                          p.pos_x, p.pos_y).coords)
            self.add_pin(nm, points=pts, width=tw, gap=tg, input_as_norm=True)
