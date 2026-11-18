r"""Wire-bond launch pad with a separate (short) back gap.

Sketch (orientation = 0; the chip edge is to the left)::

     back_gap  pad_gap
       |<->|  _____________________________
       |   | |                             \_______________________
       |   | |                                                     |
       |   | |  pad (pad_width x                taper    lead ===  | tie (pin, +x)
       |   | |       pad_height)            (taper_height)         |
       |   | |                              _______________________|
       |   | |_____________________________/
       |___|______ etched pocket: pad_gap at the sides, back_gap behind the pad ______
           ^ 'bond' pin on the pad back edge (HFSS lumped-port plane)

Same drawing as qiskit_metal LaunchpadWirebond (pad - linear taper - lead, origin at the taper / lead junction, pin
'tie' at +lead_length), except that the etch behind the pad (toward the chip edge) is `back_gap` instead of `pad_gap`:
the side gaps set the pad impedance, the back gap only adds bond length.
"""
from math import cos, radians, sin

from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond


LP = dict(chip='main', lead_length='30um', pad_width='400um', pad_height='300um',
          pad_gap='200um', taper_height='300um', back_gap='50um')
LP_INSET = 0.25  # mm inward from the nominal PCB anchor


def add_bond_pad(design, name, placement, *, inset=LP_INSET, **options):
    """Place the standard pad on a PCB anchor, preserving its azimuth and CPW."""
    theta = radians(float(placement['orientation']))
    x = design.parse_value(placement['pos_x']) + inset * cos(theta)
    y = design.parse_value(placement['pos_y']) + inset * sin(theta)
    return BondPad(design, name, options=dict(
        placement, pos_x=f'{x:.9f}mm', pos_y=f'{y:.9f}mm', **options))


class BondPad(LaunchpadWirebond):
    default_options = Dict(LaunchpadWirebond.default_options, **LP)

    def make(self):
        p = self.p
        tw, tg, lead = p.trace_width / 2, p.trace_gap, p.lead_length
        pw, ph, pg, th, bg = p.pad_width / 2, p.pad_height, p.pad_gap, p.taper_height, p.back_gap
        pad = draw.Polygon([(0, tw), (-th, pw), (-(ph + th), pw), (-(ph + th), -pw), (-th, -pw), (0, -tw),
                            (lead, -tw), (lead, tw)])
        pocket = draw.Polygon([(0, tw + tg), (-th, pw + pg), (-(ph + th), pw + pg), (-(ph + th + bg), pw + pg),
                               (-(ph + th + bg), -(pw + pg)), (-(ph + th), -(pw + pg)), (-th, -(pw + pg)),
                               (0, -(tw + tg)), (lead, -(tw + tg)), (lead, tw + tg)])
        pin = draw.LineString([(lead, tw), (lead, -tw)])
        x_back = -(ph + th)
        bond = draw.LineString([(x_back, -pw), (x_back, pw)])     # pad back edge, normal toward the chip edge
        items = draw.translate(draw.rotate([pin, pad, pocket, bond], p.orientation, origin=(0, 0)), p.pos_x, p.pos_y)
        pin, pad, pocket, bond = items
        self.add_qgeometry('poly', dict(launch_pad=pad), layer=p.layer)
        self.add_qgeometry('poly', dict(pocket=pocket), subtract=True, layer=p.layer)
        self.add_pin('tie', pin.coords, p.trace_width, gap=p.trace_gap)
        # 'bond': lumped-port reference plane across the back gap (HFSS pad S11, simulation/driven/hfss_pad.py)
        self.add_pin('bond', bond.coords, p.pad_width, gap=bg)
