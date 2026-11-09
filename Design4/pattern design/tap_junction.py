"""CPW T-junction on the feedline: a short straight prime line with a stub leaving downward.

The stub ends at pin 'second_end' where a lumped coupler (InterdigitalCap, rotated -90 deg) is
attached; the resonator continues behind the cap. Point-like replacement of the long
parallel-coupled CoupledLineTee (coupling capacitance now sits in one finger capacitor)."""
import numpy as np
from shapely.ops import unary_union
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import QComponent


class TapJunction(QComponent):
    """Pins: prime_start (left), prime_end (right), second_end (bottom, normal pointing down)."""
    component_metadata = Dict(short_name='tap', _qgeometry_table_poly='True')
    default_options = Dict(prime_length='300um', stub_length='120um', stub_up='False',
                           trace_width='10um', trace_gap='6um', layer='1')

    def make(self):
        p = self.p
        L, S, tw, tg = p.prime_length, p.stub_length, p.trace_width, p.trace_gap
        prime = draw.LineString([(-L / 2, 0), (L / 2, 0)])
        sg = 1 if str(p.stub_up) == 'True' else -1       # stub direction (+y when stub_up)
        stub = draw.LineString([(0, 0), (0, sg * S)])
        prime, stub = draw.rotate([prime, stub], p.orientation, origin=(0, 0))
        prime, stub = draw.translate([prime, stub], p.pos_x, p.pos_y)
        # one merged polygon for the metal and one for the ground cut (a T), so the HFSS subtract has no
        # overlapping tool bodies
        trace = unary_union([prime.buffer(tw / 2, cap_style=2), stub.buffer(tw / 2, cap_style=2)])
        cut = unary_union([prime.buffer(tw / 2 + tg, cap_style=2), stub.buffer(tw / 2 + tg, cap_style=2)])
        self.add_qgeometry('poly', {'tap_trace': trace}, layer=p.layer)
        self.add_qgeometry('poly', {'tap_cut': cut}, layer=p.layer, subtract=True)

        def pin(pts):
            g = draw.translate(draw.rotate(draw.LineString(pts), p.orientation, origin=(0, 0)),
                               p.pos_x, p.pos_y)
            return np.array(g.coords)

        self.add_pin('prime_start', points=pin([(-L / 2 + tw, 0), (-L / 2, 0)]), width=tw, gap=tg,
                     input_as_norm=True)
        self.add_pin('prime_end', points=pin([(L / 2 - tw, 0), (L / 2, 0)]), width=tw, gap=tg,
                     input_as_norm=True)
        self.add_pin('second_end', points=pin([(0, sg * (S - tw)), (0, sg * S)]), width=tw, gap=tg,
                     input_as_norm=True)
