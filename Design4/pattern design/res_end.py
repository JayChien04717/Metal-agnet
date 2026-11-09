"""Open end of a readout resonator where the qubit would be (no qubit in this device).

A component that only carries a pin: the resonator route ends open at this pin. The pin
normal points along `orientation` (deg), i.e. toward the resonator."""
import numpy as np
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import QComponent


class ResEnd(QComponent):
    """Pin-only termination ('end')."""
    component_metadata = Dict(short_name='resend')
    default_options = Dict(trace_width='10um', trace_gap='6um', layer='1')

    def make(self):
        p = self.p
        pts = draw.LineString([(-p.trace_width, 0), (0, 0)])
        pts = draw.translate(draw.rotate(pts, p.orientation, origin=(0, 0)), p.pos_x, p.pos_y)
        self.add_pin('end', points=np.array(pts.coords), width=p.trace_width, gap=p.trace_gap,
                     input_as_norm=True)
