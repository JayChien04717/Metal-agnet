"""PCB launch anchors for the eight-port wafer layouts with a 4 mm side pitch.

Same arrangement as the pcbport.py sketch (P1/P2 left, P3/P4 bottom, P5/P6 right, P7/P8 top),
but the four side ports sit at y = +/- 2 mm instead of +/- 1.6 mm.
"""

from pcbport import pcb_ports


wafer_pcb_ports = {name: options.copy() for name, options in pcb_ports.items()}
for name, y in (('P1', 2), ('P2', -2), ('P5', -2), ('P6', 2)):
    wafer_pcb_ports[name]['pos_y'] = f'{y}mm'
