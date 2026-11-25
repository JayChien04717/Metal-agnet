"""PCB launch anchors for the symmetric eight-port wafer layouts (4 mm pitch on all four sides).

Same port numbering and inward orientations as pcbport.py (P1/P2 left, P3/P4 bottom,
P5/P6 right, P7/P8 top), matching a 10.4 x 10.4 mm package cutout with a centred 10 x 10 mm chip:
  - Left / right (P1, P2, P5, P6): x = +/- 4.5 mm, y = +/- 2.0 mm (4 mm vertical pitch)
  - Top / bottom (P3, P4, P7, P8): x = +/- 2.0 mm, y = +/- 4.5 mm (4 mm horizontal pitch)
"""

from pcbport import pcb_ports


wafer_pcb_ports = {name: options.copy() for name, options in pcb_ports.items()}
for name, y in (('P1', 2), ('P2', -2), ('P5', -2), ('P6', 2)):
    wafer_pcb_ports[name]['pos_y'] = f'{y}mm'
for name, x in (('P3', -2), ('P4', 2), ('P7', 2), ('P8', -2)):
    wafer_pcb_ports[name]['pos_x'] = f'{x}mm'
