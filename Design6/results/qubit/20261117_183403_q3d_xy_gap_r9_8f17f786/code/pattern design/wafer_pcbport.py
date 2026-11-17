"""PCB launch anchors for the eight-port wafer layouts with a 4 mm side pitch."""

from pcbport import pcb_ports


wafer_pcb_ports = {name: options.copy() for name, options in pcb_ports.items()}
for name, y in (('P1', 2), ('P2', -2), ('P5', -2), ('P6', 2)):
    wafer_pcb_ports[name]['pos_y'] = f'{y}mm'
