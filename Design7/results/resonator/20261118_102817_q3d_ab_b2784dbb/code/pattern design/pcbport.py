"""Snapshot of the user's eight port coordinates used for this simulation.

The original custom-fluxonium directory has priority on sys.path when present.
This copy makes the notebook reproducible if that external file is moved.
"""
pcb_ports={
    'P1':dict(pos_x='-4.5mm',pos_y='1.6mm',orientation='0'),
    'P2':dict(pos_x='-4.5mm',pos_y='-1.6mm',orientation='0'),
    'P3':dict(pos_x='-1.1mm',pos_y='-4.5mm',orientation='90'),
    'P4':dict(pos_x='1.1mm',pos_y='-4.5mm',orientation='90'),
    'P5':dict(pos_x='4.5mm',pos_y='-1.6mm',orientation='180'),
    'P6':dict(pos_x='4.5mm',pos_y='1.6mm',orientation='180'),
    'P7':dict(pos_x='1.1mm',pos_y='4.5mm',orientation='270'),
    'P8':dict(pos_x='-1.1mm',pos_y='4.5mm',orientation='270'),
}
