"""Snapshot of the user's eight port coordinates used for this simulation.

Sketch (10 mm chip; each port's orientation points into the chip)::

                 P8      P7
          _______v_______v_______
         |                       |
     P1 >|                       |< P6
         |                       |
     P2 >|                       |< P5
         |_______^_______^_______|
                 P3      P4


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
