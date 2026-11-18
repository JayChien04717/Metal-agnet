"""Snapshot of the user's eight port coordinates used for this simulation.

The original custom-fluxonium directory has priority on sys.path when present.
This copy makes the notebook reproducible if that external file is moved.
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
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
