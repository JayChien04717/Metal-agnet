"""Fluxonium base metal with recesses for the fabrication junction leads."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
from qiskit_metal import Dict, draw

from my_fluxonium import MyFluxonium


class MyFluxoniumFab(MyFluxonium):
    """Cut a 2 µm wide, 4 µm deep notch into each facing stem end."""

    default_options = Dict(jj_metal='False', contact_notch_width='2um',
                           contact_notch_depth='4um')

    def _pads(self):
        top, bottom = super()._pads()
        p = self.p
        half_width = p.contact_notch_width / 2
        left, right = p.lead_x - half_width, p.lead_x + half_width
        depth = p.contact_notch_depth
        top_edge, bottom_edge = top.bounds[1], bottom.bounds[3]
        return (top.difference(draw.box(left, top_edge, right, top_edge + depth)),
                bottom.difference(draw.box(left, bottom_edge - depth, right, bottom_edge)))
