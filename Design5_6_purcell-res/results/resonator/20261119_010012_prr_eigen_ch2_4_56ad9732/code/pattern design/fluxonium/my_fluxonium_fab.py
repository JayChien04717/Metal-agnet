r"""Fluxonium base metal with recesses for the fabrication junction leads.

Sketch (orientation = 0, centre of MyFluxonium; jj_metal='False': no leads / array)::

         ________________________
        |        pad_top         |
        |_________      _________|
                  |stem|
                  |_  _|   <- notch: contact_notch_width x contact_notch_depth at x = lead_x

                           stem_gap: the fab draws the JJ + array here
                   _  _
                  | || |   <- notch
         _________|stem|_________
        |        pad_bot         |
        |________________________|
"""
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
