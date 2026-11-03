# -*- coding: utf-8 -*-
r"""Test junction: an etched pocket with two pads; the junction between the pads is left
for the fab to draw (same as the lab's layout3 ftf.gds corner structures).

    ground ____________________________
          |   ________      ________   |
          |  |  pad_l |    |  pad_r |  |   <- pocket (etched)
          |  |________|    |________|  |
          |__________________________  |
                      pad_gap

Origin = pocket centre; the pad pair is centred on it (pad_l left, pad_r right,
both vertically centred). Use it for resistance / junction-quality checks next to the chip.
"""
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import QComponent


class TestJunction(QComponent):
    """Pocket + two pads for a test junction.

    Default options:
        * pocket_width, pocket_height: etched pocket
        * pad_l_width, pad_l_height / pad_r_width, pad_r_height: the two pads
        * pad_gap: distance between the pads (junction site)
    """
    default_options = Dict(
        pocket_width='800um', pocket_height='300um',
        pad_l_width='225um', pad_l_height='129um',
        pad_r_width='295um', pad_r_height='129um',
        pad_gap='211um',
        layer='1',
    )
    component_metadata = Dict(short_name='testjj')

    def make(self):
        p = self.p
        total = p.pad_l_width + p.pad_gap + p.pad_r_width
        x_l = -total / 2 + p.pad_l_width / 2
        x_r = total / 2 - p.pad_r_width / 2
        pocket = draw.rectangle(p.pocket_width, p.pocket_height, 0, 0)
        pad_l = draw.rectangle(p.pad_l_width, p.pad_l_height, x_l, 0)
        pad_r = draw.rectangle(p.pad_r_width, p.pad_r_height, x_r, 0)
        pocket, pad_l, pad_r = draw.rotate([pocket, pad_l, pad_r], p.orientation, origin=(0, 0))
        pocket, pad_l, pad_r = draw.translate([pocket, pad_l, pad_r], p.pos_x, p.pos_y)
        self.add_qgeometry('poly', dict(pocket=pocket), subtract=True, layer=p.layer)
        self.add_qgeometry('poly', dict(pad_l=pad_l, pad_r=pad_r), layer=p.layer)
