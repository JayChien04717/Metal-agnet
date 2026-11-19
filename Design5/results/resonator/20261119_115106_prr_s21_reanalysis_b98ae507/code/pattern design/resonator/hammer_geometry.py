"""Hammer capacitor head / pocket geometry shared by PurcellReadout (pattern design/resonator/purcell_readout.py)
and the archived standalone HammerCap (pattern design/archive/hammer_cap.py).

Sketch (orientation 0, centre at the origin)::

     ________________________________
    |  ground_gap   ____   ____      |   <- etched pocket (rounded corners)
    |              |    | |    |     |
    |              |head| |head|     |   head_height across the line, head_width along it,
    |              | a  | | b  |     |   head_gap between the two plates
    |              |____| |____|     |
    |________________________________|
"""
from qiskit_metal import draw


def hammer_shapes(head_width, head_height, head_gap, ground_gap, corner_radius):
    """Local (orientation 0, centre at the origin) heads and etched pocket of the hammer.

    Returns (head_a, head_b, pocket, x_head) with x_head the outer edge of a head. Shared with
    purcell_readout.PurcellReadout so both draw the identical hammer."""
    hw, hh, g, gg, r = head_width, head_height, head_gap, ground_gap, corner_radius
    left = draw.rectangle(hw, hh, -g / 2 - hw / 2, 0)
    right = draw.rectangle(hw, hh, g / 2 + hw / 2, 0)
    if r > 0:                                        # no sharp metal corners (skill/chip_drc.md C1)
        left = left.buffer(-r, join_style=1).buffer(r, join_style=1)
        right = right.buffer(-r, join_style=1).buffer(r, join_style=1)
    x_head = g / 2 + hw                              # outer edge of a head
    pocket = draw.rectangle(2 * (x_head + gg), hh + 2 * gg, 0, 0)
    if r > 0:
        pocket = pocket.buffer(-r, join_style=1).buffer(r, join_style=1)
    return left, right, pocket, x_head
