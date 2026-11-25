# -*- coding: utf-8 -*-
r"""Grounded tunable transmon coupler with straight ㄈ bracket arms and a Manhattan SQUID.

Combines the straight ㄈ bracket arms of ``BracketTmonCouplerDolan`` with the
2-JJ Manhattan SQUID and ``manhattan_sites()`` of ``TmonCouplerManhattan``.

Sketch (orientation = 0, full coupler, stem points to -y)::

    x_tip_L                                                           x_tip_R
     ___________                                                 ___________
    |________  |  <-- upper bracket arms (gap above FX pads) --> |  ________|
    FX1 pad  | |======+-----------------------------------+======| |  FX2 pad
    _________| |      |     bar (pad_arm_length, width)   |      | |_________
    |___________|     +---------------+   +---------------+      |___________|
                                      |   |  stem (pad_head_length)
                                   ___|___|___   island crossbar
                                    \ \   \ \    2 translated Manhattan lead pairs
                                    JJ1   JJ2
                                   \ \   \ \
                                   |_________|   ground crossbar -> chip ground
                                    ______|      flux bias line (L / T / U)
"""
from qiskit_metal import Dict

from bracket_tmon_dolan import BracketTmonCouplerDolan
from tmon_coupler_manhattan import TmonCouplerManhattan


class BracketTmonCouplerManhattan(BracketTmonCouplerDolan, TmonCouplerManhattan):
    """Tunable coupler with straight ㄈ bracket arms and a 2-JJ Manhattan SQUID."""

    component_metadata = Dict(
        short_name='BracketTCManhattan',
        _qgeometry_table_path='True',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )

    default_options = Dict(
        TmonCouplerManhattan.default_options,
        bracket_gap_left='4um',
        bracket_gap_right='4um',
        bracket_length_left='6.7um',
        bracket_length_right='18.9um',
        fx_pad_edge='264um',
        fx_pad_half_height='25um',
        bracket_etch_margin='16um',
    )

    TOOLTIP = """Tunable coupler with straight ㄈ bracket arms and a Manhattan SQUID."""
