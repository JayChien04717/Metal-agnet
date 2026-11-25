# -*- coding: utf-8 -*-
r"""Grounded tunable transmon coupler with U-shaped claw ends and a Manhattan SQUID.

Combines the U-shaped claw ends of ``ClawCouplerDolan`` with the 2-JJ Manhattan
SQUID and ``manhattan_sites()`` of ``TmonCouplerManhattan``.

Sketch (orientation = 0, stem points to -y)::

        ____                                    ____
       |  __|__________________________________|__  |     claws (open toward the qubits)
  Q1 ->|_|  |__________________________________|  |_|<- Q2       bar (pad_arm_width)
       |____|               |   |              |____|
                            |   |  stem (pad_head_length)
                         ___|___|___   island crossbar
                          \ \   \ \    2 translated Manhattan lead pairs (no mirror)
                          JJ1   JJ2
                         \ \   \ \
                         |_________|   ground crossbar -> chip ground
                          ______|      L-shaped flux line
"""
from qiskit_metal import Dict

from claw_coupler_dolan import ClawCouplerDolan
from tmon_coupler_manhattan import TmonCouplerManhattan


class ClawCouplerManhattan(ClawCouplerDolan, TmonCouplerManhattan):
    """Grounded tunable transmon coupler with U-shaped claw ends and a Manhattan SQUID."""

    component_metadata = Dict(
        short_name='ClawTCManhattan',
        _qgeometry_table_path='True',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )

    default_options = Dict(
        TmonCouplerManhattan.default_options,
        qubit_arm_width='30um',
        claw_width='10um',
        claw_length='60um',
        claw_sep='8um',
        claw_gap='12um',
        flux_bias_line_options=Dict(
            TmonCouplerManhattan.default_options.flux_bias_line_options,
            fbl_type='L',
            mirror='False',
        ),
    )

    TOOLTIP = """Grounded tunable transmon coupler with U-shaped claw ends and a Manhattan SQUID."""
