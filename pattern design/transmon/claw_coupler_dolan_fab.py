# -*- coding: utf-8 -*-
r"""Fabrication mask variant of ``ClawCouplerDolan`` with Dolan contact notches.

Combines the claw ends of ``ClawCouplerDolan`` with the Dolan contact notches and
``contact_hooks()`` of ``TmonCouplerDolanNthuFab`` (2 x 4 um notches, island -1 um /
ground +1 um in x, keeping the +2 um relative Dolan offset while centring the pair
on each 10 um SQUID leg).
"""

from qiskit_metal import Dict

from claw_coupler_dolan import ClawCouplerDolan
from tmon_coupler_dolan_nthu_fab import TmonCouplerDolanNthuFab


class ClawCouplerDolanFab(ClawCouplerDolan, TmonCouplerDolanNthuFab):
    """Fab mask version of ``ClawCouplerDolan`` with Dolan contact notches."""

    default_options = Dict(
        ClawCouplerDolan.default_options,
        **{
            **TmonCouplerDolanNthuFab.default_options,
            'island_notch_offset': '-1um',
            'ground_notch_offset': '1um',
        },
    )
