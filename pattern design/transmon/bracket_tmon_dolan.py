# -*- coding: utf-8 -*-
r"""Tunable coupler with straight ㄈ bracket arms around the FX bottom pads.

Subclasses ``TmonCouplerDolanNthu`` / ``TmonCouplerDolanNthuFab`` and overrides
only ``_end_geometry()``: stem, SQUID, Dolan fab notches, flux bias line, charge
line, readout pad, and 4 um DRC checks are all inherited from the base class.

Sketch (orientation = 0, full coupler, stem points to -y)::

    x_tip_L                                                           x_tip_R
     ___________                                                 ___________
    |________  |  <-- upper bracket arms (gap above FX pads) --> |  ________|
    FX1 pad  | |======+-----------------------------------+======| |  FX2 pad
    _________| |      |     bar (pad_arm_length, width)   |      | |_________
    |___________|     +---------------+   +---------------+      |___________|
      (bracket_gap_left,              |   |  stem (pad_head_length)
       bracket_length_left)          _|___|_
                                    |x / \ x|  SQUID (2 parallel Dolan JJs)
                                    |__\__/__|  ground side
                                     ______|   flux bias line (L / T / U)

    Left bracket mirrors the right bracket; neck joins bar end (+-pad_arm_length/2)
    to each vertical spine (x_v = fx_pad_edge - gap - pad_arm_width/2).
"""

from shapely.affinity import scale
from shapely.geometry import LineString
from shapely.ops import unary_union

from qiskit_metal import Dict, draw

from tmon_coupler_dolan_nthu import TmonCouplerDolanNthu
from tmon_coupler_dolan_nthu_fab import TmonCouplerDolanNthuFab


def bracket_arm(pad_arm_length, pad_arm_width, pad_edge, pad_half_height, gap, length):
    """Straight ㄈ around the inner end of one FX bottom pad (local mm, pad at +x)."""
    half = pad_arm_width / 2.0
    x_v = pad_edge - gap - half
    y_a = pad_half_height + gap + half
    x_tip = pad_edge + length - half
    if x_tip <= x_v or x_v <= pad_arm_length / 2.0:
        raise ValueError('Bracket arm needs room between the bar end and the pad.')
    bracket = LineString([(x_tip, y_a), (x_v, y_a), (x_v, -y_a), (x_tip, -y_a)]).buffer(
        half, cap_style=1, join_style=2, mitre_limit=2.0
    )
    neck = LineString([(pad_arm_length / 2.0 - 0.018, 0.0), (x_v, 0.0)]).buffer(half, cap_style=2)
    return unary_union([bracket, neck])


class BracketTmonCouplerDolan(TmonCouplerDolanNthu):
    """TmonCouplerDolanNthu with straight ㄈ bracket arms (overrides ``_end_geometry()`` only)."""

    default_options = Dict(
        TmonCouplerDolanNthu.default_options,
        bracket_gap_left='4um',
        bracket_gap_right='4um',
        bracket_length_left='6.7um',
        bracket_length_right='18.9um',
        fx_pad_edge='264um',          # FX_X (340 um) - pad_width/2 (76 um)
        fx_pad_half_height='25um',    # pad_height/2 (25 um)
        bracket_etch_margin='16um',
    )

    def _end_geometry(self):
        p = self.p
        gaps = (p.bracket_gap_left, p.bracket_gap_right)
        lengths = (p.bracket_length_left, p.bracket_length_right)
        if any(g < 0.002 or l <= 0.0 for g, l in zip(gaps, lengths)):
            raise ValueError('Bracket arms require gap >= 2 um and positive length on both sides.')
        bar = draw.rectangle(p.pad_arm_length, p.pad_arm_width, 0, 0)
        arms = []
        for side, gap, length in zip((-1, 1), gaps, lengths):
            part = bracket_arm(
                p.pad_arm_length, p.pad_arm_width,
                p.fx_pad_edge, p.fx_pad_half_height,
                gap, length,
            )
            arms.append(part if side == 1 else scale(part, -1, 1, origin=(0, 0)))
        # Keep the base paddle-etch boxes as well so pad_etch matches the pre-refactor mask identically.
        _, _, base_etches = super()._end_geometry()
        etches = [*base_etches, *(g.buffer(p.bracket_etch_margin, join_style=1) for g in arms)]
        return bar, arms, etches


class BracketTmonCouplerDolanFab(BracketTmonCouplerDolan, TmonCouplerDolanNthuFab):
    """Fabrication tunable coupler with Dolan contact hooks and straight ㄈ bracket arms."""

    default_options = Dict(
        BracketTmonCouplerDolan.default_options,
        **TmonCouplerDolanNthuFab.default_options,
    )
