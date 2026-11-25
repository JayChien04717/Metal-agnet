# -*- coding: utf-8 -*-
r"""Design-specific tunable coupler with straight ㄈ bracket arms around the FX bottom pads.

Subclasses the shared ``MyTmonCouplerDolan`` / ``MyTmonCouplerDolanFab`` directly (no runtime
class monkey-patching).

Sketch (local coordinates, orientation 0, right half; left half is mirrored):

         ___________ x_tip
        |  ________   <-- upper bracket arm (gap above FX bottom pad)
  ======| |  FX pad   <-- neck joining coupler bar (pad_arm_length/2) to vertical
        | |________
        |___________  <-- lower bracket arm (gap below FX bottom pad)
"""
from shapely.affinity import scale
from shapely.geometry import LineString, box
from shapely.ops import unary_union

from qiskit_metal import Dict
from my_tmon_dolan import MyTmonCouplerDolan
from my_tmon_dolan_fab import MyTmonCouplerDolanFab


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


class _BracketMixin:
    """Adds straight ㄈ bracket arms and their etched pocket margin to a transmon coupler."""

    default_options = Dict(
        bracket_gap_left='4um',
        bracket_gap_right='4um',
        bracket_length_left='6.7um',
        bracket_length_right='18.9um',
        fx_pad_edge='264um',          # FX_X (340 um) - pad_width/2 (76 um)
        fx_pad_half_height='25um',    # pad_height/2 (25 um)
        bracket_etch_margin='16um',
    )

    def make_pocket(self):
        super().make_pocket()
        p = self.p
        gaps = (p.bracket_gap_left, p.bracket_gap_right)
        lengths = (p.bracket_length_left, p.bracket_length_right)
        if any(g < 0.002 or l <= 0.0 for g, l in zip(gaps, lengths)):
            raise ValueError('Bracket arms require gap >= 2 um and positive length on both sides.')

        fingers = []
        for side, gap, length in zip((-1, 1), gaps, lengths):
            part = bracket_arm(
                p.pad_arm_length, p.pad_arm_width,
                p.fx_pad_edge, p.fx_pad_half_height,
                gap, length,
            )
            fingers.append(part if side == 1 else scale(part, -1, 1, origin=(0, 0)))
        placed = self._place(fingers)

        table = self.design.qgeometry.tables['poly']
        for name in ('tmon', 'pad_etch'):
            idx, = table.index[(table.component == self.id) & (table.name == name)]
            original_geometry = table.at[idx, 'geometry']
            if name == 'tmon':
                # Drop the original flat paddles, retaining the bar and SQUID stem.
                x, y = p.pos_x, p.pos_y
                half_length = p.pad_arm_length / 2.0
                original_geometry = original_geometry.intersection(
                    box(x - half_length, y - 1.0, x + half_length, y + 1.0)
                ).buffer(0)
            addition = placed if name == 'tmon' else [
                g.buffer(p.bracket_etch_margin, join_style=1) for g in placed
            ]
            table.at[idx, 'geometry'] = unary_union([original_geometry, *addition]).simplify(1e-9)


class BracketTmonCoupler(_BracketMixin, MyTmonCouplerDolan):
    """Simulation / standard mask tunable coupler with straight ㄈ bracket arms."""
    default_options = Dict(**MyTmonCouplerDolan.default_options, **_BracketMixin.default_options)


class BracketTmonCouplerFab(_BracketMixin, MyTmonCouplerDolanFab):
    """Fabrication tunable coupler with Dolan contact hooks and straight ㄈ bracket arms."""
    default_options = Dict(**MyTmonCouplerDolanFab.default_options, **_BracketMixin.default_options)
