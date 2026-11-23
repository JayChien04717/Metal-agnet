# -*- coding: utf-8 -*-
r"""Grounded tunable transmon coupler whose two ends are U-shaped claws (Design6_TXT).

Same stem / SQUID / flux line / charge line as ``MyTmonCoupler`` (shared
``pattern design/my_tmon.py``), but the flat paddles are replaced by two claws
that wrap the arm tips of the neighbouring Xmon qubits (Q1 on the left, Q2 on the
right). The claw sets the qubit-coupler capacitance C_qc; the bar length is fixed by
the qubit spacing.

Sketch (orientation = 0, stem points to -y)::

        ____                                    ____
       |  __|__________________________________|__  |     claws (open toward the qubits)
  Q1 ->|_|  |__________________________________|  |_|<- Q2       bar (pad_arm_width)
       |____|               |   |              |____|
                            |   |  stem (pad_head_length)
                           _|___|_
                          |x /  \ x|  SQUID (2 parallel JJs, lumped Lj each in HFSS)
                          |__\__/__|  ground side
                           ______|   L-shaped flux line

Claw geometry (left claw, local frame before the y-flip of ``_place``):
    back plate: x in [-L/2 - claw_width, -L/2]  (faces the qubit arm tip, claw_sep away)
    fingers:    along -x, length claw_length (measured from the back plate's inner edge),
                at |y| in [arm/2 + claw_sep, arm/2 + claw_sep + claw_width]
    etch:       claw buffered by claw_gap (merges with the qubit cross etch).
``qubit_arm_width`` must equal the Xmon ``cross_width``.
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
from project_paths import pattern_design_paths as _pattern_paths
import os
import sys

from qiskit_metal import draw, Dict

HERE = _S('Design6_TXT/layout')
sys.path[:0] = [p for p in _pattern_paths() if p not in sys.path]
from my_tmon import MyTmonCoupler, _clean  # noqa: E402
from my_tmon_fab import MyTmonCouplerFab  # noqa: E402


class ClawCoupler(MyTmonCoupler):
    """MyTmonCoupler with claw ends. Extra options (see module doc):
    qubit_arm_width, claw_width, claw_length, claw_sep, claw_gap."""

    default_options = Dict(
        MyTmonCoupler.default_options,
        qubit_arm_width='30um',
        claw_width='10um',
        claw_length='60um',
        claw_sep='8um',
        claw_gap='12um',
        flux_bias_line_options=Dict(MyTmonCoupler.default_options.flux_bias_line_options, fbl_type='L', mirror='False'),
    )
    component_metadata = Dict(short_name='ClawTC', _qgeometry_table_path='True',
                              _qgeometry_table_poly='True', _qgeometry_table_junction='True')
    TOOLTIP = """Grounded tunable transmon coupler with claw ends (Design6_TXT)."""

    def _claw(self, sign):
        """Claw at x = sign * L/2 (sign -1: left, opening toward -x)."""
        p = self.p
        L, cw, cl, cs = p.pad_arm_length, p.claw_width, p.claw_length, p.claw_sep
        half_in = p.qubit_arm_width / 2 + cs           # inner edge of a finger
        half_out = half_in + cw
        x_back_in = sign * L / 2                       # back plate inner edge (faces the qubit)
        x_back_out = sign * (L / 2 - cw)
        x_tip = sign * (L / 2 + cl)                    # finger tips
        back = draw.box(min(x_back_in, x_back_out), -half_out, max(x_back_in, x_back_out), half_out)
        f_top = draw.box(min(x_back_in, x_tip), half_in, max(x_back_in, x_tip), half_out)
        f_bot = draw.box(min(x_back_in, x_tip), -half_out, max(x_back_in, x_tip), -half_in)
        return draw.union(back, f_top, f_bot)

    def make_pocket(self):
        p = self.p
        hw, hl = p.pad_head_width, p.pad_head_length
        jg, jw = p.jj_gap, p.jj_width
        L, aw, pg = p.pad_arm_length, p.pad_arm_width, p.paddle_gap
        squid_island, jj_pad_ground, _, rect_jj1, rect_jj2 = self._squid_geometry()

        stem = draw.rectangle(hw, hl, 0, hl / 2)
        # the bar runs between the two back plates (their outer edges at +-(L/2 - claw_width))
        bar = draw.rectangle(L - 2 * p.claw_width + 2e-6, aw, 0, 0)
        claw_l, claw_r = self._claw(-1), self._claw(1)
        jj_metal = str(p.jj_metal).lower() in ('true', '1')
        squid_parts = [_clean(squid_island)] if jj_metal else []
        tmon = _clean(draw.union(stem, bar, claw_l, claw_r, *squid_parts))
        jj_pad_ground = _clean(jj_pad_ground)

        ground_edge = hl + jg      # top of the SQUID pocket etch (pre-flip frame)
        if p.squid_ground_tail not in ('', None):
            ground_edge = hl + jg / 2 + p.squid_loop_height / 2 - p.squid_shift_toward_stem + p.squid_ground_tail
        self._squid_etch_edge = ground_edge
        etch = draw.union(
            draw.box(-hw / 2 - p.pad_gap, -jg, hw / 2 + p.pad_gap, ground_edge),
            draw.rectangle(L, aw + 2 * pg, 0, 0),
            claw_l.buffer(p.claw_gap, join_style=2),
            claw_r.buffer(p.claw_gap, join_style=2))

        rect_jj1, rect_jj2, tmon, jj_pad_ground, etch = self._place(
            [rect_jj1, rect_jj2, tmon, jj_pad_ground, etch])

        self.add_qgeometry('poly', dict(tmon=tmon))
        if jj_metal:
            self.add_qgeometry('poly', dict(squid_ground=jj_pad_ground))
        self.add_qgeometry('poly', dict(pad_etch=etch), subtract=True)
        self.add_qgeometry('junction', dict(rect_jj1=rect_jj1, rect_jj2=rect_jj2),
                           width=jw, gds_cell_name=p.gds_cell_name, hfss_inductance=p.Lj)


class ClawCouplerFab(ClawCoupler, MyTmonCouplerFab):
    """Fab mask version: Dolan contact notches of MyTmonCouplerFab (2 x 4 um, island -1 um / ground +1 um in x,
    keeping the +2 um relative Dolan offset while centring the pair on each 10 um SQUID leg)."""
    default_options = Dict(ClawCoupler.default_options,
                           **{**MyTmonCouplerFab.default_options,
                              'island_notch_offset': '-1um', 'ground_notch_offset': '1um'})
