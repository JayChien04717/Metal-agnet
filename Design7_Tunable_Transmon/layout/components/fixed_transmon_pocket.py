# -*- coding: utf-8 -*-
r"""Floating double-pad fixed-frequency transmon pocket with one Josephson junction between the pads.

Design-specific component (``Design7_Tunable_Transmon/layout/components/``). Same capacitor pads, pocket and readout
coupling pad as ``TunableTransmonPocket``; the two-leg SQUID is replaced by a single centred landing leg pair and one
junction, so there is no flux loop and the qubit frequency is fixed by EJ.

Sketch (orientation = 0)::

                      readout (pin, +y when loc_H = +1)
                              | |
                        +-----+ +-----+   <- readout coupling pad (optional)
         +--------------|-------------|--------------+
         |              +-------------+              |
         |    +---------------------------------+    |
         |    |           pad_top (+y)          |    |  <- rounded corners
         |    +---------------+  +--------------+    |
         |                    |  |                   |  <- upper landing leg (jj_leg_width)
         |                    |__|                   |
         |                     JJ                    |  <- jj_cut (5 um for TransmonJunction, squid=False)
         |                    |__|                   |
         |                    |  |                   |  <- lower landing leg
         |    +---------------+  +--------------+    |
         |    |           pad_bot (-y)          |    |
         |    +---------------------------------+    |
         +-------------------------------------------+
                       pocket (subtract=True)

Pins: one pin per entry in ``connection_pads`` (e.g. ``readout``), outward normal.
Simulation geometry only (Stage 1): Dolan contact hooks / ``contact_hooks()`` for ``TransmonJunction`` are not drawn yet.
"""

from qiskit_metal import Dict, draw
from metal_width import min_metal_width, require_strip
from tunable_transmon_pocket import TunableTransmonPocket, _clean, _rounded_rect, _is_true


class FixedTransmonPocket(TunableTransmonPocket):
    """Fixed-frequency floating double-pad transmon: one junction on the pocket centre line.

    Options (in addition to the pad / pocket / connection-pad options of ``TunableTransmonPocket``):
        * jj_leg_width: width of the upper and lower landing legs (>= 4 um metal rule)
        * jj_cut: vertical gap between the landing faces (5 um for ``TransmonJunction``)
        * jj_width: lumped-junction sheet width in the QGeometry ``junction`` table
        * Lj: HFSS lumped junction inductance (input, from the quantum-model EJ)
        * jj_metal: ``'False'`` omits the landing legs (junction-free base pads only)
    The inherited ``squid_*`` / ``Lj_a`` / ``Lj_b`` / ``fab_contacts`` options are ignored.
    """

    component_metadata = Dict(
        short_name='FTPocket',
        _qgeometry_table_path='True',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )

    default_options = Dict(
        TunableTransmonPocket.default_options,
        jj_leg_width='16um',
        Lj='15nH',
        gds_cell_name='junction_f0',
    )

    TOOLTIP = """Floating double-pad fixed-frequency transmon pocket with one centred junction."""

    def _junction_geometry(self):
        """Return (upper_leg, lower_leg, jj) in local coordinates centred at (0, 0)."""
        p = self.p
        unit = self.design.parse_value('1um')
        require_strip(self.name, 'junction landing leg width', p.jj_leg_width, min_metal_width(unit))
        cut_half = p.jj_cut / 2.0
        overlap = min(p.pad_height / 4.0, 5.0 * unit)   # clean overlap into pad_top / pad_bot
        leg_h = (p.pad_gap / 2.0 + overlap) - cut_half
        if leg_h <= 0:
            raise ValueError(f'{self.name}: pad_gap ({p.pad_gap}) must exceed jj_cut ({p.jj_cut})')
        upper = draw.rectangle(p.jj_leg_width, leg_h, 0.0, +(cut_half + leg_h / 2.0))
        lower = draw.rectangle(p.jj_leg_width, leg_h, 0.0, -(cut_half + leg_h / 2.0))
        # Junction line runs from the upper electrode (island) to the lower electrode, like the SQUID branches
        jj = draw.LineString([(0.0, +cut_half), (0.0, -cut_half)])
        return upper, lower, jj

    def make_pocket(self):
        p = self.p
        pw, ph, pg = p.pad_width, p.pad_height, p.pad_gap
        yc_pad = pg / 2.0 + ph / 2.0

        pad_top = _rounded_rect(pw, ph, p.pad_corner_radius, 0.0, +yc_pad)
        pad_bot = _rounded_rect(pw, ph, p.pad_corner_radius, 0.0, -yc_pad)
        pocket = _rounded_rect(p.pocket_width, p.pocket_height, p.pocket_corner_radius, 0.0, 0.0)

        upper, lower, jj = self._junction_geometry()
        if _is_true(p.jj_metal):
            pad_top = _clean(draw.union(pad_top, upper))
            pad_bot = _clean(draw.union(pad_bot, lower))

        pad_top, pad_bot, pocket, jj = self._place([pad_top, pad_bot, pocket, jj])

        self.add_qgeometry('poly', dict(pad_top=pad_top, pad_bot=pad_bot), chip=p.chip, layer=p.layer)
        self.add_qgeometry('poly', dict(pocket=pocket), subtract=True, chip=p.chip, layer=p.layer)
        self.add_qgeometry('junction', dict(rect_jj=jj), width=p.jj_width,
                           gds_cell_name=p.gds_cell_name, hfss_inductance=p.Lj, chip=p.chip)

    def contact_hooks(self):
        raise NotImplementedError(f'{self.name}: Dolan contact hooks are a Stage 3 (fab) addition, not drawn yet')
