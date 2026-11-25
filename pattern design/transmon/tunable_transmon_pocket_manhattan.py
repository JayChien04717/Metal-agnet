# -*- coding: utf-8 -*-
r"""Floating double-pad tunable transmon pocket with a 2-JJ Manhattan SQUID.

Same capacitor pads, pocket and readout coupling pad as ``TunableTransmonPocketDolan``;
the two Dolan SQUID legs are replaced by two parallel tapered + 45-degree-bent
Manhattan base-metal lead pairs translated by ``squid_spacing`` (never mirrored),
compatible with ``TransmonJunctionManhattan`` (``squid='True'``).

Sketch (orientation = 0)::

                      readout (pin, +y when loc_H = +1)
                              | |
                        +-----+ +-----+   <- readout coupling pad (optional)
         +--------------|-------------|--------------+
         |              +-------------+              |
         |    +---------------------------------+    |
         |    |           pad_top (+y)          |    |  <- rounded corners
         |    +------\   /---------\   /--------+    |
         |            \ \           \ \              |  <- 2 translated Manhattan
         |            JJ1   spacing JJ2              |     lead pairs (no mirror)
         |           \ \           \ \               |
         |    +------/   \---------/   \--------+    |
         |    |           pad_bot (-y)          |    |
         |    +---------------------------------+    |
         +-------------------------------------------+
                       pocket (subtract=True)

Pins: one pin per entry in ``connection_pads`` (e.g. ``readout``), outward normal.
Call ``manhattan_sites()`` after construction to get placed site centres for
``TransmonJunctionManhattan`` / ``add_ebeam_junctions()``.
"""
from qiskit_metal import Dict, draw
from shapely.geometry import Point

from manhattan_leads import DEFAULT_MANHATTAN_LEAD_OPTIONS, _clean, manhattan_lead_pair
from tunable_transmon_pocket_dolan import TunableTransmonPocketDolan


class TunableTransmonPocketManhattan(TunableTransmonPocketDolan):
    """Floating double-pad tunable transmon pocket with a 2-JJ Manhattan SQUID."""

    component_metadata = Dict(
        short_name='TTManhattan',
        _qgeometry_table_path='True',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )

    default_options = Dict(
        TunableTransmonPocketDolan.default_options,
        **DEFAULT_MANHATTAN_LEAD_OPTIONS,
        squid_spacing='30um',
        jj_width='3um',
        gds_cell_name='junction_tm0',
    )

    TOOLTIP = """Floating double-pad tunable transmon pocket with 2-JJ Manhattan SQUID."""

    def _squid_geometry(self):
        """Return (upper_leads, lower_leads, jj1, jj2) in local coordinates centred at (0, 0)."""
        p = self.p
        unit = self.design.parse_value('1um')
        x_jj = p.squid_spacing / 2.0
        if p.squid_spacing < p.lead_base_width:
            raise ValueError(
                f'{self.name}: squid_spacing ({p.squid_spacing}) must be >= lead_base_width ({p.lead_base_width})'
            )
        overlap = min(p.pad_height / 4.0, 5.0 * unit)
        y_root = p.pad_gap / 2.0 + overlap
        u1, l1, jj1 = manhattan_lead_pair(y_root, p, x_offset=-x_jj)
        u2, l2, jj2 = manhattan_lead_pair(y_root, p, x_offset=+x_jj)
        return _clean(draw.union(u1, u2)), _clean(draw.union(l1, l2)), jj1, jj2

    def manhattan_sites(self):
        """Return placed Manhattan SQUID site metadata (in design units, mm) for e-beam placement."""
        p = self.p
        x_jj = p.squid_spacing / 2.0
        c, p1, p2 = self._place([Point(0.0, 0.0), Point(-x_jj, 0.0), Point(+x_jj, 0.0)])
        ori = float(p.orientation)
        return dict(
            squid=True,
            squid_spacing=p.squid_spacing,
            center=(c.x, c.y),
            orientation=ori,
            sites=[
                dict(name='jj1', center=(p1.x, p1.y), orientation=ori),
                dict(name='jj2', center=(p2.x, p2.y), orientation=ori),
            ],
        )

    def contact_hooks(self):
        raise NotImplementedError(
            f'{self.name}: Manhattan transmons use manhattan_sites() instead of Dolan contact_hooks()'
        )
