# -*- coding: utf-8 -*-
r"""Floating double-pad fixed-frequency transmon pocket with a Manhattan junction lead pair.

Same capacitor pads, pocket and readout coupling pad as ``FixedTransmonPocketDolan``;
the straight Dolan landing legs are replaced by a tapered + 45-degree-bent rounded-tip
base-metal lead pair compatible with ``TransmonJunctionManhattan`` (``squid='False'``).

Sketch (orientation = 0)::

                      readout (pin, +y when loc_H = +1)
                              | |
                        +-----+ +-----+   <- readout coupling pad (optional)
         +--------------|-------------|--------------+
         |              +-------------+              |
         |    +---------------------------------+    |
         |    |           pad_top (+y)          |    |  <- rounded corners
         |    +-----------\       /-------------+    |
         |                 \     /                   |  <- tapered upper lead
         |                   \  \                    |  <- 45-deg bent tip (3 um)
         |                   (0,0)  JJ               |  <- Manhattan JJ site
         |                 \  \                      |  <- 45-deg bent tip (180-deg rot)
         |                 /     \                   |  <- tapered lower lead
         |    +-----------/       \-------------+    |
         |    |           pad_bot (-y)          |    |
         |    +---------------------------------+    |
         +-------------------------------------------+
                       pocket (subtract=True)

Pins: one pin per entry in ``connection_pads`` (e.g. ``readout``), outward normal.
Call ``manhattan_sites()`` after construction to get the placed site centre for
``TransmonJunctionManhattan`` / ``add_ebeam_junctions()``.
"""
from qiskit_metal import Dict
from shapely.geometry import Point

from fixed_transmon_pocket_dolan import FixedTransmonPocketDolan
from manhattan_leads import DEFAULT_MANHATTAN_LEAD_OPTIONS, manhattan_lead_pair


class FixedTransmonPocketManhattan(FixedTransmonPocketDolan):
    """Fixed-frequency floating double-pad transmon with one Manhattan junction site."""

    component_metadata = Dict(
        short_name='FTManhattan',
        _qgeometry_table_path='True',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )

    default_options = Dict(
        FixedTransmonPocketDolan.default_options,
        **DEFAULT_MANHATTAN_LEAD_OPTIONS,
        jj_width='3um',
        gds_cell_name='junction_fm0',
    )

    TOOLTIP = """Floating double-pad fixed-frequency transmon pocket with one Manhattan junction site."""

    def _junction_geometry(self):
        """Return (upper_lead, lower_lead, jj) in local coordinates centred at (0, 0)."""
        p = self.p
        unit = self.design.parse_value('1um')
        overlap = min(p.pad_height / 4.0, 5.0 * unit)
        y_root = p.pad_gap / 2.0 + overlap
        return manhattan_lead_pair(y_root, p, x_offset=0.0)

    def manhattan_sites(self):
        """Return placed Manhattan JJ site metadata (in design units, mm) for e-beam placement."""
        p = self.p
        (pt,) = self._place([Point(0.0, 0.0)])
        return dict(
            squid=False,
            center=(pt.x, pt.y),
            orientation=float(p.orientation),
            sites=[dict(name='jj1', center=(pt.x, pt.y), orientation=float(p.orientation))],
        )

    def contact_hooks(self):
        raise NotImplementedError(
            f'{self.name}: Manhattan transmons use manhattan_sites() instead of Dolan contact_hooks()'
        )
