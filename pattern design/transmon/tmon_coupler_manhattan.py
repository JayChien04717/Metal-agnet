# -*- coding: utf-8 -*-
r"""Grounded tunable transmon coupler with flat coupling paddles and a Manhattan SQUID.

Subclasses ``TmonCouplerDolanNthu`` and replaces only ``_squid_geometry()``:
the stem, paddles, flux bias line, charge line, readout pad, and 4 um ground-strip
DRC checks are all inherited from the base class.

Sketch (orientation = 0, stem points to -y)::

     ____                                   ____
    |    |_________________________________|    |   <- paddles + bar (arms)
    |    |_______________     _____________|    |
    |____|               |   |             |____|
                         |   |  stem (head)
                      ___|___|___   island crossbar
                       \ \   \ \    2 translated Manhattan lead pairs (no mirror)
                       JJ1   JJ2
                      \ \   \ \
                      |_________|   ground crossbar -> chip ground
                       ______|      L-shaped flux line (shorted at its left end)
                         | |
                    flux_bias_line (pin, -y)

Pins: 'flux_bias_line', 'charge_line', 'readout_line', 'fake_flux_bias_line' (LOM only).
Call ``manhattan_sites()`` after construction to get placed site centres for
``TransmonJunctionManhattan`` / ``add_ebeam_junctions()``.
"""
from qiskit_metal import Dict, draw
from shapely.geometry import Point

from manhattan_leads import DEFAULT_MANHATTAN_LEAD_OPTIONS, _clean, manhattan_lead_pair
from metal_width import min_metal_width, require_strip
from tmon_coupler_dolan_nthu import TmonCouplerDolanNthu


class TmonCouplerManhattan(TmonCouplerDolanNthu):
    """Grounded tunable transmon coupler with flat paddles and a 2-JJ Manhattan SQUID."""

    component_metadata = Dict(
        short_name='TmonManhattan',
        _qgeometry_table_path='True',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )

    default_options = Dict(
        TmonCouplerDolanNthu.default_options,
        **DEFAULT_MANHATTAN_LEAD_OPTIONS,
        squid_spacing='30um',
        squid_bar_height='6um',
        jj_gap='48um',
        jj_width='3um',
        squid_ground_tail='',
        gds_cell_name='junction_tmc0',
    )

    TOOLTIP = """Grounded tunable transmon coupler with flat paddles and a Manhattan SQUID."""

    def _squid_geometry(self):
        """Return (squid_island, jj_pad_ground, hole, rect_jj1, rect_jj2) in pre-flip local coords.

        ``TmonCouplerDolanNthu._place()`` flips ``y -> -y`` before rotating/translating.
        Each Manhattan lead pair is pre-flipped here so the placed leads on chip keep the
        canonical orientation (upper-right pad at +45 deg, fingers along x and y).
        """
        p = self.p
        unit = self.design.parse_value('1um')
        min_w = min_metal_width(unit)
        require_strip(self.name, 'Manhattan SQUID crossbar height', p.squid_bar_height, min_w)

        hl, jg, overlap = p.pad_head_length, p.jj_gap, p.squid_contact_overlap
        sp = p.squid_spacing
        bar_h = p.squid_bar_height
        if sp < p.lead_base_width:
            raise ValueError(f'{self.name}: squid_spacing ({sp}) must be >= lead_base_width ({p.lead_base_width})')

        yc = hl + jg / 2.0 - p.squid_shift_toward_stem
        # Available half-span between the inner faces of the island and ground crossbars
        y_root = jg / 2.0 - bar_h / 2.0
        if y_root <= p.lead_elbow_y:
            raise ValueError(
                f'{self.name}: jj_gap ({jg}) too small for Manhattan leads '
                f'(need jj_gap/2 - squid_bar_height/2 > lead_elbow_y)'
            )

        x_jj = sp / 2.0
        bar_w = max(sp + p.lead_base_width, p.pad_head_width)
        # Pre-flip each Manhattan lead pair so _place()'s y-flip restores the canonical +45 deg orientation
        u1, l1, jj1 = manhattan_lead_pair(y_root, p, x_offset=-x_jj)
        u2, l2, jj2 = manhattan_lead_pair(y_root, p, x_offset=+x_jj)
        u1, l1, jj1, u2, l2, jj2 = draw.scale([u1, l1, jj1, u2, l2, jj2], 1, -1, origin=(0, 0))
        u1, l1, jj1, u2, l2, jj2 = draw.translate([u1, l1, jj1, u2, l2, jj2], 0.0, yc)

        # Island crossbar (at y = hl - overlap .. yc - y_root) joins stem to upper leads
        isl_y0 = min(hl - overlap, yc - y_root - bar_h)
        isl_y1 = yc - y_root + 0.5 * unit
        island_bar = draw.box(-bar_w / 2.0, isl_y0, bar_w / 2.0, isl_y1)
        squid_island = _clean(draw.union(island_bar, u1, u2))

        # Ground crossbar (at y = yc + y_root .. hl + jg + overlap) joins lower leads to ground
        gnd_y0 = yc + y_root - 0.5 * unit
        gnd_y1 = max(hl + jg + overlap, yc + y_root + bar_h)
        ground_bar = draw.box(-bar_w / 2.0, gnd_y0, bar_w / 2.0, gnd_y1)
        jj_pad_ground = _clean(draw.union(ground_bar, l1, l2))

        hole = draw.rectangle(max(sp - p.lead_base_width, unit), 2.0 * y_root, 0.0, yc)
        return squid_island, jj_pad_ground, hole, jj1, jj2

    def manhattan_sites(self):
        """Return placed Manhattan SQUID site metadata (in design units, mm) for e-beam placement."""
        p = self.p
        yc = p.pad_head_length + p.jj_gap / 2.0 - p.squid_shift_toward_stem
        x_jj = p.squid_spacing / 2.0
        c, p1, p2 = self._place([Point(0.0, yc), Point(-x_jj, yc), Point(+x_jj, yc)])
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
