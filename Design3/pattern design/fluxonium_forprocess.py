# -*- coding: utf-8 -*-
"""Simulation version of MyFluxonium.

Same pads, pocket, U-loop flux line, readout and charge line as
``MyFluxonium`` (so the EM environment is identical), but the JJ and the
hairpin JJ array are removed. The two leads from the pads meet a short gap in
the middle, bridged by one lumped inductor (junction table) of length
``ind_length``. Use it for Q3D / HFSS / LOM runs.

The U-loop flux line is shorted to ground, so here it is drawn as part of the
ground plane (left out of the etch) instead of a separate metal polygon; Q3D
otherwise reports overlapping nets.

Sketch of the centre (orientation = 0)::

            pad_top
         _____|_____
              |   lead
              |
              #   <- lumped inductor (ind_length, l_inductance)
              |
         _____|_____
            pad_bot
"""

from qiskit_metal import draw, Dict
from my_fluxonium import MyFluxonium


class FluxoniumForProcess(MyFluxonium):
    """MyFluxonium with the JJ + array replaced by one lumped inductor.

    Extra options:
        * ind_length: length of the lumped-inductor gap between the two leads
        * l_inductance: inductance given to that element (hfss / q3d)
    """

    component_metadata = Dict(short_name='FluxoniumSim',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True')

    default_options = Dict(
        ind_length='40um',
    )

    TOOLTIP = """MyFluxonium for simulation: one lumped inductor, no JJ / array."""

    def make_pocket(self):
        p = self.p
        g2 = p.pad_gap / 2
        lw = p.lead_width
        il = p.ind_length
        x0 = p.lead_x

        pad_top = draw.rectangle(p.pad_width, p.pad_height, 0, g2 + p.pad_height / 2)
        pad_bot = draw.rectangle(p.pad_width, p.pad_height, 0, -g2 - p.pad_height / 2)
        lead_top = draw.rectangle(lw, g2 - il / 2, x0, (g2 + il / 2) / 2)
        lead_bot = draw.rectangle(lw, g2 - il / 2, x0, -(g2 + il / 2) / 2)
        pad_top = draw.union(pad_top, lead_top)
        pad_bot = draw.union(pad_bot, lead_bot)

        inductor = draw.LineString([(x0, il / 2), (x0, -il / 2)])

        # same pocket as MyFluxonium
        rect_pk = self._pocket_geom()
        if p.flux_line_options.make_fl:
            rect_pk = draw.subtract(rect_pk, self._flux_line_geom())

        pad_top, pad_bot, rect_pk, inductor = self._place([pad_top, pad_bot, rect_pk, inductor])

        self.add_qgeometry('poly', dict(pad_top=pad_top, pad_bot=pad_bot))
        self.add_qgeometry('poly', dict(rect_pk=rect_pk), subtract=True)
        self.add_qgeometry('junction', dict(inductor=inductor), width=lw,
                           hfss_inductance=p.l_inductance, q3d_inductance=p.l_inductance)

    def make_flux_line(self):
        """Flux U-loop is already part of the ground plane (see make_pocket); pin only."""
        self._flux_pin()
