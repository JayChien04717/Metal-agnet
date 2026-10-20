# -*- coding: utf-8 -*-
r"""My fluxonium: two rectangular pads in a ground pocket.

Sketch (orientation = 0)::

                         readout_line (pin, +y)
                               | |
               _______________ _|_|_ ______________
              |               |_____|              |   <- readout coupling pad
              |      ________________________      |
              |     |        pad_top         |     |
              |     |________________________|     |
   flux  _____|________          |                 |
   line ===============|         |                 |   U-loop flux line from the
   (pin)      |        |   ______|_                |   left wall; bottom leg shorted
              |________|  |######|_|x              |   to ground at the wall; its tip
              |           |######|_|               |   points at the array's end
              |                  |                 |
              |      ________________________      |   #### = JJ array (hairpin, ~50 um)
              |     |        pad_bot         |     |   x    = main JJ, at the pocket centre
              |     |________________________|     |
              |____________________________________|
                               _
                              | |  open-ended charge line
                         charge_line (pin, -y)

The JJ sits at the centre between the pads; the short hairpin array sits next to
it, fully inside the pocket, so there is no ground close to the array. The flux
line is a U-loop entering from the pocket's left wall (pin 'flux_line'); its
bottom leg is shorted to ground at the same wall and its tip sits fl_sep left of
the array's closed end. Set ``mirror='True'`` to put the flux line on the right.

Pins: 'flux_line', 'readout_line', 'charge_line'.
"""

import numpy as np
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import BaseQubit


class MyFluxonium(BaseQubit):
    """Rectangular-pad fluxonium, JJ + short hairpin array, U-loop flux line from the side.

    Default Options:
        * pocket_width / pocket_height: size of the main ground cut-out
        * pocket_top / pocket_bottom: optional distances from the qubit centre to the
          top / bottom pocket edge (asymmetric pocket); '' = pocket_height / 2
        * pad_width / pad_height: size of each capacitor pad
        * pad_gap: distance between the inner edges of the two pads
        * lead_x: x position of the two vertical leads from the pads
          (default -jj_offset, which puts the JJ at the pocket centre)
        * lead_width: width of the thin wires from the pads to the JJ / array
        * jj_length, jj_width: main junction (junction table, gds_cell_name)
        * jj_offset: x distance of the JJ from the leads (away from the flux line)
        * array_length, array_sep, array_width: hairpin array towards the flux line;
          two arms (junction table, gds_cell_inductor, l_inductance/2 each) of
          length array_length, array_sep apart. The loop is mostly the hairpin interior.
        * flux_line_options: make_fl, fl_height (spacing of the two U legs),
          fl_sep (gap between the U tip and the array's end), cpw_width, cpw_gap
        * readout_line_options: make_rol, pad_width, pad_height, pad_sep (gap to
          pad_top), cpw_width, cpw_gap
        * charge_line_options: make_cl, cl_sep (distance of the open end below the
          pocket edge, in the ground plane), cl_length, cpw_width, cpw_gap
    """

    component_metadata = Dict(short_name='MyFluxonium',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True')

    default_options = Dict(
        pocket_width='300um',
        pocket_height='480um',
        pocket_top='',
        pocket_bottom='',
        pad_width='130um',
        pad_height='50um',
        pad_gap='120um',
        lead_x='-12um',
        lead_width='5um',
        jj_length='10um',
        jj_width='7um',
        jj_offset='12um',
        array_length='50um',
        array_sep='32um',
        array_width='7um',
        L_j='34.38nH',
        l_inductance='200nH',
        gds_cell_name='my_jj',
        gds_cell_inductor='my_array',
        mirror='False',
        orientation='0',
        flux_line_options=Dict(
            make_fl=True,
            fl_height='50um',
            fl_sep='20um',
            cpw_width='cpw_width',
            cpw_gap='cpw_gap',
        ),
        readout_line_options=Dict(
            make_rol=True,
            pad_width='200um',
            pad_height='40um',
            pad_sep='30um',
            cpw_width='cpw_width',
            cpw_gap='cpw_gap',
        ),
        charge_line_options=Dict(
            make_cl=True,
            cl_sep='15um',
            cl_length='100um',
            cpw_width='cpw_width',
            cpw_gap='cpw_gap',
        ))

    TOOLTIP = """Rectangular-pad fluxonium, short hairpin array, U-loop flux line."""

    def _place(self, objs):
        """Mirror (optional), rotate and translate a list of geometries."""
        p = self.p
        if str(p.mirror).lower() in ('true', '1'):
            objs = draw.scale(objs, -1, 1, origin=(0, 0))
        objs = draw.rotate(objs, p.orientation, origin=(0, 0))
        return draw.translate(objs, p.pos_x, p.pos_y)

    def _pin(self, name, a, b, width):
        """Pin from a->b; order is swapped when mirrored so the normal stays outward."""
        if str(self.p.mirror).lower() in ('true', '1'):
            a, b = b, a
        line = self._place([draw.LineString([a, b])])[0]
        self.add_pin(name, np.array(line.coords), width)

    def _pocket_geom(self):
        """Ground cut-out in local coordinates (unplaced)."""
        top, bot = self._pk_top(), self._pk_bot()
        return draw.rectangle(self.p.pocket_width, top + bot, 0, (top - bot) / 2)

    def _pk_top(self):
        p = self.p
        return p.pocket_height / 2 if p.pocket_top in ('', None) else p.pocket_top

    def _pk_bot(self):
        p = self.p
        return p.pocket_height / 2 if p.pocket_bottom in ('', None) else p.pocket_bottom

    def make(self):
        self.make_pocket()
        if self.p.flux_line_options.make_fl:
            self.make_flux_line()
        if self.p.readout_line_options.make_rol:
            self.make_readout_line()
        if self.p.charge_line_options.make_cl:
            self.make_charge_line()

    def make_pocket(self):
        p = self.p
        g2 = p.pad_gap / 2
        lw = p.lead_width
        jl = p.jj_length
        al = p.array_length
        a2 = p.array_sep / 2
        jo = p.jj_offset
        x0 = p.lead_x

        pad_top = draw.rectangle(p.pad_width, p.pad_height, 0, g2 + p.pad_height / 2)
        pad_bot = draw.rectangle(p.pad_width, p.pad_height, 0, -g2 - p.pad_height / 2)

        # vertical leads from each pad down/up to the array arms (x = lead_x)
        lead_top = draw.rectangle(lw, g2 - a2 + lw / 2, x0, (g2 + a2 - lw / 2) / 2)
        lead_bot = draw.rectangle(lw, g2 - a2 + lw / 2, x0, -(g2 + a2 - lw / 2) / 2)
        # JJ branch away from the flux line: short arms to x0 + jj_offset, stubs to the JJ
        jj_arm_top = draw.rectangle(jo + lw, lw, x0 + jo / 2, a2)
        jj_arm_bot = draw.rectangle(jo + lw, lw, x0 + jo / 2, -a2)
        jj_stub_top = draw.rectangle(lw, a2 - jl / 2, x0 + jo, (a2 + jl / 2) / 2)
        jj_stub_bot = draw.rectangle(lw, a2 - jl / 2, x0 + jo, -(a2 + jl / 2) / 2)
        pad_top = draw.union(pad_top, lead_top, jj_arm_top, jj_stub_top)
        pad_bot = draw.union(pad_bot, lead_bot, jj_arm_bot, jj_stub_bot)
        # closing piece at the far end of the hairpin array
        array_end = draw.rectangle(lw, 2 * a2 + lw, x0 - al, 0)

        rect_jj = draw.LineString([(x0 + jo, jl / 2), (x0 + jo, -jl / 2)])
        # hairpin array: two arms, drawn in the same direction (parallel junctions)
        array_top = draw.LineString([(x0, a2), (x0 - al, a2)])
        array_bot = draw.LineString([(x0, -a2), (x0 - al, -a2)])

        rect_pk = self._pocket_geom()

        pad_top, pad_bot, array_end, rect_pk, rect_jj, array_top, array_bot = self._place(
            [pad_top, pad_bot, array_end, rect_pk, rect_jj, array_top, array_bot])

        self.add_qgeometry('poly', dict(pad_top=pad_top, pad_bot=pad_bot, array_end=array_end))
        self.add_qgeometry('poly', dict(rect_pk=rect_pk), subtract=True)
        self.add_qgeometry('junction', dict(rect_jj=rect_jj), width=p.jj_width,
                           hfss_inductance=p.L_j, gds_cell_name=p.gds_cell_name)
        l_half = f"{float(str(p.l_inductance).replace('nH', '')) / 2}nH"
        self.add_qgeometry('junction', dict(array_top=array_top, array_bot=array_bot),
                           width=p.array_width, hfss_inductance=l_half,
                           gds_cell_name=p.gds_cell_inductor)

    def _flux_tip_x(self):
        """x of the U tip: fl_sep left of the array's closed end."""
        p = self.p
        return p.lead_x - p.array_length - p.lead_width / 2 - p.flux_line_options.fl_sep

    def _flux_line_geom(self):
        """U-loop polygon in local coordinates (unplaced): two legs from the left wall
        at y = +/- fl_height/2, joined at the tip next to the array."""
        p = self.p
        pf = p.flux_line_options
        w = pf.cpw_width
        h2 = pf.fl_height / 2
        x_wall = -p.pocket_width / 2
        x_tip = self._flux_tip_x()
        length = x_tip - x_wall
        top_leg = draw.rectangle(length, w, (x_wall + x_tip) / 2, h2)
        bot_leg = draw.rectangle(length + 1e-3, w, (x_wall + x_tip) / 2 - 5e-4, -h2)
        tip = draw.rectangle(w, 2 * h2 + w, x_tip - w / 2, 0)
        return draw.union(top_leg, bot_leg, tip)

    def _flux_pin(self):
        p = self.p
        w = p.flux_line_options.cpw_width
        h2 = p.flux_line_options.fl_height / 2
        x_wall = -p.pocket_width / 2
        self._pin('flux_line', (x_wall, h2 - w / 2), (x_wall, h2 + w / 2), w)

    def make_flux_line(self):
        """U-loop from the left pocket wall: top leg = pin, bottom leg shorted to
        ground at the wall, tip fl_sep away from the array's closed end."""
        flux_line, = self._place([self._flux_line_geom()])
        self.add_qgeometry('poly', dict(flux_line=flux_line))
        self._flux_pin()

    def make_readout_line(self):
        """Coupling pad above pad_top with a trace to the top pocket edge."""
        p = self.p
        pr = p.readout_line_options
        y_pad = p.pad_gap / 2 + p.pad_height + pr.pad_sep + pr.pad_height / 2
        y_edge = self._pk_top()
        pad = draw.rectangle(pr.pad_width, pr.pad_height, 0, y_pad)
        y0 = y_pad + pr.pad_height / 2
        trace = draw.rectangle(pr.cpw_width, y_edge - y0, 0, (y_edge + y0) / 2)
        readout, = self._place([draw.union(pad, trace)])
        self.add_qgeometry('poly', dict(readout_line=readout))
        self._pin('readout_line', (-pr.cpw_width / 2, y_edge), (pr.cpw_width / 2, y_edge), pr.cpw_width)

    def make_charge_line(self):
        """Open-ended CPW that stops in the ground plane, cl_sep below the pocket
        (a thin ground strip stays between its end gap and the pocket -> weak coupling)."""
        p = self.p
        pc = p.charge_line_options
        w, g = pc.cpw_width, pc.cpw_gap
        y_tip = -(self._pk_bot() + pc.cl_sep)
        y_end = y_tip - pc.cl_length
        stub = draw.rectangle(w, pc.cl_length, 0, (y_tip + y_end) / 2)
        gap = draw.rectangle(w + 2 * g, pc.cl_length + g, 0, (y_tip + g + y_end) / 2)
        stub, gap = self._place([stub, gap])
        self.add_qgeometry('poly', dict(charge_line=stub))
        self.add_qgeometry('poly', dict(charge_line_gap=gap), subtract=True)
        self._pin('charge_line', (w / 2, y_end), (-w / 2, y_end), w)
