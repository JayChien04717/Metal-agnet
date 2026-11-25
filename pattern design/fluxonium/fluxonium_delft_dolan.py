# -*- coding: utf-8 -*-
r"""Delft (Andersen lab) fluxonium: full lab pattern, NTHU-style JJ + JJ array.

Standalone QComponent (imports only qiskit_metal / numpy). The base metal reproduces the
fluxonium of the lab layout3 FTF mask (Design2_FTF_Andersen Q1, measured from
``Design2_2026-11-18_wafer_r11_reference.gds``): round pads on 15 um stems, 900 x 550 um
pocket, U-loop flux line from the left wall (10 / 8.409 um), pill readout pad (15 / 10 um
line) and an open-ended charge line 30 um left of the centre (15 / 10 um). The lab Manhattan
junction drawing (JJ pads, inductor arm stubs on the stems) is NOT drawn. The loop is the
NTHU hairpin of ``MyFluxoniumDolan._loop`` (fluxonium_nthu_dolan.py), re-implemented here:

Sketch (orientation = 0, the Design2 Q1 frame)::

                               | |  readout_line (pin, +y)
                    ___________| |___________
                   (       readout pad       )
                    -------------------------
         _______________________________________________
        |                        ___                    |
        |                       /   \  pad_top          |
        |                       \___/                   |
        |                        | |   stem             |
  flux  |__________________      |_|                    |
  line ====================|  #===|                     |   #### = JJ array hairpin, end towards the flux tip
  (pin) |                  |  #   x  <- rect_jj         |   x    = main JJ (pocket centre)
        |__________________|  #===|                     |   |    = 2 um lead line
        |                  ^  ^  ___                    |
        |            U tip       | |                    |
        |             array_end  | |                    |
        |                       /   \                   |
        |                       \___/  pad_bot          |
        |_______________________________________________|
                              _
                             | |  charge_line (pin, -y), x = cl_x

* Flux line: signal enters one leg (``signal_leg``, sketch: top) through a channel in the left
  pocket wall; the other leg is shorted to the wall. The U tip tracks the hairpin end using the
  NTHU ``fl_sep`` convention (``x_tip = lead_x - array_length - lead_width/2 - fl_sep``; with
  the 53 um array and 16 um ``fl_sep`` default, the hairpin end sits at -54 um and the U tip
  at -70 um). Setting ``fl_x`` explicitly overrides ``fl_sep`` with a fixed tip coordinate.
* ``rect_jj``: short ``jj_length`` gap in the middle of the 2 um lead line (``L_j``).
* ``array_top`` / ``array_bot``: the two hairpin sides, lumped L (``l_inductance / 2`` each),
  joined by the metal ``array_end``; the hairpin points towards the flux line (-x).
* ``jj_metal='False'``: no leads and no ``array_end`` (junction-free mask). The junction
  table is always written (simulation annotations; e-beam drawing lives in
  ``pattern design/ebeam_junction/``).

Pins: 'flux_line', 'readout_line', 'charge_line'.
"""

import numpy as np
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import BaseQubit


def _true(value):
    return str(value).lower() in ('true', '1')


class FluxoniumDelftDolan(BaseQubit):
    """Lab (Delft) round-pad fluxonium pocket with the NTHU hairpin JJ + array.

    Default Options (defaults = Design2_FTF_Andersen Q1 of the lab mask):
        * pocket_width / pocket_height: ground cut-out (900 x 550 um)
        * pad_radius: radius of the round pads (60 um)
        * stem_width / stem_length: stem from each round pad towards the centre (15 / 100 um);
          pad_overlap: length of the stem inside the disc (disc centre at
          pad_gap/2 + stem_length + pad_radius - pad_overlap = 160 um)
        * pad_gap: gap between the two stem ends (30 um); the leads and JJ sit in it
        * lead_x, lead_width: x and width of the lead line (JJ at lead_x)
        * jj_length, jj_width: main JJ gap and its width
        * array_length: lead line -> centre line of the hairpin end (towards -x)
        * array_sep: arm centre to arm centre; array_width: arm / end width
        * L_j, l_inductance: main JJ and total array inductance (simulation inputs)
        * gds_cell_name, gds_cell_inductor: junction-table GDS cell names
        * jj_metal: 'False' leaves out leads and array_end (junction-free mask)
        * mirror: 'True' mirrors the component in x before rotation
        * flux_line_options: make_fl, signal_leg ('top' / 'bottom': the leg fed through the
          wall channel, the other is shorted to the wall), fl_sep (gap between the array's closed
          end and the U tip face, same as FluxoniumNTHU; optional fl_x overrides with a fixed
          centre -> U tip distance), fl_height (inner spacing of the two legs), cpw_width,
          cpw_gap (wall channel), channel_length (wall -> pin)
        * readout_line_options: make_rol, pad_width (straight part), pad_height (= 2 x end
          radius), pad_sep (pocket edge -> pad), line_length (pad top -> pin), cpw_width, cpw_gap
        * charge_line_options: make_cl, cl_x, cl_sep (ground between pocket edge and the
          gap end), end_gap (gap beyond the tip), cl_length (tip -> pin), cpw_width, cpw_gap
    """

    component_metadata = Dict(short_name='FluxoniumDelftDolan',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True')

    default_options = Dict(
        pocket_width='900um',
        pocket_height='550um',
        pad_radius='60um',
        stem_width='15um',
        stem_length='100um',
        pad_overlap='15um',  # stem end inside the disc (lab: disc centre at pad_gap/2 + 145 um)
        pad_gap='30um',
        # NTHU hairpin (MyFluxoniumDolan defaults): 55 x 10 um outer, 2 um wires and leads
        lead_x='0um',
        lead_width='2um',
        jj_length='2um',
        jj_width='2um',
        array_length='53um',
        array_sep='8um',
        array_width='2um',
        L_j='34.38nH',
        l_inductance='200nH',
        gds_cell_name='my_jj',
        gds_cell_inductor='my_array',
        jj_metal='True',
        mirror='False',
        orientation='0',
        flux_line_options=Dict(
            make_fl=True,
            signal_leg='top',   # 'top' (lab Q1) or 'bottom' (lab Q2): leg fed through the wall channel
            fl_sep='16um',      # gap from array_end to U tip (NTHU convention: 53+1+16 = 70 um at defaults)
            fl_x='',            # optional fixed centre -> U tip distance; '' uses fl_sep
            fl_height='50um',
            cpw_width='10um',
            cpw_gap='8.409um',
            channel_length='40um',
        ),
        readout_line_options=Dict(
            make_rol=True,
            pad_width='400um',
            pad_height='120um',
            pad_sep='25um',
            line_length='100um',
            cpw_width='15um',
            cpw_gap='10um',
        ),
        charge_line_options=Dict(
            make_cl=True,
            cl_x='-30um',
            cl_sep='12.5um',
            end_gap='15um',
            cl_length='150um',
            cpw_width='15um',
            cpw_gap='10um',
        ))

    TOOLTIP = """Lab (Delft) round-pad fluxonium, NTHU hairpin JJ + array."""

    # ------------------------------------------------------------------ helpers
    def _place(self, objs):
        """Mirror (optional), rotate and translate a list of geometries."""
        p = self.p
        if _true(p.mirror):
            objs = draw.scale(objs, -1, 1, origin=(0, 0))
        objs = draw.rotate(objs, p.orientation, origin=(0, 0))
        return draw.translate(objs, p.pos_x, p.pos_y)

    def _pin(self, name, a, b, width):
        """Pin from a -> b; the normal is (b - a) rotated +90 deg, so order a, b to make it
        point outward (unmirrored). The order is swapped when mirrored."""
        if _true(self.p.mirror):
            a, b = b, a
        line = self._place([draw.LineString([a, b])])[0]
        self.add_pin(name, np.array(line.coords), width)

    def make(self):
        self.make_pocket()
        if self.p.flux_line_options.make_fl:
            self.make_flux_line()
        if self.p.readout_line_options.make_rol:
            self.make_readout_line()
        if self.p.charge_line_options.make_cl:
            self.make_charge_line()

    # ------------------------------------------------------------------ qubit
    def _pads(self):
        """Round pad + stem, top and bottom (local, unplaced)."""
        p = self.p
        g2, sl, r = p.pad_gap / 2, p.stem_length, p.pad_radius
        pads = []
        for s in (1, -1):
            stem = draw.rectangle(p.stem_width, sl, 0, s * (g2 + sl / 2))
            disc = draw.Point(0, s * (g2 + sl + r - p.pad_overlap)).buffer(r)
            pads.append(draw.union(stem, disc))
        return pads

    def _loop(self):
        """NTHU hairpin (as MyFluxoniumDolan._loop): lead line with the JJ gap, two lumped-L sides
        towards -x, metal hairpin end (local, unplaced)."""
        p = self.p
        g2 = p.pad_gap / 2
        lw, aw = p.lead_width, p.array_width
        a2, al, jl, x0 = p.array_sep / 2, p.array_length, p.jj_length, p.lead_x
        xe = x0 - al                                           # centre line of the hairpin end
        return dict(
            lead_top=draw.rectangle(lw, g2 - jl / 2, x0, (g2 + jl / 2) / 2),
            lead_bot=draw.rectangle(lw, g2 - jl / 2, x0, -(g2 + jl / 2) / 2),
            array_end=draw.rectangle(lw, 2 * a2 + aw, xe, 0),
            rect_jj=draw.LineString([(x0, jl / 2), (x0, -jl / 2)]),
            array_top=draw.LineString([(x0 - lw / 2, a2), (xe + lw / 2, a2)]),
            array_bot=draw.LineString([(x0 - lw / 2, -a2), (xe + lw / 2, -a2)]))

    def make_pocket(self):
        p = self.p
        jj_metal = _true(p.jj_metal)
        pad_top, pad_bot = self._pads()
        lp = self._loop()
        if jj_metal:
            pad_top = draw.union(pad_top, lp['lead_top'])
            pad_bot = draw.union(pad_bot, lp['lead_bot'])
        rect_pk = draw.rectangle(p.pocket_width, p.pocket_height)
        (pad_top, pad_bot, rect_pk, array_end, rect_jj, array_top,
         array_bot) = self._place([pad_top, pad_bot, rect_pk, lp['array_end'], lp['rect_jj'],
                                   lp['array_top'], lp['array_bot']])

        self.add_qgeometry('poly', dict(pad_top=pad_top, pad_bot=pad_bot))
        if jj_metal:
            self.add_qgeometry('poly', dict(array_end=array_end))
        self.add_qgeometry('poly', dict(rect_pk=rect_pk), subtract=True)
        self.add_qgeometry('junction', dict(rect_jj=rect_jj), width=p.jj_width,
                           hfss_inductance=p.L_j, q3d_inductance=p.L_j,
                           gds_cell_name=p.gds_cell_name)
        l_half = f"{float(str(p.l_inductance).replace('nH', '')) / 2}nH"
        self.add_qgeometry('junction', dict(array_top=array_top, array_bot=array_bot),
                           width=p.array_width, hfss_inductance=l_half, q3d_inductance=l_half,
                           gds_cell_name=p.gds_cell_inductor)

    def _flux_tip_x(self):
        """x of the U tip: fl_sep left of the array's closed end (same convention as FluxoniumNTHU)."""
        p, pf = self.p, self.p.flux_line_options
        if pf.fl_x not in ('', None):
            return -pf.fl_x
        return p.lead_x - p.array_length - p.lead_width / 2 - pf.fl_sep

    # ------------------------------------------------------------------ lines
    def make_flux_line(self):
        """U-loop flux line from the left pocket wall; signal on the signal_leg leg, the other
        leg shorted to the wall; outer corners rounded with radius cpw_width."""
        p, pf = self.p, self.p.flux_line_options
        W, w, g = p.pocket_width / 2, pf.cpw_width, pf.cpw_gap
        h, xt, ch = pf.fl_height / 2, self._flux_tip_x(), pf.channel_length
        xi = xt - w                                            # inner face of the vertical
        eps = 1e-6                                             # 1 nm overlap so the parts merge
        top_leg = draw.box(-W - ch, h, xi + eps, h + w)
        bot_leg = draw.box(-W, -h - w, xi + eps, -h)
        vertical = draw.box(xi, -h - eps, xt, h + eps)
        corners = [draw.Point(xi, s * h).buffer(w).intersection(
            draw.box(xi - eps, min(s * (h - eps), s * (h + w)), xt, max(s * (h - eps), s * (h + w))))
            for s in (1, -1)]
        flux = draw.union(top_leg, bot_leg, vertical, *corners)
        if flux.geom_type != 'Polygon':
            raise ValueError('flux line must be one polygon')
        channel = draw.rectangle(ch, w + 2 * g, -W - ch / 2, h + w / 2)
        a, b = (-W - ch, h), (-W - ch, h + w)
        if pf.signal_leg == 'bottom':                      # signal on the bottom leg
            flux, channel = draw.scale([flux, channel], 1, -1, origin=(0, 0))
            a, b = (-W - ch, -h - w), (-W - ch, -h)
        elif pf.signal_leg != 'top':
            raise ValueError("flux_line_options.signal_leg must be 'top' or 'bottom'")
        flux, channel = self._place([flux, channel])
        self.add_qgeometry('poly', dict(flux_line=flux))
        self.add_qgeometry('poly', dict(flux_line_gap=channel), subtract=True)
        self._pin('flux_line', a, b, w)

    def make_readout_line(self):
        """Pill readout pad above the pocket with a CPW line to the pin (+y)."""
        p, pr = self.p, self.p.readout_line_options
        H, w, g = p.pocket_height / 2, pr.cpw_width, pr.cpw_gap
        pw, ph = pr.pad_width, pr.pad_height
        yc = H + pr.pad_sep + ph / 2
        y_pin = yc + ph / 2 + pr.line_length

        def pill(grow):
            # lab layout: the gap's round ends sit at +/-(pad_width/2 + gap), i.e. the gap is
            # cpw_gap above / below the pad and 2 x cpw_gap at its round ends
            r, xe = ph / 2 + grow, pw / 2 + grow
            return draw.union(draw.rectangle(2 * xe, 2 * r, 0, yc),
                              draw.Point(-xe, yc).buffer(r), draw.Point(xe, yc).buffer(r))

        pad = draw.union(pill(0), draw.rectangle(w, y_pin - yc, 0, (yc + y_pin) / 2))
        gap = draw.union(pill(g), draw.rectangle(w + 2 * g, y_pin - yc, 0, (yc + y_pin) / 2))
        pad, gap = self._place([pad, gap])
        self.add_qgeometry('poly', dict(readout_pad=pad))
        self.add_qgeometry('poly', dict(readout_pad_gap=gap), subtract=True)
        self._pin('readout_line', (-w / 2, y_pin), (w / 2, y_pin), w)

    def make_charge_line(self):
        """Open-ended CPW below the pocket, rounded tip, at x = cl_x (pin -y)."""
        p, pc = self.p, self.p.charge_line_options
        H, w, g, x = p.pocket_height / 2, pc.cpw_width, pc.cpw_gap, pc.cl_x
        r_gap = w / 2 + g
        y_gap_end = -H - pc.cl_sep                              # gap edge closest to the pocket
        y_tip = y_gap_end - pc.end_gap                          # metal tip
        y_pin = y_tip - pc.cl_length
        yc_m, yc_g = y_tip - w / 2, y_gap_end - r_gap           # centres of the round ends
        line = draw.union(draw.rectangle(w, yc_m - y_pin, x, (yc_m + y_pin) / 2),
                          draw.Point(x, yc_m).buffer(w / 2))
        gap = draw.union(draw.rectangle(2 * r_gap, yc_g - y_pin, x, (yc_g + y_pin) / 2),
                         draw.Point(x, yc_g).buffer(r_gap))
        line, gap = self._place([line, gap])
        self.add_qgeometry('poly', dict(charge_line=line))
        self.add_qgeometry('poly', dict(charge_line_gap=gap), subtract=True)
        self._pin('charge_line', (x + w / 2, y_pin), (x - w / 2, y_pin), w)
