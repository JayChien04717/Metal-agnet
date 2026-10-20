# -*- coding: utf-8 -*-
r"""Grounded tunable transmon coupler with flat coupling paddles.

Same SQUID / flux line as ``Tmon_c`` but the arm ends are flat rectangular
paddles (no half-moon palms), so each paddle faces the straight side edge of a
fluxonium pad.

Sketch (orientation = 0)::

     ____                                   ____
    |    |_________________________________|    |   <- paddles + bar (arms)
    |    |_______________     _____________|    |
    |____|               |   |             |____|
                         |   |  stem (head)
                        _|___|_
                       |x /  \ x|  SQUID: loop hole, 2 parallel JJs (Dolan)
                       |__\__/__|  ground side
                        ______|  L-shaped flux line (shorted at its left end)
                          | |
                     flux_bias_line (pin, -y)

Pins: 'flux_bias_line', 'charge_line', 'readout_line', 'fake_flux_bias_line' (LOM only).
"""

from shapely.geometry import MultiPolygon
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import BaseQubit


def _clean(geom, min_area=1e-10):
    """Drop zero-area slivers left by boolean ops (Ansys cannot draw them)."""
    parts = [g for g in getattr(geom, 'geoms', [geom]) if g.area > min_area]
    return parts[0] if len(parts) == 1 else MultiPolygon(parts)


class MyTmonCoupler(BaseQubit):
    """Grounded tunable transmon coupler with flat paddles.

    Default Options:
        * pad_arm_length / pad_arm_width: the horizontal bar between the paddles
        * paddle_width / paddle_height: flat coupling paddle at each end
        * paddle_gap: etch around the bar and the paddles
        * pad_head_width / pad_head_length: the stem from the bar down to the SQUID
        * pad_gap: etch around the stem
        * jj_gap: height of the SQUID block between stem end and ground
        * squid_width, squid_loop_width, squid_loop_height: SQUID block and its hexagonal loop hole
        * jj_cut, jj_width: cut between the two electrodes and junction width; the
          two junctions are drawn parallel (same direction) for Dolan bridges
        * Lj: inductance of each SQUID junction in HFSS (lumped element)
        * flux_as_ground: 'True' draws the (shorted) flux line as part of the ground
          plane, for Q3D
        * flux_bias_line_options: make_fbl, fbl_sep, fbl_height, fbl_length, cpw_width, cpw_gap
        * readout_line_options: make_rol, ro_width / ro_height (coupling pad parallel to
          the bar), ro_sep (pad-bar gap), ro_gap (etch around the pad), cpw_width, cpw_gap
        * charge_line_options: make_cl, cl_x (x offset from the stem), cl_sep (ground strip
          between its open end and the bar etch), cl_length, cpw_width, cpw_gap
    """

    component_metadata = Dict(short_name='MyTmon',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True')

    default_options = Dict(
        chip='main',
        jj_width='7um',
        jj_gap='40um',
        squid_width='40um',
        squid_loop_width='16um',
        squid_loop_height='24um',
        jj_cut='5um',
        pad_head_width='15um',
        pad_head_length='400um',
        pad_arm_width='25um',
        pad_arm_length='840um',
        paddle_width='30um',
        paddle_height='110um',
        paddle_gap='30um',
        pad_gap='30um',
        gds_cell_name='junction_t0',
        Lj='10nH',
        flux_as_ground='False',
        orientation='0',
        flux_bias_line_options=Dict(
            make_fbl=True,
            fbl_sep='10um',
            fbl_height='80um',
            fbl_length='100um',
            cpw_width='cpw_width',
            cpw_gap='cpw_gap',
        ),
        charge_line_options=Dict(
            make_cl=True,
            cl_x='200um',
            cl_sep='15um',
            cl_length='100um',
            cpw_width='cpw_width',
            cpw_gap='cpw_gap',
        ),
        readout_line_options=Dict(
            make_rol=False,
            ro_width='300um',
            ro_height='20um',
            ro_sep='5um',
            ro_gap='20um',
            cpw_width='cpw_width',
            cpw_gap='cpw_gap',
        ))

    TOOLTIP = """Grounded tunable transmon coupler with flat paddles."""

    def _place(self, objs):
        """Geometry is drawn with the stem along +y (Tmon_c style), then flipped
        so the stem points to -y, rotated and translated."""
        p = self.p
        objs = draw.scale(objs, 1, -1, origin=(0, 0))
        objs = draw.rotate(objs, p.orientation, origin=(0, 0))
        return draw.translate(objs, p.pos_x, p.pos_y)

    def make(self):
        self.make_pocket()
        if self.p.flux_bias_line_options.make_fbl:
            self.make_flux_bias_line()
        if self.p.charge_line_options.make_cl:
            self.make_charge_line()
        if self.p.readout_line_options.make_rol:
            self.make_readout_line()

    def make_readout_line(self):
        """Readout coupling pad parallel to the bar, on the side away from the stem,
        ro_sep from the bar, in its own etch window (ro_gap around it), with a trace
        to the 'readout_line' pin at the window edge."""
        p = self.p
        pr = p.readout_line_options
        w, rw, rh = pr.cpw_width, pr.ro_width, pr.ro_height
        # drawn with the stem along +y (see _place), so 'away from the stem' is -y
        y_pad = -(p.pad_arm_width / 2 + pr.ro_sep + rh / 2)
        y_top = y_pad - rh / 2 - pr.ro_gap                      # outer edge of the etch window
        pad = draw.rectangle(rw, rh, 0, y_pad)
        trace = draw.rectangle(w, y_pad - rh / 2 - y_top, 0, (y_pad - rh / 2 + y_top) / 2)
        window = draw.rectangle(rw + 2 * pr.ro_gap, p.pad_arm_width / 2 + pr.ro_sep + rh + pr.ro_gap,
                                0, (y_top + 0) / 2)
        port = draw.LineString([(-w / 2, y_top), (w / 2, y_top)])
        ro, window, port = self._place([draw.union(pad, trace), window, port])
        self.add_qgeometry('poly', dict(readout_line=ro))
        self.add_qgeometry('poly', dict(readout_window=window), subtract=True)
        self.add_pin('readout_line', list(port.coords), w)

    def make_charge_line(self):
        """Open-ended CPW below the bar (stem side) at x = cl_x. Its end stops in the
        ground plane cl_sep away from the bar etch -> weak capacitive drive."""
        p = self.p
        pc = p.charge_line_options
        w, g, x = pc.cpw_width, pc.cpw_gap, pc.cl_x
        y_tip = p.pad_arm_width / 2 + p.paddle_gap + pc.cl_sep
        y_end = y_tip + pc.cl_length
        stub = draw.rectangle(w, pc.cl_length, x, (y_tip + y_end) / 2)
        gap = draw.rectangle(w + 2 * g, pc.cl_length + g, x, (y_tip - g + y_end) / 2)
        port = draw.LineString([(x + w / 2, y_end), (x - w / 2, y_end)])
        stub, gap, port = self._place([stub, gap, port])
        self.add_qgeometry('poly', dict(charge_line=stub))
        self.add_qgeometry('poly', dict(charge_line_gap=gap), subtract=True)
        self.add_pin('charge_line', list(port.coords), w)

    def make_pocket(self):
        p = self.p
        hw, hl = p.pad_head_width, p.pad_head_length
        jg, jw = p.jj_gap, p.jj_width
        L, aw = p.pad_arm_length, p.pad_arm_width
        pw, ph, pg = p.paddle_width, p.paddle_height, p.paddle_gap

        # SQUID (Dolan style): a block between the stem end (y = hl) and ground
        # (y = hl + jg) with a hexagonal loop hole and a horizontal cut. The two
        # junctions bridge the cut on the left and right arms. Both junction lines
        # point the same way, so the junctions are parallel (same Dolan bridge
        # orientation), not mirrored.
        sw, lw, lh, cut = p.squid_width, p.squid_loop_width, p.squid_loop_height, p.jj_cut
        yc = hl + jg / 2
        block = draw.rectangle(sw, jg, 0, yc)
        hole = draw.Polygon([(0, yc + lh / 2), (lw / 2, yc + lh / 4), (lw / 2, yc - lh / 4),
                             (0, yc - lh / 2), (-lw / 2, yc - lh / 4), (-lw / 2, yc + lh / 4)])
        cut_strip = draw.rectangle(sw, cut, 0, yc)
        squid = draw.subtract(block, draw.union(hole, cut_strip))
        squid_island = draw.subtract(squid, draw.box(-sw, yc, sw, hl + jg))   # joins the stem
        jj_pad_ground = draw.subtract(squid, draw.box(-sw, hl, sw, yc))      # joins ground
        x_jj = (lw / 2 + sw / 2) / 2                                         # middle of each arm
        rect_jj1 = draw.LineString([(-x_jj, yc - cut / 2), (-x_jj, yc + cut / 2)])
        rect_jj2 = draw.LineString([(x_jj, yc - cut / 2), (x_jj, yc + cut / 2)])

        stem = draw.rectangle(hw, hl, 0, hl / 2)
        bar = draw.rectangle(L, aw, 0, 0)
        paddle_l = draw.rectangle(pw, ph, -(L / 2 + pw / 2), 0)
        paddle_r = draw.rectangle(pw, ph, L / 2 + pw / 2, 0)
        tmon = _clean(draw.union(stem, bar, paddle_l, paddle_r, _clean(squid_island)))
        jj_pad_ground = _clean(jj_pad_ground)

        etch = draw.union(
            draw.rectangle(hw + 2 * p.pad_gap, hl + 2 * jg, 0, hl / 2),
            draw.rectangle(L, aw + 2 * pg, 0, 0),
            draw.rectangle(pw + 2 * pg, ph + 2 * pg, -(L / 2 + pw / 2), 0),
            draw.rectangle(pw + 2 * pg, ph + 2 * pg, L / 2 + pw / 2, 0))

        rect_jj1, rect_jj2, tmon, jj_pad_ground, etch = self._place(
            [rect_jj1, rect_jj2, tmon, jj_pad_ground, etch])

        self.add_qgeometry('poly', dict(tmon=tmon, squid_ground=jj_pad_ground))
        self.add_qgeometry('poly', dict(pad_etch=etch), subtract=True)
        self.add_qgeometry('junction', dict(rect_jj1=rect_jj1, rect_jj2=rect_jj2),
                           width=jw, gds_cell_name=p.gds_cell_name, hfss_inductance=p.Lj)

    def make_flux_bias_line(self):
        """L-shaped flux line under the SQUID: comes up on the +x side (pin), runs along
        under the SQUID and is shorted straight into the ground at its -x end (no
        second leg bending back down)."""
        p = self.p
        pf = p.flux_bias_line_options
        s, fl, fh, w, g = pf.fbl_sep, pf.fbl_length, pf.fbl_height, pf.cpw_width, pf.cpw_gap
        d = p.pad_head_length + p.jj_gap

        y_end = d + 2 * (fh + s)                       # pin end of the leg
        x0, x1 = -fl / 2 - w, fl / 2 + w               # bar under the SQUID; leg at x = fl/2 .. fl/2 + w
        line = draw.union(draw.box(x0, d + s, x1, d + s + w),          # bar
                          draw.box(fl / 2, d + s, x1, y_end))          # leg to the pin
        # gaps: above the bar, below the bar, both sides of the leg; nothing at the
        # bar's -x end, so that end merges into the ground (short)
        gap = draw.union(draw.box(x0, d + s - g, x1 + g, d + s),
                         draw.box(x0, d + s + w, fl / 2, d + s + w + g),
                         draw.box(fl / 2 - g, d + s + w, fl / 2, y_end),
                         draw.box(x1, d + s - g, x1 + g, y_end))
        # pin order reversed w.r.t. Tmon_c because of the y-flip in _place
        port = draw.LineString([(x1, y_end), (fl / 2, y_end)])
        fake = draw.LineString([(x0, d + s + w), (x0, d + s)])        # shorted end (LOM only)

        line, gap, port, fake = self._place([line, gap, port, fake])
        if str(p.flux_as_ground).lower() in ('true', '1'):
            # simulation: the line is shorted to ground anyway -> leave it in the ground plane
            self.add_qgeometry('poly', dict(flux_bias_line_gap=draw.subtract(gap, line)), subtract=True)
        else:
            self.add_qgeometry('poly', dict(flux_bias_line=line))
            self.add_qgeometry('poly', dict(flux_bias_line_gap=gap), subtract=True)
        self.add_pin('flux_bias_line', list(port.coords), w)
        self.add_pin('fake_flux_bias_line', list(fake.coords), w)
