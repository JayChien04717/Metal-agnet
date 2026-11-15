# -*- coding: utf-8 -*-
"""Flux-tunable Xmon qubit for Design6, based on qiskit-metal ``TransmonCross``.

Changes w.r.t. TransmonCross:
* South arm ends in a SQUID with the Design3 / Design3_2 junction structure (``my_tmon.py`` block: rectangular loop
  20 x 24 um, 40 um block, 5 um cut, two parallel Dolan JJs, electrodes overlapping their base metal by 3 um and
  shifted 2 um toward the arm; ``squid.fab_contacts='True'``: 2 x 4 um contact notches, ground side offset +2 um in x,
  as ``my_tmon_fab.MyTmonCouplerFab``) to ground, with an L-shaped flux line under it (shorted at one end, pin at the other).
  ``Lj_a`` / ``Lj_b``: the two junction inductances (asymmetric SQUID) for HFSS / pyEPR.
  ``squid=Dict(make='False')`` falls back to a single junction ``Lj_a`` across the cross gap (fixed frequency).
* ``jj_metal='False'``: no SQUID electrodes (junction-free mask, gds_export.md); the etch window stays.
* ``flux_as_ground='True'``: the shorted flux line is drawn as part of the ground (Q3D, skill.md 8.2).
* Readout claw: a U claw that wraps the arm tip at ``claw_sep`` from the arm (same shape as the coupler claws of
  ``claw_coupler.py``) instead of the TransmonCross claw that sits outside the cross etch (Q3D 2026-11-14: only
  1.8 fF / g ~ 30 MHz). The claw back plate continues as a CPW stub of ``cpw_length`` ending in the pin.
  connection_pads[name]: connector_location ('0' west, '90' north, '180' east, '270' south, as TransmonCross),
  claw_length, claw_width, claw_sep, claw_gap (etch around the claw), cpw_width, cpw_gap, cpw_length.
"""
from qiskit_metal import Dict, draw
from fluxline import B2_FLUX_OPTIONS, geometry as flux_geometry
from qiskit_metal.qlibrary.qubits.transmon_cross import TransmonCross
from shapely.geometry import MultiPolygon

LOC_ROT = {'0': 0.0, '90': -90.0, '180': 180.0, '270': 90.0}   # local claw drawn on the west arm (-x)


def _clean(geom, min_area=1e-10):
    parts = [g for g in getattr(geom, 'geoms', [geom]) if g.area > min_area]
    return parts[0] if len(parts) == 1 else MultiPolygon(parts)


def _true(v):
    return str(v).lower() in ('true', '1')


class Xmon(TransmonCross):
    default_options = Dict(
        TransmonCross.default_options, Lj_a='10nH', Lj_b='20nH', jj_width='7um', jj_metal='True',
        flux_as_ground='False',
        squid=Dict(make='True', jj_gap='40um', squid_width='40um', loop_width='20um', loop_height='24um',
                   jj_cut='5um', contact_overlap='3um', shift_toward_arm='2um',      # Design3 TC SQUID values
                   fab_contacts='False', contact_notch_width='2um', contact_notch_depth='4um',
                   ground_notch_offset='2um'),
        flux_bias_line_options=Dict(B2_FLUX_OPTIONS, make_fbl='True', mirror='False'),
        _default_connection_pads=Dict(connector_location='0', claw_length='60um', claw_width='10um',
                                      claw_sep='8um', claw_gap='12um', cpw_width='16um', cpw_gap='9um',
                                      cpw_length='40um'))
    component_metadata = Dict(TransmonCross.component_metadata, short_name='Xmon')

    def _place(self, objs, flip=True):
        """Local frame -> chip. flip: drawn with the south arm along +y (my_tmon convention), mirrored to -y."""
        p = self.p
        if flip:
            objs = draw.scale(objs, 1, -1, origin=(0, 0))
        return draw.translate(draw.rotate(objs, p.orientation, origin=(0, 0)), p.pos_x, p.pos_y)

    def make(self):
        self.make_pocket()
        self.make_connection_pads()
        if _true(self.p.squid.make) and _true(self.p.flux_bias_line_options.make_fbl):
            self.make_flux_bias_line()

    def make_pocket(self):
        p = self.p
        cl, cw, cg = p.cross_length, p.cross_width, p.cross_gap
        cross_line = draw.shapely.ops.unary_union([
            draw.LineString([(0, cl), (0, -cl)]), draw.LineString([(cl, 0), (-cl, 0)])])
        cross = cross_line.buffer(cw / 2, cap_style=2)                 # symmetric: frame independent
        cross_etch = cross.buffer(cg, cap_style=3, join_style=2)
        sq = p.squid
        if not _true(sq.make):                                         # fixed-frequency: one junction
            rect_jj = draw.LineString([(0, -cl), (0, -cl - cg)])
            cross, cross_etch, rect_jj = self._place([cross, cross_etch, rect_jj], flip=False)
            self.add_qgeometry('poly', dict(cross=cross), chip=p.chip)
            self.add_qgeometry('poly', dict(cross_etch=cross_etch), subtract=True, chip=p.chip)
            self.add_qgeometry('junction', dict(rect_jj=rect_jj), width=p.jj_width, chip=p.chip,
                               hfss_inductance=p.Lj_a)
            return
        # SQUID block between the south arm tip (y = cl in the flipped frame) and ground (y = cl + jj_gap)
        jg, ov, sw = sq.jj_gap, sq.contact_overlap, sq.squid_width
        lw, lh, cut = sq.loop_width, sq.loop_height, sq.jj_cut
        yc = cl + jg / 2
        block = draw.rectangle(sw, jg + 2 * ov, 0, yc)
        hole = draw.rectangle(lw, lh, 0, yc)
        squid = draw.subtract(block, draw.union(hole, draw.rectangle(sw, cut, 0, yc)))
        island = draw.subtract(squid, draw.box(-sw, yc, sw, cl + jg + ov))
        ground = draw.subtract(squid, draw.box(-sw, cl - ov, sw, yc))
        x_jj = (lw / 2 + sw / 2) / 2
        jj1 = draw.LineString([(-x_jj, yc - cut / 2), (-x_jj, yc + cut / 2)])
        jj2 = draw.LineString([(x_jj, yc - cut / 2), (x_jj, yc + cut / 2)])
        island, ground, jj1, jj2 = draw.translate([island, ground, jj1, jj2], 0, -sq.shift_toward_arm)
        if _true(sq.fab_contacts):    # Dolan contact recesses, as Design3 MyTmonCouplerFab (ground side offset in x)
            hw_n, dep, off = sq.contact_notch_width / 2, sq.contact_notch_depth, sq.ground_notch_offset
            for jj in (jj1, jj2):
                (x, island_edge), (_, ground_edge) = jj.coords
                island = island.difference(draw.box(x - hw_n, island_edge - dep, x + hw_n, island_edge))
                ground = ground.difference(draw.box(x - hw_n + off, ground_edge, x + hw_n + off, ground_edge + dep))
        window = draw.box(-(cw / 2 + cg), 0, cw / 2 + cg, cl + jg)     # etch down to the SQUID ground edge
        # everything below in the flipped frame (the cross and its etch are symmetric under y -> -y)
        jj_metal = _true(p.jj_metal)
        if jj_metal:
            cross = _clean(draw.union(cross, _clean(island)))
        etch = draw.union(cross_etch, window)
        cross, etch, ground, jj1, jj2 = self._place([cross, etch, _clean(ground), jj1, jj2])
        self.add_qgeometry('poly', dict(cross=cross), chip=p.chip)
        self.add_qgeometry('poly', dict(cross_etch=etch), subtract=True, chip=p.chip)
        if jj_metal:
            self.add_qgeometry('poly', dict(squid_ground=ground), chip=p.chip)
        self.add_qgeometry('junction', dict(rect_jj1=jj1), width=p.jj_width, chip=p.chip, hfss_inductance=p.Lj_a)
        self.add_qgeometry('junction', dict(rect_jj2=jj2), width=p.jj_width, chip=p.chip, hfss_inductance=p.Lj_b)

    def make_flux_bias_line(self):
        """Build the shared local flux line under the SQUID."""
        p = self.p
        pf = p.flux_bias_line_options
        d = p.cross_length + p.squid.jj_gap
        line, gap, port, _ = flux_geometry(d, pf)
        line, gap, port = self._place([line, gap, port])
        if _true(p.flux_as_ground):
            self.add_qgeometry('poly', dict(flux_bias_line_gap=draw.subtract(gap, line)), subtract=True)
        else:
            self.add_qgeometry('poly', dict(flux_bias_line=line))
            self.add_qgeometry('poly', dict(flux_bias_line_gap=gap), subtract=True)
        self.add_pin('flux_bias_line', list(port.coords), pf.cpw_width, gap=pf.cpw_gap)

    def make_connection_pad(self, name):
        p = self.p
        pc = p.connection_pads[name]
        arm, aw = p.cross_length, p.cross_width
        cl, cw, cs, cg = pc.claw_length, pc.claw_width, pc.claw_sep, pc.claw_gap
        w, g, L = pc.cpw_width, pc.cpw_gap, pc.cpw_length
        x_in = -(arm + cs)                    # back plate inner edge (faces the arm tip)
        x_out = x_in - cw
        y_in, y_out = aw / 2 + cs, aw / 2 + cs + cw
        back = draw.box(x_out, -y_out, x_in, y_out)
        f_top = draw.box(x_in, y_in, x_in + cl, y_out)
        f_bot = draw.box(x_in, -y_out, x_in + cl, -y_in)
        stub = draw.box(x_out - L, -w / 2, x_out + 1e-6, w / 2)
        claw = draw.union(back, f_top, f_bot, stub)
        etch = draw.union(draw.union(back, f_top, f_bot).buffer(cg, join_style=2),
                          draw.box(x_out - L, -w / 2 - g, x_out, w / 2 + g))
        port = draw.LineString([(x_out - L, -w / 2), (x_out - L, w / 2)])     # normal points outward (-x)
        rot = LOC_ROT[str(int(float(pc.connector_location)))]
        items = draw.rotate([claw, etch, port], rot, origin=(0, 0))
        claw, etch, port = self._place(items, flip=False)
        self.add_qgeometry('poly', {f'{name}_connector_arm': claw}, chip=p.chip)
        self.add_qgeometry('poly', {f'{name}_connector_etcher': etch}, subtract=True, chip=p.chip)
        self.add_pin(name, port.coords, w, gap=g)
