# -*- coding: utf-8 -*-
r"""Circular wide-gap fixed-frequency (single-JJ) transmon with a keyhole readout pad (Design13).

- Circular island (minimal perimeter-to-area, no corners) in a concentric wide-gap pocket.
- **Keyhole readout:** a circular readout pad sits inside a circular keyhole cut into the island at the readout
  diagonal; its CPW stem leaves radially through a slot. The island wraps the pad, which gives a large qubit-readout
  capacitance in a small footprint and shields the pad from the four coupler C-arcs around the island.
- **Manhattan JJ:** one vertical JJ (fingers along x / y) at a free diagonal corner, spanning pad rim <-> pocket rim
  along y. South corners use the lead pattern as drawn, north corners the same pattern rotated by 180 deg (identical
  Manhattan geometry). Never rotated by 45 deg.

ASCII sketch (orientation 0, readout NE, JJ SW; not to scale)::

                       N coupler C-arc
                                   / readout pin (45 deg)
                     .--------.   /
                   /      (o)--' <- keyhole readout pad
       W C-arc    |   island   |    E C-arc
                   \          /
                    |JJ------'
                       S coupler C-arc
"""
from __future__ import annotations

import math

from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent
from shapely.geometry import LineString, Point, box
from shapely.ops import unary_union

from manhattan_leads import DEFAULT_MANHATTAN_LEAD_OPTIONS, manhattan_lead_pair
from metal_width import min_metal_width, require_strip

CORNER_ANG = {'NE': 45.0, 'NW': 135.0, 'SW': -135.0, 'SE': -45.0}


class ScalablePadQubit(QComponent):
    """Circular wide-gap single-island fixed-frequency transmon with keyhole readout and a Manhattan JJ."""

    default_options = Dict(
        pos_x='0mm',
        pos_y='0mm',
        orientation='0',
        pad_diameter='280um',      # circular island diameter
        pocket_gap='46um',         # uniform island -> ground gap
        jj_corner='auto',          # 'NE'/'NW'/'SE'/'SW' or 'auto' (SW unless readout uses it, then SE)
        jj_gap='36um',             # Manhattan JJ gap (along y)
        jj_metal='True',
        Lj='10.8nH',
        make_readout_pad='True',
        ro_corner='NE',
        ro_exit_mode='diag',       # 'diag' (45 deg radial stem) or 'ortho' (+-y stem, unit cell)
        ro_key_pad_radius='13um',  # readout pad disc radius
        ro_key_gap='20um',         # uniform pad <-> island gap
        ro_key_rim='15um',         # island metal rim outside the keyhole
        cpw_width='10um',
        cpw_gap='6um',
        cpw_length='42um',
        **DEFAULT_MANHATTAN_LEAD_OPTIONS,
        layer='1',
        chip='main',
    )

    component_metadata = Dict(
        short_name='CircPadQ',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )

    def make(self):
        p = self.p
        unit = self.design.parse_value('1um')
        r_pad = 0.5 * float(p.pad_diameter)
        r_pkt = r_pad + float(p.pocket_gap)
        ro_corner = str(p.ro_corner).upper()
        cw, cg, clen = float(p.cpw_width), float(p.cpw_gap), float(p.cpw_length)

        pad_metal = Point(0.0, 0.0).buffer(r_pad, resolution=32)
        pocket = Point(0.0, 0.0).buffer(r_pkt, resolution=32)

        # 1. Vertical Manhattan JJ at a diagonal corner (pad side +y in the JJ frame; 180 deg on north corners)
        jj_gap = float(p.jj_gap)
        jj_corner = str(p.jj_corner).upper()
        if jj_corner == 'AUTO':
            jj_corner = 'SW' if ro_corner != 'SW' else 'SE'
        theta_jj = math.radians(CORNER_ANG[jj_corner])
        r_mid = 0.5 * (r_pad + r_pkt)
        sgn_y = 1.0 if math.sin(theta_jj) > 0 else -1.0
        x_jj = r_mid * math.cos(theta_jj)
        y_circ = sgn_y * math.sqrt(max(r_pad * r_pad - x_jj * x_jj, (10.0 * unit) ** 2))
        y_pkt = sgn_y * math.sqrt(max(r_pkt * r_pkt - x_jj * x_jj, (10.0 * unit) ** 2))
        cx, cy = x_jj, 0.5 * (y_circ + y_pkt)
        stem_hw = 8.0 * unit
        pad_stem = box(-stem_hw, 0.5 * jj_gap, stem_hw, abs(cy - y_circ) + 14.0 * unit)
        gnd_stem = box(-stem_hw, -(abs(y_pkt - cy) + 14.0 * unit), stem_hw, -0.5 * jj_gap)
        up_lead, lo_lead, jj_line = manhattan_lead_pair(0.5 * jj_gap + 4.0 * unit, p, x_offset=0.0)
        phi = 0.0 if sgn_y < 0 else 180.0

        def _jj_xf(g):
            return draw.translate(draw.rotate(g, phi, origin=(0, 0)), cx, cy)

        self._jj_phi = phi
        self._jj_site_local = (cx, cy)
        pad_metal = unary_union([pad_metal, _jj_xf(pad_stem)])
        jj_line = _jj_xf(jj_line)
        jj_lead_dict = {'jj_gnd_stem': _jj_xf(gnd_stem)}
        if str(p.jj_metal).lower() == 'true':
            jj_lead_dict['jj_lead_up'] = _jj_xf(up_lead)
            jj_lead_dict['jj_lead_lo'] = _jj_xf(lo_lead)

        # 2. Keyhole readout pad
        ro_metal = None
        pin_pts = None
        if str(p.make_readout_pad).lower() == 'true':
            theta_c = math.radians(CORNER_ANG[ro_corner])
            a_p = float(p.ro_key_pad_radius)
            a_h = a_p + float(p.ro_key_gap)
            rim = float(p.ro_key_rim)
            require_strip(self.name, 'island rim around readout keyhole', rim, min_metal_width(unit))
            r_kc = r_pad - a_h - rim
            if r_kc - a_h < 0.25 * r_pad:
                raise ValueError(f'{self.name}: readout keyhole too large for the island')
            ux, uy = math.cos(theta_c), math.sin(theta_c)
            c_k = (r_kc * ux, r_kc * uy)
            edge = ((r_pad + 2.0 * unit) * ux, (r_pad + 2.0 * unit) * uy)
            pad_metal = pad_metal.difference(unary_union([
                Point(c_k).buffer(a_h, resolution=32),
                LineString([c_k, edge]).buffer(0.5 * cw + cg, cap_style=2),
            ]))
            if str(p.ro_exit_mode).lower() == 'ortho':
                r_out = r_pkt + 4.0 * unit
                p_out = (r_out * ux, r_out * uy)
                sy_dir = +1.0 if 'N' in ro_corner else -1.0
                path = [c_k, p_out, (p_out[0], p_out[1] + sy_dir * clen)]
            else:
                r_end = r_pkt + clen
                path = [c_k, (r_end * ux, r_end * uy)]
            stem_line = LineString(path)
            ro_metal = unary_union([Point(c_k).buffer(a_p, resolution=32),
                                    stem_line.buffer(0.5 * cw, cap_style=2, join_style=2)])
            pocket = unary_union([pocket, stem_line.buffer(0.5 * cw + cg, cap_style=2, join_style=2)])
            pin_pts = [path[-2], path[-1]]

        # 3. Transform and add geometry
        ang = float(p.orientation)
        px, py = float(p.pos_x), float(p.pos_y)

        def xform(g):
            return draw.translate(draw.rotate(g, ang, origin=(0, 0)), px, py)

        self.add_qgeometry('poly', {'pocket': xform(pocket)}, subtract=True, layer=p.layer, chip=p.chip)
        metal_dict = {'pad': xform(pad_metal)}
        if ro_metal is not None:
            metal_dict['readout_pad'] = xform(ro_metal)
        for lk, lg in jj_lead_dict.items():
            metal_dict[lk] = xform(lg)
        self.add_qgeometry('poly', metal_dict, subtract=False, layer=p.layer, chip=p.chip)
        if str(p.jj_metal).lower() == 'true':
            self.add_qgeometry('junction', {'jj': xform(jj_line)}, width=float(p.lead_tip_width),
                               hfss_inductance=p.Lj, layer='50', chip=p.chip)
        if pin_pts is not None:
            p0, p1 = xform(LineString(pin_pts)).coords
            self.add_pin('readout', points=[p0, p1], width=p.cpw_width, gap=p.cpw_gap, input_as_norm=True)

    def manhattan_sites(self):
        ang = float(self.p.orientation)
        px, py = float(self.p.pos_x), float(self.p.pos_y)
        pt = draw.translate(draw.rotate(draw.Point(self._jj_site_local), ang, origin=(0, 0)), px, py)
        return [dict(name=f'{self.name}_jj', center=(float(pt.x), float(pt.y)),
                     orientation=(ang + self._jj_phi) % 360.0)]
