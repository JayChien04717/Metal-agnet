# -*- coding: utf-8 -*-
r"""Floating differential tunable coupler: C-arc pads + rounded-rectangle rail-SQUID hub + straight T flux line.

- **C-shaped coupling arcs** at both ends, concentric with the neighbouring circular qubits (uniform coupling gap,
  no corners facing the qubit).
- **10/6 um CPW necks** from the arcs to the hub, so control lines can cross a coupler edge with a standard 90-deg
  airbridge over the neck.
- **Rounded-rectangle hub:** island A (left neck) and island B (right neck) side by side; the long vertical A-B gap
  sets the differential capacitance.
- **Rail SQUID:** island B sends a thin rail under island A's flat bottom edge; two Manhattan JJs bridge A's bottom
  edge and the rail. The SQUID loop is a horizontal strip (squid_spacing x squid_gap) parallel to the flux line.
  JJ fingers stay along x / y: the pattern orientation equals the coupler orientation (0/90/180/270 deg, allowed by
  the Manhattan process); both JJs of the SQUID are translated copies (never mirrored).
- **Straight T flux line** below the hub's flat bottom edge (local -y). It is placed automatically so its left arm
  covers the whole loop and the feed sits just right of the loop; the right arm runs under island B. The two shorted
  arms return the current symmetrically (quadrupole-like far field -> low coupler-to-coupler crosstalk).

ASCII sketch (local frame, orientation 0)::

     qubit L                                                                         qubit R
      (  )  )=========== neck ===========[  A  | |  B  ]=========== neck ===========(  (  )
       C-arc                             [_____| |     ]                              C-arc
                                         [JJ--JJ |     ]   <- rail of B under A
                                         ====T bar======
                                               ||  flux_bias_line (-y)
"""
from __future__ import annotations

import math

from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent
from shapely.geometry import LineString, box
from shapely.ops import unary_union

import fluxline
from fluxline import T_FLUX_OPTIONS
from manhattan_leads import DEFAULT_MANHATTAN_LEAD_OPTIONS, manhattan_lead_pair
from metal_width import min_metal_width, require_strip


def _rrect(x0, y0, x1, y1, r):
    return box(x0, y0, x1, y1).buffer(-r, join_style=1).buffer(r, join_style=1, resolution=16)


class FloatingTunableCoupler(QComponent):
    """C-arc floating coupler with a rounded-rectangle rail-SQUID hub and a straight T flux line."""

    default_options = Dict(
        pos_x='0mm',
        pos_y='0mm',
        orientation='0',          # 0/180: necks along global x; 90/270: necks along global y
        # neighbouring circular qubits (centres at local x = -+qubit_distance)
        qubit_distance='1000um',
        qubit_pad_radius='140um',
        # C-shaped coupling arcs
        c_sep='26um',             # island -> arc gap
        c_width='24um',
        c_span_deg='28',
        c_gap='18um',             # etch around the arc (merges into the qubit pocket)
        neck_width='10um',
        neck_gap='6um',
        # rounded-rectangle hub
        rect_width='230um',
        rect_height='180um',
        rect_corner='24um',       # pocket corner radius
        rect_pocket_gap='12um',   # island edge -> pocket edge
        rect_ab_gap='6um',        # long A-B gap
        rect_island_corner='12um',
        rect_rail_width='8um',    # rail of island B under island A
        squid_gap='20um',         # A bottom edge -> rail (loop height)
        squid_spacing='60um',     # JJ spacing (loop width)
        jj_metal='True',
        Lj='6.5nH',
        flux_line_options=Dict(make_fbl='True', mirror='False', **T_FLUX_OPTIONS),
        flux_as_ground='False',
        **DEFAULT_MANHATTAN_LEAD_OPTIONS,
        layer='1',
        chip='main',
    )

    component_metadata = Dict(
        short_name='FloatTC',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )

    def make(self):
        p = self.p
        unit = self.design.parse_value('1um')
        ang = float(p.orientation)
        nw, ng = float(p.neck_width), float(p.neck_gap)
        W, H = float(p.rect_width), float(p.rect_height)
        gp, gab = float(p.rect_pocket_gap), float(p.rect_ab_gap)
        Ri = float(p.rect_island_corner)
        w_r = float(p.rect_rail_width)
        sq_gap, sq_sp = float(p.squid_gap), float(p.squid_spacing)
        lb = float(p.lead_base_width)
        require_strip(self.name, 'SQUID rail width', w_r, min_metal_width(unit))
        if abs(((ang + 45.0) % 90.0) - 45.0) > 1e-6:
            raise ValueError(f'{self.name}: orientation must be a multiple of 90 deg (Manhattan JJ fingers along x/y)')

        # ---------------- 1. rounded-rectangle hub with a rail SQUID (local frame) ----------------
        x_in, y_in = 0.5 * W - gp, 0.5 * H - gp
        hub_pocket = _rrect(-0.5 * W, -0.5 * H, 0.5 * W, 0.5 * H, float(p.rect_corner))
        y_rail_top = -y_in + w_r
        y_a = y_rail_top + sq_gap                                  # island A's flat bottom edge
        x_s = -(0.5 * gab + 6.0 * unit + 0.5 * lb + 0.5 * sq_sp)   # loop centre, under A next to the A-B gap
        x_left = x_s - 0.5 * sq_sp - 0.5 * lb - 8.0 * unit
        if x_left < -x_in + Ri:
            raise ValueError(f'{self.name}: rect_width too small for squid_spacing {sq_sp}')
        island_a = _rrect(-x_in, y_a, -0.5 * gab, y_in, Ri)
        island_b = unary_union([_rrect(0.5 * gab, -y_in, x_in, y_in, Ri),
                                box(x_left, -y_in, 0.5 * gab + Ri, y_rail_top)])
        if island_a.distance(island_b) < min(gab, sq_gap) - 1e-6:
            raise ValueError(f'{self.name}: islands closer than the A-B gap')
        y_j = y_rail_top + 0.5 * sq_gap
        jj_leads, jj_lines = {}, []
        self._jj_sites_local = []
        for idx, xj in enumerate((x_s - 0.5 * sq_sp, x_s + 0.5 * sq_sp), start=1):
            up, lo, jl = manhattan_lead_pair(0.5 * sq_gap + 2.0 * unit, p, x_offset=xj)
            jj_leads[f'jj{idx}_up'] = draw.translate(up, 0.0, y_j)
            jj_leads[f'jj{idx}_lo'] = draw.translate(lo, 0.0, y_j)
            jj_lines.append((f'jj{idx}', draw.translate(jl, 0.0, y_j)))
            self._jj_sites_local.append((f'{self.name}_jj{idx}', (xj, y_j)))

        # ---------------- 2. C-arcs + necks (local frame) ----------------
        qd = float(p.qubit_distance)
        c_w, c_gap = float(p.c_width), float(p.c_gap)
        c_half = math.radians(0.5 * float(p.c_span_deg))
        r_c_in = float(p.qubit_pad_radius) + float(p.c_sep)
        r_c_mid, r_c_out = r_c_in + 0.5 * c_w, r_c_in + c_w
        pads, pkt = {}, [hub_pocket]
        for side, key, island in ((-1.0, 'pad_l', island_a), (+1.0, 'pad_r', island_b)):
            cx = side * qd
            th0 = 0.0 if side < 0 else math.pi
            angs = [th0 - c_half + 2.0 * c_half * k / 24 for k in range(25)]
            arc = LineString([(cx + r_c_mid * math.cos(a), r_c_mid * math.sin(a)) for a in angs])
            x_out = cx - side * r_c_out
            neck = LineString([(x_out + side * 2.0 * unit, 0.0), (side * (x_in - 4.0 * unit), 0.0)])
            pads[key] = unary_union([arc.buffer(0.5 * c_w, cap_style=1, join_style=1, resolution=12),
                                     neck.buffer(0.5 * nw, cap_style=2, join_style=2), island])
            pkt += [arc.buffer(0.5 * c_w + c_gap, cap_style=1, join_style=1, resolution=12),
                    neck.buffer(0.5 * nw + ng, cap_style=2, join_style=2)]
        if pads['pad_l'].intersects(pads['pad_r']):
            raise ValueError(f'{self.name}: islands A and B touch')
        pocket = unary_union(pkt)

        # ---------------- 3. straight T flux line toward local -y, placed on the loop ----------------
        fbl_opt = p.flux_line_options
        fbl_metal = fbl_gap = fbl_port = None
        if str(fbl_opt.make_fbl).lower() == 'true' and str(p.flux_as_ground).lower() != 'true':
            fbl_x = x_s + 0.5 * sq_sp + 15.0 * unit
            half = min(fbl_x - (x_s - 0.5 * sq_sp - 15.0 * unit), 0.5 * W - float(p.rect_corner) - fbl_x)
            fo = Dict({k: fbl_opt[k] for k in fbl_opt})
            fo['fbl_type'], fo['mirror'] = 'T', 'False'
            fo['fbl_x'], fo['fbl_length'] = fbl_x, 2.0 * half
            d_pre = 0.5 * H
            line_g, gap_g, port_g, _ = fluxline.geometry(d_pre, fo, etch_edge=d_pre, unit=unit, owner=self.name)

            def _flip(g):
                return draw.scale(g, xfact=1.0, yfact=-1.0, origin=(0, 0))

            fbl_metal, fbl_gap, fbl_port = _flip(line_g), _flip(gap_g), _flip(port_g)

        # ---------------- 4. add geometry ----------------
        px, py = float(p.pos_x), float(p.pos_y)

        def xform(g):
            return draw.translate(draw.rotate(g, ang, origin=(0, 0)), px, py)

        sub_dict = {'pocket': xform(pocket)}
        if fbl_gap is not None:
            sub_dict['fbl_gap'] = xform(fbl_gap)
        self.add_qgeometry('poly', sub_dict, subtract=True, layer=p.layer, chip=p.chip)
        met_dict = {'pad_l': xform(pads['pad_l']), 'pad_r': xform(pads['pad_r'])}
        if str(p.jj_metal).lower() == 'true':
            for lk, lg in jj_leads.items():
                met_dict[lk] = xform(lg)
        if fbl_metal is not None:
            met_dict['fbl_metal'] = xform(fbl_metal)
        self.add_qgeometry('poly', met_dict, subtract=False, layer=p.layer, chip=p.chip)
        if str(p.jj_metal).lower() == 'true':
            for jname, jline in jj_lines:
                self.add_qgeometry('junction', {jname: xform(jline)}, width=float(p.lead_tip_width),
                                   hfss_inductance=p.Lj, layer='50', chip=p.chip)
        if fbl_port is not None:
            p_a, p_b = list(xform(fbl_port).coords)[:2]
            mid = (0.5 * (p_a[0] + p_b[0]), 0.5 * (p_a[1] + p_b[1]))
            n_pt, o_pt = xform(draw.Point(0.0, -1.0)), xform(draw.Point(0.0, 0.0))
            nx, ny = n_pt.x - o_pt.x, n_pt.y - o_pt.y
            step = 10.0 * unit
            self.add_pin('flux_bias_line', points=[(mid[0] - step * nx, mid[1] - step * ny), mid],
                         width=fbl_opt.cpw_width, gap=fbl_opt.cpw_gap, input_as_norm=True)

    def manhattan_sites(self):
        """JJ sites in chip coordinates; Manhattan pattern orientation = coupler orientation (multiple of 90 deg)."""
        ang = float(self.p.orientation)
        px, py = float(self.p.pos_x), float(self.p.pos_y)
        out = []
        for nm, xy in self._jj_sites_local:
            pt = draw.translate(draw.rotate(draw.Point(xy), ang, origin=(0, 0)), px, py)
            out.append(dict(name=nm, center=(float(pt.x), float(pt.y)), orientation=ang % 360.0))
        return out
