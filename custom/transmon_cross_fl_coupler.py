# -*- coding: utf-8 -*-

"""TransmonCrossFL-based FTF coupler with independently sized X/Y arms."""

from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.qubits.transmon_cross_fl import TransmonCrossFL


class TransmonCrossFLCoupler(TransmonCrossFL):
    """TransmonCrossFL variant for direct fluxonium-transmon-fluxonium coupling.

    This keeps the Metal cross-flux-line component model, but lets the horizontal
    and vertical cross arms be sized independently.  For FTF layouts, make the
    horizontal arm long enough to provide the left/right transmon coupling metal
    directly, instead of adding fluxonium-side coupling pads and bridges.
    """

    component_metadata = Dict(
        short_name='TransmonCrossFLCoupler',
        _qgeometry_table_poly='True',
        _qgeometry_table_path='True',
        _qgeometry_table_junction='True',
    )

    default_options = Dict(
        cross_length_x='220um',
        cross_length_y='110um',
        make_fl=True,
        fl_options=Dict(
            t_top='45um',
            t_offset='0um',
            t_length='180um',
            t_inductive_gap='3um',
            t_width='5um',
            t_gap='3um',
        ),
        charge_line_options=Dict(
            make=True,
            side='right',
            line_sep='24um',
            length='180um',
            t_width='5um',
            t_gap='3um',
        ),
        t_coupler_options=Dict(
            make=True,
            head_width='28um',
            head_height='95um',
        ),
        connection_pads=Dict(),
    )

    def make(self):
        """Build the cross, side coupling stubs, flux line, and charge line."""
        self.make_pocket()
        self.make_connection_pads()
        if self.p.make_fl:
            self.make_flux_line()
        if self.p.charge_line_options.make:
            self.make_charge_line()

    def make_pocket(self):
        """Make the cross island using independent horizontal and vertical arms."""
        p = self.p
        cross_width = p.cross_width
        cross_gap = p.cross_gap
        cross_length_x = p.cross_length_x
        cross_length_y = p.cross_length_y
        chip = p.chip

        cross_line = draw.shapely.ops.unary_union([
            draw.LineString([(0, cross_length_y), (0, -cross_length_y)]),
            draw.LineString([(cross_length_x, 0), (-cross_length_x, 0)]),
        ])
        cross = cross_line.buffer(cross_width / 2, cap_style=2)
        tc = p.t_coupler_options
        if tc.make:
            head_left = draw.rectangle(
                tc.head_width,
                tc.head_height,
                -cross_length_x - tc.head_width / 2,
                0,
            )
            head_right = draw.rectangle(
                tc.head_width,
                tc.head_height,
                cross_length_x + tc.head_width / 2,
                0,
            )
            cross = draw.shapely.ops.unary_union([cross, head_left, head_right])
        cross_etch = cross.buffer(cross_gap, cap_style=3, join_style=2)
        rect_jj = draw.LineString([
            (0, -cross_length_y),
            (0, -cross_length_y - cross_gap),
        ])

        cross, cross_etch, rect_jj = self._transform([cross, cross_etch, rect_jj])

        self.add_qgeometry('poly', {'cross': cross}, chip=chip)
        self.add_qgeometry('poly', {'cross_etch': cross_etch}, subtract=True, chip=chip)
        self.add_qgeometry('junction', {'rect_jj': rect_jj}, width=cross_width, chip=chip)

    def make_connection_pad(self, name: str):
        """Make a connector pad using the matching X or Y arm length."""
        p = self.p
        cross_width = p.cross_width
        cross_gap = p.cross_gap
        chip = p.chip

        pc = self.p.connection_pads[name]
        c_g = pc.claw_gap
        c_l = pc.claw_length
        c_w = pc.claw_width
        g_s = pc.ground_spacing
        con_loc = float(pc.connector_location)

        claw_cpw = draw.box(0, -c_w / 2, -4 * c_w, c_w / 2)

        if int(float(pc.connector_type)) == 0:
            t_claw_height = 2 * c_g + 2 * c_w + 2 * g_s + 2 * cross_gap + cross_width
            claw_base = draw.box(-c_w, -t_claw_height / 2, c_l, t_claw_height / 2)
            claw_subtract = draw.box(0, -t_claw_height / 2 + c_w, c_l, t_claw_height / 2 - c_w)
            claw_base = claw_base.difference(claw_subtract)
            connector_arm = draw.shapely.ops.unary_union([claw_base, claw_cpw])
            connector_etcher = draw.buffer(connector_arm, c_g)
        else:
            connector_arm = claw_cpw
            connector_etcher = draw.buffer(connector_arm, c_g)

        port_line = draw.LineString([(-4 * c_w, -c_w / 2), (-4 * c_w, c_w / 2)])

        claw_rotate = 0
        arm_length = p.cross_length_x
        if con_loc > 135:
            claw_rotate = 180
            arm_length = p.cross_length_x
        elif con_loc > 45:
            claw_rotate = -90
            arm_length = p.cross_length_y

        polys = [connector_arm, connector_etcher, port_line]
        polys = draw.translate(polys, -(arm_length + cross_gap + g_s + c_g), 0)
        polys = draw.rotate(polys, claw_rotate, origin=(0, 0))
        polys = self._transform(polys)
        connector_arm, connector_etcher, port_line = polys

        self.add_qgeometry('poly', {f'{name}_connector_arm': connector_arm}, chip=chip)
        self.add_qgeometry(
            'poly',
            {f'{name}_connector_etcher': connector_etcher},
            subtract=True,
            chip=chip,
        )
        self.add_pin(name, port_line.coords, c_w)

    def make_flux_line(self):
        """Make the TransmonCrossFL T flux line relative to the Y arm."""
        p = self.p
        pf = self.p.fl_options
        t_length = getattr(pf, 't_length', 0.03)

        h_line = draw.LineString([(-pf.t_top / 2, 0), (pf.t_top / 2, 0)])
        v_line = draw.LineString([(pf.t_offset, 0), (pf.t_offset, -t_length)])
        parts = [h_line, v_line]

        parts = draw.translate(
            parts,
            0,
            -(p.cross_length_y + p.cross_gap + pf.t_inductive_gap + pf.t_width / 2 + pf.t_gap),
        )
        h_line, v_line = self._transform(parts)

        self.add_qgeometry(
            'path',
            {'h_line': h_line, 'v_line': v_line},
            width=pf.t_width,
            layer=p.layer,
        )
        self.add_qgeometry(
            'path',
            {'h_line_sub': h_line, 'v_line_sub': v_line},
            width=pf.t_width + 2 * pf.t_gap,
            subtract=True,
            layer=p.layer,
        )
        self.add_pin('flux', points=v_line.coords, width=pf.t_width, gap=pf.t_gap, input_as_norm=True)

    def make_charge_line(self):
        """Make a drive/charge line next to the transmon flux line."""
        p = self.p
        pf = self.p.fl_options
        pc = self.p.charge_line_options
        side = -1 if str(pc.side).strip().lower() == 'left' else 1
        x_offset = pf.t_offset + side * (pf.t_top / 2 + pc.line_sep)
        y_top = -(p.cross_length_y + p.cross_gap + pf.t_inductive_gap + pf.t_width / 2 + pf.t_gap)
        y_bottom = y_top - pc.length

        charge_line = draw.LineString([(x_offset, y_top), (x_offset, y_bottom)])
        charge_line = self._transform(charge_line)

        self.add_qgeometry(
            'path',
            {'charge_line': charge_line},
            width=pc.t_width,
            layer=p.layer,
        )
        self.add_qgeometry(
            'path',
            {'charge_line_sub': charge_line},
            width=pc.t_width + 2 * pc.t_gap,
            subtract=True,
            layer=p.layer,
        )
        self.add_pin('drive', points=charge_line.coords, width=pc.t_width, gap=pc.t_gap, input_as_norm=True)
        self.add_pin('charge', points=charge_line.coords, width=pc.t_width, gap=pc.t_gap, input_as_norm=True)

    def _transform(self, geometries):
        p = self.p
        geometries = draw.rotate(geometries, p.orientation, origin=(0, 0))
        return draw.translate(geometries, p.pos_x, p.pos_y)
