# -*- coding: utf-8 -*-

"""Package-level planar fluxonium module.

This subclass adds routable interfaces around the fluxonium core:
readout, charge drive, and inductive flux bias interfaces.
"""

from qiskit_metal import Dict, draw

from .fluxonium import Fluxonium


class FluxoniumModule(Fluxonium):
    """U-array fluxonium core plus package pins for FTF layouts."""

    component_metadata = Dict(
        short_name='FluxoniumModule',
        _qgeometry_table_path='True',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )

    default_options = Dict(
        **Fluxonium.default_options,
        readout_line_options=Dict(
            make=True,
            pad_sep='24um',
            pad_width='78um',
            pad_height='18um',
            pad_fillet='8um',
            launch_length='95um',
            cpw_width='10um',
            cpw_gap='6um',
            loc_H='+1',
            loc_W='+1',
            launch_axis='vertical',
        ),
        charge_line_options=Dict(
            make=True,
            side='requested',
            pad_sep='82um',
            pad_width='42um',
            pad_height='14um',
            pad_fillet='6um',
            launch_length='95um',
            cpw_width='4um',
            cpw_gap='3um',
            loc_H='-1',
            loc_W='-1',
        ),
        flux_line_options=Dict(
            make=True,
            side='outer',
            line_sep='32um',
            line_height='160um',
            pickup_length='120um',
            launch_length='80um',
            cpw_width='5um',
            cpw_gap='3um',
            fillet='20um',
            loc_W='+1',
        ),
        transmon_coupler_options=Dict(
            make=False,
            pad_sep='24um',
            pad_width='38um',
            pad_height='150um',
            pad_fillet='12um',
            neck_width='18um',
            neck_length='50um',
            cpw_width='12um',
            cpw_gap='6um',
            loc_W='-1',
        ),
    )

    def make(self):
        self.make_fluxonium()
        if self.p.readout_line_options.make:
            self.make_readout_line()
        if self.p.charge_line_options.make:
            self.make_charge_line()
        if self.p.flux_line_options.make:
            self.make_flux_line()
        if self.p.transmon_coupler_options.make:
            self.make_transmon_coupler()

    def _add_capacitive_port(self, name, opts, x_center, y_center, direction):
        """Create a coupling pad with a short CPW launch and an external pin."""
        pad = self._rounded_rectangle(
            opts.pad_width, opts.pad_height, x_center, y_center, opts.pad_fillet
        )
        pad_cut = self._rounded_rectangle(
            opts.pad_width + 2 * opts.cpw_gap,
            opts.pad_height + 2 * opts.cpw_gap,
            x_center,
            y_center,
            opts.pad_fillet + opts.cpw_gap,
        )

        if abs(direction[1]) > 0:
            start_y = y_center + direction[1] * opts.pad_height / 2
            end_y = start_y + direction[1] * opts.launch_length
            launch = draw.LineString([(x_center, start_y), (x_center, end_y)])
            pin = draw.LineString(
                [
                    (x_center - opts.cpw_width / 2, end_y),
                    (x_center + opts.cpw_width / 2, end_y),
                ]
            )
            reverse_pin = direction[1] < 0
        else:
            start_x = x_center + direction[0] * opts.pad_width / 2
            end_x = start_x + direction[0] * opts.launch_length
            launch = draw.LineString([(start_x, y_center), (end_x, y_center)])
            pin = draw.LineString(
                [
                    (end_x, y_center + opts.cpw_width / 2),
                    (end_x, y_center - opts.cpw_width / 2),
                ]
            )
            reverse_pin = direction[0] < 0

        pad, pad_cut, launch, pin = self._transform([pad, pad_cut, launch, pin])
        self.add_qgeometry('poly', {f'{name}_pad': pad})
        self.add_qgeometry('poly', {f'{name}_pad_cut': pad_cut}, subtract=True)
        self.add_qgeometry('path', {f'{name}_launch': launch}, width=opts.cpw_width)
        self.add_qgeometry(
            'path',
            {f'{name}_launch_cut': launch},
            width=opts.cpw_width + 2 * opts.cpw_gap,
            subtract=True,
        )
        pin_coords = list(draw.shapely.geometry.shape(pin).coords)
        self.add_pin(name, pin_coords[::-1] if reverse_pin else pin_coords, opts.cpw_width)

    def make_readout_line(self):
        p = self.p
        opts = p.readout_line_options
        loc_h = float(opts.loc_H)
        loc_w = float(opts.loc_W)
        pad_x = loc_w * (p.pad_gap / 2 + p.pad_width / 2)
        pad_y = p.pad_y + loc_h * (p.pad_height / 2 + opts.pad_sep + opts.pad_height / 2)
        launch_axis = str(opts.launch_axis).strip().lower()
        direction = (loc_w, 0) if launch_axis == 'horizontal' else (0, loc_h)
        self._add_capacitive_port('readout', opts, pad_x, pad_y, direction)

    def make_charge_line(self):
        p = self.p
        opts = p.charge_line_options
        loc_h = float(opts.loc_H)
        loc_w = float(opts.loc_W)
        if p.flux_line_options.make:
            flux_loc_w = self._resolved_flux_loc_w()
            side_policy = str(getattr(opts, 'side', 'requested')).strip().lower()
            if (
                side_policy == 'avoid_flux'
                and abs(loc_w) > 1e-12
                and abs(flux_loc_w) > 1e-12
                and loc_w * flux_loc_w > 0
            ):
                loc_w = -loc_w
        x_center = loc_w * (p.pad_gap / 2 + p.pad_width / 2)
        pocket_center_y = p.pad_y - p.pad_height / 2 - p.array_leg_length / 2
        pocket_edge_y = pocket_center_y + loc_h * p.pocket_height / 2
        pad_y = pocket_edge_y + loc_h * (opts.pad_sep + opts.pad_height / 2)
        self._add_capacitive_port('charge', opts, x_center, pad_y, (0, loc_h))

    def _resolved_flux_loc_w(self):
        p = self.p
        opts = p.flux_line_options
        coupler_opts = p.transmon_coupler_options
        loc_w = float(opts.loc_W)
        coupler_w = float(coupler_opts.loc_W)
        if str(opts.side).strip().lower() == 'outer':
            loc_w = -coupler_w
        return loc_w

    def make_flux_line(self):
        p = self.p
        opts = p.flux_line_options
        loc_w = self._resolved_flux_loc_w()
        _, _, left_node_x, right_node_x = self._node_coordinates()
        coupled_leg_x = right_node_x if loc_w > 0 else left_node_x
        inner_x = coupled_leg_x + loc_w * opts.line_sep
        edge_x = loc_w * p.pocket_width / 2
        outer_x = edge_x + loc_w * opts.launch_length
        pad_bottom_y = p.pad_y - p.pad_height / 2
        y_top = min(
            p.array_top_y - 2 * (opts.cpw_width + opts.cpw_gap),
            pad_bottom_y - opts.cpw_width - 2 * opts.cpw_gap,
        )
        y_bottom = y_top - opts.line_height

        # Andersen-style pocket flux-bias line: a U-shaped CPW segment enters
        # from the pocket edge, runs beside the inductor leg, and returns to an
        # external side pin without touching the fluxonium metal.
        flux_path = draw.LineString(
            [
                (edge_x, y_top),
                (inner_x, y_top),
                (inner_x, y_bottom),
                (edge_x, y_bottom),
                (outer_x, y_bottom),
            ]
        )
        pin = draw.LineString(
            [
                (outer_x, y_bottom + opts.cpw_width / 2),
                (outer_x, y_bottom - opts.cpw_width / 2),
            ]
        )

        flux_path, pin = self._transform([flux_path, pin])
        self.add_qgeometry('path', {'flux_bias_line': flux_path}, width=opts.cpw_width, fillet=opts.fillet)
        self.add_qgeometry(
            'path',
            {'flux_bias_line_cut': flux_path},
            width=opts.cpw_width + 2 * opts.cpw_gap,
            fillet=opts.fillet,
            subtract=True,
        )
        pin_coords = list(draw.shapely.geometry.shape(pin).coords)
        self.add_pin('flux', pin_coords[::-1] if loc_w < 0 else pin_coords, opts.cpw_width)

    def make_transmon_coupler(self):
        p = self.p
        opts = p.transmon_coupler_options
        loc_w = float(opts.loc_W)
        x_pad = loc_w * (p.pad_gap / 2 + p.pad_width + opts.pad_sep + opts.pad_width / 2)
        y_pad = p.pad_y
        x_neck = x_pad + loc_w * (opts.pad_width / 2 + opts.neck_length / 2)

        coupler = self._rounded_rectangle(
            opts.pad_width, opts.pad_height, x_pad, y_pad, opts.pad_fillet
        )
        neck = self._rounded_rectangle(
            opts.neck_length, opts.neck_width, x_neck, y_pad, opts.pad_fillet / 2
        )
        coupler_cut = coupler.union(neck).buffer(opts.cpw_gap, join_style=1)
        metal = coupler.union(neck)
        pin_x = x_pad + loc_w * (opts.pad_width / 2 + opts.neck_length)
        pin = draw.LineString(
            [
                (pin_x, y_pad + opts.cpw_width / 2),
                (pin_x, y_pad - opts.cpw_width / 2),
            ]
        )

        metal, coupler_cut, pin = self._transform([metal, coupler_cut, pin])
        self.add_qgeometry('poly', {'ftf_transmon_coupler': metal})
        self.add_qgeometry('poly', {'ftf_transmon_coupler_cut': coupler_cut}, subtract=True)
        pin_coords = list(draw.shapely.geometry.shape(pin).coords)
        self.add_pin(
            'transmon_coupler',
            pin_coords[::-1] if loc_w < 0 else pin_coords,
            opts.cpw_width,
        )
