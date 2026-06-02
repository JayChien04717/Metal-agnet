# -*- coding: utf-8 -*-

"""Simplified grounded transmon coupler component for FTF layouts."""

from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent


class GroundedTransmonCoupler(QComponent):
    """Compact grounded transmon coupler for fluxonium-transmon-fluxonium cells."""

    component_metadata = Dict(
        short_name='GroundedTmonCoupler',
        _qgeometry_table_path='True',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )

    default_options = Dict(
        island_width='420um',
        island_height='12um',
        coupler_pad_width='90um',
        coupler_pad_height='150um',
        coupler_pad_gap='45um',
        ground_stem_width='12um',
        ground_stem_length='120um',
        junction_length='8um',
        junction_width='4um',
        flux_line_options=Dict(
            make=True,
            loop_width='130um',
            loop_height='90um',
            loop_sep='24um',
            launch_length='220um',
            cpw_width='5um',
            cpw_gap='3um',
            loc_W='-1',
            loc_H='-1',
            fillet='20um',
        ),
        drive_pad_options=Dict(
            make=True,
            pad_width='64um',
            pad_height='22um',
            pad_sep='38um',
            pad_fillet='6um',
            launch_length='120um',
            cpw_width='10um',
            cpw_gap='6um',
            loc_W='-1',
            loc_H='-1',
        ),
        pocket_width='760um',
        pocket_height='520um',
        pocket_fillet='20um',
        pad_fillet='12um',
        Lj='12nH',
        Cj='2fF',
        orientation='0',
    )

    @staticmethod
    def _rounded_rectangle(width, height, cx, cy, radius):
        if 0 < radius < min(width, height) / 2:
            rect = draw.rectangle(width - 2 * radius, height - 2 * radius, cx, cy)
            return rect.buffer(radius, join_style=1)
        return draw.rectangle(width, height, cx, cy)

    def _transform(self, geometries):
        p = self.p
        transformed = draw.rotate(geometries, p.orientation, origin=(0, 0))
        return draw.translate(transformed, p.pos_x, p.pos_y)

    def make(self):
        p = self.p

        t_bar_width = p.island_width + 2 * (p.coupler_pad_gap + p.coupler_pad_width)
        transmon_bar = self._rounded_rectangle(
            t_bar_width, p.island_height, 0, 0, p.pad_fillet
        )
        stem_y = -(p.island_height / 2 + p.ground_stem_length / 2)
        ground_stem = draw.rectangle(p.ground_stem_width, p.ground_stem_length, 0, stem_y)
        transmon_island = transmon_bar.union(ground_stem)
        junction_y = -(p.island_height / 2 + p.junction_length / 2)
        junction = draw.LineString([(0, junction_y + p.junction_length / 2), (0, junction_y - p.junction_length / 2)])
        pocket = self._rounded_rectangle(p.pocket_width, p.pocket_height, 0, -p.ground_stem_length / 4, p.pocket_fillet)

        left_pin = draw.LineString(
            [
                (-t_bar_width / 2, -p.drive_pad_options.cpw_width / 2),
                (-t_bar_width / 2, p.drive_pad_options.cpw_width / 2),
            ]
        )
        right_pin = draw.LineString(
            [
                (t_bar_width / 2, p.drive_pad_options.cpw_width / 2),
                (t_bar_width / 2, -p.drive_pad_options.cpw_width / 2),
            ]
        )

        objects = [transmon_island, ground_stem, junction, pocket, left_pin, right_pin]
        transmon_island, ground_stem, junction, pocket, left_pin, right_pin = self._transform(objects)

        self.add_qgeometry('poly', {'transmon_t_island': transmon_island})
        self.add_qgeometry('poly', {'transmon_coupler_pocket': pocket}, subtract=True)
        self.add_qgeometry(
            'junction',
            {'transmon_ground_junction': junction},
            width=p.junction_width,
            hfss_inductance=p.Lj,
            hfss_capacitance=p.Cj,
        )

        self.add_pin('left_couple', list(draw.shapely.geometry.shape(left_pin).coords), p.drive_pad_options.cpw_width)
        self.add_pin('right_couple', list(draw.shapely.geometry.shape(right_pin).coords), p.drive_pad_options.cpw_width)

        if p.drive_pad_options.make:
            self.make_drive_pad()
        if p.flux_line_options.make:
            self.make_flux_line()

    def make_drive_pad(self):
        p = self.p
        opts = p.drive_pad_options
        loc_w = float(opts.loc_W)
        loc_h = float(opts.loc_H)

        x_center = loc_w * (p.island_width / 2 + opts.pad_sep + opts.pad_width / 2)
        y_center = loc_h * (p.island_height / 2 + opts.pad_sep + opts.pad_height / 2)
        pad = self._rounded_rectangle(opts.pad_width, opts.pad_height, x_center, y_center, opts.pad_fillet)
        pad_cut = self._rounded_rectangle(
            opts.pad_width + 2 * opts.cpw_gap,
            opts.pad_height + 2 * opts.cpw_gap,
            x_center,
            y_center,
            opts.pad_fillet + opts.cpw_gap,
        )
        launch_start_y = y_center + loc_h * opts.pad_height / 2
        launch_end_y = launch_start_y + loc_h * opts.launch_length
        launch = draw.LineString([(x_center, launch_start_y), (x_center, launch_end_y)])
        pin = draw.LineString(
            [
                (x_center - opts.cpw_width / 2, launch_end_y),
                (x_center + opts.cpw_width / 2, launch_end_y),
            ]
        )

        pad, pad_cut, launch, pin = self._transform([pad, pad_cut, launch, pin])
        self.add_qgeometry('poly', {'drive_coupling_pad': pad})
        self.add_qgeometry('poly', {'drive_coupling_pad_cut': pad_cut}, subtract=True)
        self.add_qgeometry('path', {'drive_launch': launch}, width=opts.cpw_width)
        self.add_qgeometry(
            'path',
            {'drive_launch_cut': launch},
            width=opts.cpw_width + 2 * opts.cpw_gap,
            subtract=True,
        )
        pin_coords = list(draw.shapely.geometry.shape(pin).coords)
        self.add_pin('drive', pin_coords[::-1] if loc_h > 0 else pin_coords, opts.cpw_width)

    def make_flux_line(self):
        p = self.p
        opts = p.flux_line_options
        loc_w = float(opts.loc_W)
        loc_h = float(opts.loc_H)
        left_x = -opts.loop_width / 2
        right_x = opts.loop_width / 2
        top_y = -(p.island_height / 2 + opts.loop_sep)
        bottom_y = top_y - opts.loop_height
        launch_x = left_x if loc_w < 0 else right_x
        launch_start_y = bottom_y
        pin_y = launch_start_y + loc_h * opts.launch_length

        # Match the fluxonium U geometry: two vertical pickup legs joined at
        # the bottom, with the launch tail exiting toward bottom ground.
        flux_pickup = draw.LineString(
            [
                (left_x, top_y),
                (left_x, bottom_y),
                (right_x, bottom_y),
                (right_x, top_y),
            ]
        )
        flux_launch = draw.LineString(
            [
                (launch_x, launch_start_y),
                (launch_x, pin_y),
            ]
        )
        pin = draw.LineString(
            [
                (launch_x - opts.cpw_width / 2, pin_y),
                (launch_x + opts.cpw_width / 2, pin_y),
            ]
        )

        flux_pickup, flux_launch, pin = self._transform([flux_pickup, flux_launch, pin])
        flux_paths = {
            'coupler_flux_bias_line': flux_pickup,
            'coupler_flux_bias_launch': flux_launch,
        }
        self.add_qgeometry('path', flux_paths, width=opts.cpw_width, fillet=opts.fillet)
        self.add_qgeometry(
            'path',
            {f'{name}_cut': geom for name, geom in flux_paths.items()},
            width=opts.cpw_width + 2 * opts.cpw_gap,
            fillet=opts.fillet,
            subtract=True,
        )
        pin_coords = list(draw.shapely.geometry.shape(pin).coords)
        self.add_pin('flux', pin_coords[::-1] if loc_h > 0 else pin_coords, opts.cpw_width)
