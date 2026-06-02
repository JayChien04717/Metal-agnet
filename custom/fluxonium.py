# -*- coding: utf-8 -*-

"""Planar fluxonium core component for Qiskit Metal.

This geometry follows the common optical-micrograph layout: two capacitor pads
at the top, a small junction bridging the two nodes, and a U-shaped junction
array superinductor that descends from the top nodes and returns at the bottom.
"""

from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import BaseQubit


class Fluxonium(BaseQubit):
    """Fluxonium core with a top junction and U-shaped array superinductor."""

    component_metadata = Dict(
        short_name='Fluxonium',
        _qgeometry_table_path='True',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )

    default_options = Dict(
        pad_width='92um',
        pad_height='64um',
        pad_gap='42um',
        pad_fillet='12um',
        pad_y='36um',
        pocket_width='360um',
        pocket_height='360um',
        pocket_fillet='18um',
        trace_width='3um',
        trace_gap='3um',
        lead_width='6um',
        lead_gap='4um',
        junction_gap='8um',
        junction_y='36um',
        junction_lead_length='18um',
        array_leg_gap='24um',
        array_leg_length='215um',
        array_top_y='18um',
        array_bottom_fillet='6um',
        array_mode='lumped',
        array_lumped_gap='6um',
        array_junctions_per_leg='30',
        array_island_width='9um',
        array_island_height='3.2um',
        array_junction_length='2.4um',
        array_junction_width='4.2um',
        L_j='40nH',
        C_j='0.7fF',
        l_inductance='220nH',
        orientation='0',
    )

    def make(self):
        self.make_fluxonium()

    @staticmethod
    def _rounded_rectangle(width, height, cx, cy, radius):
        if 0 < radius < min(width, height) / 2:
            rect = draw.rectangle(width - 2 * radius, height - 2 * radius, cx, cy)
            return rect.buffer(radius, join_style=1)
        return draw.rectangle(width, height, cx, cy)

    @staticmethod
    def _as_nh(value):
        text = str(value).strip().lower()
        if text.endswith('nh'):
            return float(text[:-2])
        if text.endswith('uh'):
            return float(text[:-2]) * 1000
        if text.endswith('h'):
            return float(text[:-1]) * 1e9
        return float(value)

    def _transform(self, geometries):
        p = self.p
        transformed = draw.rotate(geometries, p.orientation, origin=(0, 0))
        return draw.translate(transformed, p.pos_x, p.pos_y)

    def _node_coordinates(self):
        p = self.p
        left_pad_x = -(p.pad_gap / 2 + p.pad_width / 2)
        right_pad_x = p.pad_gap / 2 + p.pad_width / 2
        left_node_x = -p.array_leg_gap / 2
        right_node_x = p.array_leg_gap / 2
        return left_pad_x, right_pad_x, left_node_x, right_node_x

    def make_fluxonium(self):
        """Build the pads, top junction, and U-shaped junction array."""
        p = self.p
        left_pad_x, right_pad_x, left_node_x, right_node_x = self._node_coordinates()
        top_y = p.array_top_y
        bottom_y = p.array_top_y - p.array_leg_length
        junction_y = p.junction_y

        left_pad = self._rounded_rectangle(
            p.pad_width, p.pad_height, left_pad_x, p.pad_y, p.pad_fillet
        )
        right_pad = self._rounded_rectangle(
            p.pad_width, p.pad_height, right_pad_x, p.pad_y, p.pad_fillet
        )
        pocket = self._rounded_rectangle(
            p.pocket_width,
            p.pocket_height,
            0,
            p.pad_y - p.pad_height / 2 - p.array_leg_length / 2,
            p.pocket_fillet,
        )

        left_pad_inner = left_pad_x + p.pad_width / 2
        right_pad_inner = right_pad_x - p.pad_width / 2
        top_paths = {
            'left_pad_to_node': draw.LineString(
                [(left_pad_inner, junction_y), (left_node_x, junction_y)]
            ),
            'right_pad_to_node': draw.LineString(
                [(right_node_x, junction_y), (right_pad_inner, junction_y)]
            ),
            'left_array_feed': draw.LineString([(left_node_x, junction_y), (left_node_x, top_y)]),
            'right_array_feed': draw.LineString([(right_node_x, junction_y), (right_node_x, top_y)]),
        }
        top_junction = draw.LineString([(left_node_x, junction_y), (right_node_x, junction_y)])

        array_mode = str(p.array_mode).strip().lower()
        bottom_return_y = bottom_y - p.array_bottom_fillet
        array_islands = {}
        array_junctions = {}
        if array_mode == 'explicit':
            bottom_return = draw.LineString(
                [
                    (left_node_x, bottom_y),
                    (left_node_x, bottom_return_y),
                    (right_node_x, bottom_return_y),
                    (right_node_x, bottom_y),
                ]
            )
            n_per_leg = max(1, int(float(p.array_junctions_per_leg)))
            total_array_junctions = 2 * n_per_leg
            junction_inductance = f'{self._as_nh(p.l_inductance) / total_array_junctions}nH'
            pitch = (top_y - bottom_y) / n_per_leg

            for leg_name, x_pos in [('left', left_node_x), ('right', right_node_x)]:
                for idx in range(n_per_leg):
                    island_y = top_y - (idx + 0.5) * pitch
                    array_islands[f'{leg_name}_array_island_{idx:02d}'] = draw.rectangle(
                        p.array_island_width,
                        p.array_island_height,
                        x_pos,
                        island_y,
                    )
                    j_y = top_y - (idx + 1) * pitch
                    array_junctions[f'{leg_name}_array_jj_{idx:02d}'] = draw.LineString(
                        [
                            (x_pos, j_y + p.array_junction_length / 2),
                            (x_pos, j_y - p.array_junction_length / 2),
                        ]
                    )
        else:
            half_gap = p.array_lumped_gap / 2
            bottom_return = {
                'left_array_leg': draw.LineString([(left_node_x, top_y), (left_node_x, bottom_y)]),
                'right_array_leg': draw.LineString([(right_node_x, top_y), (right_node_x, bottom_y)]),
                'left_bottom_return': draw.LineString(
                    [
                        (left_node_x, bottom_y),
                        (left_node_x, bottom_return_y),
                        (-half_gap, bottom_return_y),
                    ]
                ),
                'right_bottom_return': draw.LineString(
                    [
                        (half_gap, bottom_return_y),
                        (right_node_x, bottom_return_y),
                        (right_node_x, bottom_y),
                    ]
                ),
            }
            array_junctions['array_lumped_inductor'] = draw.LineString(
                [(-half_gap, bottom_return_y), (half_gap, bottom_return_y)]
            )
            junction_inductance = p.l_inductance

        objects = [
            left_pad,
            right_pad,
            pocket,
            *top_paths.values(),
            *(bottom_return.values() if isinstance(bottom_return, dict) else [bottom_return]),
            top_junction,
            *array_islands.values(),
            *array_junctions.values(),
        ]
        transformed = self._transform(objects)
        left_pad, right_pad, pocket = transformed[:3]
        offset = 3
        top_paths = dict(zip(top_paths.keys(), transformed[offset:offset + len(top_paths)]))
        offset += len(top_paths)
        if isinstance(bottom_return, dict):
            count = len(bottom_return)
            bottom_return = dict(zip(bottom_return.keys(), transformed[offset:offset + count]))
            offset += count
        else:
            bottom_return = transformed[offset]
            offset += 1
        top_junction = transformed[offset]
        offset += 1
        array_islands = dict(
            zip(array_islands.keys(), transformed[offset:offset + len(array_islands)])
        )
        offset += len(array_islands)
        array_junctions = dict(zip(array_junctions.keys(), transformed[offset:]))

        bottom_paths = bottom_return if isinstance(bottom_return, dict) else {'array_bottom_return': bottom_return}
        metal_paths = {**top_paths, **bottom_paths}
        self.add_qgeometry('poly', {'pad_left': left_pad, 'pad_right': right_pad})
        self.add_qgeometry('poly', {'pocket': pocket}, subtract=True)
        self.add_qgeometry('path', metal_paths, width=p.lead_width, fillet=p.array_bottom_fillet)
        self.add_qgeometry(
            'path',
            {f'{name}_cut': geom for name, geom in metal_paths.items()},
            width=p.lead_width + 2 * p.lead_gap,
            fillet=p.array_bottom_fillet,
            subtract=True,
        )
        if array_islands:
            self.add_qgeometry('poly', array_islands)
            self.add_qgeometry(
                'poly',
                {
                    f'{name}_cut': geom.buffer(p.trace_gap, join_style=1)
                    for name, geom in array_islands.items()
                },
                subtract=True,
            )
        self.add_qgeometry(
            'junction',
            {'main_junction': top_junction},
            width=p.trace_width,
            hfss_inductance=p.L_j,
            hfss_capacitance=p.C_j,
        )
        self.add_qgeometry(
            'junction',
            array_junctions,
            width=p.array_junction_width,
            hfss_inductance=junction_inductance,
        )
