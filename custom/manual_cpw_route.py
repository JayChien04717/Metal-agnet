# -*- coding: utf-8 -*-

"""Small absolute-waypoint CPW route for package-aware first-pass layouts."""

from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent


class ManualCPWRoute(QComponent):
    """Draw a CPW path through explicit waypoints and optional component pins."""

    component_metadata = Dict(
        short_name='ManualCPWRoute',
        _qgeometry_table_path='True',
    )

    default_options = Dict(
        start_component='',
        start_pin='',
        end_component='',
        end_pin='',
        points=[],
        lead=Dict(
            start_straight='0um',
            end_straight='0um',
        ),
        trace_width='10um',
        trace_gap='6um',
        fillet='99um',
    )

    def _pin(self, component_name, pin_name):
        component = self.design.components[component_name]
        return component.pins[pin_name]

    def _pin_middle(self, component_name, pin_name):
        return tuple(self._pin(component_name, pin_name)['middle'])

    def _pin_normal(self, component_name, pin_name):
        pin = self._pin(component_name, pin_name)
        normal = pin.get('normal', None)
        if normal is not None:
            return tuple(normal)

        points = pin.get('points', None)
        if points is None or len(points) < 2:
            return (0, 0)
        (x0, y0), (x1, y1) = points[:2]
        dx = y1 - y0
        dy = -(x1 - x0)
        length = (dx * dx + dy * dy) ** 0.5
        if length == 0:
            return (0, 0)
        return (dx / length, dy / length)

    def _lead_point(self, component_name, pin_name, distance):
        middle = self._pin_middle(component_name, pin_name)
        normal = self._pin_normal(component_name, pin_name)
        return (
            middle[0] + normal[0] * distance,
            middle[1] + normal[1] * distance,
        )

    def _parse_point(self, point):
        return (self.parse_value(point[0]), self.parse_value(point[1]))

    @staticmethod
    def _append_manhattan(points, point):
        if not points:
            points.append(point)
            return

        last_x, last_y = points[-1]
        next_x, next_y = point
        if abs(last_x - next_x) > 1e-15 and abs(last_y - next_y) > 1e-15:
            points.append((next_x, last_y))
        if points[-1] != point:
            points.append(point)

    def make(self):
        p = self.p
        points = []
        if p.start_component and p.start_pin:
            start_middle = self._pin_middle(p.start_component, p.start_pin)
            self._append_manhattan(points, start_middle)
            if p.lead.start_straight:
                self._append_manhattan(
                    points,
                    self._lead_point(
                        p.start_component, p.start_pin, p.lead.start_straight
                    ),
                )
        for point in p.points:
            self._append_manhattan(points, self._parse_point(point))
        if p.end_component and p.end_pin:
            if p.lead.end_straight:
                self._append_manhattan(
                    points,
                    self._lead_point(
                        p.end_component, p.end_pin, p.lead.end_straight
                    ),
                )
            self._append_manhattan(points, self._pin_middle(p.end_component, p.end_pin))

        route = draw.LineString(points)
        self.add_qgeometry('path', {'trace': route}, width=p.trace_width, fillet=p.fillet)
        self.add_qgeometry(
            'path',
            {'cut': route},
            width=p.trace_width + 2 * p.trace_gap,
            fillet=p.fillet,
            subtract=True,
        )
