"""Inward jogs on a Design5 CPW tap, with the original outer line axis."""

from math import atan2, sqrt

import numpy as np
from qiskit_metal import Dict, draw
from shapely.ops import unary_union
from tap_junction import TapJunction


def transition_span(inset, radius):
    """Longitudinal span of two opposite circular arcs joining parallel lines."""
    if not 0 < inset < 2 * radius:
        raise ValueError('Tap inset requires two arcs below 90 degrees')
    return sqrt(4 * radius * inset - inset**2)


class InsetTap(TapJunction):
    """CPW T with a straight inset centre and tangent circular transitions."""

    default_options = Dict(inset='105um', bend_radius='99um', end_straight='5um')

    def make(self):
        p = self.p
        half = p.prime_length / 2
        span = transition_span(p.inset, p.bend_radius)
        if span + p.end_straight >= half:
            raise ValueError('Tap must have a straight centre segment')
        sign = 1 if str(p.stub_up).lower() == 'true' else -1
        theta = np.linspace(0, 2 * atan2(p.inset, span), 65)
        advance = p.bend_radius * np.sin(theta)
        shift = p.bend_radius * (1 - np.cos(theta))
        first = np.column_stack((advance, shift))
        second = (span, p.inset) - first[-2::-1]
        left = np.concatenate((first, second))
        left[:, 0] += p.end_straight - half
        left[:, 1] *= sign
        left = np.concatenate(([[-half, 0]], left))
        right = left[::-1] * (-1, 1)
        prime = draw.LineString(np.concatenate((left, right)))
        junction = sign * p.inset
        tip = junction + sign * p.stub_length
        stub = draw.LineString([(0, junction), (0, tip)])

        def place(geometry):
            rotated = draw.rotate(geometry, p.orientation, origin=(0, 0))
            return draw.translate(rotated, p.pos_x, p.pos_y)

        for name, distance, subtract in (
            ('tap_trace', p.trace_width / 2, False),
            ('tap_cut', p.trace_width / 2 + p.trace_gap, True),
        ):
            outline = unary_union([path.buffer(distance, cap_style=2) for path in (prime, stub)])
            self.add_qgeometry('poly', {name: place(outline)}, layer=p.layer, subtract=subtract)

        width = p.trace_width
        for name, points in (
            ('prime_start', [(-half + width, 0), (-half, 0)]),
            ('prime_end', [(half - width, 0), (half, 0)]),
            ('second_end', [(0, tip - sign * width), (0, tip)]),
        ):
            pin = np.array(place(draw.LineString(points)).coords)
            self.add_pin(name, points=pin, width=width, gap=p.trace_gap, input_as_norm=True)
