"""Parameterized transmon (SQUID) junction-mask geometry from ``layout3 ftf.gds``.

This is a two-layer fabrication drawing, not a lumped HFSS junction model.
It follows the same construction and option names as FluxoniumJunction:
a pad with a small vertical finger (``finger_width``), a horizontal junction bar
(``junction_width``) that joins a lower lead, and recesses held at the shared
``junction_gap``. See transmon_junction.md and transmon_junction_parameters.png.
All dimensions are parsed in the design's units (normally mm in Qiskit Metal).
"""
from math import isfinite

import numpy as np
from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent


class TransmonJunction(QComponent):
    """Dolan-bridge junction, or a mirrored SQUID pair, for a transmon.

    Local origin is the SQUID centre of the source (0, -1578) um. The left
    junction is built from parameters; with ``squid=True`` (source default)
    it is mirrored about x=0 to form the second junction.

    ``finger_width``: x width of the vertical finger (source 0.1um). Its right
    edge stays flush with the top pad edge, so it grows away from the bar.
    ``junction_width``: y width of the horizontal junction bar (source
    0.25um). Its bottom edge stays on the top of the lower lead so the bar
    and lead remain joined; it grows upwards.
    ``junction_gap``: shared recess clearance at the pad landing, both sides
    and the tip of the finger, and around the bridge between finger and bar.
    ``contact_extension``: lengthens the top pad upwards (source end y=3.0um)
    and the lower lead downwards (source end y=-3.183um); pins move with them.
    """

    # Source left-junction coordinates (um), relative to the SQUID centre.
    PAD_EDGE_X = -4.325        # finger right edge = top pad right edge
    PAD_BOTTOM_Y = 0.5         # finger top = top pad bottom
    FINGER_TIP_Y = -0.491      # finger bottom
    BAR_X = (-4.076, -2.308)   # junction bar x extent (right end over lead)
    BAR_BOTTOM_Y = -0.054      # bar bottom = lower lead top

    default_options = Dict(
        junction_gap='0.075um', junction_width='0.25um', finger_width='0.1um',
        contact_extension='0um', squid='True', layer='1', layer_recess='2', jj_metal='True',
    )
    component_metadata = Dict(short_name='transmon_junction',
                              _qgeometry_table_poly='True')

    @staticmethod
    def reference():
        """Return the immutable, device-local source polygon snapshot (um)."""
        return {'source_file': 'layout3 ftf.gds',
         'source_sha256': '6546a6202ce78dc716b100501d88884f8ed6f5854e17b0e59d962bce10ea2b26',
         'source_cell': 'TOP',
         'unit': 'um',
         'selection_window_um': [[-40, -1620], [40, -1540]],
         'translation_um': [0, 1578],
         'note': 'Transmon SQUID junction pair; only layers 1 and 2 have polygons in this window. Created '
                 '2026-11-18.',
         'polygons_um': {'1': [[[-5.008, 0.5], [-5.0, 3.0], [-4.325, 3.0], [-4.325, 0.5]],
                               [[-2.809, -3.183], [-2.809, -0.054], [-2.308, -0.054], [-2.309, -3.183]],
                               [[4.325, 0.5], [4.325, 3.0], [5.0, 3.0], [5.008, 0.5]],
                               [[2.309, -3.183], [2.308, -0.054], [2.809, -0.054], [2.809, -3.183]],
                               [[-4.425, -0.491], [-4.425, 0.5], [-4.325, 0.5], [-4.325, -0.491]],
                               [[-4.076, -0.054], [-4.076, 0.196], [-2.308, 0.196], [-2.308, -0.054]],
                               [[4.325, -0.491], [4.325, 0.5], [4.425, 0.5], [4.425, -0.491]],
                               [[2.308, -0.054], [2.308, 0.196], [4.076, 0.196], [4.076, -0.054]]],
                         '2': [[[-5.095, -0.741], [-5.095, 0.45], [-4.501, 0.45], [-4.501, -0.741]],
                               [[-4.501, -0.741], [-4.501, -0.566], [-4.251, -0.566], [-4.251, -0.741]],
                               [[-4.251, -0.741], [-4.251, -0.129], [-3.677, -0.129], [-3.677, -0.741]],
                               [[-4.251, -0.129], [-4.251, 0.271], [-4.151, 0.271], [-4.151, -0.129]],
                               [[-4.251, 0.271], [-4.251, 1.312], [-3.677, 1.312], [-3.677, 0.271]],
                               [[4.251, -0.741], [4.251, -0.566], [4.501, -0.566], [4.501, -0.741]],
                               [[3.677, -0.741], [3.677, -0.129], [4.251, -0.129], [4.251, -0.741]],
                               [[4.151, -0.129], [4.151, 0.271], [4.251, 0.271], [4.251, -0.129]],
                               [[3.677, 0.271], [3.677, 1.312], [4.251, 1.312], [4.251, 0.271]],
                               [[4.501, -0.741], [4.501, 0.45], [5.095, 0.45], [5.095, -0.741]]]}}

    @staticmethod
    def _flag(value, name):
        text = str(value).lower()
        if text not in ('true', 'false', '1', '0'):
            raise ValueError(f'{name} must be True or False')
        return text in ('true', '1')

    def make(self):
        p = self.p
        u = self.parse_value('1um')
        gap, jj_width, sj = (float(v) / u for v in
                             (p.junction_gap, p.junction_width, p.finger_width))
        if not all(isfinite(v) and v > 0 for v in (gap, jj_width, sj)):
            raise ValueError('junction_gap, junction_width and finger_width must be positive')
        ext = float(p.contact_extension) / u
        if not isfinite(ext) or ext < 0:
            raise ValueError('contact_extension must be >= 0')
        squid = self._flag(p.squid, 'squid')
        if not self._flag(p.jj_metal, 'jj_metal'):
            return

        def box(x0, y0, x1, y1, what):
            # draw.box silently reorders inverted corners, so check explicitly.
            if not (x1 - x0 > 1e-9 and y1 - y0 > 1e-9):
                raise ValueError(f'{what} has non-positive size; reduce '
                                 'junction_gap, junction_width or finger_width')
            return draw.box(x0, y0, x1, y1)

        source = self.reference()['polygons_um']
        fx1 = self.PAD_EDGE_X
        fx0 = fx1 - sj
        bx0, bx1 = self.BAR_X
        y0 = self.BAR_BOTTOM_Y
        y1 = y0 + jj_width
        g = gap
        # Left junction; source layer-1 indices 0 (pad) and 1 (lead) are kept.
        # Outer ends move by contact_extension (pad top y=3.0, lead end y=-3.183).
        pad = [(x, y + ext if y >= 3.0 - 1e-9 else y) for x, y in source['1'][0]]
        lead = [(x, y - ext if y <= -3.183 + 1e-9 else y) for x, y in source['1'][1]]
        metal = {
            'pad_L': draw.Polygon(pad),
            'lead_L': draw.Polygon(lead),
            'finger_L': box(fx0, self.FINGER_TIP_Y, fx1, self.PAD_BOTTOM_Y, 'finger'),
            'bar_L': box(bx0, y0, bx1, y1, 'junction bar'),
        }
        tip = self.FINGER_TIP_Y - g
        recess = {
            'landing_L': box(-5.095, -0.741, fx0-g, self.PAD_BOTTOM_Y-g, 'recess landing'),
            'tip_L': box(fx0-g, -0.741, fx1+g, tip, 'recess tip'),
            'below_bar_L': box(fx1+g, -0.741, -3.677, y0-g, 'recess below bar'),
            'bridge_L': box(fx1+g, y0-g, bx0-g, y1+g, 'recess bridge'),
            'above_bar_L': box(fx1+g, y1+g, -3.677, 1.312, 'recess above bar'),
        }
        if squid:
            for group in (metal, recess):
                for name in [k for k in group if k.endswith('_L')]:
                    group[name[:-1] + 'R'] = draw.scale(group[name], -1, 1, origin=(0, 0))

        for shapes, layer in ((metal, p.layer), (recess, p.layer_recess)):
            placed = {}
            for name, shape in shapes.items():
                if not shape.is_valid or shape.is_empty or shape.area <= 0:
                    raise ValueError(f'Invalid mask geometry: {name}')
                shape = draw.scale(shape, u, u, origin=(0, 0))
                shape = draw.rotate(shape, p.orientation, origin=(0, 0))
                placed[name] = draw.translate(shape, p.pos_x, p.pos_y)
            self.add_qgeometry('poly', placed, layer=layer, chip=p.chip, subtract=False)

        # Pins at the outer ends of the pads (top) and leads (bottom).
        pads = [('L', 1)] + ([('R', -1)] if squid else [])
        for side, sx in pads:
            for name, x, y, w, direction in (
                    ('top', -4.6665, 3.0 + ext, 0.683, 1),
                    ('bottom', -2.5585, -3.183 - ext, 0.501, -1)):
                line = draw.LineString([(sx*x, y-direction), (sx*x, y)])
                line = draw.scale(line, u, u, origin=(0, 0))
                line = draw.rotate(line, p.orientation, origin=(0, 0))
                line = draw.translate(line, p.pos_x, p.pos_y)
                self.add_pin(f'{name}_{side}', np.asarray(line.coords), width=w*u,
                             input_as_norm=True, chip=p.chip)


transmon_junction = TransmonJunction
