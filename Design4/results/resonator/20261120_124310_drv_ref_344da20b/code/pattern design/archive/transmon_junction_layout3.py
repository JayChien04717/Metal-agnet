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

    ``hook_contact=True`` (development option) replaces the straight pad and
    lead ends by tapered contacts that fill the two Dolan hooks of the base
    metal. Each contact keeps the source pad/lead as a narrow neck of
    ``neck_length`` next to the junction, then tapers to ``hook_width + 2 *
    hook_margin`` at ``taper_inset`` inside the hook mouth and runs to
    ``hook_overlap`` past the hook bottom. The hook mouths are ``hook_gap``
    apart in y, centred on the source pad/lead midpoint (y=-0.0915um). The top
    hook is centred on the pad (x=-4.6665um), the bottom hook ``hook_offset``
    further along +x. ``tip_width`` narrows each neck when it is below the
    source width (pad 0.683um, lead 0.501um), keeping the edge that meets the
    finger/bar fixed. ``layer_contact`` ('' = ``layer``) puts the tapered part
    on its own layer. ``contact_extension`` is ignored in this mode.
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
        hook_contact='False', hook_gap='5um', hook_width='2um', hook_depth='4um',
        hook_offset='2um', hook_margin='1um', hook_overlap='1um', taper_inset='1um',
        neck_length='0.9um', tip_width='0.5um', layer_contact='',
    )
    # Source pad/lead midpoint in y; the hook mouths are centred on it.
    HOOK_CENTRE_Y = (3.0 + -3.183) / 2
    PAD_CENTRE_X = -4.6665
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
        hook = self._flag(p.hook_contact, 'hook_contact')
        contacts = {}
        pin_ends = None
        if hook:
            hk = {k: float(p[k]) / u for k in (
                'hook_gap', 'hook_width', 'hook_depth', 'hook_offset', 'hook_margin',
                'hook_overlap', 'taper_inset', 'neck_length', 'tip_width')}
            tip = hk['tip_width']
            if not tip > 0:
                raise ValueError('tip_width must be > 0')
            if not all(isfinite(v) for v in hk.values()):
                raise ValueError('hook_* options must be finite')
            if min(hk['hook_gap'], hk['hook_width'], hk['hook_depth'], hk['neck_length']) <= 0 \
                    or min(hk['hook_margin'], hk['hook_overlap'], hk['taper_inset']) < 0:
                raise ValueError('hook_gap/width/depth and neck_length must be > 0; '
                                 'hook_margin, hook_overlap and taper_inset >= 0')
            w = hk['hook_width'] + 2 * hk['hook_margin']
            inset, far = hk['taper_inset'], hk['hook_depth'] + hk['hook_overlap']
            pad_poly, lead_poly = draw.Polygon(source['1'][0]), draw.Polygon(source['1'][1])
            # Top: neck = source pad up to y_neck, then taper into the top hook.
            mouth = self.HOOK_CENTRE_Y + hk['hook_gap'] / 2
            y_n = self.PAD_BOTTOM_Y + hk['neck_length']
            hx = self.PAD_CENTRE_X
            if not (y_n < mouth + inset and y_n < 3.0):
                raise ValueError('neck_length too long for hook_gap / taper_inset')
            neck = pad_poly.intersection(box(-6, self.PAD_BOTTOM_Y - 1, -4, y_n, 'top neck'))
            xl, _, xr, _ = neck.intersection(draw.LineString([(-6, y_n), (-4, y_n)])).bounds
            if tip < xr - xl - 0.005:   # narrow the neck, finger-side edge fixed
                xl = xr - tip
                neck = box(xl, self.PAD_BOTTOM_Y, xr, y_n, 'top neck')
            top =draw.Polygon([(xl, y_n), (xr, y_n), (hx + w/2, mouth + inset),
                                (hx + w/2, mouth + far), (hx - w/2, mouth + far),
                                (hx - w/2, mouth + inset)])
            # Bottom: neck = source lead down to y_neck, then taper into the
            # bottom hook, which sits hook_offset further along +x.
            mouth_b = self.HOOK_CENTRE_Y - hk['hook_gap'] / 2
            y_nb = self.BAR_BOTTOM_Y - hk['neck_length']
            hxb = hx + hk['hook_offset']
            if not (y_nb > mouth_b - inset and y_nb > -3.183):
                raise ValueError('neck_length too long for hook_gap / taper_inset')
            neck_b = lead_poly.intersection(box(-3.5, y_nb, -2, self.BAR_BOTTOM_Y + 1, 'bottom neck'))
            bl, _, br, _ = neck_b.intersection(draw.LineString([(-3.5, y_nb), (-2, y_nb)])).bounds
            if tip < br - bl - 0.005:   # narrow the neck, bar-side edge fixed
                bl, br = self.BAR_X[1] - tip, self.BAR_X[1]
                neck_b = box(bl, y_nb, br, self.BAR_BOTTOM_Y, 'bottom neck')
            bottom = draw.Polygon([(br, y_nb), (bl, y_nb), (hxb - w/2, mouth_b - inset),
                                   (hxb - w/2, mouth_b - far), (hxb + w/2, mouth_b - far),
                                   (hxb + w/2, mouth_b - inset)])
            for name, shape in (('top', top), ('bottom', bottom)):
                if not shape.is_valid:
                    raise ValueError(f'{name} hook contact is self-intersecting; '
                                     'widen hook_width/hook_margin or move the hook')
                for rname, r in recess.items():
                    if shape.intersection(r).area > 1e-9:
                        raise ValueError(f'{name} hook contact overlaps recess {rname}; '
                                         'increase neck_length')
            if str(p.layer_contact).strip() in ('', str(p.layer)):
                metal['pad_L'] = neck.union(top)
                metal['lead_L'] = neck_b.union(bottom)
            else:
                metal['pad_L'], metal['lead_L'] = neck, neck_b
                contacts['hook_top_L'], contacts['hook_bottom_L'] = top, bottom
            pin_ends = ((hx, mouth + far, w), (hxb, mouth_b - far, w))

        if squid:
            for group in (metal, recess, contacts):
                for name in [k for k in group if k.endswith('_L')]:
                    group[name[:-1] + 'R'] = draw.scale(group[name], -1, 1, origin=(0, 0))

        for shapes, layer in ((metal, p.layer), (recess, p.layer_recess),
                              (contacts, p.layer_contact or p.layer)):
            if not shapes:
                continue
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
        ends = (('top', -4.6665, 3.0 + ext, 0.683, 1),
                ('bottom', -2.5585, -3.183 - ext, 0.501, -1))
        if pin_ends:
            ends = tuple((n, x, y, w, d) for (n, _, _, _, d), (x, y, w) in zip(ends, pin_ends))
        for side, sx in pads:
            for name, x, y, w, direction in ends:
                line = draw.LineString([(sx*x, y-direction), (sx*x, y)])
                line = draw.scale(line, u, u, origin=(0, 0))
                line = draw.rotate(line, p.orientation, origin=(0, 0))
                line = draw.translate(line, p.pos_x, p.pos_y)
                self.add_pin(f'{name}_{side}', np.asarray(line.coords), width=w*u,
                             input_as_norm=True, chip=p.chip)


transmon_junction = TransmonJunction
