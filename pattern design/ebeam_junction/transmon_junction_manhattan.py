"""Transmon / coupler e-beam mask: Manhattan-style JJs (single JJ or SQUID).

Two-layer fabrication drawing (main JJ dose + undercut / bandage layer), not a
lumped HFSS junction model. Geometry matches the 45-degree-pad Manhattan JJ
standard (4 x 8 um rotated landing pads, orthogonal vertical/horizontal fingers
crossing at (+2, -2) um with 1 um stubs, 5 x 9 um pad undercut/bandage boxes,
and 1 um tip undercut boxes).

Sketch of one Manhattan junction, local frame of JJ1 (orientation 0, um)::

                        /\\   upper-right pad (4 x 8 um, +45 deg)
                       /  \\  with 5 x 9 um undercut/bandage box
                      /    \\
                      \\    /
                       \\  /
                        ||   vertical finger (jjN_finger_y_width, x = +2)
                        ||
          /\\            ||
         /  \\===========++====[ ]  horizontal finger (jjN_finger_x_width, y = -2)
        /    \\          ||         + 1 um stub + tip undercut box
        \\    /          ||
         \\  /          [  ]  1 um vertical stub + tip undercut box
          \\/
       lower-left pad (4 x 8 um, +45 deg)

* **Orientation rule.** Every Manhattan JJ has its fingers strictly along x and y
  so one two-angle shadow evaporation covers the entire chip. In a SQUID
  (``squid='True'``), JJ2 is JJ1 translated by ``squid_spacing`` along +x —
  never mirrored.
* Per-branch finger widths (``jj1_finger_x_width``, ``jj1_finger_y_width``,
  ``jj2_finger_x_width``, ``jj2_finger_y_width``) support both symmetric and
  asymmetric SQUIDs; the crossing centre ``(x_cross, y_cross)`` stays fixed at
  ``(+2, -2) um`` as finger widths change.
* Layers: ``layer_jj`` (default ``'1'``; main JJ dose, single merged 20-vertex
  polygon per JJ) and ``layer_jj_undercut`` (default ``'2'``; pad undercut /
  bandage boxes and finger-tip undercut boxes). On a chip with base metal on
  layer 1, set e.g. ``layer_jj='11'``, ``layer_jj_undercut='12'`` (or
  ``'13'`` / ``'18'``).

Origin: JJ centre (or midpoint between JJ1 and JJ2 when ``squid='True'``).
All dimensions are parsed in the design's units (normally mm).
"""
from math import isfinite, sqrt

import numpy as np
from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent


class TransmonJunctionManhattan(QComponent):
    """Manhattan-style e-beam JJ (single JJ or 2-JJ SQUID) for transmons and couplers.

    Pins, per drawn junction N = 1, 2: ``top_jjN`` at the upper-right pad centre
    (facing +45 deg outward along +u) and ``bottom_jjN`` at the lower-left pad
    centre (facing -135 deg outward along -u).
    """

    SQRT2_INV = 1.0 / sqrt(2.0)

    default_options = Dict(
        jj1_finger_x_width='190nm',
        jj1_finger_y_width='170nm',
        jj2_finger_x_width='190nm',
        jj2_finger_y_width='170nm',
        finger_stub_length='1.0um',
        cross_x='2.0um',
        cross_y='-2.0um',
        pad_length='8.0um',
        pad_width='4.0um',
        pad_inner_gap='4.456um',
        pad_undercut_margin='0.5um',
        tip_undercut_length='1.0um',
        tip_undercut_overlap='0.1um',
        tip_undercut_margin='0.125um',
        squid='True',
        squid_spacing='30um',
        layer_jj='1',
        layer_jj_undercut='2',
        jj_metal='True',
    )
    component_metadata = Dict(
        short_name='transmon_junction_manhattan',
        _qgeometry_table_poly='True',
    )

    def make(self):
        p = self.p
        u = self.parse_value('1um')
        positive = (
            'jj1_finger_x_width', 'jj1_finger_y_width',
            'jj2_finger_x_width', 'jj2_finger_y_width',
            'finger_stub_length', 'pad_length', 'pad_width',
            'pad_inner_gap', 'squid_spacing',
            'tip_undercut_length',
        )
        non_negative = ('pad_undercut_margin', 'tip_undercut_overlap', 'tip_undercut_margin')
        signed = ('cross_x', 'cross_y')
        v = {k: float(p[k]) / u for k in positive + non_negative + signed}
        if not all(isfinite(x) for x in v.values()):
            raise ValueError('all dimensions must be finite')
        bad = [k for k in positive if not v[k] > 0]
        if bad:
            raise ValueError(f'{", ".join(bad)} must be positive')
        bad_nn = [k for k in non_negative if v[k] < 0]
        if bad_nn:
            raise ValueError(f'{", ".join(bad_nn)} must be non-negative')
        if v['tip_undercut_overlap'] >= v['tip_undercut_length']:
            raise ValueError('tip_undercut_overlap must be smaller than tip_undercut_length')

        squid = _flag(p.squid, 'squid')
        if not _flag(p.jj_metal, 'jj_metal'):
            return

        x1 = -v['squid_spacing'] / 2.0 if squid else 0.0
        drawn = [('jj1', x1, v['jj1_finger_x_width'], v['jj1_finger_y_width'])]
        if squid:
            drawn.append(('jj2', x1 + v['squid_spacing'], v['jj2_finger_x_width'], v['jj2_finger_y_width']))

        metal, undercut, pins = {}, {}, []
        for tag, xc, wx, wy in drawn:
            m, uc, pn = self._junction(v, wx, wy)
            for name, shape in m.items():
                metal[f'{name}_{tag}'] = draw.translate(shape, xc, 0.0)
            for name, shape in uc.items():
                undercut[f'{name}_{tag}'] = draw.translate(shape, xc, 0.0)
            pins += [(f'{n}_{tag}', x + xc, y, dx, dy, w) for n, x, y, dx, dy, w in pn]

        if squid:
            for k, a in metal.items():
                if not k.endswith('jj1'):
                    continue
                for j, b in metal.items():
                    if j.endswith('jj2') and a.intersection(b).area > 1e-9:
                        raise ValueError('squid_spacing too small: the two Manhattan junctions overlap')
            for k, a in undercut.items():
                if not k.endswith('jj1'):
                    continue
                for j, b in undercut.items():
                    if j.endswith('jj2') and a.intersection(b).area > 1e-9:
                        raise ValueError('squid_spacing too small: the two Manhattan undercut regions overlap')

        for shapes, layer in ((metal, p.layer_jj), (undercut, p.layer_jj_undercut)):
            placed = {}
            for name, shape in shapes.items():
                if not shape.is_valid or shape.is_empty or shape.area <= 0:
                    raise ValueError(f'Invalid mask geometry: {name}')
                shape = draw.scale(shape, u, u, origin=(0, 0))
                shape = draw.rotate(shape, p.orientation, origin=(0, 0))
                placed[name] = draw.translate(shape, p.pos_x, p.pos_y)
            self.add_qgeometry('poly', placed, layer=layer, chip=p.chip, subtract=False)

        for name, x, y, dx, dy, w in pins:
            line = draw.LineString([(x - dx, y - dy), (x, y)])
            line = draw.scale(line, u, u, origin=(0, 0))
            line = draw.rotate(line, p.orientation, origin=(0, 0))
            line = draw.translate(line, p.pos_x, p.pos_y)
            self.add_pin(name, np.asarray(line.coords), width=w * u, input_as_norm=True, chip=p.chip)

    def _junction(self, v, wx, wy):
        """One Manhattan JJ in its local frame (um), centred at (0, 0)."""
        pl, pw, gap = v['pad_length'], v['pad_width'], v['pad_inner_gap']
        cx, cy = v['cross_x'], v['cross_y']
        stub = v['finger_stub_length']
        u_in = gap / 2.0
        u_out = u_in + pl
        u_mid = (u_in + u_out) / 2.0

        # Upper-right (+u) and lower-left (-u) 45-degree landing pads
        pad_top = draw.rotate(draw.box(u_in, -pw / 2.0, u_out, pw / 2.0), 45.0, origin=(0, 0))
        pad_bot = draw.rotate(draw.box(-u_out, -pw / 2.0, -u_in, pw / 2.0), 45.0, origin=(0, 0))

        # Check that each finger starts cleanly inside its pad
        # Upper pad inner edge is line x + y = u_in * sqrt(2); vertical finger at x = cx meets it at y = u_in*sqrt(2) - cx
        y_top_entry = u_in * sqrt(2.0) - cx
        x_bot_entry = -u_in * sqrt(2.0) - cy
        y_top_anchor = y_top_entry + min(1.0, pl / 4.0)
        x_bot_anchor = x_bot_entry - min(1.0, pl / 4.0)

        y_v_tip = cy - wx / 2.0 - stub
        x_h_tip = cx + wy / 2.0 + stub

        if y_top_entry <= cy + wx / 2.0:
            raise ValueError('vertical finger crossing is too high to clear the upper pad')
        if x_bot_entry >= cx - wy / 2.0:
            raise ValueError('horizontal finger crossing is too far left to clear the lower pad')
        if not pad_top.contains(draw.Point(cx, y_top_anchor)):
            raise ValueError('vertical finger does not intersect the upper-right pad; check cross_x / pad_width')
        if not pad_bot.contains(draw.Point(x_bot_anchor, cy)):
            raise ValueError('horizontal finger does not intersect the lower-left pad; check cross_y / pad_width')

        finger_v = draw.box(cx - wy / 2.0, y_v_tip, cx + wy / 2.0, y_top_anchor)
        finger_h = draw.box(x_bot_anchor, cy - wx / 2.0, x_h_tip, cy + wx / 2.0)

        jj_poly = pad_top.union(pad_bot).union(finger_v).union(finger_h).simplify(0)
        if jj_poly.geom_type != 'Polygon':
            raise ValueError('Manhattan JJ main dose failed to merge into a single polygon')
        metal = {'jj_main': draw.Polygon(jj_poly.exterior.coords)}

        # Undercut / bandage layer
        m = v['pad_undercut_margin']
        uc_top = draw.rotate(draw.box(u_in - m, -(pw / 2.0 + m), u_out + m, pw / 2.0 + m), 45.0, origin=(0, 0))
        uc_bot = draw.rotate(draw.box(-(u_out + m), -(pw / 2.0 + m), -(u_in - m), pw / 2.0 + m), 45.0, origin=(0, 0))

        tl, tov, tm = v['tip_undercut_length'], v['tip_undercut_overlap'], v['tip_undercut_margin']
        uc_tip_h = draw.box(x_h_tip - tov, cy - wx / 2.0 - tm, x_h_tip - tov + tl, cy + wx / 2.0 + tm)
        uc_tip_v = draw.box(cx - wy / 2.0 - tm, y_v_tip + tov - tl, cx + wy / 2.0 + tm, y_v_tip + tov)

        undercut = {
            'undercut_pad_top': uc_top,
            'undercut_pad_bot': uc_bot,
            'undercut_tip_h': uc_tip_h,
            'undercut_tip_v': uc_tip_v,
        }

        c_xy = u_mid * self.SQRT2_INV
        s = self.SQRT2_INV
        pins = [
            ('top', c_xy, c_xy, s, s, pw),
            ('bottom', -c_xy, -c_xy, -s, -s, pw),
        ]
        return metal, undercut, pins


def _flag(value, name):
    text = str(value).lower()
    if text not in ('true', 'false', '1', '0'):
        raise ValueError(f'{name} must be True or False')
    return text in ('true', '1')


transmon_junction_manhattan = TransmonJunctionManhattan
