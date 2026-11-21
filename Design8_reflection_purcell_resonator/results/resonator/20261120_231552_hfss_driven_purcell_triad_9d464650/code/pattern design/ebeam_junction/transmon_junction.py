"""Transmon / coupler SQUID e-beam mask: Dolan JJs, NCKU EBL standard geometry.

Two-layer fabrication drawing, not a lumped HFSS junction model. Rules,
layers and the EJ -> width model: skill/ebeam_junction.md. The geometry and
option names were adopted from the NCKU fabrication team's EBL standard
(2026-11); defaults are its Design2_FTF_Andersen coupler SQUID. Checks:
code/general/check_ebeam_junctions.py. The previous mirrored-pair version is
``archive/transmon_junction_layout3.py``.

Sketch of one junction, local frame of JJ1 (orientation 0, not to scale, um)::

            upper_contact_body_width (4; 7.6 for a rounded bay)
            |<------------>|              pin 'top_jj1'
            |   upper      |  6 um (upper_contact_body_length)
            |   contact    |
            |______________|  y = 1.5 (1 um below the bay mouth at y = 2.5)
                  ||  finger (jj_finger_width, 80 nm), right edge x = 0.338
                  ||                   ____
      bridge -> ||  ==========cross===|    |   y = 0   lower_contact_lead_width (2)
      (jj_bridge_gap)   electrode     |    |   cross electrode top edge y = -0.304,
                  ||  (jjN_cross_     |    |   width jjN_cross_electrode_width
                      electrode_width)|    |___
                                      |        |  lower_contact_body_width (4),
                                      |________|  y = -2.3 .. -7.5, centred at
                                        pin 'bottom_jj1'   x = lower_contact_offset

* Both junctions have the same orientation: JJ2 is JJ1 translated by
  ``squid_spacing`` along +x, never mirrored (NCKU rule; the asymmetric
  base-metal hooks accommodate one evaporation direction). The finger is on
  the left and the cross electrode on the right.
* ``jj1_cross_electrode_width`` / ``jj2_cross_electrode_width`` allow an
  asymmetric SQUID (Design6_TXT Q1: 384 / 157 nm). Resizing keeps the finger right
  edge and the electrode top edge fixed; the shadow openings follow the lower
  edge with ``shadow_clearance`` (75 nm).
* ``lower_contact_lead_offset`` moves the 2 um stem within the lower body
  (default 0 = centred, NCKU). Our Design3_FTF_NTHU / ``my_tmon_fab`` ground hooks sit
  only +2 um from the island hook, so use ``lower_contact_offset='2um'`` with
  ``lower_contact_lead_offset='1um'``; the stem must stay ``shadow_clearance``
  away from the shadow openings.
* ``squid=False`` draws only JJ1 (test keys). ``upper_contact_offset`` moves
  the upper body along x (NCKU T keys: 3 um); a 0.5 um ledge then joins the
  finger to the body.

Origin: SQUID centre, x midway between the two upper-body centres, y = 0 at
the top of the lower lead. Use ``code/qubit/jj_resistance_model.py``
(``squid_options(ej1, ej2)``) to turn EJ1 / EJ2 into the electrode widths.
Layers: 1 main dose, 2 shadow. All dimensions are parsed
in the design's units (normally mm).
"""
from math import isfinite

import numpy as np
from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent


class TransmonJunction(QComponent):
    """Dolan-bridge SQUID (or one JJ) for a transmon / tunable coupler, NCKU style.

    Pins, per drawn junction N = 1, 2: ``top_jjN`` at the upper body top
    (y = 1.5 + upper_contact_body_length, width upper_contact_body_width,
    facing +y) and ``bottom_jjN`` at the lower body bottom (y = -7.5, width
    lower_contact_body_width, facing -y).
    """

    FINGER_RIGHT = 0.338        # from the upper-body centre (NCKU)
    FINGER_TIP = -0.991
    UPPER_BOTTOM = 1.5          # upper body bottom = finger top
    BAR_TOP = -0.304            # cross electrode top edge (fixed)
    LOWER_STEM_BOTTOM, LOWER_END = -2.3, -7.5
    SHADOW_FLOOR, SHADOW_TOP = -1.241, 0.0
    LANDING_GAP, FINGER_SIDE_GAP, SHADOW_REACH, LANDING_WIDTH = 0.176, 0.074, 0.648, 0.594
    LEDGE = 0.5

    default_options = Dict(
        jj_finger_width='80nm', jj1_cross_electrode_width='271nm',
        jj2_cross_electrode_width='271nm', jj_bridge_gap='0.249um',
        shadow_clearance='0.075um', squid='True', squid_spacing='30um',
        upper_contact_body_width='4um', upper_contact_body_length='6um',
        upper_contact_offset='0um', lower_contact_lead_width='2um',
        lower_contact_body_width='4um', lower_contact_offset='4um',
        lower_contact_lead_offset='0um',
        layer_jj='1', layer_jj_shadow='2', jj_metal='True',
    )
    component_metadata = Dict(short_name='transmon_junction',
                              _qgeometry_table_poly='True')

    def make(self):
        p = self.p
        u = self.parse_value('1um')
        positive = ('jj_finger_width', 'jj1_cross_electrode_width',
                    'jj2_cross_electrode_width', 'jj_bridge_gap', 'shadow_clearance',
                    'squid_spacing', 'upper_contact_body_width', 'upper_contact_body_length',
                    'lower_contact_lead_width', 'lower_contact_body_width')
        v = {k: float(p[k]) / u for k in positive + ('upper_contact_offset',
                                                      'lower_contact_offset',
                                                      'lower_contact_lead_offset')}
        if not all(isfinite(x) for x in v.values()):
            raise ValueError('all dimensions must be finite')
        bad = [k for k in positive if not v[k] > 0]
        if bad:
            raise ValueError(f'{", ".join(bad)} must be positive')
        squid = _flag(p.squid, 'squid')
        if not _flag(p.jj_metal, 'jj_metal'):
            return

        x1 = -v['squid_spacing'] / 2 if squid else 0.0
        drawn = [('jj1', x1, v['jj1_cross_electrode_width'])]
        if squid:
            drawn.append(('jj2', x1 + v['squid_spacing'], v['jj2_cross_electrode_width']))
        metal, shadow, pins = {}, {}, []
        for tag, xc, cw in drawn:
            m, s, pn = self._junction(v, cw)
            for name, shape in m.items():
                metal[f'{name}_{tag}'] = draw.translate(shape, xc, 0)
            for name, shape in s.items():
                shadow[f'{name}_{tag}'] = draw.translate(shape, xc, 0)
            pins += [(f'{n}_{tag}', x + xc, y, w, d) for n, x, y, w, d in pn]
        if squid and any(a.intersection(b).area > 1e-9
                         for k, a in metal.items() if k.endswith('jj1')
                         for j, b in metal.items() if j.endswith('jj2')):
            raise ValueError('squid_spacing too small: the two junctions overlap')

        for shapes, layer in ((metal, p.layer_jj), (shadow, p.layer_jj_shadow)):
            placed = {}
            for name, shape in shapes.items():
                if not shape.is_valid or shape.is_empty or shape.area <= 0:
                    raise ValueError(f'Invalid mask geometry: {name}')
                shape = draw.scale(shape, u, u, origin=(0, 0))
                shape = draw.rotate(shape, p.orientation, origin=(0, 0))
                placed[name] = draw.translate(shape, p.pos_x, p.pos_y)
            self.add_qgeometry('poly', placed, layer=layer, chip=p.chip, subtract=False)

        for name, x, y, w, direction in pins:
            line = draw.LineString([(x, y - direction), (x, y)])
            line = draw.scale(line, u, u, origin=(0, 0))
            line = draw.rotate(line, p.orientation, origin=(0, 0))
            line = draw.translate(line, p.pos_x, p.pos_y)
            self.add_pin(name, np.asarray(line.coords), width=w*u,
                         input_as_norm=True, chip=p.chip)

    def _junction(self, v, cw):
        """One JJ in its local frame (origin: upper-body centre x, lower-lead top y)."""
        fw, bg, sc = v['jj_finger_width'], v['jj_bridge_gap'], v['shadow_clearance']
        bu, lu, ou = (v['upper_contact_body_width'], v['upper_contact_body_length'],
                      v['upper_contact_offset'])
        ls, bl, ol = (v['lower_contact_lead_width'], v['lower_contact_body_width'],
                      v['lower_contact_offset'])
        fr, ub = self.FINGER_RIGHT, self.UPPER_BOTTOM
        et, eb, el = self.BAR_TOP, self.BAR_TOP - cw, fr + bg
        sx = ol + v['lower_contact_lead_offset']      # stem centre
        stem_left = sx - ls / 2
        if eb < self.FINGER_TIP:
            raise ValueError('cross electrode extends below the finger tip')
        if abs(sx - ol) + ls / 2 > bl / 2 + 1e-9:
            raise ValueError('the lower stem must lie within lower_contact_body_width')
        if stem_left < fr + self.SHADOW_REACH + sc - 1e-9:
            raise ValueError('lower stem closer than shadow_clearance to the shadow openings; '
                             'increase lower_contact_offset or lower_contact_lead_offset')
        upper_left = ou - bu / 2
        if ou + bu / 2 < fr:
            raise ValueError('upper contact body must reach the finger right edge; '
                             'increase upper_contact_offset or its width')
        metal = {
            'upper_contact': _box(upper_left, ub, ou + bu / 2, ub + lu, 'upper contact'),
            'finger': _box(fr - fw, self.FINGER_TIP, fr, ub, 'finger'),
            'cross_electrode': _box(el, eb, stem_left, et, 'cross electrode'),
            'lower_contact': _merge(
                _box(stem_left, self.LOWER_STEM_BOTTOM, sx + ls / 2, 0, 'lower stem'),
                _box(ol - bl / 2, self.LOWER_END, ol + bl / 2, self.LOWER_STEM_BOTTOM,
                     'lower body')),
        }
        if upper_left > fr - fw:
            # Body moved off the finger (test-key bay): 0.5 um ledge joins them.
            ledge = _box(min(fr - 0.1, fr - fw), ub, upper_left, ub + self.LEDGE, 'ledge')
            metal['upper_contact'] = metal['upper_contact'].union(ledge)
        lr = min(fr - self.LANDING_GAP, fr - fw - sc)
        fs, reach = fr + self.FINGER_SIDE_GAP, fr + self.SHADOW_REACH
        lo, hi = self.SHADOW_FLOOR, self.SHADOW_TOP
        shadow = {
            'shadow_landing': _box(lr - self.LANDING_WIDTH, lo, lr, hi, 'landing opening'),
            'shadow_tip': _box(lr, lo, fs, self.FINGER_TIP - sc, 'tip opening'),
            'shadow_below': _box(fs, lo, reach, eb - sc, 'opening below electrode'),
            'shadow_bridge': _box(fs, eb - sc, el - sc, et + sc, 'bridge opening'),
            'shadow_above': _box(fs, et + sc, reach, hi, 'opening above electrode'),
        }
        pins = [('top', ou, ub + lu, bu, 1), ('bottom', ol, self.LOWER_END, bl, -1)]
        return metal, shadow, pins


def _flag(value, name):
    text = str(value).lower()
    if text not in ('true', 'false', '1', '0'):
        raise ValueError(f'{name} must be True or False')
    return text in ('true', '1')


def _box(x0, y0, x1, y1, what):
    if not (x1 - x0 > 1e-9 and y1 - y0 > 1e-9):
        raise ValueError(f'{what} has non-positive size for these options')
    return draw.box(x0, y0, x1, y1)


def _merge(a, b):
    # Stem and body share an edge; merge to one simple polygon without
    # collinear vertices so the drawing equals the NCKU 8-vertex outline.
    return draw.Polygon(a.union(b).simplify(0).exterior.coords)


transmon_junction = TransmonJunction
