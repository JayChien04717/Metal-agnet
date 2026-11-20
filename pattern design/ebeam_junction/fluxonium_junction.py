"""Fluxonium e-beam mask: small Dolan JJ + JJ array, NCKU EBL standard geometry.

Multilayer fabrication drawing, not a lumped HFSS junction model. Rules,
layers and the EJ -> width model: skill/ebeam_junction.md. The geometry and
option names were adopted from the NCKU fabrication team's EBL standard
(2026-11); defaults are its Design2_FTF_Andersen fluxonium. Checks:
code/general/check_ebeam_junctions.py. The previous layout3-based version is
``archive/fluxonium_junction_layout3.py``.

Sketch, ``side='left'`` (orientation 0, not to scale, um)::

                      array_contact_body_width (4)
                         |<------>|            pin 'top' (0, +array_contact_end)
                         |        |
                          \      /   taper  y = array_contact_taper_end .. _start
                           |    |    funnel mouth = array_contact_lead_width (2)
   return link            _|____|____  finger (jj_finger_width, 80 nm) + bridge
   (array_return_link_    |_||_||_|..==|__|- -===  cross electrode (jj_cross_
    width, 2)             ......        |  |       electrode_width) -> return
   |  num_array per arm, array_element_length x array_row_width,   |  |  lead
   |  array_bridge_gap apart                                       |  | (jj_return_
   |_ |_||_||_|.............................|_||_|  ____ ________ |__|  lead_width)
                                           |    |
                          bottom contact    \  /     pin 'bottom' (0, -array_contact_end)

* Every small JJ has the vertical finger on the LEFT and the horizontal cross
  electrode on the RIGHT, for one shared evaporation direction (NCKU rule).
  ``side='right'`` is the NCKU ``FLUX_RIGHT`` drawing: the array and contacts
  are rotated by 180 deg (array towards +x, JJ on the lower arm), and the JJ
  is re-arranged (finger hangs from the return lead, cross electrode joins the
  mouth bar) so that the finger stays on the left. It is not a mirror image.
* Resizing keeps the finger right edge and the cross electrode top edge
  fixed; the electrode's lower edge sets ``jj_cross_electrode_width``. The
  local shadow openings follow that edge with ``shadow_clearance`` (75 nm).
* The funnel mouth keeps its array-facing edge (x = -0.501 um for ``left``)
  and the whole small JJ moves with the mouth's other edge.
* The contact body is centred on x = 0 and should cover the base-metal bay
  (NCKU: 4 um for a 2 um rectangular hook, 7.6 um for the 5.6 um rounded bay).

Use ``code/qubit/jj_resistance_model.py`` to turn EJ into
``jj_cross_electrode_width`` (``size_from_ej(4.0).options()``).
Layers: 1 small-JJ main dose, 2 small-JJ/contact shadow, 3 array main dose
(array, return, contacts), 4 array shadow; base metal is not drawn here.
All dimensions are parsed in the design's units (normally mm).
"""
from math import isfinite

import numpy as np
from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import QComponent
from shapely.affinity import rotate as _rotate


class FluxoniumJunction(QComponent):
    """Fluxonium small JJ + array e-beam mask (NCKU EBL standard).

    Origin: contact centre (array-row centre lines at y = +/-4 um). Pins
    ``top`` / ``bottom`` at (0, +/-``array_contact_end``), width
    ``array_contact_body_width``, facing outwards.

    Array: ``num_array`` cells per arm (cells, not an established number of
    electrical junctions), ``array_element_length`` along x,
    ``array_row_width`` along y, ``array_bridge_gap`` between cells and from
    the last cell to the return block. ``array_return_link_width`` is the x
    width of the vertical link joining the two arms (its array-side edge is
    fixed). Contacts: funnel mouth ``array_contact_lead_width`` up to
    ``array_contact_taper_start``, taper to ``array_contact_body_width`` at
    ``array_contact_taper_end``, body to ``array_contact_end``; shadow strips
    ``array_contact_shadow_width`` wide beside the body (0.504 um neck strips).
    Small JJ: ``jj_finger_width`` (x of the vertical finger),
    ``jj_cross_electrode_width`` (y of the horizontal electrode),
    ``jj_bridge_gap`` (finger to electrode end), ``jj_return_lead_width``.
    """

    ARM_Y = 4.0                 # array-row centre line (um)
    MOUTH_EDGE = -0.501         # array-facing funnel edge = first cell edge
    # Small JJ, side='left', relative to the mouth's JJ-side edge xm.
    FINGER_RIGHT = 2.983        # finger right edge = top bar end
    RETURN_LEFT = 5.229         # return lead left edge = cross electrode end
    BAR_TOP = 2.696             # cross electrode top edge (fixed)
    FINGER_TIP = 2.009
    # Small JJ, side='right' (after the 180 deg turn), relative to xm.
    R_ELECTRODE_RIGHT = -2.885  # cross electrode joins the mouth bar here
    R_ELECTRODE_LENGTH = 1.998
    R_BAR_TOP = -3.0
    R_FINGER = (-3.687, -2.142)  # finger tip, return-lead bottom
    # Finger-side shadow offsets from the finger right edge (NCKU, fixed).
    LANDING_GAP, FINGER_SIDE_GAP, SHADOW_REACH = 0.176, 0.074, 0.647
    NECK_STRIP = 0.504

    default_options = Dict(
        num_array='43', array_element_length='0.998um', array_row_width='2um',
        array_bridge_gap='0.202um', array_return_link_width='2um',
        array_contact_lead_width='2um', array_contact_body_width='4um',
        array_contact_taper_start='6um', array_contact_taper_end='8um',
        array_contact_end='20um', array_contact_shadow_width='0.504um',
        jj_finger_width='80nm', jj_cross_electrode_width='92nm',
        jj_bridge_gap='0.248um', jj_return_lead_width='2um',
        shadow_clearance='0.075um', side='left',
        layer_jj='1', layer_jj_shadow='2', layer_array='3', layer_array_shadow='4',
        jj_metal='True',
    )
    component_metadata = Dict(short_name='fluxonium_junction',
                              _qgeometry_table_poly='True')

    def make(self):
        p = self.p
        u = self.parse_value('1um')
        n = float(p.num_array)
        if not isfinite(n) or not n.is_integer() or n < 1:
            raise ValueError('num_array must be a positive integer per arm')
        n = int(n)
        keys = ('array_element_length', 'array_row_width', 'array_bridge_gap',
                'array_return_link_width', 'array_contact_lead_width',
                'array_contact_body_width', 'array_contact_taper_start',
                'array_contact_taper_end', 'array_contact_end', 'array_contact_shadow_width',
                'jj_finger_width', 'jj_cross_electrode_width', 'jj_bridge_gap',
                'jj_return_lead_width', 'shadow_clearance')
        v = {k: float(p[k]) / u for k in keys}
        bad = [k for k, x in v.items() if not (isfinite(x) and x > 0)]
        if bad:
            raise ValueError(f'{", ".join(bad)} must be positive')
        side = str(p.side).strip().lower()
        if side not in ('left', 'right'):
            raise ValueError("side must be 'left' or 'right'")
        if v['array_row_width'] > 2:
            raise ValueError('array_row_width must be <= 2um for the 8um arm spacing')
        if not _flag(p.jj_metal, 'jj_metal'):
            return

        array, array_shadow, contact_shadow = self._array_and_contacts(n, v)
        if side == 'left':
            jj, jj_shadow = self._jj_left(v)
        else:
            turn = lambda g: {_turned_name(k): _rotate(s, 180, origin=(0, 0)) for k, s in g.items()}
            array, array_shadow, contact_shadow = map(turn, (array, array_shadow, contact_shadow))
            jj, jj_shadow = self._jj_right(v)

        groups = ((jj, p.layer_jj), ({**jj_shadow, **contact_shadow}, p.layer_jj_shadow),
                  (array, p.layer_array), (array_shadow, p.layer_array_shadow))
        for shapes, layer in groups:
            placed = {}
            for name, shape in shapes.items():
                if not shape.is_valid or shape.is_empty or shape.area <= 0:
                    raise ValueError(f'Invalid mask geometry: {name}')
                shape = draw.scale(shape, u, u, origin=(0, 0))
                shape = draw.rotate(shape, p.orientation, origin=(0, 0))
                placed[name] = draw.translate(shape, p.pos_x, p.pos_y)
            self.add_qgeometry('poly', placed, layer=layer, chip=p.chip, subtract=False)

        end, width = v['array_contact_end'], v['array_contact_body_width']
        for name, direction in (('top', 1), ('bottom', -1)):
            line = draw.LineString([(0, direction*(end - 1)), (0, direction*end)])
            line = draw.scale(line, u, u, origin=(0, 0))
            line = draw.rotate(line, p.orientation, origin=(0, 0))
            line = draw.translate(line, p.pos_x, p.pos_y)
            self.add_pin(name, np.asarray(line.coords), width=width*u,
                         input_as_norm=True, chip=p.chip)

    # ---- geometry in the side='left' frame (um) --------------------------
    def _array_and_contacts(self, n, v):
        L, W, gap = v['array_element_length'], v['array_row_width'], v['array_bridge_gap']
        ml = self.MOUTH_EDGE
        xm = ml + v['array_contact_lead_width']
        half = v['array_contact_body_width'] / 2
        y0, y1, ye = (v['array_contact_taper_start'], v['array_contact_taper_end'],
                      v['array_contact_end'])
        sw = v['array_contact_shadow_width']
        if not (5.1 < y0 < y1 < ye):
            raise ValueError('need 5.1um < array_contact_taper_start < _taper_end < array_contact_end')
        if half < xm or -half > ml:
            raise ValueError('array_contact_body_width must cover the funnel mouth')

        array = {}
        pitch = L + gap
        for arm, y in (('top', self.ARM_Y), ('bottom', -self.ARM_Y)):
            for i in range(n):
                right = ml - i*pitch
                array[f'{arm}_array_{i:04d}'] = draw.box(right - L, y - W/2, right, y + W/2)
        return_right = ml - L - (n - 1)*pitch - gap
        link_right = return_right - 1.008
        link_left = link_right - v['array_return_link_width']
        array['return_top'] = draw.box(link_left, 3, return_right, 5)
        array['return_link'] = draw.box(link_left, -3, link_right, 3)
        array['return_bottom'] = draw.box(link_left, -5, return_right, -3)
        top = _poly([(ml, 3), (xm, 3), (xm, y0), (half, y1), (half, ye),
                     (-half, ye), (-half, y1), (ml, y0)], 'top contact')
        array['contact_top'] = top
        array['contact_bottom'] = draw.scale(top, 1, -1, origin=(0, 0))

        bar_x0 = return_right - 2
        if bar_x0 >= -1.069:
            raise ValueError('num_array too small for the array shadow bars')
        array_shadow = {
            f'array_shadow_{arm}': draw.box(bar_x0, s*self.ARM_Y - W/2 - 0.05,
                                            -1.069, s*self.ARM_Y + W/2 + 0.05)
            for arm, s in (('top', 1), ('bottom', -1))}

        k = self.NECK_STRIP
        contact_shadow = {
            'return_shadow_outer': draw.box(link_left - 0.697, -5.05, link_left - 0.011, 5.049),
            'return_shadow_inner': draw.box(link_right - 0.010, -2.95, link_right + 0.496, 2.95),
        }
        for name, s, start_l, start_r in (('top', 1, 5.1, 3.6), ('bottom', -1, 5.1, 3.7)):
            left = [(ml, start_l), (ml, y0), (-half, y1), (-half, ye), (-half - sw, ye),
                    (-half - sw, y1), (ml - k, y0), (ml - k, start_l)]
            right = [(xm, start_r), (xm, y0), (half, y1), (half, ye), (half + sw, ye),
                     (half + sw, y1), (xm + k, y0), (xm + k, start_r)]
            for tag, pts in (('l', left), ('r', right)):
                contact_shadow[f'contact_shadow_{name}_{tag}'] = _poly(
                    [(x, s*y) for x, y in pts], f'{name} contact shadow')
        return array, array_shadow, contact_shadow

    def _finger_shadows(self, fr, fw, tip_y, el, eb, et, y_lo, y_hi, sc):
        """Landing, tip, below-electrode and bridge openings around one finger.

        Openings span y_lo .. y_hi (landing); the finger-side offsets are the
        NCKU values relative to the fixed finger right edge ``fr``, and the
        landing keeps at least ``sc`` from a wider finger.
        """
        lr = min(fr - self.LANDING_GAP, fr - fw - sc)
        fs = fr + self.FINGER_SIDE_GAP
        reach = fr + self.SHADOW_REACH
        return {
            'shadow_landing': _box(lr - 0.595, y_lo, lr, y_hi, 'landing opening'),
            'shadow_tip': _box(lr, y_lo, fs, tip_y - sc, 'tip opening'),
            'shadow_below': _box(fs, y_lo, reach, eb - sc, 'opening below electrode'),
            'shadow_bridge': _box(fs, eb - sc, el - sc, et + sc, 'bridge opening'),
        }, fs, reach

    def _jj_left(self, v):
        fw, cw, bg = v['jj_finger_width'], v['jj_cross_electrode_width'], v['jj_bridge_gap']
        rw, sc = v['jj_return_lead_width'], v['shadow_clearance']
        xm = self.MOUTH_EDGE + v['array_contact_lead_width']
        fr, xr = xm + self.FINGER_RIGHT, xm + self.RETURN_LEFT
        et, eb, el = self.BAR_TOP, self.BAR_TOP - cw, fr + bg
        if eb < self.FINGER_TIP:
            raise ValueError('jj_cross_electrode_width extends below the finger tip')
        jj = {
            'jj_top_bar': _box(xm, 3.0, fr, 3.5, 'top bar'),
            'jj_finger': _box(fr - fw, self.FINGER_TIP, fr, 3.0, 'finger'),
            'jj_cross_electrode': _box(el, eb, xr, et, 'cross electrode'),
            'jj_return_lead': _box(xr, -3.0, xr + rw, 3.412, 'return lead'),
            'jj_bottom_bar': _box(xm, -3.6, xr + rw, -3.0, 'bottom bar'),
        }
        shadow, fs, reach = self._finger_shadows(fr, fw, self.FINGER_TIP, el, eb, et,
                                                 1.759, 2.95, sc)
        shadow['shadow_above'] = _box(fs, et + sc, reach, 3.812, 'opening above electrode')
        shadow['shadow_return_below'] = _box(xr - 0.549, -2.95, xr - 0.049, eb - sc,
                                             'return strip below electrode')
        shadow['shadow_return_above'] = _box(xr - 0.549, et + sc, xr - 0.049, 3.812,
                                             'return strip above electrode')
        shadow['shadow_return_outer'] = _box(xr + rw + 0.054, -4.086, xr + rw + 0.554, 3.912,
                                             'outer return strip')
        return jj, shadow

    def _jj_right(self, v):
        fw, cw, bg = v['jj_finger_width'], v['jj_cross_electrode_width'], v['jj_bridge_gap']
        rw, sc = v['jj_return_lead_width'], v['shadow_clearance']
        xm = -(self.MOUTH_EDGE + v['array_contact_lead_width'])
        er = xm + self.R_ELECTRODE_RIGHT
        el = er - self.R_ELECTRODE_LENGTH
        fr = el - bg
        et, eb = self.R_BAR_TOP, self.R_BAR_TOP - cw
        tip, lead_bottom = self.R_FINGER
        if eb < tip:
            raise ValueError('jj_cross_electrode_width extends below the finger tip')
        jj = {
            'jj_top_bar': _box(fr - rw, 3.0, xm, 3.6, 'top bar'),
            'jj_return_lead': _box(fr - rw, lead_bottom, fr, 3.0, 'return lead'),
            'jj_finger': _box(fr - fw, tip, fr, lead_bottom, 'finger'),
            'jj_cross_electrode': _box(el, eb, er, et, 'cross electrode'),
            'jj_bottom_bar': _box(er, -3.5, xm, -3.0, 'bottom bar'),
        }
        shadow, fs, reach = self._finger_shadows(fr, fw, tip, el, eb, et,
                                                 -3.937, lead_bottom - 0.05, sc)
        shadow['shadow_above'] = _box(fs, et + sc, reach, 2.95, 'opening above electrode')
        shadow['shadow_return_outer'] = _box(fr - rw - 0.554, -3.912, fr - rw - 0.054, 4.086,
                                             'outer return strip')
        return jj, shadow


def _turned_name(name):
    """Name of a shape after the 180 deg turn: top <-> bottom, left (_l) <-> right (_r)."""
    swap = {'top': 'bottom', 'bottom': 'top', 'l': 'r', 'r': 'l'}
    parts = name.split('_')
    parts = [swap.get(t, t) if t in ('top', 'bottom') else t for t in parts]
    if name.startswith('contact_shadow_'):
        parts[-1] = swap[parts[-1]]
    return '_'.join(parts)


def _flag(value, name):
    text = str(value).lower()
    if text not in ('true', 'false', '1', '0'):
        raise ValueError(f'{name} must be True or False')
    return text in ('true', '1')


def _box(x0, y0, x1, y1, what):
    # draw.box silently reorders inverted corners, so check explicitly.
    if not (x1 - x0 > 1e-9 and y1 - y0 > 1e-9):
        raise ValueError(f'{what} has non-positive size for these options')
    return draw.box(x0, y0, x1, y1)


def _poly(points, what):
    shape = draw.Polygon(points)
    if not shape.is_valid or shape.area <= 0:
        raise ValueError(f'{what} is invalid for these options')
    return shape


fluxonium_junction = FluxoniumJunction
