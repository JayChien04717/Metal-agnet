# -*- coding: utf-8 -*-
r"""Junction test coupon: two floating probe pads with facing landing stems in an etched pocket.

The fab draws the junction across the landing gap. The two styles (geometry inputs in
DesignTestJunction/parameters/geometry.json, placed by DesignTestJunction/layout/build.py
``add_test_junction``):

* F (fluxonium style): the two Dolan recesses face each other (``hook_offset`` 0);
  30 um landing gap with room for an array.
* T (transmon style): the right recess is shifted along y (``hook_offset``); 5 um gap.

Sketch (local frame, orientation 0; pads and stems along x)::

    ________________________ pocket (etched) ________________________
   |    ____________                              ____________       |
   |   |   pad_l    |=====[hook]  gap  [hook]=====|   pad_r    |     |
   |   |____________|     stem              stem  |____________|     |
   |__________________________________________________________________|

``fab_contacts=True`` cuts the 2 x ``hook_depth`` Dolan recesses into the stem tips.
``label='auto'`` (default) etches an upright ``F`` or ``T`` label at the upper-right
of the pocket (``T`` when ``hook_offset != 0`` or ``landing_gap <= 15 um``, else ``F``;
pass explicit ``'F'``, ``'T'``, or ``''`` to disable).
Neither pad is connected to ground. Moved unchanged on 2026-11-20 from the retired
``pattern design/test_junction.py`` (whose old corner-pocket ``TestJunction`` was removed).
"""
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import QComponent
from xmark import text_geometry


class HookedTestJunction(QComponent):
    """Two floating probe pads with facing landing stems, optional Dolan recesses,
    and an etched upright F/T style label at the upper-right of the pocket.

    In local coordinates the pad length and stems run along x. The left hook
    is centred at y=0; the right hook can be shifted along y. The gap stays open
    for the fab's junction and neither landing is connected to ground.
    """

    default_options = Dict(
        fab_contacts='False',
        hook_offset='0um',
        label='auto',
        label_size='150um',
        label_gap='40um',
    )
    component_metadata = Dict(short_name='hooktestjj', _qgeometry_table_poly='True')

    def _resolved_label(self):
        p = self.p
        raw = '' if p.label is None else str(p.label).strip()
        if raw.lower() in ('', 'none', 'false'):
            return ''
        if raw.lower() == 'auto':
            return 'T' if (abs(p.hook_offset) > 1e-9 or p.landing_gap <= 0.015 + 1e-9) else 'F'
        return raw

    def _label_geom(self, pocket):
        text = self._resolved_label()
        if not text:
            return None
        p = self.p
        g = text_geometry(text, p.label_size)
        gx0, _, _, gy1 = g.bounds
        _, _, px1, py1 = pocket.bounds
        return draw.translate(g, px1 + p.label_gap - gx0, py1 - gy1)

    def make(self):
        p = self.p
        if not (
            0 < p.pad_corner_radius < min(p.pad_width, p.pad_length) / 2
            and 0 < p.stem_width < p.pad_width - 2 * p.pad_corner_radius
            and 0 < p.hook_depth < (p.pad_body_gap - p.landing_gap) / 2
            and p.landing_gap > 0
            and 0 < p.hook_width and abs(p.hook_offset) + p.hook_width / 2 < p.stem_width / 2
            and p.pocket_width > 2 * p.pad_length + p.pad_body_gap
            and p.pocket_height > p.pad_width
        ):
            raise ValueError('Invalid pad, pocket, stem or hook dimensions for HookedTestJunction')

        pocket = draw.rectangle(p.pocket_width, p.pocket_height, 0, 0)
        y = p.stem_width / 2
        body_edge, tip = p.pad_body_gap / 2, p.landing_gap / 2
        radius = p.pad_corner_radius
        pad_r = draw.rectangle(p.pad_length - 2 * radius, p.pad_width - 2 * radius,
                               body_edge + p.pad_length / 2, 0).buffer(radius)
        pad_r = pad_r.union(draw.box(tip, -y, body_edge + radius, y))
        pad_l = draw.scale(pad_r, -1, 1, origin=(0, 0))

        if str(p.fab_contacts).lower() in ('true', '1'):
            w = p.hook_width / 2
            pad_l = pad_l.difference(draw.box(-tip - p.hook_depth, -w, -tip, w))
            pad_r = pad_r.difference(draw.box(tip, p.hook_offset - w,
                                              tip + p.hook_depth, p.hook_offset + w))

        pocket, pad_l, pad_r = draw.rotate([pocket, pad_l, pad_r], p.orientation, origin=(0, 0))
        pocket, pad_l, pad_r = draw.translate([pocket, pad_l, pad_r], p.pos_x, p.pos_y)
        subtract_geoms = dict(pocket=pocket)
        lbl = self._label_geom(pocket)
        if lbl is not None:
            subtract_geoms['label'] = lbl
        self.add_qgeometry('poly', subtract_geoms, subtract=True, layer=p.layer)
        self.add_qgeometry('poly', dict(pad_l=pad_l, pad_r=pad_r), layer=p.layer)

    def keepouts(self):
        """Etched pocket (and F/T label bbox, if enabled) keep-out polygons in design units (mm)."""
        p = self.p
        pocket = draw.rectangle(p.pocket_width, p.pocket_height, 0, 0)
        pocket = draw.rotate(pocket, p.orientation, origin=(0, 0))
        pocket = draw.translate(pocket, p.pos_x, p.pos_y)
        out = {'pocket': pocket}
        lbl = self._label_geom(pocket)
        if lbl is not None:
            out['label'] = draw.box(*lbl.bounds)
        return out
