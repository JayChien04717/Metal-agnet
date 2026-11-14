# -*- coding: utf-8 -*-
r"""XMark: metal cross in an etched square window (die corner / alignment mark), with an
optional die-ID label next to it, at the lower right.

     ground (metal)
      ______________
     |      |       |
     |      |       |   <- window (etched, window_size square)
     |------+-------|   <- metal cross, arm_length x arm_width bars
     |      |       |   _________
     |______|_______|  | FTF07   |  <- label, etched in the ground, label_gap from
                       |_________|     the window edge, top label_gap below the cross

Origin = cross centre.

    from xmark import XMark
    XMark(design, 'xm_UL', options=dict(pos_x='-4.65mm', pos_y='4.65mm', label='FTF07'))
"""
import gdspy
import shapely
from shapely.affinity import translate
from shapely.geometry import Polygon

from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import QComponent


def text_geometry(text, size, spacing=None):
    """Text as a shapely (Multi)Polygon, gdspy's built-in font, character height `size`
    (any unit, the result has the same unit). Bounding box starts at (0, 0).
    The characters are packed tightly: `spacing` between neighbouring glyphs (default one
    font pixel, size / 9) instead of gdspy's wide fixed advance."""
    spacing = size / 9 if spacing is None else spacing
    glyphs, x = [], 0.0
    for ch in text:
        if ch == ' ':
            x += size / 3
            continue
        g = shapely.union_all([Polygon(p) for p in gdspy.Text(ch, size, (0, 0)).polygons])
        gx0, _, gx1, _ = g.bounds
        glyphs.append(translate(g, x - gx0, 0))       # same baseline for all characters
        x += gx1 - gx0 + spacing
    g = shapely.union_all(glyphs)
    x0, y0, _, _ = g.bounds
    return translate(g, -x0, -y0)


def label_offset(arm_width, window_size, label_gap):
    """(x_left, y_top) of the label relative to the cross centre: right of the window,
    label_gap from its edge, top label_gap below the horizontal arm."""
    return max(window_size, arm_width) / 2 + label_gap, -(arm_width / 2 + label_gap)


def label_geometry(text, size, arm_width, window_size, label_gap):
    """Label polygon placed relative to the cross centre (before rotation / translation)."""
    g = text_geometry(text, size)
    x, y = label_offset(arm_width, window_size, label_gap)
    return translate(g, x, y - g.bounds[3])


CORNERS = (('UL', -1, 1), ('UR', 1, 1), ('LL', -1, -1), ('LR', 1, -1))


def corner_mark_centres(chip_x, chip_y, inset):
    """{'UL': (x, y), ...} cross centres `inset` in from each chip corner (chip centred at 0)."""
    return {n: (sx * (chip_x / 2 - inset), sy * (chip_y / 2 - inset)) for n, sx, sy in CORNERS}


def add_corner_marks(design, die_id='', inset=0.35, label_corner='UL', chip='main', **options):
    """Four XMarks `inset` mm in from the chip corners; the die ID goes next to the
    `label_corner` cross. Returns {corner: (x, y)} (mm) of the cross centres."""
    size = design.chips[chip].size
    cx = float(str(size.size_x).replace('mm', ''))
    cy = float(str(size.size_y).replace('mm', ''))
    centres = corner_mark_centres(cx, cy, inset)
    for n, (x, y) in centres.items():
        XMark(design, f'xm_{n}', options=dict(
            chip=chip, pos_x=f'{x}mm', pos_y=f'{y}mm',
            label=die_id if n == label_corner else '', **options))
    return centres


class XMark(QComponent):
    """Metal cross in an etched square window, with optional etched label.

    Default options:
        * arm_length: length of each bar of the cross (250 um)
        * arm_width: width of each bar (10 um)
        * window_size: side of the etched square around the cross (400 um)
        * label: text right of the window, below the cross line ('' = no label)
        * label_size: character height passed to gdspy.Text (capitals / digits are ~0.78 of it)
        * label_gap: distance from the label to the window edge and to the cross line
        * label_left: place the label on the left side of the window, keeping it upright
        * label_subtract: True = label etched in the ground, False = metal
    """
    default_options = Dict(
        arm_length='250um', arm_width='10um', window_size='400um',
        label='', label_size='100um', label_gap='30um', label_subtract='True',
        label_left='False',
        layer='1',
    )
    component_metadata = Dict(short_name='xm')

    def make(self):
        p = self.p
        window = draw.rectangle(p.window_size, p.window_size, 0, 0)
        cross = shapely.union_all([draw.rectangle(p.arm_length, p.arm_width, 0, 0),
                                   draw.rectangle(p.arm_width, p.arm_length, 0, 0)])
        geoms = [window, cross]
        if str(p.label):
            label = label_geometry(str(p.label), p.label_size, p.arm_width,
                                   p.window_size, p.label_gap)
            if str(p.label_left).lower() in ('true', '1'):
                label = translate(label, -label.bounds[0] - label.bounds[2], 0)
            geoms.append(label)
        geoms = draw.rotate(geoms, p.orientation, origin=(0, 0))
        geoms = draw.translate(geoms, p.pos_x, p.pos_y)
        self.add_qgeometry('poly', dict(window=geoms[0]), subtract=True, layer=p.layer)
        self.add_qgeometry('poly', dict(cross=geoms[1]), layer=p.layer)
        if len(geoms) > 2:
            self.add_qgeometry('poly', dict(label=geoms[2]), subtract=p.label_subtract, layer=p.layer)
