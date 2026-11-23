# -*- coding: utf-8 -*-
r"""TSRI A12 e-beam alignment marks: one 12 x 12 mm field, reproduced from TSRI's
`Fabrication output/Template/TSRI ebeam gds/A12mark_TSRI.gds` (layer 0/0, 76 polygons).

Field (origin = field centre, um; not to scale):

    +----------------------------------------------------------+ (6000, 6000)
    |                                                    ___|  <- corner L (UR only):
    |        coarse +                     + coarse          |     2 bars 500 x 50 um
    |     (-3600, 4200)                 (3600, 4200)        |
    |       ++++                           ++++             |  <- fine group: 4 x 4 crosses,
    |       ++++  x,y in -3400..-3100      ++++ 3100..3400  |     pitch 100 um, 50 x 4 um arms
    |       ++++                           ++++             |
    |       ++++                           ++++             |
    |                          o (0, 0)                     |
    |       ++++                           ++++             |
    |       ++++                           ++++             |
    |     coarse + (-3600, -4200)       coarse + (3600, -4200)
    | |___                                                   |  <- corner L (LL only)
    +----------------------------------------------------------+
  (-6000, -6000)

* fine crosses: centres at x, y in +-{3100, 3200, 3300, 3400} um (64 crosses), arms
  `fine_arm` (50 um) x `fine_width` (4 um)
* coarse crosses: centres (+-3600, +-4200) um, arms `coarse_arm` (1200 um) x `coarse_width`
  (4 um)
* corner L marks at (6000, 6000) and (-6000, -6000) um, legs `corner_leg` (500 um) x
  `corner_width` (50 um), pointing into the field

Polarity, target layer, per-die vs. per-field placement and the bare-rectangle margin are still
pending with TSRI (skill/fabrication_output.md); `layer` defaults to the reference layer 0 and
no ground opening is cut unless `ground_opening='True'`.

    from tsri_a12_mark import TSRIA12Mark
    TSRIA12Mark(design, 'a12', options=dict(pos_x='0mm', pos_y='0mm'))
"""
import shapely
from shapely.affinity import rotate, scale, translate
from shapely.geometry import box

from qiskit_metal import Dict
from qiskit_metal.qlibrary.core import QComponent

FINE_CENTRES_UM = (3100.0, 3200.0, 3300.0, 3400.0)
COARSE_CENTRE_UM = (3600.0, 4200.0)
CORNER_UM = 6000.0
REFERENCE_LAYER = (0, 0)


def cross(cx, cy, arm, width):
    """Plus-shaped cross (one polygon), centre (cx, cy), bar length `arm`, bar width `width`."""
    return shapely.union_all([box(cx - arm / 2, cy - width / 2, cx + arm / 2, cy + width / 2),
                              box(cx - width / 2, cy - arm / 2, cx + width / 2, cy + arm / 2)])


def a12_groups(fine_arm=50.0, fine_width=4.0, coarse_arm=1200.0, coarse_width=4.0,
               corner_leg=500.0, corner_width=50.0):
    """{group name: [shapely polygons]} of one A12 field in um, field centred at (0, 0).
    Groups: fine_<UL|UR|LL|LR>, coarse_<..>, corner_<UR|LL>."""
    quads = (('UL', -1, 1), ('UR', 1, 1), ('LL', -1, -1), ('LR', 1, -1))
    groups = {}
    for name, sx, sy in quads:
        groups[f'fine_{name}'] = [cross(sx * x, sy * y, fine_arm, fine_width)
                                  for y in FINE_CENTRES_UM for x in FINE_CENTRES_UM]
    for name, sx, sy in quads:
        cx, cy = sx * COARSE_CENTRE_UM[0], sy * COARSE_CENTRE_UM[1]
        groups[f'coarse_{name}'] = [cross(cx, cy, coarse_arm, coarse_width)]
    c, leg, w = CORNER_UM, corner_leg, corner_width
    ur = shapely.union_all([box(c - w, c - leg, c, c), box(c - leg, c - w, c, c)])
    groups['corner_UR'] = [ur]
    groups['corner_LL'] = [rotate(ur, 180, origin=(0, 0))]
    return groups


def a12_polygons(**kwargs):
    """All A12 polygons (um, field centred at (0, 0)) as a flat list."""
    return [g for gs in a12_groups(**kwargs).values() for g in gs]


class TSRIA12Mark(QComponent):
    """One TSRI A12 mark field (12 x 12 mm), see the module docstring for the sketch.

    Default options:
        * fine_arm, fine_width: fine cross bar length / width (50 um, 4 um)
        * coarse_arm, coarse_width: coarse cross bar length / width (1200 um, 4 um)
        * corner_leg, corner_width: corner L leg length / width (500 um, 50 um)
        * layer: marker layer (0 = reference layer; final layer pending with TSRI)
        * ground_opening: 'True' cuts a bare rectangle (group bbox + opening_margin) in
          ground_layer around every group; off by default while the margin is pending
        * opening_margin, ground_layer: bare-rectangle margin per side, layer to cut
    """
    default_options = Dict(
        fine_arm='50um', fine_width='4um', coarse_arm='1200um', coarse_width='4um',
        corner_leg='500um', corner_width='50um',
        layer='0', ground_opening='False', opening_margin='50um', ground_layer='1',
    )
    component_metadata = Dict(short_name='a12')

    def make(self):
        p = self.p
        um = 1e3     # parsed options are in mm
        groups = a12_groups(p.fine_arm * um, p.fine_width * um, p.coarse_arm * um,
                            p.coarse_width * um, p.corner_leg * um, p.corner_width * um)

        def place(g):
            g = scale(g, 1e-3, 1e-3, origin=(0, 0))
            g = rotate(g, p.orientation, origin=(0, 0))
            return translate(g, p.pos_x, p.pos_y)

        marks = {f'{name}_{i}': place(g) for name, gs in groups.items() for i, g in enumerate(gs)}
        self.add_qgeometry('poly', marks, layer=p.layer)
        if str(p.ground_opening).lower() in ('true', '1'):
            m = p.opening_margin * um
            openings = {}
            for name, gs in groups.items():
                x0, y0, x1, y1 = shapely.union_all(gs).bounds
                openings[f'opening_{name}'] = place(box(x0 - m, y0 - m, x1 + m, y1 + m))
            self.add_qgeometry('poly', openings, subtract=True, layer=p.ground_layer)
