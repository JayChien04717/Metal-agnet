# -*- coding: utf-8 -*-
r"""TSRI A12 e-beam alignment marks inside one chip: the fine-cross groups and the four coarse
crosses of TSRI's A12 field (`Fabrication output/Template/TSRI ebeam gds/A12mark_TSRI.gds`).

The A12 reference also has two corner L marks at (+-6000, +-6000) um. They mark the dicing
street, outside a 10 mm chip, and are drawn with the street artwork, not by this component
(`a12_groups(include_corners=True)` still returns them for the reference check).

Chip (origin = chip / field centre, um; not to scale; everything within +-4800 um):

    +--------------------------------------------------+ (5000, 5000) chip edge
    |      |                                    |      |
    |   ---+--- coarse cross            coarse ---+--- |  <- 4 coarse crosses, centres
    |      |    (-3600, 4200)    (3600, 4200)   |      |     (+-3600, +-4200), 1200 x 4 um bars
    |    ++++                                ++++      |  <- 4 fine groups: 4 x 4 crosses,
    |    ++++  x, y in +-{3100 ... 3400}     ++++      |     pitch 100 um, 50 x 4 um bars
    |    ++++                                ++++      |
    |    ++++                                ++++      |
    |                       o (0, 0)                   |
    |    ++++                                ++++      |
    |    ++++                                ++++      |
    |   ---+---                             ---+---    |
    |      |                                    |      |
    +--------------------------------------------------+ (-5000, -5000)

* fine crosses: centres at x, y in +-{3100, 3200, 3300, 3400} um (64 crosses), bars
  `fine_arm` (50 um) x `fine_width` (4 um)
* coarse crosses: centres (+-3600, +-4200) um, bars `coarse_arm` (1200 um) x `coarse_width` (4 um)
* bare rectangle in `ground_layer` around every group: group bbox + `opening_margin` (50 um):
  450 x 450 um per fine group, 1300 x 1300 um per coarse cross

Layer 2 and the 50 um margin follow the TSRI plan in skill/fabrication_output.md (not yet
confirmed by TSRI).

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
               include_corners=False, corner_leg=500.0, corner_width=50.0):
    """{group name: [shapely polygons]} of the A12 field in um, centred at (0, 0).
    Groups: fine_<UL|UR|LL|LR>, coarse_<..>; with include_corners also the dicing-street corner
    L marks corner_<UR|LL> of the reference (not part of the chip marks)."""
    quads = (('UL', -1, 1), ('UR', 1, 1), ('LL', -1, -1), ('LR', 1, -1))
    groups = {}
    for name, sx, sy in quads:
        groups[f'fine_{name}'] = [cross(sx * x, sy * y, fine_arm, fine_width)
                                  for y in FINE_CENTRES_UM for x in FINE_CENTRES_UM]
    for name, sx, sy in quads:
        cx, cy = sx * COARSE_CENTRE_UM[0], sy * COARSE_CENTRE_UM[1]
        groups[f'coarse_{name}'] = [cross(cx, cy, coarse_arm, coarse_width)]
    if include_corners:
        c, leg, w = CORNER_UM, corner_leg, corner_width
        ur = shapely.union_all([box(c - w, c - leg, c, c), box(c - leg, c - w, c, c)])
        groups['corner_UR'] = [ur]
        groups['corner_LL'] = [rotate(ur, 180, origin=(0, 0))]
    return groups


def a12_polygons(**kwargs):
    """A12 polygons (um, centred at (0, 0)) as a flat list; kwargs as a12_groups()."""
    return [g for gs in a12_groups(**kwargs).values() for g in gs]


class TSRIA12Mark(QComponent):
    """TSRI A12 chip marks (fine groups + coarse crosses), see the module docstring.

    Default options:
        * fine_arm, fine_width: fine cross bar length / width (50 um, 4 um)
        * coarse_arm, coarse_width: coarse cross bar length / width (1200 um, 4 um)
        * layer: marker layer (2, the TSRI chip-mark layer of the e-beam GDS)
        * include_marks: draw the cross artwork ('False' keeps only the ground openings)
        * ground_opening: cut a bare rectangle (group bbox + opening_margin) in ground_layer
          around every group
        * opening_margin, ground_layer: bare-rectangle margin per side (50 um), layer to cut (1)
    """
    default_options = Dict(
        fine_arm='50um', fine_width='4um', coarse_arm='1200um', coarse_width='4um',
        layer='2', include_marks='True', ground_opening='True', opening_margin='50um', ground_layer='1',
    )
    component_metadata = Dict(short_name='a12')

    def _place(self, g):
        p = self.p
        g = scale(g, 1e-3, 1e-3, origin=(0, 0))
        g = rotate(g, p.orientation, origin=(0, 0))
        return translate(g, p.pos_x, p.pos_y)

    def _groups(self):
        p = self.p
        um = 1e3     # parsed options are in mm
        return a12_groups(p.fine_arm * um, p.fine_width * um, p.coarse_arm * um, p.coarse_width * um)

    def keepouts(self):
        """Bare-rectangle keep-out polygons {group_name: Polygon} in design units (mm)."""
        p = self.p
        m = p.opening_margin * 1e3
        out = {}
        for name, gs in self._groups().items():
            x0, y0, x1, y1 = shapely.union_all(gs).bounds
            out[name] = self._place(box(x0 - m, y0 - m, x1 + m, y1 + m))
        return out

    def make(self):
        p = self.p
        groups = self._groups()
        if str(p.include_marks).lower() in ('true', '1'):
            marks = {f'{name}_{i}': self._place(g) for name, gs in groups.items() for i, g in enumerate(gs)}
            self.add_qgeometry('poly', marks, layer=p.layer)
        if str(p.ground_opening).lower() in ('true', '1'):
            openings = {f'opening_{name}': poly for name, poly in self.keepouts().items()}
            self.add_qgeometry('poly', openings, subtract=True, layer=p.ground_layer)
