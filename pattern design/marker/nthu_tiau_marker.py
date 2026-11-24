# -*- coding: utf-8 -*-
r"""NTHU TiAu e-beam alignment markers (pre-fabricated NTHU TiAu photomask), reproduced from
`Fabrication output/Template/NTHU ebeam gds/`:

* `TiAuChipMarker`: the 4 chip-marker sites of one 10 x 10 mm die (layer 2/0, 36 polygons),
  reference `TiAu_marker_chip_single_die.gds`
* `TiAuWaferMarker`: the 8 wafer-marker sites of the 4-inch wafer (layer 3/0, 24 polygons),
  reference `TiAu_marker_wafer_4inch.gds`

The photomask already exists, so these coordinates are fixed: use them to draw / check a layout
around the markers, never to change the marker artwork. Grid and site map:
note/tiau_marker_alignment.md.

Chip-marker site, UR prototype (um, die centre at (0, 0); 84 x 68 um, centre (2507.25, 4100)):

      y=4134 |  .  .  /|  .  .  .  .  .  .  .  .  .  .
             |  _|_  /_|              ___           <- triangle 5 um, 20 um cross (3 um bars)
             |   |               ____|   |  <- L vernier squares (5 um)
      y=4100 |                --------+--------      <- 40 um cross (5 um bars),
             |  [] _|_ []             |                 centre (2529.25, 4102.5)
             |      |
      y=4066 |     []                                <- 5 um squares
            x=2465.25                          x=2549.25

The other sites are mirror images: UL (-x, y), LL (-x, -y), LR (x, -y).

Wafer-marker site (um, relative to the site centre; 2300 x 3155 um):

             +---------+  label block 430 x 308 um with 16 slots (top 1269.5..1577.5)
             |  ::  :: |
             +---------+
                 _|_     550 x 450 um cross, 150 um bars, centre (0, 922.5)
                  |
            ______|______
           |      |      |   2300 x 1800 um cross, 500 um bars, centre (0, -677.5)
           |______|______|
                  |

All 8 sites are translated (not mirrored) copies, centres (+-22500, 37547.5),
(+-22500, -37432.5), (+-37350, 23297.5), (+-37350, -22702.5) um from the wafer centre.

Each site needs a bare-substrate opening in base metal, site bbox + 50 um per side
(`ground_opening`, default on).

    from nthu_tiau_marker import TiAuChipMarker
    TiAuChipMarker(design, 'tiau', options=dict(pos_x='0mm', pos_y='0mm'))
"""
import shapely
from shapely.affinity import rotate, scale, translate
from shapely.geometry import Polygon, box

from qiskit_metal import Dict
from qiskit_metal.qlibrary.core import QComponent

CHIP_LAYER = (2, 0)
WAFER_LAYER = (3, 0)
QUADRANTS = (('UR', 1, 1), ('UL', -1, 1), ('LL', -1, -1), ('LR', 1, -1))
CHIP_SITE_CENTRE_UM = (2507.25, 4100.0)
WAFER_SITE_CENTRES_UM = {
    'top_L': (-22500.0, 37547.5), 'top_R': (22500.0, 37547.5),
    'bottom_L': (-22500.0, -37432.5), 'bottom_R': (22500.0, -37432.5),
    'upper_L': (-37350.0, 23297.5), 'upper_R': (37350.0, 23297.5),
    'lower_L': (-37350.0, -22702.5), 'lower_R': (37350.0, -22702.5),
}
MARGIN_UM = 50.0


def cross(cx, cy, length_x, length_y, width):
    """Plus-shaped cross (one polygon): horizontal bar length_x x width, vertical bar
    width x length_y, centre (cx, cy)."""
    return shapely.union_all([box(cx - length_x / 2, cy - width / 2, cx + length_x / 2, cy + width / 2),
                              box(cx - width / 2, cy - length_y / 2, cx + width / 2, cy + length_y / 2)])


def chip_site_ur():
    """The 9 polygons of the UR chip-marker site (um, die centre at (0, 0))."""
    return [
        cross(2529.25, 4102.5, 40.0, 40.0, 5.0),                       # 40 um cross
        cross(2476.75, 4077.5, 20.0, 20.0, 3.0),                       # 20 um crosses
        cross(2476.75, 4122.5, 20.0, 20.0, 3.0),
        Polygon([(2483.25, 4129.0), (2483.25, 4134.0), (2488.25, 4129.0)]),   # triangle
        box(2539.25, 4117.5, 2544.25, 4122.5),                         # L vernier
        box(2539.25, 4112.5, 2549.25, 4117.5),
        box(2465.25, 4084.0, 2470.25, 4089.0),                         # squares
        box(2483.25, 4084.0, 2488.25, 4089.0),
        box(2483.25, 4066.0, 2488.25, 4071.0),
    ]


def chip_marker_sites():
    """{'UR'|'UL'|'LL'|'LR': [9 polygons]} of one die (um, die centre at (0, 0))."""
    proto = chip_site_ur()
    return {n: [scale(g, sx, sy, origin=(0, 0)) for g in proto] for n, sx, sy in QUADRANTS}


def wafer_site_proto():
    """The 3 polygons of one wafer-marker site (um, relative to the site centre)."""
    big_y = -677.5
    label = box(-215.0, big_y + 1947.0, 215.0, big_y + 2255.0)
    slots = []
    for y0, y1, xs in ((1967, 1998, ((-173, -132), (-118, -77), (77, 118), (132, 173))),
                       (2012, 2043, ((-173, -132), (-118, -77), (77, 118), (132, 173))),
                       (2125, 2165, ((-195, -140), (-110, -55), (55, 110), (140, 195))),
                       (2195, 2235, ((-195, -140), (-110, -55), (55, 110), (140, 195)))):
        slots += [box(x0, big_y + y0, x1, big_y + y1) for x0, x1 in xs]
    return [
        cross(0.0, big_y, 2300.0, 1800.0, 500.0),                      # large cross
        cross(0.0, big_y + 1600.0, 550.0, 450.0, 150.0),               # small cross
        label.difference(shapely.union_all(slots)),                    # label block
    ]


def wafer_marker_sites():
    """{site name: [3 polygons]} of the 4-inch wafer (um, wafer centre at (0, 0))."""
    proto = wafer_site_proto()
    return {n: [translate(g, x, y) for g in proto] for n, (x, y) in WAFER_SITE_CENTRES_UM.items()}


def site_opening(polys, margin=MARGIN_UM):
    """Bare-substrate rectangle: bbox of the site polygons + margin per side."""
    x0, y0, x1, y1 = shapely.union_all(polys).bounds
    return box(x0 - margin, y0 - margin, x1 + margin, y1 + margin)


class _SiteMarker(QComponent):
    """Shared make(): place {site: [polygons in um]} at pos / orientation, plus openings."""
    _sites = staticmethod(dict)

    def _selected_sites(self):
        p = self.p
        wanted = [s.strip() for s in str(p.sites).split(',') if s.strip()]
        sites = self._sites()
        unknown = set(wanted) - set(sites) - {'all'}
        if unknown:
            raise ValueError(f'{self.name}: unknown sites {sorted(unknown)}; choose from {list(sites)}')
        if 'all' not in wanted:
            sites = {n: gs for n, gs in sites.items() if n in wanted}
        return sites

    def _place(self, g):
        p = self.p
        g = scale(g, 1e-3, 1e-3, origin=(0, 0))
        g = rotate(g, p.orientation, origin=(0, 0))
        return translate(g, p.pos_x, p.pos_y)

    def keepouts(self):
        """Bare-rectangle keep-out polygons {site_name: Polygon} in design units (mm)."""
        m = self.p.opening_margin * 1e3
        return {n: self._place(site_opening(gs, m)) for n, gs in self._selected_sites().items()}

    def make(self):
        p = self.p
        sites = self._selected_sites()
        self.add_qgeometry('poly', {f'{n}_{i}': self._place(g) for n, gs in sites.items()
                                    for i, g in enumerate(gs)}, layer=p.layer)
        if str(p.ground_opening).lower() in ('true', '1'):
            self.add_qgeometry('poly', {f'opening_{n}': poly for n, poly in self.keepouts().items()},
                               subtract=True, layer=p.ground_layer)


class TiAuChipMarker(_SiteMarker):
    """NTHU TiAu chip markers of one die; pos = die centre. See the module docstring.

    Default options:
        * sites: 'all' or a comma list of UR, UL, LL, LR
        * layer: marker layer (2 = NTHU TiAu mask chip-marker layer)
        * ground_opening: cut site bbox + opening_margin (50 um) in ground_layer (1)
    """
    default_options = Dict(sites='all', layer='2', ground_opening='True',
                           opening_margin=f'{MARGIN_UM:g}um', ground_layer='1')
    component_metadata = Dict(short_name='tiauc')
    _sites = staticmethod(chip_marker_sites)


class TiAuWaferMarker(_SiteMarker):
    """NTHU TiAu wafer markers of the 4-inch wafer; pos = wafer centre (centre die at the
    origin of the 10.2 mm grid). See the module docstring.

    Default options:
        * sites: 'all' or a comma list of top_L, top_R, bottom_L, bottom_R, upper_L, upper_R,
          lower_L, lower_R
        * layer: marker layer (3 = NTHU TiAu mask wafer-marker layer)
        * ground_opening: cut site bbox + opening_margin (50 um) in ground_layer (1)
    """
    default_options = Dict(sites='all', layer='3', ground_opening='True',
                           opening_margin=f'{MARGIN_UM:g}um', ground_layer='1')
    component_metadata = Dict(short_name='tiauw')
    _sites = staticmethod(wafer_marker_sites)
