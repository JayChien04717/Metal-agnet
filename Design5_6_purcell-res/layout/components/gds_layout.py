# -*- coding: utf-8 -*-
r"""GdsLayout: rebuild a qiskit-metal planar layout from a saved GDS.

Two source formats:
  * ``source_format='sim'``: GDS written by the qiskit-metal GDS renderer (negative ground plane on ``layer/0``,
    positive signal metal on ``layer/10`` polygons and ``layer/11`` paths);
  * ``source_format='metal'``: a mask with all positive metal on one layer/datatype (e.g. a delivered fab GDS, ``1/0``);
  * ``source_format='all'``: every datatype of the layer is positive metal (reader only, e.g. a fab render 1/10, 1/11, 1/100).
Metal = union of the source shapes; with ``fill_holes_um2 > 0`` interior holes smaller than that area (cheese) are
filled first. The component adds
  * ``metal_<i>``: every connected metal piece except the one touching the chip edge (the main ground), subtract=False;
  * ``etch_<i>``: each connected piece of chip box - metal, subtract=True.
Rendering it therefore gives the same metal as the source. There are no pins or per-part names (the source layout
scripts no longer exist, Design5_6 README 2026-11-25). GDS units are read from the file (um); design units are mm.
"""
from __future__ import annotations

from pathlib import Path

from qiskit_metal import Dict
from qiskit_metal.qlibrary.core import QComponent
from shapely.geometry import Polygon, box
from shapely.ops import unary_union


def _pieces(geom):
    return [g for g in getattr(geom, 'geoms', [geom]) if not g.is_empty and g.area > 0]


def fill_small_holes(geom, max_area):
    return unary_union([Polygon(g.exterior, [h for h in g.interiors if Polygon(h).area > max_area])
                        for g in _pieces(geom)])


def read_gds_metal(gds_file, source_format='sim', layer=1, datatype=0, cell=None, fill_holes_um2=0.0):
    """Return (chip_box, metal) shapely geometries in mm."""
    import gdstk
    lib = gdstk.read_gds(str(gds_file))
    top = next(c for c in lib.cells if c.name == cell) if cell else lib.top_level()[0]
    scale = lib.unit / 1e-3                       # GDS user unit -> mm
    polys = top.get_polygons()                    # flattened, all references applied
    shp = lambda sel: unary_union([Polygon(p.points * scale).buffer(0) for p in polys if sel(p)])
    if source_format == 'sim':
        metal = shp(lambda p: p.layer == layer and p.datatype in (0, 10, 11))
    elif source_format == 'metal':
        metal = shp(lambda p: p.layer == layer and p.datatype == datatype)
    elif source_format == 'all':          # every datatype of the layer is positive metal (e.g. a fab render)
        metal = shp(lambda p: p.layer == layer)
    else:
        raise ValueError(f'source_format {source_format!r}')
    if fill_holes_um2 > 0:
        metal = fill_small_holes(metal, fill_holes_um2 * 1e-6)     # um^2 -> mm^2
    (x0, y0), (x1, y1) = top.bounding_box()
    return box(x0 * scale, y0 * scale, x1 * scale, y1 * scale), metal


def read_gds_regions(gds_file, layer=1, cell=None, ground_dt=0, metal_dts=(10, 11)):
    """Backward-compatible reader for a qiskit-metal sim GDS: (chip_box, ground, signal metal) in mm."""
    import gdstk
    lib = gdstk.read_gds(str(gds_file))
    top = next(c for c in lib.cells if c.name == cell) if cell else lib.top_level()[0]
    scale = lib.unit / 1e-3
    polys = top.get_polygons()
    to_shp = lambda p: Polygon(p.points * scale).buffer(0)
    ground = unary_union([to_shp(p) for p in polys if p.layer == layer and p.datatype == ground_dt])
    metal = unary_union([to_shp(p) for p in polys if p.layer == layer and p.datatype in metal_dts])
    (x0, y0), (x1, y1) = top.bounding_box()
    return box(x0 * scale, y0 * scale, x1 * scale, y1 * scale), ground, metal


class GdsLayout(QComponent):
    """Whole-chip layout imported from a saved GDS (see module docstring)."""

    default_options = Dict(gds_file='', cell='', source_format='sim', layer='1', datatype='0',
                           fill_holes_um2='0', chip='main')
    component_metadata = Dict(short_name='gds', _qgeometry_table_poly='True')

    def make(self):
        p = self.options
        if not p.gds_file or not Path(p.gds_file).exists():
            raise FileNotFoundError(f'GdsLayout: gds_file not found: {p.gds_file!r}')
        chip, metal = read_gds_metal(p.gds_file, source_format=p.source_format, layer=int(p.layer),
                                     datatype=int(p.datatype), cell=p.cell or None,
                                     fill_holes_um2=float(p.fill_holes_um2))
        edge = chip.exterior
        islands = [g for g in _pieces(metal) if not g.intersects(edge)]
        etch = chip.difference(metal)
        self.add_qgeometry('poly', {f'metal_{i}': g for i, g in enumerate(islands)}, subtract=False,
                           chip=p.chip, layer=int(p.layer))
        self.add_qgeometry('poly', {f'etch_{i}': g for i, g in enumerate(_pieces(etch))}, subtract=True,
                           chip=p.chip, layer=int(p.layer))
