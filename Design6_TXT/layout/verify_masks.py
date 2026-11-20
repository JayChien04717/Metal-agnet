"""Geometry checks for the Design6_TXT simulation drawing and Dolan fab mask."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
from collections import Counter

import gdstk
import numpy as np
from shapely import get_parts


def check_masks(sim_file, fab_file, junctions, unnotched, fab_design):
    def layers(path):
        lib = gdstk.read_gds(str(path))
        assert lib.unit == 1e-6 and lib.precision == 1e-9
        top, = lib.top_level()
        assert np.allclose(top.bounding_box(), [(-5000, -5000), (5000, 5000)])
        return top.get_polygons()

    sim = layers(sim_file)
    fab = layers(fab_file)
    sim_metal = [p for p in sim if p.layer == 1]
    symbols = [p for p in sim if p.layer == 50]
    assert len(symbols) == len(junctions) == 6
    assert {(p.layer, p.datatype) for p in fab} == {(1, 0)}

    def contains(polygons, xy):
        return any(gdstk.inside([xy], polygons))

    def shapes(design, name):
        table = design.qgeometry.tables['poly']
        return table.loc[table.component.eq(design.components[name].id)].set_index('name').geometry

    cuts = Counter()
    for name in ('Q1', 'TC', 'Q2'):
        plain = shapes(unnotched, name)
        notched = shapes(fab_design, name)
        island = 'tmon' if name == 'TC' else 'cross'
        for electrode in (island, 'squid_ground'):
            for cut in get_parts(plain[electrode].difference(notched[electrode])):
                assert abs(cut.area * 1e6 - 8) < .01, (name, electrode, cut.area)
                xy = np.asarray(cut.representative_point().coords[0]) * 1000
                assert not contains(fab, xy), (name, electrode, xy)
                cuts[name] += 1
    assert cuts == {'Q1': 4, 'TC': 4, 'Q2': 4}
    for j in junctions:
        xy = np.asarray(j.geometry.interpolate(.5, normalized=True).coords[0]) * 1000
        assert contains(symbols, xy)
        assert not contains(sim_metal, xy) and not contains(fab, xy)
    for name in ('Q1', 'TC', 'Q2'):
        geometry = shapes(fab_design, name)
        gaps = geometry[geometry.index.str.startswith('flux_bias_line_gap')]
        assert len(gaps) == 2, (name, 'unexpected local gap geometry')
        for gap in gaps:
            xy = np.asarray(gap.representative_point().coords[0]) * 1000
            assert not contains(fab, xy) and not contains(sim_metal, xy), (name, 'local flux-line gap', xy)
    print('Verified six simulation annotations, twelve fab hooks and three local flux-line gaps.')
