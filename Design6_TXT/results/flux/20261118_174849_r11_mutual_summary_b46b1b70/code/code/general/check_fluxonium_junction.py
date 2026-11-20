"""Check fixed-cell repetition, common gaps, placement and source retention.

Run with the repository qmetal Python. Optional --preview-dir writes draft
Qiskit Metal GUI screenshots; no fabrication GDS or solver run is produced.
"""
import argparse
from pathlib import Path
import sys

import numpy as np
from qiskit_metal import designs, draw

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'pattern design'))
from fluxonium_junction import FluxoniumJunction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preview-dir', type=Path)
    args = parser.parse_args()
    for count in (1, 12, 43, 80):
        for gap in ('0.05um', '0.075um', '0.1um'):
            d = designs.DesignPlanar()
            c = FluxoniumJunction(d, 'sample', options=dict(num_cells=count, junction_gap=gap))
            assert c.status == 'good', c.status
            rows = d.qgeometry.tables['poly']
            assert all(rows.geometry.is_valid) and all(rows.geometry.area > 0)
            cells = rows[rows.name.str.contains('_cell_')]
            assert len(cells) == 2*count
            for shape in cells.geometry:
                x0, y0, x1, y1 = shape.bounds
                assert np.allclose([x1-x0, y1-y0], [0.000998, 0.002], atol=1e-12)
            geom = {(row.layer, row['name']): row.geometry for _, row in rows.iterrows()}
            expected = c.p.junction_gap
            # Arrows at top landing, both sides of the finger and right bridge.
            for metal, recess in ((0, 0), (1, 0), (1, 9), (2, 11), (2, 9), (2, 10)):
                measured = geom[1, f'reference_{metal}'].distance(geom[2, f'reference_{recess}'])
                assert abs(measured-expected) < 1e-12, (count, gap, metal, recess, measured)
            assert np.allclose(c.pins.top.normal, [0, 1])
            assert np.allclose(c.pins.bottom.normal, [0, -1])
    # Preserve every red source shape exactly at default placement.
    d = designs.DesignPlanar()
    c = FluxoniumJunction(d, 'default')
    rows = d.qgeometry.tables['poly']
    for i, coords in enumerate(c.reference()['polygons_um']['1']):
        expected = draw.scale(draw.Polygon(coords), 0.001, 0.001, origin=(0, 0))
        actual = rows[(rows.layer == 1) & (rows.name == f'reference_{i}')].iloc[0].geometry
        assert actual.equals_exact(expected, 1e-12)
    shifted = FluxoniumJunction(d, 'shifted', options=dict(
        orientation=90, pos_x='100um', pos_y='200um'))
    assert np.allclose(shifted.pins.top.normal, [-1, 0])
    assert np.allclose(shifted.pins.top.middle, [0.084, 0.2])
    off = FluxoniumJunction(d, 'off', options=dict(jj_metal=False))
    assert len(d.qgeometry.tables['poly'].query('component == @off.id')) == 0
    # QComponent records build failures; directly call make to verify rejection.
    for option in (dict(num_cells='2.5'), dict(num_cells=0), dict(junction_gap='0um'),
                   dict(cell_length='-1um'), dict(junction_gap='0.25um')):
        invalid = FluxoniumJunction(d, 'invalid', options=option, make=False)
        try:
            invalid.make()
        except ValueError:
            pass
        else:
            raise AssertionError(f'Accepted invalid dimensions: {option}')
        d.delete_component('invalid')
    print('PASS: 12 count/gap combinations; fixed cell dimensions; six equal gaps; '
          'source red geometry; transformed pins; jj_metal=False; invalid inputs')
    if args.preview_dir:
        from qiskit_metal import MetalGUI
        args.preview_dir.mkdir(parents=True, exist_ok=True)
        preview = designs.DesignPlanar()
        FluxoniumJunction(preview, 'fluxonium_junction')
        gui = MetalGUI(preview)
        gui.rebuild()
        gui.zoom_on_components(['fluxonium_junction'])
        gui.screenshot(name=str(args.preview_dir / 'fluxonium_junction_draft'), display=False)
        gui.canvas.zoom_to_rectangle((0.0023, 0.0014, 0.0045, 0.0041))
        gui.screenshot(name=str(args.preview_dir / 'fluxonium_junction_detail'), display=False)
        gui.main_window.force_close = True
        gui.main_window.close()


if __name__ == '__main__':
    main()
