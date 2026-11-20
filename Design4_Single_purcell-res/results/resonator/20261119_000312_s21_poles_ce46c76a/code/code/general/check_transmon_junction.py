"""Check TransmonJunction parameters, shared gaps, mirroring and source retention.

Run with the repository qmetal Python. No GDS or solver run is produced.
"""
from pathlib import Path
import sys

import numpy as np
from qiskit_metal import designs, draw

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'pattern design' / g) for g in ('', 'ebeam_junction')]
from transmon_junction import TransmonJunction

UM = 0.001  # mm per um in the default design units

# (metal, recess) pairs whose distance must equal junction_gap.
GAP_PAIRS = (('pad', 'landing'), ('finger', 'landing'), ('finger', 'below_bar'),
             ('finger', 'bridge'), ('bar', 'below_bar'), ('bar', 'bridge'),
             ('bar', 'above_bar'))


def build(**options):
    d = designs.DesignPlanar()
    c = TransmonJunction(d, 'tj', options=options)
    assert c.status == 'good', c.status
    rows = d.qgeometry.tables['poly']
    return c, {(int(r.layer), r['name']): r.geometry for _, r in rows.iterrows()}


def main():
    ref = TransmonJunction.reference()['polygons_um']
    # Defaults reproduce every red source polygon; recesses within 2 nm
    # (the source has 0.074-0.076um clearances, uniformized to junction_gap).
    c, geom = build()
    red = [g for (layer, _), g in geom.items() if layer == 1]
    purple = [g for (layer, _), g in geom.items() if layer == 2]
    assert len(red) == 8 and len(purple) == 10
    for coords in ref['1']:
        src = draw.scale(draw.Polygon(coords), UM, UM, origin=(0, 0))
        assert min(src.hausdorff_distance(g) for g in red) < 1e-12
    # The pad landing source clearance is 0.050um; like FluxoniumJunction it is
    # set by junction_gap, so the two landing recesses differ by 0.025um.
    deviations = sorted(min(draw.scale(draw.Polygon(coords), UM, UM, origin=(0, 0))
                            .hausdorff_distance(g) for g in purple) / UM for coords in ref['2'])
    assert all(v < 0.002 + 1e-9 for v in deviations[:-2]), deviations
    assert all(abs(v - 0.025) < 1e-9 for v in deviations[-2:]), deviations

    for gap in ('0.05um', '0.075um', '0.1um'):
        for jw in ('0.15um', '0.25um', '0.4um'):
            for sj in ('0.06um', '0.1um', '0.2um'):
                c, geom = build(junction_gap=gap, junction_width=jw, finger_width=sj)
                for side, sx in (('L', 1), ('R', -1)):
                    for metal, recess in GAP_PAIRS:
                        measured = geom[1, f'{metal}_{side}'].distance(geom[2, f'{recess}_{side}'])
                        assert abs(measured-c.p.junction_gap) < 1e-12, (gap, jw, sj, side, metal, recess, measured)
                    x0, y0, x1, y1 = geom[1, f'finger_{side}'].bounds
                    assert abs((x1-x0)-c.p.finger_width) < 1e-12
                    assert abs((x1 if sx > 0 else -x0) - (-4.325*UM)) < 1e-12
                    x0, y0, x1, y1 = geom[1, f'bar_{side}'].bounds
                    assert abs((y1-y0)-c.p.junction_width) < 1e-12 and abs(y0+0.054*UM) < 1e-12
                    # Bar stays joined to its lower lead.
                    assert geom[1, f'bar_{side}'].intersects(geom[1, f'lead_{side}'])
                # Mirror symmetry of the SQUID pair.
                for (layer, name), g in geom.items():
                    if name.endswith('_L'):
                        mirror = draw.scale(g, -1, 1, origin=(0, 0))
                        assert mirror.hausdorff_distance(geom[layer, name[:-1]+'R']) < 1e-12

    c, geom = build(squid=False)
    assert len(geom) == 9 and set(c.pins) == {'top_L', 'bottom_L'}
    c, geom = build()
    assert set(c.pins) == {'top_L', 'bottom_L', 'top_R', 'bottom_R'}
    assert np.allclose(c.pins.top_L.normal, [0, 1]) and np.allclose(c.pins.bottom_R.normal, [0, -1])
    assert np.allclose(c.pins.top_L.middle, [-4.6665*UM, 3*UM])
    d = designs.DesignPlanar()
    shifted = TransmonJunction(d, 'shifted', options=dict(orientation=90, pos_x='100um', pos_y='200um'))
    assert np.allclose(shifted.pins.top_L.normal, [-1, 0])
    assert np.allclose(shifted.pins.top_L.middle, [0.1-3*UM, 0.2-4.6665*UM])
    off = TransmonJunction(d, 'off', options=dict(jj_metal=False))
    assert len(d.qgeometry.tables['poly'].query('component == @off.id')) == 0
    for option in (dict(junction_gap='0um'), dict(junction_width='-1um'), dict(finger_width='0um'),
                   dict(junction_gap='0.2um'), dict(finger_width='0.8um'),
                   dict(junction_width='1.3um'), dict(squid='maybe')):
        invalid = TransmonJunction(d, 'invalid', options=option, make=False)
        try:
            invalid.make()
        except ValueError:
            pass
        else:
            raise AssertionError(f'Accepted invalid input: {option}')
        d.delete_component('invalid')
    print('PASS: source red geometry; recesses within 2 nm of source except landing (+25 nm); 27 gap/width/finger '
          'combinations with 7 equal gaps per junction; bar-lead contact; SQUID mirror; '
          'single junction; pins; jj_metal=False; invalid inputs')


if __name__ == '__main__':
    main()
