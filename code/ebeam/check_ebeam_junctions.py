"""Self-contained checks for the e-beam junction components and the EJ model.

Measures the drawn geometry of ``FluxoniumJunctionDolan`` and ``TransmonJunctionDolan``
(pattern design/ebeam_junction) and checks it against the rules in
skill/ebeam_junction.md:

* drawn finger / cross-electrode widths equal the options (also via EJ ->
  options from code/qubit/jj_resistance_model.py);
* every small JJ has the finger on the left of the cross electrode, with the
  set bridge gap, the finger spanning the electrode, and fixed finger right /
  electrode top edges when widths change;
* shadow openings keep ``shadow_clearance`` from the electrode; no small-JJ
  main / shadow overlap; small-JJ and array main doses only touch;
* SQUID branch 2 is branch 1 translated by ``squid_spacing`` (no mirror);
* pins, rotation, ``num_array``, ``jj_metal=False`` and invalid inputs.

No external reference data and no solver. Run with the qmetal Python:

    python -I -X utf8 code/general/check_ebeam_junctions.py
"""
from pathlib import Path
import sys

from qiskit_metal import designs
from shapely import union_all
from shapely.affinity import scale, translate

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / 'pattern design'),
    str(ROOT / 'pattern design' / 'ebeam_junction'),
    str(ROOT / 'pattern design' / 'transmon'),
    str(ROOT / 'code' / 'qubit'),
    str(ROOT / 'code' / 'ebeam'),
]
from fluxonium_junction_dolan import FluxoniumJunctionDolan  # noqa: E402
from transmon_junction_dolan import TransmonJunctionDolan  # noqa: E402
from transmon_junction_manhattan import TransmonJunctionManhattan  # noqa: E402
from fixed_transmon_pocket_manhattan import FixedTransmonPocketManhattan  # noqa: E402
from tunable_transmon_pocket_manhattan import TunableTransmonPocketManhattan  # noqa: E402
from xmon_manhattan import XmonManhattan  # noqa: E402
from tmon_coupler_manhattan import TmonCouplerManhattan  # noqa: E402
from claw_coupler_manhattan import ClawCouplerManhattan  # noqa: E402
from bracket_tmon_manhattan import BracketTmonCouplerManhattan  # noqa: E402
from ebeam_junction_placement import add_ebeam_junctions  # noqa: E402
import jj_resistance_model as jjm  # noqa: E402

TOL = 1e-6   # um


def build(cls, **options):
    d = designs.DesignPlanar()
    c = cls(d, 'c', options=options)
    assert c.status == 'good', c.status
    rows = d.qgeometry.tables['poly']
    # Design units are mm; return {name: (layer, geometry in um)}.
    return c, {r['name']: (int(r.layer), scale(r.geometry, 1000, 1000, origin=(0, 0)))
               for _, r in rows.iterrows()}


def width(g, axis):
    x0, y0, x1, y1 = g.bounds
    return (x1 - x0) if axis == 'x' else (y1 - y0)


def check_junction(shapes, finger, electrode, fw, cw, gap, tag):
    """One small JJ: widths, finger left of electrode, bridge gap, span, clearances."""
    f, e = shapes[finger][1], shapes[electrode][1]
    assert abs(width(f, 'x') - fw) < TOL, (tag, 'finger width', width(f, 'x'))
    assert abs(width(e, 'y') - cw) < TOL, (tag, 'electrode width', width(e, 'y'))
    assert f.bounds[2] < e.bounds[0], (tag, 'finger must be left of the electrode')
    assert abs(e.bounds[0] - f.bounds[2] - gap) < TOL, (tag, 'bridge gap')
    assert f.bounds[1] <= e.bounds[1] + TOL and f.bounds[3] >= e.bounds[3] - TOL, (tag, 'finger spans electrode')
    return f.bounds[2], e.bounds[3]          # fixed edges: finger right, electrode top


def check_layers(shapes, tag):
    def u(layers):
        return union_all([g for l, g in shapes.values() if l in layers])
    m, s, a = u((1,)), u((2,)), u((3,))
    assert m.intersection(s).area < 1e-9, (tag, 'small-JJ main overlaps shadow')
    assert m.intersection(a).area < 1e-9, (tag, 'small-JJ main overlaps array main')


def electrode_clearance(shapes, electrode, sc, tag):
    e = shapes[electrode][1]
    shadow = union_all([g for l, g in shapes.values() if l == 2])
    gap = e.distance(shadow)
    assert abs(gap - sc) < 1e-4, (tag, 'shadow clearance', gap)


def fluxonium():
    fixed = None
    for side in ('left', 'right'):
        for fw, cw in ((80, 92), (80, 300), (100, 200), (120, 600)):
            c, s = build(FluxoniumJunctionDolan, side=side, jj_finger_width=f'{fw}nm',
                         jj_cross_electrode_width=f'{cw}nm')
            tag = f'FX {side} {fw}/{cw}'
            edges = check_junction(s, 'jj_finger', 'jj_cross_electrode', fw / 1000, cw / 1000, 0.248, tag)
            if fixed is None or fixed[0] != side:
                fixed = (side, edges)
            assert max(abs(a - b) for a, b in zip(edges, fixed[1])) < TOL, (tag, 'fixed edges moved')
            electrode_clearance(s, 'jj_cross_electrode', 0.075, tag)
            check_layers(s, tag)
            # Both contacts and the small JJ meet the array only at shared edges.
            arr = union_all([g for l, g in s.values() if l == 3])
            jj = union_all([g for l, g in s.values() if l == 1])
            assert jj.intersects(arr) and jj.intersection(arr).area < 1e-9, (tag, 'JJ-array contact')
    # Shape names follow the real position after the 180 deg turn.
    _, s = build(FluxoniumJunctionDolan, side='right')
    assert s['contact_top'][1].centroid.y > 0 and s['top_array_0000'][1].centroid.y > 0
    assert s['top_array_0000'][1].centroid.x > 0, 'right device array runs towards +x'
    # Rounded-bay variant.
    _, s = build(FluxoniumJunctionDolan, array_contact_body_width='7.6um', array_contact_shadow_width='0.958um')
    assert abs(width(s['contact_top'][1], 'x') - 7.6) < TOL
    # num_array: cell count and fixed junction.
    for n in (1, 20, 80):
        _, s = build(FluxoniumJunctionDolan, num_array=n)
        assert sum('_array_' in k for k in s) == 2 * n


def transmon():
    fixed = None
    for jj1, jj2 in ((271, 271), (384, 157), (594, 199), (100, 600)):
        c, s = build(TransmonJunctionDolan, jj1_cross_electrode_width=f'{jj1}nm', jj2_cross_electrode_width=f'{jj2}nm')
        for tag, cw in (('jj1', jj1), ('jj2', jj2)):
            edges = check_junction(s, f'finger_{tag}', f'cross_electrode_{tag}', 0.080, cw / 1000, 0.249,
                                   f'TC {jj1}/{jj2} {tag}')
            electrode_clearance({k: v for k, v in s.items() if k.endswith(tag)}, f'cross_electrode_{tag}', 0.075, tag)
            if tag == 'jj1':
                fixed = fixed or edges
                assert max(abs(a - b) for a, b in zip(edges, fixed)) < TOL, 'fixed edges moved'
        check_layers(s, f'TC {jj1}/{jj2}')
    # Translation-only SQUID (equal widths): every jj2 shape = jj1 shape + squid_spacing.
    _, s = build(TransmonJunctionDolan, squid_spacing='30um')
    for name, (layer, g) in s.items():
        if name.endswith('_jj1'):
            other = s[name[:-1] + '2'][1]
            assert translate(g, 30, 0).symmetric_difference(other).area < 1e-9, name
    # Design3_FTF_NTHU-style hooks (+2 um) and test-key body offset.
    build(TransmonJunctionDolan, lower_contact_offset='2um', lower_contact_lead_offset='1um')
    _, s = build(TransmonJunctionDolan, squid=False, upper_contact_offset='3um')
    assert s['upper_contact_jj1'][1].intersects(s['finger_jj1'][1]), 'test-key ledge joins the finger'


def model_to_component():
    jjm.self_check()
    _, s = build(FluxoniumJunctionDolan, **jjm.size_from_ej(4.0).options())
    assert abs(width(s['jj_cross_electrode'][1], 'y') - 0.092) < TOL
    sizes = jjm.squid_sizes(23.056, jjm.asymmetry_from_flux_points(23.056, 9.746))
    _, s = build(TransmonJunctionDolan, **sizes['options'])
    a1 = (0.080 + 0.030) * (width(s['cross_electrode_jj1'][1], 'y') + 0.070)
    a2 = (0.080 + 0.030) * (width(s['cross_electrode_jj2'][1], 'y') + 0.070)
    assert abs(a1 / a2 - sizes['area_ratio']) / sizes['area_ratio'] < 0.01, (a1 / a2, sizes['area_ratio'])


def misc():
    d = designs.DesignPlanar()
    f = FluxoniumJunctionDolan(d, 'f', options=dict(orientation=90, pos_x='100um', pos_y='200um'))
    top = f.pins['top'].middle * 1000
    assert abs(top[0] - 80) < 1e-9 and abs(top[1] - 200) < 1e-9, top
    t = TransmonJunctionDolan(d, 't', options=dict(squid=False))
    assert set(t.pins) == {'top_jj1', 'bottom_jj1'}
    d = designs.DesignPlanar()
    FluxoniumJunctionDolan(d, 'off', options=dict(jj_metal=False))
    TransmonJunctionDolan(d, 'off2', options=dict(jj_metal=False))
    assert len(d.qgeometry.tables['poly']) == 0
    for cls, bad in ((FluxoniumJunctionDolan, dict(jj_cross_electrode_width='0.7um')),
                     (FluxoniumJunctionDolan, dict(side='middle')),
                     (FluxoniumJunctionDolan, dict(jj_bridge_gap='0.1um')),
                     (FluxoniumJunctionDolan, dict(array_contact_body_width='1um')),
                     (TransmonJunctionDolan, dict(jj1_cross_electrode_width='0.7um')),
                     (TransmonJunctionDolan, dict(squid_spacing='3um')),
                     (TransmonJunctionDolan, dict(upper_contact_offset='-3um')),
                     (TransmonJunctionDolan, dict(lower_contact_offset='2um'))):
        d = designs.DesignPlanar()
        c = cls(d, 'bad', options=bad, make=False)
        try:
            c.make()
        except ValueError:
            continue
        raise AssertionError(f'{cls.__name__} accepted {bad}')


def manhattan():
    # Single Manhattan JJ geometry & vertex count (20-vertex merged polygon, 4 undercut boxes)
    c, s = build(TransmonJunctionManhattan, squid='False',
                 jj1_finger_x_width='190nm', jj1_finger_y_width='170nm')
    assert set(c.pins) == {'top_jj1', 'bottom_jj1'}
    jj1 = s['jj_main_jj1'][1]
    assert jj1.geom_type == 'Polygon' and len(jj1.exterior.coords) - 1 == 20, len(jj1.exterior.coords) - 1
    assert abs(s['undercut_pad_top_jj1'][1].area - 45.0) < 1e-4
    assert abs(s['undercut_pad_bot_jj1'][1].area - 45.0) < 1e-4
    # Tip undercut boxes: 1.0 x (0.19 + 0.25) = 0.44 and 1.0 x (0.17 + 0.25) = 0.42
    assert abs(s['undercut_tip_h_jj1'][1].area - 0.44) < 1e-4
    assert abs(s['undercut_tip_v_jj1'][1].area - 0.42) < 1e-4

    # Translation-only SQUID (never mirrored)
    _, sq = build(TransmonJunctionManhattan, squid='True', squid_spacing='30um')
    for name, (layer, g) in sq.items():
        if name.endswith('_jj1'):
            other = sq[name[:-1] + '2'][1]
            assert translate(g, 30.0, 0.0).symmetric_difference(other).area < 1e-6, name

    # Asymmetric SQUID widths & jj_metal=False
    build(TransmonJunctionManhattan, squid='True',
          jj1_finger_x_width='220nm', jj1_finger_y_width='200nm',
          jj2_finger_x_width='120nm', jj2_finger_y_width='100nm')
    d_off = designs.DesignPlanar()
    TransmonJunctionManhattan(d_off, 'off', options=dict(jj_metal='False'))
    assert len(d_off.qgeometry.tables['poly']) == 0

    # Invalid inputs rejected
    for bad in (dict(squid_spacing='5um'), dict(jj1_finger_x_width='0nm'),
                dict(tip_undercut_overlap='1.5um'), dict(cross_x='10um')):
        d = designs.DesignPlanar()
        comp = TransmonJunctionManhattan(d, 'bad', options=bad, make=False)
        try:
            comp.make()
        except ValueError:
            continue
        raise AssertionError(f'TransmonJunctionManhattan accepted {bad}')

    # All 6 Manhattan transmon components + e-beam placement overlap check
    cases = [
        (FixedTransmonPocketManhattan, 'FTM', {}),
        (TunableTransmonPocketManhattan, 'TTM', {}),
        (XmonManhattan, 'XMS', dict(squid=dict(make='True'))),
        (XmonManhattan, 'XMF', dict(squid=dict(make='False'))),
        (TmonCouplerManhattan, 'TCM', {}),
        (ClawCouplerManhattan, 'CLM', {}),
        (BracketTmonCouplerManhattan, 'BKM', dict(pad_arm_length='400um')),
    ]
    for cls, name, opts in cases:
        d = designs.DesignPlanar()
        q = cls(d, name, options=opts)
        assert q.status == 'good', (name, q.status)
        rec = add_ebeam_junctions(d, manhattan={name: dict(jj1_finger_x_width='190nm', jj1_finger_y_width='170nm')})
        assert name + '_EBJ' in rec['placements']
        # Verify each e-beam Manhattan JJ pad overlaps the base-metal 45-deg landing arms cleanly
        polys = d.qgeometry.tables['poly']
        base = union_all([r.geometry for _, r in polys.iterrows() if not r['subtract'] and int(r['layer']) == 1])
        ebj = union_all([r.geometry for _, r in polys.iterrows() if int(r['layer']) == 11])
        overlap_um2 = base.intersection(ebj).area * 1e6
        expected_sites = 2 if q.manhattan_sites()['squid'] else 1
        # Each site has 2 pads crossing 3 um x 4 um = 12 um^2 each -> ~24 um^2 per site
        assert overlap_um2 > 20.0 * expected_sites, (name, 'pad-to-lead overlap too small', overlap_um2)


def main():
    fluxonium()
    transmon()
    manhattan()
    model_to_component()
    misc()
    print('PASS: drawn widths, finger-left orientation, bridge gaps, fixed edges, 75 nm shadow clearance, '
          'no main/shadow overlap, translated SQUID, Manhattan JJ + 6 transmons + placement, '
          'EJ -> options, pins, num_array, jj_metal, invalid inputs')


if __name__ == '__main__':
    main()
