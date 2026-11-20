"""Design4_Single_purcell-res: test and validation chips for the shared Purcell + readout module (pattern design/resonator/purcell_readout.py).

Three builds share one channel (feed line - n0 tee - C_p EnclosingCap - PurcellReadout):
  build('open')       S21 / eigenmode fixture: resonator open end = OpenToGround (no qubit), two 50 Ohm feed ports
  build('fluxonium')  MyFluxonium readout pad abutted to the module's 'qubit' pin
  build('transmon')   MyTmonCoupler (readout_line_options.make_rol) abutted to the 'qubit' pin

Geometry inputs: parameters/geometry.json (loaded here; lengths per build come from the frequency model or from an
explicit length file). Run from the repository root in the qmetal env:

    python -I -B Design4_Single_purcell-res/layout/build.py            # three GDS (sim geometry, no cheese) + figures + checks
"""
from pathlib import Path
import copy
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'code' / 'general')]
import project_paths  # noqa: E402
project_paths.configure(__file__)

import numpy as np  # noqa: E402
from qiskit_metal import designs, Dict  # noqa: E402
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight  # noqa: E402
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround  # noqa: E402

from enclosing_cap import EnclosingCap  # noqa: E402
from tap_junction import TapJunction  # noqa: E402
from res_end import ResEnd  # noqa: E402
from purcell_readout import PurcellReadout  # noqa: E402
import purcell_readout_model as M  # noqa: E402

DESIGN = ROOT / 'Design4_Single_purcell-res'
GEOMETRY = json.loads((DESIGN / 'parameters' / 'geometry.json').read_text(encoding='utf-8'))
PORTS = [('feed_in', 'start', 50), ('feed_out', 'end', 50)]


def um(v):
    return f'{v:.6f}um'


def mm(v):
    return f'{v:.9f}mm'


def model_lengths(kind, caps=None):
    """Filter / resonator lengths (mm) for the bare f_target of geometry.json, from the frequency model.
    caps: Q3D values dict(cp_mutual, cp_plate_ground, cj, head_excess) overriding the estimates in geometry.json."""
    g = GEOMETRY
    cpw = g['cpw']['model_preset']
    est = dict(g['model_loads'])
    if caps:
        est.update(caps)
    h = g['hammer']
    stub = g['module']['line_spacing_um'] / 2 / 1000                  # tap -> hammer centre (deembedded plane)
    c_stub = est['head_excess'] + est['cj']
    c_end_f = est['cp_mutual'] + est['cp_plate_ground']
    c_end_r = dict(open=M.open_end_c(cpw), fluxonium=est['qubit_pad_fluxonium'],
                   transmon=est['qubit_pad_transmon'])[kind]
    f = g['f_target_GHz']
    lf = M.length_for_f(f, g['module']['filter_tap_mm'], c_end_f, stub, c_stub, cpw)
    lr = M.length_for_f(f, g['module']['resonator_tap_mm'], c_end_r, stub, c_stub, cpw)
    return dict(filter_mm=lf, resonator_mm=lr, c_end_filter_fF=c_end_f, c_end_resonator_fF=c_end_r,
                c_stub_fF=c_stub, stub_mm=stub, cpw=cpw, f_target_GHz=f, loads=est, hammer=h)


def module_options(lengths):
    g = GEOMETRY
    m, h, c = g['module'], g['hammer'], g['cpw']
    return dict(filter_length=mm(lengths['filter_mm']), resonator_length=mm(lengths['resonator_mm']),
                filter_tap=mm(m['filter_tap_mm']), resonator_tap=mm(m['resonator_tap_mm']),
                line_spacing=um(m['line_spacing_um']), tee_half=um(m['tee_half_um']), lead=um(m['lead_um']),
                jog=um(m['jog_um']), meander_pitch=um(m['meander_pitch_um']),
                filter_meander_height=um(m['filter_meander_height_um']),
                resonator_meander_height=um(m['resonator_meander_height_um']), fillet=um(m['fillet_um']),
                trace_width=um(c['width_um']), trace_gap=um(c['gap_um']),
                head_width=um(h['head_width_um']), head_height=um(h['head_height_um']), head_gap=um(h['head_gap_um']),
                hammer_ground_gap=um(h['ground_gap_um']), corner_radius=um(h['corner_radius_um']))


def _place_qubit_at_pin(design, qubit, pin_name, target):
    """Translate an already built qubit so its readout pin middle sits on `target` (x, y in mm)."""
    middle = np.asarray(qubit.pins[pin_name]['middle'], float)
    x = design.parse_value(qubit.options.pos_x) + target[0] - middle[0]
    y = design.parse_value(qubit.options.pos_y) + target[1] - middle[1]
    qubit.options.pos_x, qubit.options.pos_y = mm(x), mm(y)
    qubit.rebuild()
    assert np.allclose(qubit.pins[pin_name]['middle'], target, atol=1e-9)


def build(kind='open', lengths=None, jj_metal=False, wire_bonds=True):
    """Return (design, info). kind: 'open' | 'fluxonium' | 'transmon'. lengths: dict(filter_mm, resonator_mm) or None
    (model lengths for this end load)."""
    g = GEOMETRY
    if lengths is None:
        lengths = model_lengths(kind)
    c = g['cpw']
    cpw = dict(trace_width=um(c['width_um']), trace_gap=um(c['gap_um']))
    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = g['substrate']['material']
    chip = g['chip']
    design.chips.main.size.update(size_x=mm(chip['size_x_mm']), size_y=mm(chip['size_y_mm']), size_z='-500um',
                                  center_x='0mm', center_y=mm(chip['center_y_mm']))
    design.variables.update(cpw_width=cpw['trace_width'], cpw_gap=cpw['trace_gap'])
    m = g['module']
    s = m['line_spacing_um'] / 1000
    xf = -s / 2
    feed_y = g['feed']['y_mm']
    enc = {k: v for k, v in g['cp_cap'].items() if k != 'source'}
    enc_half = (enc['plate_length'] / 2 + enc['leg_ext'] + enc['lead_length']) / 1000      # centre -> pin
    # feed line: pin-only ends (ports at the route ends), tee n0 above the filter, stub down to C_p
    for name, x, ori in (('feed_l', -g['feed']['half_length_mm'], '0'), ('feed_r', g['feed']['half_length_mm'], '180')):
        ResEnd(design, name, options=dict(pos_x=mm(x), pos_y=mm(feed_y), orientation=ori, **cpw))
    TapJunction(design, 'n0', options=dict(pos_x=mm(xf), pos_y=mm(feed_y), stub_up='False',
                                           stub_length=um(g['feed']['stub_um']), **cpw))
    cp_y = feed_y - g['feed']['stub_um'] / 1000 - g['feed']['cp_lead_mm'] - enc_half
    EnclosingCap(design, 'cp', options=dict({k: um(v) for k, v in enc.items()}, **cpw,
                                            pos_x=mm(xf), pos_y=mm(cp_y), orientation='270'))
    for name, a, b in (('feed_in', ('feed_l', 'end'), ('n0', 'prime_start')),
                       ('feed_out', ('n0', 'prime_end'), ('feed_r', 'end')),
                       ('stub_p', ('n0', 'second_end'), ('cp', 'a'))):
        RouteStraight(design, name, Dict(pin_inputs=Dict(start_pin=Dict(component=a[0], pin=a[1]),
                                                         end_pin=Dict(component=b[0], pin=b[1])), **cpw))
    # module: its 'cp' pin abuts C_p pin b
    pb = np.asarray(design.components['cp'].pins['b']['middle'], float)
    opts = module_options(lengths)
    ft = m['filter_tap_mm']
    purcell_resonator_pair = PurcellReadout(design, 'purcell_resonator_pair', options=dict(opts, pos_x=mm(pb[0] - xf), pos_y=mm(pb[1] - ft)))
    qpin = np.asarray(purcell_resonator_pair.pins['qubit']['middle'], float)
    qubit = None
    if kind == 'open':
        OpenToGround(design, 'readout_open', options=dict(pos_x=mm(qpin[0]), pos_y=mm(qpin[1]), orientation='270',
                                                          width=cpw['trace_width'], gap=cpw['trace_gap']))
    elif kind == 'fluxonium':
        from my_fluxonium import MyFluxonium
        fx = copy.deepcopy(g['qubits']['fluxonium'])
        qubit = MyFluxonium(design, 'qubit', options=dict(fx, jj_metal=str(jj_metal)))
        _place_qubit_at_pin(design, qubit, 'readout_line', qpin)
    elif kind == 'transmon':
        from my_tmon import MyTmonCoupler
        tm = copy.deepcopy(g['qubits']['transmon'])
        tm['readout_line_options'] = Dict(tm['readout_line_options'], make_rol=True)
        qubit = MyTmonCoupler(design, 'qubit', options=dict(tm, jj_metal=str(jj_metal)))
        _place_qubit_at_pin(design, qubit, 'readout_line', qpin)
    else:
        raise ValueError(kind)
    if wire_bonds:
        t = design.qgeometry.tables['path']
        t['hfss_wire_bonds'] = t['subtract'].astype(bool)
    info = dict(kind=kind, lengths=lengths, length_report=purcell_resonator_pair.length_report, feed_y_mm=feed_y)
    return design, info


# ---- checks ---------------------------------------------------------------------------------------------------------
def pin_check(design):
    """Every connected pin pair coincides with opposite normals; the module pins are met by a component."""
    out = []
    purcell_resonator_pair = design.components['purcell_resonator_pair']
    partners = dict(cp=('cp', 'b'), qubit=None)
    qubit = design.components['qubit'] if 'qubit' in design.components else None
    if qubit is not None:
        partners['qubit'] = ('qubit', 'readout_line')
    elif 'readout_open' in design.components:
        partners['qubit'] = ('readout_open', 'open')
    for pin, (comp, other) in partners.items():
        a, b = purcell_resonator_pair.pins[pin], design.components[comp].pins[other]
        d = float(np.linalg.norm(np.asarray(a['middle']) - np.asarray(b['middle']))) * 1000
        n = float(np.dot(a['normal'], b['normal']))
        wa, wb = float(a['width']) * 1000, float(b['width']) * 1000
        out.append(dict(module_pin=pin, partner=f'{comp}.{other}', offset_um=d, normal_dot=n,
                        width_um=(wa, wb), ok=bool(d < 1e-3 and n < -0.999 and abs(wa - wb) < 1e-6)))
    return out


def spacing_check(design, min_um=400.0):
    """Side-by-side CPW segment pairs of different lines closer than min_um (skill.md section 11). Exempt: segments of
    the same component, segments shorter than 2 x fillet (bends), perpendicular pairs."""
    t = design.qgeometry.tables['path']
    rows = t[~t['subtract'].astype(bool)]
    segs = []
    for r in rows.itertuples():
        pts = np.asarray(r.geometry.coords)
        for a, b in zip(pts[:-1], pts[1:]):
            segs.append((int(r.component), r.name, a, b))
    bad = []
    for i, (ca, na, a0, a1) in enumerate(segs):
        ua = (a1 - a0) / np.linalg.norm(a1 - a0)
        for cb, nb, b0, b1 in segs[i + 1:]:
            if ca == cb and na == nb:
                continue
            ub = (b1 - b0) / np.linalg.norm(b1 - b0)
            if abs(abs(ua @ ub) - 1) > 1e-6:
                continue
            nrm = np.array([-ua[1], ua[0]])
            dist = abs((b0 - a0) @ nrm) * 1000
            ta, tb = sorted([0, np.linalg.norm(a1 - a0)]), sorted([(b0 - a0) @ ua, (b1 - a0) @ ua])
            overlap = (min(ta[1], tb[1]) - max(ta[0], tb[0])) * 1000
            if dist < 1e-6 or overlap <= 1.0 or dist >= min_um - 1e-6:
                continue
            names = {design._components[ca].name, design._components[cb].name}
            bad.append(dict(a=f'{design._components[ca].name}.{na}', b=f'{design._components[cb].name}.{nb}',
                            spacing_um=round(dist, 2), overlap_um=round(overlap, 2), components=sorted(names)))
    return bad


def export(kind, out_gds, out_png_prefix, lengths=None):
    """Junction-free simulation-geometry GDS (no cheese), layout figures and the 4 um metal check."""
    import gds_tools
    import check_metal_width as cmw
    design, info = build(kind, lengths=lengths, jj_metal=False)
    exempt = cmw.design_exemptions(design)
    pins = pin_check(design)
    spacing = spacing_check(design)
    gds_tools.export_gds(design, str(out_gds), fab=False)
    found = cmw.check(out_gds, exempt=exempt)
    chip = GEOMETRY['chip']
    cy, hx, hy = chip['center_y_mm'], chip['size_x_mm'] / 2, chip['size_y_mm'] / 2
    gds_tools.plot_gds(str(out_gds), f'{out_png_prefix}_layout.png', window=(-hx, hx, cy - hy, cy + hy),
                       figsize=(10, 10 * hy / hx))
    # zooms: hammer / T region and the qubit coupling
    hx0 = design.parse_value(design.components['purcell_resonator_pair'].options.pos_x)
    hy0 = design.parse_value(design.components['purcell_resonator_pair'].options.pos_y)
    gds_tools.plot_gds(str(out_gds), f'{out_png_prefix}_hammer.png', window=(hx0 - .35, hx0 + .35, hy0 - .25, hy0 + .25),
                       figsize=(9, 6.5), dpi=150)
    q = np.asarray(design.components['purcell_resonator_pair'].pins['qubit']['middle'])
    gds_tools.plot_gds(str(out_gds), f'{out_png_prefix}_qubit.png', window=(q[0] - .45, q[0] + .45, q[1] - .7, q[1] + .2),
                       figsize=(8, 8), dpi=150)
    return dict(info, pins=pins, spacing_violations=spacing,
                metal_width=dict(regions_below_4um=[dict(width_um=w, bounds=g.bounds) for g, w in found],
                                 exempt='design_exemptions: pattern design/resonator classes incl. PurcellReadout'),
                gds=str(Path(out_gds).relative_to(ROOT)), figures=[f'{out_png_prefix}_{v}.png' for v in ('layout', 'hammer', 'qubit')])


def main():
    import datetime
    rev = GEOMETRY['revision']
    out = {}
    for kind in ('open', 'fluxonium', 'transmon'):
        gds = DESIGN / 'gds' / f'Design8_{rev}_{kind}_sim.gds'
        png = DESIGN / 'fig' / f'Design8_{rev}_{kind}'
        r = export(kind, gds, png)
        r['figures'] = [str(Path(f).relative_to(ROOT)) for f in r['figures']]
        out[kind] = r
        print(kind, 'lengths', {k: round(v, 6) for k, v in r['length_report'].items() if 'mm' in k},
              '| pins ok', all(p['ok'] for p in r['pins']), '| spacing violations', len(r['spacing_violations']),
              '| <4um regions', len(r['metal_width']['regions_below_4um']))
    out['created'] = datetime.datetime.now().isoformat(timespec='seconds')
    (DESIGN / 'parameters' / f'build_checks_{rev}.json').write_text(json.dumps(out, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
