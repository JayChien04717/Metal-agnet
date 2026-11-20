"""Design4: PurcellReadout hammer-stub test, open fixture (no qubit), two variants of the same channel.

feed line (two 50 Ohm ports) - n0 tee - C_p (Cai T-type EnclosingCap) - PurcellReadout - OpenToGround
  variant 'ref'   line_spacing 400 um: metal stub 153.5 um (Design5 / Design6 module as drawn)
  variant 'half'  line_spacing 243 um: metal stub 75 um
Everything else (16/9 um CPW, eps_r 11.5, hammer 25 x 65 um heads, 27 um gap, taps 1.2 / 1.3 mm, f_target 7 GHz)
is shared; see parameters/geometry.json. Lengths come from the frequency model with the Q3D loads of the variant's
q3d_module run (or the estimates in geometry.json for dry runs).

    python -I -B Design4/layout/build.py        # both variants: simulation GDS + figures + checks (estimate lengths)
"""
from pathlib import Path
import json
import os
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

DESIGN = ROOT / 'Design4'
GEOMETRY = json.loads((DESIGN / 'parameters' / 'geometry.json').read_text(encoding='utf-8'))
PORTS = [('feed_in', 'start', 50), ('feed_out', 'end', 50)]
os.environ['SI_EPS'] = f"{GEOMETRY['substrate']['eps_r']:g}"     # save_aedt.ensure_si_material reads it


def um(v):
    return f'{v:.6f}um'


def mm(v):
    return f'{v:.9f}mm'


def spacing_um(variant):
    return GEOMETRY['variants'][variant]['line_spacing_um']


def stub_metal_um(variant):
    g = GEOMETRY
    h, w = g['hammer'], g['cpw']['width_um']
    return spacing_um(variant) / 2 - w / 2 - (h['head_gap_um'] / 2 + h['head_width_um'])


def model_lengths(variant, caps=None):
    """Filter / resonator lengths (mm) for the bare f_target. caps: Q3D dict(cp_mutual, cp_plate_ground, cj,
    head_excess); without it the hammer loads come from the purcell_readout_model fits (dry runs only)."""
    g = GEOMETRY
    cpw = g['cpw']['model_preset']
    h = g['hammer']
    est = dict(cp_mutual=g['model_loads']['cp_mutual'], cp_plate_ground=g['model_loads']['cp_plate_ground'],
               cj=M.hammer_cj(h['head_gap_um'], h['head_height_um']),
               head_excess=M.hammer_excess_ground(h['head_gap_um'], h['head_height_um']))
    if caps:
        est.update(caps)
    stub = spacing_um(variant) / 2 / 1000                  # tap -> hammer centre (deembedded reference plane)
    c_stub = est['head_excess'] + est['cj']
    c_end_f = est['cp_mutual'] + est['cp_plate_ground']
    c_end_r = M.open_end_c(cpw)
    f = g['f_target_GHz']
    lf = M.length_for_f(f, g['module']['filter_tap_mm'], c_end_f, stub, c_stub, cpw)
    lr = M.length_for_f(f, g['module']['resonator_tap_mm'], c_end_r, stub, c_stub, cpw)
    return dict(variant=variant, filter_mm=lf, resonator_mm=lr, c_end_filter_fF=c_end_f, c_end_resonator_fF=c_end_r,
                c_stub_fF=c_stub, stub_mm=stub, cpw=cpw, f_target_GHz=f, loads=est,
                loads_source='q3d' if caps else 'estimates')


def module_options(variant, lengths):
    g = GEOMETRY
    m, h, c = g['module'], g['hammer'], g['cpw']
    return dict(filter_length=mm(lengths['filter_mm']), resonator_length=mm(lengths['resonator_mm']),
                filter_tap=mm(m['filter_tap_mm']), resonator_tap=mm(m['resonator_tap_mm']),
                line_spacing=um(spacing_um(variant)), tee_half=um(m['tee_half_um']), lead=um(m['lead_um']),
                jog=um(m['jog_um']), meander_pitch=um(m['meander_pitch_um']),
                filter_meander_height=um(m['filter_meander_height_um']),
                resonator_meander_height=um(m['resonator_meander_height_um']), fillet=um(m['fillet_um']),
                trace_width=um(c['width_um']), trace_gap=um(c['gap_um']),
                head_width=um(h['head_width_um']), head_height=um(h['head_height_um']), head_gap=um(h['head_gap_um']),
                hammer_ground_gap=um(h['ground_gap_um']), corner_radius=um(h['corner_radius_um']))


def cp_options():
    return {k: um(v) for k, v in GEOMETRY['cp_cap'].items() if k != 'source'}


def build(variant, lengths=None, wire_bonds=True):
    """Return (design, info) of the open fixture for `variant` ('ref' | 'half')."""
    g = GEOMETRY
    if lengths is None:
        lengths = model_lengths(variant)
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
    s = spacing_um(variant) / 1000
    xf = -s / 2
    feed_y = g['feed']['y_mm']
    enc = g['cp_cap']
    enc_half = (enc['plate_length'] / 2 + enc['leg_ext'] + enc['lead_length']) / 1000      # centre -> pin
    for name, x, ori in (('feed_l', -g['feed']['half_length_mm'], '0'), ('feed_r', g['feed']['half_length_mm'], '180')):
        ResEnd(design, name, options=dict(pos_x=mm(x), pos_y=mm(feed_y), orientation=ori, **cpw))
    TapJunction(design, 'n0', options=dict(pos_x=mm(xf), pos_y=mm(feed_y), stub_up='False',
                                           stub_length=um(g['feed']['stub_um']), **cpw))
    cp_y = feed_y - g['feed']['stub_um'] / 1000 - g['feed']['cp_lead_mm'] - enc_half
    EnclosingCap(design, 'cp', options=dict(cp_options(), **cpw, pos_x=mm(xf), pos_y=mm(cp_y), orientation='270'))
    for name, a, b in (('feed_in', ('feed_l', 'end'), ('n0', 'prime_start')),
                       ('feed_out', ('n0', 'prime_end'), ('feed_r', 'end')),
                       ('stub_p', ('n0', 'second_end'), ('cp', 'a'))):
        RouteStraight(design, name, Dict(pin_inputs=Dict(start_pin=Dict(component=a[0], pin=a[1]),
                                                         end_pin=Dict(component=b[0], pin=b[1])), **cpw))
    pb = np.asarray(design.components['cp'].pins['b']['middle'], float)
    ft = m['filter_tap_mm']
    pr = PurcellReadout(design, 'purcell_resonator_pair',
                        options=dict(module_options(variant, lengths), pos_x=mm(pb[0] - xf), pos_y=mm(pb[1] - ft)))
    qpin = np.asarray(pr.pins['qubit']['middle'], float)
    OpenToGround(design, 'readout_open', options=dict(pos_x=mm(qpin[0]), pos_y=mm(qpin[1]), orientation='270',
                                                      width=cpw['trace_width'], gap=cpw['trace_gap']))
    if wire_bonds:
        t = design.qgeometry.tables['path']
        t['hfss_wire_bonds'] = t['subtract'].astype(bool)
    # drawn geometry of the circuit-model sections (mm)
    n0 = design.components['n0']
    tee_x = float(np.mean([n0.pins['prime_start']['middle'][0], n0.pins['prime_end']['middle'][0]]))
    pa = np.asarray(design.components['cp'].pins['a']['middle'], float)
    ports = [float(design.components['feed_l'].pins['end']['middle'][0]),
             float(design.components['feed_r'].pins['end']['middle'][0])]
    info = dict(variant=variant, lengths=lengths, length_report=pr.length_report, line_spacing_um=spacing_um(variant),
                stub_metal_um=stub_metal_um(variant), feed_y_mm=feed_y,
                circuit_sections_mm=dict(port_x=ports, tee_x=tee_x, tee_to_cp_pin_a=float(feed_y - pa[1])))
    return design, info


def pin_check(design):
    pr = design.components['purcell_resonator_pair']
    out = []
    for pin, (comp, other) in dict(cp=('cp', 'b'), qubit=('readout_open', 'open')).items():
        a, b = pr.pins[pin], design.components[comp].pins[other]
        d = float(np.linalg.norm(np.asarray(a['middle']) - np.asarray(b['middle']))) * 1000
        n = float(np.dot(a['normal'], b['normal']))
        wa, wb = float(a['width']) * 1000, float(b['width']) * 1000
        out.append(dict(module_pin=pin, partner=f'{comp}.{other}', offset_um=d, normal_dot=n, width_um=(wa, wb),
                        ok=bool(d < 1e-3 and n < -0.999 and abs(wa - wb) < 1e-6)))
    return out


def export(variant, out_gds, out_png_prefix, lengths=None):
    """Junction-free simulation-geometry GDS (no cheese), layout figures and the 4 um metal check."""
    import gds_tools
    import check_metal_width as cmw
    design, info = build(variant, lengths=lengths)
    exempt = cmw.design_exemptions(design)
    pins = pin_check(design)
    gds_tools.export_gds(design, str(out_gds), fab=False)
    found = cmw.check(out_gds, exempt=exempt)
    chip = GEOMETRY['chip']
    cy, hx, hy = chip['center_y_mm'], chip['size_x_mm'] / 2, chip['size_y_mm'] / 2
    gds_tools.plot_gds(str(out_gds), f'{out_png_prefix}_layout.png', window=(-hx, hx, cy - hy, cy + hy),
                       figsize=(10, 10 * hy / hx))
    p = design.components['purcell_resonator_pair'].options
    hx0, hy0 = design.parse_value(p.pos_x), design.parse_value(p.pos_y)
    gds_tools.plot_gds(str(out_gds), f'{out_png_prefix}_hammer.png', window=(hx0 - .35, hx0 + .35, hy0 - .25, hy0 + .25),
                       figsize=(9, 6.5), dpi=150)
    return dict(info, pins=pins,
                metal_width=dict(regions_below_4um=[dict(width_um=w, bounds=gm.bounds) for gm, w in found]),
                gds=str(Path(out_gds).relative_to(ROOT)),
                figures=[str(Path(f'{out_png_prefix}_{v}.png').relative_to(ROOT)) for v in ('layout', 'hammer')])


def main():
    import datetime
    rev = GEOMETRY['revision']
    out = {}
    for variant in GEOMETRY['variants']:
        r = export(variant, DESIGN / 'gds' / f'Design4_{rev}_{variant}_sim.gds', DESIGN / 'fig' / f'Design4_{rev}_{variant}')
        out[variant] = r
        print(variant, 'stub', round(r['stub_metal_um'], 3), 'um | lengths',
              {k: round(v, 6) for k, v in r['length_report'].items() if 'mm' in k},
              '| pins ok', all(p['ok'] for p in r['pins']), '| <4um regions', len(r['metal_width']['regions_below_4um']))
    out['created'] = datetime.datetime.now().isoformat(timespec='seconds')
    out['note'] = 'estimate lengths (purcell_readout_model fits); HFSS runs use their own Q3D loads'
    (DESIGN / 'parameters' / f'build_checks_{rev}.json').write_text(json.dumps(out, indent=1, default=float), encoding='utf-8')


if __name__ == '__main__':
    main()
