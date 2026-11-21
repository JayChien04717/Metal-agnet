# -*- coding: utf-8 -*-
"""Design8_reflection_purcell_resonator: single-port reflection lambda/2 Purcell filter
coupled to three downward-meandering lambda/4 readout resonators (without qubits).

Topology (6.0 x 5.0 mm chip, 10/6 um CPW):
    1. Top launch pad (`LP_IN`) at the upper edge pointing downward (-y).
    2. Short CPW elbow route (`feed_in`) into pin `b` of a horizontal hexagon-tapered
       interdigital capacitor (`IDC_in`, `layout/components/interdigital_cap.py`).
    3. Reflection lambda/2 Purcell filter (`PF1`): leaves `IDC_in` pin `a` toward -x,
       drops vertically along x_left to y_bus, runs horizontally across above the three
       readout resonators to x_right, and folds into a 5-column vertical serpentine on
       the right, terminating in `PF1_open` (`OpenToGround`).
    4. Three lambda/4 readout resonators (`RR1`, `RR2`, `RR3`): each starts with a
       horizontal coupling section (`ShortToGround` at the left tip) parallel to the
       Purcell bus, meanders downward through 8 horizontal bars, and exits at the bottom
       via an inward dogleg tail terminating in `RRi_open` (`OpenToGround`).

Usage (from repository root, qmetal / gds_only Python):
    python -X utf8 -B Design8_reflection_purcell_resonator/layout/build.py
"""

from collections import OrderedDict
from pathlib import Path
import json
import math
import sys
import warnings

import numpy as np

warnings.filterwarnings('ignore')

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO / 'code' / 'general')]
import project_paths  # noqa: E402

project_paths.configure(__file__)
DESIGN = project_paths.design_dir(__file__)
sys.path[:0] = [str(DESIGN / 'layout' / 'components'), *project_paths.pattern_design_paths()]

from qiskit_metal import Dict, designs  # noqa: E402
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround  # noqa: E402
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround  # noqa: E402
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors  # noqa: E402

from bondpad import add_bond_pad  # noqa: E402
from interdigital_cap import InterdigitalCap, alley_cap_fF  # noqa: E402
from gds_tools import export_gds  # noqa: E402
import check_metal_width as cmw  # noqa: E402


def load_config():
    return json.loads((DESIGN / 'parameters' / 'geometry.json').read_text(encoding='utf-8'))


def _um(val_um):
    return f'{val_um:g}um'


def _mm(val_mm):
    return f'{val_mm:.9f}mm'


def _assert_segments_ge_2fillet(name, pts, r_mm):
    """Assert every consecutive straight segment in ``pts`` is >= 2 * r_mm (DRC rule B2)."""
    min_seg = 2.0 * r_mm
    for i in range(len(pts) - 1):
        seg_len = float(np.linalg.norm(np.asarray(pts[i + 1]) - np.asarray(pts[i])))
        if seg_len < min_seg - 1e-9:
            raise ValueError(
                f'{name}: segment {i}->{i+1} length {seg_len*1e3:.1f} um < 2*fillet ({min_seg*1e3:.1f} um)'
            )


def _rendered_filleted_length_mm(design, route_comp, n_bends, r_mm):
    ptable = design.qgeometry.tables['path']
    row = ptable[(ptable['component'] == route_comp.id) & (~ptable['subtract'])].iloc[0]
    poly_len = float(row.geometry.length)
    fillet_loss = n_bends * r_mm * (2.0 - math.pi / 2.0)
    return poly_len - fillet_loss


def add_reflection_purcell_filter(design, idc_comp, geo):
    """Build the L-bus + right-side 5-column vertical meander lambda/2 Purcell filter."""
    cpw_cfg = geo['cpw']
    pf = geo['purcell_filter']
    r_mm = cpw_cfg['fillet_um'] / 1000.0
    v_ph = cpw_cfg['assumed_v_ph_m_s']
    c_per_m = cpw_cfg['assumed_c_per_m_f']
    z0 = cpw_cfg['assumed_z0_ohm']

    # Analytic + Q3D pi-network admittance of IDC_in at pin a
    idc_cfg = geo['interdigital_cap']
    c_in_alley_ff = alley_cap_fF(
        idc_cfg['finger_count'],
        idc_cfg['overlap_um'],
        finger_width_um=idc_cfg['finger_width_um'],
        h_um=abs(geo['chip']['size_z_um']),
        eps_r=11.9,
    )
    q3d_loads = pf.get('q3d_loads') or {}
    c_in_ff = float(q3d_loads.get('c_in_q3d_fF', c_in_alley_ff))
    c_a_gnd_ff = float(q3d_loads.get('c_a_gnd_deembedded_fF', 0.0))
    c_b_gnd_ff = float(q3d_loads.get('c_b_gnd_raw_fF', c_a_gnd_ff))
    omega = 2.0 * math.pi * pf['f_target_GHz'] * 1e9
    y0 = 1.0 / z0
    yb = y0 + 1j * omega * (c_b_gnd_ff * 1e-15)
    jw_cin = 1j * omega * (c_in_ff * 1e-15)
    y_load_a = (jw_cin * yb) / (jw_cin + yb) + 1j * omega * (c_a_gnd_ff * 1e-15)
    theta_load = math.atan(float(y_load_a.imag) / y0)
    c_eff_f = float(y_load_a.imag) / omega
    # Exact lambda/2 length with open-end phase shift theta_load: beta * L = pi - theta_load
    target_len_mm = (v_ph / (2.0 * pf['f_target_GHz'] * 1e9)) * (1.0 - theta_load / math.pi) * 1e3

    pin_a = idc_comp.pins['a']
    x_a, y_idc = float(pin_a['middle'][0]), float(pin_a['middle'][1])
    x_left = float(pf['x_left_mm'])
    x_right = float(pf['x_right_mm'])
    y_bus = float(pf['y_bus_mm'])
    dy_clear = pf['meander_clear_y_um'] / 1000.0
    pitch = pf['meander_pitch_um'] / 1000.0
    dy_drop = pf['open_tail_drop_um'] / 1000.0

    # 11 ninety-degree bends across the L-bus and 5-column vertical serpentine
    n_bends = 11
    fillet_loss = n_bends * r_mm * (2.0 - math.pi / 2.0)
    fixed_poly = (
        (x_a - x_left)
        + (y_idc - y_bus)
        + (x_right - x_left)
        + dy_clear
        + 4.0 * pitch
        - dy_drop
    )
    h_meander = (target_len_mm + fillet_loss - fixed_poly) / 5.0

    y_bot = y_bus + dy_clear
    y_top = y_bot + h_meander
    y_open = y_top - dy_drop

    x5 = x_right
    x4 = x_right - 1.0 * pitch
    x3 = x_right - 2.0 * pitch
    x2 = x_right - 3.0 * pitch
    x1 = x_right - 4.0 * pitch

    anchor_list = [
        np.array([x_left, y_idc]),
        np.array([x_left, y_bus]),
        np.array([x5, y_bus]),
        np.array([x5, y_top]),
        np.array([x4, y_top]),
        np.array([x4, y_bot]),
        np.array([x3, y_bot]),
        np.array([x3, y_top]),
        np.array([x2, y_top]),
        np.array([x2, y_bot]),
        np.array([x1, y_bot]),
    ]
    all_pts = [np.array([x_a, y_idc]), *anchor_list, np.array([x1, y_open])]
    _assert_segments_ge_2fillet(pf['name'], all_pts, r_mm)

    cpw = dict(trace_width=_um(cpw_cfg['trace_width_um']),
               trace_gap=_um(cpw_cfg['trace_gap_um']))

    open_name = f"{pf['name']}_open"
    OpenToGround(design, open_name, options=dict(
        pos_x=_mm(x1),
        pos_y=_mm(y_open),
        orientation='90',  # incoming vertical column approaches from -y -> pin normal (0, -1)
        width=cpw['trace_width'],
        gap=cpw['trace_gap'],
        termination_gap=cpw['trace_gap'],
    ))

    route = RouteAnchors(design, pf['name'], options=Dict(
        chip='main',
        fillet=_um(cpw_cfg['fillet_um']),
        lead=Dict(start_straight='0um', end_straight='0um'),
        anchors=OrderedDict(enumerate(anchor_list)),
        pin_inputs=Dict(
            start_pin=Dict(component=idc_comp.name, pin='a'),
            end_pin=Dict(component=open_name, pin='open'),
        ),
        **cpw,
    ))

    actual_len_mm = _rendered_filleted_length_mm(design, route, n_bends, r_mm)
    return dict(
        name=pf['name'],
        mode=pf['mode'],
        f_target_GHz=pf['f_target_GHz'],
        c_in_alley_fF=round(c_in_ff, 3),
        c_in_eff_shunt_fF=round(c_eff_f * 1e15, 3),
        target_len_mm=target_len_mm,
        actual_filleted_mm=actual_len_mm,
        meander_height_mm=h_meander,
        y_bus_mm=y_bus,
        y_top_mm=y_top,
        open_pos_mm=[x1, y_open],
    )


def add_downward_readout_resonator(design, ch_cfg, geo):
    """Build one downward-meandering lambda/4 readout resonator coupled below the Purcell bus."""
    cpw_cfg = geo['cpw']
    pf = geo['purcell_filter']
    rr = geo['readout_resonators']

    w_mm = cpw_cfg['trace_width_um'] / 1000.0
    g_mm = cpw_cfg['trace_gap_um'] / 1000.0
    r_mm = cpw_cfg['fillet_um'] / 1000.0
    v_ph = cpw_cfg['assumed_v_ph_m_s']

    strip_mm = rr['coupling_ground_strip_um'] / 1000.0
    l_coup = rr['coupling_length_um'] / 1000.0
    pitch = rr['meander_pitch_um'] / 1000.0
    dy_tail1 = rr['tail_drop1_um'] / 1000.0
    dy_tail2 = rr['tail_drop2_um'] / 1000.0

    name = ch_cfg['name']
    x_c = float(ch_cfg['x_c_mm'])
    f_ghz = float(ch_cfg['f_target_GHz'])
    d_turn8 = ch_cfg['bar8_turn_offset_um'] / 1000.0
    dx_jog = ch_cfg['tail_jog_x_um'] / 1000.0

    target_len_mm = (v_ph / (4.0 * f_ghz * 1e9)) * 1e3

    # 19 ninety-degree bends: 1 top coupling turn + 8*2 meander turns + 2 dogleg tail turns
    n_bends = 19
    fillet_loss = n_bends * r_mm * (2.0 - math.pi / 2.0)
    fixed_poly = l_coup + 8.0 * pitch + d_turn8 + dy_tail1 + abs(dx_jog) + dy_tail2
    w_meander = (target_len_mm + fillet_loss - fixed_poly) / 7.5

    y_bus = float(pf['y_bus_mm'])
    y_coup = y_bus - (w_mm + 2.0 * g_mm + strip_mm)
    x_l = x_c - w_meander / 2.0
    x_r = x_c + w_meander / 2.0
    x_short = x_r - l_coup
    x_turn8 = x_c + d_turn8
    x_end = x_turn8 + dx_jog

    y_bars = [y_coup - k * pitch for k in range(1, 9)]  # y_1 .. y_8
    y_tail_mid = y_bars[7] - dy_tail1
    y_end = y_tail_mid - dy_tail2

    anchor_list = [
        np.array([x_r, y_coup]),
        np.array([x_r, y_bars[0]]),
        np.array([x_l, y_bars[0]]),
        np.array([x_l, y_bars[1]]),
        np.array([x_r, y_bars[1]]),
        np.array([x_r, y_bars[2]]),
        np.array([x_l, y_bars[2]]),
        np.array([x_l, y_bars[3]]),
        np.array([x_r, y_bars[3]]),
        np.array([x_r, y_bars[4]]),
        np.array([x_l, y_bars[4]]),
        np.array([x_l, y_bars[5]]),
        np.array([x_r, y_bars[5]]),
        np.array([x_r, y_bars[6]]),
        np.array([x_l, y_bars[6]]),
        np.array([x_l, y_bars[7]]),
        np.array([x_turn8, y_bars[7]]),
        np.array([x_turn8, y_tail_mid]),
        np.array([x_end, y_tail_mid]),
    ]
    all_pts = [np.array([x_short, y_coup]), *anchor_list, np.array([x_end, y_end])]
    _assert_segments_ge_2fillet(name, all_pts, r_mm)

    cpw = dict(trace_width=_um(cpw_cfg['trace_width_um']),
               trace_gap=_um(cpw_cfg['trace_gap_um']))

    short_name = f'{name}_short'
    ShortToGround(design, short_name, options=dict(
        pos_x=_mm(x_short),
        pos_y=_mm(y_coup),
        orientation='180',  # outward pin normal points +x
        width=cpw['trace_width'],
        gap=cpw['trace_gap'],
    ))

    open_name = f'{name}_open'
    OpenToGround(design, open_name, options=dict(
        pos_x=_mm(x_end),
        pos_y=_mm(y_end),
        orientation='270',  # incoming route approaches from +y -> outward pin normal (0, +1)
        width=cpw['trace_width'],
        gap=cpw['trace_gap'],
        termination_gap=cpw['trace_gap'],
    ))

    route = RouteAnchors(design, name, options=Dict(
        chip='main',
        fillet=_um(cpw_cfg['fillet_um']),
        lead=Dict(start_straight='0um', end_straight='0um'),
        anchors=OrderedDict(enumerate(anchor_list)),
        pin_inputs=Dict(
            start_pin=Dict(component=short_name, pin='short'),
            end_pin=Dict(component=open_name, pin='open'),
        ),
        **cpw,
    ))

    actual_len_mm = _rendered_filleted_length_mm(design, route, n_bends, r_mm)
    return dict(
        name=name,
        f_target_GHz=f_ghz,
        target_len_mm=target_len_mm,
        actual_filleted_mm=actual_len_mm,
        meander_width_mm=w_meander,
        coupling_length_mm=l_coup,
        x_c_mm=x_c,
        short_pos_mm=[x_short, y_coup],
        open_pos_mm=[x_end, y_end],
    )


def build_design(*, overrides=None, wire_bonds=True):
    """Build the ``Design8_reflection_purcell_resonator`` Qiskit Metal design."""
    geo = load_config()
    if overrides:
        geo = json.loads(json.dumps(geo))
        geo.update(overrides)

    chip = geo['chip']
    cpw_cfg = geo['cpw']
    r_mm = cpw_cfg['fillet_um'] / 1000.0
    cpw = dict(trace_width=_um(cpw_cfg['trace_width_um']),
               trace_gap=_um(cpw_cfg['trace_gap_um']))

    design = designs.DesignPlanar(overwrite_enabled=True)
    design.chips.main.material = chip['material']
    design.chips.main.size.update(
        size_x=_mm(chip['size_x_mm']),
        size_y=_mm(chip['size_y_mm']),
        size_z=_um(chip['size_z_um']),
        center_x=_mm(chip['center_x_mm']),
        center_y=_mm(chip['center_y_mm']),
    )
    design.variables.update(cpw_width=cpw['trace_width'], cpw_gap=cpw['trace_gap'])

    # 1. Top input wire-bond launch pad pointing downward (-y)
    lp_cfg = geo['launch_pad']
    lp_comp = add_bond_pad(
        design,
        lp_cfg['name'],
        dict(
            pos_x=_mm(lp_cfg['pos_x_mm']),
            pos_y=_mm(lp_cfg['pos_y_mm']),
            orientation=str(lp_cfg['orientation_deg']),
        ),
        **cpw,
    )

    # 2. Hexagon-tapered interdigital capacitor (IDC_in)
    idc_cfg = geo['interdigital_cap']
    idc_comp = InterdigitalCap(design, idc_cfg['name'], options=Dict(
        chip='main',
        pos_x=_mm(idc_cfg['pos_x_mm']),
        pos_y=_mm(idc_cfg['pos_y_mm']),
        orientation=str(idc_cfg['orientation_deg']),
        finger_count=str(idc_cfg['finger_count']),
        finger_width=_um(idc_cfg['finger_width_um']),
        finger_gap=_um(idc_cfg['finger_gap_um']),
        overlap=_um(idc_cfg['overlap_um']),
        end_gap=_um(idc_cfg['end_gap_um']),
        bus_width=_um(idc_cfg['bus_width_um']),
        taper_length=_um(idc_cfg['taper_length_um']),
        lead_length=_um(idc_cfg['lead_length_um']),
        ground_gap=_um(idc_cfg['ground_gap_um']),
        **cpw,
    ))

    # 3. Input elbow CPW from LP_IN.tie -> IDC_in.b
    p_tie = np.asarray(lp_comp.pins['tie']['middle'], float)
    p_idc_b = np.asarray(idc_comp.pins['b']['middle'], float)
    elbow = np.array([p_tie[0], p_idc_b[1]])
    _assert_segments_ge_2fillet('feed_in', [p_tie, elbow, p_idc_b], r_mm)
    RouteAnchors(design, 'feed_in', options=Dict(
        chip='main',
        fillet=_um(cpw_cfg['fillet_um']),
        lead=Dict(start_straight='0um', end_straight='0um'),
        anchors=OrderedDict([(0, elbow)]),
        pin_inputs=Dict(
            start_pin=Dict(component=lp_comp.name, pin='tie'),
            end_pin=Dict(component=idc_comp.name, pin='b'),
        ),
        **cpw,
    ))

    # 4. Reflection lambda/2 Purcell filter (L-bus + right 5-column vertical meander)
    pf_report = add_reflection_purcell_filter(design, idc_comp, geo)

    # 5. Three downward-meandering lambda/4 readout resonators coupled to the Purcell bus
    rr_reports = [
        add_downward_readout_resonator(design, ch_cfg, geo)
        for ch_cfg in geo['readout_resonators']['channels']
    ]

    if wire_bonds:
        t = design.qgeometry.tables['path']
        t['hfss_wire_bonds'] = t['subtract'].astype(bool)

    design.metadata['geometry_revision'] = geo['geometry_revision']
    design.metadata['purcell_filter'] = pf_report
    design.metadata['readout_resonators'] = rr_reports
    return design


def verify_layout(design):
    pf = design.metadata['purcell_filter']
    err_pf_um = (pf['actual_filleted_mm'] - pf['target_len_mm']) * 1e3
    print('\n=== Reflection Purcell Filter (lambda/2 open) ===')
    print(f"  {pf['name']} ({pf['f_target_GHz']:.2f} GHz): "
          f"C_in(Alley)={pf['c_in_alley_fF']:.1f} fF, "
          f"target={pf['target_len_mm']:.4f} mm, actual={pf['actual_filleted_mm']:.4f} mm, "
          f"err={err_pf_um:+.3f} um, H_meander={pf['meander_height_mm']*1e3:.1f} um")
    assert abs(err_pf_um) < 0.5, f"Purcell filter length error {err_pf_um:.3f} um"

    print('\n=== Readout Resonators (lambda/4 coupled to Purcell bus) ===')
    for rep in design.metadata['readout_resonators']:
        err_um = (rep['actual_filleted_mm'] - rep['target_len_mm']) * 1e3
        print(f"  {rep['name']} ({rep['f_target_GHz']:.2f} GHz): "
              f"target={rep['target_len_mm']:.4f} mm, actual={rep['actual_filleted_mm']:.4f} mm, "
              f"err={err_um:+.3f} um, W_meander={rep['meander_width_mm']*1e3:.1f} um")
        assert abs(err_um) < 0.5, f"{rep['name']} length error {err_um:.3f} um"


def main():
    geo = load_config()
    rev = geo['geometry_revision']
    gds_dir = DESIGN / 'gds'
    gds_dir.mkdir(parents=True, exist_ok=True)

    design = build_design()
    verify_layout(design)

    sim_gds = gds_dir / f'{DESIGN.name}_{rev}_sim.gds'
    export_gds(design, str(sim_gds), fab=False)
    print(f'\nWrote simulation companion GDS: {sim_gds}')

    exempt = cmw.design_exemptions(design)
    thin = cmw.check(sim_gds, exempt=exempt)
    print(f'Metal-width check (sim GDS, >= 4 um): {len(thin)} violation(s)')
    if thin:
        for reg, w in thin[:10]:
            print(f'  VIOLATION: width={w:.2f} um at bounds={reg.bounds}')


if __name__ == '__main__':
    main()
