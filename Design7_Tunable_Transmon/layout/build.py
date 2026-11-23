"""Design7_Tunable_Transmon: 10 x 10 mm 8-port chip — half bare resonators, half fixed-frequency transmons.

Revision 2026-11-23_r3 replaces the tunable SQUID transmons with single-JJ fixed-frequency transmons
(Q1 3.5, Q3 3.8, Q2 4.1, Q4 4.4 GHz targets). The folder name is kept for path stability.

Topology (wafer_pcb_ports, 10 x 10 mm die):
    Top half    (P1 <-> P6, y = +2.0 mm): 4 bare lambda/4 hanger resonators (R1-R4, 6.30-6.90 GHz).
    Bottom half (P2 <-> P5, y = -2.0 mm): 4 fixed-frequency transmon pockets (Q1-Q4, no flux/charge lines)
                                          read out by 4 lambda/4 hanger resonators (RR1-RR4, 6.30-6.90 GHz).
    Spare ports (P3, P4, P7, P8):         Standard BondPads grounded at their neck.

All 4 transmons use ``FixedTransmonPocket`` (layout/components/fixed_transmon_pocket.py) with
``orientation='0'`` so every junction shares the vertical finger-left orientation required by
``TransmonJunction`` (single JJ, ``squid=False``). Stage 1 draft: only the simulation GDS is written;
Dolan hooks / e-beam layers are a Stage 3 addition.

Usage (from repository root, qmetal Python):
    python -X utf8 -B Design7_Tunable_Transmon/layout/build.py
"""

from collections import OrderedDict
from pathlib import Path
import argparse
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
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight  # noqa: E402

from bondpad import add_bond_pad  # noqa: E402
from wafer_pcbport import wafer_pcb_ports  # noqa: E402
from fixed_transmon_pocket import FixedTransmonPocket  # noqa: E402
from gds_tools import export_gds, write_wafer_gds  # noqa: E402
import check_metal_width as cmw  # noqa: E402


# -----------------------------------------------------------------------------
# Module-level knobs (mm unless noted); loaded from parameters/geometry.json
# -----------------------------------------------------------------------------
def load_config():
    geo = json.loads((DESIGN / 'parameters' / 'geometry.json').read_text(encoding='utf-8'))
    qtm = json.loads((DESIGN / 'parameters' / 'quantum.json').read_text(encoding='utf-8'))
    return geo, qtm


def _um(val_um):
    return f'{val_um:g}um'


def _mm(val_mm):
    return f'{val_mm:.9f}mm'


def quarter_wave_mm(f_ghz, v_ph=1.2016e8, load_shorten_mm=0.0):
    """Analytic lambda/4 length in mm for a CPW with phase velocity v_ph."""
    return (v_ph / (4.0 * f_ghz * 1e9)) * 1e3 - load_shorten_mm


def add_hanger_resonator(design, name, x_c, y_feed, side, target_len_mm, geo,
                         end_pin=None):
    """Build a deterministic filleted lambda/4 hanger resonator along y_feed.

    Parameters:
        x_c: channel center x (mm)
        y_feed: feedline center y (mm)
        side: +1 (above feedline) or -1 (below feedline)
        target_len_mm: desired filleted centerline length (mm)
        geo: geometry config dict
        end_pin: None -> creates OpenToGround at end; or (comp_name, pin_name) to connect to a qubit.

    Returns:
        dict with actual filleted length, end coordinates, and meander width.
    """
    cpw_cfg = geo['cpw']
    hr = geo['hanger_resonator']
    w_mm = cpw_cfg['trace_width_um'] / 1000.0
    g_mm = cpw_cfg['trace_gap_um'] / 1000.0
    r_mm = cpw_cfg['fillet_um'] / 1000.0
    strip_mm = hr['coupling_ground_strip_um'] / 1000.0

    l_coup = hr['coupling_length_um'] / 1000.0
    dy_lead = hr['lead_y_um'] / 1000.0
    pitch = hr['meander_pitch_um'] / 1000.0
    dy_exit = hr.get('exit_y_um', hr['lead_y_um']) / 1000.0
    min_trench_gap = hr.get('meander_trench_gap_um', 300.0) / 1000.0
    trench_w = w_mm + 2.0 * g_mm

    # Assert DRC B2: every fixed straight segment >= 2 * fillet
    min_seg = 2.0 * r_mm
    for seg_name, seg_val in (('coupling_length', l_coup), ('dy_lead', dy_lead),
                              ('meander_pitch', pitch), ('dy_exit', dy_exit)):
        if seg_val < min_seg - 1e-9:
            raise ValueError(f'{name}: {seg_name}={seg_val*1e3:.1f} um < 2*fillet ({min_seg*1e3:.1f} um)')

    # Assert DRC B7: trench-to-trench distance between parallel meander segments >= 300 um
    for gap_name, c2c in (('lead_y', dy_lead), ('meander_pitch', pitch)):
        trench_gap = c2c - trench_w
        if trench_gap < min_trench_gap - 1e-9:
            raise ValueError(
                f'{name}: {gap_name} trench-to-trench gap {trench_gap*1e3:.1f} um '
                f'< {min_trench_gap*1e3:.1f} um (DRC B7)'
            )

    # 9 ninety-degree bends in the 4-leg horizontal serpentine
    n_bends = 9
    fillet_loss = n_bends * r_mm * (2.0 - math.pi / 2.0)
    fixed_poly = 1.5 * l_coup + dy_lead + 3.0 * pitch + dy_exit
    w_meander = (target_len_mm + fillet_loss - fixed_poly) / 3.0
    if w_meander / 2.0 < min_seg - 1e-9:
        raise ValueError(f'{name}: half meander width {w_meander*500:.1f} um < 2*fillet ({min_seg*1e3:.1f} um)')

    cpw = dict(trace_width=_um(cpw_cfg['trace_width_um']),
               trace_gap=_um(cpw_cfg['trace_gap_um']))

    y_coup = y_feed + side * (w_mm + 2.0 * g_mm + strip_mm)
    x_start = x_c - l_coup / 2.0
    x_turn = x_c + l_coup / 2.0
    x_left = x_c - w_meander / 2.0
    x_right = x_c + w_meander / 2.0

    y0 = y_coup + side * dy_lead
    y1 = y0 + side * pitch
    y2 = y1 + side * pitch
    y3 = y2 + side * pitch
    y_end = y3 + side * dy_exit

    # Shorted coupling end pointing +x along feedline (ShortToGround orientation='180' -> normal (+1, 0))
    short_name = f'{name}_short'
    ShortToGround(design, short_name, options=dict(
        pos_x=_mm(x_start), pos_y=_mm(y_coup), orientation='180',
        width=cpw['trace_width'], gap=cpw['trace_gap'],
    ))

    if end_pin is None:
        open_name = f'{name}_open'
        # OpenToGround pin faces back toward the incoming route (-side * y): 90 -> (0, -1), 270 -> (0, +1)
        open_ori = '90' if side > 0 else '270'
        OpenToGround(design, open_name, options=dict(
            pos_x=_mm(x_c), pos_y=_mm(y_end), orientation=open_ori,
            width=cpw['trace_width'], gap=cpw['trace_gap'],
            termination_gap=_um(cpw_cfg['trace_gap_um']),
        ))
        end_pin = (open_name, 'open')

    anchors = OrderedDict(enumerate([
        np.array([x_turn, y_coup]),
        np.array([x_turn, y0]),
        np.array([x_left, y0]),
        np.array([x_left, y1]),
        np.array([x_right, y1]),
        np.array([x_right, y2]),
        np.array([x_left, y2]),
        np.array([x_left, y3]),
        np.array([x_c, y3]),
    ]))

    route = RouteAnchors(design, name, options=Dict(
        chip='main',
        fillet=_um(cpw_cfg['fillet_um']),
        lead=Dict(start_straight='0um', end_straight='0um'),
        anchors=anchors,
        pin_inputs=Dict(
            start_pin=Dict(component=short_name, pin='short'),
            end_pin=Dict(component=end_pin[0], pin=end_pin[1]),
        ),
        **cpw,
    ))

    # Measure rendered filleted centerline length from QGeometry
    ptable = design.qgeometry.tables['path']
    row = ptable[(ptable['component'] == route.id) & (~ptable['subtract'])].iloc[0]
    rendered_poly_len = float(row.geometry.length)
    actual_filleted = rendered_poly_len - fillet_loss

    return dict(
        name=name,
        target_len_mm=target_len_mm,
        actual_filleted_mm=actual_filleted,
        meander_width_mm=w_meander,
        x_c_mm=x_c,
        y_end_mm=y_end,
        side=side,
    )


def build_design(*, jj_metal=True, overrides=None):
    """Build the 10 x 10 mm Design7 Qiskit Metal design (fixed-frequency transmons, simulation geometry)."""

    geo, qtm = load_config()
    if overrides:
        geo = json.loads(json.dumps(geo))
        geo.update(overrides)

    design = designs.DesignPlanar(overwrite_enabled=True)
    design.chips.main.size.size_x = f"{geo['chip_size_mm'][0]}mm"
    design.chips.main.size.size_y = f"{geo['chip_size_mm'][1]}mm"
    design.chips.main.size.size_z = '-500um'
    design.chips.main.material = 'silicon'

    cpw_cfg = geo['cpw']
    cpw = dict(trace_width=_um(cpw_cfg['trace_width_um']),
               trace_gap=_um(cpw_cfg['trace_gap_um']))

    # 1. All 8 standard wafer launch pads (active: P1/P6, P2/P5; grounded spares: P3/P4/P7/P8)
    for p_name, placement in wafer_pcb_ports.items():
        add_bond_pad(design, f'LP{p_name}', placement, **cpw)

    # 2. Two straight horizontal feedlines
    RouteStraight(design, 'feed_res', options=Dict(
        chip='main',
        pin_inputs=Dict(
            start_pin=Dict(component='LPP1', pin='tie'),
            end_pin=Dict(component='LPP6', pin='tie'),
        ),
        **cpw,
    ))
    RouteStraight(design, 'feed_qbt', options=Dict(
        chip='main',
        pin_inputs=Dict(
            start_pin=Dict(component='LPP2', pin='tie'),
            end_pin=Dict(component='LPP5', pin='tie'),
        ),
        **cpw,
    ))

    # 3. Top half (y = +2.0 mm): 4 bare lambda/4 hanger resonators
    v_ph = cpw_cfg['assumed_v_ph_m_s']
    res_reports = []
    for ch in geo['bare_resonators']:
        l_target = quarter_wave_mm(ch['f_target_GHz'], v_ph=v_ph, load_shorten_mm=0.0)
        rep = add_hanger_resonator(
            design, ch['name'], ch['x_mm'], +2.0, ch['side'], l_target, geo, end_pin=None
        )
        rep['f_target_GHz'] = ch['f_target_GHz']
        res_reports.append(rep)

    # 4. Bottom half (y = -2.0 mm): 4 fixed-frequency transmon pockets + readout hanger resonators
    tp = geo['transmon_pocket']
    hr = geo['hanger_resonator']
    w_mm = cpw_cfg['trace_width_um'] / 1000.0
    g_mm = cpw_cfg['trace_gap_um'] / 1000.0
    strip_mm = hr['coupling_ground_strip_um'] / 1000.0
    exit_y_um = hr.get('exit_y_um', hr['lead_y_um'])
    dy_total = (w_mm + 2.0 * g_mm + strip_mm) + (hr['lead_y_um'] + exit_y_um + 3.0 * hr['meander_pitch_um']) / 1000.0
    q_pin_offset = (tp['pocket_height_um'] / 2.0 + tp['readout_cpw_extend_um']) / 1000.0

    qbt_reports = []

    for ch in geo['transmons']:
        qname = ch['name']
        rname = ch['res_name']
        side = ch['side']
        y_end = -2.0 + side * dy_total
        y_q = y_end + side * q_pin_offset
        # Readout pad faces back toward the feedline (-side): qubit below feed -> loc_H=+1; above -> loc_H=-1
        loc_h = '+1' if side < 0 else '-1'
        q_in = qtm['qubits'][qname]

        FixedTransmonPocket(design, qname, options=Dict(
            chip='main',
            pos_x=_mm(ch['x_mm']),
            pos_y=_mm(y_q),
            orientation='0',  # MUST stay 0 so all junctions share the finger-left Dolan direction
            pad_width=_um(tp['pad_width_um']),
            pad_height=_um(tp['pad_height_um']),
            pad_gap=_um(tp['pad_gap_um']),
            pad_corner_radius=_um(tp['pad_corner_radius_um']),
            pocket_width=_um(tp['pocket_width_um']),
            pocket_height=_um(tp['pocket_height_um']),
            pocket_corner_radius=_um(tp['pocket_corner_radius_um']),
            jj_leg_width=_um(tp['jj_leg_width_um']),
            jj_cut=_um(tp['jj_cut_um']),
            jj_width=_um(tp['jj_width_um']),
            Lj=f"{q_in['Lj_nH']}nH",
            jj_metal=str(bool(jj_metal)),
            connection_pads=Dict(
                readout=Dict(
                    pad_gap=_um(tp['readout_pad_gap_um']),
                    pad_width=_um(tp['readout_pad_width_um']),
                    pad_height=_um(tp['readout_pad_height_um']),
                    cpw_width=cpw['trace_width'],
                    cpw_gap=cpw['trace_gap'],
                    cpw_extend=_um(tp['readout_cpw_extend_um']),
                    loc_W='0',
                    loc_H=loc_h,
                )
            ),
        ))

        load_shorten = float(hr.get('readout_load_shorten_mm', 0.220))
        l_target = quarter_wave_mm(ch['f_r_target_GHz'], v_ph=v_ph, load_shorten_mm=load_shorten)
        rep = add_hanger_resonator(
            design, rname, ch['x_mm'], -2.0, side, l_target, geo, end_pin=(qname, 'readout')
        )
        rep['qubit'] = qname
        rep['qubit_pos_mm'] = [ch['x_mm'], y_q]
        rep['f_target_GHz'] = ch['f_r_target_GHz']
        rep['f01_target_GHz'] = ch['f01_target_GHz']
        qbt_reports.append(rep)

    design.metadata['geometry_revision'] = geo['geometry_revision']
    design.metadata['bare_resonators'] = res_reports
    design.metadata['transmons'] = qbt_reports
    return design


def verify_pins_and_lengths(design):
    print('\n=== Resonator length check (target vs. filleted) ===')
    for group in ('bare_resonators', 'transmons'):
        for rep in design.metadata[group]:
            err_um = (rep['actual_filleted_mm'] - rep['target_len_mm']) * 1e3
            print(f"  {rep['name']:4s} ({rep['f_target_GHz']:.2f} GHz): "
                  f"target={rep['target_len_mm']:.4f} mm, actual={rep['actual_filleted_mm']:.4f} mm, "
                  f"err={err_um:+.3f} um, W_meander={rep['meander_width_mm']*1e3:.1f} um")
            assert abs(err_um) < 0.5, f"{rep['name']} length error {err_um:.3f} um"

    print('\n=== Qubit readout pin normals ===')
    for q in ('Q1', 'Q2', 'Q3', 'Q4'):
        pin = design.components[q].pins['readout']
        print(f"  {q}.readout: middle=({pin.middle[0]:+.4f}, {pin.middle[1]:+.4f}) mm, "
              f"normal=({pin.normal[0]:+.2f}, {pin.normal[1]:+.2f})")


def _annotate_junctions_layer50(design, gds_path):
    """Add lumped junction sheets from design.qgeometry.tables['junction'] to layer 50/0."""
    import gdstk
    jt = design.qgeometry.tables['junction']
    if jt.empty:
        return
    lib = gdstk.read_gds(str(gds_path))
    top, = lib.top_level()
    unit_um = 1.0 / design.parse_value('1um')
    for row in jt.itertuples():
        geom = row.geometry.buffer(row.width / 2.0, cap_style=2)
        pts = [(x * unit_um, y * unit_um) for x, y in geom.exterior.coords[:-1]]
        top.add(gdstk.Polygon(pts, layer=50, datatype=0))
    lib.write_gds(str(gds_path), max_points=199)


def main():
    parser = argparse.ArgumentParser(description='Build the Design7 layout and export the simulation GDS.')
    parser.parse_args()

    geo, _ = load_config()
    rev = geo['geometry_revision']
    gds_dir = DESIGN / 'gds'
    gds_dir.mkdir(parents=True, exist_ok=True)

    # 1. Clean companion simulation GDS (jj_metal=True, unnotched, layer 50/0 junction annotations)
    d_sim = build_design(jj_metal=True)
    verify_pins_and_lengths(d_sim)
    # Save junction table before export_gds empties it
    jt_backup = d_sim.qgeometry.tables['junction'].copy()
    sim_raw = gds_dir / f'_tmp_Design7_{rev}_sim.gds'
    sim_gds = gds_dir / f'Design7_{rev}_sim.gds'
    export_gds(d_sim, str(sim_raw), fab=True, cheese=False)
    write_wafer_gds(sim_raw, sim_gds)
    sim_raw.unlink(missing_ok=True)
    d_sim.qgeometry.tables['junction'] = jt_backup
    _annotate_junctions_layer50(d_sim, sim_gds)
    print(f'\nWrote simulation companion GDS: {sim_gds}', flush=True)

    # Windowed metal-width DRC check (>= 4 um) on all 4 qubit regions
    exempt = cmw.design_exemptions(d_sim)
    total_viol = 0
    for rep in d_sim.metadata['transmons']:
        qx, qy = rep['qubit_pos_mm'][0] * 1e3, rep['qubit_pos_mm'][1] * 1e3
        win = (qx - 450, qy - 450, qx + 450, qy + 450)
        thin = cmw.check(sim_gds, window=win, exempt=exempt)
        total_viol += len(thin)
        for reg, w in thin[:5]:
            print(f"  VIOLATION in {rep['qubit']}: width={w:.2f} um at bounds={reg.bounds}")
    print(f'Metal-width check (sim GDS qubits, >= 4 um): {total_viol} violation(s)', flush=True)


if __name__ == '__main__':
    main()
