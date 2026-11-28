# -*- coding: utf-8 -*-
"""Design14_FlipChipSurface17: Dual-Tier Flip-Chip Rotated d=3 Surface-17 Layout (Stage 1).

Architecture:
  - Q-tier (Top Qubit Chip, 8 x 8 mm):
      ONLY 17 pure circular wide-gap ScalablePadQubits (256 um pad, 46 um gap, NO coplanar readout pad)
      + 24 FloatingTunableCouplers (flux_as_ground=True, NO on-chip flux lines, NO airbridges!)
      + Indium micro-bump ground-stitching grid (Layer 30).
  - C-tier (Bottom Control & Readout Interposer Chip, 16 x 16 mm):
      56 perimeter launch pads + 4 shared-bus lambda/2 Purcell filters (6.85-7.18 GHz)
      + 17 VerticalReadoutCouplers feeding 17 exact lambda/4 DRC-B7 readout resonators (zero series airbridges)
      + 17 vertical XY drive stubs + 24 VerticalFluxLoops directly beneath every coupler SQUID
      + matching Indium micro-bump grid (Layer 30).

Run (qmetal env):
    python -I -X utf8 -B Design14_FlipChipSurface17/layout/build.py
"""
from __future__ import annotations

import math
import sys
from collections import OrderedDict
from pathlib import Path

import gdspy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly, Rectangle as MplRect
import numpy as np
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

HERE = Path(__file__).resolve().parent
DESIGN_DIR = HERE.parent
REPO_ROOT = DESIGN_DIR.parent
D13_LAYOUT = REPO_ROOT / 'Design13_ScalableFloatingCoupler2Q' / 'layout'
for _p in (
    str(REPO_ROOT),
    str(HERE),
    str(D13_LAYOUT),
    str(D13_LAYOUT / 'components'),
    str(REPO_ROOT / 'pattern design' / 'ebeam_junction'),
    str(REPO_ROOT / 'pattern design' / 'fluxline'),
    str(REPO_ROOT / 'pattern design' / 'qubit'),
):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from qiskit_metal import Dict, designs
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight

import importlib.util
def _load_mod(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_d13 = _load_mod(D13_LAYOUT / 'build.py', 'd13_build_helpers')
BondPad = _d13.BondPad
EnclosingCap = _d13.EnclosingCap

from scalable_pad_qubit import ScalablePadQubit
from floating_tunable_coupler import FloatingTunableCoupler
from airbridge_crossover import AirbridgeSpan
from components.flipchip_couplers import VerticalReadoutCoupler, VerticalFluxLoop

# Reuse verified geometry helpers and Surface-17 lattice spec from Design13
SURFACE17_QUBITS = _d13.SURFACE17_QUBITS
PITCH_MM = _d13.PITCH_MM
CPW_RO = _d13.CPW_RO
CPW_CTL = _d13.CPW_CTL
ENC_OPT = _d13.ENC_OPT
QUBIT_OPT = _d13.QUBIT_OPT
COUPLER_OPT = _d13.COUPLER_OPT
FILLET_S17 = _d13.FILLET_S17
V_PH_MM_GHZ = _d13.V_PH_MM_GHZ
_pos = _d13._pos
_surface17_edges = _d13._surface17_edges
_perimeter_pads_14mm = _d13._perimeter_pads_14mm
_coupler_orientation = _d13._coupler_orientation
_octilinear_fanout = _d13._octilinear_fanout
_xy_to_diag_uv = _d13._xy_to_diag_uv
_build_diag_lambda4_resonator = _d13._build_diag_lambda4_resonator
_segment_intersections_general = _d13._segment_intersections_general
_draw_overpass_with_bridges = _d13._draw_overpass_with_bridges
_assert_no_acute_bends = _d13._assert_no_acute_bends
_assert_zero_unbridged_cpw_crossings = _d13._assert_zero_unbridged_cpw_crossings
_assert_max_airbridge_length = _d13._assert_max_airbridge_length
_place_ground_airbridges = _d13._place_ground_airbridges

REVISION = 'r1'
DATE_TAG = '2026-11-27'
QTIER_DIE_MM = 8.0
CTIER_DIE_MM = 16.0
BUMP_PITCH_MM = 0.35
BUMP_SIZE_MM = 0.024
BUMP_LAYER = 30


def _assert_ctier_zero_pad_collisions(design):
    """Hard Stage-1 DRC: verify no C-tier CPW route intersects any VerticalReadoutCoupler, VerticalFluxLoop, or EnclosingCap."""
    id_to_name = {c.id: n for n, c in design.components.items()}
    poly_tbl = design.qgeometry.tables['poly']
    path_tbl = design.qgeometry.tables['path']

    keepouts = []
    for row in poly_tbl.itertuples():
        cname = id_to_name.get(row.component, '')
        if int(row.layer) != 1:
            continue
        if not row.subtract and cname.startswith(('VRO_', 'VFL_')):
            keepouts.append((cname, row.name, row.geometry))
        elif row.subtract and cname.startswith('PF_c'):
            keepouts.append((cname, 'cap_pocket', row.geometry.buffer(0.018)))

    collisions = []
    for row in path_tbl.itertuples():
        if row.subtract or int(row.layer) != 1:
            continue
        rname = id_to_name.get(row.component, str(row.component))
        trace_poly = row.geometry.buffer(float(row.width) / 2.0, cap_style=2)
        comp = design.components[rname] if rname in design.components else None
        attached = set()
        if comp is not None and 'pin_inputs' in comp.options:
            for side_k in ('start_pin', 'end_pin'):
                pi = comp.options.pin_inputs.get(side_k, {})
                if pi and pi.get('component'):
                    attached.add(pi['component'])
        for cname, pname, geom in keepouts:
            if cname != rname and cname not in attached and trace_poly.intersects(geom):
                collisions.append((rname, cname, pname))
    if collisions:
        sample = ', '.join(f'{r}x{c}.{p}' for r, c, p in collisions[:12])
        raise AssertionError(f'Found {len(collisions)} C-tier route-to-component collisions: {sample}')


def _collect_keepout_union(design, buffer_mm: float = 0.11):
    """Union of all Layer-1 etched pockets & CPW cuts buffered by buffer_mm for clean Indium bump placement."""
    geoms = []
    for row in design.qgeometry.tables['poly'].itertuples():
        if int(row.layer) == 1:
            geoms.append(row.geometry.buffer(buffer_mm))
    for row in design.qgeometry.tables['path'].itertuples():
        if int(row.layer) == 1:
            w = float(row.width) * 0.5 + buffer_mm
            geoms.append(row.geometry.buffer(w))
    return unary_union(geoms) if geoms else Polygon()


def compute_shared_indium_bumps(q_design, c_design, pitch_mm=BUMP_PITCH_MM, clearance_mm=0.105):
    """Compute (x, y) coordinates inside the 8x8 mm Q-tier footprint that are clear ground on BOTH tiers."""
    ku_q = _collect_keepout_union(q_design, buffer_mm=clearance_mm)
    ku_c = _collect_keepout_union(c_design, buffer_mm=clearance_mm)
    ku_all = ku_q.union(ku_c)
    half = 0.5 * QTIER_DIE_MM - 0.25
    xs = np.arange(-half, half + 1e-6, pitch_mm)
    ys = np.arange(-half, half + 1e-6, pitch_mm)
    pts = []
    for x in xs:
        for y in ys:
            if not ku_all.contains(Point(float(x), float(y))):
                pts.append((round(float(x), 4), round(float(y), 4)))
    return pts


def build_qtier(variant='sim'):
    """Build Top Qubit Chip (8 x 8 mm): ONLY 17 circular qubits + 24 floating tunable couplers!"""
    jjm = 'True' if variant == 'sim' else 'False'
    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon_eps11p45'
    design.chips.main.size.update(dict(
        size_x=f'{QTIER_DIE_MM:g}mm', size_y=f'{QTIER_DIE_MM:g}mm',
        size_z='-500um', center_x='0mm', center_y='0mm',
    ))

    # 1. Place 17 PURE circular wide-gap ScalablePadQubits (NO coplanar readout pad!)
    q_opt = dict(QUBIT_OPT, make_readout_pad='False')
    for qname, qinfo in SURFACE17_QUBITS.items():
        ScalablePadQubit(design, qname, options=dict(
            **_pos(qinfo['i'] * PITCH_MM, qinfo['j'] * PITCH_MM),
            Lj=qinfo['Lj'], jj_metal=jjm, **q_opt,
        ))

    # 2. Place 24 FloatingTunableCouplers with NO on-chip flux line (flux_as_ground=True -> pure floating coupler!)
    edges = _surface17_edges()
    for na, nb in edges:
        qa, qb = SURFACE17_QUBITS[na], SURFACE17_QUBITS[nb]
        xm = 0.5 * (qa['i'] + qb['i']) * PITCH_MM
        ym = 0.5 * (qa['j'] + qb['j']) * PITCH_MM
        is_horiz = (qa['j'] == qb['j'])
        orient = _coupler_orientation(na, nb, xm, ym, is_horiz)
        FloatingTunableCoupler(design, f'TC_{na}_{nb}', options=dict(
            **_pos(xm, ym), orientation=orient, jj_metal=jjm, flux_as_ground='True', **COUPLER_OPT,
        ))

    design._design14_summary = dict(
        tier='Q-tier', revision=REVISION, die_mm=(QTIER_DIE_MM, QTIER_DIE_MM),
        n_qubits=len(SURFACE17_QUBITS), n_couplers=len(edges),
        n_control_lines=0, n_resonators=0, n_airbridges=0,
    )
    return design


def build_ctier(variant='sim'):
    """Build Bottom Control & Readout Interposer Chip (16 x 16 mm) aligned to Q-tier."""
    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon_eps11p45'
    design.chips.main.size.update(dict(
        size_x=f'{CTIER_DIE_MM:g}mm', size_y=f'{CTIER_DIE_MM:g}mm',
        size_z='-500um', center_x='0mm', center_y='0mm',
    ))
    design.variables.update(cpw_width=CPW_RO['trace_width'], cpw_gap=CPW_RO['trace_gap'])

    # 1. 56 perimeter launch pads on 16 x 16 mm interposer die
    lp_opt = dict(pad_width='240um', pad_height='260um', pad_gap='110um',
                  taper_height='200um', back_gap='50um', lead_length='30um', **CPW_CTL)
    for pname, pdict in _perimeter_pads_14mm().items():
        BondPad(design, 'LP_' + pname, options=dict(pdict, **lp_opt))

    # 2. 17 VerticalReadoutCouplers positioned directly beneath each Q-tier qubit center
    for qname, qinfo in SURFACE17_QUBITS.items():
        qx, qy = qinfo['i'] * PITCH_MM, qinfo['j'] * PITCH_MM
        VerticalReadoutCoupler(design, f'VRO_{qname}', options=dict(
            **_pos(qx, qy), pad_diameter='100um', pocket_gap='28um',
            ro_corner=qinfo['ro_corner'], lead_length='42um', **CPW_RO,
        ))

    # 3. 24 VerticalFluxLoops positioned directly beneath each Q-tier coupler SQUID center!
    edges = _surface17_edges()
    vfl_list = []
    for na, nb in edges:
        qa, qb = SURFACE17_QUBITS[na], SURFACE17_QUBITS[nb]
        xm = 0.5 * (qa['i'] + qb['i']) * PITCH_MM
        ym = 0.5 * (qa['j'] + qb['j']) * PITCH_MM
        is_horiz = (qa['j'] == qb['j'])
        orient = _coupler_orientation(na, nb, xm, ym, is_horiz)
        vfl_name = f'VFL_{na}_{nb}'
        vfl_list.append((vfl_name, na, nb, xm, ym, orient))
        VerticalFluxLoop(design, vfl_name, options=dict(
            **_pos(xm, ym), orientation=orient, bar_length='72um', bar_width='6um',
            bar_gap='6um', lead_length='55um', **CPW_CTL,
        ))

    # 4. 4 Diagonal (45-deg) Shared-Bus lambda/2 Purcell Filters
    u_purcell = 4.95
    purcell_quads = OrderedDict([
        ('NW', dict(sx=-1, sy=+1, p_in='LP_N00', p_out='LP_W13', v_span=(-2.35, +2.35), cin_ori='-90', cout_ori='180')),
        ('NE', dict(sx=+1, sy=+1, p_in='LP_E13', p_out='LP_N13', v_span=(-2.65, +2.65), cin_ori='180', cout_ori='90')),
        ('SW', dict(sx=-1, sy=-1, p_in='LP_W00', p_out='LP_S00', v_span=(-2.35, +2.35), cin_ori='0',   cout_ori='-90')),
        ('SE', dict(sx=+1, sy=-1, p_in='LP_S13', p_out='LP_E00', v_span=(-2.35, +2.35), cin_ori='90',  cout_ori='0')),
    ])
    planned_routes = []

    for qk, qs in purcell_quads.items():
        sx, sy = qs['sx'], qs['sy']
        pin_in = design.components[qs['p_in']].pins['tie']
        pin_out = design.components[qs['p_out']].pins['tie']
        c_in_xy = tuple(np.asarray(pin_in.middle, float) + 0.45 * np.asarray(pin_in.normal, float))
        c_out_xy = tuple(np.asarray(pin_out.middle, float) + 0.45 * np.asarray(pin_out.normal, float))
        EnclosingCap(design, f'PF_cin_{qk}', options=dict(**_pos(*c_in_xy), orientation=qs['cin_ori'], **ENC_OPT, **CPW_RO))
        EnclosingCap(design, f'PF_cout_{qk}', options=dict(**_pos(*c_out_xy), orientation=qs['cout_ori'], **ENC_OPT, **CPW_RO))
        RouteStraight(design, f'pf_feed_in_{qk}', Dict(
            pin_inputs=Dict(start_pin=Dict(component=qs['p_in'], pin='tie'),
                            end_pin=Dict(component=f'PF_cin_{qk}', pin='a')), **CPW_RO,
        ))
        RouteStraight(design, f'pf_feed_out_{qk}', Dict(
            pin_inputs=Dict(start_pin=Dict(component=f'PF_cout_{qk}', pin='b'),
                            end_pin=Dict(component=qs['p_out'], pin='tie')), **CPW_RO,
        ))
        S = u_purcell * math.sqrt(2.0)
        x_c = 6.05
        y_c = S - x_c
        x_outer = x_c + 0.42
        cin_b = design.components[f'PF_cin_{qk}'].pins['b']
        cout_a = design.components[f'PF_cout_{qk}'].pins['a']
        b_xy = tuple(np.asarray(cin_b.middle, float))
        a_xy = tuple(np.asarray(cout_a.middle, float))
        if abs(cin_b.normal[0]) > 0.5:
            pf_pts = [
                (sx * x_outer, b_xy[1]), (sx * x_outer, sy * y_c), (sx * x_c, sy * y_c),
                (sx * y_c, sy * x_c), (sx * y_c, sy * x_outer), (a_xy[0], sy * x_outer),
            ]
        else:
            pf_pts = [
                (b_xy[0], sy * x_outer), (sx * y_c, sy * x_outer), (sx * y_c, sy * x_c),
                (sx * x_c, sy * y_c), (sx * x_outer, sy * y_c), (sx * x_outer, a_xy[1]),
            ]
        planned_routes.append(dict(
            name=f'PF_body_{qk}', priority=0,
            start=(f'PF_cin_{qk}', 'b'), end=(f'PF_cout_{qk}', 'a'),
            anchors=pf_pts, cpw=CPW_RO, fillet='70um',
        ))

    # 5. 17 EXACT lambda/4 DRC-B7 Readout Resonators (ALL 17 straight out of VRO_* with ZERO coupler-neck bridges!)
    quad_lists = {k: [] for k in purcell_quads}
    for qname, qinfo in SURFACE17_QUBITS.items():
        quad_lists[qinfo['ro_quad']].append(qname)

    res_report = {}
    for quad, qnames in quad_lists.items():
        sx, sy = purcell_quads[quad]['sx'], purcell_quads[quad]['sy']
        qnames.sort(key=lambda qn: _xy_to_diag_uv(
            SURFACE17_QUBITS[qn]['i'] * PITCH_MM, SURFACE17_QUBITS[qn]['j'] * PITCH_MM, sx, sy
        )[1])
        v_lo, v_hi = purcell_quads[quad]['v_span']
        v_slots = np.linspace(v_lo + 0.48, v_hi - 1.05, len(qnames))
        tang_deg = math.degrees(math.atan2(sx, -sy))
        for qi, qname in enumerate(qnames):
            qinfo = SURFACE17_QUBITS[qname]
            L_target = V_PH_MM_GHZ / (4.0 * qinfo['f_ro_GHz'])
            p_pin = tuple(np.asarray(design.components[f'VRO_{qname}'].pins['readout'].middle, float))
            if qname == 'D5':
                # Jog around VRO_D6 at (1.20, 1.20) on the NE diagonal cleanly (+0.38 mm lateral offset!)
                p_jog1 = (0.52, 0.52)
                p_jog2 = (0.85, 0.52)
                p_jog3 = (1.55, 1.22)
                m_anc, short_xy, L_tot = _build_diag_lambda4_resonator(
                    p_jog3, sx, sy, float(v_slots[qi]), u_purcell, L_target - 1.55,
                    fillet_mm=0.055, pitch_mm=0.35, coup_leg_mm=0.25,
                )
                m_anc = [p_jog1, p_jog2, p_jog3] + m_anc
                L_tot += 1.55
            else:
                m_anc, short_xy, L_tot = _build_diag_lambda4_resonator(
                    p_pin, sx, sy, float(v_slots[qi]), u_purcell, L_target,
                    fillet_mm=0.055, pitch_mm=0.35, coup_leg_mm=0.25,
                )
            ShortToGround(design, f'ro_short_{qname}', options=dict(
                **_pos(*short_xy), orientation=f'{tang_deg:.1f}', width=CPW_RO['trace_width'],
            ))
            planned_routes.append(dict(
                name=f'ro_{qname}', priority=0,
                start=(f'VRO_{qname}', 'readout'), end=(f'ro_short_{qname}', 'short'),
                anchors=m_anc, cpw=CPW_RO, fillet=FILLET_S17,
            ))
            res_report[qname] = dict(f_GHz=qinfo['f_ro_GHz'], L_target_mm=round(L_target, 3), L_actual_mm=round(L_tot, 3))

    # 6 & 7. MONOTONIC PERIMETER PAD SORTING (17 vertical XY + 24 vertical Flux sorted TOGETHER per side!)
    pad_pools = {
        'N': [f'N{k:02d}' for k in range(1, 13)],
        'S': [f'S{k:02d}' for k in range(1, 13)],
        'W': [f'W{k:02d}' for k in range(1, 13)],
        'E': [f'E{k:02d}' for k in range(1, 13)],
    }
    side_exits = {'N': [], 'S': [], 'W': [], 'E': []}
    P = PITCH_MM

    for qname, qinfo in SURFACE17_QUBITS.items():
        side = qinfo['xy_side']
        qx, qy = qinfo['i'] * PITCH_MM, qinfo['j'] * PITCH_MM
        if side == 'N':
            if qname == 'X1':
                sx, sy = qx - 0.18, qy + 0.04
            elif qname == 'Z3':
                sx, sy = qx - 0.16, qy + 0.10
            else:
                sx, sy = qx + 0.16, qy + 0.10
            ori = '-90'
        elif side == 'S':
            if qname == 'X2':
                sx, sy = qx + 0.18, qy - 0.04
            elif qname == 'Z4':
                sx, sy = qx + 0.16, qy - 0.10
            else:
                sx, sy = qx - 0.16, qy - 0.10
            ori = '90'
        elif side == 'W':
            if qname == 'D5':
                sx, sy = qx - 0.14, qy + 0.14
            elif qname == 'Z1':
                sx, sy = qx - 0.04, qy - 0.18
            elif qname == 'X3':
                sx, sy = qx - 0.10, qy - 0.16
            else:
                sx, sy = qx - 0.10, qy + 0.16
            ori = '0'
        else:
            if qname == 'Z2':
                sx, sy = qx + 0.04, qy + 0.18
            elif qname == 'X4':
                sx, sy = qx + 0.10, qy + 0.16
            else:
                sx, sy = qx + 0.10, qy - 0.16
            ori = '180'
        OpenToGround(design, f'vxy_{qname}', options=dict(
            **_pos(sx, sy), orientation=ori, width=CPW_CTL['trace_width'], gap=CPW_CTL['trace_gap'],
        ))
        p0 = tuple(np.asarray(design.components[f'vxy_{qname}'].pins['open'].middle, float))
        lead = 2.15 * P - abs(p0[1]) if (side in ('N', 'S') and abs(p0[1]) < 2.0 * P) else (
               2.15 * P - abs(p0[0]) if (side in ('E', 'W') and abs(p0[0]) < 2.0 * P) else 0.35)
        side_exits[side].append(dict(
            sort_key=p0[0] if side in ('N', 'S') else p0[1],
            name=f'xy_{qname}', start=(f'vxy_{qname}', 'open'),
            p_exit=p0, pre=[], lead=max(0.35, lead),
        ))

    orient_to_side = {'180': 'N', '0': 'S', '270': 'W', '90': 'E'}
    for vfl_name, na, nb, xm, ym, orient in vfl_list:
        side = orient_to_side[orient]
        pin_xy = tuple(np.asarray(design.components[vfl_name].pins['flux'].middle, float))
        jog = []
        if side in ('N', 'S') and abs(xm) < 1.5 * P and abs(ym) < 1.5 * P:
            sign_x = +1.0 if xm >= 0 else -1.0
            dx_off = (+0.26 if abs(ym) < 0.8 * P else -0.26) * sign_x
            sy_dir = +1.0 if side == 'N' else -1.0
            j1 = (pin_xy[0], pin_xy[1] + sy_dir * 0.14)
            j2 = (pin_xy[0] + dx_off, pin_xy[1] + sy_dir * (0.14 + abs(dx_off)))
            jog = [j1, j2]
            p_start = j2
        elif side in ('E', 'W') and abs(ym) < 1.5 * P and abs(xm) < 1.5 * P:
            sign_y = +1.0 if ym >= 0 else -1.0
            dy_off = (+0.26 if abs(xm) < 0.8 * P else -0.26) * sign_y
            sx_dir = +1.0 if side == 'E' else -1.0
            j1 = (pin_xy[0] + sx_dir * 0.14, pin_xy[1])
            j2 = (pin_xy[0] + sx_dir * (0.14 + abs(dy_off)), pin_xy[1] + dy_off)
            jog = [j1, j2]
            p_start = j2
        else:
            p_start = pin_xy
        lead = 2.15 * P - abs(p_start[1]) if (side in ('N', 'S') and abs(p_start[1]) < 2.0 * P) else (
               2.15 * P - abs(p_start[0]) if (side in ('E', 'W') and abs(p_start[0]) < 2.0 * P) else 0.35)
        side_exits[side].append(dict(
            sort_key=p_start[0] if side in ('N', 'S') else p_start[1],
            name=f'fl_{na}_{nb}', start=(vfl_name, 'flux'),
            p_exit=p_start, pre=jog, lead=max(0.35, lead),
        ))

    for side, items in side_exits.items():
        items.sort(key=lambda it: it['sort_key'])
        pool = pad_pools[side]
        n_it = len(items)
        idxs = [int(round(k * (len(pool) - 1) / max(1, n_it - 1))) for k in range(n_it)] if n_it > 1 else [len(pool) // 2]
        for it, pi in zip(items, idxs):
            pad_id = pool[pi]
            pad_xy = tuple(np.asarray(design.components['LP_' + pad_id].pins['tie'].middle, float))
            pts = it['pre'] + _octilinear_fanout(it['p_exit'], pad_xy, side, lead_mm=it['lead'])
            planned_routes.append(dict(
                name=it['name'], priority=1,
                start=it['start'], end=('LP_' + pad_id, 'tie'),
                anchors=pts, cpw=CPW_CTL, fillet=FILLET_S17,
            ))

    # 8. Automatic AirbridgeSpan insertion at EVERY C-tier CPW-over-CPW crossing (0 on Purcell/readout!)
    polylines = []
    for r in planned_routes:
        p0 = tuple(np.asarray(design.components[r['start'][0]].pins[r['start'][1]].middle, float))
        p1 = tuple(np.asarray(design.components[r['end'][0]].pins[r['end'][1]].middle, float))
        polylines.append([p0] + [tuple(pt) for pt in r['anchors']] + [p1])

    raw_by_route_seg = {}
    for i in range(len(planned_routes)):
        for j in range(i + 1, len(planned_routes)):
            xs = _segment_intersections_general(polylines[i], polylines[j])
            for xc, yc, seg_i, seg_j, ang_i, ang_j, cross_ang in xs:
                pi, pj = planned_routes[i]['priority'], planned_routes[j]['priority']
                if pi == 0 and pj == 0:
                    continue
                if pi < pj:
                    o_idx, o_seg, o_ang = j, seg_j, ang_j
                elif pj < pi:
                    o_idx, o_seg, o_ang = i, seg_i, ang_i
                else:
                    ai0, ai1 = polylines[i][seg_i], polylines[i][seg_i + 1]
                    bj0, bj1 = polylines[j][seg_j], polylines[j][seg_j + 1]
                    mi = min(math.hypot(xc - ai0[0], yc - ai0[1]), math.hypot(xc - ai1[0], yc - ai1[1]))
                    mj = min(math.hypot(xc - bj0[0], yc - bj0[1]), math.hypot(xc - bj1[0], yc - bj1[1]))
                    if mj >= mi:
                        o_idx, o_seg, o_ang = j, seg_j, ang_j
                    else:
                        o_idx, o_seg, o_ang = i, seg_i, ang_i
                a0 = polylines[o_idx][o_seg]
                t_dist = math.hypot(xc - a0[0], yc - a0[1])
                raw_by_route_seg.setdefault((o_idx, o_seg), []).append((t_dist, xc, yc, o_ang, cross_ang))

    bridges_by_route = {i: [] for i in range(len(planned_routes))}
    ab_count = 0
    for (o_idx, o_seg), hits in raw_by_route_seg.items():
        hits.sort(key=lambda h: h[0])
        a0, a1 = polylines[o_idx][o_seg], polylines[o_idx][o_seg + 1]
        seg_len = math.hypot(a1[0] - a0[0], a1[1] - a0[1])
        ux = (a1[0] - a0[0]) / max(seg_len, 1e-9)
        uy = (a1[1] - a0[1]) / max(seg_len, 1e-9)
        # Non-chaining cluster guard: merge ONLY a tight pair < 92 um apart (NEVER chains into mega-bridges!)
        clusters = []
        for h in hits:
            if clusters and (h[0] - clusters[-1][0][0]) < 0.092:
                clusters[-1].append(h)
            else:
                clusters.append([h])
        for cl in clusters:
            t_lo, t_hi = cl[0][0], cl[-1][0]
            t_mid = 0.5 * (t_lo + t_hi)
            min_cross = min(h[4] for h in cl)
            base_span = 32.0 if min_cross > 70.0 else 54.0
            span_um = base_span + (t_hi - t_lo) * 1000.0
            tot_um = span_um + 34.0
            pad_half = 0.5 * tot_um / 1000.0
            t_mid = max(pad_half + 0.018, min(seg_len - pad_half - 0.018, t_mid))
            xc = a0[0] + ux * t_mid
            yc = a0[1] + uy * t_mid
            o_ang = cl[0][3]
            u_rel = 0.0 if min_cross > 70.0 else (90.0 - min_cross)
            ab_count += 1
            ab_name = f'ABX_{ab_count:03d}'
            AirbridgeSpan(design, ab_name, options=dict(
                **_pos(xc, yc), orientation=f'{o_ang - 90.0:.2f}',
                clear_span=f'{span_um:.1f}um', total_length=f'{tot_um:.1f}um',
                pad_width='10um', pad_length='8um', bridge_width='8um', pier_size='6um',
                under_angle=f'{u_rel:.1f}', under_keepout='16um',
                **CPW_CTL,
            ))
            bridges_by_route[o_idx].append((o_seg, xc, yc, ab_name, 'over_bot', 'over_top'))

    for idx, r in enumerate(planned_routes):
        _draw_overpass_with_bridges(
            design, r['name'], r['start'], r['end'], r['anchors'],
            bridges_by_route[idx], r['cpw'], fillet=r['fillet'],
        )

    _assert_ctier_zero_pad_collisions(design)
    _assert_no_acute_bends(design, min_interior_deg=88.0)
    _assert_zero_unbridged_cpw_crossings(design)
    _assert_max_airbridge_length(design, max_l20_um=115.0)
    n_gnd_ab = _place_ground_airbridges(design, spacing_mm=0.42)
    design._design14_summary = dict(
        tier='C-tier', revision=REVISION, die_mm=(CTIER_DIE_MM, CTIER_DIE_MM),
        n_purcell_filters=4, n_readout_resonators=len(res_report),
        n_xy_lines=17, n_flux_lines=len(vfl_list),
        n_ctier_airbridges=ab_count, n_ground_airbridges=n_gnd_ab,
        pad_collisions=0, unbridged_crossings=0,
    )
    return design


def _export_gds_with_bumps(design, gds_path: Path, bump_pts: list[tuple[float, float]]):
    gds_path.parent.mkdir(parents=True, exist_ok=True)
    gds = design.renderers.gds
    gds.options['no_cheese']['buffer'] = '50um'
    gds.options['cheese']['cheese_0_x'] = '0um'
    gds.options['cheese']['cheese_0_y'] = '0um'
    gds.options['cheese']['view_in_file'] = {'main': {1: False}}
    gds.options['no_cheese']['view_in_file'] = {'main': {1: False}}
    gds.options['short_segments_to_not_fillet'] = 'False'
    gds.export_to_gds(str(gds_path))

    # Inject shared Indium micro-bump UBM grid on Layer 30 (in mm, matching Qiskit Metal's 1e-3 GDS unit!)
    lib = gdspy.GdsLibrary(infile=str(gds_path))
    top = lib.top_level()[0]
    half_mm = 0.5 * BUMP_SIZE_MM
    for bx, by in bump_pts:
        top.add(gdspy.Rectangle((bx - half_mm, by - half_mm), (bx + half_mm, by + half_mm), layer=BUMP_LAYER))
    lib.write_gds(str(gds_path))


def _polys_from_gds(gds_path: Path):
    lib = gdspy.GdsLibrary(infile=str(gds_path))
    top = lib.top_level()[0]
    return top.get_polygons(by_spec=True)


def _draw_layer(ax, polys_by_spec, layer_id, facecolor, edgecolor='none', lw=0.0, alpha=1.0, zorder=2):
    for (lay, _dt), plist in polys_by_spec.items():
        if lay != layer_id:
            continue
        for pts in plist:
            arr = np.asarray(pts, float)
            ax.add_patch(MplPoly(arr, closed=True, facecolor=facecolor, edgecolor=edgecolor,
                                 linewidth=lw, alpha=alpha, zorder=zorder))


def render_figures(qtier_gds: Path, ctier_gds: Path, bump_pts: list[tuple[float, float]]):
    fig_dir = DESIGN_DIR / 'fig'
    fig_dir.mkdir(parents=True, exist_ok=True)
    qp = _polys_from_gds(qtier_gds)
    cp = _polys_from_gds(ctier_gds)

    # 1. Q-tier ONLY (8 x 8 mm — pure Qubits + Floating Couplers + Indium Bumps!)
    fig, ax = plt.subplots(figsize=(8.5, 8.5), dpi=220)
    ax.set_facecolor('#ffffff')
    _draw_layer(ax, qp, 1, facecolor='#3b5b7a', zorder=2)
    _draw_layer(ax, qp, BUMP_LAYER, facecolor='#00a896', alpha=0.85, zorder=4)
    ax.set_xlim(-4.2, 4.2); ax.set_ylim(-4.2, 4.2); ax.set_aspect('equal')
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    ax.set_title(f'Design14 {REVISION} Q-Tier (8x8 mm Top Chip): ONLY 17 Qubits + 24 Floating Couplers + {len(bump_pts)} In Bumps',
                 fontsize=10.5, fontweight='bold')
    fig.tight_layout(); fig.savefig(fig_dir / 'qtier.png', dpi=220); plt.close(fig)

    from matplotlib.patches import Patch as _Patch
    ctier_legend = [
        _Patch(facecolor='#f59e0b', edgecolor='#7c2d12', label='Crossover Airbridge (L20)'),
        _Patch(facecolor='#9333ea', edgecolor='#581c87', label='Ground Airbridge (L22)'),
        _Patch(facecolor='#10b981', edgecolor='#065f46', label='Indium Bump (L30)'),
    ]

    # 2. C-tier ONLY (16 x 16 mm — Interposer Control & Readout Chip)
    fig, ax = plt.subplots(figsize=(10.5, 10.5), dpi=220)
    ax.set_facecolor('#ffffff')
    _draw_layer(ax, cp, 1, facecolor='#6d7a89', zorder=2)
    _draw_layer(ax, cp, 22, facecolor='#9333ea', edgecolor='#581c87', lw=0.25, zorder=4)
    _draw_layer(ax, cp, 20, facecolor='#f59e0b', edgecolor='#7c2d12', lw=0.25, zorder=5)
    _draw_layer(ax, cp, BUMP_LAYER, facecolor='#10b981', alpha=0.85, zorder=4)
    ax.add_patch(MplRect((-4.0, -4.0), 8.0, 8.0, fill=False, edgecolor='#0284c7', lw=1.4, ls='--', zorder=6))
    ax.legend(handles=ctier_legend, loc='upper right', fontsize=8.5, framealpha=0.92)
    ax.set_xlim(-8.3, 8.3); ax.set_ylim(-8.3, 8.3); ax.set_aspect('equal')
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    ax.set_title(f'Design14 {REVISION} C-Tier (16x16 mm Interposer): 4 Purcell + 17 Vertical Readout + 17 XY + 24 Vertical Flux',
                 fontsize=10.5, fontweight='bold')
    fig.tight_layout(); fig.savefig(fig_dir / 'ctier.png', dpi=220); plt.close(fig)

    # 3. Aligned Flip-Chip Stack Overlay (chip.png)
    fig, ax = plt.subplots(figsize=(10.5, 10.5), dpi=220)
    ax.set_facecolor('#ffffff')
    _draw_layer(ax, cp, 1, facecolor='#8b98a8', alpha=0.85, zorder=2)
    _draw_layer(ax, cp, 22, facecolor='#9333ea', edgecolor='#581c87', lw=0.25, zorder=4)
    _draw_layer(ax, cp, 20, facecolor='#f59e0b', edgecolor='#7c2d12', lw=0.25, zorder=5)
    # Q-tier etched pockets highlighted in crisp blue-navy overlay
    _draw_layer(ax, qp, 1, facecolor='#1d4ed8', alpha=0.35, zorder=3)
    _draw_layer(ax, qp, BUMP_LAYER, facecolor='#10b981', alpha=0.90, zorder=6)
    ax.add_patch(MplRect((-4.0, -4.0), 8.0, 8.0, fill=False, edgecolor='#1d4ed8', lw=1.8, ls='--', zorder=7))
    ax.legend(handles=ctier_legend, loc='upper right', fontsize=8.5, framealpha=0.92)
    ax.set_xlim(-8.3, 8.3); ax.set_ylim(-8.3, 8.3); ax.set_aspect('equal')
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    ax.set_title(f'Design14 {REVISION} Flip-Chip Stack Overlay (8x8 mm Q-Tier on 16x16 mm C-Tier Interposer)',
                 fontsize=11, fontweight='bold')
    fig.tight_layout(); fig.savefig(fig_dir / 'chip.png', dpi=220); plt.close(fig)

    # 4. Core 2Q+1TC Flip-Chip Alignment Zoom (core.png: D5 - TC_D5_Z2 - Z2)
    fig, ax = plt.subplots(figsize=(9.5, 6.2), dpi=240)
    ax.set_facecolor('#f8fafc')
    # C-tier metal in soft slate, Q-tier metal in translucent sapphire blue so vertical alignment pops!
    _draw_layer(ax, cp, 1, facecolor='#cbd5e1', edgecolor='#475569', lw=0.35, alpha=0.90, zorder=2)
    _draw_layer(ax, cp, 22, facecolor='#9333ea', edgecolor='#581c87', lw=0.45, zorder=5)
    _draw_layer(ax, cp, 20, facecolor='#f59e0b', edgecolor='#7c2d12', lw=0.5, zorder=6)
    _draw_layer(ax, qp, 1, facecolor='#2563eb', edgecolor='#1e3a8a', lw=0.45, alpha=0.38, zorder=4)
    _draw_layer(ax, qp, BUMP_LAYER, facecolor='#10b981', edgecolor='#065f46', lw=0.3, alpha=0.95, zorder=7)
    ax.legend(handles=ctier_legend, loc='upper right', fontsize=8.0, framealpha=0.92)
    ax.set_xlim(-0.42, 1.62); ax.set_ylim(-0.68, 0.72); ax.set_aspect('equal')
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    ax.set_title('Design14 r1 Flip-Chip Alignment Zoom (Blue = Q-Tier Qubits/Coupler | Grey = C-Tier Vertical Readout/XY/Flux | Green = In Bumps)',
                 fontsize=9.2, fontweight='bold')
    fig.tight_layout(); fig.savefig(fig_dir / 'core.png', dpi=240); plt.close(fig)


if __name__ == '__main__':
    d_q = build_qtier(variant='sim')
    d_c = build_ctier(variant='sim')
    bumps = compute_shared_indium_bumps(d_q, d_c)
    gds_dir = DESIGN_DIR / 'gds'
    q_gds = gds_dir / f'Design14_{DATE_TAG}_{REVISION}_qtier_sim.gds'
    c_gds = gds_dir / f'Design14_{DATE_TAG}_{REVISION}_ctier_sim.gds'
    _export_gds_with_bumps(d_q, q_gds, bumps)
    _export_gds_with_bumps(d_c, c_gds, bumps)
    render_figures(q_gds, c_gds, bumps)
    print('DESIGN14_QTIER_OK:', d_q._design14_summary)
    print('DESIGN14_CTIER_OK:', d_c._design14_summary)
    print('INDIUM_BUMPS_OK:', len(bumps))
