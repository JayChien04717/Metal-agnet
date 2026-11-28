# -*- coding: utf-8 -*-
r"""Stage 1 draft layout entry point for Design13_ScalableFloatingCoupler2Q (r9).

14 x 14 mm octilinear (135-deg obtuse-bend / 45-deg diagonal) rotated d=3
surface code (Surface-17: 17 circular wide-gap fixed-frequency qubits, 24
floating tunable couplers) on a 1.00 mm square-lattice pitch.

Key physics & layout improvements over r3:
1. Upgraded Qubit (ScalablePadQubit v2): Pure circular 256 um island with wide
   46 um uniform pocket gap (minimal perimeter-to-area ratio, zero corners) and
   smooth radial 45-deg arc readout pad (eliminates sharp tabs/claws).
2. Exact lambda/4 Readout Resonators (L_total = 4.17-4.36 mm for 6.85-7.17 GHz):
   The 4 continuous layer-1 shared-bus lambda/2 Purcell filters sit diagonally
   at +-45 deg IMMEDIATELY outside the 4 diamond edges of Surface-17 (only
   1.0-2.0 mm from the qubits!), so EVERY readout resonator's TOTAL filleted
   length from qubit readout pad to ShortToGround equals exact lambda/4!
   Only the broadband 50-ohm feedlines extend to the perimeter bond pads.
   ZERO airbridges on any Purcell filter or readout resonator meander!
3. 135-deg Obtuse-Bend (45-deg Diagonal) Control Fanout:
   XY and TC flux lines use 135-deg obtuse bends (lower reflection, shorter
   paths, fewer crossings than 90-deg staircases) on a compact 14 x 14 mm die.
4. Dual-mode entry point:
   - mode='surface17' (default): 14 x 14 mm octilinear Surface-17.
   - mode='unit_cell': 10 x 10 mm 8-port 2Q + 1TC EM unit cell (identical
     circular ScalablePadQubit + FloatingTunableCoupler + PurcellReadout).

Stage 1 draft only: EM simulations have NOT been run yet.
"""
from __future__ import annotations

from collections import OrderedDict
import math
import re
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'code' / 'general'))
import project_paths
project_paths.configure(__file__)
DESIGN = project_paths.design_dir(__file__)

COMPONENTS_DIR = Path(__file__).resolve().parent / 'components'
if str(COMPONENTS_DIR) not in sys.path:
    sys.path.insert(0, str(COMPONENTS_DIR))

import numpy as np
from qiskit_metal import Dict, designs
from qiskit_metal.qlibrary.core import QRoute
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight
from shapely.geometry import LineString, MultiPoint, Point
from shapely.ops import unary_union

from airbridge_crossover import AirbridgeCrossover, AirbridgeGround, AirbridgeSpan
from fluxline import T_FLUX_OPTIONS
from bondpad import BondPad, add_bond_pad
from check_keepouts import check_keepouts
from enclosing_cap import EnclosingCap
from floating_tunable_coupler import FloatingTunableCoupler
from purcell_readout import PurcellReadout
from scalable_pad_qubit import ScalablePadQubit
from tap_junction import TapJunction
from wafer_pcbport import wafer_pcb_ports

REVISION = 'r9'
CPW_RO = dict(trace_width='10um', trace_gap='6um')
CPW_CTL = dict(trace_width='10um', trace_gap='6um')
FILLET_S17 = '55um'
FILLET_UC = '99um'
PITCH_MM = 2.00
# r8: 20 x 20 mm die. Everything that depends on the die size is derived from these.
DIE_MM = 30.0   # r9: back to 30 x 30 mm (r8 tried 20 mm)
DIE_HALF = 0.5 * DIE_MM
PAD_INSET_MM = 0.65              # launch-pad centre inset from the die edge (ground encloses every port)
PAD_SPAN_MM = 11.05              # outermost pad coordinate along an edge (r7 pad ring on the 30 mm die)
U_PURCELL_MM = 5.60              # diagonal Purcell filters: outermost qubits sit at u = 4.24 mm (2 mm pitch)
PF_XC_MM = 6.90                  # Purcell diagonal corner coordinate (outer legs at PF_XC_MM + 0.28)
V_PH_MM_GHZ = 119.5  # lambda/4 (mm) = V_PH_MM_GHZ / (4 * f_GHz)
INV_SQRT2 = 1.0 / math.sqrt(2.0)

QUBIT_OPT = dict(
    pad_diameter='280um',      # Ec 193 MHz (small-model Q3D, 2026-11-28)
    pocket_gap='46um',
    jj_gap='36um',             # vertical Manhattan JJ at a free diagonal corner (180 deg pattern on N corners)
    ro_key_pad_radius='13um',  # keyhole readout: C_qr ~10 fF, 2chi ~1.4 MHz, readout<->coupler parasitic 1.7 %
    ro_key_gap='20um',
    cpw_width='10um',
    cpw_gap='6um',
    cpw_length='42um',
)

COUPLER_OPT = dict(
    # r9 (Q3D + HFSS 2026-11-28): C-arcs concentric with the qubits, 10/6 um necks, rounded-rectangle hub with a
    # rail SQUID (island B's thin rail under island A's flat bottom edge, two Manhattan JJs between them -> loop
    # parallel to a straight T flux line). JJs at multiples of 90 deg (Manhattan-compliant). Small-model results:
    # Ec(TC) 199 MHz, g_qc 80 MHz, own M 0.73 pH, coupler-to-coupler crosstalk <= 0.2 %.
    qubit_distance=f'{0.5 * PITCH_MM * 1000.0:.1f}um',
    qubit_pad_radius='140um',
    c_sep='26um',
    c_width='24um',
    c_span_deg='28',
    c_gap='18um',
    neck_width='10um',
    neck_gap='6um',
    rect_width='230um',
    rect_height='180um',
    rect_ab_gap='6um',
    squid_gap='20um',
    squid_spacing='60um',
    flux_line_options=Dict(make_fbl='True', mirror='False', **T_FLUX_OPTIONS),   # T placed automatically on the loop
    Lj='6.5nH',
)


def _xy_corner(qname: str) -> str:
    """Diagonal corner used by each qubit's XY stub (mirrors the XY placement in _build_surface17)."""
    fixed = {'D5': 'NW', 'X1': 'NE', 'X2': 'SW', 'Z1': 'NW', 'Z2': 'SE'}
    if qname in fixed:
        return fixed[qname]
    side = SURFACE17_QUBITS[qname]['xy_side']
    return {'N': 'NW' if qname == 'Z3' else 'NE', 'S': 'SE' if qname == 'Z4' else 'SW',
            'W': 'SW' if qname == 'X3' else 'NW', 'E': 'NE' if qname == 'X4' else 'SE'}[side]


def _jj_corner(qname: str) -> str:
    """Free diagonal corner (not readout / XY) for the vertical Manhattan JJ; N corners use the 180-deg pattern."""
    taken = {SURFACE17_QUBITS[qname]['ro_corner'], _xy_corner(qname)}
    for c in ('SW', 'SE', 'NW', 'NE'):
        if c not in taken:
            return c
    raise AssertionError(f'{qname}: both NE and SW are taken by readout / XY')

SURFACE17_QUBITS = OrderedDict([
    # Data qubits (9) — exact C4-symmetric open 45-deg diagonal readout corners!
    ('D1', dict(role='data', i=0,  j=-2, f_GHz=4.65, Lj='10.8nH', ro_quad='SE', ro_corner='SE', xy_side='S', f_ro_GHz=7.05)),
    ('D2', dict(role='data', i=1,  j=-1, f_GHz=4.60, Lj='11.0nH', ro_quad='SE', ro_corner='SE', xy_side='S', f_ro_GHz=6.95)),
    ('D3', dict(role='data', i=2,  j=0,  f_GHz=4.65, Lj='10.8nH', ro_quad='NE', ro_corner='NE', xy_side='E', f_ro_GHz=6.85)),
    ('D4', dict(role='data', i=-1, j=-1, f_GHz=4.60, Lj='11.0nH', ro_quad='SW', ro_corner='SW', xy_side='W', f_ro_GHz=6.95)),
    ('D5', dict(role='data', i=0,  j=0,  f_GHz=4.65, Lj='10.8nH', ro_quad='NE', ro_corner='NE', xy_side='W', f_ro_GHz=7.17)),
    ('D6', dict(role='data', i=1,  j=1,  f_GHz=4.60, Lj='11.0nH', ro_quad='NE', ro_corner='NE', xy_side='E', f_ro_GHz=7.01)),
    ('D7', dict(role='data', i=-2, j=0,  f_GHz=4.65, Lj='10.8nH', ro_quad='SW', ro_corner='SW', xy_side='W', f_ro_GHz=6.85)),
    ('D8', dict(role='data', i=-1, j=1,  f_GHz=4.60, Lj='11.0nH', ro_quad='NW', ro_corner='NW', xy_side='N', f_ro_GHz=6.95)),
    ('D9', dict(role='data', i=0,  j=2,  f_GHz=4.65, Lj='10.8nH', ro_quad='NW', ro_corner='NW', xy_side='N', f_ro_GHz=7.05)),
    # X-ancillas (4)
    ('X1', dict(role='X',    i=0,  j=1,  f_GHz=4.50, Lj='11.5nH', ro_quad='NW', ro_corner='NW', xy_side='N', f_ro_GHz=7.15)),
    ('X2', dict(role='X',    i=0,  j=-1, f_GHz=4.50, Lj='11.5nH', ro_quad='SE', ro_corner='SE', xy_side='S', f_ro_GHz=7.15)),
    ('X3', dict(role='X',    i=-2, j=1,  f_GHz=4.50, Lj='11.5nH', ro_quad='NW', ro_corner='NW', xy_side='W', f_ro_GHz=6.85)),
    ('X4', dict(role='X',    i=2,  j=-1, f_GHz=4.50, Lj='11.5nH', ro_quad='SE', ro_corner='SE', xy_side='E', f_ro_GHz=6.85)),
    # Z-ancillas (4)
    ('Z1', dict(role='Z',    i=-1, j=0,  f_GHz=4.50, Lj='11.5nH', ro_quad='SW', ro_corner='SW', xy_side='W', f_ro_GHz=7.15)),
    ('Z2', dict(role='Z',    i=1,  j=0,  f_GHz=4.50, Lj='11.5nH', ro_quad='NE', ro_corner='NE', xy_side='E', f_ro_GHz=7.09)),
    ('Z3', dict(role='Z',    i=1,  j=2,  f_GHz=4.50, Lj='11.5nH', ro_quad='NE', ro_corner='NE', xy_side='N', f_ro_GHz=6.93)),
    ('Z4', dict(role='Z',    i=-1, j=-2, f_GHz=4.50, Lj='11.5nH', ro_quad='SW', ro_corner='SW', xy_side='S', f_ro_GHz=7.05)),
])

ENC_OPT = dict(plate_width='220um', plate_length='30um', gap='6um',
               enclosing_width='10um', enclosing_gap='6um', lead='40um')
PURCELL_CHANNELS = OrderedDict([
    ('Q1', dict(side=-1, f_GHz=7.05, L_filter=3.80, L_res=3.95, head_gap_um=48.0, head_height_um=40.0, cp_plate_um=220.0)),
    ('Q2', dict(side=+1, f_GHz=6.85, L_filter=3.95, L_res=4.10, head_gap_um=48.0, head_height_um=40.0, cp_plate_um=160.0)),
])
HAMMER_LINE_SPACING_UM = 200.0
HAMMER_HEAD_WIDTH_UM = 40.0
PORT_FEED = dict(y_top=3.55, x_side=2.70)


class _PolyRoute(QRoute):
    """Deterministic filleted CPW route supporting arbitrary 45/90/135-deg anchor polylines."""

    default_options = Dict(anchors=OrderedDict())

    def make(self):
        self.set_pin('start')
        self.set_pin('end')
        self.set_lead('start')
        self.set_lead('end')
        pts = [np.asarray(v, float) for v in self.options.anchors.values()]
        self.intermediate_pts = np.asarray(pts, float) if pts else np.empty((0, 2), float)
        raw = [np.asarray(p, float) for p in self.get_points()]
        clean = [raw[0]]
        for p in raw[1:-1]:
            if np.linalg.norm(p - clean[-1]) > 0.008 and np.linalg.norm(p - raw[-1]) > 0.008:
                clean.append(p)
        clean.append(raw[-1])
        self.make_elements(np.asarray(clean, float))


def _pos(x_mm: float, y_mm: float) -> dict:
    return dict(pos_x=f'{x_mm:.6f}mm', pos_y=f'{y_mm:.6f}mm')


def _route(design, name, start, end, pts, cpw, fillet=FILLET_UC):
    p_s = np.asarray(design.components[start[0]].pins[start[1]].middle, float)
    p_e = np.asarray(design.components[end[0]].pins[end[1]].middle, float)
    clean_pts = []
    for pt in pts:
        arr = np.array(pt, float)
        if not clean_pts or np.linalg.norm(arr - clean_pts[-1]) > 1e-4:
            clean_pts.append(arr)
    # Remove any nearly-collinear redundant anchor (> 178.5 deg) so gdstk/fillet never emits a miter spike!
    full = [p_s] + clean_pts + [p_e]
    pruned = []
    for k in range(1, len(full) - 1):
        v1 = full[k - 1] - full[k]
        v2 = full[k + 1] - full[k]
        n1, n2 = float(np.linalg.norm(v1)), float(np.linalg.norm(v2))
        if n1 > 1e-4 and n2 > 1e-4:
            ang = math.degrees(math.acos(np.clip(float(np.dot(v1, v2) / (n1 * n2)), -1.0, 1.0)))
            if ang > 178.5:
                continue
        pruned.append(full[k])
    clean_pts = pruned
    if not clean_pts:
        return RouteStraight(design, name, Dict(
            pin_inputs=Dict(
                start_pin=Dict(component=start[0], pin=start[1]),
                end_pin=Dict(component=end[0], pin=end[1]),
            ),
            **cpw,
        ))
    return _PolyRoute(design, name, Dict(
        pin_inputs=Dict(
            start_pin=Dict(component=start[0], pin=start[1]),
            end_pin=Dict(component=end[0], pin=end[1]),
        ),
        anchors=OrderedDict((i, arr) for i, arr in enumerate(clean_pts)),
        fillet=fillet, **cpw,
    ))


def _surface17_edges():
    edges = []
    names = list(SURFACE17_QUBITS.keys())
    for idx_a, na in enumerate(names):
        qa = SURFACE17_QUBITS[na]
        for nb in names[idx_a + 1:]:
            qb = SURFACE17_QUBITS[nb]
            if abs(qa['i'] - qb['i']) + abs(qa['j'] - qb['j']) == 1:
                edges.append((na, nb) if (qa['i'], qa['j']) <= (qb['i'], qb['j']) else (nb, na))
    return edges


def _perimeter_pads_14mm():
    """56 symmetric launch pads (14 per side) on the DIE_MM die, inset PAD_INSET_MM so solid ground encloses every port."""
    coords = np.linspace(-PAD_SPAN_MM, PAD_SPAN_MM, 14)
    edge = DIE_HALF - PAD_INSET_MM
    pads = OrderedDict()
    for k, x in enumerate(coords):
        pads[f'N{k:02d}'] = dict(pos_x=f'{x:.4f}mm', pos_y=f'{+edge:.4f}mm', orientation='-90')
    for k, x in enumerate(coords):
        pads[f'S{k:02d}'] = dict(pos_x=f'{x:.4f}mm', pos_y=f'{-edge:.4f}mm', orientation='90')
    for k, y in enumerate(coords):
        pads[f'W{k:02d}'] = dict(pos_x=f'{-edge:.4f}mm', pos_y=f'{y:.4f}mm', orientation='0')
    for k, y in enumerate(coords):
        pads[f'E{k:02d}'] = dict(pos_x=f'{+edge:.4f}mm', pos_y=f'{y:.4f}mm', orientation='180')
    return pads


def _coupler_orientation(na: str, nb: str, xm: float, ym: float, is_horiz: bool) -> str:
    pair = frozenset((na, nb))
    overrides = {
        frozenset(('Z2', 'D3')): '180',
        frozenset(('X1', 'D9')): '270',
        frozenset(('D1', 'X2')): '90',
        frozenset(('D7', 'Z1')): '0',
        frozenset(('D5', 'X1')): '90',
        frozenset(('X2', 'D5')): '90',
        frozenset(('Z1', 'D5')): '0',
        frozenset(('D5', 'Z2')): '0',
    }
    if pair in overrides:
        return overrides[pair]
    if is_horiz:
        return '180' if ym >= 0.5 else '0'
    return '270' if xm <= -0.5 else '90'


def _poly_filleted_length(verts: list[tuple[float, float]], fillet_mm: float) -> float:
    """Exact filleted centre-line length for any polyline (including 90-deg and 135-deg bends!)."""
    total = 0.0
    for k in range(len(verts) - 1):
        total += math.hypot(verts[k + 1][0] - verts[k][0], verts[k + 1][1] - verts[k][1])
    for k in range(1, len(verts) - 1):
        v1 = np.array(verts[k]) - np.array(verts[k - 1])
        v2 = np.array(verts[k + 1]) - np.array(verts[k])
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 < 1e-6 or n2 < 1e-6:
            continue
        cos_t = np.clip(float(np.dot(v1, v2) / (n1 * n2)), -1.0, 1.0)
        turn = math.acos(cos_t)  # exterior turn angle (pi/2 for 90 deg, pi/4 for 135 deg)
        if turn > 1e-3:
            total -= fillet_mm * (2.0 * math.tan(0.5 * turn) - turn)
    return total


def _diag_uv_to_xy(u: float, v: float, sx: int, sy: int) -> tuple[float, float]:
    """Map quadrant diagonal frame (u = radial outward along (sx, sy), v = tangent along (-sy, sx)) to (x, y)."""
    # e_u = (sx, sy)/sqrt(2), e_v = (-sy, sx)/sqrt(2) -> right-handed when sx*sy > 0
    return (
        (sx * u - sy * v) * INV_SQRT2,
        (sy * u + sx * v) * INV_SQRT2,
    )


def _xy_to_diag_uv(x: float, y: float, sx: int, sy: int) -> tuple[float, float]:
    return (
        (sx * x + sy * y) * INV_SQRT2,
        (-sy * x + sx * y) * INV_SQRT2,
    )


def _build_diag_lambda4_resonator(
    p_pin: tuple[float, float], sx: int, sy: int,
    v_slot: float, u_purcell: float, target_total_mm: float,
    fillet_mm: float, pitch_mm: float = 0.35, coup_leg_mm: float = 0.28, early_turn: bool = False, u_clear: float = 0.0,
):
    """Build an EXACT total-length lambda/4 DRC-B7 readout resonator in the 45-deg diagonal quadrant!

    - Starts right at the qubit's 45-deg radial readout pin `p_pin`
    - Runs radially outward into the 45-deg meander bay at tangent `v_slot` (just outside the diamond!)
    - Solves meander depth `h` so the ENTIRE filleted path (pin -> lead -> meander -> coupling leg)
      equals `target_total_mm` (4.17-4.36 mm for 6.85-7.17 GHz) to < 5 um!
    - Coupling leg runs parallel to the 45-deg Purcell filter at `u_coup = u_purcell - 0.032 mm` (10 um GND strip).
    """
    u_pin, v_pin = _xy_to_diag_uv(p_pin[0], p_pin[1], sx, sy)
    u_coup = u_purcell - 0.032

    def _uv_verts(h: float):
        # Shoot 100% STRAIGHT radially at 45 deg past ALL control lines all the way to u_lo (outer meander bay)!
        u_lo = u_coup - h
        u_hi = u_coup - 0.22
        # early_turn: shift sideways right after leaving the qubit, so the lead never runs through the v-range of
        # a neighbouring resonator's meander (dense readout banks on small dies)
        u_turn = max(u_pin + 0.25, u_clear) if early_turn else max(u_pin + 0.25, u_lo)
        v0 = v_slot
        v1 = v_slot + pitch_mm
        v2 = v_slot + 2.0 * pitch_mm
        v_end = v2 + coup_leg_mm
        pts_uv = [(u_pin, v_pin)]
        if abs(v0 - v_pin) > 0.14:
            pts_uv += [(u_turn, v_pin), (u_turn, v0)]
        else:
            v0 = v_pin
            v1 = v0 + pitch_mm
            v2 = v0 + 2.0 * pitch_mm
            v_end = v2 + coup_leg_mm
        pts_uv += [
            (u_hi, v0),
            (u_hi, v1),
            (u_lo, v1),
            (u_lo, v2),
            (u_coup, v2),
            (u_coup, v_end),
        ]
        return [_diag_uv_to_xy(u, v, sx, sy) for u, v in pts_uv]

    # Solve meander depth h so TOTAL filleted length (from qubit pin to ShortToGround!) == target_total_mm
    lo, hi = 0.35, 1.35
    for _ in range(36):
        mid = 0.5 * (lo + hi)
        if _poly_filleted_length(_uv_verts(mid), fillet_mm) < target_total_mm:
            lo = mid
        else:
            hi = mid
    verts = _uv_verts(0.5 * (lo + hi))
    return verts[1:-1], verts[-1], _poly_filleted_length(verts, fillet_mm)


def _segment_intersections_general(poly_a: list[tuple[float, float]], poly_b: list[tuple[float, float]]):
    """Return list of (xc, yc, ia, ib, angle_a_deg, angle_b_deg, cross_angle_deg) for interior segment crossings."""
    from shapely.geometry import Point as _Pt
    out = []
    for ia in range(len(poly_a) - 1):
        a0, a1 = poly_a[ia], poly_a[ia + 1]
        la = math.hypot(a1[0] - a0[0], a1[1] - a0[1])
        if la < 0.008:
            continue
        sa = LineString([a0, a1])
        ang_a = math.degrees(math.atan2(a1[1] - a0[1], a1[0] - a0[0]))
        for ib in range(len(poly_b) - 1):
            b0, b1 = poly_b[ib], poly_b[ib + 1]
            lb = math.hypot(b1[0] - b0[0], b1[1] - b0[1])
            if lb < 0.008:
                continue
            sb = LineString([b0, b1])
            inter = sa.intersection(sb)
            if not isinstance(inter, _Pt):
                continue
            xc, yc = float(inter.x), float(inter.y)
            # Only skip literal shared pin endpoints (< 15 um from BOTH route endpoints)
            if (min(math.hypot(xc - a0[0], yc - a0[1]), math.hypot(xc - a1[0], yc - a1[1])) < 0.015 and
                    min(math.hypot(xc - b0[0], yc - b0[1]), math.hypot(xc - b1[0], yc - b1[1])) < 0.015):
                continue
            ang_b = math.degrees(math.atan2(b1[1] - b0[1], b1[0] - b0[0]))
            diff = abs((ang_a - ang_b + 90.0) % 180.0 - 90.0)
            if diff < 15.0:  # nearly parallel, not a transverse crossing
                continue
            out.append((xc, yc, ia, ib, ang_a, ang_b, diff))
    return out


def _draw_overpass_with_bridges(design, base_name, start_pin, end_pin, anchors, bridges_on_me, cpw, fillet):
    """Route a polyline, splitting cleanly at AirbridgeSpan pins (over_bot/over_top) so layer 1 has zero shorts."""
    p_start = tuple(np.asarray(design.components[start_pin[0]].pins[start_pin[1]].middle, float))
    p_end = tuple(np.asarray(design.components[end_pin[0]].pins[end_pin[1]].middle, float))
    verts = [p_start] + [tuple(pt) for pt in anchors] + [p_end]

    def _key(c):
        seg_idx, xc, yc, *_ = c
        vx, vy = verts[seg_idx]
        return (seg_idx, math.hypot(xc - vx, yc - vy))

    def _safe_fillet(p_a, p_b, anc_list):
        """Clamp fillet on short sub-routes while ALWAYS preserving smooth circular arc bends (never 0 um sharp corners!)."""
        f_req = float(design.parse_value(fillet))
        pts_all = [p_a] + [tuple(p) for p in anc_list] + [p_b]
        if len(pts_all) < 3:
            return fillet
        min_seg = min(math.hypot(pts_all[k + 1][0] - pts_all[k][0], pts_all[k + 1][1] - pts_all[k][1])
                      for k in range(len(pts_all) - 1))
        if min_seg < 3.0 * f_req:
            max_r_um = 0.28 * min_seg * 1000.0
            return f'{max(4.0, min(f_req * 1000.0, max_r_um)):.1f}um'
        return fillet

    if not bridges_on_me:
        _route(design, base_name, start_pin, end_pin, anchors, cpw, fillet=_safe_fillet(p_start, p_end, anchors))
        return

    ordered = sorted(bridges_on_me, key=_key)
    cur_pin = start_pin
    cur_seg_start = 0
    sub_idx = 0
    for seg_idx, xc, yc, ab_name, pin_in, pin_out in ordered:
        sub_anchors = [verts[k] for k in range(cur_seg_start + 1, seg_idx + 1)]
        pa = tuple(np.asarray(design.components[cur_pin[0]].pins[cur_pin[1]].middle, float))
        pb = tuple(np.asarray(design.components[ab_name].pins[pin_in].middle, float))
        _route(design, f'{base_name}_s{sub_idx}', cur_pin, (ab_name, pin_in),
               sub_anchors, cpw, fillet=_safe_fillet(pa, pb, sub_anchors))
        cur_pin = (ab_name, pin_out)
        cur_seg_start = seg_idx
        sub_idx += 1

    tail_anchors = [verts[k] for k in range(cur_seg_start + 1, len(verts) - 1)]
    pa = tuple(np.asarray(design.components[cur_pin[0]].pins[cur_pin[1]].middle, float))
    _route(design, f'{base_name}_s{sub_idx}', cur_pin, end_pin,
           tail_anchors, cpw, fillet=_safe_fillet(pa, p_end, tail_anchors))


def _place_ground_airbridges(design, spacing_mm: float = 0.42, corner_clear_mm: float = 0.075,
                             end_clear_mm: float = 0.10, bridge_clear_mm: float = 0.095, step_mm: float = 0.01):
    """Place ground-stitching AirbridgeGround (L22/L23) along every 50-ohm control CPW (skill/airbridge.md §4).

    Arc-length based: candidate sites every `step_mm` along each drawn sub-route are valid when they keep
    `corner_clear_mm` from bends, `end_clear_mm` from sub-route ends (crossover-bridge pins / component pins),
    `bridge_clear_mm` from any L20 crossover bridge, and both piers land on ground (outside every component
    pocket / pad and every other CPW trench). A greedy pass then takes the farthest valid site within
    `spacing_mm` of the previous bridge, so a bridge is only missing where no valid site exists at all
    (reported by _assert_ground_bridge_coverage)."""
    from shapely.geometry import Point as _Pt
    from shapely.ops import unary_union as _union

    id_to_name = {c.id: n for n, c in design.components.items()}
    poly_tbl = design.qgeometry.tables['poly']
    path_tbl = design.qgeometry.tables['path']

    l20, comp = [], []
    for row in poly_tbl.itertuples():
        lay = int(row.layer)
        if lay == 20:
            l20.append(row.geometry.buffer(bridge_clear_mm))
        elif lay == 1:
            comp.append(row.geometry.buffer(0.050))
    ko = _union(l20 + comp) if (l20 or comp) else _Pt(999, 999).buffer(0.001)

    cpws = []
    for row in path_tbl.itertuples():
        if not row.subtract and int(row.layer) == 1:
            cpws.append((id_to_name.get(row.component, ''), row.geometry))

    gnd_count = 0
    for rname, geom in cpws:
        if not rname.startswith(('xy_', 'fl_', 'pf_feed_')):
            continue
        others = _union([g.buffer(0.016) for nm, g in cpws if nm != rname]) if len(cpws) > 1 else None
        coords = list(geom.coords)
        L = geom.length
        if L < 2.0 * end_clear_mm:
            continue
        # arc positions of interior vertices
        s_acc, v_s = 0.0, []
        for k in range(1, len(coords)):
            s_acc += math.hypot(coords[k][0] - coords[k - 1][0], coords[k][1] - coords[k - 1][1])
            if k < len(coords) - 1:
                v_s.append(s_acc)
        cand = []
        for t in np.arange(end_clear_mm, L - end_clear_mm + 1e-9, step_mm):
            if any(abs(t - vs) < corner_clear_mm for vs in v_s):
                continue
            p = geom.interpolate(float(t))
            q = geom.interpolate(float(min(L, t + 0.005)))
            q0 = geom.interpolate(float(max(0.0, t - 0.005)))
            dx, dy = q.x - q0.x, q.y - q0.y
            nn = math.hypot(dx, dy)
            if nn < 1e-9:
                continue
            ux, uy = dx / nn, dy / nn
            nx, ny = -uy, ux
            pa = _Pt(p.x + nx * 0.0175, p.y + ny * 0.0175)
            pb = _Pt(p.x - nx * 0.0175, p.y - ny * 0.0175)
            if ko.contains(p) or ko.contains(pa) or ko.contains(pb):
                continue
            if others is not None and (others.contains(pa) or others.contains(pb) or others.contains(p)):
                continue
            cand.append((float(t), p.x, p.y, math.degrees(math.atan2(uy, ux))))
        if not cand:
            continue
        chosen, cur = [], 0.0
        while True:
            if L - cur <= spacing_mm:
                break
            win = [c for c in cand if cur + 0.12 <= c[0] <= cur + spacing_mm]
            if not win:
                win = [c for c in cand if c[0] > cur + 0.12][:1]
                if not win:
                    break
            c = win[-1]
            chosen.append(c)
            cur = c[0]
        for t, cx, cy, ang in chosen:
            gnd_count += 1
            AirbridgeGround(design, f'ABG_{gnd_count:03d}', options=dict(
                **_pos(cx, cy), orientation=f'{ang:.2f}',
                span_length='35um', bridge_width='10um', pier_size='5um',
                bridge_layer='22', pier_layer='23',
            ))
    return gnd_count


def _assert_max_airbridge_length(design, max_l20_um: float = 95.0):
    """Hard Stage-1 DRC: verify NO crossover airbridge (L20) exceeds fab arch-length limit."""
    id_to_name = {c.id: n for n, c in design.components.items()}
    bad = []
    for row in design.qgeometry.tables['poly'].itertuples():
        if int(row.layer) != 20:
            continue
        pts = np.asarray(row.geometry.exterior.coords)
        L_um = max(math.hypot(pts[k + 1][0] - pts[k][0], pts[k + 1][1] - pts[k][1])
                   for k in range(len(pts) - 1)) * 1000.0
        if L_um > max_l20_um:
            c = pts.mean(axis=0)
            bad.append((id_to_name.get(row.component, ''), round(L_um, 1), round(float(c[0]), 2), round(float(c[1]), 2)))
    if bad:
        bad.sort(key=lambda x: -x[1])
        sample = ', '.join(f'{n}({l}um@{x},{y})' for n, l, x, y in bad[:8])
        raise AssertionError(f'Found {len(bad)} oversized L20 airbridges (>{max_l20_um:.0f} um): {sample}')


def _assert_zero_unbridged_cpw_crossings(design):
    """Hard Stage-1 DRC: verify ZERO unbridged Layer-1 CPW-to-CPW crossings anywhere on the chip."""
    id_to_name = {c.id: n for n, c in design.components.items()}
    path_tbl = design.qgeometry.tables['path']
    rows = [r for r in path_tbl.itertuples() if not r.subtract and int(r.layer) == 1]
    bad = []
    for i in range(len(rows)):
        gi = rows[i].geometry
        ni = id_to_name.get(rows[i].component, str(rows[i].component))
        ends_i = (gi.coords[0], gi.coords[-1])
        for j in range(i + 1, len(rows)):
            nj = id_to_name.get(rows[j].component, str(rows[j].component))
            if ni == nj:
                continue
            gj = rows[j].geometry
            if not gi.intersects(gj):
                continue
            inter = gi.intersection(gj)
            if inter.geom_type == 'Point':
                p = (float(inter.x), float(inter.y))
                ends_j = (gj.coords[0], gj.coords[-1])
                if (any(math.hypot(p[0] - e[0], p[1] - e[1]) < 0.015 for e in ends_i) and
                        any(math.hypot(p[0] - e[0], p[1] - e[1]) < 0.015 for e in ends_j)):
                    continue
                bad.append((ni, nj, round(p[0], 3), round(p[1], 3)))
            else:
                bad.append((ni, nj, inter.geom_type, round(float(inter.centroid.x), 3), round(float(inter.centroid.y), 3)))
    if bad:
        sample = ', '.join(str(b) for b in bad[:10])
        raise AssertionError(f'Found {len(bad)} unbridged Layer-1 CPW crossings: {sample}')


def _assert_no_acute_bends(design, min_interior_deg: float = 88.0):
    """Hard Stage-1 DRC: verify every CPW route bend has interior angle >= 90 deg (no acute <90 corners)."""
    id_to_name = {c.id: n for n, c in design.components.items()}
    path_tbl = design.qgeometry.tables['path']
    bad = []
    for row in path_tbl.itertuples():
        if row.subtract or int(row.layer) != 1:
            continue
        coords = list(row.geometry.coords)
        for k in range(1, len(coords) - 1):
            v1 = np.array(coords[k - 1]) - np.array(coords[k])
            v2 = np.array(coords[k + 1]) - np.array(coords[k])
            n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
            if n1 < 1e-3 or n2 < 1e-3:
                continue
            cos_a = np.clip(float(np.dot(v1, v2) / (n1 * n2)), -1.0, 1.0)
            interior_deg = math.degrees(math.acos(cos_a))
            if interior_deg < min_interior_deg:
                bad.append((id_to_name.get(row.component, str(row.component)), round(interior_deg, 1), coords[k]))
    if bad:
        sample = ', '.join(f'{n}({a}deg@{p[0]:.2f},{p[1]:.2f})' for n, a, p in bad[:6])
        raise AssertionError(f'Found {len(bad)} acute (<90 deg) CPW bends: {sample}')


def _assert_zero_pad_collisions(design):
    """Hard Stage-1 DRC: verify no CPW route trace intersects any qubit, coupler, or EnclosingCap."""
    id_to_name = {c.id: n for n, c in design.components.items()}
    poly_tbl = design.qgeometry.tables['poly']
    path_tbl = design.qgeometry.tables['path']

    islands = []
    cap_boxes = []
    for row in poly_tbl.itertuples():
        cname = id_to_name.get(row.component, '')
        if int(row.layer) != 1:
            continue
        if not row.subtract and (cname in SURFACE17_QUBITS or cname.startswith('TC_')):
            if row.name in ('pad', 'pad_l', 'pad_r'):
                islands.append((cname, row.name, row.geometry))
        elif row.subtract and cname.startswith('PF_c'):
            cap_boxes.append((cname, row.geometry.buffer(0.02)))  # 20 um clearance around EnclosingCap pocket

    collisions = []
    for row in path_tbl.itertuples():
        if row.subtract or int(row.layer) != 1:
            continue
        rname = id_to_name.get(row.component, str(row.component))
        trace_poly = row.geometry.buffer(float(row.width) / 2.0, cap_style=2)
        for cname, pname, geom in islands:
            if trace_poly.intersects(geom):
                collisions.append((rname, cname, pname))
        # Check EnclosingCap collisions for every route NOT directly connected to that capacitor
        comp = design.components[rname] if rname in design.components else None
        attached = set()
        if comp is not None and 'pin_inputs' in comp.options:
            for side_k in ('start_pin', 'end_pin'):
                pi = comp.options.pin_inputs.get(side_k, {})
                if pi and pi.get('component'):
                    attached.add(pi['component'])
        for cname, cgeom in cap_boxes:
            if cname != rname and cname not in attached and trace_poly.intersects(cgeom):
                collisions.append((rname, cname, 'cap_pocket'))
    if collisions:
        sample = ', '.join(f'{r}x{c}.{p}' for r, c, p in collisions[:12])
        raise AssertionError(f'Found {len(collisions)} layer-1 route-to-component collisions: {sample}')


def _net_key(name: str) -> str:
    """Electrical net of a drawn/planned route: strip bridge split (_sN) and SPAN-joined (_in/_out/_tap) suffixes."""
    base = re.sub(r'_s\d+$', '', name)
    return re.sub(r'_(in|out|tap)$', '', base)


def _pin_xy_dir(design, comp, pin):
    p = design.components[comp].pins[pin]
    return tuple(float(v) for v in p.middle), tuple(float(v) for v in p.normal)


def _route_control_lines(design, planned_routes, side_exits, pad_pools, h=0.05,
                         d_min=0.18, d_pref=0.25):
    """Assign perimeter pads monotonically, then route every control line with the clearance-aware
    octilinear A* router (code/general/octilinear_router.py). Appends the routed nets to planned_routes."""
    from octilinear_router import OctiGrid, RouteFailed, dir_index, polyline_ok

    r_hard = d_min + 0.01          # interior runs sit on cell centres; exact check is _assert_cpw_clearance
    r_soft = d_pref + 0.01
    grid = OctiGrid(-DIE_HALF, DIE_HALF, -DIE_HALF, DIE_HALF, h=h)

    # --- hard obstacles: every placed (non-planned) component, buffered by trench + gap clearance
    id_to_name = {c.id: n for n, c in design.components.items()}
    geoms = {}
    for tbl in ('poly', 'path', 'junction'):
        if tbl not in design.qgeometry.tables:
            continue
        for row in design.qgeometry.tables[tbl].itertuples():
            if int(row.layer) != 1:
                continue
            nm = id_to_name.get(row.component, '')
            g = row.geometry
            if tbl == 'path':
                g = g.buffer(max(float(row.width), 0.0) / 2.0 + 0.006)
            geoms.setdefault(nm, []).append(g)
    for nm, gl in geoms.items():
        grid.add_hard(unary_union(gl).buffer(0.065), owner=f'comp:{nm}')

    # --- die margin just inside the launch-pad tie pins
    tie_abs = min(max(abs(v) for v in _pin_xy_dir(design, 'LP_' + pid, 'tie')[0])
                  for pool in pad_pools.values() for pid in pool)
    grid.add_hard_rect_outside(-tie_abs + 0.02, tie_abs - 0.02, -tie_abs + 0.02, tie_abs - 0.02, owner='die_margin')

    # --- fixed planned routes: readout / Purcell (oriented buffers + hard meander hulls) and short control stubs
    for r in planned_routes:
        pa, _ = _pin_xy_dir(design, *r['start'])
        pb, _ = _pin_xy_dir(design, *r['end'])
        poly = [pa] + [tuple(map(float, p)) for p in r['anchors']] + [pb]
        grid.add_line(poly, owner=_net_key(r['name']), radius=r_hard, soft_radius=r_soft)
        if r['name'].startswith('ro_') and not r['name'].endswith('_tap') and len(poly) >= 6:
            # meander only: (u_hi, v0) ... short; never the lead's sideways turn near the qubit
            hull = MultiPoint(poly[-6:]).convex_hull.buffer(r_hard)
            grid.add_hard(hull, owner=f'hull:{r["name"]}')

    # --- monotonic perimeter pad assignment (XY + flux sorted together per side)
    nets = []
    for side, items in side_exits.items():
        items.sort(key=lambda it: it['sort_key'])
        pool = pad_pools[side]
        n_it = len(items)
        idxs = ([int(round(k * (len(pool) - 1) / max(1, n_it - 1))) for k in range(n_it)]
                if n_it > 1 else [len(pool) // 2])
        for it, pi in zip(items, idxs):
            pin_xy, pin_n = _pin_xy_dir(design, *it['start'])
            prefix = [pin_xy] + [tuple(map(float, p)) for p in it['pre']]
            if len(prefix) >= 2:
                v = (prefix[-1][0] - prefix[-2][0], prefix[-1][1] - prefix[-2][1])
            else:
                v = pin_n
            d_s = dir_index(v)
            if d_s is None:
                raise AssertionError(f"{it['name']}: start direction {v} is not octilinear")
            pad = 'LP_' + pool[pi]
            tie_xy, tie_n = _pin_xy_dir(design, pad, 'tie')
            d_e = (dir_index(tie_n) + 4) % 8
            nets.append(dict(name=it['name'], key=_net_key(it['name']), start=it['start'], pad=pad,
                             prefix=prefix, d_s=d_s, tie=tie_xy, tie_n=tie_n, d_e=d_e))

    # --- reserve every pending net's fixed prefix + pad approach so earlier nets cannot block them
    from octilinear_router import DX, DY
    for n in nets:
        p_last = n['prefix'][-1]
        nrm = math.hypot(DX[n['d_s']], DY[n['d_s']])
        # axis-aligned pin exits get a long protected escape corridor; diagonal XY stubs a short one so
        # they can still turn away from neighbouring parallel exits
        # (a diagonal XY stub reserves a 0.40 mm escape toward its pad side when that is <= 45 deg away)
        d_res = n['d_s']
        if n['d_s'] % 2 == 1 and min((n['d_s'] - n['d_e']) % 8, (n['d_e'] - n['d_s']) % 8) == 1:
            d_res = n['d_e']
        L_res = 0.40 if d_res % 2 == 0 else 0.15
        nrm = math.hypot(DX[d_res], DY[d_res])
        stub = (p_last[0] + L_res * DX[d_res] / nrm, p_last[1] + L_res * DY[d_res] / nrm)
        grid.add_line(n['prefix'] + [stub], owner=n['key'], radius=r_hard, soft_radius=r_soft)
        tie = n['tie']
        appr = (tie[0] + 0.30 * n['tie_n'][0], tie[1] + 0.30 * n['tie_n'][1])
        grid.add_line([tie, appr], owner=n['key'], radius=r_hard, soft_radius=r_soft)

    # --- route innermost nets first; rip up the nets that block a failed net's pin escapes, route the
    #     failed net first, then re-route the ripped ones (bounded number of rounds)
    nets.sort(key=lambda n: math.hypot(*n['prefix'][-1]))
    by_name = {n['name']: n for n in nets}
    routed, failed = {}, []

    def _try(n):
        poly = grid.route(n['key'], n['prefix'][-1], n['d_s'], n['tie'], n['d_e'],
                          ignore_owners=(('route', n['key']),), start_is_anchor=len(n['prefix']) >= 2)
        if not polyline_ok(poly, min_seg=0.02):
            print(f'WARNING {n["name"]}: degenerate polyline after re-intersection: '
                  f'{[(round(x, 3), round(y, 3)) for x, y in poly]}')
        grid.add_line(poly, owner=('route', n['key']), radius=r_hard, soft_radius=r_soft)
        routed[n['name']] = n['prefix'][1:] + poly[1:-1]

    def _blockers(n, radius=0.9):
        """Routed nets whose clearance buffers occupy the escape zones around n's start and end."""
        cells = set()
        for (px, py) in (n['prefix'][-1], n['tie']):
            ix0, iy0 = grid.cell(px, py)
            rr = int(math.ceil(radius / grid.h))
            for iy in range(max(0, iy0 - rr), min(grid.ny, iy0 + rr + 1)):
                for ix in range(max(0, ix0 - rr), min(grid.nx, ix0 + rr + 1)):
                    if math.hypot(ix - ix0, iy - iy0) * grid.h <= radius:
                        cells.add(iy * grid.nx + ix)
        cells = np.fromiter(cells, np.int64)
        out = []
        for owner, flat, _bits in grid.items:
            if isinstance(owner, tuple) and owner[0] == 'route' and owner[1] != n['key'] and owner not in out:
                if np.isin(flat, cells, assume_unique=False).any():
                    out.append(owner)
        return [o[1] for o in out]

    queue = list(nets)
    rounds = {n['name']: 0 for n in nets}
    while queue:
        n = queue.pop(0)
        try:
            _try(n)
            continue
        except RouteFailed as exc:
            err = str(exc)
        rounds[n['name']] += 1
        if rounds[n['name']] > 4:
            failed.append(err)
            continue
        rip = [b for b in _blockers(n) if b != n['key']]
        if not rip:
            failed.append(err)
            continue
        ripped = []
        for key in rip:
            grid.remove_owner(('route', key))
            for nm in [m['name'] for m in nets if m['key'] == key and m['name'] in routed]:
                routed.pop(nm)
                ripped.append(by_name[nm])
        print(f'rip-up for {n["name"]}: {[m["name"] for m in ripped]}')
        queue = [n] + ripped + [q for q in queue if q not in ripped]
    _route_control_lines.debug = dict(grid=grid, nets=nets, routed=routed, failed=failed)
    if failed:
        raise AssertionError(f'A* router failed for {len(failed)} control nets: {failed}')

    for n in nets:
        planned_routes.append(dict(
            name=n['name'], priority=1, start=n['start'], end=(n['pad'], 'tie'),
            anchors=routed[n['name']], cpw=CPW_CTL, fillet=FILLET_S17,
        ))


def _control_and_readout_nets(design):
    """Group drawn layer-1 CPW path rows into electrical nets for the clearance checker."""
    id_to_name = {c.id: n for n, c in design.components.items()}
    nets = {}
    for row in design.qgeometry.tables['path'].itertuples():
        if row.subtract or int(row.layer) != 1:
            continue
        nm = id_to_name.get(row.component, '')
        if not nm.startswith(('xy_', 'fl_', 'ro_', 'PF_body_', 'pf_feed_')):
            continue
        key = _net_key(nm)
        pr = 1 if nm.startswith(('xy_', 'fl_')) else 0
        nets.setdefault(key, dict(polys=[], priority=pr))['polys'].append(list(row.geometry.coords))
    return nets


def _assert_cpw_clearance(design, d_min_mm=0.18):
    """Hard Stage-1 DRC (skill/airbridge.md §1): true centerline clearance >= d_min between distinct nets,
    zero oblique crossings; reports consecutive-bridge spacing < 500 um along control lines."""
    from cpw_clearance import bridge_spacing_violations, clearance_violations, crossing_report

    nets = _control_and_readout_nets(design)
    bridges = [(float(design.parse_value(c.options.pos_x)), float(design.parse_value(c.options.pos_y)))
               for n, c in design.components.items() if n.startswith(('ABX_', 'SPAN_'))]
    bad = clearance_violations(nets, d_min=d_min_mm, min_cross_deg=85.0, bridges=bridges)
    allc, obl = crossing_report(nets, min_cross_deg=85.0)
    sp = bridge_spacing_violations(nets, allc, min_gap=0.50)
    # min clearance excluding crossing discs (report value)
    min_um = min([b[0] for b in bad], default=None)
    if obl:
        raise AssertionError(f'Found {len(obl)} oblique (<85 deg) CPW crossings: {obl[:10]}')
    if bad:
        raise AssertionError(f'Found {len(bad)} CPW clearance violations (<{d_min_mm*1000:.0f} um): {bad[:12]}')
    if sp:
        print(f'NOTE: {len(sp)} consecutive crossover pairs < 500 um along a control line: {sp[:12]}')
    return dict(min_um=f'>= {d_min_mm*1000:.0f}' if min_um is None else min_um, bridge_spacing=len(sp),
                n_crossings=len(allc))


def _assert_ground_bridge_coverage(design, max_gap_mm=0.55, end_allow_mm=0.12):
    """Hard Stage-1 DRC (skill/airbridge.md §4): every control-line run > max_gap carries a ground bridge."""
    centers = []
    for n, c in design.components.items():
        if n.startswith('ABG_'):
            centers.append((float(design.parse_value(c.options.pos_x)), float(design.parse_value(c.options.pos_y))))
    centers = np.asarray(centers, float) if centers else np.zeros((0, 2))
    id_to_name = {c.id: n for n, c in design.components.items()}
    bad, worst = [], 0.0
    for row in design.qgeometry.tables['path'].itertuples():
        if row.subtract or int(row.layer) != 1:
            continue
        nm = id_to_name.get(row.component, '')
        if not nm.startswith(('xy_', 'fl_')):
            continue
        ls = row.geometry
        L = ls.length
        pos = []
        if len(centers):
            near = centers[(centers[:, 0] > ls.bounds[0] - 0.02) & (centers[:, 0] < ls.bounds[2] + 0.02) &
                           (centers[:, 1] > ls.bounds[1] - 0.02) & (centers[:, 1] < ls.bounds[3] + 0.02)]
            for x, y in near:
                pt = Point(x, y)
                if ls.distance(pt) < 0.005:
                    pos.append(ls.project(pt))
        anchors = [0.0] + sorted(pos) + [L]
        for k in range(len(anchors) - 1):
            gap = anchors[k + 1] - anchors[k]
            if k == 0:
                gap -= end_allow_mm
            if k == len(anchors) - 2:
                gap -= end_allow_mm
            worst = max(worst, gap)
            if gap > max_gap_mm:
                p = ls.interpolate(0.5 * (anchors[k] + anchors[k + 1]))
                bad.append((nm, round(gap * 1000.0), round(p.x, 2), round(p.y, 2)))
    if bad:
        bad.sort(key=lambda b: -b[1])
        raise AssertionError(f'Found {len(bad)} control-line runs without ground bridge (> {max_gap_mm*1000:.0f} um): {bad[:12]}')
    return round(worst * 1000.0)


def _build_surface17(variant='sim', overrides=None, flux_as_ground=False):
    jjm = 'True' if variant == 'sim' else 'False'
    fag = 'True' if flux_as_ground else 'False'
    ov = overrides or {}

    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon_eps11p45'
    design.chips.main.size.update(dict(size_x=f'{DIE_MM:g}mm', size_y=f'{DIE_MM:g}mm', size_z='-500um', center_x='0mm', center_y='0mm'))
    design.variables.update(cpw_width=CPW_RO['trace_width'], cpw_gap=CPW_RO['trace_gap'])

    # 1. Place 17 circular wide-gap ScalablePadQubits (45-deg radial readout exit!)
    q_base = dict(QUBIT_OPT, make_readout_pad='True', ro_exit_mode='diag')
    q_base.update(ov.get('qubit', {}))
    for qname, qinfo in SURFACE17_QUBITS.items():
        ScalablePadQubit(design, qname, options=dict(
            **_pos(qinfo['i'] * PITCH_MM, qinfo['j'] * PITCH_MM),
            ro_corner=qinfo['ro_corner'], jj_corner=_jj_corner(qname),
            Lj=qinfo['Lj'], jj_metal=jjm, **q_base,
        ))

    # 2. Place 24 FloatingTunableCouplers on all nearest-neighbor data-ancilla edges
    tc_base = dict(COUPLER_OPT)
    tc_base.update(ov.get('coupler', {}))
    edges = _surface17_edges()
    tc_list = []
    for na, nb in edges:
        qa, qb = SURFACE17_QUBITS[na], SURFACE17_QUBITS[nb]
        xm = 0.5 * (qa['i'] + qb['i']) * PITCH_MM
        ym = 0.5 * (qa['j'] + qb['j']) * PITCH_MM
        is_horiz = (qa['j'] == qb['j'])
        orient = _coupler_orientation(na, nb, xm, ym, is_horiz)
        tc_name = f'TC_{na}_{nb}'
        tc_list.append((tc_name, na, nb, xm, ym, orient))
        FloatingTunableCoupler(design, tc_name, options=dict(
            **_pos(xm, ym), orientation=orient, jj_metal=jjm, flux_as_ground=fag, **tc_base,
        ))

    # 3. Place 56 perimeter launch pads on 14 x 14 mm die (14 per side, 0.75 mm pitch)
    lp_opt = dict(pad_width='240um', pad_height='260um', pad_gap='110um',
                  taper_height='200um', back_gap='50um', lead_length='30um', **CPW_CTL)
    for pname, pdict in _perimeter_pads_14mm().items():
        BondPad(design, 'LP_' + pname, options=dict(pdict, **lp_opt))

    # 4. Place 10 standard <= 45 um-span CPW airbridges STRICTLY over the 10/6 um CPW necks of the couplers!
    #    NEVER crosses any wide qubit or coupler pad pattern!
    P = PITCH_MM
    span_specs = OrderedDict([
        ('SPAN_ro_D5',   dict(x=+1.0 * P,          y=+0.5 * P + 0.380, orient='90')),
        ('SPAN_fl_D5X1', dict(x=+1.0 * P,          y=+0.5 * P - 0.380, orient='90')),
        ('SPAN_xy_D5',   dict(x=-0.5 * P - 0.380,  y=+1.0 * P,         orient='0')),
        ('SPAN_fl_D5Z2', dict(x=+0.5 * P + 0.380,  y=-1.0 * P,         orient='0')),
        ('SPAN_fl_X2D5', dict(x=+0.5 * P - 0.380,  y=-1.0 * P,         orient='0')),
        ('SPAN_fl_Z1D5', dict(x=-1.0 * P,          y=-0.5 * P - 0.380, orient='90')),
        ('SPAN_fl_Z1D8', dict(x=-2.0 * P,          y=+0.5 * P - 0.500, orient='90')),
        ('SPAN_fl_X1D6', dict(x=+0.5 * P - 0.500,  y=+2.0 * P,         orient='0')),
        ('SPAN_fl_D2Z2', dict(x=+2.0 * P,          y=-0.5 * P + 0.500, orient='90')),
        ('SPAN_fl_D4X2', dict(x=-0.5 * P + 0.500,  y=-2.0 * P,         orient='0')),
        ('SPAN_xy_X1',   dict(x=+0.5 * P - 0.290,  y=+2.0 * P,         orient='0')),
        ('SPAN_xy_X2',   dict(x=-0.5 * P + 0.290,  y=-2.0 * P,         orient='0')),
        ('SPAN_xy_Z1',   dict(x=-2.0 * P,          y=+0.5 * P - 0.290, orient='90')),
        ('SPAN_xy_Z2',   dict(x=+2.0 * P,          y=-0.5 * P + 0.290, orient='90')),
    ])
    for sname, sp in span_specs.items():
        AirbridgeSpan(design, sname, options=dict(
            **_pos(sp['x'], sp['y']), orientation=sp['orient'],
            clear_span='30um', total_length='58um',
            pad_width='10um', pad_length='6um', bridge_width='8um', pier_size='5um',
            under_angle='0', under_keepout='15um', end_gap='2.5um',
            **CPW_CTL,
        ))

    # 5. Place 4 DIAGONAL (45-deg) SHARED-BUS lambda/2 PURCELL FILTERS IMMEDIATELY OUTSIDE THE DIAMOND!
    u_purcell = U_PURCELL_MM
    purcell_quads = OrderedDict([
        ('NW', dict(sx=-1, sy=+1, p_in='LP_N00', p_out='LP_W13', v_span=(-3.50, +3.50), cin_ori='-90', cout_ori='180')),
        ('NE', dict(sx=+1, sy=+1, p_in='LP_E13', p_out='LP_N13', v_span=(-4.10, +3.30), cin_ori='180', cout_ori='90')),
        ('SW', dict(sx=-1, sy=-1, p_in='LP_W00', p_out='LP_S00', v_span=(-3.50, +3.50), cin_ori='0',   cout_ori='-90')),
        ('SE', dict(sx=+1, sy=-1, p_in='LP_S13', p_out='LP_E00', v_span=(-3.50, +3.50), cin_ori='90',  cout_ori='0')),
    ])

    # Collect all routes so every CPW-over-CPW crossing (45 deg or 90 deg) gets an explicit AirbridgeSpan
    # on the 50-ohm overpass ONLY (priority 0 = Purcell/readout ALWAYS continuous in layer 1!).
    planned_routes = []

    for qk, qs in purcell_quads.items():
        sx, sy = qs['sx'], qs['sy']
        v_lo, v_hi = qs['v_span']
        pin_in = design.components[qs['p_in']].pins['tie']
        pin_out = design.components[qs['p_out']].pins['tie']
        c_in_xy = tuple(np.asarray(pin_in.middle, float) + 0.45 * np.asarray(pin_in.normal, float))
        c_out_xy = tuple(np.asarray(pin_out.middle, float) + 0.45 * np.asarray(pin_out.normal, float))
        EnclosingCap(design, f'PF_cin_{qk}', options=dict(
            **_pos(*c_in_xy), orientation=qs['cin_ori'], **ENC_OPT, **CPW_RO,
        ))
        EnclosingCap(design, f'PF_cout_{qk}', options=dict(
            **_pos(*c_out_xy), orientation=qs['cout_ori'], **ENC_OPT, **CPW_RO,
        ))
        RouteStraight(design, f'pf_feed_in_{qk}', Dict(
            pin_inputs=Dict(start_pin=Dict(component=qs['p_in'], pin='tie'),
                            end_pin=Dict(component=f'PF_cin_{qk}', pin='a')), **CPW_RO,
        ))
        RouteStraight(design, f'pf_feed_out_{qk}', Dict(
            pin_inputs=Dict(start_pin=Dict(component=f'PF_cout_{qk}', pin='b'),
                            end_pin=Dict(component=qs['p_out'], pin='tie')), **CPW_RO,
        ))
        S = u_purcell * math.sqrt(2.0)
        x_c = PF_XC_MM
        y_c = S - x_c
        x_outer = x_c + 0.28
        cin_b = design.components[f'PF_cin_{qk}'].pins['b']
        cout_a = design.components[f'PF_cout_{qk}'].pins['a']
        b_xy = tuple(np.asarray(cin_b.middle, float))
        a_xy = tuple(np.asarray(cout_a.middle, float))
        if abs(cin_b.normal[0]) > 0.5:
            pf_pts = [
                (sx * x_outer, b_xy[1]),
                (sx * x_outer, sy * y_c),
                (sx * x_c, sy * y_c),
                (sx * y_c, sy * x_c),
                (sx * y_c, sy * x_outer),
                (a_xy[0], sy * x_outer),
            ]
        else:
            pf_pts = [
                (b_xy[0], sy * x_outer),
                (sx * y_c, sy * x_outer),
                (sx * y_c, sy * x_c),
                (sx * x_c, sy * y_c),
                (sx * x_outer, sy * y_c),
                (sx * x_outer, a_xy[1]),
            ]
        planned_routes.append(dict(
            name=f'PF_body_{qk}', priority=0,
            start=(f'PF_cin_{qk}', 'b'), end=(f'PF_cout_{qk}', 'a'),
            anchors=pf_pts, cpw=CPW_RO, fillet='70um',
        ))

    # 6. Build ALL 17 EXACT lambda/4 (4.17-4.36 mm TOTAL!) DRC-B7 READOUT RESONATORS at 45 deg!
    #    ZERO airbridges on any readout resonator meander!
    quad_lists = {k: [] for k in purcell_quads}
    for qname, qinfo in SURFACE17_QUBITS.items():
        quad_lists[qinfo['ro_quad']].append(qname)

    res_report = {}
    for quad, qnames in quad_lists.items():
        sx = purcell_quads[quad]['sx']
        sy = purcell_quads[quad]['sy']
        # Sort by actual readout start-pin tangent coordinate v so zero readout resonators cross one another!
        def _ro_v(qn):
            p = (design.components['SPAN_ro_D5'].pins['over_bot'].middle if qn == 'D5'
                 else design.components[qn].pins['readout'].middle)
            return _xy_to_diag_uv(float(p[0]), float(p[1]), sx, sy)[1]
        qnames.sort(key=_ro_v)
        n_q = len(qnames)
        v_lo, v_hi = purcell_quads[quad]['v_span']
        v_slots = np.linspace(v_lo + 0.48, v_hi - 1.05, n_q)
        meander_w = 2.0 * 0.35 + 0.25                       # 2 meander pitches + coupling leg (v extent)
        # sideways shifts happen only outside the lattice (outermost qubit u + pocket + margin)
        u_clear = max(_xy_to_diag_uv(SURFACE17_QUBITS[q]['i'] * PITCH_MM, SURFACE17_QUBITS[q]['j'] * PITCH_MM,
                                     sx, sy)[0] for q in SURFACE17_QUBITS) + 0.19 + 0.20
        pin_v = [_ro_v(q) for q in qnames]
        early = [any(k != qi and v_slots[k] - 0.15 <= pin_v[qi] <= v_slots[k] + meander_w + 0.15
                     for k in range(n_q)) for qi in range(n_q)]
        tang_deg = math.degrees(math.atan2(sx, -sy))

        for qi, qname in enumerate(qnames):
            qinfo = SURFACE17_QUBITS[qname]
            L_target = V_PH_MM_GHZ / (4.0 * qinfo['f_ro_GHz'])
            if qname == 'D5':
                d5_pin = design.components['D5'].pins['readout']
                d5_lead = tuple(np.asarray(d5_pin.middle, float) + 0.22 * np.asarray(d5_pin.normal, float))
                sp_in = design.components['SPAN_ro_D5'].pins['over_top']
                sp_in_lead = tuple(np.asarray(sp_in.middle, float) + 0.22 * np.asarray(sp_in.normal, float))
                # Octilinear tap: continue the 45-deg readout ray, then run along the span axis into SPAN_ro_D5
                # (a straight non-octilinear tap would cross fl_TC_D5_X1_in obliquely).
                sp_n = np.asarray(sp_in.normal, float)
                if abs(sp_n[0]) > 0.5:   # span pin faces +-x -> horizontal final run at y = sp_in_lead.y
                    t_corner = (sp_in_lead[1] - d5_lead[1]) / float(d5_pin.normal[1])
                else:                    # span pin faces +-y -> vertical final run at x = sp_in_lead.x
                    t_corner = (sp_in_lead[0] - d5_lead[0]) / float(d5_pin.normal[0])
                tap_corner = tuple(np.asarray(d5_lead) + t_corner * np.asarray(d5_pin.normal, float))
                planned_routes.append(dict(
                    name='ro_D5_tap', priority=0,
                    start=('D5', 'readout'), end=('SPAN_ro_D5', 'over_top'),
                    anchors=[d5_lead, tap_corner, sp_in_lead], cpw=CPW_RO, fillet=FILLET_S17,
                ))
                sp_out_pin = design.components['SPAN_ro_D5'].pins['over_bot']
                sp_out_lead = tuple(np.asarray(sp_out_pin.middle, float) + 0.22 * np.asarray(sp_out_pin.normal, float))
                m_anc, short_xy, L_tot = _build_diag_lambda4_resonator(
                    sp_out_lead, sx, sy, float(v_slots[qi]), u_purcell, L_target - 0.90,
                    fillet_mm=0.055, pitch_mm=0.35, coup_leg_mm=0.25, early_turn=early[qi], u_clear=u_clear,
                )
                ShortToGround(design, 'ro_short_D5', options=dict(
                    **_pos(*short_xy), orientation=f'{tang_deg:.1f}', width=CPW_RO['trace_width'],
                ))
                planned_routes.append(dict(
                    name='ro_D5', priority=0,
                    start=('SPAN_ro_D5', 'over_bot'), end=('ro_short_D5', 'short'),
                    anchors=[sp_out_lead] + m_anc, cpw=CPW_RO, fillet=FILLET_S17,
                ))
                res_report[qname] = dict(f_GHz=qinfo['f_ro_GHz'], L_target_mm=round(L_target, 3), L_actual_mm=round(L_tot + 0.90, 3))
                continue

            p_pin = tuple(np.asarray(design.components[qname].pins['readout'].middle, float))
            m_anc, short_xy, L_tot = _build_diag_lambda4_resonator(
                p_pin, sx, sy, float(v_slots[qi]), u_purcell, L_target,
                fillet_mm=0.055, pitch_mm=0.35, coup_leg_mm=0.25, early_turn=early[qi], u_clear=u_clear,
            )
            ShortToGround(design, f'ro_short_{qname}', options=dict(
                **_pos(*short_xy), orientation=f'{tang_deg:.1f}', width=CPW_RO['trace_width'],
            ))
            planned_routes.append(dict(
                name=f'ro_{qname}', priority=0,
                start=(qname, 'readout'), end=(f'ro_short_{qname}', 'short'),
                anchors=m_anc, cpw=CPW_RO, fillet=FILLET_S17,
            ))
            res_report[qname] = dict(f_GHz=qinfo['f_ro_GHz'], L_target_mm=round(L_target, 3), L_actual_mm=round(L_tot, 3))

    # 7 & 8. MONOTONIC PERIMETER PAD SORTING (17 XY + 24 TC flux lines sorted TOGETHER per side!)
    #        Sorting XY + Flux monotonically by exit coordinate (x for N/S, y for E/W) eliminates
    #        ALL same-side criss-crossing and minimizes series airbridges per signal line!
    pad_pools = {
        'N': [f'N{k:02d}' for k in range(1, 13)],
        'S': [f'S{k:02d}' for k in range(1, 13)],
        'W': [f'W{k:02d}' for k in range(1, 13)],
        'E': [f'E{k:02d}' for k in range(1, 13)],
    }
    side_exits = {'N': [], 'S': [], 'W': [], 'E': []}

    for qname, qinfo in SURFACE17_QUBITS.items():
        side = qinfo['xy_side']
        qx, qy = qinfo['i'] * PITCH_MM, qinfo['j'] * PITCH_MM
        if qname == 'D5':
            # Aim D5 XY charge line diagonally at 135 deg directly toward D5 center (far from both adjacent couplers!)
            OpenToGround(design, 'xy_end_D5', options=dict(
                **_pos(-0.173, +0.173), orientation='-45', width=CPW_CTL['trace_width'], gap=CPW_CTL['trace_gap'],
            ))
            sp_in_xy = tuple(np.asarray(design.components['SPAN_xy_D5'].pins['over_bot'].middle, float))
            planned_routes.append(dict(
                name='xy_D5_in', priority=1,
                start=('xy_end_D5', 'open'), end=('SPAN_xy_D5', 'over_bot'),
                anchors=[(-0.34, +0.34), (sp_in_xy[0], +0.34)], cpw=CPW_CTL, fillet=FILLET_S17,
            ))
            sp_out = tuple(np.asarray(design.components['SPAN_xy_D5'].pins['over_top'].middle, float))
            p_exit = (sp_out[0], 1.36 * P)
            side_exits['W'].append(dict(
                sort_key=p_exit[1], name='xy_D5_out',
                start=('SPAN_xy_D5', 'over_top'), p_exit=p_exit,
                pre=[], lead=3.20,
            ))
            continue

        # Aim EVERY XY charge line radially at 45 deg toward its circular qubit center from an EMPTY diagonal corner
        # (>= 90 deg away from its readout resonator AND 45 deg away from all 4 cardinal coupler pads!).
        r_xy = 0.245
        if qname in ('X1', 'X2', 'Z1', 'Z2'):
            # X1 ro=NW -> xy=NE(+45); X2 ro=SE -> xy=SW(-135); Z1 ro=SW -> xy=NW(+135); Z2 ro=NE -> xy=SE(-45)
            # Route cleanly over the outer 585 um coupler neck via SPAN_xy_{qname} (90 deg away from readout!).
            diag_map = {'X1': (+1.0, +1.0, 'over_bot', 'over_top'),
                        'X2': (-1.0, -1.0, 'over_top', 'over_bot'),
                        'Z1': (-1.0, +1.0, 'over_bot', 'over_top'),
                        'Z2': (+1.0, -1.0, 'over_top', 'over_bot')}
            dx_s, dy_s, pin_in, pin_out = diag_map[qname]
            assert _xy_corner(qname) == ('N' if dy_s > 0 else 'S') + ('E' if dx_s > 0 else 'W'), qname
            sx = qx + dx_s * r_xy * INV_SQRT2
            sy = qy + dy_s * r_xy * INV_SQRT2
            ori_in = math.degrees(math.atan2(-dy_s, -dx_s))
            OpenToGround(design, f'xy_end_{qname}', options=dict(
                **_pos(sx, sy), orientation=f'{ori_in:.1f}', width=CPW_CTL['trace_width'], gap=CPW_CTL['trace_gap'],
            ))
            sname = f'SPAN_xy_{qname}'
            sp_in = tuple(np.asarray(design.components[sname].pins[pin_in].middle, float))
            sp_out = tuple(np.asarray(design.components[sname].pins[pin_out].middle, float))
            if side in ('N', 'S'):
                d_diag = abs(sp_in[0] - qx)
                p_turn = (sp_in[0], qy + dy_s * d_diag)
            else:
                d_diag = abs(sp_in[1] - qy)
                p_turn = (qx + dx_s * d_diag, sp_in[1])
            planned_routes.append(dict(
                name=f'xy_{qname}_in', priority=1,
                start=(f'xy_end_{qname}', 'open'), end=(sname, pin_in),
                anchors=[p_turn], cpw=CPW_CTL, fillet=FILLET_S17,
            ))
            side_exits[side].append(dict(
                sort_key=sp_out[0] if side in ('N', 'S') else sp_out[1],
                name=f'xy_{qname}', start=(sname, pin_out),
                p_exit=sp_out, pre=[], lead=0.80,
            ))
            continue
        else:
            if side == 'N':
                dx_s, dy_s = (-1.0, +1.0) if qname == 'Z3' else (+1.0, +1.0)
            elif side == 'S':
                dx_s, dy_s = (+1.0, -1.0) if qname == 'Z4' else (-1.0, -1.0)
            elif side == 'W':
                dx_s, dy_s = (-1.0, -1.0) if qname == 'X3' else (-1.0, +1.0)
            else:
                dx_s, dy_s = (+1.0, +1.0) if qname == 'X4' else (+1.0, -1.0)
            assert _xy_corner(qname) == ('N' if dy_s > 0 else 'S') + ('E' if dx_s > 0 else 'W'), qname
            sx = qx + dx_s * r_xy * INV_SQRT2
            sy = qy + dy_s * r_xy * INV_SQRT2
            ori_in = math.degrees(math.atan2(-dy_s, -dx_s))
            OpenToGround(design, f'xy_end_{qname}', options=dict(
                **_pos(sx, sy), orientation=f'{ori_in:.1f}', width=CPW_CTL['trace_width'], gap=CPW_CTL['trace_gap'],
            ))
            d_step = 0.285 if qname in ('D1', 'D3', 'D4', 'D7', 'D9') else 0.315
            p_diag = (qx + dx_s * d_step, qy + dy_s * d_step)
        lead = 2.45 * P - abs(p_diag[1]) if (side in ('N', 'S') and abs(p_diag[1]) < 2.1 * P) else (
               2.45 * P - abs(p_diag[0]) if (side in ('E', 'W') and abs(p_diag[0]) < 2.1 * P) else 0.35)
        # Shift xy_D2, xy_Z3, xy_Z4 fanout diagonals to avoid triple junctions and keep 100% 90-deg bridges!
        if qname == 'D2':
            lead += 0.55
        elif qname == 'Z3':
            lead += 1.15
        elif qname == 'Z4':
            lead += 0.55
        side_exits[side].append(dict(
            sort_key=p_diag[0] if side in ('N', 'S') else p_diag[1],
            name=f'xy_{qname}', start=(f'xy_end_{qname}', 'open'),
            p_exit=p_diag, pre=[p_diag], lead=max(0.35, lead),
        ))

    if not flux_as_ground:
        span_tc_map = {
            frozenset(('D5', 'X1')): ('SPAN_fl_D5X1', 'over_top', 'over_bot', 'E'),
            frozenset(('D5', 'Z2')): ('SPAN_fl_D5Z2', 'over_top', 'over_bot', 'S'),
            frozenset(('X2', 'D5')): ('SPAN_fl_X2D5', 'over_top', 'over_bot', 'S'),
            frozenset(('Z1', 'D5')): ('SPAN_fl_Z1D5', 'over_bot', 'over_top', 'W'),
            frozenset(('Z1', 'D8')): ('SPAN_fl_Z1D8', 'over_bot', 'over_top', 'W'),
            frozenset(('X1', 'D6')): ('SPAN_fl_X1D6', 'over_bot', 'over_top', 'N'),
            frozenset(('D2', 'Z2')): ('SPAN_fl_D2Z2', 'over_top', 'over_bot', 'E'),
            frozenset(('D4', 'X2')): ('SPAN_fl_D4X2', 'over_top', 'over_bot', 'S'),
        }
        orient_to_side = {'180': 'N', '0': 'S', '270': 'W', '90': 'E'}
        for tc_name, na, nb, xm, ym, orient in tc_list:
            pair = frozenset((na, nb))
            pin_xy = tuple(np.asarray(design.components[tc_name].pins['flux_bias_line'].middle, float))
            if pair in span_tc_map:
                sname, pin_in, pin_out, side = span_tc_map[pair]
                sp_in = tuple(np.asarray(design.components[sname].pins[pin_in].middle, float))
                sp_out = tuple(np.asarray(design.components[sname].pins[pin_out].middle, float))
                pin_n = np.asarray(design.components[tc_name].pins['flux_bias_line'].normal, float)
                sp_n = np.asarray(design.components[sname].pins[pin_in].normal, float)
                if abs(pin_xy[0] - sp_in[0]) < 0.01 or abs(pin_xy[1] - sp_in[1]) < 0.01:
                    in_pts = []
                elif abs(float(np.dot(pin_n, sp_n))) < 0.5:
                    in_pts = [(sp_in[0], pin_xy[1])] if abs(pin_n[0]) > 0.5 else [(pin_xy[0], sp_in[1])]
                elif orient in ('0', '180'):
                    y_mid = 0.5 * (pin_xy[1] + sp_in[1])
                    in_pts = [(pin_xy[0], y_mid), (sp_in[0], y_mid)]
                else:
                    if pair == frozenset(('D5', 'X1')):
                        dy = abs(pin_xy[1] - sp_in[1])
                        in_pts = [(0.75, pin_xy[1]), (0.75 + dy, sp_in[1])]
                    else:
                        x_mid = 0.5 * (pin_xy[0] + sp_in[0])
                        in_pts = [(x_mid, pin_xy[1]), (x_mid, sp_in[1])]
                planned_routes.append(dict(
                    name=f'fl_{tc_name}_in', priority=1,
                    start=(tc_name, 'flux_bias_line'), end=(sname, pin_in),
                    anchors=in_pts, cpw=CPW_CTL, fillet=FILLET_S17,
                ))
                # Give fl_TC_Z1_D5_out / fl_TC_D5_Z2_out a clean 160 um 135-deg shift to clear the SW/SE meander turn!
                pre_pts = []
                p_ex = sp_out
                if pair == frozenset(('Z1', 'D5')):
                    # Keep fl_TC_Z1_D5_out straight at y=-0.775 past x=-1.74 (clear of xy_D4 at y=-0.885!)
                    pre_pts, p_ex, lead = [], sp_out, 1.35
                elif pair == frozenset(('D5', 'Z2')):
                    # Keep fl_TC_D5_Z2_out straight at x=0.775 past ro_D1 (-3.175) before fanning out!
                    pre_pts, p_ex, lead = [], sp_out, 2.35
                elif pair in (frozenset(('D2', 'Z2')), frozenset(('D4', 'X2'))):
                    lead = 3.60
                elif pair == frozenset(('X1', 'D6')):
                    lead = 1.10
                else:
                    lead = 2.45 * P - abs(sp_out[1]) if (side in ('N', 'S') and abs(sp_out[1]) < 2.1 * P) else (
                           2.45 * P - abs(sp_out[0]) if (side in ('E', 'W') and abs(sp_out[0]) < 2.1 * P) else 0.35)
                side_exits[side].append(dict(
                    sort_key=p_ex[0] if side in ('N', 'S') else p_ex[1],
                    name=f'fl_{tc_name}_out', start=(sname, pin_out),
                    p_exit=p_ex, pre=pre_pts, lead=max(0.35, lead),
                ))
                continue

            side = orient_to_side[orient]
            # r7: with the rim-SQUID arc-T pin, TC_D1_X2's flux pin sits right next to the S-side bundle; sending it
            # to an S pad (parallel to the bundle) avoids crossing that bundle where its diagonals are.
            side = {frozenset(('D1', 'X2')): 'S'}.get(pair, side)
            base_k = (2.65 if pair == frozenset(('D6', 'Z3')) else
                      2.35 if pair == frozenset(('D7', 'Z1')) else
                      2.90 if pair == frozenset(('X1', 'D9')) else 3.08)
            lead = base_k * P - abs(pin_xy[1]) if (side in ('N', 'S') and abs(pin_xy[1]) < 2.9 * P) else (
                   base_k * P - abs(pin_xy[0]) if (side in ('E', 'W') and abs(pin_xy[0]) < 2.9 * P) else 0.35)
            pre_fl = []
            p_ex_fl = pin_xy
            if pair == frozenset(('D1', 'X2')):
                # reserve the narrow corridor beside the D1-X2 neck before the S-side bundle is routed
                # (otherwise fl_TC_D1_X2 / fl_TC_X2_D5_out / fl_TC_X2_D2 block each other in a rip-up cycle)
                pre_fl = [(pin_xy[0] + 0.12, pin_xy[1]), (pin_xy[0] + 0.12, pin_xy[1] - 0.60)]
                p_ex_fl = pre_fl[-1]
            # No hand-tuned detours: the clearance-aware A* router (skill/airbridge.md §4a) finds the path.
            side_exits[side].append(dict(
                sort_key=p_ex_fl[0] if side in ('N', 'S') else p_ex_fl[1],
                name=f'fl_{tc_name}', start=(tc_name, 'flux_bias_line'),
                p_exit=p_ex_fl, pre=pre_fl, lead=max(0.35, lead),
            ))

    # 9. Clearance-aware octilinear A* routing of ALL control lines (skill/airbridge.md §4a):
    #    >= 180 um centerline clearance by construction, crossings only straight at 90 deg.
    _route_control_lines(design, planned_routes, side_exits, pad_pools)

    def _build_polylines():
        out = []
        for r in planned_routes:
            p0 = tuple(np.asarray(design.components[r['start'][0]].pins[r['start'][1]].middle, float))
            p1 = tuple(np.asarray(design.components[r['end'][0]].pins[r['end'][1]].middle, float))
            out.append([p0] + [tuple(pt) for pt in r['anchors']] + [p1])
        return out

    # 9b. General-angle AirbridgeSpan insertion at EVERY CPW-over-CPW crossing!
    #     Purcell filters & readout resonators (priority 0) stay 100% continuous in layer 1 (NEVER broken!).
    polylines = _build_polylines()

    # Collect EVERY crossing per overpass route/segment first, then cluster closely-spaced crossings on the
    # same segment into one clean wider-span AirbridgeSpan so ZERO crossings are ever silently skipped!
    raw_by_route_seg = {}
    for i in range(len(planned_routes)):
        for j in range(i + 1, len(planned_routes)):
            xs = _segment_intersections_general(polylines[i], polylines[j])
            for xc, yc, seg_i, seg_j, ang_i, ang_j, cross_ang in xs:
                # Pick overpass: strictly higher priority (control line over readout/Purcell);
                # between two control lines, prefer whichever crossing is farther from its own segment endpoints!
                pi, pj = planned_routes[i]['priority'], planned_routes[j]['priority']
                if pi == 0 and pj == 0:
                    continue  # Never break a Purcell filter or readout resonator!
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
    oblique_bridges = []
    for (o_idx, o_seg), hits in raw_by_route_seg.items():
        hits.sort(key=lambda h: h[0])
        a0, a1 = polylines[o_idx][o_seg], polylines[o_idx][o_seg + 1]
        seg_len = math.hypot(a1[0] - a0[0], a1[1] - a0[1])
        ux = (a1[0] - a0[0]) / max(seg_len, 1e-9)
        uy = (a1[1] - a0[1]) / max(seg_len, 1e-9)
        # Non-chaining cluster guard: merge ONLY a tight pair < 105 um apart so bridge pins never overlap!
        clusters = []
        for h in hits:
            if clusters and (h[0] - clusters[-1][0][0]) < 0.105:
                clusters[-1].append(h)
            else:
                clusters.append([h])
        for cl in clusters:
            t_lo, t_hi = cl[0][0], cl[-1][0]
            t_mid = 0.5 * (t_lo + t_hi)
            min_cross = min(h[4] for h in cl)
            # Flush 10 um pad_width + >=4 um intact ground strip on both sides of underpass CPW:
            # 30 um clear_span at 90 deg (suspended pier-to-pier arch = 39 um, total L20 = 49 um!)
            if min_cross < 85.0 or len(cl) > 1:
                oblique_bridges.append((planned_routes[o_idx]['name'], round(cl[0][1], 3),
                                        round(cl[0][2], 3), round(min_cross, 1), len(cl)))
            base_span = 30.0 if min_cross > 70.0 else 50.0
            span_um = base_span + (t_hi - t_lo) * 1000.0
            tot_um = span_um + 28.0
            pad_half = 0.5 * tot_um / 1000.0
            t_mid = max(pad_half + 0.025, min(seg_len - pad_half - 0.025, t_mid))
            xc = a0[0] + ux * t_mid
            yc = a0[1] + uy * t_mid
            o_ang = cl[0][3]
            u_rel = 0.0 if min_cross > 70.0 else (90.0 - min_cross)
            ab_count += 1
            ab_name = f'ABX_{ab_count:03d}'
            AirbridgeSpan(design, ab_name, options=dict(
                **_pos(xc, yc), orientation=f'{o_ang - 90.0:.2f}',
                clear_span=f'{span_um:.1f}um', total_length=f'{tot_um:.1f}um',
                pad_width='10um', pad_length='6um', bridge_width='8um', pier_size='5um',
                under_angle=f'{u_rel:.1f}', under_keepout='15um', end_gap='2.5um',
                **CPW_CTL,
            ))
            bridges_by_route[o_idx].append((o_seg, xc, yc, ab_name, 'over_bot', 'over_top'))

    for idx, r in enumerate(planned_routes):
        _draw_overpass_with_bridges(
            design, r['name'], r['start'], r['end'], r['anchors'],
            bridges_by_route[idx], r['cpw'], fillet=r['fillet'],
        )

    _assert_zero_pad_collisions(design)
    _assert_no_acute_bends(design, min_interior_deg=88.0)
    _assert_zero_unbridged_cpw_crossings(design)
    _assert_max_airbridge_length(design, max_l20_um=46.0)
    if oblique_bridges:
        raise AssertionError(f'Found {len(oblique_bridges)} oblique (<85 deg) or merged crossover bridges: '
                             f'{oblique_bridges[:10]}')
    clr = _assert_cpw_clearance(design, d_min_mm=0.18)
    n_gnd_ab = _place_ground_airbridges(design, spacing_mm=0.42)
    gnd_cov = _assert_ground_bridge_coverage(design, max_gap_mm=0.55)
    design._design13_summary = dict(
        mode='surface17', revision=REVISION, die_mm=(DIE_MM, DIE_MM),
        n_qubits=len(SURFACE17_QUBITS), n_couplers=len(edges),
        n_purcell_filters=4, n_readout_resonators=len(res_report),
        n_cpw_airbridges=ab_count + len(span_specs),
        n_ground_airbridges=n_gnd_ab,
        resonator_lengths_mm=res_report, pad_collisions=0,
        oblique_bridges=0, min_cpw_clearance_um=clr['min_um'],
        bridge_spacing_lt_500um=clr['bridge_spacing'], ground_bridge_max_gap_um=gnd_cov,
    )
    return design


def _build_unit_cell(variant='sim', overrides=None, channels=True, flux_as_ground=False, check_keepout_um=100.0):
    jjm = 'True' if variant == 'sim' else 'False'
    fag = 'True' if flux_as_ground else 'False'
    ov = overrides or {}

    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon_eps11p45'
    design.chips.main.size.update(dict(size_x='10mm', size_y='10mm', size_z='-500um', center_x='0mm', center_y='0mm'))
    design.variables.update(cpw_width=CPW_RO['trace_width'], cpw_gap=CPW_RO['trace_gap'])

    qx = 0.5 * PITCH_MM
    yq = 0.0
    q_opt = dict(QUBIT_OPT, make_readout_pad=str(channels), ro_exit_mode='ortho')
    q_opt.update(ov.get('qubit', {}))
    for name, sx, lj, rc in (('Q1', -1, '10.8nH', 'NW'), ('Q2', +1, '11.5nH', 'NE')):
        ScalablePadQubit(design, name, options=dict(
            **_pos(sx * qx, yq), ro_corner=rc, Lj=lj, jj_metal=jjm, **q_opt,
        ))

    tc_opt = dict(COUPLER_OPT, jj_metal=jjm, flux_as_ground=fag)
    tc_opt.update(ov.get('coupler', {}))
    FloatingTunableCoupler(design, 'TC', options=dict(**_pos(0.0, yq), **tc_opt))

    for pn in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'):
        if not channels and pn in ('P7', 'P8'):
            continue
        cpw = CPW_RO if pn in ('P7', 'P8') else CPW_CTL
        add_bond_pad(design, 'LP' + pn, wafer_pcb_ports[pn], **cpw)

    if channels:
        cap_pin = (30.0 / 2.0 + 40.0 + 16.0) / 1000.0
        fillet_loss = design.parse_value(FILLET_UC) * (2.0 - np.pi / 2.0)
        claw_q1 = np.asarray(design.components['Q1'].pins['readout'].middle, float)
        xf_abs = 1.80
        xr_abs = xf_abs - HAMMER_LINE_SPACING_UM / 1000.0
        y_jog = claw_q1[1] + 0.28
        horiz = abs(xr_abs - abs(claw_q1[0]))
        rt_mm = 0.32
        dy_top = 0.25
        ext_mm = dy_top + horiz + 0.28 - 2.0 * fillet_loss
        y_pin_target = y_jog + dy_top
        y_tap = y_pin_target + rt_mm
        feed_y = y_tap + 1.10 + 0.150 + 2.0 * cap_pin

        for name, ch in PURCELL_CHANNELS.items():
            s = '_' + name
            side = ch['side']
            claw = np.asarray(design.components[name].pins['readout'].middle, float)
            xf = side * xf_abs
            xr = side * xr_abs
            TapJunction(design, 'n0' + s, options=dict(
                **_pos(xf, feed_y), orientation='0',
                prime_length='300um', stub_length='100um', stub_up='False', **CPW_RO,
            ))
            y_cp = feed_y - 0.150 - cap_pin
            EnclosingCap(design, 'cp' + s, options=dict(
                **_pos(xf, y_cp), orientation='-90', **dict(ENC_OPT, plate_width=f"{ch['cp_plate_um']:g}um"),
                **CPW_RO,
            ))
            RouteStraight(design, 'stub_p' + s, Dict(
                pin_inputs=Dict(start_pin=Dict(component='n0' + s, pin='second_end'),
                                end_pin=Dict(component='cp' + s, pin='a')), **CPW_RO,
            ))
            cp_b = np.asarray(design.components['cp' + s].pins['b'].middle, float)
            y_tap_ch = cp_b[1] - 1.10
            x0 = 0.5 * (xf + xr)
            PurcellReadout(design, 'prr' + s, options=dict(
                **_pos(x0, y_tap_ch), orientation='0', mirror=str(side > 0),
                filter_length=f"{ch['L_filter']:.6f}mm",
                resonator_length=f"{ch['L_res']:.6f}mm",
                filter_tap='1.10mm', resonator_tap=f'{rt_mm:.6f}mm',
                resonator_open_extension=f'{ext_mm:.6f}mm',
                line_spacing=f'{HAMMER_LINE_SPACING_UM:g}um', tee_half='100um', lead='100um',
                jog='400um', meander_pitch='400um', fillet=FILLET_UC,
                filter_meander_height='1.05mm', resonator_meander_height='1.05mm',
                filter_first_run='back', resonator_first_run='along',
                head_width=f'{HAMMER_HEAD_WIDTH_UM:g}um', head_height=f"{ch['head_height_um']:g}um",
                head_gap=f"{ch['head_gap_um']:.1f}um",
                hammer_ground_gap='16um', corner_radius='5um', **CPW_RO,
            ))
            q_pin = np.asarray(design.components['prr' + s].pins['qubit'].middle, float)
            _route(design, 'res_open' + s, ('prr' + s, 'qubit'), (name, 'readout'),
                   [(q_pin[0], y_jog), (claw[0], y_jog)], CPW_RO, fillet=FILLET_UC)

        p8 = np.asarray(design.components['LPP8'].pins['tie'].middle, float)
        p7 = np.asarray(design.components['LPP7'].pins['tie'].middle, float)
        y_cin = p8[1] - 0.15 - cap_pin
        EnclosingCap(design, 'cin', options=dict(
            **_pos(p8[0], y_cin), orientation='-90', **ENC_OPT, **CPW_RO,
        ))
        RouteStraight(design, 'feed_cin', Dict(
            pin_inputs=Dict(start_pin=Dict(component='LPP8', pin='tie'),
                            end_pin=Dict(component='cin', pin='a')), **CPW_RO,
        ))
        cin_b = np.asarray(design.components['cin'].pins['b'].middle, float)
        yt = min(PORT_FEED['y_top'], cin_b[1] - 0.25)
        xs = PORT_FEED['x_side']
        _route(design, 'feed_in', ('cin', 'b'), ('n0_Q1', 'prime_start'),
               [(cin_b[0], yt), (-xs, yt), (-xs, feed_y)], CPW_RO, fillet=FILLET_UC)
        RouteStraight(design, 'feed_mid', Dict(
            pin_inputs=Dict(start_pin=Dict(component='n0_Q1', pin='prime_end'),
                            end_pin=Dict(component='n0_Q2', pin='prime_start')), **CPW_RO,
        ))
        _route(design, 'feed_out', ('LPP7', 'tie'), ('n0_Q2', 'prime_end'),
               [(p7[0], yt), (xs, yt), (xs, feed_y)], CPW_RO, fillet=FILLET_UC)

    for name, sx, xy_port in (('Q1', -1, 'LPP1'), ('Q2', +1, 'LPP6')):
        xe = sx * (qx + 0.128 + 0.046 + 0.080)
        OpenToGround(design, 'xy_end_' + name, options=dict(
            **_pos(xe, yq), orientation='0' if sx < 0 else '180',
            width=CPW_CTL['trace_width'], gap=CPW_CTL['trace_gap'],
        ))
        xy_y = float(design.components[xy_port].pins['tie'].middle[1])
        _route(design, 'xy_' + name, (xy_port, 'tie'), ('xy_end_' + name, 'open'),
               [(sx * 3.45, xy_y), (sx * 3.45, yq)], CPW_CTL, fillet=FILLET_UC)

    if not flux_as_ground:
        xf_tc = float(design.components['TC'].pins['flux_bias_line'].middle[0])
        xp4 = float(design.components['LPP4'].pins['tie'].middle[0])
        _route(design, 'fl_TC', ('TC', 'flux_bias_line'), ('LPP4', 'tie'),
               [(xf_tc, -3.80), (xp4, -3.80)], CPW_CTL, fillet=FILLET_UC)

    if check_keepout_um is not None:
        check_keepouts(design, min_um=float(check_keepout_um), include_site_markers=('NTHU', 'TSRI'))
    return design


def build_design(variant='sim', overrides=None, mode='surface17', channels=True,
                 flux_as_ground=False, check_keepout_um=100.0):
    if variant not in ('sim', 'fab'):
        raise ValueError(f'variant must be "sim" or "fab"; got {variant!r}')
    if mode == 'surface17':
        return _build_surface17(variant=variant, overrides=overrides, flux_as_ground=flux_as_ground)
    if mode == 'unit_cell':
        return _build_unit_cell(
            variant=variant, overrides=overrides, channels=channels,
            flux_as_ground=flux_as_ground, check_keepout_um=check_keepout_um,
        )
    raise ValueError(f'mode must be "surface17" or "unit_cell"; got {mode!r}')


def render_multilayer_gds(gds_path: Path, png_path: Path, window=None, title='', figsize=(9, 9), dpi=180):
    import gdstk
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    top, = gdstk.read_gds(str(gds_path)).top_level()
    all_polys = top.get_polygons()
    l1 = [p.points / 1000.0 for p in all_polys if p.layer == 1]
    l20 = [p.points / 1000.0 for p in all_polys if p.layer == 20]
    l21 = [p.points / 1000.0 for p in all_polys if p.layer == 21]
    l22 = [p.points / 1000.0 for p in all_polys if p.layer == 22]
    l23 = [p.points / 1000.0 for p in all_polys if p.layer == 23]
    l50 = [p.points / 1000.0 for p in all_polys if p.layer == 50]

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_facecolor('white')
    if l1:
        ax.add_collection(PolyCollection(l1, fc='#7d8a98', ec='#7d8a98', lw=0.2, zorder=1))
    if l23:
        ax.add_collection(PolyCollection(l23, fc='#581c87', ec='none', alpha=0.9, zorder=3))
    if l22:
        ax.add_collection(PolyCollection(l22, fc='#9333ea', ec='#581c87', lw=0.35, alpha=0.92, zorder=4))
    if l21:
        ax.add_collection(PolyCollection(l21, fc='#1d4ed8', ec='none', alpha=0.9, zorder=5))
    if l20:
        ax.add_collection(PolyCollection(l20, fc='#f59e0b', ec='#b45309', lw=0.35, alpha=0.95, zorder=6))
    if l50:
        ax.add_collection(PolyCollection(l50, fc='#d7322b', ec='none', zorder=7))

    from matplotlib.patches import Patch as _Patch
    handles = []
    if l20:
        handles.append(_Patch(facecolor='#f59e0b', edgecolor='#b45309', label='Crossover Airbridge (L20)'))
    if l22:
        handles.append(_Patch(facecolor='#9333ea', edgecolor='#581c87', label='Ground Airbridge (L22)'))
    if handles:
        ax.legend(handles=handles, loc='upper right', fontsize=8.5, framealpha=0.92)

    if window is None:
        bb = np.asarray(top.bounding_box()) / 1000.0
        ax.set(xlim=(bb[0, 0], bb[1, 0]), ylim=(bb[0, 1], bb[1, 1]), aspect='equal')
    else:
        ax.set(xlim=(window[0], window[1]), ylim=(window[2], window[3]), aspect='equal')

    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    if title:
        ax.set_title(title, fontsize=11, weight='bold')
    fig.tight_layout()
    fig.savefig(str(png_path), dpi=dpi, facecolor='white')
    plt.close(fig)


if __name__ == '__main__':
    from gds_tools import export_gds

    d_s17 = build_design(variant='sim', mode='surface17')
    gds_s17 = DESIGN / 'gds' / f'Design13_2026-11-27_{REVISION}_surface17_sim.gds'
    gds_s17.parent.mkdir(parents=True, exist_ok=True)
    export_gds(d_s17, str(gds_s17), fab=False, cheese=False)

    d_uc = build_design(variant='sim', mode='unit_cell', channels=True, check_keepout_um=100.0)
    gds_uc = DESIGN / 'gds' / f'Design13_2026-11-27_{REVISION}_unitcell_sim.gds'
    export_gds(d_uc, str(gds_uc), fab=False, cheese=False)

    fig_dir = DESIGN / 'fig'
    fig_dir.mkdir(parents=True, exist_ok=True)
    render_multilayer_gds(
        gds_s17, fig_dir / 'chip.png',
        title=f'Design13 {REVISION} - {DIE_MM:g}x{DIE_MM:g} mm Surface-17, A*-routed control lines (>=180 um clearance, 90-deg bridges)',
        figsize=(9.5, 9.5), dpi=180,
    )
    render_multilayer_gds(
        gds_s17, fig_dir / 'readout_bank.png', window=(0.0, DIE_HALF, 0.0, DIE_HALF),
        title=f'Design13 {REVISION} - NE Diagonal Readout Bank ({DIE_MM:g} mm die)',
        figsize=(8.5, 8.5), dpi=180,
    )
    render_multilayer_gds(
        gds_s17, fig_dir / 'core.png', window=(-0.45, 2.45, -0.75, 0.85),
        title=f'Design13 {REVISION} - Circular 46 um Gap + 585 um CPW-Neck Airbridge Unit Cell (D5 - FloatTC - Z2)',
        figsize=(9.5, 6.2), dpi=180,
    )
    render_multilayer_gds(
        gds_s17, fig_dir / 'qubit_c_arcs.png', window=(-0.42, 0.42, -0.42, 0.42),
        title=f'Design13 {REVISION} - D5: circular island, keyhole readout, 4 coupler C-arcs, XY, radial JJ',
        figsize=(8.0, 8.0), dpi=220,
    )
    render_multilayer_gds(
        gds_s17, fig_dir / 'coupler_hub.png', window=(0.84, 1.16, -0.22, 0.10),
        title=f'Design13 {REVISION} - Rounded-rect hub, rail SQUID (Manhattan JJs), straight T flux line',
        figsize=(8.0, 8.0), dpi=260,
    )
    render_multilayer_gds(
        gds_s17, fig_dir / 'crossing_zoom.png', window=(1.85, 3.65, -0.05, 1.75),
        title=f'Design13 {REVISION} - 90-deg L20 crossovers, >=180 um CPW spacing, L22 stitching',
        figsize=(8.5, 8.5), dpi=200,
    )
    render_multilayer_gds(
        gds_uc, fig_dir / 'unit_cell_chip.png',
        title=f'Design13 {REVISION} - Extracted 2Q + 1TC EM Unit Cell (10 x 10 mm, Circular Pads)',
        figsize=(8.5, 8.5), dpi=160,
    )
    s = d_s17._design13_summary
    print(f'DESIGN13_{REVISION.upper()}_STAGE1_OK:', {k: v for k, v in s.items() if k != 'resonator_lengths_mm'})
    for qn, info in list(s['resonator_lengths_mm'].items())[:6]:
        print('  RES', qn, info)
