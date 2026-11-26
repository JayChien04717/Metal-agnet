# -*- coding: utf-8 -*-
r"""Single layout entry point for Design11_Transmon2Q_Manhattan (Stage 1 simulation geometry).

2 flux-tunable Manhattan-junction Xmons (Q1, Q2) + grounded tunable coupler (TC) + two <6.5 GHz
PurcellReadout channels (|2chi|/kappa_r = 1.0) on a 10 x 10 mm 8-port symmetric die (wafer_pcb_ports).

Layout topology (per user sketch 2026-11-25):
  - Readout claws on the NORTH ('90') arm of Q1 and Q2 -> straight up into PurcellReadout meanders.
  - Feedline raised to y ~ 2.75 mm (P8 -> C_in -> n0_Q1 -> n0_Q2 -> P7).
  - XY drive lines (P1 -> Q1 west arm, P6 -> Q2 east arm) enter horizontally at y = 0.
  - Q1 / Q2 flux lines mirrored outward (Q1 enters left, Q2 enters right) with 45-degree fan-out to P2 / P5.
  - TC flux line goes straight down to y = -3.90 mm -> P4 (P3 grounded spare).
"""
from collections import OrderedDict
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'code' / 'general'))
import project_paths
project_paths.configure(__file__)

COMPONENTS_DIR = Path(__file__).resolve().parent / 'components'
if str(COMPONENTS_DIR) not in sys.path:
    sys.path.insert(0, str(COMPONENTS_DIR))

import numpy as np
from qiskit_metal import Dict, designs
from qiskit_metal.qlibrary.core.qroute import QRoute
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight

from bondpad import add_bond_pad
from check_keepouts import check_keepouts
from claw_coupler_manhattan import ClawCouplerManhattan
from enclosing_cap import EnclosingCap
from fluxline import B2_FLUX_OPTIONS
from purcell_readout import PurcellReadout
from round_corners import round_component_polys, square_merged_etch_corners
from tap_junction import TapJunction
from wafer_pcbport import wafer_pcb_ports
from xmon_manhattan import XmonManhattan

CPW_RO = dict(trace_width='10um', trace_gap='6um')
CPW_CTL = dict(trace_width='10um', trace_gap='6um')
FILLET = '99um'

# Core geometry r8 (10/6 CPW, 14 um Xmon pocket + 3 um ground-strip readout claw, TC claw_sep=12.5 um;
# Design11 local min metal width override: > 2 um per user fab verification >= 1.5 um)
QX = 0.40       # mm, Q1 at -QX, Q2 at +QX
YQ = 0.0        # mm
XMON_OPT = dict(
    cross_length='142um', cross_width='28um', cross_gap='14um',
)
# r9 (readout option B, |2chi| = kappa_r = 3 MHz at fixed E_C): no ground strip between the Xmon pocket and the
# readout claw (claw_sep 18 um - cross_gap 14 um < claw_gap 6 um -> claw etch merges with the pocket); C_qr is set
# per qubit by claw_length, C_sigma restored by cross_length (Q3D core sweeps optB_sweepA..E).
RO_CLAW = dict(
    claw_length='122um', claw_width='10um', claw_sep='18um', claw_gap='6um',
    cpw_width='10um', cpw_gap='6um', cpw_length='40um',
)
# Per-qubit overrides merged on top of XMON_OPT / RO_CLAW (readout option B retarget needs different C_qr per qubit
# at fixed E_C). Q1/Q2 cross_length may differ; TC length and centre follow both inner arm tips.
XMON_BY_QUBIT = dict(Q1=dict(cross_length='143.4um'), Q2=dict(cross_length='143.5um'))
RO_CLAW_BY_QUBIT = dict(Q1=dict(claw_length='100um'), Q2=dict(claw_length='114um'))
COUPLER_OPT = dict(
    pad_arm_width='15um', paddle_width='30um', paddle_height='78um', paddle_gap='55um',
    pad_head_width='15um', pad_head_length='285um', pad_gap='30um', jj_gap='48um',
    claw_width='10um', claw_length='17um', claw_sep='12.5um', claw_gap='12um',
    squid_spacing='30um', squid_bar_height='6um', squid_contact_overlap='3um',
    flux_bias_line_options=Dict(B2_FLUX_OPTIONS, make_fbl=True),
    charge_line_options=Dict(make_cl=False),
    readout_line_options=Dict(make_rol=False),
)
# Per-branch junction inductances (TC symmetric SQUID: two 14.03 nH branches in parallel -> 7.015 nH)
LJ = dict(Q1a='17.39nH', Q1b='34.79nH', Q2a='20.55nH', Q2b='41.11nH', TC='14.03nH')
XY_SEP_MM = 0.14  # horizontal separation from outer cross arm etch to XY open end

# PurcellReadout channels (10/6 CPW, option B: kappa_r = |2chi| = 3 MHz, kappa_p = 30 MHz (user 2026-11-26, was 20),
# J = 4.74 MHz). r16: hammer line spacing 200 um (fixed rule), head 20 x 40 um, head_gap and lengths from
# optB_channel_solve 20261126_091830_optB_channel_solve_4e2b6737 (200 um coupons q3d_coupons_cj_s200_h40_1be966a4,
# r15 HFSS - model offsets carried over from 20261126_090714_optB_channel_solve_f5449433, kp-correction 0.718 / 2.093);
# C_p plates kept at r15. r15 = r14 with L_res from simulation/resonator/optB_channel_solve.py --kappa-p 30 --kp-correction (Q3D coupon
# C_j / C_p tables incl. q3d_coupons_cp_r13_9a35735e, Kron claw end load, full-feed kappa_p model scaled per channel
# by the HFSS-fitted / model kappa_p of r13: Q1 0.723, Q2 2.121) with the per-line anchor on the complex-fit poles of
# the r13 3-pass S21 (fit run 20261126_071122_fit_s21_r13_kp30_3p_fine4_30189ea9; solve run
# 20261126_071206_optB_channel_solve_58ba120b). tap_r = resonator_tap + open extension = 1.735 mm.
# History (L_filter / L_res mm, head_gap um, cp_plate um):
#   r9  (kappa_p 20): Q1 4.2420 / 4.3645, 39.8, 168.5; Q2 4.4033 / 4.4957, 38.5, 160.0
#   r10 (kappa_p 30, r9 anchor): Q1 4.1639 / 4.3410, 34.7, 218.5; Q2 4.3205 / 4.4402, 33.5, 209.5
#   r11 (r10 anchor): Q1 4.1599 / 4.3433; Q2 4.3217 / 4.4661 (r10 gaps / plates)
#   r12 (r11 anchor): Q1 4.1589 / 4.3484; Q2 4.3234 / 4.4599 (r10 gaps / plates); HFSS kappa_p 22.9 / 59.8 MHz
#   r13 (r12 fit anchor, kp-correction 0.762 / 1.994): Q1 4.1275 / 4.3552, 34.4, 259.0; Q2 4.4096 / 4.4484, 34.2, 129.0;
#       HFSS kappa_p 28.5 / 31.9, f_r - f_p -18.8 / -2.1 MHz
#   r14 (r13 fit anchor, kp-correction 0.723 / 2.121): Q1 4.1142 / 4.3385, 34.3, 267.0; Q2 4.4246 / 4.4543, 34.2, 123.0;
#       HFSS (fit 20261126_074256_..._705e5a1f) f_p 6.3496 / 6.1506, f_r 6.3648 / 6.1581, kappa_p 29.7 / 29.8, J 4.89 / 4.81
#   r15: filters, gaps and plates kept at r14 (on target); L_res from solve 20261126_074353_optB_channel_solve_d30bb89c
#       (r14 fit anchor, kp-correction 0.716 / 2.107); 400 um spacing. HFSS fit: f_r - f_p -0.9 / -1.8 MHz, kappa_p 29.8 / 29.6,
#       J 4.69 / 4.73 MHz (accepted 2026-11-26, superseded by the 200 um spacing rule)
#   r15 values: Q1 4.1142 / 4.3490, 34.3, 267.0; Q2 4.4246 / 4.4604, 34.2, 123.0 (head 20 x 40)
#   r16 (200 um spacing, r15 offsets carried over): Q1 4.1828 / 4.4104, 38.9, 267.0; Q2 4.4984 / 4.5209, 38.8, 123.0;
#       HFSS eigen pass 2 ~6.344/6.358 (Q1), ~6.143/6.160 (Q2); driven f_r not converged (skill/em_simulation.md section 4)
#   r17: Q3D-only solve 20261126_155318_q3d_length_solve_11ac23ec (no HFSS offsets): straight-CPW Q3D line constants
#       (eps_eff 6.365) + full-fixture in-situ Q3D loads and C_j (model check 20261126_144132_..._fa91743e); head_gap from
#       the in-situ C_j scaled by the 200 um / 40 um coupon C_j(gap) shape; C_p plates kept (driven kappa_p 28.6 / 30.6).
# Hammer: filter / resonator main lines fixed at 200 um centre to centre (skill/purcell_resonator_design.md section 1,
# user 2026-11-26; r9-r15 used 400 um). J is tuned only with head_height / head_gap; head_width fixed.
HAMMER_LINE_SPACING_UM = 200.0
HAMMER_HEAD_WIDTH_UM = 20.0
PURCELL_CHANNELS = {
    'Q1': dict(side=-1, f_GHz=6.35, L_filter=4.1577, L_res=4.4060, head_height_um=40.0, head_gap_um=45.2, cp_plate_um=267.0),
    'Q2': dict(side=+1, f_GHz=6.15, L_filter=4.5024, L_res=4.5370, head_height_um=40.0, head_gap_um=44.6, cp_plate_um=123.0),
}
ENC_OPT = dict(
    plate_width='140um', plate_length='41um', gap='8um',
    wall='52um', wall_side='57um', leg_ext='72um', lead_length='40um', ground_gap='12um',
)
PORT_FEED = dict(y_top=3.55, x_side=2.75)


class RoutePolyline(QRoute):
    """Exact waypoint polyline route (supports 45-degree diagonals with native fillet)."""
    default_options = Dict(anchors=OrderedDict())

    def make(self):
        self.set_pin('start')
        self.set_pin('end')
        p = self.p
        sp = self.design.components[p.pin_inputs.start_pin.component].pins[p.pin_inputs.start_pin.pin]
        ep = self.design.components[p.pin_inputs.end_pin.component].pins[p.pin_inputs.end_pin.pin]
        pts = ([np.asarray(sp.middle, dtype=float)]
               + [np.asarray(v, dtype=float) for v in self.options.anchors.values()]
               + [np.asarray(ep.middle, dtype=float)])
        self.make_elements(np.array(pts, dtype=float))


def _pos(x, y):
    return dict(pos_x=f'{x:.6f}mm', pos_y=f'{y:.6f}mm')


def _route(design, name, start, end, anchors, cpw):
    return RouteAnchors(design, name, Dict(
        chip='main', fillet=FILLET, hfss_wire_bonds=False,
        lead=Dict(start_straight='0um', end_straight='0um'),
        anchors=OrderedDict({i: np.array(a, float) for i, a in enumerate(anchors)}),
        pin_inputs=Dict(start_pin=Dict(component=start[0], pin=start[1]),
                        end_pin=Dict(component=end[0], pin=end[1])),
        **cpw,
    ))


def _polyroute(design, name, start, end, anchors, cpw):
    return RoutePolyline(design, name, Dict(
        chip='main', fillet=FILLET, hfss_wire_bonds=False,
        anchors=OrderedDict({i: np.array(a, float) for i, a in enumerate(anchors)}),
        pin_inputs=Dict(start_pin=Dict(component=start[0], pin=start[1]),
                        end_pin=Dict(component=end[0], pin=end[1])),
        **cpw,
    ))


def _poly_geom(design, comp_name, poly_name):
    t = design.qgeometry.tables['poly']
    row = t[(t['component'] == design.components[comp_name].id) & (t['name'] == poly_name)]
    return row.iloc[0]['geometry'] if len(row) else None


def build_design(variant='sim', overrides=None, channels=True, flux_as_ground=False, check_keepout_um=100.0):
    """Build Design11_Transmon2Q_Manhattan (Stage 1 simulation geometry)."""
    if variant != 'sim':
        raise ValueError(f'Stage 1/2 supports variant="sim" only; got {variant!r}')
    ov = overrides or {}
    xmon_opt = dict(XMON_OPT, **ov.get('xmon', {}))
    ro_claw = dict(RO_CLAW, **ov.get('ro_claw', {}))
    coupler_opt = dict(COUPLER_OPT, **ov.get('coupler', {}))
    qx = float(ov.get('qx_mm', QX))
    fag = 'True' if flux_as_ground else 'False'

    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon_eps11p45'
    design.chips.main.size.update(dict(size_x='10mm', size_y='10mm', size_z='-500um', center_x='0mm', center_y='0mm'))
    design.variables.update(cpw_width=CPW_RO['trace_width'], cpw_gap=CPW_RO['trace_gap'])

    # Qubits: Q1 left (-qx), Q2 right (+qx); readout claw on NORTH ('90') arm;
    # flux lines mirrored outward (Q1: mirror=True -> enters left; Q2: mirror=False -> enters right)
    xmon_q, arms = {}, {}
    for name, sx in (('Q1', -1), ('Q2', +1)):
        fbl = Dict(B2_FLUX_OPTIONS, make_fbl=True, mirror=str(sx < 0))
        q_ro = {**RO_CLAW, **RO_CLAW_BY_QUBIT[name], **ov.get('ro_claw', {}), **ov.get('ro_claw_' + name, {})}
        xmon_q[name] = {**xmon_opt, **XMON_BY_QUBIT[name], **ov.get('xmon_' + name, {})}
        arms[name] = design.parse_value(xmon_q[name]['cross_length'])
        XmonManhattan(design, name, options=dict(
            **_pos(sx * qx, YQ), **xmon_q[name],
            Lj_a=LJ[name + 'a'], Lj_b=LJ[name + 'b'],
            jj_metal='True', flux_as_ground=fag,
            flux_bias_line_options=fbl,
            connection_pads=dict(readout=dict(q_ro, connector_location='90')),
        ))
        raw_etcher = _poly_geom(design, name, 'readout_connector_etcher')
        round_component_polys(design.components[name], 5e-3)
        # r15: the claw etch merges with the Xmon pocket (no ground strip since r9); its R = 5 um rounded corners
        # there left 0.58 um ground cusps (GDS check 2026-11-25). Restore the square etch corner within 5 um of the
        # pocket only (about 3 um^2 of ground removed per corner; no Q3D rerun, user 2026-11-26).
        if ov.get('square_claw_corners', True):   # override False: keep the plain R = 5 um rounding (r9-r15)
            square_merged_etch_corners(design.components[name], 'readout_connector_etcher', raw_etcher,
                                       _poly_geom(design, name, 'cross_etch'), pad=5e-3)

    # Grounded tunable coupler (TC) spanning between Q1 and Q2 inner arms
    # Inner arm tips at -qx + arm_Q1 and qx - arm_Q2 -> TC bar spans between them, centred on their midpoint
    claw_sep = design.parse_value(coupler_opt['claw_sep'])
    L_bar = (2.0 * qx - arms['Q1'] - arms['Q2']) - 2.0 * claw_sep
    x_tc = 0.5 * (arms['Q1'] - arms['Q2'])
    ClawCouplerManhattan(design, 'TC', options=dict(
        **_pos(x_tc, YQ), pad_arm_length=f'{L_bar * 1e3:.4f}um',
        qubit_arm_width=xmon_opt['cross_width'], Lj=LJ['TC'],
        jj_metal='True', flux_as_ground=fag, **coupler_opt,
    ))
    round_component_polys(design.components['TC'], 5e-3)

    # 8 symmetric wafer launch pads (4 mm pitch on all four sides)
    for pn in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'):
        if not channels and pn in ('P7', 'P8'):
            continue
        cpw = CPW_RO if pn in ('P7', 'P8') else CPW_CTL
        add_bond_pad(design, 'LP' + pn, wafer_pcb_ports[pn], **cpw)

    # Readout / Purcell channels (feedline raised to y ~ 2.75 mm; resonator drops straight into north claw)
    if channels:
        cap_pin = (41.0 / 2.0 + 72.0 + 40.0) / 1000.0
        fillet_loss = design.parse_value(FILLET) * (2.0 - np.pi / 2.0)
        claw_q1 = np.asarray(design.components['Q1'].pins['readout'].middle, float)
        # Filter line and its feed tee stay at xf = side * 1.80 (feed standing wave / kappa_p unchanged);
        # resonator line at xr = xf - HAMMER_LINE_SPACING (1.60 mm from r16; 1.40 mm with the 400 um spacing of r9-r15).
        # Open extension: straight down from qubit pin to y_jog, jog horizontally to claw_x, straight down to claw
        xf_abs = 1.80
        xr_abs = xf_abs - HAMMER_LINE_SPACING_UM / 1000.0
        y_jog = claw_q1[1] + 0.25
        horiz = abs(xr_abs - abs(claw_q1[0]))
        rt_mm = 0.32
        # Let qubit pin sit at y_jog + 0.25 -> total extension = 0.25 + horiz + 0.25 - 2*fillet_loss
        # Adjust y_pin_target so total extension = 1.30 - rt_mm = 0.98 mm
        ext_mm = 1.30 - rt_mm
        dy_top = ext_mm - horiz - 0.25 + 2.0 * fillet_loss
        if dy_top < 0.22:
            # Put qubit pin directly above claw_x? No, raise y_pin_target slightly and adjust rt_mm
            dy_top = 0.25
            ext_mm = dy_top + horiz + 0.25 - 2.0 * fillet_loss
            rt_mm = max(0.32, 1.30 - ext_mm)
        y_pin_target = y_jog + dy_top
        y_tap = y_pin_target + rt_mm
        feed_y = y_tap + 1.20 + 0.150 + 2.0 * cap_pin

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
            round_component_polys(design.components['cp' + s], 5e-3, skip_names=())
            RouteStraight(design, 'stub_p' + s, Dict(
                pin_inputs=Dict(start_pin=Dict(component='n0' + s, pin='second_end'),
                                end_pin=Dict(component='cp' + s, pin='a')), **CPW_RO,
            ))
            cp_b = np.asarray(design.components['cp' + s].pins['b'].middle, float)
            y_tap_ch = cp_b[1] - 1.20
            x0 = 0.5 * (xf + xr)
            PurcellReadout(design, 'prr' + s, options=dict(
                **_pos(x0, y_tap_ch), orientation='0', mirror=str(side > 0),
                filter_length=f"{ch['L_filter']:.6f}mm",
                resonator_length=f"{ch['L_res']:.6f}mm",
                filter_tap='1.20mm', resonator_tap=f'{rt_mm:.6f}mm',
                resonator_open_extension=f'{ext_mm:.6f}mm',
                line_spacing=f'{HAMMER_LINE_SPACING_UM:g}um', tee_half='100um', lead='100um',
                jog='400um', meander_pitch='400um', fillet=FILLET,
                filter_meander_height='1.10mm', resonator_meander_height='1.10mm',
                filter_first_run='back', resonator_first_run='along',
                head_width=f'{HAMMER_HEAD_WIDTH_UM:g}um', head_height=f"{ch['head_height_um']:g}um", head_gap=f"{ch['head_gap_um']:.1f}um",
                hammer_ground_gap='16um', corner_radius='5um', **CPW_RO,
            ))
            q_pin = np.asarray(design.components['prr' + s].pins['qubit'].middle, float)
            _route(design, 'res_open' + s, ('prr' + s, 'qubit'), (name, 'readout'),
                   [(q_pin[0], y_jog), (claw[0], y_jog)], CPW_RO)

        # Feedline P8 -> C_in -> n0_Q1 -> n0_Q2 -> P7
        p8 = np.asarray(design.components['LPP8'].pins['tie'].middle, float)
        p7 = np.asarray(design.components['LPP7'].pins['tie'].middle, float)
        y_cin = p8[1] - 0.15 - cap_pin
        EnclosingCap(design, 'cin', options=dict(
            **_pos(p8[0], y_cin), orientation='-90', **ENC_OPT, **CPW_RO,
        ))
        round_component_polys(design.components['cin'], 5e-3, skip_names=())
        RouteStraight(design, 'feed_cin', Dict(
            pin_inputs=Dict(start_pin=Dict(component='LPP8', pin='tie'),
                            end_pin=Dict(component='cin', pin='a')), **CPW_RO,
        ))
        cin_b = np.asarray(design.components['cin'].pins['b'].middle, float)
        yt = min(PORT_FEED['y_top'], cin_b[1] - 0.25)
        xs = PORT_FEED['x_side']
        _route(design, 'feed_in', ('cin', 'b'), ('n0_Q1', 'prime_start'),
               [(cin_b[0], yt), (-xs, yt), (-xs, feed_y)], CPW_RO)
        RouteStraight(design, 'feed_mid', Dict(
            pin_inputs=Dict(start_pin=Dict(component='n0_Q1', pin='prime_end'),
                            end_pin=Dict(component='n0_Q2', pin='prime_start')), **CPW_RO,
        ))
        _route(design, 'feed_out', ('LPP7', 'tie'), ('n0_Q2', 'prime_end'),
               [(p7[0], yt), (xs, yt), (xs, feed_y)], CPW_RO)

    # Control lines:
    # - XY lines (P1 -> Q1 west arm, P6 -> Q2 east arm) enter horizontally at y = 0
    # - Q1 / Q2 flux lines fan out at 45 degrees down to y = -2.0 mm -> P2 / P5
    for name, sx, xy_port, fl_port in (('Q1', -1, 'LPP1', 'LPP2'), ('Q2', +1, 'LPP6', 'LPP5')):
        # Horizontal XY open end facing outer cross arm tip at y = YQ
        cg = design.parse_value(xmon_q[name]['cross_gap'])
        xe = sx * (qx + arms[name] + cg + XY_SEP_MM)
        OpenToGround(design, 'xy_end_' + name, options=dict(
            **_pos(xe, YQ), orientation='0' if sx < 0 else '180',
            width=CPW_CTL['trace_width'], gap=CPW_CTL['trace_gap'],
        ))
        xy_y = float(design.components[xy_port].pins['tie'].middle[1])
        _route(design, 'xy_' + name, (xy_port, 'tie'), ('xy_end_' + name, 'open'),
               [(sx * 3.45, xy_y), (sx * 3.45, YQ)], CPW_CTL)

        # 45-degree diagonal flux fan-out to P2 / P5
        xf, yf = [float(v) for v in design.components[name].pins['flux_bias_line'].middle]
        fl_y = float(design.components[fl_port].pins['tie'].middle[1])
        y_bend = yf - 0.20
        dx_diag = abs(y_bend - fl_y)
        _polyroute(design, 'fl_' + name, (name, 'flux_bias_line'), (fl_port, 'tie'),
                   [(xf, y_bend), (xf + sx * dx_diag, fl_y)], CPW_CTL)

    xf_tc = float(design.components['TC'].pins['flux_bias_line'].middle[0])
    xp4 = float(design.components['LPP4'].pins['tie'].middle[0])
    _route(design, 'fl_TC', ('TC', 'flux_bias_line'), ('LPP4', 'tie'),
           [(xf_tc, -3.90), (xp4, -3.90)], CPW_CTL)

    if check_keepout_um is not None:
        check_keepouts(design, min_um=float(check_keepout_um), include_site_markers=('NTHU', 'TSRI'))
    return design


if __name__ == '__main__':
    from gds_tools import export_gds, plot_gds
    d = build_design(variant='sim', channels=True, check_keepout_um=100.0)
    gds_out = Path(__file__).resolve().parents[1] / 'gds' / 'Design11_2026-11-26_r17_sim.gds'
    fig_out = Path(__file__).resolve().parents[1] / 'fig' / 'chip.png'
    fig_core = Path(__file__).resolve().parents[1] / 'fig' / 'core.png'
    fig_squid = Path(__file__).resolve().parents[1] / 'fig' / 'squid_q1.png'
    gds_out.parent.mkdir(parents=True, exist_ok=True)
    fig_out.parent.mkdir(parents=True, exist_ok=True)
    export_gds(d, str(gds_out), fab=False, cheese=False)
    plot_gds(str(gds_out), str(fig_out), ports=wafer_pcb_ports)
    plot_gds(str(gds_out), str(fig_core), window=(-0.75, 0.75, -0.45, 0.45), dpi=160, figsize=(10, 6))
    plot_gds(str(gds_out), str(fig_squid), window=(-0.45, -0.35, -0.24, -0.14), dpi=160, figsize=(6, 6))
    print('BUILD_AND_FIG_OK:', gds_out, fig_out, fig_core, fig_squid)
