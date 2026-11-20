"""Design6: three flux-tunable transmons - two Xmon qubits (asymmetric SQUID) + one grounded SQUID tunable coupler
(CZ gate) - each qubit read out through its own Purcell filter (Design4 / Cai 2024 channel), 10 x 10 mm chip.

    Q1 (Xmon, SQUID) -- claw -- TC (grounded SQUID transmon, claws on both ends) -- claw -- Q2 (Xmon, SQUID)

Ports (PCB positions of Design1-3, pattern design/pcbport.py):
    P8 -> P7 : readout feedline (top). Two hanger channels (Design4 channel: n0 tee - C_p enclosing cap - lambda/4
               filter, hammer C_j at the 1.2 / 1.3 mm taps - lambda/4 resonator whose open end is the qubit claw)
    P1 / P6  : XY (charge) lines of Q1 / Q2
    P2 / P5  : flux lines of Q1 / Q2
    P4       : TC flux line
    P3       : grounded spare on wafer layouts; no TC charge line

Readout claws sit on the OUTER arms (Q1 west, Q2 east) so that the resonator open run keeps the Design4
1.3 mm tap distance; coupler claws on the inner arms; SQUID + L flux line on the south arm
(flux leg toward the coupler side); XY line end in the outer-south quadrant.

Usage:
    from design6_layout import build_design
    d = build_design()                     # simulation geometry (JJs in the junction table, Lj in HFSS)
    d = build_design(jj_metal=False)       # mask (junction-free) version
    python design6_layout.py               # GDS (gds/Design6_<date>_{sim,fab}.gds), GUI screenshots, checks
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I, glob_files as _glob_files, glob_paths as _glob_paths
_configure(__file__)
from project_paths import pattern_design_paths as _pattern_paths
import glob
import json
import os
import sys
import time
import warnings
from collections import OrderedDict

import numpy as np

warnings.filterwarnings('ignore')
HERE = _S('Design6/layout')
ROOT = _S('Design6')                       # Design6
REPO = _S('.')
sys.path[:0] = [HERE, _I('Design4/layout'), *_pattern_paths()]
import purcell_6500 as P          # noqa: E402  (sets SI_EPS = 11.5 before design5_layout is imported)
import design5_layout as d5       # noqa: E402
from qiskit_metal import Dict      # noqa: E402
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround   # noqa: E402
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors  # noqa: E402
from wafer_pcbport import wafer_pcb_ports      # noqa: E402
from bondpad import LP, add_bond_pad      # noqa: E402
from fluxline import B2_FLUX_OPTIONS
from xmon import Xmon              # noqa: E402
from claw_coupler import ClawCoupler, ClawCouplerFab   # noqa: E402

TAG = f"Design6_{time.strftime('%Y-%m-%d')}"
GDS_DIR, FIG_DIR, PAR_DIR = (_relocate(os.path.join(ROOT, k)) for k in ('gds', 'fig', 'parameters'))

# ---------------------------------------------------------------------------------------------
# Knobs (mm unless noted). Every number used by the simulations lives here.
# ---------------------------------------------------------------------------------------------
FILLET = '99um'
CPW_RO = dict(trace_width='16um', trace_gap='9um')        # readout system (Design4 / Cai: 16 / 9 um)
CPW_CTL = dict(trace_width='10um', trace_gap='6um')       # control lines (Design1-3)
QX, YQ = 0.45, 0.0          # Q1 at (-QX, YQ), Q2 at (+QX, YQ), coupler bar at y = YQ
# Hamiltonian targets (simulation/cz, 2026-11-14): Q1 / Q2 idle at their upper sweet spots 4.80 / 4.40 GHz (asymmetric
# SQUIDs, EJa:EJb = 2:1, f_min ~ 0.55 f_max), Ec 0.20 GHz; coupler f_max 6.3 GHz (idle, sweet spot), Ec_c 0.21 GHz;
# g_qc(n) ~ 57 MHz (C_qc 3.3 fF, Q3D v4). Idle Delta = 400 MHz -> small static ZZ; CZ: Q1 flux pulse to the
# |11> <-> |02> condition (~4.18 GHz) + coupler pulse. Readout g ~ 128 MHz (C_qr ~ 8 fF), chi ~ 0.7 MHz.
XMON = dict(cross_width='30um', cross_length='185um', cross_gap='30um')
# Junction part = Design3 / Design3_2 TC SQUID (my_tmon.py) for all three SQUIDs: rectangular loop 20 x 24 um, 40 um
# block, 5 um cut, jj_width 7 um, two parallel Dolan JJs. Simulation geometry: electrodes flush with their base metal
# (Q3D rejects the overlap: 'squid_ground overlaps ground_main_plane'); mask: SQUID_MASK (3 um overlap, shifted 2 um
# toward the arm / stem -> 1 um left on the ground side), landing electrodes kept, fab mask with contact notches.
QSQUID = dict(make='True', jj_gap='40um', squid_width='40um', loop_width='20um', loop_height='24um', jj_cut='5um',
              contact_overlap='0um', shift_toward_arm='0um')
SQUID_MASK = dict(overlap='3um', shift='2um')
# flux lines: variant B of the 2026-11-14 study (README 7.3), shortened 2026-11-15 (B2): L type, sep 9 um (3 um ground
# strip), bar 45 um; shorted end at the SQUID's outer edge (x = -20.5 um), feed leg just outside the other edge (34.5-44.5 um)
QFLUX = dict(B2_FLUX_OPTIONS, make_fbl='True')
# junction inductances (fab targets), re-targeted from the HFSS Ec (hfss_core_v8 fit: 0.185 / 0.188 / 0.185 GHz) for
# f_max Q1 4.80, Q2 4.40, TC 6.30 GHz (cz_scqubits.py --retarget); asymmetric 2:1 (L_b = 2 L_a). The HFSS core sweep
# and readout runs of 2026-11-14 used the earlier Q3D-based values 15.89 / 31.77, 18.78 / 37.56, 12.85 nH.
LJ = dict(Q1a='14.56nH', Q1b='29.11nH', Q2a='17.21nH', Q2b='34.41nH', TC='11.64nH')
READOUT_CLAW = dict(claw_length='55um', claw_width='10um', claw_sep='8um', claw_gap='12um',    # xmon.py claw
                    cpw_width='16um', cpw_gap='9um', cpw_length='40um')
COUPLER = dict(
    pad_arm_width='15um', paddle_gap='30um',               # bar and its etch
    pad_head_width='15um', pad_head_length='180um', pad_gap='30um', jj_gap='40um',
    claw_width='10um', claw_length='10um', claw_sep='15um', claw_gap='12um',
    squid_loop_shape='rectangle', squid_loop_width='20um', squid_loop_height='24um',   # Design3_2 SQUID (M ~ 1.07 pH)
    squid_contact_overlap='0um', squid_shift_toward_stem='0um', gds_cell_name='junction_t0',   # mask: SQUID_MASK
    flux_bias_line_options=Dict(B2_FLUX_OPTIONS, make_fbl=True),
    # no TC charge line (2026-11-14: its stub coupled 3.3x more to Q1 than to TC); coupler spectroscopy via flux line
    charge_line_options=Dict(make_cl=False),
    readout_line_options=Dict(make_rol=False),
)
XY_END = dict(dx=0.16, sep=0.25)   # XY open end: dx outward of the qubit centre, sep below the arm etch (outer-south quadrant)

# readout channels (Design4 channel; f_p = f_r)
READOUT = OrderedDict([('Q1', dict(f=7.00, side=-1)), ('Q2', dict(f=6.80, side=+1))])   # >= 0.5 GHz above the coupler f_max
HFSS_CORR = 6.08 / 5.975         # Design4 eigenmode vs circuit model (chip6_d4)
CLAW_DL = 0.100                  # mm, resonator shortening for the claw end load (C_end ~ 17 fF / 167 fF/mm); HFSS trims it
# A simulation WITHOUT the claw must use L_res + CLAW_DL (skill/purcell_resonator_design.md section 3).
HAMMER = dict(P.PAPER_HAMMER, gap_um=45)      # 25 x 65 um plates, gap 45 um: Q3D C_j ~ 0.60 fF (Design4 q3d400 sweep)
CH_SPACING = P.MIN_CPW_SPACING   # filter - resonator centre lines 400 um
JR_X = 1.30                      # |x| of the resonator tap jr; jf at |x| = JR_X + 0.4 (filter outside)
RES_HEIGHT = 0.90                # resonator serpentine run height (2 steps, ends above the qubit)
PORT_FEED = dict(y_top=3.40, x_side=2.90)
MARK = dict(arm_length='250um', arm_width='10um', window_size='400um', label_size='150um', label_gap='30um')   # LABEL_SIZE_UM in xmark.py
TEST_JJ = dict(x=3.95, y=4.25, pocket=(0.800, 0.300), pad_left=(0.225, 0.129), pad_right=(0.295, 0.129), pad_gap=0.211)
DIE_ID = 'TXT00'     # lower right, like Design3 FTF00 (any ID is safe since the gds_tools union fix of 2026-11-15)


def _length_file():
    files = sorted(_glob_files(_relocate(os.path.join(PAR_DIR, 'Design6_*.json'))))
    return files[-1] if files else None


def tune_lengths(save=True):
    """Filter / resonator lengths for each channel: purcell_6500.tune_lengths (Design4 circuit model, Q3D caps)
    at f * HFSS_CORR, resonator shortened by CLAW_DL for the claw. Returns {name: dict(L_filter, L_res, J, kp)}."""
    out = {}
    for name, ch in READOUT.items():
        P.F0 = ch['f'] * HFSS_CORR
        p, r = P.tune_lengths()
        out[name] = dict(f_target_GHz=ch['f'], f_model_GHz=P.F0, L_filter=float(p['L_filt']),
                         L_res=float(p['L_res']) - CLAW_DL, J_model_MHz=float(r['J']), kappa_p_model_MHz=float(r['kp']))
        print(f"{name} {ch['f']:.2f} GHz: L_filter {out[name]['L_filter']:.4f}, L_res {out[name]['L_res']:.4f} mm "
              f"(model J {r['J']:.1f}, kappa_p {r['kp']:.1f} MHz at C_j of purcell_6500)", flush=True)
    return out


def lengths_from_params():
    f = _length_file()
    if f:
        par = json.load(open(f))
        if 'channels' in par:
            return {c['name']: dict(L_filter=c['L_filter_mm'], L_res=c['L_res_mm']) for c in par['channels']}
    return tune_lengths()


def route(design, name, start, end, anchors, cpw, ss='0um', es='0um'):
    return RouteAnchors(design, name, Dict(
        chip='main', fillet=FILLET, lead=Dict(start_straight=ss, end_straight=es),
        anchors=OrderedDict({i: np.array(a, float) for i, a in enumerate(anchors)}),
        pin_inputs=Dict(start_pin=Dict(component=start[0], pin=start[1]),
                        end_pin=Dict(component=end[0], pin=end[1])), **cpw))


def pos(x, y):
    return dict(pos_x=f'{x:.9f}mm', pos_y=f'{y:.9f}mm')


def add_channel(design, name, side, Lf, Lr, feed_y, claw_pin, readout_hammer=None):
    """Design4 channel hanging below the feed line at n0 = (jf_x, feed_y). side -1: Q1 (filter on the left),
    +1: Q2 (mirrored in x). The resonator open run goes from jr straight down to the claw height, then
    horizontally into the readout claw of the qubit (its open end); its filleted length is RES_TAP - TEE_HALF."""
    s = '_' + name
    H = HAMMER if readout_hammer is None else readout_hammer
    spacing = H.get('cpw_centres_um', CH_SPACING * 1000) / 1000
    jr_x = side * JR_X
    jf_x = jr_x + side * spacing
    x0 = jf_x
    d5.TapJunction(design, 'n0' + s, options=dict(**CPW_RO, **pos(x0, feed_y), orientation='0',
                                                     stub_length='100um', stub_up='False'))
    opts = {k: f'{v}um' for k, v in dict(d5.ENC, **P.ENC_CAI).items()}
    d5.EnclosingCap(design, 'cp' + s, options=dict(opts, **CPW_RO, **pos(x0, feed_y - P.CP_Y), orientation='-90'))
    P.straight(design, 'stub_p' + s, ('n0' + s, 'second_end'), ('cp' + s, 'a'))
    stub = spacing / 2 - (H['gap_um'] / 2 + H['width_um'] + d5.CJ_LEAD) / 1000 - 0.02
    if stub <= 0:
        raise ValueError('Hammer geometry does not fit the selected CPW spacing')
    y_tap = feed_y - P.TAP_Y
    # orientation 90: prime_start points down (away from the feed), prime_end up; stub_up False -> +x
    d5.TapJunction(design, 'jf' + s, options=dict(**CPW_RO, **pos(jf_x, y_tap), orientation='90',
                                                     stub_up=str(side > 0), stub_length=f'{stub:.9f}mm'))
    d5.TapJunction(design, 'jr' + s, options=dict(**CPW_RO, **pos(jr_x, y_tap), orientation='90',
                                                     stub_up=str(side < 0), stub_length=f'{stub:.9f}mm'))
    d5.HammerCap(design, 'cj' + s, options=dict(**CPW_RO, **pos((jf_x + jr_x) / 2, y_tap),
                                                head_width=f"{H['width_um']}um", head_height=f"{H['height_um']}um",
                                                head_gap=f"{H['gap_um']}um", lead_length=f'{d5.CJ_LEAD}um',
                                                ground_gap=f"{H['ground_um']}um"))
    a_f, a_r = ('a', 'b') if side < 0 else ('b', 'a')        # hammer pin on the filter / resonator side
    P.straight(design, 'stub_f' + s, ('jf' + s, 'second_end'), ('cj' + s, a_f))
    P.straight(design, 'stub_r' + s, ('cj' + s, a_r), ('jr' + s, 'second_end'))
    P.straight(design, 'feed_jf' + s, ('cp' + s, 'b'), ('jf' + s, 'prime_end'))

    # resonator open run: jr -> down to the claw height -> horizontally into the claw (pin normal points outward)
    route(design, 'res_open' + s, ('jr' + s, 'prime_start'), (name, 'readout'),
          [(jr_x, claw_pin[1])], CPW_RO)
    open_len = d5.f3.path_length(design, 'res_open' + s)
    target = P.RES_TAP - P.TEE_HALF
    if abs(open_len - target) > 0.002:
        raise ValueError(f'res_open{s}: {open_len:.4f} mm vs tap target {target:.4f} mm (move the feed line)')

    # filter tail: lead down past the tap (length absorbs the rest), one step outward, run up to feed - 0.65
    y_bottom = P.WIDE['y_bottom']
    Lft = Lf - P.FILTER_TAP - P.TEE_HALF
    loss = design.parse_value(FILLET) * (2 - np.pi / 2)
    top = feed_y - y_bottom
    start_y = y_tap - P.TEE_HALF
    # Lft = lead + step + (top - (start_y - lead)) - 2 loss  ->  lead
    f_lead = (Lft - CH_SPACING - (top - start_y) + 2 * loss) / 2
    if f_lead < 0.2:
        raise ValueError(f'filter {name}: lead {f_lead:.3f} mm too short')
    P.serpentine_tail(design, 'filter_tail' + s, ('jf' + s, 'prime_start'), Lft, f_lead, side, -1,
                      height=top - (start_y - f_lead) + 1e-6)
    # resonator tail: up to feed - 0.65, then a 2-step serpentine toward the centre (runs RES_HEIGHT)
    r_lead = (P.TAP_Y - P.TEE_HALF) - y_bottom
    P.serpentine_tail(design, 'res_tail' + s, ('jr' + s, 'prime_end'), Lr - P.RES_TAP - P.TEE_HALF,
                      r_lead, -side, 1, height=RES_HEIGHT)
    for nm, tgt, fixed in [('filter_tail' + s, Lf, P.FILTER_TAP), ('res_tail' + s, Lr, P.RES_TAP)]:
        actual = d5.f3.path_length(design, nm) + fixed + P.TEE_HALF
        if abs(actual - tgt) > 0.001:
            raise ValueError(f'{nm}: drawn {actual:.5f} mm vs {tgt:.5f} mm')
    radius = design.parse_value(FILLET)
    overlap = (-P.TEE_HALF - f_lead + radius, P.TEE_HALF + r_lead - radius)
    return dict(n0='n0' + s, jf_x=jf_x, jr_x=jr_x, y_tap=y_tap, open_len=open_len,
                parallel_limits_from_tap_mm=overlap)


def build_design(jj_metal=True, lengths=None, xmon=None, readout_claw=None, coupler=None, lj=None, qsquid=None, qflux=None,
                 marks=None, test_junctions=False, channels=True, control_lines=True, landing_pads=True,
                 fab_contacts=False, wafer_ports=False, port_options=None, readout_hammer=None):
    """Full chip. jj_metal=False: mask version (Design3 convention, gds_export.md): the junction table is emptied on
    export; landing_pads=True keeps the SQUID landing electrodes of all three SQUIDs with the Design3 mask overlap /
    shift (SQUID_MASK); fab_contacts=True adds the Dolan contact notches (fab GDS only).
    xmon / readout_claw / coupler / lj: overrides of XMON / READOUT_CLAW / COUPLER / LJ (tuning).
    channels=False: no feed line / Purcell channels (core-only model; the claws keep a 40 um CPW stub)."""
    xm = dict(XMON, **(xmon or {}))
    rc = dict(READOUT_CLAW, **(readout_claw or {}))
    cp = dict(COUPLER, **(coupler or {}))
    lj = dict(LJ, **(lj or {}))
    lengths = lengths or lengths_from_params()
    if fab_contacts and jj_metal:
        raise ValueError('fab contact notches only on the junction-free mask (jj_metal=False)')
    sq = dict(QSQUID, **(qsquid or {}))
    squid_metal = jj_metal or landing_pads
    if not jj_metal and landing_pads:        # mask: Design3 landing electrodes
        sq.update(contact_overlap=SQUID_MASK['overlap'], shift_toward_arm=SQUID_MASK['shift'],
                  fab_contacts=str(bool(fab_contacts)))
        cp.update(squid_contact_overlap=SQUID_MASK['overlap'], squid_shift_toward_stem=SQUID_MASK['shift'])
    TCcls = ClawCouplerFab if fab_contacts else ClawCoupler

    design = d5.designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = d5.SI_MATERIAL
    design.chips.main.size.update(dict(size_x='10mm', size_y='10mm', size_z='-500um', center_x='0mm', center_y='0mm'))
    design.variables.update(cpw_width=CPW_RO['trace_width'], cpw_gap=CPW_RO['trace_gap'])

    # ---- qubits and coupler ----
    for name, sx, loc in (('Q1', -1, '0'), ('Q2', 1, '180')):       # readout claw on the outer arm
        fbl = dict(QFLUX, **(qflux or {}), mirror=str(sx > 0))
        Xmon(design, name, options=dict(**pos(sx * QX, YQ), **xm, Lj_a=lj[name + 'a'], Lj_b=lj[name + 'b'],
                                        squid=sq, flux_bias_line_options=fbl,
                                        jj_metal='True' if squid_metal else 'False',
                                        flux_as_ground=cp.get('flux_as_ground', 'False'),
                                        connection_pads=dict(readout=dict(rc, connector_location=loc))))
    arm = design.parse_value(xm['cross_length'])
    L_bar = 2 * (QX - arm) - 2 * design.parse_value(cp['claw_sep'])
    TCcls(design, 'TC', options=dict(**pos(0, YQ), pad_arm_length=f'{L_bar * 1e3:.4f}um',
                                     qubit_arm_width=xm['cross_width'], Lj=lj['TC'],
                                     jj_metal='True' if squid_metal else 'False', **cp))

    # ---- launch pads ----
    ports = wafer_pcb_ports if port_options is None else port_options
    used = ('P1', 'P2', 'P4', 'P5', 'P6') + (('P7', 'P8') if channels else ())      # P3 free
    for pn in used:
        cpw = CPW_RO if pn in ('P7', 'P8') else CPW_CTL
        add_bond_pad(design, 'LP' + pn, ports[pn], **cpw)
    if wafer_ports:
        # The pocket and centre conductor end together at the tie pin, so
        # the unused launch joins the ground there without a control route.
        add_bond_pad(design, 'LPP3', ports['P3'], **CPW_CTL)

    # ---- readout: feed line + two Purcell channels ----
    info = {}
    if channels:
        claw = design.components['Q1'].pins['readout'].middle
        # feed height so that the open run (vertical + horizontal, one corner) is exactly RES_TAP - TEE_HALF
        loss = design.parse_value(FILLET) * (2 - np.pi / 2)
        horiz = abs(claw[0] - (-JR_X))
        y_tap = claw[1] + (P.RES_TAP - P.TEE_HALF) - horiz + loss + P.TEE_HALF
        feed_y = y_tap + P.TAP_Y
        for name, ch in READOUT.items():
            pin = design.components[name].pins['readout'].middle
            info[name] = add_channel(design, name, ch['side'], lengths[name]['L_filter'], lengths[name]['L_res'],
                                     feed_y, pin, readout_hammer)
        fy, xs, yt = feed_y, PORT_FEED['x_side'], PORT_FEED['y_top']
        route(design, 'feed_in', ('LPP8', 'tie'), (info['Q1']['n0'], 'prime_start'),
              [(-1.1, yt), (-xs, yt), (-xs, fy)], CPW_RO)
        P.straight(design, 'feed_mid', (info['Q1']['n0'], 'prime_end'), (info['Q2']['n0'], 'prime_start'))
        route(design, 'feed_out', ('LPP7', 'tie'), (info['Q2']['n0'], 'prime_end'),     # drawn from the pad (mirror of feed_in)
              [(1.1, yt), (xs, yt), (xs, fy)], CPW_RO)
        design.metadata['feed_y_mm'] = feed_y

    # ---- control lines ----
    if control_lines:
        cw, cg = design.parse_value(xm['cross_width']), design.parse_value(xm['cross_gap'])
        for name, sx, xy_port, fl_port in (('Q1', -1, 'LPP1', 'LPP2'), ('Q2', 1, 'LPP6', 'LPP5')):
            xe, ye = sx * (QX + XY_END['dx']), YQ - cw / 2 - cg - XY_END['sep']
            OpenToGround(design, 'xy_end_' + name, options=dict(**pos(xe, ye), orientation='90', **{
                'width': CPW_CTL['trace_width'], 'gap': CPW_CTL['trace_gap']}))
            route(design, 'xy_' + name, (xy_port, 'tie'), ('xy_end_' + name, 'open'),
                  [(sx * 3.3, design.components[xy_port].pins['tie'].middle[1]),
                   (sx * 3.3, -0.6), (xe, -0.6)], CPW_CTL)
            xf, yf = design.components[name].pins['flux_bias_line'].middle
            route(design, 'fl_' + name, (name, 'flux_bias_line'), (fl_port, 'tie'),
                  [(xf, -1.0), (sx * 1.6, -1.0),
                   (sx * 1.6, design.components[fl_port].pins['tie'].middle[1])], CPW_CTL)
        tc = design.components['TC']
        xf, yf = tc.pins['flux_bias_line'].middle
        route(design, 'fl_TC', ('TC', 'flux_bias_line'), ('LPP4', 'tie'), [(xf, -3.9), (1.1, -3.9)], CPW_CTL)

    # ---- marks / test junctions (mask only) ----
    if (not jj_metal) if marks is None else marks:
        from xmark import add_corner_marks
        add_corner_marks(design, die_id=DIE_ID, inset=0.35, label_corner='LR', label_left='True', **MARK)
    if test_junctions:
        from test_junction import TestJunction
        t = TEST_JJ
        (wl, hl), (wr, hr) = t['pad_left'], t['pad_right']
        for nm, sx, sy in (('LL', -1, -1), ('LR', 1, -1), ('UL', -1, 1), ('UR', 1, 1)):
            TestJunction(design, f'testjj_{nm}', options=dict(
                chip='main', pos_x=f"{sx * t['x']}mm", pos_y=f"{sy * t['y']}mm",
                pocket_width=f"{t['pocket'][0]}mm", pocket_height=f"{t['pocket'][1]}mm",
                pad_l_width=f'{wl}mm', pad_l_height=f'{hl}mm', pad_r_width=f'{wr}mm', pad_r_height=f'{hr}mm',
                pad_gap=f"{t['pad_gap']}mm"))
    design.metadata['channels'] = info
    return design


def save_params(lengths, extra=None):
    os.makedirs(PAR_DIR, exist_ok=True)
    fn = _relocate(os.path.join(PAR_DIR, TAG + '.json'))
    old = json.load(open(fn)) if os.path.exists(fn) else {}
    par = dict(old, design='Design6', date=time.strftime('%Y-%m-%d'),
               description='3 flux-tunable transmons: 2 Xmon (asymmetric SQUID) + grounded SQUID coupler (CZ), Design4 Purcell readout on P8-P7',
               substrate=dict(material='Si', eps_r=float(d5.SI_EPS), thickness_um=500),
               chip=dict(size_mm=[10, 10], qubit_x_mm=QX, qubit_y_mm=YQ),
               cpw=dict(readout=CPW_RO, control=CPW_CTL), xmon=XMON, readout_claw=READOUT_CLAW,
               coupler={k: (dict(v) if isinstance(v, dict) else v) for k, v in COUPLER.items()}, Lj=LJ,
               qubit_squid=QSQUID, qubit_flux_line=QFLUX, xy_end_mm=XY_END,
               C_p=dict(type='T-type enclosing', **P.ENC_CAI), C_j_hammer=HAMMER,
               routing=dict(cpw_min_spacing_um=400, filter_tap_mm=P.FILTER_TAP, res_tap_mm=P.RES_TAP,
                            jr_x_mm=JR_X, fillet=FILLET),
               channels=[dict(name=k, L_filter_mm=round(v['L_filter'], 5),
                              L_res_mm=round(v['L_res'], 5),
                              **{kk: vv for kk, vv in v.items() if kk not in ('L_filter', 'L_res')})
                         for k, v in lengths.items()],
               solver_rules=dict(q3d=dict(passes=8), eigenmode=dict(passes=6, mesh_um=10), s21=dict(passes=3, mesh_um=10)))
    par.update(extra or {})
    json.dump(par, open(fn, 'w'), indent=1, default=str)
    return fn


def check(design):
    """CPW spacing (>= 400 um rule, Design4 checker) and pins."""
    try:
        import cpw_spacing_check as S
        S.report(design, 400)
    except Exception as e:
        print('spacing check skipped:', repr(e)[:200])
    for n in ('Q1', 'Q2', 'TC'):
        for pn, pin in design.components[n].pins.items():
            print(f'  {n}.{pn}: middle {np.round(pin.middle, 4)}, normal {np.round(pin.normal, 3)}')


def export_simulation_gds(design, filename):
    """Export solver electrodes and the layer 50/0 lumped-junction annotations."""
    import gdstk
    from gds_tools import export_gds

    junctions = list(design.qgeometry.tables['junction'].itertuples())
    export_gds(design, str(filename))
    library = gdstk.read_gds(str(filename))
    top, = library.top_level()
    for junction in junctions:
        marker, = gdstk.FlexPath(np.asarray(junction.geometry.coords) * 1000, junction.width * 1000,
                                 layer=50).to_polygons()
        top.add(marker)
    library.write_gds(str(filename))
    return junctions


def export_masks(f_sim, f_fab, *, test_junctions, wafer_ports=False, port_options=None,
                 readout_hammer=None, lengths=None):
    """Export and verify both variants, without changing the saved model inputs."""
    from pathlib import Path
    from tempfile import TemporaryDirectory
    from gds_tools import export_gds
    round_mask_dir = _I('Design3/layout/components')
    sys.path.insert(0, round_mask_dir)
    try:
        from round_mask import round_mask
    finally:
        sys.path.remove(round_mask_dir)
    from verify_masks import check_masks

    # The simulator consumes the original junction table. In the GDS, preserve its
    # locations and widths as separate annotation polygons on layer 50/0.
    settings = dict(wafer_ports=wafer_ports, port_options=port_options,
                    readout_hammer=readout_hammer, lengths=lengths)
    d = build_design(**settings)
    fab = build_design(jj_metal=False, test_junctions=test_junctions,
                       fab_contacts=True, **settings)
    unnotched = build_design(jj_metal=False, **settings)
    with TemporaryDirectory(dir=f_sim.parent) as temp:
        staging = Path(temp)
        sim = staging / f_sim.name
        mask = staging / f_fab.name
        junctions = export_simulation_gds(d, sim)
        raw_mask = staging / 'raw_fab.gds'
        export_gds(fab, str(raw_mask), fab=True)
        round_mask(raw_mask, mask, 5.0)
        check_masks(sim, mask, junctions, unnotched, fab)
        sim.replace(f_sim)
        mask.replace(f_fab)


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')    # round_mask prints 'µm'
    import geopandas as gpd
    if not hasattr(gpd.GeoDataFrame, 'append'):
        gpd.GeoDataFrame.append = gpd.GeoDataFrame._append
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--date', default=time.strftime('%Y-%m-%d'))
    parser.add_argument('--retune', action='store_true')
    parser.add_argument('--no-gui', action='store_true')
    parser.add_argument('--test-junctions', choices=('yes', 'no'), nargs='?', const='yes')
    args = parser.parse_args()
    TAG = 'Design6_' + args.date
    if args.retune or _length_file() is None:
        save_params(tune_lengths())
    d = build_design()
    print('feed y =', round(d.metadata['feed_y_mm'], 4), 'mm; channels:', d.metadata['channels'])
    check(d)
    os.makedirs(GDS_DIR, exist_ok=True); os.makedirs(FIG_DIR, exist_ok=True)
    tj = args.test_junctions
    while tj not in ('yes', 'no'):
        tj = input('Include corner test-junction pads in the fab mask? [yes/no] ').strip().lower()
    tj = tj == 'yes'
    f_sim = Path(GDS_DIR) / (TAG + '_sim.gds')
    f_fab = Path(GDS_DIR) / (TAG + '_fab.gds')
    export_masks(f_sim, f_fab, test_junctions=tj)
    print('wrote', f_sim, f_fab)
    if not args.no_gui:
        from qiskit_metal import MetalGUI
        gui = MetalGUI(build_design(jj_metal=False, test_junctions=tj, fab_contacts=True))
        gui.rebuild(); gui.autoscale()
        gui.screenshot(name=_relocate(os.path.join(FIG_DIR, TAG + '_layout')), display=False)
        gui.zoom_on_components(['Q1', 'TC', 'Q2'])
        gui.screenshot(name=_relocate(os.path.join(FIG_DIR, TAG + '_core')), display=False)
        gui.zoom_on_components(['cp_Q1', 'cj_Q1', 'jf_Q1', 'jr_Q1'])
        gui.screenshot(name=_relocate(os.path.join(FIG_DIR, TAG + '_channel_q1')), display=False)

        gui.main_window.close()
