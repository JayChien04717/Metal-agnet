# -*- coding: utf-8 -*-
r"""Single layout entry point for Design12_Transmon2Q_TransmissionPurcell (Stage 1 simulation geometry).

Same transmon core as Design11_Transmon2Q_Manhattan r9 (two flux-tunable Manhattan Xmons Q1 / Q2, grounded
tunable coupler TC, XY / flux lines, 10 x 10 mm 8-port die), but the two individual PurcellReadout channels are
replaced by ONE SHARED TRANSMISSION-TYPE lambda/2 Purcell filter:

    P8 -> IDC_in -> PF (lambda/2, open-open, capacitively coupled at both ends) -> IDC_out -> P7

Both lambda/4 readout resonators (RR_Q1, RR_Q2) couple to the central horizontal bus of PF through a short
side-by-side coupling section next to their shorted end (Design8_reflection_purcell_resonator style, 6 um ground
strip between the two CPW gaps); their open ends drop straight into the Q1 / Q2 north readout claws.

Sketch (not to scale, y up)::

            LP P8                                   LP P7
              |                                       |
           [IDC_in]                               [IDC_out]
              |                                       |
              |       ____________PF bus___________   |
              |      |  ==coup==        ==coup==   |  |
              |______|  |RR_Q1 |        | RR_Q2|   |__|      (filter dips to y_lo on both sides
                        |meandr|        |meandr|              to set its lambda/2 length)
                        |      |        |      |
                        Q1 ----- TC ----- Q2
                     (XY / flux lines as Design11)

All readout lengths are Stage 1 draft estimates from the CPW line model (code/resonator/purcell_readout_model.py)
with Design11 Q3D claw end loads and Design8 Q3D IDC loads; nothing here is EM-validated for this geometry.
"""
from collections import OrderedDict
from pathlib import Path
import math
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
from qiskit_metal.qlibrary.core.qroute import QRoute
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight

from bondpad import add_bond_pad
from check_keepouts import check_keepouts
from claw_coupler_manhattan import ClawCouplerManhattan
from fluxline import B2_FLUX_OPTIONS, T_FLUX_OPTIONS
from interdigital_cap import InterdigitalCap, alley_cap_fF
from junction_test_coupon import ManhattanTestJunction
from logo import SquidFlowerLogo
from tap_junction import TapJunction
from purcell_readout_model import CPW as CPW_MODEL, open_end_c
from round_corners import round_component_polys, square_merged_etch_corners
from wafer_pcbport import wafer_pcb_ports
from xmon_manhattan import XmonManhattan

REVISION = 'r7'
CPW_RO = dict(trace_width='10um', trace_gap='6um')
CPW_CTL = dict(trace_width='10um', trace_gap='6um')
FILLET = '99um'

# ---- Transmon core: Design11 r9 core + Design12 r5 flux-line ('T') & compact SQUID (36x30 um) update ----------
QX = 0.40       # mm, Q1 at -QX, Q2 at +QX
YQ = 0.0        # mm
R5_FLUX_OPTIONS = dict(T_FLUX_OPTIONS)
R5_SQUID_OPTIONS = dict(jj_gap='36um', squid_spacing='36um', squid_bar_height='6um')
XMON_OPT = dict(cross_length='142um', cross_width='28um', cross_gap='14um')
RO_CLAW = dict(
    claw_length='122um', claw_width='10um', claw_sep='18um', claw_gap='6um',
    cpw_width='10um', cpw_gap='6um', cpw_length='40um',
)
XMON_BY_QUBIT = dict(Q1=dict(cross_length='150.4um'), Q2=dict(cross_length='150.5um'))
RO_CLAW_BY_QUBIT = dict(Q1=dict(claw_length='100um'), Q2=dict(claw_length='114um'))
COUPLER_OPT = dict(
    pad_arm_width='15um', paddle_width='30um', paddle_height='78um', paddle_gap='55um',
    pad_head_width='15um', pad_head_length='285um', pad_gap='30um',
    jj_gap=R5_SQUID_OPTIONS['jj_gap'],
    claw_width='10um', claw_length='17um', claw_sep='12.5um', claw_gap='12um',
    squid_spacing=R5_SQUID_OPTIONS['squid_spacing'],
    squid_bar_height=R5_SQUID_OPTIONS['squid_bar_height'],
    squid_contact_overlap='3um',
    flux_bias_line_options=Dict(R5_FLUX_OPTIONS, make_fbl=True),
    charge_line_options=Dict(make_cl=False),
    readout_line_options=Dict(make_rol=False),
)
LJ = dict(Q1a='17.39nH', Q1b='34.79nH', Q2a='20.55nH', Q2b='41.11nH', TC='14.03nH')
XY_SEP_MM = 0.14

# ---- Shared transmission Purcell filter (Stage 1 draft targets) -------------------------------------------------
# Readout targets kept from Design11 option B: f_r 6.35 / 6.15 GHz, kappa_r = |2chi| = 3 MHz.
# Shared filter: lambda/4 (open at the weak input IDC, shorted at the far end), f_PF = 6.25 GHz, kappa_PF ~ 300 MHz
# (Q_F ~ 21), strongly asymmetric so the output collects kappa_out / kappa_PF ~ 95 % of the emitted readout signal:
#   input  P7 -> weak IDC_in at the open end, kappa_in ~ 15 MHz
#   output P8 <- galvanic tap (TapJunction) at tap_mm from the short, kappa_out = (4 w / pi) sin^2(beta l) ~ 285 MHz
# Needed resonator-filter coupling: J = sqrt(kappa_r kappa_PF (1 + (2 delta / kappa_PF)^2)) / 2 ~ 17.4 MHz.
LINE = CPW_MODEL['10/6']            # eps_r 11.45 line model (v, Z0, C'), matches chip material silicon_eps11p45
# r3 (user 2026-11-26): filter DESIGN (model) value f_PF = 6.5 GHz, no HFSS anchor on the filter. With the r1
# model / HFSS ratio (5.9987 / 6.25) it is expected near 6.24 GHz in HFSS, between the two resonators.
# History (filter model f): r1 6.25; r2 6.772 (= 6.5 GHz HFSS target via the r1 anchor; misread request, superseded).
PF = dict(f_GHz=6.5, kappa_in_MHz=15.0, kappa_out_MHz=285.0,
          x_in_mm=2.0, x_tap_mm=-1.55, y_bus_mm=2.60, tap_prime_um=300.0, tap_stub_um=150.0, min_drop_mm=0.30)
# HFSS anchor (r1 3-pass driven S21, re-analysis results/resonator/20261126_165813_fit_s21_r1_3p_d6c72cc5, parent
# 20261126_163208_purcell_s21_driven_r1_3p_2b2c127c): model target f / HFSS fitted f of the same geometry. Each line
# is drawn for f_model = f_target x f_model,r1 / f_HFSS,r1. 3-pass driven f is not converged evidence
# (skill/em_simulation.md section 4); the r1 values are the only reference so far.
# r4: resonator anchors moved to the r3 3-pass S21 (run 20261126_180157_purcell_s21_driven_r3_3p_20b82629; f_model =
# the r3 drawn model target, f_hfss = r3 3-mode fit). The r1 Q2 anchor did not carry over because most of its -152 MHz
# offset came from the 460 um coupling section (r3: 244 um, offset -30 MHz). PF entry kept for reference only (r3+
# draws the filter at its design value). History: r1 anchor Q1 6.35 -> 6.31642, Q2 6.15 -> 5.99805 (used in r2, r3).
HFSS_ANCHOR = dict(run='20261126_180157_purcell_s21_driven_r3_3p_20b82629',
                   PF=dict(f_model=6.25, f_hfss=5.99870, run='20261126_165813_fit_s21_r1_3p_d6c72cc5'),
                   Q1=dict(f_model=6.38384, f_hfss=6.32485),
                   Q2=dict(f_model=6.30581, f_hfss=6.27625))


def model_target(key, f_target):
    a = HFSS_ANCHOR[key]
    return f_target * a['f_model'] / a['f_hfss']
# IDC_in: Design8 hexagon-tapered comb style (6 / 6 um fingers, 20 um bus, 100 um taper, 24 um pocket) with fewer
# fingers. Loads scaled from the Design8 Q3D run 20261120_231033_q3d_purcell_idc_3b7c138a (N = 10, overlap 200 um:
# C_in 98.67 fF = 0.938 x Alley, comb ground 54.84 fF de-embedded): C_in via Alley x 0.938, ground C x comb height.
IDC_OPT = dict(
    finger_count='6', finger_width='6um', finger_gap='6um', overlap='75um', end_gap='6um',
    bus_width='20um', taper_length='100um', lead_length='45um', ground_gap='24um',
)
IDC_Q3D_REF = dict(n=10, overlap_um=200.0, c_in_fF=98.67, c_gnd_fF=54.84)
# Readout resonators: lambda/4, short at the coupling-section tip, open end into the qubit north claw.
# End load = Design11 r9 Kron claw load (41.75 / 45.14 fF). Coupling length: Design8 (380 um -> J ~ 31 MHz near the
# lambda/2 bus current antinode, i.e. ~0.082 MHz/um) scaled by cos(beta z) of the filter current at the section
# (z = distance from the filter short): RR_Q1 sits near the short, RR_Q2 near the open end. Rough draft only;
# Stage 2 Q3D / HFSS sets the final value.
RR = dict(strip_um=6.0, pitch_mm=0.40, bars=3)
# r4 coupling lengths: r3 HFSS J = 22.5 / 20.4 MHz at 122 / 244 um, f_PF 6.233 GHz, kappa_PF 357 MHz. With the
# resonators back on target (delta = +117 / -83 MHz) kappa_r,eff would be 3.96 / 3.84 MHz; scale J by sqrt(3 / that)
# -> 106 / 216 um.
# r3 coupling lengths: r1 HFSS fit J = 35.9 / 33.8 MHz at 230 / 460 um (Design8 scaling underestimated J ~2x).
# Needed J for kappa_r,eff = 3 MHz with the filter expected at 6.5 x 5.9987 / 6.25 = 6.239 GHz in HFSS and
# kappa_PF = 333 MHz (r1 fit): J^2 = kappa_r ((kappa_PF/2)^2 + delta^2) / kappa_PF, delta = +111 / -89 MHz
# -> J = 19.0 / 17.9 MHz; lengths scaled linearly in J. Resonator lengths use the r1 HFSS anchor (model_target).
# History (coupling um): r1 230 / 460; r2 136 / 501 (superseded); r3 122 / 244.
RR_CHANNELS = {
    'Q1': dict(side=-1, f_GHz=6.35, c_claw_fF=41.75, coupling_um=106.0),
    'Q2': dict(side=+1, f_GHz=6.15, c_claw_fF=45.14, coupling_um=216.0),
}
CORNER_LOSS = lambda r_mm: r_mm * (2.0 - math.pi / 2.0)     # filleted length lost per 90-degree bend


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


def _min_segment_check(name, pts, r_mm):
    """DRC B2: every inner straight segment >= 2 fillet, end segments >= fillet."""
    pts = [np.asarray(p, float) for p in pts]
    for i in range(len(pts) - 1):
        seg = float(np.linalg.norm(pts[i + 1] - pts[i]))
        need = r_mm if i in (0, len(pts) - 2) else 2.0 * r_mm
        if seg < need - 1e-9:
            raise ValueError(f'{name}: segment {i} is {seg * 1e3:.1f} um < {need * 1e3:.1f} um')


def route_filleted_mm(design, comp_name, r_mm):
    """Filleted centre-line length of a route as drawn (poly length minus the per-corner fillet saving)."""
    t = design.qgeometry.tables['path']
    row = t[(t['component'] == design.components[comp_name].id) & (~t['subtract'])].iloc[0]
    c = np.asarray(row.geometry.coords, float)
    n_corner = 0
    for p0, p1, p2 in zip(c[:-2], c[1:-1], c[2:]):
        a, b = p1 - p0, p2 - p1
        if np.linalg.norm(a) > 1e-12 and np.linalg.norm(b) > 1e-12 and abs(a @ b) < 1e-9:
            n_corner += 1
    return float(row.geometry.length) - n_corner * CORNER_LOSS(r_mm)


# ---- draft circuit estimates (line model; not EM results) -------------------------------------------------------
def idc_loads(n_fingers, overlap_um):
    """Series C_in and comb ground C [fF] of the IDC, scaled from the Design8 Q3D coupon (see IDC_Q3D_REF)."""
    ref = IDC_Q3D_REF
    q3d_over_alley = ref['c_in_fF'] / alley_cap_fF(ref['n'], ref['overlap_um'])
    c_in = alley_cap_fF(n_fingers, overlap_um) * q3d_over_alley
    c_gnd = ref['c_gnd_fF'] * (2 * n_fingers - 1) / (2 * ref['n'] - 1)    # comb height N w + (N-1) g, w = g
    return c_in, c_gnd


def port_load(f_ghz, n_fingers, overlap_um, z_port=50.0):
    """Admittance seen from the filter open end into IDC -> 50 Ohm port (comb grounds on both sides)."""
    c_in, c_gnd = idc_loads(n_fingers, overlap_um)
    w = 2 * math.pi * f_ghz * 1e9
    y_port = 1.0 / z_port + 1j * w * c_gnd * 1e-15
    y_ser = 1j * w * c_in * 1e-15
    return y_ser * y_port / (y_ser + y_port) + 1j * w * c_gnd * 1e-15


def filter_length_mm(f_ghz, n_fingers, overlap_um):
    """lambda/4 length (open end loaded by the input IDC, galvanic tap assumed reactance-free):
    beta L = pi / 2 - atan(B / Y0)."""
    y = port_load(f_ghz, n_fingers, overlap_um)
    theta = math.atan(y.imag * LINE['z0'])
    return LINE['v'] / (4 * f_ghz * 1e9) * (1 - 2 * theta / math.pi) * 1e3


def filter_c_eff_f(length_mm):
    return LINE['c_line'] * 1e-15 * length_mm / 2.0                 # lambda/4 mode, c_line in fF/mm


def kappa_in_mhz(f_ghz, n_fingers, overlap_um, length_mm):
    """kappa_in = Re Y_load / C_eff at the open end."""
    return float(port_load(f_ghz, n_fingers, overlap_um).real) / filter_c_eff_f(length_mm) / (2 * math.pi) / 1e6


def kappa_tap_mhz(f_ghz, tap_mm, length_mm):
    """Galvanic tap at tap_mm from the short into a matched line: V_tap = V0 sin(beta l), P = V_tap^2 / 2 Z0,
    kappa = sin^2(beta l) / (Z0 C_eff)."""
    beta = 2 * math.pi * f_ghz * 1e9 / LINE['v'] * 1e-3              # rad / mm
    return math.sin(beta * tap_mm) ** 2 / (LINE['z0'] * filter_c_eff_f(length_mm)) / (2 * math.pi) / 1e6


def tap_for_kappa_mm(f_ghz, kappa_mhz, length_mm):
    beta = 2 * math.pi * f_ghz * 1e9 / LINE['v'] * 1e-3
    s2 = kappa_mhz * 1e6 * 2 * math.pi * LINE['z0'] * filter_c_eff_f(length_mm)
    return math.asin(math.sqrt(s2)) / beta


def resonator_length_mm(f_ghz, c_end_ff):
    """Open-short lambda/4 with the claw end load (L_eff = L + C_end / C')."""
    return LINE['v'] / (4 * f_ghz * 1e9) * 1e3 - c_end_ff / LINE['c_line']


def build_design(variant='sim', overrides=None, channels=True, flux_as_ground=False, check_keepout_um=100.0,
                 readout_only=False):
    """Build Design12_Transmon2Q_TransmissionPurcell (Stage 1 simulation geometry).

    readout_only=True: EM fixture variant for the filter + resonators without qubits. Each resonator then ends in an
    etched open end (OpenToGround 'RR_<Q>_open') below the claw pin instead of at the qubit claw, extended by
    (C_claw - C_open) / C' so that the bare lambda/4 keeps its target f_r; the meander is identical to the chip
    layout. Render only the readout components (qubits are still built but must not be rendered)."""
    if variant not in ('sim', 'fab'):
        raise ValueError(f'variant must be "sim" or "fab"; got {variant!r}')
    if variant == 'fab' and (readout_only or flux_as_ground):
        raise ValueError('fab variant: readout_only / flux_as_ground are simulation-only options')
    # fab (Stage 3, NCKU single die, user 2026-11-26): only jj_metal='False' differs from sim. The Manhattan Xmons and
    # coupler have no Dolan contact hooks, and the junction sites stay empty (NCKU draws the junctions).
    jjm = 'True' if variant == 'sim' else 'False'
    ov = overrides or {}
    xmon_opt = dict(XMON_OPT, **ov.get('xmon', {}))
    coupler_opt = dict(COUPLER_OPT, **ov.get('coupler', {}))
    pf = dict(PF, **ov.get('pf', {}))
    idc_opt = dict(IDC_OPT, **ov.get('idc', {}))
    rr = dict(RR, **ov.get('rr', {}))
    qx = float(ov.get('qx_mm', QX))
    fag = 'True' if flux_as_ground else 'False'
    r_mm = float(FILLET.rstrip('um')) / 1000.0

    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon_eps11p45'
    design.chips.main.size.update(dict(size_x='10mm', size_y='10mm', size_z='-500um', center_x='0mm', center_y='0mm'))
    design.variables.update(cpw_width=CPW_RO['trace_width'], cpw_gap=CPW_RO['trace_gap'])

    # ---- Qubits + coupler (r5: 'T' flux lines + compact 36x30 um SQUID loops) ----
    xmon_q, arms = {}, {}
    for name, sx in (('Q1', -1), ('Q2', +1)):
        fbl = Dict(R5_FLUX_OPTIONS, **ov.get('fbl', {}), **ov.get('fbl_' + name, {}),
                   make_fbl=True, mirror=str(sx < 0))
        q_ro = {**RO_CLAW, **RO_CLAW_BY_QUBIT[name], **ov.get('ro_claw', {}), **ov.get('ro_claw_' + name, {})}
        xmon_q[name] = {**xmon_opt, **XMON_BY_QUBIT[name], **ov.get('xmon_' + name, {})}
        sq_ov = dict(XmonManhattan.default_options.squid, **R5_SQUID_OPTIONS,
                     **ov.get('xmon', {}).get('squid', {}), **ov.get('squid', {}), **ov.get('squid_' + name, {}))
        xmon_q[name]['squid'] = Dict(**sq_ov)
        arms[name] = design.parse_value(xmon_q[name]['cross_length'])
        XmonManhattan(design, name, options=dict(
            **_pos(sx * qx, YQ), **xmon_q[name],
            Lj_a=LJ[name + 'a'], Lj_b=LJ[name + 'b'],
            jj_metal=jjm, flux_as_ground=fag,
            flux_bias_line_options=fbl,
            connection_pads=dict(readout=dict(q_ro, connector_location='90')),
        ))
        raw_etcher = _poly_geom(design, name, 'readout_connector_etcher')
        round_component_polys(design.components[name], 5e-3)
        square_merged_etch_corners(design.components[name], 'readout_connector_etcher', raw_etcher,
                                   _poly_geom(design, name, 'cross_etch'), pad=5e-3)

    claw_sep = design.parse_value(coupler_opt['claw_sep'])
    L_bar = (2.0 * qx - arms['Q1'] - arms['Q2']) - 2.0 * claw_sep
    x_tc = 0.5 * (arms['Q1'] - arms['Q2'])
    ClawCouplerManhattan(design, 'TC', options=dict(
        **_pos(x_tc, YQ), pad_arm_length=f'{L_bar * 1e3:.4f}um',
        qubit_arm_width=xmon_opt['cross_width'], Lj=LJ['TC'],
        jj_metal=jjm, flux_as_ground=fag, **coupler_opt,
    ))
    round_component_polys(design.components['TC'], 5e-3)

    for pn in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'):
        if not channels and pn in ('P7', 'P8'):
            continue
        add_bond_pad(design, 'LP' + pn, wafer_pcb_ports[pn], **(CPW_RO if pn in ('P7', 'P8') else CPW_CTL))

    report = {}
    if channels:
        # ---- lambda/4 shared filter: P7 -> IDC_in (open end) -> bus -> output tap T -> short; T stub -> P8 ----
        n_f = int(idc_opt['finger_count'])
        overlap_um = design.parse_value(idc_opt['overlap']) * 1e3
        f_pf = pf['f_GHz']                                                     # design value, no anchor (r3)
        L_pf = filter_length_mm(f_pf, n_f, overlap_um)
        tap = tap_for_kappa_mm(f_pf, pf['kappa_out_MHz'], L_pf)               # tap centre -> short, mm
        xin, xt, yb = pf['x_in_mm'], pf['x_tap_mm'], pf['y_bus_mm']
        half_prime = pf['tap_prime_um'] / 2000.0
        TapJunction(design, 'PF_tap', options=dict(
            **_pos(xt, yb), orientation='0', prime_length=f"{pf['tap_prime_um']:g}um",
            stub_length=f"{pf['tap_stub_um']:g}um", stub_up='True', **CPW_RO))
        # short arm: T prime_start -> short, straight along -x
        ShortToGround(design, 'PF_short', options=dict(
            **_pos(xt - tap, yb), orientation='180', width=CPW_RO['trace_width'], gap=CPW_RO['trace_gap']))
        RouteStraight(design, 'PF_short_arm', Dict(
            pin_inputs=Dict(start_pin=Dict(component='PF_tap', pin='prime_start'),
                            end_pin=Dict(component='PF_short', pin='short')), **CPW_RO))
        if tap - half_prime < r_mm:
            raise ValueError(f'PF: short arm {1e3 * (tap - half_prime):.0f} um shorter than the fillet')
        # open side: IDC_in.a -> straight down at x_in -> bus at y_bus -> T prime_end
        L_open_route = L_pf - tap - half_prime
        tie7 = np.asarray(design.components['LPP7'].pins['tie'].middle, float)
        dy = L_open_route - (xin - (xt + half_prime)) + CORNER_LOSS(r_mm)
        if dy < pf['min_drop_mm']:        # short filter: move IDC_in inward, feed_in jogs from the P7 tie
            xin = xt + half_prime + L_open_route + CORNER_LOSS(r_mm) - pf['min_drop_mm']
            xin = min(xin, tie7[0] - 2 * r_mm - 1e-3)
            dy = L_open_route - (xin - (xt + half_prime)) + CORNER_LOSS(r_mm)
        y_a = yb + dy
        tmp = InterdigitalCap(design, 'IDC_in', options=dict(**_pos(0, 0), orientation='90', **idc_opt, **CPW_RO))
        half = tmp.center_to_pin_mm()
        if y_a + 2 * half > tie7[1] - r_mm:
            raise ValueError(f'IDC_in top {y_a + 2 * half:.3f} mm too close to the P7 tie {tie7[1]:.3f} mm')
        InterdigitalCap(design, 'IDC_in', options=dict(
            **_pos(xin, y_a + half), orientation='90', **idc_opt, **CPW_RO))
        if abs(xin - tie7[0]) < 1e-9:
            RouteStraight(design, 'feed_in', Dict(
                pin_inputs=Dict(start_pin=Dict(component='LPP7', pin='tie'), end_pin=Dict(component='IDC_in', pin='b')),
                **CPW_RO))
        else:
            b_in = np.asarray(design.components['IDC_in'].pins['b'].middle, float)
            y_j = 0.5 * (tie7[1] + b_in[1])
            in_pts = [tuple(tie7), (tie7[0], y_j), (xin, y_j), tuple(b_in)]
            _min_segment_check('feed_in', in_pts, r_mm)
            _route(design, 'feed_in', ('LPP7', 'tie'), ('IDC_in', 'b'), in_pts[1:-1], CPW_RO)
        pf_pts = [(xin, y_a), (xin, yb), (xt + half_prime, yb)]
        _min_segment_check('PF', pf_pts, r_mm)
        _route(design, 'PF', ('IDC_in', 'a'), ('PF_tap', 'prime_end'), pf_pts[1:-1], CPW_RO)
        L_pf_drawn = (route_filleted_mm(design, 'PF', r_mm) + 2 * half_prime
                      + route_filleted_mm(design, 'PF_short_arm', r_mm))
        # output: T stub (second_end) -> up -> over to P8 -> up into the pad
        se = np.asarray(design.components['PF_tap'].pins['second_end'].middle, float)
        tie8 = np.asarray(design.components['LPP8'].pins['tie'].middle, float)
        y_mid = 0.5 * (se[1] + tie8[1])
        out_pts = [tuple(se), (se[0], y_mid), (tie8[0], y_mid), tuple(tie8)]
        _min_segment_check('feed_out', out_pts, r_mm)
        _route(design, 'feed_out', ('PF_tap', 'second_end'), ('LPP8', 'tie'), out_pts[1:-1], CPW_RO)
        report['PF'] = dict(
            f_target_GHz=pf['f_GHz'], f_model_GHz=f_pf, x_in_mm=xin, target_mm=L_pf, drawn_mm=L_pf_drawn, tap_from_short_mm=tap,
            idc_fingers=n_f, idc_overlap_um=overlap_um, c_in_fF=idc_loads(n_f, overlap_um)[0],
            c_comb_gnd_fF=idc_loads(n_f, overlap_um)[1],
            kappa_in_est_MHz=kappa_in_mhz(f_pf, n_f, overlap_um, L_pf),
            kappa_out_est_MHz=kappa_tap_mhz(f_pf, tap, L_pf))
        x_short_pf = xt - tap

        # ---- lambda/4 readout resonators side-coupled under the bus ----
        w_um = design.parse_value(CPW_RO['trace_width']) * 1e3
        g_um = design.parse_value(CPW_RO['trace_gap']) * 1e3
        y_c = yb - (w_um + 2 * g_um + rr['strip_um']) / 1000.0     # 6 um ground strip between the two gaps
        p, nb = rr['pitch_mm'], int(rr['bars'])
        if nb % 2 == 0:
            raise ValueError('RR bars must be odd so the last bar ends above the claw')
        for name, ch in RR_CHANNELS.items():
            s = ch['side']
            lc = ch['coupling_um'] / 1000.0 + r_mm      # straight parallel run = coupling_um; the bend fillet follows
            claw = np.asarray(design.components[name].pins['readout'].middle, float)
            xq = claw[0]
            f_rm = model_target(name, ch['f_GHz'])                         # HFSS-anchored model target
            L_r = resonator_length_mm(f_rm, ch['c_claw_fF'])
            n_corner = 1 + 2 * nb
            # L_poly = lc + nb p + nb W + (y_c - nb p - y_claw)  ->  W
            W = (L_r + n_corner * CORNER_LOSS(r_mm) - lc - (y_c - claw[1])) / nb
            if W < lc:                                   # short tip must stay above the meander block
                raise ValueError(f'RR_{name}: meander width {W * 1e3:.0f} um too short for the coupling section')
            x_out = xq + s * W                               # outer edge of the meander (away from the centre)
            x_short = x_out - s * lc                         # coupling section runs outward from the short
            pts = [(x_short, y_c), (x_out, y_c)]
            y = y_c
            for k in range(nb):
                y -= p
                pts += [(pts[-1][0], y), (xq if k % 2 == 0 else x_out, y)]
            ext = 0.0
            if readout_only:   # no claw load: extend the drop so the bare open end keeps f_r (layout meander unchanged)
                ext = (ch['c_claw_fF'] - open_end_c('10/6')) / LINE['c_line']
                OpenToGround(design, 'RR_' + name + '_open', options=dict(
                    **_pos(xq, claw[1] - ext), orientation='270',
                    width=CPW_RO['trace_width'], gap=CPW_RO['trace_gap'], termination_gap=CPW_RO['trace_gap']))
            pts.append((xq, claw[1] - ext))
            _min_segment_check('RR_' + name, pts, r_mm)
            if y - claw[1] < 2 * r_mm:
                raise ValueError(f'RR_{name}: drop to claw {1e3 * (y - claw[1]):.0f} um too short')
            ShortToGround(design, 'RR_' + name + '_short', options=dict(
                **_pos(x_short, y_c), orientation='0' if s < 0 else '180',   # pin normal points along the coupling run
                width=CPW_RO['trace_width'], gap=CPW_RO['trace_gap']))
            end = ('RR_' + name + '_open', 'open') if readout_only else (name, 'readout')
            _route(design, 'RR_' + name, ('RR_' + name + '_short', 'short'), end, pts[1:-1], CPW_RO)
            report['RR_' + name] = dict(f_target_GHz=ch['f_GHz'], f_model_GHz=f_rm, target_mm=L_r + ext, open_extension_mm=ext,
                                        drawn_mm=route_filleted_mm(design, 'RR_' + name, r_mm),
                                        meander_width_um=W * 1e3, coupling_um=ch['coupling_um'],
                                        coupling_from_pf_short_mm=sorted(abs(x - x_short_pf) for x in (x_short, x_out)),
                                        coupling_x_mm=sorted((x_short, x_out)))

    # ---- control lines (Design11 r9, unchanged) ----
    for name, sx, xy_port, fl_port in (('Q1', -1, 'LPP1', 'LPP2'), ('Q2', +1, 'LPP6', 'LPP5')):
        cg = design.parse_value(xmon_q[name]['cross_gap'])
        xe = sx * (qx + arms[name] + cg + XY_SEP_MM)
        OpenToGround(design, 'xy_end_' + name, options=dict(
            **_pos(xe, YQ), orientation='0' if sx < 0 else '180',
            width=CPW_CTL['trace_width'], gap=CPW_CTL['trace_gap'],
        ))
        xy_y = float(design.components[xy_port].pins['tie'].middle[1])
        _route(design, 'xy_' + name, (xy_port, 'tie'), ('xy_end_' + name, 'open'),
               [(sx * 3.45, xy_y), (sx * 3.45, YQ)], CPW_CTL)
        xf, yf = [float(v) for v in design.components[name].pins['flux_bias_line'].middle]
        fl_y = float(design.components[fl_port].pins['tie'].middle[1])
        y_bend = yf - 0.20
        dx_diag = abs(y_bend - fl_y)
        _polyroute(design, 'fl_' + name, (name, 'flux_bias_line'), (fl_port, 'tie'),
                   [(xf, y_bend), (xf + sx * dx_diag, fl_y)], CPW_CTL)

    xf_tc = float(design.components['TC'].pins['flux_bias_line'].middle[0])
    xp4 = float(design.components['LPP4'].pins['tie'].middle[0])
    _route(design, 'fl_TC', ('TC', 'flux_bias_line'), ('LPP4', 'tie'), [(xf_tc, -3.90), (xp4, -3.90)], CPW_CTL)

    # Two 4-probe Manhattan test-junction coupons in open lower ground (>= 600 um from routes, >= 299 um from markers)
    if channels:
        for tname, tx in (('tj_m1', -1.50), ('tj_m2', +1.50)):
            ManhattanTestJunction(design, tname, options=dict(
                **_pos(tx, -3.10), label='M', jj_metal='True',
            ))
        # Decorative 6-petal SQUID-flower logo in upper-center open ground
        SquidFlowerLogo(design, 'chip_logo', options=dict(
            **_pos(0.0, 3.90), size='1400um', stroke_width='35um', sparkle_stroke='18um', tone='etched',
        ))

    design.metadata['readout_report'] = report
    if check_keepout_um is not None:
        check_keepouts(design, min_um=float(check_keepout_um), include_site_markers=('NTHU', 'TSRI'))
    return design


def print_report(design):
    rep = design.metadata['readout_report']
    for k, v in rep.items():
        err = (v['drawn_mm'] - v['target_mm']) * 1e3
        extra = (f"tap {1e3 * v['tap_from_short_mm']:.0f} um from short, C_in {v['c_in_fF']:.1f} fF, "
                 f"kappa_in / out est {v['kappa_in_est_MHz']:.1f} / {v['kappa_out_est_MHz']:.0f} MHz"
                 if k == 'PF' else f"W {v['meander_width_um']:.0f} um, coupling {v['coupling_um']:g} um at "
                 f"{v['coupling_from_pf_short_mm'][0]:.2f}-{v['coupling_from_pf_short_mm'][1]:.2f} mm from the PF short")
        print(f"{k}: f {v['f_target_GHz']} GHz, target {v['target_mm']:.4f} mm, drawn {v['drawn_mm']:.4f} mm "
              f"({err:+.2f} um); {extra}")
        assert abs(err) < 0.5, f'{k}: drawn length off by {err:.2f} um'


if __name__ == '__main__':
    from gds_tools import export_gds
    d = build_design(variant='sim', channels=True, check_keepout_um=100.0)
    print_report(d)
    for nm, c in d.components.items():
        for pn, pin in c.pins.items():
            if np.linalg.norm(pin.normal) < 1e-9:
                raise ValueError(f'{nm}.{pn}: zero pin normal')
    gds_out = DESIGN / 'gds' / f'{DESIGN.name}_2026-11-27_{REVISION}_sim.gds'
    gds_out.parent.mkdir(parents=True, exist_ok=True)
    export_gds(d, str(gds_out), fab=False, cheese=False)
    print('BUILD_OK:', gds_out)
