"""One FTF set (fluxonium - tunable transmon coupler - fluxonium) on a 10 x 10 mm chip.

Launch pads sit at the PCB port positions in pcbport.py. Port use:
    P7 / P8 : readout feedline (U shape at the top, three coupled-line tees:
              resonators of FX1 6.3 GHz, coupler 6.0 GHz, FX2 6.5 GHz)
    P1 / P6 : flux lines of FX1 / FX2
    P2 / P5 : charge lines of FX1 / FX2
    P3 / P4 : flux line / charge line of the transmon coupler

All CPW routes use a 99 um fillet.

Usage (notebook):  from ftf_layout import build_design
                   design = build_design()                  # full layout (Q3D / HFSS)
                   design = build_design(for_process=True)  # simulation: lumped inductor
                   design = build_design(jj_metal=False)    # junction-free layout (GDS)

GDS (python ftf_layout.py): gds/ftf_layout.gds and gds/ftf_layout_fab.gds (cheesed), both
junction-free: keeps the TC SQUID landing electrodes; omits the JJs and FX arrays.
"""
import os
import sys
from collections import OrderedDict

import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from qiskit_metal.qlibrary.couplers.coupled_line_tee import CoupledLineTee
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.tlines.mixed_path import RouteMixed
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight
from qiskit_metal.qlibrary.sample_shapes.rectangle import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))           # <repo>/pattern design
ROOT = os.path.dirname(HERE)                                # <repo>
GDS_DIR = os.path.join(ROOT, 'gds')
FIG_DIR = os.path.join(ROOT, 'fig', 'layout')   # GDS layout plots
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'pattern design'))   # shared QComponents
from pcbport import pcb_ports
from my_fluxonium import MyFluxonium
from fluxonium import FluxoniumArrayJJ
from my_tmon import MyTmonCoupler
from xmark import add_corner_marks
from gds_tools import export_gds, plot_gds  # shared, skill/gds_export.md, skill/layout_figures.md
from test_junction import TestJunction

# ---- main knobs (mm unless noted) -------------------------------------------
Y0 = -0.8           # y of the FTF row
FX_X = 0.34         # FX1 at -FX_X, FX2 at +FX_X  (Q3D-tuned: Ec_TC = 0.20 GHz)
FEED_Y = 3.0        # horizontal part of the U-shaped feedline
FILLET = '99um'

# ---- readout resonators (lambda/2, open at the tee and at the qubit pad) -----
# f [GHz], kappa [MHz], C_pad [fF] = Maxwell self capacitance of the readout pad + trace
# (Q3D, simulation/results/previous/v4_pads_only_cmat_fF.csv, 15 / 5.909 um CPW)
# C_extra [fF] = resonator-feedline coupling that is not in the isolated-tee fit (the rest
# of the resonator near the feedline): full-circuit Q3D C_1 (12 passes, 0.2 %) minus the fit,
# kappa_v4c_conv.csv (15 / 5.909 um CPW, in the git history, commit ce5674d)
# Design3: FX1 / FX2 scale additionally times f_HFSS / f_target of the 10 um readout pads (eigen_resonators.csv: 6.3472 / 6.5465 GHz, lengths 9.268 / 8.983 mm); not re-checked with HFSS.
# hfss_scale = f_HFSS / f_target of the analytic length (eigenmode, 17 passes,
# simulation/results/previous/eigen_resonators.csv, 15 / 5.909 um CPW); the length is multiplied by it to put f_r on target
# C_extra: full-circuit coupling minus the tee fit, from q3d_kappa d3_kappa_fit19p9 (2026-11-12)
RESONATORS = OrderedDict([
    ('FX1', dict(f=6.3, kappa=2.0, C_pad=24.5, C_extra=-0.328, hfss_scale=6.2892 / 6.3 * 6.3472 / 6.3)),
    ('TC', dict(f=6.0, kappa=3.0, C_pad=48.8, C_extra=-0.745, hfss_scale=6.0259 / 6.0)),
    ('FX2', dict(f=6.5, kappa=2.0, C_pad=24.1, C_extra=-0.244, hfss_scale=6.4892 / 6.5 * 6.5465 / 6.5)),
])
# all CPWs (feedline, resonators, control lines, tees): 10 um centre / 6 um gap
CPW_WIDTH, CPW_GAP = '10 um', '6 um'
CPW_Z0 = 50.9           # Ohm, 10 / 6 um CPW on Si (conformal mapping, eps_eff = 6.225)
CPW_C = 163.4e-12       # F/m
CPW_V = 1.2016e8        # m/s, phase velocity
# tee coupling vs coupling length, Q3D fit (simulation/q3d_tee.py, simulation/results/previous/tee_sweep.csv)
# NOTE: fitted for the old 15 / 5.909 um CPW; re-run q3d_tee.py for 10 / 6 um
CC_SLOPE, CC_OFFSET = 19.896, -1.717     # Cc [fF] = CC_SLOPE * L [mm] + CC_OFFSET (q3d_tee 0.6 / 0.9 / 1.2 mm, CPW 10 / 6 um, 2026-11-12)
TEE_DOWN = 0.15         # mm, tee second-line drop
TEE_GAP = 0.25          # mm, straight feedline between neighbouring tees


def filleted_length(coords, r):
    """Length of a polyline after rounding every corner with radius r (as in the GDS):
    each corner of turning angle t is shortened by 2 r tan(t/2) - r t."""
    c = np.asarray(coords, float)
    L = np.linalg.norm(np.diff(c, axis=0), axis=1).sum()
    for p0, p1, p2 in zip(c[:-2], c[1:-1], c[2:]):
        a, b = p1 - p0, p2 - p1
        t = np.arccos(np.clip(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1))
        L -= 2 * r * np.tan(t / 2) - r * t
    return L


def path_length(design, comp, name=None):
    """Filleted centre-line length [mm] of a route / tee line in the design."""
    t = design.qgeometry.tables['path']
    rows = t[(t['component'] == design.components[comp].id) & (~t['subtract'])]
    if name:
        rows = rows[rows['name'] == name]
    r = rows.iloc[0]
    fil = r['fillet'] if r['fillet'] == r['fillet'] and r['fillet'] else 0.0
    return filleted_length(r.geometry.coords, float(fil))


def coupling_cap(f_ghz, kappa_mhz, Z0=CPW_Z0, Zr=CPW_Z0):
    """Cc [fF] for kappa (lab formula, FXTFX24 cell 42): Qc = pi / (Z0 Zr (Cc w)^2)."""
    w = 2 * np.pi * f_ghz * 1e9
    Qc = f_ghz * 1e3 / kappa_mhz
    return np.sqrt(np.pi / (Qc * Z0 * Zr)) / w * 1e15


def res_length(f_ghz, C_load_ff):
    """lambda/2 length [mm] with end loading: f = v / 2L * (1 - C_load / C_line)."""
    f = f_ghz * 1e9
    L = CPW_V / (2 * f)
    for _ in range(50):
        L = CPW_V / (2 * f) * (1 - C_load_ff * 1e-15 / (CPW_C * L))
    return L * 1e3


for _k, _r in RESONATORS.items():
    _r['Cc'] = coupling_cap(_r['f'], _r['kappa'])                   # fF
    _r['coupling_length'] = round((_r['Cc'] - _r['C_extra'] - CC_OFFSET) / CC_SLOPE, 3)   # mm
    _r['length'] = round(res_length(_r['f'], _r['C_pad'] + _r['Cc']) * _r.get('hfss_scale', 1.0), 3)   # mm, total
RES_FREQ = tuple(RESONATORS[k]['f'] for k in ('FX1', 'FX2'))        # GHz (used by the simulation)

# Q3D-tuned (simulation/results, run s3_x034): Ec_FX = 1.31 GHz, Ec_TC = 0.20 GHz
FX_OPTIONS = dict(
    # smaller, asymmetric pocket: 40 um below the bottom pad, room for the readout pad on top
    pocket_width='260um', pocket_height='360um', pocket_top='210um', pocket_bottom='150um',
    pad_width='152um', pad_height='50um', pad_gap='120um',     # Q3D d3_pw152_ro21 (with stems): Ec_FX 1.0 GHz
    # capacitor stems from both pads to the centre, 10 um wide (was 30 um), 30 um gap between
    # their ends (JJ + array sit in the gap) - added 2026-11-12
    stem_width='10um', stem_gap='30um',
    # 'jj_array_50um' is not in the JJ library yet: add your 50 um array cell under this name
    gds_cell_name='junction_1', gds_cell_inductor='jj_array_50um',
    # JJ array (55 x 10 um, 2 um wires and leads, 2 um JJ gap) and the flux-line tip
    # (15 um from the array end) are the MyFluxonium defaults
    # charge line stops in the ground, 15 um below the pocket edge (~9 um ground strip
    # after its end gap); the pocket edge is only 40 um from the bottom pad now
    charge_line_options=Dict(make_cl=True, cl_sep='15um', cl_length='100um'),
    # readout coupling pad: 21 x 40 um (was 10 um wide: Ec 1 GHz needs a larger C_sigma, so a
    # larger pad keeps beta ~0.1); pad_sep is set per qubit in FX_PER_QUBIT for g = 80 MHz
    readout_line_options=Dict(make_rol=True, pad_width='21um', pad_height='40um', pad_sep='10um'),
)
# per-qubit overrides (merged into FX_OPTIONS). The readout gap differs because the
# resonators sit at different frequencies (V_zpf grows with f_r): g = 80 MHz for both (Q3D d3_pw152_ro21: 81.7 / 80.2 MHz).
FX_PER_QUBIT = {
    'FX1': dict(readout_line_options=Dict(make_rol=True, pad_width='21um', pad_height='40um', pad_sep='10um')),
    'FX2': dict(readout_line_options=Dict(make_rol=True, pad_width='21um', pad_height='40um', pad_sep='11.45um')),
}
PADDLE_SEP = 0.020      # mm, gap between a coupler paddle and the FX bottom pad
TMON_OPTIONS = dict(
    # pad_arm_length and the coupler height are set in build_design from the FX pads
    pad_arm_width='15um',
    paddle_width='30um', paddle_gap='30um',
    pad_head_length='158um', pad_gap='30um', jj_gap='40um',      # Q3D-tuned (v4): Ec_TC 0.2 GHz
    squid_loop_shape='rectangle', squid_loop_width='20um', squid_loop_height='14um',
    squid_contact_overlap='3um',  # nominal overlap onto stem and ground for mask preview
    squid_shift_toward_stem='2um',  # leaves 1 um of overlap onto the ground plane
    gds_cell_name='junction_t0',
    # L-shaped flux line shifted 25 um to the left (leg at x = +32.5 um)
    flux_bias_line_options=Dict(make_fbl=True, fbl_sep='10um', fbl_length='100um', fbl_x='-25um'),
    # right of the stem, between the flux line and FX2's pocket (x = +160 um was too close to
    # FX2). Tip 250 um below the bar etch: Q3D (tc5_x120_cl250) 0.051 fF -> T1 limit 0.93 ms
    charge_line_options=Dict(make_cl=True, cl_x='120um', cl_sep='250um'),
    # Initial 100 MHz target for the bottom-pad layout (transmon convention, EJ/Ec = 50).
    readout_line_options=Dict(make_rol=True, ro_width='150um', ro_height='20um', ro_sep='5um'),
)
LP = dict(chip='main', lead_length='30um', pad_width='200um', pad_height='200um', pad_gap='76.5um')
# test-junction structures in the four chip corners (as in the lab's layout3 ftf.gds), mm:
# an etched pocket with two pads; the junction between them is left for the fab to draw
TEST_JJ = dict(x=3.95, y=4.25,                       # centres at (+-x, +-y)
               pocket=(0.800, 0.300),                # pocket width, height
               pad_left=(0.225, 0.129), pad_right=(0.295, 0.129), pad_gap=0.211)
# die marks (shared XMark, pattern design/xmark.py): metal crosses (250 x 10 um bars) in etched 400 um
# square windows, centred MARK_INSET mm
# in from the four chip corners; the die ID (etched) sits left of the LR window, below the cross line.
# Die ID: build_design(die_id='FTF07'), FTF00 ... FTF99
MARK_INSET = 0.35
MARK = dict(arm_length='250um', arm_width='10um', window_size='400um', label_size='100um', label_gap='30um')
DIE_ID = 'FTF00'


def route(design, name, start, end, anchors, ss='0um', es='0um'):
    """CPW through explicit anchor points (mm)."""
    return RouteAnchors(design, name, Dict(
        chip='main', fillet=FILLET,
        lead=Dict(start_straight=ss, end_straight=es),
        anchors=OrderedDict({i: np.array(a, float) for i, a in enumerate(anchors)}),
        pin_inputs=Dict(start_pin=Dict(component=start[0], pin=start[1]),
                        end_pin=Dict(component=end[0], pin=end[1]))))


def _mm(v):
    """'240um' / '0.24mm' / 0.24 -> 0.24 (mm)."""
    v = str(v).strip()
    if v.endswith('um'):
        return float(v[:-2]) / 1000
    if v.endswith('mm'):
        return float(v[:-2])
    return float(v)


def build_design(for_process=False, fx_options=None, tmon_options=None, fx_x=None,
                 fx1_options=None, fx2_options=None, marks=None, die_id=DIE_ID, jj_metal=True,
                 include_test_junctions=True, include_tc_landing_pads=False):
    """for_process=True uses FluxoniumArrayJJ (simulation model: two lumped array arms, a lumped main JJ, and a grounded flux line)
    for Q3D / HFSS; False gives the fabrication layout.
    fx_options / tmon_options override FX_OPTIONS / TMON_OPTIONS, fx1_options /
    fx2_options override a single fluxonium, and fx_x overrides FX_X (for tuning).
    include_test_junctions: add the four corner test structures when True.
    include_tc_landing_pads: retain the base-metal SQUID electrodes when jj_metal=False.
    marks: corner crosses + die ID (default: only in the fabrication layout, not for_process).
    jj_metal=False: junction-free layout for the GDS - the FX pads without leads / JJ stubs /
    array end and the coupler stem without the SQUID block (see export_gds)."""
    FX_X_ = FX_X if fx_x is None else fx_x
    FX = FluxoniumArrayJJ if for_process else MyFluxonium
    fx_opt = dict(FX_OPTIONS, **(fx_options or {}))
    tm_opt = dict(TMON_OPTIONS, **(tmon_options or {}))
    if for_process:
        tm_opt.setdefault('flux_as_ground', 'True')
        fx_opt.setdefault('flux_as_ground', 'True')
    if not jj_metal:
        fx_opt['jj_metal'] = 'False'
        tm_opt['jj_metal'] = 'True' if include_tc_landing_pads else 'False'
    # coupler geometry follows the FX pads: paddles PADDLE_SEP from the bottom pads,
    # bar at the height of the bottom-pad centre
    pw, ph, pg = _mm(fx_opt['pad_width']), _mm(fx_opt['pad_height']), _mm(fx_opt['pad_gap'])
    paddle_w = _mm(tm_opt['paddle_width'])
    tm_opt.setdefault('pad_arm_length', f'{2 * (FX_X_ - pw / 2 - PADDLE_SEP - paddle_w) * 1000:.3f}um')
    tm_opt.setdefault('paddle_height', fx_opt['pad_height'])
    y_tc = Y0 - pg / 2 - ph / 2
    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon'
    design.variables['cpw_width'] = CPW_WIDTH
    design.variables['cpw_gap'] = CPW_GAP
    design.chips.main.size.size_x = '10mm'
    design.chips.main.size.size_y = '10mm'
    design.chips.main.size.size_z = '-525um'

    # ---- launch pads at the PCB ports ----
    for name, opt in pcb_ports.items():
        LaunchpadWirebond(design, 'LP' + name, options=dict(LP, **opt))

    # ---- test-junction pockets in the chip corners (pads only, junction site empty) ----
    if include_test_junctions:
        t = TEST_JJ
        (wl, hl), (wr, hr) = t['pad_left'], t['pad_right']
        for name, sx, sy in (('LL', -1, -1), ('LR', 1, -1), ('UL', -1, 1), ('UR', 1, 1)):
            TestJunction(design, f'testjj_{name}', options=dict(
                chip='main', pos_x=f"{sx * t['x']}mm", pos_y=f"{sy * t['y']}mm",
                pocket_width=f"{t['pocket'][0]}mm", pocket_height=f"{t['pocket'][1]}mm",
                pad_l_width=f'{wl}mm', pad_l_height=f'{hl}mm', pad_r_width=f'{wr}mm', pad_r_height=f'{hr}mm',
                pad_gap=f"{t['pad_gap']}mm"))

    # ---- die marks: corner crosses + die ID ----
    if (not for_process) if marks is None else marks:
        add_corner_marks(design, die_id=die_id, inset=MARK_INSET, label_corner='LR',
                         label_left='True', **MARK)

    # ---- qubits ----
    fx1 = {**fx_opt, **FX_PER_QUBIT['FX1'], **(fx1_options or {})}
    fx2 = {**fx_opt, **FX_PER_QUBIT['FX2'], **(fx2_options or {})}
    FX(design, 'FX1', options=dict(pos_x=f'{-FX_X_}mm', pos_y=f'{Y0}mm', **fx1))
    FX(design, 'FX2', options=dict(pos_x=f'{FX_X_}mm', pos_y=f'{Y0}mm', mirror='True', **fx2))
    # coupler bar at the height of the FX bottom pads, flat paddles face the pad edges
    MyTmonCoupler(design, 'TC', options=dict(pos_x='0mm', pos_y=f'{y_tc}mm', **tm_opt))

    # ---- feedline with three coupled-line tees (FX1, TC, FX2 from left to right) ----
    tee = dict(chip='main', prime_width='cpw_width', prime_gap='cpw_gap',
               second_width='cpw_width', second_gap='cpw_gap', coupling_space='3um',
               down_length=f'{TEE_DOWN}mm', fillet=FILLET, open_termination=True)
    widths = [2 * r['coupling_length'] for r in RESONATORS.values()]      # prime length of each tee
    x = -(sum(widths) + TEE_GAP * (len(widths) - 1)) / 2
    drops = {}
    for (name, r), wdt in zip(RESONATORS.items(), widths):
        cx, cl = x + wdt / 2, r['coupling_length']
        # second line drops at cx + cl/2 (mirror False) or cx - cl/2 (True): take the one nearer 0
        mirror = abs(cx - cl / 2) < abs(cx + cl / 2)
        drops[name] = cx - cl / 2 if mirror else cx + cl / 2
        CoupledLineTee(design, f'tee_{name}', options=dict(
            pos_x=f'{cx:.4f}mm', pos_y=f'{FEED_Y}mm', mirror=mirror, coupling_length=f'{cl}mm', **tee))
        x += wdt + TEE_GAP
    x_left = -(sum(widths) + TEE_GAP * (len(widths) - 1)) / 2
    feed_x = abs(x_left) + 0.3                                            # corners of the feedline U
    names = list(RESONATORS)
    route(design, 'feed_in', ('LPP8', 'tie'), (f'tee_{names[0]}', 'prime_start'),
          [(-1.1, 3.9), (-feed_x, 3.9), (-feed_x, FEED_Y)])
    for a_, b_ in zip(names[:-1], names[1:]):
        RouteStraight(design, f'feed_{a_}_{b_}', Dict(chip='main', pin_inputs=Dict(
            start_pin=Dict(component=f'tee_{a_}', pin='prime_end'),
            end_pin=Dict(component=f'tee_{b_}', pin='prime_start'))))
    route(design, 'feed_out', (f'tee_{names[-1]}', 'prime_end'), ('LPP7', 'tie'),
          [(feed_x, FEED_Y), (feed_x, 3.9), (1.1, 3.9)])

    # ---- lambda/2 readout resonators: tee (open end) -> meander -> qubit readout pad ----
    # RouteMixed: meander straight below the tee down to y_jog, straight jog sideways
    # to the qubit's x, straight down into the readout pin; total_length covers all of it
    jog_up = {'FX1': 0.30, 'FX2': 0.30, 'TC': 0.75}     # mm, height of the jog above the qubit pin
    for name, r in RESONATORS.items():
        pin = design.components[name].pins['readout_line'].middle
        y_jog = pin[1] + jog_up[name]
        RouteMixed(design, f'res_{name}', Dict(
            chip='main', total_length=f'{r["length"] - r["coupling_length"] - TEE_DOWN:.4f}mm',
            fillet=FILLET, hfss_wire_bonds=False,
            meander=Dict(spacing='230um', asymmetry='0um'),
            lead=Dict(start_straight='200um', end_straight='100um'),
            anchors=OrderedDict({0: np.array([drops[name], y_jog]), 1: np.array([pin[0], y_jog])}),
            between_anchors=OrderedDict({0: 'M', 1: 'S', 2: 'S'}),
            pin_inputs=Dict(start_pin=Dict(component=f'tee_{name}', pin='second_end'),
                            end_pin=Dict(component=name, pin='readout_line'))))
        # RouteMixed misses total_length and fillets shorten the line: iterate on the
        # real (filleted) length of tee second line + route until it hits r['length']
        res = design.components[f'res_{name}']
        tee_len = path_length(design, f'tee_{name}', 'second_cpw')
        for _ in range(8):
            err = r['length'] - (tee_len + path_length(design, f'res_{name}'))
            if abs(err) < 0.002:
                break
            res.options.total_length = f'{float(res.options.total_length[:-2]) + err:.4f}mm'
            res.rebuild()

    # ---- flux lines (FX: side U-loop, coupler: bottom) ----
    yf = design.components['FX1'].pins['flux_line'].middle[1]   # height of the FX flux pins
    route(design, 'fl_FX1', ('FX1', 'flux_line'), ('LPP1', 'tie'), [(-3.8, yf), (-3.8, 1.6)])
    route(design, 'fl_FX2', ('FX2', 'flux_line'), ('LPP6', 'tie'), [(3.8, yf), (3.8, 1.6)])
    xt = design.components['TC'].pins['flux_bias_line'].middle[0]
    route(design, 'fl_TC', ('TC', 'flux_bias_line'), ('LPP3', 'tie'), [(xt, -3.9), (-1.1, -3.9)])

    # ---- charge lines ----
    route(design, 'cl_FX1', ('FX1', 'charge_line'), ('LPP2', 'tie'), [(-FX_X_, -1.6)])
    route(design, 'cl_FX2', ('FX2', 'charge_line'), ('LPP5', 'tie'), [(FX_X_, -1.6)])
    xc = design.components['TC'].pins['charge_line'].middle[0]
    route(design, 'cl_TC', ('TC', 'charge_line'), ('LPP4', 'tie'), [(xc, -3.6), (1.1, -3.6)])
    return design


def plot_design(design, filename, window=None, dpi=130):
    """Matplotlib rendering without the Qt GUI (metal = grey, etched = white)."""
    import shapely
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from shapely.geometry import box

    names = {v.id: k for k, v in design.components.items()}
    sub, pos, jj = [], [], []
    for _, r in design.qgeometry.tables['poly'].iterrows():
        (sub if r['subtract'] else pos).append(r.geometry)
    for _, r in design.qgeometry.tables['path'].iterrows():
        g = r.geometry.buffer(r['width'] / 2, cap_style=2, join_style=2)
        (sub if r['subtract'] else pos).append(g)
    for _, r in design.qgeometry.tables['junction'].iterrows():
        jj.append(r.geometry.buffer(max(r['width'], 0.002) / 2, cap_style=2))
    sx = float(design.chips.main.size.size_x.replace('mm', '')) / 2
    sy = float(design.chips.main.size.size_y.replace('mm', '')) / 2
    metal = box(-sx, -sy, sx, sy).difference(shapely.union_all(sub)).union(shapely.union_all(pos))

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_facecolor('white')

    def fill(g, **kw):
        for p in getattr(g, 'geoms', [g]):
            if p.is_empty:
                continue
            ax.fill(*p.exterior.xy, **kw)
            for h in p.interiors:
                ax.fill(*h.xy, color='white', lw=0)
    fill(metal, color='#9aa3ad', lw=0)
    for g in jj:
        fill(g, color='red', lw=0)
    for cname, c in design.components.items():
        for pn, pin in c.pins.items():
            if cname.startswith(('res', 'fl_', 'cl_', 'feed', 'tee')):
                continue
            ax.annotate(f'{cname}.{pn}', pin.middle, fontsize=6, color='navy')
    for name, opt in pcb_ports.items():
        ax.annotate(name, (float(opt['pos_x'][:-2]), float(opt['pos_y'][:-2])), fontsize=11,
                    color='darkred', weight='bold', ha='center', va='center')
    ax.set_aspect('equal')
    if window:
        ax.set_xlim(window[0], window[1]); ax.set_ylim(window[2], window[3])
    else:
        ax.set_xlim(-sx, sx); ax.set_ylim(-sy, sy)
    ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
    fig.tight_layout()
    fig.savefig(filename, dpi=dpi)
    plt.close(fig)


if __name__ == '__main__':
    import argparse
    import warnings
    warnings.filterwarnings('ignore')
    parser = argparse.ArgumentParser(description='Export Design3 junction-free GDS.')
    parser.add_argument('--test-junctions', choices=('yes', 'no'))
    args = parser.parse_args()
    choice = args.test_junctions
    if choice is None:
        choice = input('Include corner test junction structures? [yes/no] ').strip().lower()
        if choice not in ('yes', 'no'):
            parser.error('choose yes or no before exporting')
    export_options = dict(jj_metal=False, include_tc_landing_pads=True,
                          include_test_junctions=choice == 'yes')
    gds_file = os.path.join(GDS_DIR, 'ftf_layout.gds')
    fab_file = os.path.join(GDS_DIR, 'ftf_layout_fab.gds')
    d = build_design(**export_options)
    export_gds(d, gds_file)
    export_gds(build_design(**export_options), fab_file, fab=True)
    plot_gds(gds_file, os.path.join(FIG_DIR, 'ftf_layout_chip.png'), ports=pcb_ports)
    plot_gds(gds_file, os.path.join(FIG_DIR, 'ftf_layout_zoom.png'), window=(-1.0, 1.0, Y0 - 0.55, Y0 + 0.55))
    plot_gds(gds_file, os.path.join(FIG_DIR, 'ftf_layout_fx1.png'), window=(-0.95, -0.3, Y0 - 0.3, Y0 + 0.3))
    # FX1 JJ + array site (empty in the GDS) and the flux-line tip
    plot_gds(gds_file, os.path.join(FIG_DIR, 'ftf_layout_fxjj.png'), window=(-0.48, -0.32, Y0 - 0.045, Y0 + 0.045), dpi=140)
    print('wrote', gds_file, fab_file)
