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
junction-free: no JJ, no array, no leads / JJ stubs / array end, no SQUID block.
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
from project_paths import pattern_design_paths as _pattern_paths
import os
import sys
from collections import OrderedDict

import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight
from qiskit_metal.qlibrary.sample_shapes.rectangle import Rectangle

HERE = _S('Design1_FTF_v1/layout')           # <repo>/pattern design
ROOT = _S('Design1_FTF_v1')                                # <repo>
GDS_DIR = _S('Design1_FTF_v1/gds')
FIG_DIR = _S('Design1_FTF_v1/fig/layout')   # GDS layout plots
sys.path.insert(0, HERE)
sys.path[1:1] = _pattern_paths()   # shared QComponents
from bondpad import add_bond_pad
from wafer_pcbport import wafer_pcb_ports as pcb_ports
from fluxonium_nthu_dolan import MyFluxoniumDolan, FluxoniumArrayJJDolan
from tmon_coupler_dolan_nthu import TmonCouplerDolanNthu
from fluxline import B2_FLUX_OPTIONS
from xmark import add_corner_marks
from fluxonium_nthu_dolan_fab import MyFluxoniumDolanFab
from tmon_coupler_dolan_nthu_fab import TmonCouplerDolanNthuFab
from gds_tools import export_gds, plot_gds  # shared, skill/gds_export.md, skill/layout_figures.md

# ---- main knobs (mm unless noted) -------------------------------------------
Y0 = -0.8           # y of the FTF row
FX_X = 0.34         # FX1 at -FX_X, FX2 at +FX_X  (Q3D-tuned: Ec_TC = 0.20 GHz)
FEED_Y = 1.35       # horizontal part of the U-shaped feedline
FEED_X = 2.00       # half-width of the feedline U
FILLET = '99um'     # feedline & control lines
FILLET_RES = '70um' # Design7-style lambda/4 hanger resonators

# ---- readout resonators (Design7-style lambda/4 hanger, shorted at feedline, open at qubit pad) ----
# Analytic target length: v_ph / (4 f) - load_shorten_mm (unsimulated after 2026-11-22 switch to lambda/4 hanger
# and 10 um readout coupling pads; needs Q3D + HFSS revalidation).
RESONATORS = OrderedDict([
    ('FX1', dict(f=6.3, x_c=-1.25, load_shorten_mm=0.100)),
    ('TC',  dict(f=6.0, x_c=0.00,  load_shorten_mm=0.150)),
    ('FX2', dict(f=6.5, x_c=+1.25, load_shorten_mm=0.100)),
])
# Design7 hanger-resonator geometry (10 / 6 um CPW, 8 um ground strip, 300 um trench gap -> 322 um pitch)
HANGER = dict(
    coupling_length_mm=0.300,
    coupling_ground_strip_mm=0.008,
    meander_pitch_mm=0.322,
    lead_y_mm=0.322,
    fillet_mm=0.070,
)
# all CPWs (feedline, resonators, control lines): 10 um centre / 6 um gap
CPW_WIDTH, CPW_GAP = '10 um', '6 um'
CPW_W_MM, CPW_G_MM = 0.010, 0.006
CPW_Z0 = 50.9           # Ohm, 10 / 6 um CPW on Si (conformal mapping, eps_eff = 6.225)
CPW_C = 163.4e-12       # F/m
CPW_V = 1.2016e8        # m/s, phase velocity


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
    """Filleted centre-line length [mm] of a route in the design."""
    t = design.qgeometry.tables['path']
    rows = t[(t['component'] == design.components[comp].id) & (~t['subtract'])]
    if name:
        rows = rows[rows['name'] == name]
    r = rows.iloc[0]
    fil = r['fillet'] if r['fillet'] == r['fillet'] and r['fillet'] else 0.0
    return filleted_length(r.geometry.coords, float(fil))


def quarter_wave_length(f_ghz, load_shorten_mm=0.100):
    """Analytic lambda/4 CPW length [mm] minus capacitive end-load allowance."""
    return (CPW_V / (4.0 * f_ghz * 1e9)) * 1e3 - load_shorten_mm


for _k, _r in RESONATORS.items():
    _r['coupling_length'] = HANGER['coupling_length_mm']
    _r['length'] = round(quarter_wave_length(_r['f'], _r['load_shorten_mm']), 4)   # mm, filleted target
RES_FREQ = tuple(RESONATORS[k]['f'] for k in ('FX1', 'FX2'))        # GHz (used by the simulation)

# Q3D-tuned (simulation/results, run s3_x034): Ec_FX = 1.31 GHz, Ec_TC = 0.20 GHz
FX_OPTIONS = dict(
    # smaller, asymmetric pocket: 40 um below the bottom pad, room for the readout pad on top
    pocket_width='260um', pocket_height='360um', pocket_top='210um', pocket_bottom='150um',
    pad_width='125um', pad_height='50um', pad_gap='120um',     # Q3D-tuned (v4): Ec_FX 1.3 GHz (before the stems)
    # capacitor stems from both pads to the centre, 30 um wide, 30 um gap between their ends
    # (JJ + array sit in the gap) - added 2026-11-12, Ec / J / g not re-simulated yet
    stem_width='30um', stem_gap='30um',
    # 'jj_array_50um' is not in the JJ library yet: add your 50 um array cell under this name
    gds_cell_name='junction_1', gds_cell_inductor='jj_array_50um',
    # JJ array (55 x 10 um, 2 um wires and leads, 2 um JJ gap) and the flux-line tip
    # (15 um from the array end) are the MyFluxoniumDolan defaults
    # charge line stops in the ground, 15 um below the pocket edge (~9 um ground strip
    # after its end gap); the pocket edge is only 40 um from the bottom pad now
    charge_line_options=Dict(make_cl=True, cl_sep='15um', cl_length='100um'),
    # readout coupling pad: 10 um thickness (matching Design7_Tunable_Transmon); pad_sep set per qubit in FX_PER_QUBIT
    readout_line_options=Dict(make_rol=True, pad_width='100um', pad_height='10um', pad_sep='30um'),
    # Design1_FTF_v1 flux line (M = 0.575 pH, Design1_FTF_v1/notes/mutual_inductance_method.md). Explicit here because the shared
    # MyFluxoniumDolan (SQL HFSS/pattern design/fluxonium_nthu_dolan.py) now carries the Design3_FTF_NTHU defaults (80 / 11 um)
    flux_line_options=Dict(make_fl=True, fl_height='50um', fl_sep='15um', cpw_width='10um', cpw_gap='6um'),
)
# per-qubit overrides (merged into FX_OPTIONS). The readout gap differs because the
# resonators sit at different frequencies (V_zpf grows with f_r): g = 150 MHz for both.
FX_PER_QUBIT = {
    'FX1': dict(readout_line_options=Dict(make_rol=True, pad_width='100um', pad_height='10um', pad_sep='19.75um')),
    'FX2': dict(readout_line_options=Dict(make_rol=True, pad_width='100um', pad_height='10um', pad_sep='21.2um')),
}
PADDLE_SEP = 0.020      # mm, gap between a coupler paddle and the FX top pad
TMON_OPTIONS = dict(
    # pad_arm_length and the coupler height are set in build_design from the FX pads
    pad_arm_width='15um',
    paddle_width='30um', paddle_gap='30um',
    pad_head_width='15um', pad_head_length='158um', pad_gap='30um', jj_gap='40um',
    # Design3_FTF_NTHU rectangular SQUID + 7 um ground tail
    squid_loop_shape='rectangle', squid_loop_width='20um', squid_loop_height='24um',
    squid_ground_tail='7um',
    gds_cell_name='junction_t0',
    # Design3_FTF_NTHU / Design6 B2 short L-shaped flux line (fbl_sep 10.5 um, >= 4 um ground strip)
    flux_bias_line_options=Dict(B2_FLUX_OPTIONS, make_fbl=True),
    # right of the stem, between the flux line and FX2's pocket (x = +160 um was too close to
    # FX2). Tip 250 um below the bar etch: Q3D (tc5_x120_cl250) 0.051 fF -> T1 limit 0.93 ms
    charge_line_options=Dict(make_cl=True, cl_x='120um', cl_sep='250um'),
    # readout coupling pad: 10 um thickness (matching Design7_Tunable_Transmon); needs Q3D revalidation
    readout_line_options=Dict(make_rol=True, ro_width='273um', ro_height='10um', ro_sep='5um'),
)
# test-junction structures in the four chip corners (as in the lab's layout3 ftf.gds), mm:
# an etched pocket with two pads; the junction between them is left for the fab to draw
TEST_JJ = dict(x=3.95, y=4.25,                       # centres at (+-x, +-y)
               pocket=(0.800, 0.300),                # pocket width, height
               pad_left=(0.225, 0.129), pad_right=(0.295, 0.129), pad_gap=0.211)
# die marks (shared XMark, pattern design/xmark.py): metal crosses (250 x 10 um bars) in etched 400 um
# square windows, centred MARK_INSET mm
# in from the four chip corners; the die ID (etched) sits right of the UL window, below the cross line.
# Die ID: build_design(die_id='FTF07'), FTF00 ... FTF99
MARK_INSET = 0.35
MARK = dict(arm_length='250um', arm_width='10um', window_size='400um', label_size='150um', label_gap='30um')   # LABEL_SIZE_UM in xmark.py
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
                 fab_contacts=False, include_tc_landing_pads=False, include_test_junctions=True):
    """for_process=True uses FluxoniumArrayJJDolan (simulation model: array = hairpin loop of two
    lumped inductors, JJ = short lumped segment; U-loop flux line drawn as ground)
    for Q3D / HFSS; False gives the fabrication layout.
    fx_options / tmon_options override FX_OPTIONS / TMON_OPTIONS, fx1_options /
    fx2_options override a single fluxonium, and fx_x overrides FX_X (for tuning).
    marks: corner crosses + die ID (default: only in the fabrication layout, not for_process).
    jj_metal=False: junction-free layout for the GDS - the FX pads without leads / JJ stubs /
    array end and the coupler stem without the SQUID block (see export_gds)."""
    FX_X_ = FX_X if fx_x is None else fx_x
    if fab_contacts and (for_process or jj_metal or not include_tc_landing_pads):
        raise ValueError('Fab contacts require a junction-free mask with TC landing pads.')
    FX = FluxoniumArrayJJDolan if for_process else (MyFluxoniumDolanFab if fab_contacts else MyFluxoniumDolan)
    TC = TmonCouplerDolanNthuFab if fab_contacts else TmonCouplerDolanNthu
    fx_opt = dict(FX_OPTIONS, **(fx_options or {}))
    tm_opt = dict(TMON_OPTIONS, **(tmon_options or {}))
    if for_process:
        tm_opt.setdefault('flux_as_ground', 'True')
        fx_opt.setdefault('flux_as_ground', 'True')
    if not jj_metal:
        fx_opt['jj_metal'] = 'False'
        tm_opt['jj_metal'] = 'True' if include_tc_landing_pads else 'False'
        if include_tc_landing_pads:
            tm_opt.update(squid_contact_overlap='3um', squid_shift_toward_stem='2um')
    # coupler geometry follows the FX pads: paddles PADDLE_SEP from the top pads,
    # bar at the height of the top-pad centre
    pw, ph, pg = _mm(fx_opt['pad_width']), _mm(fx_opt['pad_height']), _mm(fx_opt['pad_gap'])
    paddle_w = _mm(tm_opt['paddle_width'])
    tm_opt.setdefault('pad_arm_length', f'{2 * (FX_X_ - pw / 2 - PADDLE_SEP - paddle_w) * 1000:.3f}um')
    tm_opt.setdefault('paddle_height', fx_opt['pad_height'])
    y_tc = Y0 + pg / 2 + ph / 2
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
        add_bond_pad(design, 'LP' + name, opt)

    # ---- test-junction pockets in the chip corners (pads only, junction site empty) ----
    t = TEST_JJ
    (wl, hl), (wr, hr), gap = t['pad_left'], t['pad_right'], t['pad_gap']
    x_l = -(wl + gap + wr) / 2 + wl / 2                  # pad centres relative to the pocket centre
    x_r = (wl + gap + wr) / 2 - wr / 2
    corners = (('LL', -1, -1), ('LR', 1, -1), ('UL', -1, 1), ('UR', 1, 1))
    for name, sx, sy in corners if include_test_junctions else ():
        x0, y0 = sx * t['x'], sy * t['y']
        Rectangle(design, f'testjj_{name}_pocket', options=dict(
            chip='main', width=f"{t['pocket'][0]}mm", height=f"{t['pocket'][1]}mm",
            subtract='True', pos_x=f'{x0}mm', pos_y=f'{y0}mm'))
        Rectangle(design, f'testjj_{name}_padL', options=dict(
            chip='main', width=f'{wl}mm', height=f'{hl}mm', pos_x=f'{x0 + x_l}mm', pos_y=f'{y0}mm'))
        Rectangle(design, f'testjj_{name}_padR', options=dict(
            chip='main', width=f'{wr}mm', height=f'{hr}mm', pos_x=f'{x0 + x_r}mm', pos_y=f'{y0}mm'))

    # ---- die marks: corner crosses + die ID ----
    if (not for_process) if marks is None else marks:
        add_corner_marks(design, die_id=die_id, inset=MARK_INSET, **MARK)

    # ---- qubits ----
    fx1 = dict(fx_opt, **FX_PER_QUBIT['FX1'], **(fx1_options or {}))
    fx2 = dict(fx_opt, **FX_PER_QUBIT['FX2'], **(fx2_options or {}))
    FX(design, 'FX1', options=dict(pos_x=f'{-FX_X_}mm', pos_y=f'{Y0}mm', **fx1))
    FX(design, 'FX2', options=dict(pos_x=f'{FX_X_}mm', pos_y=f'{Y0}mm', mirror='True', **fx2))
    # coupler bar at the height of the FX top pads, flat paddles face the pad edges
    TC(design, 'TC', options=dict(pos_x='0mm', pos_y=f'{y_tc}mm', **tm_opt))

    # ---- U-shaped straight feedline (P8 -> P7, horizontal segment at y = FEED_Y) ----
    xp8, xp7 = float(pcb_ports['P8']['pos_x'][:-2]), float(pcb_ports['P7']['pos_x'][:-2])
    route(design, 'feedline', ('LPP8', 'tie'), ('LPP7', 'tie'),
          [(xp8, 3.6), (-FEED_X, 3.6), (-FEED_X, FEED_Y),
           (FEED_X, FEED_Y), (FEED_X, 3.6), (xp7, 3.6)])

    # ---- Design7-style lambda/4 hanger readout resonators: ShortToGround -> serpentine -> qubit readout_line ----
    l_coup = HANGER['coupling_length_mm']
    strip = HANGER['coupling_ground_strip_mm']
    pitch = HANGER['meander_pitch_mm']
    lead_y = HANGER['lead_y_mm']
    r_fil = HANGER['fillet_mm']
    min_seg = 2.0 * r_fil
    y_coup = FEED_Y - (CPW_W_MM + 2.0 * CPW_G_MM + strip)
    y0 = y_coup - lead_y
    y1 = y0 - pitch
    y2 = y1 - pitch
    y3 = y2 - pitch

    jog_up = {'FX1': 0.25, 'FX2': 0.25, 'TC': 0.0}
    for name, r in RESONATORS.items():
        pin = design.components[name].pins['readout_line'].middle
        px, py = float(pin[0]), float(pin[1])
        x_c = float(r['x_c'])
        has_jog = abs(x_c - px) > 1e-6
        y_jog = py + jog_up[name]
        dy_down = (y3 - y_jog) if has_jog else (y3 - py)
        if dy_down < min_seg - 1e-9:
            raise ValueError(f'res_{name}: vertical exit {dy_down*1e3:.1f} um < 2*fillet ({min_seg*1e3:.1f} um)')
        if has_jog and jog_up[name] < min_seg - 1e-9:
            raise ValueError(f'res_{name}: jog_up {jog_up[name]*1e3:.1f} um < 2*fillet ({min_seg*1e3:.1f} um)')

        n_bends = 11 if has_jog else 9
        fillet_loss = n_bends * r_fil * (2.0 - np.pi / 2.0)
        fixed_poly = 1.5 * l_coup + lead_y + 3.0 * pitch + dy_down + (abs(x_c - px) + jog_up[name] if has_jog else 0.0)
        w_m = (r['length'] + fillet_loss - fixed_poly) / 3.0
        if w_m / 2.0 < min_seg - 1e-9:
            raise ValueError(f'res_{name}: half meander width {w_m*500:.1f} um < 2*fillet ({min_seg*1e3:.1f} um)')

        x_start = x_c - l_coup / 2.0
        x_turn = x_c + l_coup / 2.0
        x_left = x_c - w_m / 2.0
        x_right = x_c + w_m / 2.0

        short_name = f'res_{name}_short'
        ShortToGround(design, short_name, options=dict(
            pos_x=f'{x_start:.9f}mm', pos_y=f'{y_coup:.9f}mm', orientation='180',
            width=CPW_WIDTH, gap=CPW_GAP,
        ))

        pts = [
            (x_turn, y_coup),
            (x_turn, y0),
            (x_left, y0),
            (x_left, y1),
            (x_right, y1),
            (x_right, y2),
            (x_left, y2),
            (x_left, y3),
            (x_c, y3),
        ]
        if has_jog:
            pts += [(x_c, y_jog), (px, y_jog)]

        RouteAnchors(design, f'res_{name}', Dict(
            chip='main', fillet=FILLET_RES, hfss_wire_bonds=False,
            lead=Dict(start_straight='0um', end_straight='0um'),
            anchors=OrderedDict({i: np.array(pt, float) for i, pt in enumerate(pts)}),
            pin_inputs=Dict(start_pin=Dict(component=short_name, pin='short'),
                            end_pin=Dict(component=name, pin='readout_line')),
        ))
        actual_len = path_length(design, f'res_{name}')
        if abs(actual_len - r['length']) > 0.002:
            raise AssertionError(f'res_{name}: actual {actual_len:.4f} mm != target {r["length"]:.4f} mm')

    # ---- flux lines (FX: side U-loop, coupler: bottom) ----
    yf = design.components['FX1'].pins['flux_line'].middle[1]   # height of the FX flux pins
    route(design, 'fl_FX1', ('FX1', 'flux_line'), ('LPP1', 'tie'), [(-3.8, yf), (-3.8, float(pcb_ports['P1']['pos_y'][:-2]))])
    route(design, 'fl_FX2', ('FX2', 'flux_line'), ('LPP6', 'tie'), [(3.8, yf), (3.8, float(pcb_ports['P6']['pos_y'][:-2]))])
    xt = design.components['TC'].pins['flux_bias_line'].middle[0]
    route(design, 'fl_TC', ('TC', 'flux_bias_line'), ('LPP3', 'tie'), [(xt, -3.9), (float(pcb_ports['P3']['pos_x'][:-2]), -3.9)])

    # ---- charge lines ----
    route(design, 'cl_FX1', ('FX1', 'charge_line'), ('LPP2', 'tie'), [(-FX_X_, float(pcb_ports['P2']['pos_y'][:-2]))])
    route(design, 'cl_FX2', ('FX2', 'charge_line'), ('LPP5', 'tie'), [(FX_X_, float(pcb_ports['P5']['pos_y'][:-2]))])
    xc = design.components['TC'].pins['charge_line'].middle[0]
    route(design, 'cl_TC', ('TC', 'charge_line'), ('LPP4', 'tie'), [(xc, -3.6), (float(pcb_ports['P4']['pos_x'][:-2]), -3.6)])
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

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test-junctions', choices=('yes', 'no'))
    args = parser.parse_args()
    test_junctions = args.test_junctions
    while test_junctions not in ('yes', 'no'):
        test_junctions = input('Include corner test-junction pads? [yes/no] ').strip().lower()
    warnings.filterwarnings('ignore')
    gds_file = _S('Design1_FTF_v1/gds/ftf_layout.gds')
    fab_file = _S('Design1_FTF_v1/gds/ftf_layout_fab.gds')
    mask_options = dict(jj_metal=False, include_tc_landing_pads=True,
                        include_test_junctions=test_junctions == 'yes')
    export_gds(build_design(**mask_options), gds_file)
    export_gds(build_design(**mask_options, fab_contacts=True), fab_file, fab=True)
    plot_gds(gds_file, _S('Design1_FTF_v1/fig/layout/ftf_layout_chip.png'), ports=pcb_ports)
    plot_gds(gds_file, _S('Design1_FTF_v1/fig/layout/ftf_layout_zoom.png'), window=(-1.0, 1.0, Y0 - 0.55, Y0 + 0.55))
    plot_gds(gds_file, _S('Design1_FTF_v1/fig/layout/ftf_layout_fx1.png'), window=(-0.95, -0.3, Y0 - 0.3, Y0 + 0.3))
    # FX1 JJ + array site (empty in the GDS) and the flux-line tip
    plot_gds(gds_file, _S('Design1_FTF_v1/fig/layout/ftf_layout_fxjj.png'), window=(-0.48, -0.32, Y0 - 0.045, Y0 + 0.045), dpi=140)
    print('wrote', gds_file, fab_file)
