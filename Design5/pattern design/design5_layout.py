"""Design5 (= Design6 variant W1, Wallraff-group parameters, thesis Table 2.4 Q35; see SQL HFSS/note/wallraff reference).
Default build = W1; D5_VARIANT=P1 gives the old Design6 c6_P1 geometry.

Original Design6 description: one readout channel in the style of Cai et al. 2024, Fig. 3(e) (Heinsoo-type individual Purcell filter).

    P1 --C_in-- node N0 ------------------------------- P2          (feedline, CPW 10/6 um)
                  |
                 C_p  (T-type enclosing capacitor)
                  |
                node J1 --- Purcell filter (lambda/4, open end at J1, far end shorted)
                  |
                 C_j  (small interdigital, ~2 fF)
                  |
                readout resonator (lambda/4, open end at C_j, far end shorted)

Targets (paper example, scaled to 6.7 GHz): C_in = C_p = 40 fF, C_j = 2 fF, f_filter = f_res.
Usage: python design5_layout.py
"""
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                                   # .../Design6
HFSS = os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(HFSS, 'Design4', 'pattern design'))
import purcell_layout as pl  # noqa: E402
sys.path.insert(1, os.path.join(HFSS, 'Design3', 'pattern design'))
import ftf_layout as f3  # noqa: E402
from qiskit_metal import designs, Dict  # noqa: E402
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond  # noqa: E402
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander  # noqa: E402
from qiskit_metal.qlibrary.tlines.mixed_path import RouteMixed  # noqa: E402
from collections import OrderedDict  # noqa: E402
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight  # noqa: E402
from interdigital_cap import InterdigitalCap  # noqa: E402  (C_in/C_p no longer use it; kept for idc helper)
from res_end import ResEnd  # noqa: E402
from tap_junction import TapJunction  # noqa: E402
from pcbport import pcb_ports  # noqa: E402
sys.path.insert(0, HERE)
from enclosing_cap import EnclosingCap  # noqa: E402
from hammer_cap import HammerCap  # noqa: E402

GDS_FILE = os.path.join(ROOT, 'gds', 'design5_P1.gds')
V, Z0 = pl.CPW_V, pl.CPW_Z0

F0 = 6.7                              # GHz (filter and resonator)
C_IN_FF, C_P_FF, C_J_FF = 40.0, 40.0, 2.0
CIN_N, CIN_RATIO = 6, 0.90            # input finger cap: fingers, Q3D / Alley guess
FINGER = dict(w=6.0, g=6.0, eg=6.0, bw=20.0, taper=100.0, lead=40.0)           # um
CIN_OVERLAP = round(pl.overlap_for(C_IN_FF / CIN_RATIO, CIN_N), 1)
CJ_HW, CJ_HH, CJ_HG, CJ_LEAD, CJ_GG = 25.0, 65.0, 26.0, 10.0, 16.0   # c6_P1 (paper values): gap 26 um = Q3D 2.00 fF (c6_T5 had 35 um, 1.48 fF);   # C_j hammer (um), read off paper Fig. 3e inset (2 fF there): two parallel plates 25 x 65, gap 15, pocket margin 16
CIN_TYPE, CP_TYPE = 'T', 'T'           # enclosing capacitor type (Cai 2024): T = opening toward the filter side, C = opening toward the transmission line
# enclosing capacitor read off the paper (Fig. 3d, 100 um scale bar): thick solid frame 324 x 172 um with a slot,
# a 193 x 41 um bar with a stem inside it, gap ~8 um, frame ~12 um from the ground plane
ENC = dict(plate_width=140.0, plate_length=41.0, gap=8.0, wall=52.0, wall_side=57.0, leg_ext=72.0, lead_length=40.0, ground_gap=12.0)
# C_p plate: c6_P1 (paper values C_in = C_p = 40 fF): 140 um = Q3D 40.6 fF (c6_T5 used 175 um / 47.3 fF from a wrong fit factor)
ENC_P = dict(ENC, plate_width=140.0)
LAMBDA4 = V / (4 * F0 * 1e9) * 1e3    # mm
LOAD = lambda c_ff: Z0 * 2 * np.pi * F0 * 1e9 * c_ff * 1e-15 * V / (2 * np.pi * F0 * 1e9) * 1e3   # open-end shortening, mm
L_FILTER = LAMBDA4 - LOAD(C_P_FF + C_J_FF) - 0.1533  # c6_P1: raw Q3D model, S21 zeros centred on 6.70 GHz (T5 -0.2975); open end -> short; -0.125 mm: HFSS c6_v3 showed the filter 4.5 % low (Cp + 16 fF to ground load the open end more than assumed)
L_RES = LAMBDA4 - LOAD(C_J_FF)                # from the C_j side to the short
Y0 = 3.0                               # feedline height
D_CJ = 1.02                            # C_j junction distance from the open ends (paper: ~1 mm)
VERT, OPEN_DY = 1.63, 0.44             # straight length of the vertical lines, resonator open end offset (mm)
XM, YM = 0.87, 1.2                     # horizontal run before the snake, snake height (mm); c6_T5: 0.9 -> 0.87 keeps legs >= 350 um
XMF, YMF = 0.8, 0.8                    # filter: run before the snake, snake height (mm)
MSP_F, MSP_R = 200, 350                      # meander spacing (um) of filter / resonator
CJ_TAPER = 0.0                         # um, C_j has no taper (compact, as in the paper)
CJ_TAPER_X = CJ_HG / 2 + CJ_HW + CJ_LEAD                 # um, centre -> pin of C_j
CJ_GAP = 0.02                          # mm, gap between a stub end and a C_j pin
L_RES_TOTAL = LAMBDA4 - 0.0807                 # open end -> shorted end; c6_T3 -0.085; c6_T5 -0.329;
                                               # c6_P1: raw Q3D model, aligned with the filter -> zeros 6.685 / 6.715 GHz
L_RES_TOTAL += float(os.environ.get('DLR_EXTRA', 0))   # diagnostic / calibration: extra resonator length (mm), default 0
L_FILTER += float(os.environ.get('DLF_EXTRA', 0))      # diagnostic / calibration: extra filter length (mm), default 0
CIN_X = -0.3025                               # c6_T5: C_in next to n0 (C_in -> n0 line 1.25 -> 0.25 mm, symmetric peak); was -1.3

# ---- variant W1: Wallraff-group parameters (thesis Table 2.4, Q35; note/wallraff reference, simulation/wallraff_model.py) ----
# D5_VARIANT=W1 (default):  hugger flipped (frame = big ground C on the FILTER side, C_pfg ~ 93 fF), C_pf ~ 45.5 fF, hammer ~ 2.8 fF,
#                 filter tap L_pf = 1.2 mm from C_p, resonator tap L_rr = 1.3 mm from the open end,
#                 lengths from the W1 circuit model (wr 6.700, wp 6.719 GHz; J 16.9, kp 37.6 MHz) -> to be re-solved with Q3D C
VARIANT = os.environ.get('D5_VARIANT', 'W1')
if VARIANT == 'W1':
    CP_TYPE = 'C'                                   # 'C': plate (pin b) toward n0, frame (pin a) toward the filter
    ENC_P = dict(ENC, plate_width=165.0)            # interpolated: 140 um -> 40.6 fF, 175 um -> 47.3 fF  => 45.5 fF (Q3D to confirm)
    CJ_HG = 19.0                                    # C = 52.6/gap - 0.02 from Q3D (26 um 2.00 fF, 35 um 1.48 fF) => 2.77 fF
    CJ_TAPER_X = CJ_HG / 2 + CJ_HW + CJ_LEAD
    D_CJ = 1.20                                     # L_pf
    OPEN_DY = 0.87                                  # resonator tap from the open end = VERT + OPEN_DY - D_CJ = 1.30 mm (L_rr)
    XMF, YMF = 0.5, 0.5                             # filter run / snake height (shorter filter, 3.63 mm; scanned: routes cleanly)
    # lengths: Wallraff-cap model gave 3.6289 / 4.4578; re-solved with the W1 Q3D pi networks (solve_w1.py:
    # C_p 45.41 fF, frame 119 fF to ground on the filter side; C_j 2.59 fF, 13 fF to ground) -> wr 6.700, wp 6.719 GHz,
    # J 14.5, kp 36.6 MHz, S21 zeros 6.6917 / 6.7256 GHz
    L_FILTER = 3.4542 + float(os.environ.get('DLF_EXTRA', 0))
    L_RES_TOTAL = 4.4070 + float(os.environ.get('DLR_EXTRA', 0))
    GDS_FILE = os.path.join(ROOT, 'gds', 'design5.gds')


def finger_xpin(overlap_um):
    return (overlap_um + 2 * FINGER['eg']) / 2 + FINGER['bw'] + FINGER['taper'] + FINGER['lead']


def enc_xpin():
    return ENC['plate_length'] / 2 + ENC['leg_ext'] + ENC['lead_length']


def build_design(filter_only=False):
    warnings.filterwarnings('ignore')
    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon'
    design.variables['cpw_width'] = f3.CPW_WIDTH
    design.variables['cpw_gap'] = f3.CPW_GAP
    design.chips.main.size.size_x = '10mm'
    design.chips.main.size.size_y = '10mm'
    design.chips.main.size.size_z = '-525um'
    for name in ('P8', 'P7'):
        LaunchpadWirebond(design, 'LP' + name, options=dict(f3.LP, **pcb_ports[name]))
    cpw = dict(trace_width=f3.CPW_WIDTH, trace_gap=f3.CPW_GAP)

    def idc(name, n, overlap, x, y, orient, taper=None, fw=None, bus=None):
        InterdigitalCap(design, name, options=dict(
            chip='main', pos_x=f'{x:.5f}mm', pos_y=f'{y:.5f}mm', orientation=str(orient), finger_count=str(n),
            finger_width=f"{FINGER['w'] if fw is None else fw}um", finger_gap=f"{FINGER['g']}um", overlap=f'{overlap}um',
            end_gap=f"{FINGER['eg']}um", bus_width=f"{FINGER['bw'] if bus is None else bus}um", taper_length=f"{FINGER['taper'] if taper is None else taper}um",
            lead_length=f"{FINGER['lead']}um", ground_gap='6um', **cpw))

    def straight(name, a, b):
        RouteStraight(design, name, Dict(chip='main', pin_inputs=Dict(
            start_pin=Dict(component=a[0], pin=a[1]), end_pin=Dict(component=b[0], pin=b[1]))))

    # ---- feedline: P1 -> C_in -> node N0 -> P2 ----
    # C_in: enclosing capacitor in line with the feedline (as in the paper, Fig. 3e): frame on the P1 side, plate toward the node
    EnclosingCap(design, 'cin', options=dict(
        chip='main', pos_x=f'{CIN_X}mm', pos_y=f'{Y0}mm', orientation='0' if CIN_TYPE == 'T' else '180',
        plate_width=f"{ENC['plate_width']}um", plate_length=f"{ENC['plate_length']}um", gap=f"{ENC['gap']}um",
        wall=f"{ENC['wall']}um", wall_side=f"{ENC['wall_side']}um", ground_gap=f"{ENC['ground_gap']}um", leg_ext=f"{ENC['leg_ext']}um", lead_length=f"{ENC['lead_length']}um", **cpw))
    f3.route(design, 'feed_in', ('LPP8', 'tie'), ('cin', 'a' if CIN_TYPE == 'T' else 'b'), [(-1.1, 3.9), (CIN_X - 0.6, 3.9), (CIN_X - 0.6, Y0)] if CIN_X < -1.0 else [(-1.1, Y0)])
    TapJunction(design, 'n0', options=dict(pos_x='0mm', pos_y=f'{Y0}mm', orientation='0', prime_length='0.3mm',
                                           stub_length='0.1mm', stub_up='False', **cpw))
    straight('feed_c', ('cin', 'b' if CIN_TYPE == 'T' else 'a'), ('n0', 'prime_start'))
    f3.route(design, 'out_line', ('n0', 'prime_end'), ('LPP7', 'tie'), [(1.1, Y0)])

    # ---- C_p (T-type enclosing capacitor), vertical, pin a on top ----
    y = Y0 - 0.1 - 0.05 - enc_xpin() / 1000
    EnclosingCap(design, 'cp', options=dict(
        chip='main', pos_x='0mm', pos_y=f'{y:.5f}mm', orientation='-90' if CP_TYPE == 'T' else '90',
        plate_width=f"{ENC_P['plate_width']}um", plate_length=f"{ENC['plate_length']}um", gap=f"{ENC['gap']}um",
        wall=f"{ENC['wall']}um", wall_side=f"{ENC['wall_side']}um", ground_gap=f"{ENC['ground_gap']}um", leg_ext=f"{ENC['leg_ext']}um", lead_length=f"{ENC['lead_length']}um", **cpw))
    CP_FPIN = 'b' if CP_TYPE == 'T' else 'a'                   # C_p pin on the filter side
    straight('stub_p', ('n0', 'second_end'), ('cp', 'a' if CP_TYPE == 'T' else 'b'))
    # ---- paper layout (Cai et al., Fig. 3e): filter and resonator run side by side, anti-parallel;
    #      C_j bridges the two straight segments ~1 mm from each open end ----
    y_c = y - enc_xpin() / 1000                                  # C_p pin b = open end of the filter
    y_j = y_c - D_CJ                                             # height of the C_j junctions
    DX = 0.20                                                    # spacing of the two vertical lines
    prime, half = 0.3, 0.15
    if filter_only:
        # Cin + Cp + filter only (no resonator, no C_j): filter line from C_p pin b (open end) to the short
        y_f = y_c - VERT
        ResEnd(design, 'end_filter', options=dict(pos_x=f'{-XMF:.4f}mm', pos_y=f'{y_f + YMF:.4f}mm', orientation='270'))
        rf = Dict(chip='main', total_length=f'{L_FILTER:.4f}mm', fillet=f3.FILLET,
                  hfss_wire_bonds=False, meander=Dict(spacing=f'{MSP_F}um', asymmetry='0um'),
                  lead=Dict(start_straight='100um', end_straight='100um'),
                  anchors=OrderedDict({0: np.array([-XMF, y_f]), 1: np.array([0.0, y_f])}),
                  between_anchors=OrderedDict({0: 'M', 1: 'S', 2: 'S'}),
                  pin_inputs=Dict(start_pin=Dict(component='end_filter', pin='end'),
                                  end_pin=Dict(component='cp', pin=CP_FPIN)))
        RouteMixed(design, 'filt', rf)
        return design
    # filter junction (vertical line, stub to the right)
    stub_len = (DX - 2 * (CJ_TAPER_X) / 1000 - 2 * CJ_GAP) / 2
    TapJunction(design, 'jf', options=dict(pos_x='0mm', pos_y=f'{y_j:.5f}mm', orientation='90', prime_length=f'{prime}mm',
                                           stub_length=f'{stub_len:.5f}mm', stub_up='False', **cpw))
    straight('feed_jf', ('cp', CP_FPIN), ('jf', 'prime_end'))
    # resonator junction (vertical line, stub to the left)
    TapJunction(design, 'jr', options=dict(pos_x=f'{DX}mm', pos_y=f'{y_j:.5f}mm', orientation='90', prime_length=f'{prime}mm',
                                           stub_length=f'{stub_len:.5f}mm', stub_up='True', **cpw))
    # C_j between the two stubs
    HammerCap(design, 'cj', options=dict(chip='main', pos_x=f'{DX / 2:.5f}mm', pos_y=f'{y_j:.5f}mm', orientation='0',
                                         head_width=f'{CJ_HW}um', head_height=f'{CJ_HH}um', head_gap=f'{CJ_HG}um',
                                         lead_length=f'{CJ_LEAD}um', ground_gap=f'{CJ_GG}um', **cpw))
    straight('stub_f', ('jf', 'second_end'), ('cj', 'a'))
    straight('stub_r', ('jr', 'second_end'), ('cj', 'b'))
    # filter: down along x = 0 to y_f, left, then a snake going up, shorted far end
    y_f = y_c - VERT
    ResEnd(design, 'end_filter', options=dict(pos_x=f'{-XMF:.4f}mm', pos_y=f'{y_f + YMF:.4f}mm', orientation='270'))
    rf = Dict(chip='main', total_length=f'{L_FILTER - (y_c - (y_j + half)) - prime:.4f}mm', fillet=f3.FILLET,
              hfss_wire_bonds=False, meander=Dict(spacing=f'{MSP_F}um', asymmetry='0um'),
              lead=Dict(start_straight='100um', end_straight='100um'),
              anchors=OrderedDict({0: np.array([-XMF, y_f]), 1: np.array([0.0, y_f])}),
              between_anchors=OrderedDict({0: 'M', 1: 'S', 2: 'S'}),
              pin_inputs=Dict(start_pin=Dict(component='end_filter', pin='end'),
                              end_pin=Dict(component='jf', pin='prime_start')))
    RouteMixed(design, 'filt', rf)
    # resonator: open end at the bottom, up along x = DX, right, then a snake going down, shorted far end
    y_ro = y_c - VERT - OPEN_DY
    ResEnd(design, 'open_res', options=dict(pos_x=f'{DX}mm', pos_y=f'{y_ro:.5f}mm', orientation='90'))
    straight('res_a', ('open_res', 'end'), ('jr', 'prime_start'))
    y_top = y_c - 0.44
    ResEnd(design, 'end_res', options=dict(pos_x=f'{DX + XM}mm', pos_y=f'{y_top - YM:.4f}mm', orientation='90'))
    rr = Dict(chip='main', total_length=f'{L_RES_TOTAL - (y_j - half - y_ro) - prime:.4f}mm', fillet=f3.FILLET,
              hfss_wire_bonds=False, meander=Dict(spacing=f'{MSP_R}um', asymmetry='0um'),
              lead=Dict(start_straight='100um', end_straight='100um'),
              anchors=OrderedDict({0: np.array([DX + XM, y_top]), 1: np.array([DX, y_top])}),
              between_anchors=OrderedDict({0: 'M', 1: 'S', 2: 'S'}),
              pin_inputs=Dict(start_pin=Dict(component='end_res', pin='end'),
                              end_pin=Dict(component='jr', pin='prime_end')))
    RouteMixed(design, 'res', rr)
    return design


if __name__ == '__main__':
    import geopandas as gpd
    if not hasattr(gpd.GeoDataFrame, 'append'):          # qiskit-metal's gds renderer still calls DataFrame.append
        gpd.GeoDataFrame.append = gpd.GeoDataFrame._append
    d = build_design()
    print(f'C_in overlap {CIN_OVERLAP} um ({C_IN_FF} fF); C_in plate {ENC["plate_width"]} um, C_p plate {ENC_P["plate_width"]}x{ENC_P["plate_length"]} um; '
          f'C_j hammer {CJ_HW}x{CJ_HH} um, gap {CJ_HG} um; filter {L_FILTER:.3f} mm, resonator {L_RES_TOTAL:.3f} mm (lambda/4 {LAMBDA4:.3f})')
    os.makedirs(os.path.dirname(GDS_FILE), exist_ok=True)
    pl.export_gds(d, GDS_FILE)
    png = os.path.join(ROOT, 'fig', 'design5_layout.png' if VARIANT == 'W1' else f'design5_{VARIANT}_layout.png')
    pl.plot_gds(png.replace('.png', '.gds') if False else GDS_FILE, png, window=(-3, 3, -0.5, 5))
    print('gds', GDS_FILE, 'png', png)
