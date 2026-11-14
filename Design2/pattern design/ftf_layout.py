"""Design2: one FTF set (Q1 - Q3 - Q2) of the lab design FXTFX24-W7-0x-V.3.1, with the CPW
set to 10 / 6 um (the default of the components in this folder).

Same positions, orientations and options as cells 7-8 of
2FX_with_Tmon/FXTFX24-W7-0x-V.3.1-ALL_Sims.ipynb, and the readout in the notebook's style
(cells 4, 11): a straight feedline down the chip centre (LPReadIn (0, 4) -> LPReadOut (0, -4)) and
lambda/2 resonators that start at an OpenToGround 32 um left of the feedline, run parallel to it
for `start_straight` (sets kappa) and meander to the qubit's readout pad. Only Q1 / Q2 (6.3 / 6.5 GHz)
of the lab's eight resonators; no control-line routes / other sets. The components are the ones in this folder (fluxonium_withJJ.py,
tmon_c.py), not the copies in qiskit_metal.

    from ftf_layout import build_design
    design = build_design()                 # junction-free (default): no JJ / SQUID pads between the stems
    design = build_design(jj_metal=True)    # old version with the lab JJ / SQUID pads
    design = build_design(for_process=True) # simulation model (Q3D / HFSS): array = U-shaped hairpin of
                                            # two lumped L + metal end, JJ = short line (as Design1 / 3)
"""
import os
import sys

from collections import OrderedDict

from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'pattern design'))   # shared QComponents

from fluxonium_withJJ import FluxoniumPocket  # noqa: E402
from tmon_c import Tmon_c  # noqa: E402
from tmon_squid_fab import TmonSQUID, TmonSQUIDFab  # noqa: E402
from xmark import add_corner_marks  # noqa: E402
from gds_tools import export_gds, plot_gds, junction_window  # noqa: E402

GDS_DIR = os.path.join(os.path.dirname(HERE), 'gds')
FIG_DIR = os.path.join(os.path.dirname(HERE), 'fig', 'layout')   # GDS layout plots (same set as Design1 / 3)

CPW_WIDTH, CPW_GAP = '10 um', '6 um'
CPW_Z0 = 50.9           # Ohm, 10 / 6 um CPW on Si (as Design1)

# positions / orientations of the FXTFX24 set (notebook cells 7, 8)
Q1 = dict(pos_x='-2.3mm', pos_y='-2.270mm', orientation=270)
Q2 = dict(pos_x='-2.3mm', pos_y='-1.30mm', orientation=90)
Q3 = dict(pos_x='-2.46mm', pos_y='-1.785mm', orientation='90')

FX_COMMON = dict(chip='main', nanowire=False, jj_length='10um', inductor_orientation='-1',
                 flux_bias_line_options=Dict(make_fbl=True))
# per qubit: the lab design mirrors the readout / charge line side and the JJ orientation
FX1_OPTIONS = dict(FX_COMMON, charge_line_options=Dict(loc_H=-1, make_cl=True),
                   readout_line_options=Dict(loc_H=+1, make_rol=True))
FX2_OPTIONS = dict(FX_COMMON, jj_orientation='-1',
                   charge_line_options=Dict(loc_H=+1, make_cl=True),
                   readout_line_options=Dict(loc_H=-1, make_rol=True))
TMON_OPTIONS = dict(jj_gap='30um', pad_arm_length='410um', pad_head_length='400um',
                    pad_arm_width='25um', palm_radius='63um', pad_gap='30um',
                    circle_gap='8um', x_gap='38.5um',
                    flux_bias_line_options=Dict(make_fbl=True, fbl_sep='10um', fbl_length='100um'))
# readout (notebook cells 4, 11): feedline on x = 0, resonators start 32 um left of it
LP = dict(chip='main', lead_length='30um', pad_width='200um', pad_height='200um', pad_gap='76.5um')
RES_X = -0.032          # mm, x of the resonator open end (OpenToGround) next to the feedline
# per resonator: design frequency, kappa target, OpenToGround y (notebook), parallel length (kappa)
# and total length (f_r); start_straight / length are tuned with Q3D (kappa) and HFSS eigenmode (f_r)
RESONATORS = OrderedDict([
    ('Q1', dict(f=6.3, kappa=2.0, otg_y=-2.12, start_straight=1.60, length=9.25, jog='R')),   # Q3D d2_kappa_ss160_155: 2.00 MHz
    ('Q2', dict(f=6.5, kappa=2.0, otg_y=-0.40, start_straight=1.52, length=8.97, jog='R')),   # 1.55 mm: 2.07 MHz -> 1.52
])
FILLET = '99um'

# die marks (as Design1 / 3): metal crosses in etched square windows, MARK_INSET mm in from the chip corners, die ID at the
# lower-right of the UL cross; build_design(die_id='FTF07'), FTF00 ... FTF99
MARK_INSET = 0.35
MARK = dict(arm_length='250um', arm_width='10um', window_size='400um', label_size='100um', label_gap='30um')
DIE_ID = 'FTF00'


def build_design(fx1_options=None, fx2_options=None, tmon_options=None, marks=True, die_id=DIE_ID,
                 jj_metal=False, for_process=False, readout=True, resonators=None,
                 include_tc_landing_pads=False, fab_contacts=False):
    """Default jj_metal=False: no JJ pads between the Q1 / Q2 stems and no SQUID pads on Q3
    (used by the existing simulation). Export enables include_tc_landing_pads, and fab_contacts adds the Dolan recesses. jj_metal=True restores the lab pads."""
    jj = {} if jj_metal else dict(jj_metal='False')
    fx_sim = dict(u_array='True') if for_process else {}
    design = designs.DesignPlanar()
    design.chips.main.material = 'silicon'
    design.variables['cpw_width'] = CPW_WIDTH
    design.variables['cpw_gap'] = CPW_GAP
    design._chips['main']['size']['size_x'] = '10mm'
    design._chips['main']['size']['size_y'] = '10mm'
    design._chips['main']['size']['size_z'] = '-525um'
    design.overwrite_enabled = True
    FluxoniumPocket(design, 'Q1', options=dict(**Q1, **dict(FX1_OPTIONS, **jj, **fx_sim, **(fx1_options or {}))))
    FluxoniumPocket(design, 'Q2', options=dict(**Q2, **dict(FX2_OPTIONS, **jj, **fx_sim, **(fx2_options or {}))))
    if fab_contacts and (jj_metal or for_process or not include_tc_landing_pads):
        raise ValueError("Fab recesses require a junction-free mask with TC landing pads.")
    tmon = Tmon_c
    options = dict(TMON_OPTIONS, **jj, **(tmon_options or {}))
    if include_tc_landing_pads:
        tmon = TmonSQUIDFab if fab_contacts else TmonSQUID
        options['jj_metal'] = 'False'
    tmon(design, 'Q3', options=dict(Q3, **options))
    if readout:
        add_readout(design, resonators)
    if marks:
        add_corner_marks(design, die_id=die_id, inset=MARK_INSET, **MARK)
    return design


def add_readout(design, resonators=None):
    """Feedline LPReadIn -> LPReadOut on x = 0 and one lambda/2 resonator per qubit (res_<Q>),
    notebook style. resonators: {name: overrides} merged into RESONATORS."""
    LaunchpadWirebond(design, 'LPReadIn', options=dict(LP, pos_x='0mm', pos_y='4.0mm', orientation='270'))
    LaunchpadWirebond(design, 'LPReadOut', options=dict(LP, pos_x='0mm', pos_y='-4.0mm', orientation='90'))
    RouteStraight(design, 'TL', options=Dict(pin_inputs=Dict(
        start_pin=Dict(component='LPReadIn', pin='tie'), end_pin=Dict(component='LPReadOut', pin='tie'))))
    for q, r in RESONATORS.items():
        r = dict(r, **(resonators or {}).get(q, {}))
        OpenToGround(design, f'otg_{q}', options=dict(chip='main', pos_x=f'{RES_X}mm', pos_y=f"{r['otg_y']}mm",
                                                      orientation='90'))
        jogs = OrderedDict([(0, [r['jog'], '200um']), (1, [r['jog'], '200um'])])
        RouteMeander(design, f'res_{q}', Dict(
            meander=Dict(spacing='200um'), total_length=f"{r['length']}mm", fillet=FILLET,
            lead=Dict(start_straight=f"{r['start_straight']}mm", start_jogged_extension=jogs),
            pin_inputs=Dict(start_pin=Dict(component=f'otg_{q}', pin='open'),
                            end_pin=Dict(component=q, pin='readout_line'))))


if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    gds_file = os.path.join(GDS_DIR, 'ftf_layout.gds')
    fab_file = os.path.join(GDS_DIR, 'ftf_layout_fab.gds')
    d = build_design(include_tc_landing_pads=True)
    for n in ('Q1', 'Q2', 'Q3'):
        print(n, [k for k in d.components[n].pins])
    fxjj_win = junction_window(d, 'Q1', half=0.12, half_y=0.045)
    tc_window = junction_window(d, 'Q3', half=0.055)
    export_gds(d, gds_file)
    export_gds(build_design(include_tc_landing_pads=True, fab_contacts=True), fab_file, fab=True)
    plot_gds(gds_file, os.path.join(FIG_DIR, 'ftf_layout_chip.png'))
    plot_gds(gds_file, os.path.join(FIG_DIR, 'ftf_layout_zoom.png'), window=(-3.4, 0.4, -3.0, 0.2))
    plot_gds(gds_file, os.path.join(FIG_DIR, 'ftf_layout_fx1.png'), window=(-2.85, -1.75, -2.82, -1.72))
    # Q1 junction site (empty in the GDS) between the two stems
    plot_gds(gds_file, os.path.join(FIG_DIR, 'ftf_layout_fxjj.png'), window=fxjj_win, dpi=140)
    plot_gds(fab_file, os.path.join(FIG_DIR, 'ftf_layout_tc_fab.png'), window=tc_window, dpi=160)
    print('wrote', gds_file, fab_file)
