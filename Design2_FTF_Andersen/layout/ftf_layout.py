"""Design2_FTF_Andersen layout build and export entry points.

build_design() loads the full layout3 design.
build_notebook_design() reproduces the previous CPW 10/6 notebook model.
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

from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight

HERE = _S('Design2_FTF_Andersen/layout')
sys.path.insert(0, HERE)
sys.path[1:1] = _pattern_paths()   # shared QComponents

sys.path.insert(0, os.path.join(HERE, 'components'))   # Design2-specific components
from fluxonium_pocket import FluxoniumPocket  # noqa: E402  (lab FluxoniumPocket, Design2 only)
from tmon_delft import TmonDelft  # noqa: E402
from tmon_squid_fab import TmonSQUID, TmonSQUIDFab  # noqa: E402
from xmark import add_corner_marks  # noqa: E402
from bondpad import add_bond_pad
from gds_tools import export_gds, plot_gds, junction_window  # noqa: E402

GDS_DIR = _S('Design2_FTF_Andersen/gds')
FIG_DIR = _S('Design2_FTF_Andersen/fig/layout')   # GDS layout plots (same set as Design1_FTF_v1 / 3)

CPW_WIDTH, CPW_GAP = '10 um', '6 um'
CPW_Z0 = 50.9           # Ohm, 10 / 6 um CPW on Si (as Design1_FTF_v1)

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
RES_X = -0.032          # mm, x of the resonator open end (OpenToGround) next to the feedline
# per resonator: design frequency, kappa target, OpenToGround y (notebook), parallel length (kappa)
# and total length (f_r); start_straight / length are tuned with Q3D (kappa) and HFSS eigenmode (f_r)
RESONATORS = OrderedDict([
    ('Q1', dict(f=6.3, kappa=2.0, otg_y=-2.12, start_straight=1.60, length=9.25, jog='R')),   # Q3D d2_kappa_ss160_155: 2.00 MHz
    ('Q2', dict(f=6.5, kappa=2.0, otg_y=-0.40, start_straight=1.52, length=8.97, jog='R')),   # 1.55 mm: 2.07 MHz -> 1.52
])
FILLET = '99um'

# die marks (as Design1_FTF_v1 / 3): metal crosses in etched square windows, MARK_INSET mm in from the chip corners, die ID at the
# lower-right of the UL cross; build_design(die_id='FTF07'), FTF00 ... FTF99
MARK_INSET = 0.35
MARK = dict(arm_length='250um', arm_width='10um', window_size='400um', label_size='150um', label_gap='30um')   # LABEL_SIZE_UM in xmark.py
DIE_ID = 'FTF00'



def build_design(source=None, *, readout_ports=False):
    """Build the restored full-chip Design2_FTF_Andersen geometry from the reference mask."""
    from pathlib import Path
    from layout3_full import build_design as build_layout3

    reference = _P('Design2_FTF_Andersen/gds/layout3 ftf.gds')
    return build_layout3(source if source is not None else reference, readout_ports)


def build_notebook_design(fx1_options=None, fx2_options=None, tmon_options=None, marks=True, die_id=DIE_ID,
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
    tmon = TmonDelft
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
    add_bond_pad(design, 'LPReadIn', dict(pos_x='0mm', pos_y='4.0mm', orientation='270'), inset=0)
    add_bond_pad(design, 'LPReadOut', dict(pos_x='0mm', pos_y='-4.0mm', orientation='90'), inset=0)
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
    from layout3_full import main

    main()
