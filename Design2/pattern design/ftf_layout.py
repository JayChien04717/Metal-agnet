"""Design2: one FTF set (Q1 - Q3 - Q2) of the lab design FXTFX24-W7-0x-V.3.1, with the CPW
set to 10 / 6 um (the default of the components in this folder).

Same positions, orientations and options as cells 7-8 of
2FX_with_Tmon/FXTFX24-W7-0x-V.3.1-ALL_Sims.ipynb, but without the launch pads, routes and the
other three sets: for the Q3D capacitance matrix only the three qubits are rendered (as in
that notebook, cell 33). The components are the ones in this folder (fluxonium_withJJ.py,
tmon_c.py), not the copies in qiskit_metal.

    from ftf_layout import build_design
    design = build_design()
"""
import os
import sys

from qiskit_metal import designs, Dict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'pattern design'))   # shared QComponents

from fluxonium_withJJ import FluxoniumPocket  # noqa: E402
from tmon_c import Tmon_c  # noqa: E402

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


def build_design(fx1_options=None, fx2_options=None, tmon_options=None):
    design = designs.DesignPlanar()
    design.chips.main.material = 'silicon'
    design.variables['cpw_width'] = CPW_WIDTH
    design.variables['cpw_gap'] = CPW_GAP
    design._chips['main']['size']['size_x'] = '10mm'
    design._chips['main']['size']['size_y'] = '10mm'
    design._chips['main']['size']['size_z'] = '-525um'
    design.overwrite_enabled = True
    FluxoniumPocket(design, 'Q1', options=dict(**Q1, **dict(FX1_OPTIONS, **(fx1_options or {}))))
    FluxoniumPocket(design, 'Q2', options=dict(**Q2, **dict(FX2_OPTIONS, **(fx2_options or {}))))
    Tmon_c(design, 'Q3', options=dict(**Q3, **dict(TMON_OPTIONS, **(tmon_options or {}))))
    return design


def plot_design(design, path, xlim=None, ylim=None, title=''):
    """Metal (poly) filled, ground cut-outs (subtract) hatched, from the qgeometry tables."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 10))
    polys = design.qgeometry.tables['poly']
    for _, row in polys.iterrows():
        geoms = getattr(row.geometry, 'geoms', [row.geometry])
        for g in geoms:
            x, y = g.exterior.xy
            if row['subtract']:
                ax.fill(x, y, fc='none', ec='tab:red', lw=0.6, hatch='///', alpha=0.5)
            else:
                ax.fill(x, y, fc='#9aa3af', ec='k', lw=0.4)
    ax.set_aspect('equal')
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title(title)
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    d = build_design()
    out = os.path.join(os.path.dirname(HERE), 'fig', 'layout')
    os.makedirs(out, exist_ok=True)
    plot_design(d, os.path.join(out, 'ftf_layout_set.png'), title='Design2: Q1 - Q3 - Q2 (CPW 10 / 6 um)')
    for n in ('Q1', 'Q2', 'Q3'):
        print(n, [k for k in d.components[n].pins])
