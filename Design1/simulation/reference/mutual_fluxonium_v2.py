"""Flux-line mutual inductance of the reference FluxoniumPocket (comparison).

Loads FluxoniumPocket from the repo file qiskit_metal/qlibrary/qubits/fluxonium_v2.py
(the installed qmetal copy differs) with its default options + make_fbl=True, and runs
the same HFSS driven-modal extraction as ../hfss_mutual_fx.py:
Port 1 = flux_bias_line pin, Port 2 = rect_jj (port only), the inductor keeps
hfss_inductance (200 nH).

Two things in fluxonium_v2.py as it stands:
* make_pocket reads p.nanowire, which is not an option (the option is
  nanowire_inductor), so the JJ-array branch is always drawn: inductor 130 um long at
  x = l_arm_length - l_arm_width = 18 um.
* the flux_bias_line pin is drawn from y = 0 to -(fbl_height + cpw_width), so its
  width is 60 um instead of cpw_width; an HFSS port on it would touch the ground next
  to the trace. FluxoniumPocketSim re-adds that pin 10 um wide, with gap = cpw_gap,
  at the end of the bottom leg. The geometry is unchanged.

Usage:  python reference/mutual_fluxonium_v2.py [tag]
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'driven'))
from qiskit_metal import designs, draw, Dict  # noqa: E402
from hfss_mutual_fx import run_design  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
V2_FILE = os.path.join(REPO, 'qiskit_metal', 'qlibrary', 'qubits', 'fluxonium_v2.py')
_spec = importlib.util.spec_from_file_location('fluxonium_v2_repo', V2_FILE)
_v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_v2)


class FluxoniumPocketSim(_v2.FluxoniumPocket):
    """FluxoniumPocket with a cpw_width-wide flux_bias_line pin (for an HFSS port)."""

    def make_flux_bias_line(self):
        super().make_flux_bias_line()
        p = self.p
        pfb = p.flux_bias_line_options
        w = pfb.cpw_width
        x_end = p.pocket_width / 2 + pfb.fbl_sep / 2      # end of the bottom leg
        y_leg = -pfb.fbl_height / 2 - w / 2                # centre of the bottom leg
        line = draw.LineString([(x_end, y_leg + w / 2), (x_end, y_leg - w / 2)])
        line = draw.translate(draw.rotate(line, p.orientation, origin=(0, 0)), p.pos_x, p.pos_y)
        self.add_pin('flux_bias_line', list(line.coords), w, gap=pfb.cpw_gap)


def build_v2_design(options=None):
    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = 'silicon'
    design.chips.main.size.size_z = '-525um'
    design.variables['cpw_width'] = '10um'
    design.variables['cpw_gap'] = '6um'
    opt = dict(flux_bias_line_options=Dict(make_fbl=True))
    opt.update(options or {})
    FluxoniumPocketSim(design, 'FXv2', options=opt)
    return design


if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    tag = (sys.argv[1:2] or ['mutual_fluxonium_v2'])[0]
    if tag == 'check':
        d = build_v2_design()
        for k, v in d.components['FXv2'].pins.items():
            print(k, v.middle, v.normal, v.width, v.gap)
    else:
        run_design(tag, build_v2_design(), 'FXv2', 'flux_bias_line', 'rect_jj',
                   info=dict(reference=V2_FILE))
