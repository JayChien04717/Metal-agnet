# -*- coding: utf-8 -*-
r"""Grounded tunable transmon coupler whose two ends are U-shaped claws (Dolan bridge SQUID).

Subclasses ``TmonCouplerDolanNthu`` and overrides only ``_end_geometry()``: stem,
SQUID, flux bias line, charge line, readout pad, and 4 um DRC checks are all
inherited from the base class.

Sketch (orientation = 0, stem points to -y)::

        ____                                    ____
       |  __|__________________________________|__  |     claws (open toward the qubits)
  Q1 ->|_|  |__________________________________|  |_|<- Q2       bar (pad_arm_width)
       |____|               |   |              |____|
                            |   |  stem (pad_head_length)
                           _|___|_
                          |x /  \ x|  SQUID (2 parallel JJs, lumped Lj each in HFSS)
                          |__\__/__|  ground side
                           ______|   L-shaped flux line

Claw geometry (left claw, local frame before the y-flip of ``_place``):
    back plate: x in [-L/2 - claw_width, -L/2]  (faces the qubit arm tip, claw_sep away)
    fingers:    along -x, length claw_length (measured from the back plate's inner edge),
                at |y| in [arm/2 + claw_sep, arm/2 + claw_sep + claw_width]
    etch:       claw buffered by claw_gap (merges with the qubit cross etch).
``qubit_arm_width`` must equal the Xmon ``cross_width``.
"""

from qiskit_metal import Dict, draw

from tmon_coupler_dolan_nthu import TmonCouplerDolanNthu


class ClawCouplerDolan(TmonCouplerDolanNthu):
    """TmonCouplerDolanNthu with U-shaped claw ends (overrides ``_end_geometry()`` only)."""

    default_options = Dict(
        TmonCouplerDolanNthu.default_options,
        qubit_arm_width='30um',
        claw_width='10um',
        claw_length='60um',
        claw_sep='8um',
        claw_gap='12um',
        flux_bias_line_options=Dict(
            TmonCouplerDolanNthu.default_options.flux_bias_line_options,
            fbl_type='L',
            mirror='False',
        ),
    )
    component_metadata = Dict(
        short_name='ClawTC',
        _qgeometry_table_path='True',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )
    TOOLTIP = """Grounded tunable transmon coupler with claw ends (Dolan bridge SQUID)."""

    def _claw(self, sign):
        """Claw at x = sign * L/2 (sign -1: left, opening toward -x)."""
        p = self.p
        L, cw, cl, cs = p.pad_arm_length, p.claw_width, p.claw_length, p.claw_sep
        half_in = p.qubit_arm_width / 2 + cs
        half_out = half_in + cw
        x_back_in = sign * L / 2
        x_back_out = sign * (L / 2 - cw)
        x_tip = sign * (L / 2 + cl)
        back = draw.box(min(x_back_in, x_back_out), -half_out, max(x_back_in, x_back_out), half_out)
        f_top = draw.box(min(x_back_in, x_tip), half_in, max(x_back_in, x_tip), half_out)
        f_bot = draw.box(min(x_back_in, x_tip), -half_out, max(x_back_in, x_tip), -half_in)
        return draw.union(back, f_top, f_bot)

    def _end_geometry(self):
        p = self.p
        bar = draw.rectangle(p.pad_arm_length - 2 * p.claw_width + 2e-6, p.pad_arm_width, 0, 0)
        claws = [self._claw(-1), self._claw(1)]
        etches = [c.buffer(p.claw_gap, join_style=2) for c in claws]
        return bar, claws, etches
