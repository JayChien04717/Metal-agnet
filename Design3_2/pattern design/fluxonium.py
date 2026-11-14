# -*- coding: utf-8 -*-
r"""Simulation models of the MyFluxonium hairpin and main JJ.

Same geometry, options and loop as ``MyFluxonium`` (``MyFluxonium._loop``): pads + stems,
pocket, U-loop flux line, readout and charge line. The loop (55 x 10 um, 2 um wires):

* ``array_lumped_model='arms'`` uses two straight L/2 boundaries and a metal
  closing end ``array_end``.
* ``array_lumped_model='closure'`` uses metal arms and one straight
  ``array_closure`` boundary with the total array inductance.
* the JJ: the short ``jj_length`` gap in the middle of the lead line, at the pocket centre,
  ``rect_jj`` with ``L_j``.

In HFSS driven modal the flux line pin becomes Port 1 and ``rect_jj`` Port 2
(the EMF of the whole fluxonium loop):

    M = Im Z21 / (2 pi f)        (quasi-static, far below any resonance)

Sketch of the arms model (orientation = 0)::

                                    | lead (from pad_top)
               #=====array_top======|
     flux      |                    |
     U-loop    | array_end          x  <- rect_jj (L_j): Port 2
     tip ->    |                    |
               #=====array_bot======|
                                    | lead (to pad_bot)
       = sides of the U = lumped L, #/| U end = metal (array_end)

The U end faces the flux-line tip (fl_sep away), as in MyFluxonium.
Use ``port_list=[(name, 'flux_line', 50)]`` and
``jj_to_port=[(name, 'rect_jj', 50, False)]`` (port only) in
``QHFSSRenderer.render_design``; the selected array model is rendered as lumped-L boundaries.
With ``True`` L_j is drawn in parallel with the port, which closes the loop and
scales the port voltage by L_j / (L_j + l_inductance), so M would come out too small.

Q3D (Ec, J, g): ``flux_as_ground='True'`` draws the shorted flux line as part of the
ground plane (Q3D otherwise reports overlapping nets). Q3D skips the junction table, so the
arms model leaves a floating ``array_end_<name>`` for ``note/cap_matrix.py``
to eliminate. The closure model connects each half of the hairpin to its
own qubit electrode and leaves a gap across the lumped-array location.
For a capacitance model that omits the entire array, set
``array_lumped_model='arms'`` and ``array_end_metal='False'``.

Pins: 'flux_line', 'readout_line', 'charge_line'.
"""

from qiskit_metal import draw, Dict
from my_fluxonium import MyFluxonium


class FluxoniumArrayJJ(MyFluxonium):
    """MyFluxonium for simulation: same loop, flux_as_ground option for Q3D.

    Options of the loop (MyFluxonium defaults = real array):
        * lead_x, lead_width: the lead line (JJ at lead_x = pocket centre)
        * array_length: lead line -> centre of the U end (53 um -> 55 um outer);
          array_sep: arm centre to arm centre (8 um -> 10 um outer); array_width
        * jj_length, jj_width: JJ gap (Port 2) and its width
        * array_width: width of the lumped-L sheets
        * array_lumped_model: 'arms' or 'closure'
        * array_closure_gap: length of the straight closure sheet
        * array_end_metal: retain the metal return in the arms model
        * l_inductance: total array inductance
        * L_j: main JJ inductance

    .. image:: fluxonium.png
    """

    component_metadata = Dict(short_name='FluxoniumArrayJJ',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True')

    default_options = Dict(
        flux_as_ground='False',   # 'True' for Q3D: U-loop flux line drawn as ground
        array_lumped_model='arms',
        array_closure_gap='2um',
        array_end_metal='True',
    )

    TOOLTIP = """Fluxonium hairpin with lumped array inductance and a lumped main JJ."""

    def _loop(self):
        geometry = super()._loop()
        p = self.p
        if p.array_lumped_model == 'arms':
            if str(p.array_end_metal).lower() not in ('true', '1'):
                geometry.pop('array_end')
            return geometry
        if p.array_lumped_model != 'closure':
            raise ValueError("array_lumped_model must be 'arms' or 'closure'")
        if not 0 < p.array_closure_gap < p.array_sep - p.array_width:
            raise ValueError('array_closure_gap must leave metal on both sides of the gap')

        end_x = p.lead_x - p.array_length
        half_gap = p.array_closure_gap / 2
        outer_y = (p.array_sep + p.array_width) / 2
        for side, sign in (('top', 1), ('bot', -1)):
            arm = geometry.pop(f'array_{side}').buffer(p.array_width / 2, cap_style=2)
            closure = draw.rectangle(p.array_width, outer_y - half_gap,
                                     end_x, sign * (outer_y + half_gap) / 2)
            # Carry the horizontal metal into the return so the corner shares an area.
            elbow = draw.rectangle((p.lead_width + p.array_width) / 2, p.array_width,
                                   end_x + (p.lead_width - p.array_width) / 4,
                                   sign * p.array_sep / 2)
            lead = f'lead_{side}'
            geometry[lead] = draw.union(geometry[lead], arm, elbow, closure)
        geometry.pop('array_end')
        geometry['array_closure'] = draw.LineString([(end_x, half_gap), (end_x, -half_gap)])
        return geometry

    def _pocket_for_make(self):
        """Same loop as MyFluxonium (MyFluxonium._loop); for Q3D the shorted U-loop flux
        line is part of the ground plane (flux_as_ground)."""
        rect_pk = self._pocket_geom()
        if self._flux_as_ground() and self.p.flux_line_options.make_fl:
            rect_pk = draw.subtract(rect_pk, self._flux_line_geom())
        return rect_pk

    def _flux_as_ground(self):
        return str(self.p.flux_as_ground).lower() in ('true', '1')

    def make_flux_line(self):
        """flux_as_ground: the U-loop is already part of the ground plane, pin only."""
        if self._flux_as_ground():
            MyFluxonium._flux_pin(self)
        else:
            super().make_flux_line()

    def _flux_pin(self):
        """Flux pin with the CPW gap, so an HFSS port at this pin spans cpw_gap."""
        p = self.p
        pf = p.flux_line_options
        w, h2 = pf.cpw_width, pf.fl_height / 2
        x_wall = -p.pocket_width / 2
        a, b = (x_wall, h2 - w / 2), (x_wall, h2 + w / 2)
        if str(p.mirror).lower() in ('true', '1'):
            a, b = b, a
        line = self._place([draw.LineString([a, b])])[0]
        self.add_pin('flux_line', list(line.coords), w, gap=pf.cpw_gap)


if __name__ == '__main__':
    # layout check: one fluxonium with the FTF options, plots in fig/mutual/
    import os
    from qiskit_metal import designs
    from ftf_layout import FX_OPTIONS, FX_PER_QUBIT, CPW_WIDTH, CPW_GAP, plot_design

    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.variables['cpw_width'] = CPW_WIDTH
    design.variables['cpw_gap'] = CPW_GAP
    design.chips.main.size.size_x = '1mm'
    design.chips.main.size.size_y = '1mm'
    FluxoniumArrayJJ(design, 'FX1', options=dict(FX_OPTIONS, **FX_PER_QUBIT['FX1']))
    for pn, pin in design.components['FX1'].pins.items():
        print(pn, pin.middle, pin.normal, pin.width, pin.gap)

    fig = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'fig', 'mutual')
    os.makedirs(fig, exist_ok=True)
    plot_design(design, os.path.join(fig, 'fluxonium_arrayjj.png'), window=(-0.15, 0.15, -0.17, 0.23))
    plot_design(design, os.path.join(fig, 'fluxonium_arrayjj_zoom.png'),
                window=(-0.135, 0.015, -0.04, 0.04), dpi=160)
