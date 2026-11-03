# -*- coding: utf-8 -*-
r"""HFSS version of MyFluxonium, for the flux-line mutual inductance.

Same geometry and options as ``MyFluxonium`` (pads, pocket, U-loop flux line,
readout and charge line, JJ array 55 x 10 um with 2 um wires and leads). Following
``FluxoniumPocket`` (qlibrary/qubits/fluxonium_v2.py), both inductive elements are
junction-table lines with an inductance, drawn so that they only touch metal at
their ends (clean lumped boundaries in HFSS):

* the two array arms: lumped inductors, ``l_inductance / 2`` each, from the lead
  edge to the metal closing end;
* the main JJ: the ``jj_length`` gap bridged by ``rect_jj`` with ``L_j``.

In HFSS driven modal the flux line pin becomes Port 1 and ``rect_jj`` Port 2
(the EMF of the whole fluxonium loop):

    M = Im Z21 / (2 pi f)        (quasi-static, far below any resonance)

Sketch of the centre (orientation = 0)::

                         pad_top
                ____________|____________
                            |  lead (lead_width)
     flux   ________________|__
     U-loop |  #####array_top## |__  ### = array arm, lumped L = l_inductance/2
     tip -> | |              |     |
            | |  loop        |     x  <- rect_jj (L_j), jj_length gap: Port 2
            | |#####array_bot## __|
            |__________________ |
           ^ array_end (metal)
                            |
                ____________|____________
                         pad_bot

The array closing end faces the flux-line tip (fl_sep away), as in MyFluxonium.
Use ``port_list=[(name, 'flux_line', 50)]`` and
``jj_to_port=[(name, 'rect_jj', 50, False)]`` (port only) in
``QHFSSRenderer.render_design``; the array arms are rendered as lumped-L boundaries.
With ``True`` L_j is drawn in parallel with the port, which closes the loop and
scales the port voltage by L_j / (L_j + l_inductance), so M would come out too small.

Pins: 'flux_line', 'readout_line', 'charge_line'.
"""

from qiskit_metal import draw, Dict
from my_fluxonium import MyFluxonium


class FluxoniumArrayJJ(MyFluxonium):
    """MyFluxonium with lumped-L array arms and a JJ port gap (same options).

    Options used for the centre (MyFluxonium defaults = real array):
        * lead_x, lead_width: leads from the pads to the inner edge of the arms;
          the JJ branch wires and the closing end are lead_width wide
        * array_length: from the leads to the centre of the closing end (53 um ->
          55 um outer); array_sep: arm centre to arm centre (8 um -> 10 um outer);
          array_width: width of the lumped-L arms
        * jj_offset: x distance from the leads to the JJ branch
        * jj_length, jj_width: JJ gap (Port 2) and its width
        * l_inductance (split l_inductance/2 per arm), L_j
    """

    component_metadata = Dict(short_name='FluxoniumArrayJJ',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True')

    TOOLTIP = """MyFluxonium for HFSS: lumped-L array arms + JJ port (mutual inductance)."""

    def make_pocket(self):
        p = self.p
        g2 = p.pad_gap / 2
        lw = p.lead_width
        aw = p.array_width
        a2 = p.array_sep / 2                        # arm centre line
        al = p.array_length
        jl = p.jj_length
        x0 = p.lead_x
        xj = x0 + p.jj_offset

        pad_top = draw.rectangle(p.pad_width, p.pad_height, 0, g2 + p.pad_height / 2)
        pad_bot = draw.rectangle(p.pad_width, p.pad_height, 0, -g2 - p.pad_height / 2)
        # leads from each pad to the inner edge of the arms (as MyFluxonium)
        y_in = a2 - aw / 2
        lead_top = draw.rectangle(lw, g2 - y_in, x0, (g2 + y_in) / 2)
        lead_bot = draw.rectangle(lw, g2 - y_in, x0, -(g2 + y_in) / 2)
        # JJ branch: arm extensions to x = xj, stubs down / up to the JJ gap
        ext = xj + lw / 2 - x0
        jj_arm_top = draw.rectangle(ext, aw, x0 + ext / 2, a2)
        jj_arm_bot = draw.rectangle(ext, aw, x0 + ext / 2, -a2)
        stub = a2 + aw / 2 - jl / 2
        jj_stub_top = draw.rectangle(lw, stub, xj, (a2 + aw / 2 + jl / 2) / 2)
        jj_stub_bot = draw.rectangle(lw, stub, xj, -(a2 + aw / 2 + jl / 2) / 2)
        pad_top = draw.union(pad_top, lead_top, jj_arm_top, jj_stub_top)
        pad_bot = draw.union(pad_bot, lead_bot, jj_arm_bot, jj_stub_bot)

        # closing end (metal, as MyFluxonium) and the two arms as lumped inductors
        # between the lead edge and the closing end, drawn in the same direction
        array_end = draw.rectangle(lw, 2 * a2 + aw, x0 - al, 0)
        x_a, x_b = x0 - lw / 2, x0 - al + lw / 2
        array_top = draw.LineString([(x_a, a2), (x_b, a2)])
        array_bot = draw.LineString([(x_a, -a2), (x_b, -a2)])

        rect_jj = draw.LineString([(xj, jl / 2), (xj, -jl / 2)])
        rect_pk = self._pocket_geom()

        pad_top, pad_bot, array_end, rect_pk, rect_jj, array_top, array_bot = self._place(
            [pad_top, pad_bot, array_end, rect_pk, rect_jj, array_top, array_bot])

        self.add_qgeometry('poly', dict(pad_top=pad_top, pad_bot=pad_bot, array_end=array_end))
        self.add_qgeometry('poly', dict(rect_pk=rect_pk), subtract=True)
        self.add_qgeometry('junction', dict(rect_jj=rect_jj), width=p.jj_width,
                           hfss_inductance=p.L_j, q3d_inductance=p.L_j,
                           gds_cell_name=p.gds_cell_name)
        l_half = f"{float(str(p.l_inductance).replace('nH', '')) / 2}nH"
        self.add_qgeometry('junction', dict(array_top=array_top, array_bot=array_bot),
                           width=aw, hfss_inductance=l_half, q3d_inductance=l_half,
                           gds_cell_name=p.gds_cell_inductor)

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
