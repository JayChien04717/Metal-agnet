r"""FluxoniumDelft base metal with recesses for the fabrication junction leads.

Same as ``FluxoniumDelft`` (fluxonium_delft.py) with ``jj_metal='False'`` (no leads, no hairpin
end) plus one Dolan contact notch in each facing stem end, as ``MyFluxoniumFab`` does for the
NTHU fluxonium. The e-beam JJ + array (``pattern design/ebeam_junction/``) land in the notches.

Sketch (orientation = 0, centre of FluxoniumDelft)::

                    ___
                   /   \    pad_top
                   \___/
                    | |   stem
                    |_ _|   <- notch: contact_notch_width x contact_notch_depth at x = lead_x

                            pad_gap: the fab draws the JJ + array here (array towards -x)
                     _ _
                    | | |   <- notch
                    | |   stem
                   /   \
                   \___/    pad_bot
"""
from qiskit_metal import Dict, draw
from shapely.geometry import Point

from fluxonium_delft import FluxoniumDelft


class FluxoniumDelftFab(FluxoniumDelft):
    """Cut a 2 µm wide, 4 µm deep notch into each facing stem end."""

    component_metadata = Dict(short_name='FluxoniumDelftFab',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True')

    default_options = Dict(jj_metal='False', contact_notch_width='2um',
                           contact_notch_depth='4um')

    def _stem_edges(self):
        g2 = self.p.pad_gap / 2
        return g2, -g2

    def _pads(self):
        top, bottom = super()._pads()
        p = self.p
        half_width = p.contact_notch_width / 2
        left, right = p.lead_x - half_width, p.lead_x + half_width
        depth = p.contact_notch_depth
        top_edge, bottom_edge = self._stem_edges()
        if not (p.stem_width / 2 > abs(p.lead_x) + half_width and depth < p.stem_length):
            raise ValueError('contact notch must lie inside the stem end')
        return (top.difference(draw.box(left, top_edge, right, top_edge + depth)),
                bottom.difference(draw.box(left, bottom_edge - depth, right, bottom_edge)))

    def contact_hooks(self):
        """Placed contact notches, in design units (mm), for drawing e-beam junctions on them.

        Returns ``{'hooks': [...], 'array_dir': (dx, dy)}``. Each hook is a dict with ``name``
        ('top' / 'bottom' in the local frame), ``mouth`` (notch opening on the stem end),
        ``end`` (closed notch end) and ``width``. ``array_dir`` is the placed unit vector of the
        local -x axis, the side the simulation array runs to. Call after the component is built.
        """
        p = self.p
        top_edge, bottom_edge = self._stem_edges()
        depth = p.contact_notch_depth
        pts = self._place([Point(p.lead_x, top_edge), Point(p.lead_x, top_edge + depth),
                           Point(p.lead_x, bottom_edge), Point(p.lead_x, bottom_edge - depth),
                           Point(0, 0), Point(-1, 0)])
        xy = [(q.x, q.y) for q in pts]
        hooks = [dict(name='top', mouth=xy[0], end=xy[1], width=p.contact_notch_width),
                 dict(name='bottom', mouth=xy[2], end=xy[3], width=p.contact_notch_width)]
        return dict(hooks=hooks, array_dir=(xy[5][0] - xy[4][0], xy[5][1] - xy[4][1]))
