# -*- coding: utf-8 -*-
r"""Floating double-pad tunable transmon pocket with SQUID landing electrodes
compatible with ``TransmonJunctionDolan`` (Dolan bridge e-beam mask).

Design-specific component (currently in ``Design7_Tunable_Transmon/layout/components/``;
candidate for promotion to ``pattern design/transmon/`` after verification).

Sketch (orientation = 0)::

                      readout (pin, +y when loc_H = +1)
                              | |
                        +-----+ +-----+   <- readout coupling pad (optional)
         +--------------|-------------|--------------+
         |              +-------------+              |
         |                                           |
         |    +---------------------------------+    |
         |    |           pad_top (+y)          |    |  <- rounded corners
         |    +----------+           +----------+    |
         |               |  |     |  |               |  <- upper SQUID legs
         |               |__|     |__|               |     (2x4 um hooks in fab)
         |              JJ1         JJ2              |  <- 5 um jj_cut, 30 um spacing
         |               |__|     |__|               |     (2x4 um hooks, +2 um x-offset)
         |               |  |     |  |               |  <- lower SQUID legs
         |    +----------+           +----------+    |
         |    |           pad_bot (-y)          |    |
         |    +---------------------------------+    |
         |                                           |
         +-------------------------------------------+
                       pocket (subtract=True)

Pins: one pin per entry in ``connection_pads`` (e.g. ``readout``), outward normal.
No flux bias line and no charge line (flux tuning via external coil; drive/readout via resonator).
"""

from shapely.geometry import Point, MultiPolygon
from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.core import BaseQubit
from metal_width import min_metal_width, require_strip


def _clean(geom, min_area=1e-10):
    """Drop zero-area slivers left by boolean ops (Ansys cannot draw them)."""
    parts = [g for g in getattr(geom, 'geoms', [geom]) if g.area > min_area]
    return parts[0] if len(parts) == 1 else MultiPolygon(parts)


def _rounded_rect(w, h, r, xoff=0.0, yoff=0.0):
    """Rectangle centered at (xoff, yoff) with corner radius r (0 -> sharp)."""
    rect = draw.rectangle(w, h, xoff, yoff)
    if r > 0:
        r_eff = min(r, w / 2 - 1e-9, h / 2 - 1e-9)
        rect = rect.buffer(-r_eff, join_style=1).buffer(r_eff, join_style=1)
    return rect


def _is_true(val):
    return str(val).lower() in ('true', '1')


class TunableTransmonPocket(BaseQubit):
    """Floating double-pad tunable transmon pocket with a 2-JJ SQUID between the pads.

    All arithmetic uses ``self.p`` (design units, mm).
    Always keep ``orientation='0'`` when placing ``TransmonJunctionDolan`` so every
    Dolan bridge shares the vertical finger-left evaporation direction.

    Default Options:
        * pad_width / pad_height / pad_gap: floating capacitor pads and vertical separation
        * pad_corner_radius: corner fillet on both capacitor pads (>= 10 um recommended)
        * pocket_width / pocket_height / pocket_corner_radius: etched ground cutout
        * squid_width / squid_loop_width: outer width and inner loop width of the two
          vertical SQUID landing legs spanning ``pad_gap``. Their centers sit at
          ``x = +- (squid_width + squid_loop_width) / 4`` -> ``squid_spacing = (squid_width + squid_loop_width) / 2``.
          With 46 um outer / 14 um loop, ``squid_spacing = 30 um`` and each leg is 16 um wide,
          leaving >= 5 um metal on both sides of the 2x4 um Dolan contact hooks (even with +2 um offset).
        * jj_cut: vertical gap between upper and lower SQUID landing faces (5 um for ``TransmonJunctionDolan``)
        * jj_width: lumped-junction sheet width in the QGeometry ``junction`` table
        * Lj_a / Lj_b: per-branch HFSS lumped inductances (nH)
        * jj_metal: ``'False'`` omits the SQUID landing legs (junction-free base pads only)
        * fab_contacts: ``'True'`` cuts 2x4 um Dolan contact hooks into the SQUID landing faces
          (lower hooks shifted by ``ground_notch_offset`` along +x), compatible with
          ``code/general/ebeam_junction_placement.py`` ``add_ebeam_junctions()``.
        * connection_pads: optional capacitive readout coupling pads above (``loc_H=+1``) or
          below (``loc_H=-1``) the transmon pads.
    """

    component_metadata = Dict(
        short_name='TTPocket',
        _qgeometry_table_path='True',
        _qgeometry_table_poly='True',
        _qgeometry_table_junction='True',
    )

    default_options = Dict(
        chip='main',
        pos_x='0um',
        pos_y='0um',
        orientation='0',
        layer='1',
        pad_width='450um',
        pad_height='120um',
        pad_gap='50um',
        pad_corner_radius='15um',
        pocket_width='650um',
        pocket_height='550um',
        pocket_corner_radius='25um',
        # SQUID landing electrodes (default: 30 um spacing, 16 um legs, 5 um cut -> TransmonJunctionDolan ready)
        squid_width='46um',
        squid_loop_width='14um',
        jj_cut='5um',
        jj_width='7um',
        Lj_a='16nH',
        Lj_b='16nH',
        gds_cell_name='junction_t0',
        jj_metal='True',
        # Dolan contact hooks for TransmonJunctionDolan
        fab_contacts='False',
        contact_notch_width='2um',
        contact_notch_depth='4um',
        ground_notch_offset='2um',
        # Optional readout coupling pads
        connection_pads=Dict(),
        _default_connection_pads=Dict(
            pad_gap='15um',
            pad_width='160um',
            pad_height='30um',
            pad_cpw_shift='0um',
            pad_cpw_extent='25um',
            cpw_width='10um',
            cpw_gap='6um',
            cpw_extend='50um',
            pocket_extent='5um',
            pocket_rise='30um',
            loc_W='0',
            loc_H='+1',
        ),
    )

    TOOLTIP = """Floating double-pad tunable transmon pocket with TransmonJunctionDolan-compatible SQUID hooks."""

    def _place(self, objs):
        p = self.p
        objs = draw.rotate(objs, p.orientation, origin=(0, 0))
        return draw.translate(objs, p.pos_x, p.pos_y)

    def _squid_params(self):
        p = self.p
        sw, lw = p.squid_width, p.squid_loop_width
        arm_w = (sw - lw) / 2.0
        x_jj = (lw / 2.0 + sw / 2.0) / 2.0
        return sw, lw, arm_w, x_jj

    def _squid_geometry(self):
        """Return (upper_legs, lower_legs, jj1, jj2) in local coordinates centered at (0, 0)."""
        p = self.p
        sw, lw, arm_w, x_jj = self._squid_params()
        unit = self.design.parse_value('1um')
        min_w = min_metal_width(unit)

        if arm_w <= 0:
            raise ValueError(f'{self.name}: squid_width ({sw}) must exceed squid_loop_width ({lw})')

        pad_edge = p.pad_gap / 2.0
        cut_half = p.jj_cut / 2.0
        overlap = min(p.pad_height / 4.0, 5.0 * unit)   # clean overlap into pad_top / pad_bot
        leg_h = (pad_edge + overlap) - cut_half
        if leg_h <= 0:
            raise ValueError(f'{self.name}: pad_gap ({p.pad_gap}) must exceed jj_cut ({p.jj_cut})')

        yc_up = cut_half + leg_h / 2.0
        yc_lo = -(cut_half + leg_h / 2.0)

        upper_legs = draw.union(
            draw.rectangle(arm_w, leg_h, -x_jj, yc_up),
            draw.rectangle(arm_w, leg_h, +x_jj, yc_up),
        )
        lower_legs = draw.union(
            draw.rectangle(arm_w, leg_h, -x_jj, yc_lo),
            draw.rectangle(arm_w, leg_h, +x_jj, yc_lo),
        )

        # Both JJ lines run from upper electrode (island) to lower electrode (ground-side in TransmonJunctionDolan frame)
        jj1 = draw.LineString([(-x_jj, +cut_half), (-x_jj, -cut_half)])
        jj2 = draw.LineString([(+x_jj, +cut_half), (+x_jj, -cut_half)])

        if _is_true(p.fab_contacts):
            nw = p.contact_notch_width
            nd = p.contact_notch_depth
            off = p.ground_notch_offset
            # Check 4 um minimum metal width on both sides of upper and lower notches
            up_margin = arm_w / 2.0 - nw / 2.0
            lo_margin = arm_w / 2.0 - nw / 2.0 - abs(off)
            require_strip(self.name, 'upper SQUID leg beside contact hook', up_margin, min_w)
            require_strip(self.name, 'lower SQUID leg beside offset contact hook', lo_margin, min_w)
            for x in (-x_jj, +x_jj):
                upper_legs = upper_legs.difference(
                    draw.box(x - nw / 2.0, cut_half, x + nw / 2.0, cut_half + nd)
                )
                lower_legs = lower_legs.difference(
                    draw.box(x + off - nw / 2.0, -(cut_half + nd), x + off + nw / 2.0, -cut_half)
                )
        else:
            require_strip(self.name, 'SQUID landing leg width', arm_w, min_w)

        return _clean(upper_legs), _clean(lower_legs), jj1, jj2

    def make(self):
        self.make_pocket()
        self.make_connection_pads()

    def make_pocket(self):
        p = self.p
        pw, ph, pg = p.pad_width, p.pad_height, p.pad_gap
        yc_pad = pg / 2.0 + ph / 2.0

        pad_top = _rounded_rect(pw, ph, p.pad_corner_radius, 0.0, +yc_pad)
        pad_bot = _rounded_rect(pw, ph, p.pad_corner_radius, 0.0, -yc_pad)
        pocket = _rounded_rect(p.pocket_width, p.pocket_height, p.pocket_corner_radius, 0.0, 0.0)

        upper_legs, lower_legs, jj1, jj2 = self._squid_geometry()

        if _is_true(p.jj_metal):
            pad_top = _clean(draw.union(pad_top, upper_legs))
            pad_bot = _clean(draw.union(pad_bot, lower_legs))

        pad_top, pad_bot, pocket, jj1, jj2 = self._place([pad_top, pad_bot, pocket, jj1, jj2])

        self.add_qgeometry('poly', dict(pad_top=pad_top, pad_bot=pad_bot), chip=p.chip, layer=p.layer)
        self.add_qgeometry('poly', dict(pocket=pocket), subtract=True, chip=p.chip, layer=p.layer)
        self.add_qgeometry('junction', dict(rect_jj1=jj1), width=p.jj_width,
                           gds_cell_name=p.gds_cell_name, hfss_inductance=p.Lj_a, chip=p.chip)
        self.add_qgeometry('junction', dict(rect_jj2=jj2), width=p.jj_width,
                           gds_cell_name=p.gds_cell_name, hfss_inductance=p.Lj_b, chip=p.chip)

    def make_connection_pads(self):
        for name in self.options.connection_pads:
            self.make_connection_pad(name)

    def make_connection_pad(self, name):
        """Capacitive coupling pad above (loc_H=+1) or below (loc_H=-1) the transmon pads."""
        p = self.p
        pc = p.connection_pads[name]
        loc_w = int(pc.loc_W)
        loc_h = int(pc.loc_H)
        if loc_h not in (-1, +1):
            raise ValueError(f'{self.name}.{name}: loc_H must be +1 or -1, got {loc_h}')

        cpw_w = pc.cpw_width
        cpw_g = pc.cpw_gap
        c_pw = pc.pad_width
        c_ph = pc.pad_height
        c_pg = pc.pad_gap

        # x position along pad width
        x_pad = loc_w * (p.pad_width / 2.0 - c_pw / 2.0) + pc.pad_cpw_shift
        # y position: outside the transmon pad by c_pg
        y_pad_outer = p.pad_gap / 2.0 + p.pad_height
        y_cpad_in = loc_h * (y_pad_outer + c_pg)
        y_cpad_out = loc_h * (y_pad_outer + c_pg + c_ph)
        y_pocket_edge = loc_h * (p.pocket_height / 2.0)
        y_pin = y_pocket_edge + loc_h * pc.cpw_extend

        c_pad = _rounded_rect(c_pw, c_ph, min(c_ph / 4.0, self.design.parse_value('5um')),
                              x_pad, (y_cpad_in + y_cpad_out) / 2.0)
        # Straight trace from coupling pad to pin outside pocket
        trace_y0 = min(y_cpad_out, y_pin)
        trace_y1 = max(y_cpad_out, y_pin)
        trace = draw.box(x_pad - cpw_w / 2.0, trace_y0, x_pad + cpw_w / 2.0, trace_y1)
        # Etch around trace from pocket edge to pin
        etch_y0 = min(y_pocket_edge, y_pin)
        etch_y1 = max(y_pocket_edge, y_pin)
        trace_etch = draw.box(x_pad - (cpw_w / 2.0 + cpw_g), etch_y0,
                              x_pad + (cpw_w / 2.0 + cpw_g), etch_y1)

        # Pin segment oriented so right-hand normal points outward (along loc_h * +y)
        if loc_h > 0:
            pin_line = draw.LineString([(x_pad - cpw_w / 2.0, y_pin), (x_pad + cpw_w / 2.0, y_pin)])
        else:
            pin_line = draw.LineString([(x_pad + cpw_w / 2.0, y_pin), (x_pad - cpw_w / 2.0, y_pin)])

        metal = _clean(draw.union(c_pad, trace))
        metal, trace_etch, pin_line = self._place([metal, trace_etch, pin_line])

        self.add_qgeometry('poly', {f'{name}_connector_pad': metal}, chip=p.chip, layer=p.layer)
        self.add_qgeometry('poly', {f'{name}_connector_etch': trace_etch}, subtract=True, chip=p.chip, layer=p.layer)
        self.add_pin(name, list(pin_line.coords), cpw_w, gap=cpw_g)

    def contact_hooks(self):
        """Return placed Dolan contact hooks in design units (mm) for ``add_ebeam_junctions()``.

        Hook names follow the shared SQUID contract in ``code/general/ebeam_junction_placement.py``:
        ``island_jj1`` / ``island_jj2`` (upper pad landing hooks) and
        ``ground_jj1`` / ``ground_jj2`` (lower pad landing hooks, offset by ``ground_notch_offset``).
        """
        p = self.p
        _, _, _, x_jj = self._squid_params()
        cut_half = p.jj_cut / 2.0
        nd = p.contact_notch_depth
        off = p.ground_notch_offset
        hooks = []
        for tag, x in (('jj1', -x_jj), ('jj2', +x_jj)):
            pts = self._place([
                Point(x, +cut_half),
                Point(x, +cut_half + nd),
                Point(x + off, -cut_half),
                Point(x + off, -(cut_half + nd)),
            ])
            xy = [(pt.x, pt.y) for pt in pts]
            hooks += [
                dict(name=f'island_{tag}', mouth=xy[0], end=xy[1], width=p.contact_notch_width),
                dict(name=f'ground_{tag}', mouth=xy[2], end=xy[3], width=p.contact_notch_width),
            ]
        return hooks


class TunableTransmonPocketFab(TunableTransmonPocket):
    """Fabrication variant of ``TunableTransmonPocket`` with Dolan contact hooks enabled."""

    default_options = Dict(TunableTransmonPocket.default_options, fab_contacts='True')
