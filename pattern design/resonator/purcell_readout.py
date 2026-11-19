r"""Purcell filter + readout resonator module joined only by a hammer capacitor C_j.

One quarter-wave Purcell filter (open at the C_p pin 'cp', shorted at its meander end) and one quarter-wave
readout resonator (open at the qubit pin 'qubit', shorted at its meander end). The two main lines pass
through one T junction each, ``line_spacing`` apart, and the hammer (``hammer_geometry.hammer_shapes``, the
same heads and pocket as the archived ``HammerCap``) joins the two T stubs. Frequencies are set by the two filleted lengths;
the filter-resonator coupling is set by the hammer: the two lines leave the T region in opposite directions,
so the only side-by-side run is the T region itself (2 x (tee_half + lead) at ``line_spacing``).

Sketch (orientation 0; F = filter, R = resonator, ## = etched CPW, [H] = hammer pocket)::

                    cp (pin, +y)  -> C_p / feed line
                     |
                     | filter open arm          R short-arm meander (runs +x, stacked
                     | (filter_tap)               up to resonator_meander_height)
                     |                          ____      ____
                     |              lead + jog |    |    |    |
                     |            +----------->+    |    |    |  short (end in ground)
                     |            |   R tee         |____|    |____
          F tee   ___|___   [H]   |  (+T..-T)
       tap  (-s/2) ==|==stub=|hh|=stub==|== tap (+s/2)      s = line_spacing
                     |                  |
         ____      __|                  | resonator open arm (resonator_tap)
        |    |    |     <-- jog + lead  |
        |    |    |    F short-arm      |
    ____|    |____|    meander (runs -x, stacked down to filter_meander_height)
                                        |
                                   qubit (pin, -y) -> qubit readout pad

* Filter: 'cp' -> straight down (filter_tap) -> tap -> lead down -> jog -x -> vertical meander (runs
  filter_meander_height, pitch meander_pitch) towards -x -> short. filter_length = filleted centre-line length
  from 'cp' to the short (stub excluded).
* Resonator: point mirror of the filter about the hammer centre: 'qubit' -> straight up (resonator_tap) ->
  tap -> lead up -> jog +x -> meander towards +x -> short. resonator_length likewise.
* Local frame: hammer centre at (0, 0), filter tap at (-s/2, 0), resonator tap at (+s/2, 0).

Options (lengths as qiskit strings):
    filter_length, filter_tap, resonator_length, resonator_tap: filleted lengths (open end -> short) and
        tap distances from the open ends
    line_spacing: centre spacing of the two main lines at the taps (default 400 um; check with parallel_j)
    tee_half: straight main line each side of a tap inside the T polygon;  lead: straight after the T
        before the first bend (>= fillet);  jog: sideways step before the meander (meander clearance)
    meander_pitch: run-to-run centre spacing (default 400 um);  filter_meander_height /
        resonator_meander_height: maximum run length; the runs are equalised so the length is exact
    filter_first_run / resonator_first_run: 'along' (first meander run continues away from the tap, as sketched)
        or 'back' (first run folds back past the tap: a tall, narrow channel, e.g. Design5's 2.6 mm pitch)
    mirror: 'True' mirrors the module in local x before rotation (filter on the +x side)
    resonator_open_extension: filleted length of a route the host draws from 'qubit' to the qubit's readout pin
        (e.g. a bend into a side claw). It is part of the resonator: resonator_length counts it and the tap sits
        resonator_tap + extension from the open end (use that tap in purcell_readout_model)
    fillet: bend radius (>= 3 W, B1);  trace_width / trace_gap: CPW 10/6 um (or 16/9 um)
    head_width, head_height, head_gap, hammer_ground_gap, corner_radius: hammer (HammerCap options)
    draw_mode: 'full' (default) or 'coupon' (Q3D coupon: only the two T regions, stubs and hammer, with
        straight main lines of coupon_lead each way from the taps and etched open ends)
    short_ends: 'short' (default, mask and HFSS) or 'open' (Q3D only: an etched trace_gap at each meander end, so
        the lines stay separate nets; method note section 3, 'render shorted ends open')

Pins: 'cp' (filter open end, normal +y) and 'qubit' (resonator open end, normal -y). Both open ends need a
component there (EnclosingCap / qubit readout pad); a bare pin ends the CPW in ground (a short).

Helpers: ``length_report`` (filleted lengths as drawn), ``strap_points()`` (ground straps / air bridges at the
hammer T junctions and on the stubs, design units), ``parallel_runs()`` (side-by-side filter-resonator runs for
``purcell_readout_model.parallel_j``), ``hammer_coupon()`` (Q3D-ready design). Circuit model (not a
QComponent): code/resonator/purcell_readout_model.py; layout strategy: note/purcell_filter_method.md section 6.
"""
import numpy as np
from shapely.ops import unary_union
from qiskit_metal import designs, draw, Dict
from qiskit_metal.qlibrary.core import QComponent

from hammer_geometry import hammer_shapes
from metal_width import min_metal_width, require_strip


def filleted_length(coords, r):
    """Length of a polyline after rounding every corner with radius r (as in the GDS): each corner of
    turning angle t is shortened by 2 r tan(t/2) - r t. Same as Design3/layout/ftf_layout.filleted_length."""
    c = np.asarray(coords, float)
    L = np.linalg.norm(np.diff(c, axis=0), axis=1).sum()
    for p0, p1, p2 in zip(c[:-2], c[1:-1], c[2:]):
        a, b = p1 - p0, p2 - p1
        t = np.arccos(np.clip(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1))
        L -= 2 * r * np.tan(t / 2) - r * t
    return float(L)


def _meander(start, u, w, length, lead, jog, pitch, height, r, first=1.0):
    """Polyline from `start` (tee end) heading `u`, then a jog along `w`, then runs parallel to u stepping by
    `pitch` along w; the first run goes along u (first=+1) or back along -u (first=-1, folds the meander back past
    the tap). Filleted length (radius r) equals `length`. Returns (points, runs)."""
    start, u, w = np.asarray(start, float), np.asarray(u, float), np.asarray(w, float)
    loss = r * (2 - np.pi / 2)                      # per 90 deg corner
    if lead < r or jog < 2 * r or pitch < 2 * r:
        raise ValueError('PurcellReadout: lead >= fillet, jog and meander_pitch >= 2 fillet required')
    for n in range(1, 200):
        corners = 2 * n
        run = (length - lead - jog - (n - 1) * pitch + corners * loss) / n
        if run <= height + 1e-12:
            break
    else:
        raise ValueError('PurcellReadout: meander needs more than 200 runs')
    if run < 2 * r:
        raise ValueError(f'PurcellReadout: run {run:.4g} shorter than 2 fillet; lower the meander height, '
                         'lead or jog, or lengthen the line')
    p = start + u * lead
    pts = [start, p]
    p = p + w * jog
    pts.append(p)
    d = float(first)
    for i in range(n):
        p = p + u * d * run
        pts.append(p)
        if i < n - 1:
            p = p + w * pitch
            pts.append(p)
        d = -d
    return np.array(pts), dict(runs=n, run=run)


class PurcellReadout(QComponent):
    """Quarter-wave Purcell filter + quarter-wave readout resonator + hammer C_j (see module docstring)."""

    component_metadata = Dict(short_name='purcell_resonator_pair', _qgeometry_table_poly='True', _qgeometry_table_path='True')
    default_options = Dict(
        filter_length='4.8mm', filter_tap='1.2mm',
        resonator_length='4.9mm', resonator_tap='1.3mm',
        line_spacing='400um', tee_half='50um', lead='100um', jog='400um',
        meander_pitch='400um', filter_meander_height='700um', resonator_meander_height='700um',
        fillet='99um', trace_width='10um', trace_gap='6um',
        head_width='25um', head_height='65um', head_gap='30um', hammer_ground_gap='16um', corner_radius='2um',
        filter_first_run='along', resonator_first_run='along', mirror='False', resonator_open_extension='0um',
        draw_mode='full', coupon_lead='400um', short_ends='short', layer='1')
    TOOLTIP = 'Purcell filter + readout resonator joined by a hammer C_j'

    def _mirrored(self):
        return str(self.p.mirror).lower() in ('true', '1')

    def _place(self, objs):
        p = self.p
        if self._mirrored():
            objs = draw.scale(objs, -1, 1, origin=(0, 0))
        objs = draw.rotate(objs, p.orientation, origin=(0, 0))
        return draw.translate(objs, p.pos_x, p.pos_y)

    def _xy(self, pts):
        """Local -> design coordinates for an (n, 2) array."""
        p = self.p
        a = np.radians(p.orientation)
        rot = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
        pts = np.asarray(pts, float) * ([-1, 1] if self._mirrored() else [1, 1])
        return pts @ rot.T + [p.pos_x, p.pos_y]

    def make(self):
        p = self.p
        tw, tg, r = p.trace_width, p.trace_gap, p.fillet
        s, T = p.line_spacing, p.tee_half
        xf, xr = -s / 2, s / 2
        head_a, head_b, pocket, x_head = hammer_shapes(p.head_width, p.head_height, p.head_gap,
                                                       p.hammer_ground_gap, p.corner_radius)
        unit = self.design.parse_value('1um')
        # ground between the T cut and the hammer pocket, and the T must hold the stub cut
        require_strip(self.name, 'ground between tee cut and hammer pocket',
                      (s / 2 - tw / 2 - tg) - (x_head + p.hammer_ground_gap), min_metal_width(unit))
        require_strip(self.name, 'tee_half - stub cut half width', T - (tw / 2 + tg), min_metal_width(unit))

        coupon = str(p.draw_mode).lower() == 'coupon'
        tee_len = 2 * p.coupon_lead if coupon else 2 * T
        tee_f = draw.LineString([(xf, -tee_len / 2), (xf, tee_len / 2)])
        tee_r = draw.LineString([(xr, -tee_len / 2), (xr, tee_len / 2)])
        stub_f = draw.LineString([(xf, 0), (-x_head + 1e-4 * unit, 0)])
        stub_r = draw.LineString([(x_head - 1e-4 * unit, 0), (xr, 0)])
        metal = [ln.buffer(tw / 2, cap_style=2) for ln in (tee_f, tee_r, stub_f, stub_r)]
        metal_f = unary_union([metal[0], metal[2], head_a])
        metal_r = unary_union([metal[1], metal[3], head_b])
        cut = unary_union([ln.buffer(tw / 2 + tg, cap_style=2) for ln in (tee_f, tee_r, stub_f, stub_r)] + [pocket])
        if coupon:          # etched open ends (end gap = trace gap)
            ends = [draw.rectangle(tw + 2 * tg, tg, x, sg * (tee_len / 2 + tg / 2))
                    for x in (xf, xr) for sg in (-1, 1)]
            cut = unary_union([cut] + ends)
        metal_f, metal_r, cut = self._place([metal_f, metal_r, cut])
        self.add_qgeometry('poly', {'filter_tee': metal_f, 'resonator_tee': metal_r}, layer=p.layer)
        self.add_qgeometry('poly', {'tee_cut': cut}, layer=p.layer, subtract=True)
        self.length_report = {}
        self._lines = {}
        if coupon:
            return

        ft, rt = p.filter_tap, p.resonator_tap
        if ft < T + p.lead + r or rt < T + p.lead + r:
            raise ValueError(f'{self.name}: filter_tap / resonator_tap must exceed tee_half + lead + fillet')
        ext = p.resonator_open_extension
        lf_short, lr_short = p.filter_length - ft, p.resonator_length - rt - ext
        # filter: short arm down (-y), meander towards -x; resonator: point mirror (+y, +x)
        f_pts, f_info = _meander((xf, -T), (0, -1), (-1, 0), lf_short - T, p.lead, p.jog, p.meander_pitch,
                                 p.filter_meander_height, r, self._first(p.filter_first_run))
        r_pts, r_info = _meander((xr, T), (0, 1), (1, 0), lr_short - T, p.lead, p.jog, p.meander_pitch,
                                 p.resonator_meander_height, r, self._first(p.resonator_first_run))
        f_open = np.array([(xf, T), (xf, ft)])
        r_open = np.array([(xr, -T), (xr, -rt)])
        if str(p.short_ends).lower() == 'open':
            gaps = []
            for pts in (f_pts, r_pts):
                end, prev = pts[-1], pts[-2]
                u = (end - prev) / np.linalg.norm(end - prev)
                gaps.append(draw.LineString([end, end + u * tg]).buffer(tw / 2 + tg, cap_style=2))
            gaps = self._place(gaps)
            self.add_qgeometry('poly', {'short_end_gaps': unary_union(gaps)}, layer=p.layer, subtract=True)
        for name, pts in (('filter_open', f_open), ('filter_short', f_pts),
                          ('resonator_open', r_open), ('resonator_short', r_pts)):
            line, = self._place([draw.LineString(pts)])
            self.add_qgeometry('path', {name: line}, width=tw, fillet=r, layer=p.layer)
            self.add_qgeometry('path', {name + '_cut': line}, width=tw + 2 * tg, fillet=r, layer=p.layer,
                               subtract=True)
        # full main lines from the open end (for lengths and side-by-side runs)
        self._lines = dict(filter=np.vstack([[(xf, ft)], [(xf, 0)], f_pts]),
                           resonator=np.vstack([[(xr, -rt - ext)], [(xr, 0)], r_pts]))   # extension drawn straight
        lf = ft + filleted_length(f_pts, r) + T
        lr = rt + filleted_length(r_pts, r) + T + ext
        self.length_report = dict(
            filter_mm=lf / self.design.parse_value('1mm'), resonator_mm=lr / self.design.parse_value('1mm'),
            filter_target_mm=p.filter_length / self.design.parse_value('1mm'),
            resonator_target_mm=p.resonator_length / self.design.parse_value('1mm'),
            filter_runs=f_info['runs'], filter_run_um=f_info['run'] / unit,
            resonator_runs=r_info['runs'], resonator_run_um=r_info['run'] / unit,
            stub_um=(s / 2 - x_head) / unit, tap_to_hammer_centre_um=s / 2 / unit,
            resonator_open_extension_mm=ext / self.design.parse_value('1mm'),
            resonator_tap_from_open_end_mm=(rt + ext) / self.design.parse_value('1mm'))
        for nm, (a, b) in (('cp', ((xf, ft - tw), (xf, ft))), ('qubit', ((xr, -rt + tw), (xr, -rt)))):
            self.add_pin(nm, points=self._xy([a, b]), width=tw, gap=tg, input_as_norm=True)

    @staticmethod
    def _first(value):
        v = str(value).lower()
        if v not in ('along', 'back'):
            raise ValueError("first_run must be 'along' or 'back'")
        return 1.0 if v == 'along' else -1.0

    # ---- helpers --------------------------------------------------------------------------------------
    def strap_points(self):
        """Ground straps / air bridges that keep the ground continuous around the hammer (method note section 5,
        Wallraff C_AB): across each stub halfway between the T cut and the pocket, and across both main lines
        just outside each T. Returns [dict(name, x, y, angle_deg, span)] in design units; angle is the strap
        direction (across the CPW), span covers the CPW etch plus 2 x 10 um landing."""
        p = self.p
        tw, tg, s, T = p.trace_width, p.trace_gap, p.line_spacing, p.tee_half
        unit = self.design.parse_value('1um')
        _, _, _, x_head = hammer_shapes(p.head_width, p.head_height, p.head_gap, p.hammer_ground_gap,
                                        p.corner_radius)
        x_stub = (s / 2 - tw / 2 - tg + x_head + p.hammer_ground_gap) / 2
        y_tee = T + tg + 10 * unit
        span = tw + 2 * tg + 20 * unit
        local = [('stub_f', -x_stub, 0, 90), ('stub_r', x_stub, 0, 90)]
        for nm, x in (('f', -s / 2), ('r', s / 2)):
            local += [(f'tee_{nm}_up', x, y_tee, 0), (f'tee_{nm}_down', x, -y_tee, 0)]
        out = []
        for nm, x, y, ang in local:
            (gx, gy), = self._xy([(x, y)])
            if self._mirrored():
                ang = 180 - ang
            out.append(dict(name=nm, x=float(gx), y=float(gy), angle_deg=(ang + float(p.orientation)) % 180,
                            span=float(span)))
        return out

    def parallel_runs(self, max_spacing_um=1000.0):
        """Side-by-side filter/resonator segment pairs (straight polyline segments, fillets ignored) closer than
        max_spacing_um: [dict(spacing_um, zp=(a, b), zr=(a, b), length_um)], z measured from each open end in
        mm (input of purcell_readout_model.parallel_j)."""
        unit_mm = self.design.parse_value('1mm')
        unit_um = self.design.parse_value('1um')
        out = []
        if not self._lines:
            return out

        def segments(pts):
            z = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
            return [(pts[i], pts[i + 1], z[i], z[i + 1]) for i in range(len(pts) - 1)]

        for a0, a1, za0, za1 in segments(self._lines['filter']):
            da = a1 - a0
            la = np.linalg.norm(da)
            ua = da / la
            for b0, b1, zb0, zb1 in segments(self._lines['resonator']):
                db = b1 - b0
                lb = np.linalg.norm(db)
                ub = db / lb
                if abs(abs(ua @ ub) - 1) > 1e-9:
                    continue
                nrm = np.array([-ua[1], ua[0]])
                spacing = abs((b0 - a0) @ nrm)
                if spacing / unit_um > max_spacing_um or spacing < 1e-12:
                    continue
                ta = sorted([0.0, la])
                tb = sorted([(b0 - a0) @ ua, (b1 - a0) @ ua])
                lo, hi = max(ta[0], tb[0]), min(ta[1], tb[1])
                if hi - lo <= 1e-12:
                    continue
                zp = (za0 + lo, za0 + hi)
                # resonator positions of the same projected interval
                tb0 = (b0 - a0) @ ua
                sgn = ua @ ub
                zr = (zb0 + (lo - tb0) * sgn, zb0 + (hi - tb0) * sgn)
                out.append(dict(spacing_um=spacing / unit_um, zp=tuple(z / unit_mm for z in zp),
                                zr=tuple(z / unit_mm for z in zr), length_um=(hi - lo) / unit_um))
        return out


def hammer_coupon(options=None, chip_size_mm=(3.0, 2.0)):
    """Q3D-ready design: the module's T regions, stubs and hammer with straight main lines of coupon_lead each way
    (etched open ends), on a chip of chip_size_mm. Use with QQ3DRenderer; C_j = -Maxwell(filter, resonator),
    extract as in Design4/simulation/resonator/hammer_deembed.py (straight-CPW reference for the excess ground)."""
    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.size.update(size_x=f'{chip_size_mm[0]}mm', size_y=f'{chip_size_mm[1]}mm', size_z='-500um')
    PurcellReadout(design, 'cj_coupon', options=dict(options or {}, draw_mode='coupon'))
    return design
