"""Single 6.5 GHz Cai-style Purcell filter and readout resonator for Design6."""
import argparse
from pathlib import Path
import sys

import os
os.environ.setdefault('SI_EPS', '11.5')      # Cai 2024: Si eps_r 11.5 (read by design5_layout / save_aedt at import)
import numpy as np
from scipy.optimize import root
from shapely.affinity import translate
from shapely.geometry import Point, box

HERE = Path(__file__).resolve().parent
DESIGN_ROOT = HERE.parent if HERE.name == 'pattern design' else Path.cwd().parent / 'Design6'
DESIGN5 = DESIGN_ROOT.parent / 'Design5'
sys.path[:0] = [str(DESIGN5 / directory) for directory in ('pattern design', 'simulation')]

import design5_layout as d5
import modes5
import theory5
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround

# Cai et al., Chin. Phys. B 33, 068501 (2024) parameters (Jay 2026-11-12): 16/9 um CPW (least fab-sensitive),
# T-type enclosing C_in = C_p = 40 fF with L = 12 um, C_j = 2 fF (two rectangular plates), f_r = f_p = 6.08 GHz.
CPW = dict(trace_width='16um', trace_gap='9um')
CHIP = dict(size_x='5.2mm', size_y='4.3mm', size_z='-500um',
            center_x='0mm', center_y='1.35mm')
PAD_X, CIN_X = 1.65, -0.66
CAP_PIN = d5.enc_xpin() / 1000
CP_Y = 0.15 + CAP_PIN
FILTER_TAP, RES_TAP, TEE_HALF = 1.20, 1.30, 0.15
TAP_Y = CP_Y + CAP_PIN + FILTER_TAP
RES_X, CJ_GAP = 0.22, 45
# Fig. 3(e) inset, dimensions estimated from the 50 um scale bar.
# C_j: Cai Fig. 3e two rectangular plates 25 x 65 um; gap 26 um = 2.00 fF (Design5 Q3D, 10/6 CPW) -> check with Q3D.
PAPER_HAMMER = dict(spacing_um=190, width_um=25, height_um=65, gap_um=42, ground_um=16)   # Q3D: 26 um 2.78, 36 um 2.20 fF -> 42 um ~2.0 fF
# C_in / C_p: T-type enclosing cap, characteristic gap L = 12 um (Cai); plate width sets C (Q3D-tuned)
ENC_CAI = dict(gap=12.0, plate_width=float(os.environ.get('ENC_PLATE_UM', 182.0)))
F0 = 6.08
LENGTHS = (4.535956670430899, 4.852581441726434)   # tune_lengths() at 6.08 GHz (circuit model, Design5 Q3D caps); check with HFSS


def tune_lengths():
    """Use the Design5 ABCD network and the actual C_in environment."""
    _, frame, plate, _ = theory5.q3d_caps('cin')
    _, ground_a, ground_b, _ = theory5.q3d_caps('cj')
    # Interpolate earlier Design6 Q3D values at 40 and 60 um in 1/gap.
    cj = np.interp(1 / CJ_GAP, [1 / 60, 1 / 40], [1.66, 2.25361])
    scale = np.sqrt(6.225 / ((1 + d5.SI_EPS) / 2))
    lead = d5.ENC['lead_length'] / 1000
    cap = (40.57409, frame, plate)
    p = dict(Z0=d5.Z0 * scale, v=d5.V * scale, cin=cap, cp=cap,
             cj=(cj, ground_a, ground_b), D_cj=FILTER_TAP, d_r=RES_TAP,
             l_feed_in=PAD_X + CIN_X - CAP_PIN + lead,
             l_feed_c=-CIN_X - CAP_PIN + lead, l_out=PAD_X, l_stub=0.15 + lead)

    def residual(lengths):
        p['L_filt'], p['L_res'] = map(float, lengths)
        result = modes5.extract(lambda f: theory5.branch_Y(f, p), p['Z0'], F0,
                                span=0.55, Yenv=lambda f: theory5.n0_env_Y(f, p))
        if result is None:
            raise RuntimeError('Could not resolve the two readout modes')
        return np.array([result['wp'], result['wr']]) - F0

    solution = root(residual, [4.19 * 6.5 / F0, 4.51 * 6.5 / F0])
    if not solution.success or max(abs(residual(solution.x))) > 1e-5:
        raise RuntimeError(f'Length tuning failed: {solution.message}')
    result = modes5.extract(lambda f: theory5.branch_Y(f, p), p['Z0'], F0,
                            span=0.55, Yenv=lambda f: theory5.n0_env_Y(f, p))
    return p, result


def position(x, y):
    return dict(pos_x=f'{x:.9f}mm', pos_y=f'{y:.9f}mm')


def straight(design, name, start, end):
    d5.RouteStraight(design, name, d5.Dict(
        pin_inputs=d5.Dict(start_pin=d5.Dict(component=start[0], pin=start[1]),
                           end_pin=d5.Dict(component=end[0], pin=end[1])), **CPW))


def enclosing(design, name, x, y, orientation):
    options = {key: f'{value}um' for key, value in dict(d5.ENC, **ENC_CAI).items()}
    d5.EnclosingCap(design, name, options=dict(
        options, **CPW, **position(x, y), orientation=str(orientation)))


def folded_tail(design, name, start, length, turn_run, side, direction, trim=0):
    """Five transverse runs and a grounded end. Dimensions are in mm."""
    x, y = design.components[start[0]].pins[start[1]].middle
    pitch, clearance, turns = 0.23, 0.25, 5
    radius = design.parse_value(d5.f3.FILLET)
    corner_loss = radius * (2 - np.pi / 2)
    span = (length - turn_run - turns * pitch - clearance + 2 * turns * corner_loss) / turns
    if min(span, pitch) < 2 * radius or turn_run < radius:
        raise ValueError(f'{name}: insufficient length for folded CPW')
    y += direction * turn_run
    anchors = [(x, y)]
    for index in range(turns):
        fold_x = x + side * (clearance + (span if index % 2 == 0 else 0))
        anchors.append((fold_x, y))
        y -= direction * pitch
        anchors.append((fold_x, y))
    end_x, end_y = anchors.pop()
    end_y -= direction * trim
    end = name + '_short'
    d5.ResEnd(design, end, options=dict(**position(end_x, end_y),
                                       orientation='90' if direction > 0 else '270'))
    d5.f3.route(design, name, start, (end, 'end'), anchors)


# qiskit-metal HFSS/Q3D auto wire bonds (renderer options): one bond per wb_threshold of every straight
# CPW segment, span wb_size x 22 um etch = 66 um, PEC, 15 um wire. The hammer stubs and tee arms are too short
# for auto bonds: add them with draw_hammer_straps(renderer, design, 'all'). Jay 2026-11-12.
WB_OPTIONS = dict(wb_threshold='350um', wb_offset='0um', wb_size=2)   # span 2 x 34 um etch = 68 um


def set_wire_bonds(design, on=True):
    """Flag every etched CPW path for qiskit-metal auto wire bonds (hfss_wire_bonds)."""
    t = design.qgeometry.tables['path']
    t['hfss_wire_bonds'] = t['subtract'].astype(bool) & on


def wire_bond_positions(design, threshold_mm=0.35):
    """Where qiskit-metal auto_wirebonds will put bonds (same rule), for plots: list of (x, y, ux, uy)."""
    out = []
    for row in design.qgeometry.tables['path'].itertuples():
        if not (row.subtract and row.hfss_wire_bonds):
            continue
        pts = np.asarray(row.geometry.coords)
        for a, b in zip(pts[:-1], pts[1:]):
            L = np.linalg.norm(b - a)
            if L > threshold_mm:
                n = int(L // threshold_mm)
                step = (b - a) / (n + 1)
                u = step / np.linalg.norm(step)
                out += [(*(a + step * (i + 1)), -u[1], u[0]) for i in range(n)]
    return out


MIN_CPW_SPACING = 0.40    # mm, Jay 2026-11-12: every pair of CPW lines >= 400 um centre to centre
WIDE = dict(step=0.40, y_bottom=0.65,
            filter=dict(lead=0.40, side=-1, direction=1),
            res=dict(lead=None, side=1, direction=-1))   # res lead: down to y_bottom


def serpentine_tail(design, name, start, length, lead, side, direction, height,
                    step=MIN_CPW_SPACING, trim=0, min_last=0.15):
    """Vertical serpentine with runs `step` apart (>= 400 um rule) and a grounded end.

    From the start pin: a vertical lead (along `direction`), then alternately a sideways step of
    `step` (towards `side`) and a vertical run of `height`, reversed each time; the last vertical
    run takes the remaining length. Length is the filleted centre-line length (mm); each 90 deg
    corner subtracts r (2 - pi/2). `trim` (mm) lengthens the last run.
    """
    x, y = design.components[start[0]].pins[start[1]].middle
    radius = design.parse_value(d5.f3.FILLET)
    loss = radius * (2 - np.pi / 2)
    d = direction
    y += d * lead
    anchors, straight = [(x, y)], lead
    while True:
        x += side * step
        anchors.append((x, y))
        straight += step
        d = -d
        remaining = length - straight + (len(anchors)) * loss    # corners = anchors after the pop
        if remaining <= height:
            break
        if remaining - height - step + 2 * loss < min_last:   # one more step would leave too little: run long
            break
        y += d * height
        anchors.append((x, y))
        straight += height
    if remaining < min_last:
        raise ValueError(f'{name}: last run {1000 * remaining:.0f} um too short, change lead/height')
    end_x, end_y = x, y + d * (remaining + trim)
    end = name + '_short'
    d5.ResEnd(design, end, options=dict(**position(end_x, end_y),
                                       orientation='270' if d > 0 else '90'))
    d5.f3.route(design, name, start, (end, 'end'), anchors)


def build(lengths=None, *, paper_hammer=False, paper_hammer_gap_um=None, res_trim=0,
          spacing_mm=None, wide=True):
    """Return a Qiskit Metal design with one filter and one readout resonator.

    spacing_mm overrides the filter-resonator centerline spacing (the hammer stays centred).
    The hammer bar (stub_f + pocket + stub_r) cuts the ground between the two lines; add
    hammer_ground_straps() in HFSS (draw_hammer_straps) or the coupling is not capacitive.
    """
    if lengths is None:
        p, _ = tune_lengths()
        lengths = p['L_filt'], p['L_res']
    filter_length, resonator_length = lengths
    if not -.12 <= res_trim <= .16:
        raise ValueError('Terminal segment trim outside validated geometry range')
    design = d5.designs.DesignPlanar()
    design.chips.main.material = d5.SI_MATERIAL
    design.chips.main.size.update(CHIP)
    design.variables.update(cpw_width=CPW['trace_width'], cpw_gap=CPW['trace_gap'])
    for name, x, orientation in [('LPPin', -PAD_X, '0'), ('LPPout', PAD_X, '180')]:
        d5.LaunchpadWirebond(design, name, options=dict(d5.f3.LP, **position(x, 0),
                                                     orientation=orientation))
    enclosing(design, 'cin', CIN_X, 0, 0)
    d5.TapJunction(design, 'n0', options=dict(**CPW, stub_length='100um', stub_up='True'))
    enclosing(design, 'cp', 0, CP_Y, 90)
    straight(design, 'feed_in', ('LPPin', 'tie'), ('cin', 'a'))
    straight(design, 'feed_mid', ('cin', 'b'), ('n0', 'prime_start'))
    straight(design, 'feed_out', ('n0', 'prime_end'), ('LPPout', 'tie'))
    straight(design, 'stub_p', ('n0', 'second_end'), ('cp', 'a'))

    spacing, width, height, gap, ground = RES_X, d5.CJ_HW, d5.CJ_HH, CJ_GAP, d5.CJ_GG
    if paper_hammer:
        spacing = PAPER_HAMMER['spacing_um'] / 1000
        width, height, gap, ground = (PAPER_HAMMER[key] for key in
                                      ('width_um', 'height_um', 'gap_um', 'ground_um'))
        if paper_hammer_gap_um is not None:
            gap = paper_hammer_gap_um
    if spacing_mm is not None:
        spacing = float(spacing_mm)
    hammer_pin = (gap / 2 + width + d5.CJ_LEAD) / 1000
    stub = spacing / 2 - hammer_pin - 0.02
    for name, x, stub_up in [('jf', 0, 'True'), ('jr', spacing, 'False')]:
        d5.TapJunction(design, name, options=dict(
            **CPW, **position(x, TAP_Y), orientation='-90',
            stub_up=stub_up, stub_length=f'{stub:.9f}mm'))
    d5.HammerCap(design, 'cj', options=dict(
        **CPW, **position(spacing / 2, TAP_Y), head_width=f'{width}um',
        head_height=f'{height}um', head_gap=f'{gap}um',
        lead_length=f'{d5.CJ_LEAD}um', ground_gap=f'{ground}um'))
    straight(design, 'stub_f', ('jf', 'second_end'), ('cj', 'a'))
    straight(design, 'stub_r', ('cj', 'b'), ('jr', 'second_end'))
    straight(design, 'feed_jf', ('cp', 'b'), ('jf', 'prime_end'))

    OpenToGround(design, 'readout_open', options=dict(
        **position(spacing, TAP_Y + RES_TAP), orientation='90',
        width=CPW['trace_width'], gap=CPW['trace_gap']))
    straight(design, 'res_open_run', ('jr', 'prime_start'), ('readout_open', 'open'))
    if wide:     # serpentines, every CPW pair >= 400 um apart (cpw_spacing_check.py)
        y_start = TAP_Y + TEE_HALF
        f_lead = WIDE['filter']['lead']
        serpentine_tail(design, 'filter_tail', ('jf', 'prime_start'),
                        filter_length - FILTER_TAP - TEE_HALF, f_lead, -1, 1,
                        height=y_start + f_lead - WIDE['y_bottom'])
        r_lead = (TAP_Y - TEE_HALF) - WIDE['y_bottom']
        serpentine_tail(design, 'res_tail', ('jr', 'prime_end'),
                        resonator_length - RES_TAP - TEE_HALF, r_lead, 1, -1,
                        height=TAP_Y + RES_TAP - WIDE['y_bottom'], trim=res_trim)
    else:        # original folded tails (230 um pitch, violates the 400 um rule)
        folded_tail(design, 'filter_tail', ('jf', 'prime_start'),
                    filter_length - FILTER_TAP - TEE_HALF, 0.52, -1, 1)
        folded_tail(design, 'res_tail', ('jr', 'prime_end'),
                    resonator_length - RES_TAP - TEE_HALF, 0.70, 1, -1, trim=res_trim)
    for name, target, fixed in [('filter_tail', filter_length, FILTER_TAP),
                                ('res_tail', resonator_length + res_trim, RES_TAP)]:
        actual = d5.f3.path_length(design, name) + fixed + TEE_HALF
        if abs(actual - target) > 0.001:
            raise ValueError(f'{name}: drawn length {actual:.6f} mm differs from {target:.6f} mm')
    return design


def _etched_region(design):
    """Union of every etched (subtract) shape, in mm."""
    from shapely.ops import unary_union
    shapes = [row.geometry.buffer(row.width / 2, cap_style=2)
              for row in design.qgeometry.tables['path'].itertuples() if row.subtract]
    shapes += [row.geometry for row in design.qgeometry.tables['poly'].itertuples()
               if row.subtract]
    return unary_union(shapes)


def hammer_ground_straps(design, span_um=60, diameter_um=15):
    """Ground straps across the two hammer stubs (Wallraff C_AB at the hammer T-junctions).

    Without them the etch of stub_f + hammer pocket + stub_r splits the ground between the
    filter and the resonator, and HFSS shows a ground-current coupling (J ~ 60 MHz at 190 um,
    hammer or not) that no capacitance describes. One strap per stub, at the middle of the
    straight etch between the main-line gap and the pocket. Returns a list of
    dict(pos=(x, y), ori=(ox, oy), span=..., diameter=...) in mm.
    """
    etch = _etched_region(design)
    xf, xr = (design.parse_value(design.components[name].options.pos_x) for name in ('jf', 'jr'))
    y = design.parse_value(design.components['jf'].options.pos_y)
    pocket = next(row.geometry for row in design.qgeometry.tables['poly'].itertuples()
                  if row.component == design.components['cj'].id and row.name == 'pocket')
    trace = design.parse_value(CPW['trace_width']) + 2 * design.parse_value(CPW['trace_gap'])
    px0, _, px1, _ = pocket.bounds
    span, radius = span_um / 1000, diameter_um / 2000
    straps = []
    # straight stub etch: from the main-line gap edge (centerline +- trace/2) to the pocket
    for x0, x1 in ((xf + trace / 2, px0), (px1, xr - trace / 2)):
        if x1 - x0 < 2 * radius + 0.004:
            raise ValueError(f'stub etch {1000 * (x1 - x0):.1f} um too short for a strap')
        x = (x0 + x1) / 2
        for foot_y in (y - span / 2, y + span / 2):
            foot = Point(x, foot_y).buffer(radius)
            if foot.intersects(etch):
                raise ValueError(f'strap foot at ({x:.4f}, {foot_y:.4f}) mm lands in an etch')
        straps.append(dict(pos=(x, y), ori=(0.0, 1.0), span=span, diameter=2 * radius))
    return straps


def junction_ground_straps(design, span_um=60, diameter_um=15, clearance_um=12):
    """Straps across the arms of the T-junctions (Wallraff C_AB at every t-junction).

    n0 (feed tee): both feed arms and the stub to the hugger; jf / jr: the two main-line arms
    (their stub arm is covered by hammer_ground_straps). qiskit-metal's auto_wirebonds skips all
    of these (segments < 400 um). Each strap sits clearance_um beyond the etch of the crossing
    line; feet are checked against the etch. Returns the hammer_ground_straps() format.
    """
    etch = _etched_region(design)
    half = (design.parse_value(CPW['trace_width']) + 2 * design.parse_value(CPW['trace_gap'])) / 2
    span, radius = span_um / 1000, diameter_um / 2000
    dist = half + radius + clearance_um / 1000          # from the junction centre along the arm
    arms = dict(n0=('prime_start', 'prime_end', 'second_end'), jf=('prime_start', 'prime_end'),
                jr=('prime_start', 'prime_end'))
    straps = []
    for name, pins in arms.items():
        comp = design.components[name]
        centre = np.array([design.parse_value(comp.options.pos_x), design.parse_value(comp.options.pos_y)])
        for pin in pins:
            u = comp.pins[pin].middle - centre
            u = u / np.linalg.norm(u)
            # the stub-arm strap must clear the feet of the two main-arm straps (HFSS: "parts intersect")
            stub_dist = max(dist + half, span / 2 + 3 * radius + clearance_um / 1000)
            pos = centre + u * (stub_dist if pin == 'second_end' else dist)
            ori = np.array([-u[1], u[0]])
            if any(Point(*foot).buffer(radius).intersects(etch)
                   for foot in (pos - ori * span / 2, pos + ori * span / 2)):
                continue      # e.g. the stub pocket of build_no_hammer; the CPW segment straps cover it
            straps.append(dict(pos=tuple(map(float, pos)), ori=tuple(map(float, ori)), span=span,
                               diameter=2 * radius, where=f'{name}.{pin}'))
    return straps


def cpw_ground_straps(design, pitch_um=200, span_um=60, diameter_um=15, end_clear_um=30):
    """Straps along every straight CPW segment (every etch path), at most pitch_um apart, at
    least one per segment where it fits. Interior vertices keep the fillet radius clear, path
    ends end_clear_um. Candidates whose feet land in an etch are skipped."""
    from shapely.geometry import LineString
    etch = _etched_region(design)
    span, radius = span_um / 1000, diameter_um / 2000
    names = {c.id: n for n, c in design.components.items()}
    straps = []
    for row in design.qgeometry.tables['path'].itertuples():
        if not row.subtract:
            continue
        pts = np.asarray(row.geometry.coords)
        fillet = row.fillet if row.fillet == row.fillet else 0.0          # nan -> 0
        for k in range(len(pts) - 1):
            a, b = pts[k], pts[k + 1]
            length = np.linalg.norm(b - a)
            u = (b - a) / length
            ca = fillet if k > 0 else end_clear_um / 1000
            cb = fillet if k < len(pts) - 2 else end_clear_um / 1000
            usable = length - ca - cb
            n = max(1, int(np.ceil(usable / (pitch_um / 1000)))) if usable > 0 else 1
            for j in range(n):
                s = ca + usable * (j + 0.5) / n if usable > 0 else length / 2
                pos, ori = a + u * s, np.array([-u[1], u[0]])
                feet = [pos - ori * span / 2, pos + ori * span / 2]
                if any(Point(*f).buffer(radius).intersects(etch) for f in feet):
                    continue
                if LineString(feet).length < span * 0.99:
                    continue
                straps.append(dict(pos=tuple(map(float, pos)), ori=tuple(map(float, ori)), span=span,
                                   diameter=2 * radius, where=f'{names[row.component]}.{row.name}[{k}]'))
    return straps


def _drop_crossing(straps, gap_um=5):
    """Keep straps in order; drop any whose footprint touches an already kept one."""
    from shapely.geometry import LineString
    kept, shapes = [], []
    for s in straps:
        p, o = np.asarray(s['pos']), np.asarray(s['ori']) * s['span'] / 2
        shape = LineString([p - o, p + o]).buffer(s['diameter'] / 2 + gap_um / 1000)
        if any(shape.intersects(x) for x in shapes):
            continue
        kept.append(s); shapes.append(shape)
    return kept


def ground_straps(design, which='all', pitch_um=200):
    """'hammer': the two hammer-stub straps; 'all': hammer + every T-junction arm;
    'cpw': hammer + T-junctions + every CPW segment (pitch_um) - use instead of auto_wirebonds."""
    cj = design.components['cj'].id
    no_hammer = not any(r.component == cj and r.name.startswith('head')
                        for r in design.qgeometry.tables['poly'].itertuples())
    # build_no_hammer keeps ground metal between the two stub pockets: no hammer straps needed there
    straps = [] if no_hammer else [dict(s, where=w) for s, w in
                                   zip(hammer_ground_straps(design), ('stub_f', 'stub_r'))]
    if which in ('all', 'cpw'):
        js = junction_ground_straps(design)
        if no_hammer:   # _widen_no_hammer moves jr geometry but not its options/pins: skip jf / jr
            js = [s for s in js if s['where'].startswith('n0.')]
        straps += js
    if which == 'cpw':
        straps += cpw_ground_straps(design, pitch_um=pitch_um)
    return _drop_crossing(straps)


def draw_hammer_straps(renderer, design, which='hammer', **kwargs):
    """Draw ground_straps(which) in an open HFSS/Q3D model, after renderer.render_design().

    Same bond-wire object as qiskit-metal's auto_wirebonds (PEC, name 'g_wb'), so q3d6.py
    still merges them into ground.
    """
    from pyEPR.ansys import parse_units      # mm -> renderer units, as in auto_wirebonds
    names = []
    for strap in ground_straps(design, which, **kwargs):
        names.append(renderer.modeler.draw_wirebond(
            pos=parse_units(np.asarray(strap['pos'])), ori=np.asarray(strap['ori']),
            width=parse_units(strap['span']), height=parse_units(strap['span']),
            z=0, wire_diameter=f"{strap['diameter']}mm", NumSides=6, name='g_wb',
            material='pec', solve_inside=False))
    return names


def build_no_hammer(lengths, *, spacing_mm=.400, res_trim=0):
    """Retain the CPW tap stubs and move the resonator and its local etch."""
    reference_spacing_mm = PAPER_HAMMER['spacing_um'] / 1000
    if not np.isfinite(spacing_mm) or spacing_mm < reference_spacing_mm:
        raise ValueError('This pattern requires at least 190 µm to preserve the tap geometry')
    design = build(lengths, paper_hammer=True, paper_hammer_gap_um=100, res_trim=res_trim)
    cj = design.components['cj']
    for kind in ('poly', 'path'):
        table = design.qgeometry.tables[kind]
        metal = (table['component'] == cj.id) & ~table['subtract'].astype(bool)
        table.drop(index=table.index[metal], inplace=True)
    if spacing_mm > reference_spacing_mm:
        _widen_no_hammer(design, spacing_mm - reference_spacing_mm)
    return design


def _widen_no_hammer(d, shift_mm):
    components = d.components
    res_ids = {components[name].id for name in (
        'jr', 'stub_r', 'res_open_run', 'readout_open', 'res_tail', 'res_tail_short')}
    cj = components['cj']
    for kind in ('path', 'poly', 'junction'):
        table = d.qgeometry.tables[kind]
        shifted = table['component'].isin(res_ids)
        if kind == 'path':
            shifted |= (table['component'] == cj.id) & table['name'].eq('lead_b_sub')
        table.loc[shifted, 'geometry'] = table.loc[shifted, 'geometry'].map(
            lambda geometry: translate(geometry, xoff=shift_mm))

    # Retain each stub's local etch, with ground metal between the two pockets.
    poly = d.qgeometry.tables['poly']
    idx, = poly.index[(poly['component'] == cj.id) & poly['name'].eq('pocket')]
    pocket = poly.at[idx, 'geometry']
    xmin, ymin, xmax, ymax = pocket.bounds
    middle = (xmin + xmax) / 2
    poly.at[idx, 'geometry'] = pocket.intersection(box(xmin, ymin, middle, ymax))
    right = pocket.intersection(box(middle, ymin, xmax, ymax))
    cj.add_qgeometry('poly', {'pocket_r': translate(right, xoff=shift_mm)},
                     subtract=True, layer=poly.at[idx, 'layer'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DESIGN_ROOT, help='parent of gds/ and fig/')
    parser.add_argument('--lengths', nargs=2, type=float, metavar=('FILTER_MM', 'RES_MM'),
                        help='export an explicit pair of centerline lengths in mm')
    parser.add_argument('--paper-hammer', action='store_true', help='use the inset dimensions')
    args = parser.parse_args()
    if args.paper_hammer and not args.lengths:
        parser.error('--paper-hammer requires explicit --lengths in mm')
    if args.lengths:
        lengths, result = args.lengths, None
        stem = 'purcell_6500_custom'
    else:
        p, result = tune_lengths()
        lengths = p['L_filt'], p['L_res']
        stem = 'purcell_6500'
    if args.paper_hammer:
        stem = 'purcell_6500_paper_hammer'
    design = build(lengths, paper_hammer=args.paper_hammer)
    gds = args.output / 'gds' / f'{stem}.gds'
    preview = args.output / 'fig' / f'{stem}.png'
    for path in (gds, preview):
        path.parent.mkdir(parents=True, exist_ok=True)
    d5.pl.export_gds(design, str(gds))
    d5.pl.plot_gds(str(gds), str(preview), window=(-2.6, 2.6, -0.8, 3.5))
    print(f'CPW 10/6 um, Si 500 um; L_filter={lengths[0]:.5f} mm, L_res={lengths[1]:.5f} mm')
    if result:
        print(f'Model: f_p={result["wp"]:.6f}, f_r={result["wr"]:.6f} GHz; '
              f'J/2pi={result["J"]:.2f}, kappa_p/2pi={result["kp"]:.2f} MHz')
    if args.paper_hammer:
        zoom = preview.with_name(f'{stem}_zoom.png')
        d5.pl.plot_gds(str(gds), str(zoom), window=(-.055, .245, TAP_Y-.10, TAP_Y+.10))
    print(gds, preview, sep='\n')


if __name__ == '__main__':
    main()





