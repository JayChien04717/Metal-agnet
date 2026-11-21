# -*- coding: utf-8 -*-
"""Design9_CPW_Bandpass_Filter: 4-8 GHz 4-stub quarter-wave short-circuited
CPW bandpass filter using uniform 10/6 um CPW everywhere.

Topology (10.0 x 6.0 mm silicon chip, 10/6 um CPW, R = 99 um fillets):
    LP_IN (x = -4.5 mm, y = 0)
      -- feed_in --> T1 (x = -3.0 mm) ==[L12, lam/4]==> T2 (x = -1.0 mm)
                                      ==[L23, lam/4]==> T3 (x = +1.0 mm)
                                      ==[L34, lam/4]==> T4 (x = +3.0 mm)
      -- feed_out -> LP_OUT (x = +4.5 mm, y = 0)

    - 4 identical quarter-wave short-circuited CPW shunt stubs (S1..S4) branch upward
      (+y) from TapJunctions T1..T4, each making one wide symmetric U-fold (650 um pitch)
      and terminating in ShortToGround at y = +350 um.
    - 3 identical quarter-wave connecting sections (L12, L23, L34) loop downward (-y)
      between adjacent T-junctions in single wide U-loops (1.20 mm width).
    - Every signal trace and stub uses the exact same 10/6 um CPW (Z0 ~ 50.9 ohm):
      no width transitions (except standard input/output BondPads), no narrow coupling
      gaps, no interdigital fingers, and no vias.

Usage (from repository root, qmetal Python):
    python -X utf8 -B Design9_CPW_Bandpass_Filter/layout/build.py
    python -X utf8 -B Design9_CPW_Bandpass_Filter/layout/build.py --export-sim-gds
"""

from collections import OrderedDict
from pathlib import Path
import argparse
import json
import math
import sys
import warnings

import numpy as np

warnings.filterwarnings('ignore')

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO / 'code' / 'general')]
import project_paths  # noqa: E402

project_paths.configure(__file__)
DESIGN = project_paths.design_dir(__file__)
sys.path[:0] = [*project_paths.pattern_design_paths()]

from shapely.geometry import LineString, MultiPolygon, Polygon, box  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

from qiskit_metal import Dict, designs  # noqa: E402
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround  # noqa: E402
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors  # noqa: E402
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight  # noqa: E402

from bondpad import add_bond_pad  # noqa: E402
from tap_junction import TapJunction  # noqa: E402

METAL_COLOR = '#7d8a98'
EDGE_COLOR = '#2b3540'
BOX_COLOR = '#d98b00'


def load_config():
    return json.loads((DESIGN / 'parameters' / 'geometry.json').read_text(encoding='utf-8'))


def _um(val_um):
    return f'{val_um:g}um'


def _mm(val_mm):
    return f'{val_mm:.9f}mm'


def quarter_wave_mm(f_ghz, v_ph=1.2016e8):
    """Analytic quarter-wave length [mm] at frequency f_ghz [GHz]."""
    return (v_ph / (4.0 * f_ghz * 1e9)) * 1e3


def _assert_segments_ge_2fillet(name, pts, r_mm):
    """Assert every consecutive straight segment in ``pts`` is >= 2 * r_mm (DRC rule B2)."""
    min_seg = 2.0 * r_mm
    for i in range(len(pts) - 1):
        seg_len = float(np.linalg.norm(np.asarray(pts[i + 1]) - np.asarray(pts[i])))
        if seg_len < min_seg - 1e-9:
            raise ValueError(
                f'{name}: segment {i}->{i + 1} length {seg_len * 1e3:.1f} um '
                f'< 2*fillet ({min_seg * 1e3:.1f} um)'
            )


def _rendered_filleted_length_mm(design, route_comp, n_bends, r_mm):
    ptable = design.qgeometry.tables['path']
    row = ptable[(ptable['component'] == route_comp.id) & (~ptable['subtract'])].iloc[0]
    poly_len = float(row.geometry.length)
    fillet_loss = n_bends * r_mm * (2.0 - math.pi / 2.0)
    return poly_len - fillet_loss


def build_design(*, overrides=None, wire_bonds=True):
    """Build the ``Design9_CPW_Bandpass_Filter`` Qiskit Metal design."""
    geo = load_config()
    if overrides:
        geo = json.loads(json.dumps(geo))
        geo.update(overrides)

    chip = geo['chip']
    cpw_cfg = geo['cpw']
    flt = geo['filter_target']
    tap_cfg = geo['tap_junction']
    stub_cfg = geo['stubs']
    conn_cfg = geo['connecting_lines']

    r_mm = cpw_cfg['fillet_um'] / 1000.0
    tw_mm = cpw_cfg['trace_width_um'] / 1000.0
    tg_mm = cpw_cfg['trace_gap_um'] / 1000.0
    assert r_mm >= 3.0 * tw_mm - 1e-12, f'DRC B1 violated: fillet {r_mm*1e3:.1f} um < 3W'

    cpw = dict(trace_width=_um(cpw_cfg['trace_width_um']),
               trace_gap=_um(cpw_cfg['trace_gap_um']))

    design = designs.DesignPlanar(overwrite_enabled=True)
    design.chips.main.material = chip['material']
    design.chips.main.size.update(
        size_x=_mm(chip['size_x_mm']),
        size_y=_mm(chip['size_y_mm']),
        size_z=_um(chip['size_z_um']),
        center_x=_mm(chip['center_x_mm']),
        center_y=_mm(chip['center_y_mm']),
    )
    design.variables.update(cpw_width=cpw['trace_width'], cpw_gap=cpw['trace_gap'])

    # 1. Input and output wire-bond launch pads on y = 0 mm
    lps = {}
    for lp_name, lp_opt in geo['launch_pads'].items():
        lps[lp_name] = add_bond_pad(
            design,
            lp_name,
            dict(
                pos_x=_mm(lp_opt['pos_x_mm']),
                pos_y=_mm(lp_opt['pos_y_mm']),
                orientation=str(lp_opt['orientation_deg']),
            ),
            **cpw,
        )

    # 2. Four TapJunctions along y = 0 mm
    l_prime_mm = tap_cfg['prime_length_um'] / 1000.0
    s_tap_mm = tap_cfg['stub_length_um'] / 1000.0
    stub_up = bool(tap_cfg['stub_up'])
    stub_sgn = 1.0 if stub_up else -1.0

    x_nodes = [float(x) for x in stub_cfg['x_mm']]
    taps = []
    for idx, x_k in enumerate(x_nodes, start=1):
        t = TapJunction(
            design,
            f'T{idx}',
            options=Dict(
                chip='main',
                pos_x=_mm(x_k),
                pos_y=_mm(0.0),
                orientation='0',
                prime_length=_um(tap_cfg['prime_length_um']),
                stub_length=_um(tap_cfg['stub_length_um']),
                stub_up='True' if stub_up else 'False',
                **cpw,
            ),
        )
        taps.append(t)

    # 3. Straight 50-ohm input and output feeds
    RouteStraight(
        design,
        'feed_in',
        options=Dict(
            chip='main',
            pin_inputs=Dict(
                start_pin=Dict(component='LP_IN', pin='tie'),
                end_pin=Dict(component=taps[0].name, pin='prime_start'),
            ),
            **cpw,
        ),
    )
    RouteStraight(
        design,
        'feed_out',
        options=Dict(
            chip='main',
            pin_inputs=Dict(
                start_pin=Dict(component=taps[-1].name, pin='prime_end'),
                end_pin=Dict(component='LP_OUT', pin='tie'),
            ),
            **cpw,
        ),
    )

    # Quarter-wave target length at f0
    ell_qw_mm = quarter_wave_mm(flt['f0_GHz'], cpw_cfg['assumed_v_ph_m_s'])
    fillet_loss_90 = r_mm * (2.0 - math.pi / 2.0)

    # 4. Three quarter-wave connecting lines (L12, L23, L34) looping downward
    loop_dir = float(conn_cfg['loop_dir'])
    lead_off_mm = conn_cfg['lead_offset_x_um'] / 1000.0
    conn_scales = [float(v) for v in conn_cfg['electrical_scale']]
    conn_reports = []
    for idx in range(len(taps) - 1):
        name = f'L{idx + 1}{idx + 2}'
        t_left, t_right = taps[idx], taps[idx + 1]
        x_l_node, x_r_node = x_nodes[idx], x_nodes[idx + 1]
        dx_nodes = x_r_node - x_l_node
        target_mm = ell_qw_mm * conn_scales[idx]

        # Node-to-node filleted length = dx_nodes + 2 * H_conn - 4 * fillet_loss_90
        n_bends = 4
        h_conn = (target_mm + n_bends * fillet_loss_90 - dx_nodes) / 2.0
        y_loop = loop_dir * h_conn
        x_leg_l = x_l_node + lead_off_mm
        x_leg_r = x_r_node - lead_off_mm

        p_start = np.asarray(t_left.pins['prime_end']['middle'], float)
        p_end = np.asarray(t_right.pins['prime_start']['middle'], float)
        anchors = [
            np.array([x_leg_l, 0.0]),
            np.array([x_leg_l, y_loop]),
            np.array([x_leg_r, y_loop]),
            np.array([x_leg_r, 0.0]),
        ]
        _assert_segments_ge_2fillet(name, [p_start, *anchors, p_end], r_mm)

        route = RouteAnchors(
            design,
            name,
            options=Dict(
                chip='main',
                fillet=_um(cpw_cfg['fillet_um']),
                lead=Dict(start_straight='0um', end_straight='0um'),
                anchors=OrderedDict(enumerate(anchors)),
                pin_inputs=Dict(
                    start_pin=Dict(component=t_left.name, pin='prime_end'),
                    end_pin=Dict(component=t_right.name, pin='prime_start'),
                ),
                **cpw,
            ),
        )
        route_len_mm = _rendered_filleted_length_mm(design, route, n_bends, r_mm)
        total_node_to_node_mm = route_len_mm + l_prime_mm
        conn_reports.append(
            dict(
                name=name,
                target_len_mm=target_mm,
                actual_filleted_mm=total_node_to_node_mm,
                route_only_mm=route_len_mm,
                loop_depth_mm=h_conn,
                loop_width_mm=x_leg_r - x_leg_l,
                y_loop_mm=y_loop,
            )
        )

    # 5. Four quarter-wave short-circuited stubs (S1..S4) folding above y = 0
    fold_dirs = [float(v) for v in stub_cfg['fold_dir']]
    fold_w_mm = stub_cfg['fold_width_um'] / 1000.0
    y_short_abs_mm = stub_cfg['short_y_um'] / 1000.0
    stub_scales = [float(v) for v in stub_cfg['electrical_scale']]
    stub_reports = []
    for idx, (t_comp, x_k, fdir, sc) in enumerate(
        zip(taps, x_nodes, fold_dirs, stub_scales), start=1
    ):
        name = f'S{idx}'
        target_mm = ell_qw_mm * sc
        n_bends = 2
        # Node-to-short filleted length = H_stub + fold_w + (H_stub - y_short_abs) - 2*fillet_loss_90
        h_stub = (target_mm + n_bends * fillet_loss_90 - fold_w_mm + y_short_abs_mm) / 2.0
        y_top = stub_sgn * h_stub
        y_short = stub_sgn * y_short_abs_mm
        x_short = x_k + fdir * fold_w_mm

        p_start = np.asarray(t_comp.pins['second_end']['middle'], float)
        anchors = [
            np.array([x_k, y_top]),
            np.array([x_short, y_top]),
        ]
        p_end = np.array([x_short, y_short])
        _assert_segments_ge_2fillet(name, [p_start, *anchors, p_end], r_mm)

        # ShortToGround outward pin normal: orientation=90 -> (0, +1) when stub approaches from +y;
        # orientation=270 -> (0, -1) when stub approaches from -y.
        short_ori = '90' if stub_sgn > 0 else '270'
        short_name = f'{name}_short'
        short_comp = ShortToGround(
            design,
            short_name,
            options=dict(
                pos_x=_mm(x_short),
                pos_y=_mm(y_short),
                orientation=short_ori,
                width=cpw['trace_width'],
                gap=cpw['trace_gap'],
            ),
        )
        # Verify short pin normal points toward the incoming vertical leg
        expected_ny = 1.0 if stub_sgn > 0 else -1.0
        actual_normal = np.asarray(short_comp.pins['short']['normal'], float)
        if np.dot(actual_normal, np.array([0.0, expected_ny])) < 0.99:
            # Flip orientation if qiskit-metal convention differs
            short_ori = '270' if short_ori == '90' else '90'
            short_comp.options.orientation = short_ori
            short_comp.rebuild()

        route = RouteAnchors(
            design,
            name,
            options=Dict(
                chip='main',
                fillet=_um(cpw_cfg['fillet_um']),
                lead=Dict(start_straight='0um', end_straight='0um'),
                anchors=OrderedDict(enumerate(anchors)),
                pin_inputs=Dict(
                    start_pin=Dict(component=t_comp.name, pin='second_end'),
                    end_pin=Dict(component=short_name, pin='short'),
                ),
                **cpw,
            ),
        )
        route_len_mm = _rendered_filleted_length_mm(design, route, n_bends, r_mm)
        total_stub_mm = route_len_mm + s_tap_mm
        stub_reports.append(
            dict(
                name=name,
                target_len_mm=target_mm,
                actual_filleted_mm=total_stub_mm,
                route_only_mm=route_len_mm,
                stub_height_mm=h_stub,
                fold_width_mm=fold_w_mm,
                x_tap_mm=x_k,
                x_short_mm=x_short,
                y_top_mm=y_top,
                y_short_mm=y_short,
            )
        )

    if wire_bonds:
        t = design.qgeometry.tables['path']
        t['hfss_wire_bonds'] = t['subtract'].astype(bool)

    design.metadata['geometry_revision'] = geo['geometry_revision']
    design.metadata['config'] = geo
    design.metadata['connecting_lines'] = conn_reports
    design.metadata['stubs'] = stub_reports
    return design


# -----------------------------------------------------------------------------
# Exact filleted Shapely rendering (no GDS written!) + Stage 1 DRC checks
# -----------------------------------------------------------------------------
def _fillet_linestring(line, radius, n_arc=24):
    """Apply exact circular-arc fillets of ``radius`` to interior vertices of a LineString."""
    coords = [np.asarray(pt, dtype=float) for pt in line.coords]
    if radius <= 0.0 or len(coords) < 3:
        return line
    out = [coords[0]]
    for i in range(1, len(coords) - 1):
        p_prev, p_cur, p_next = coords[i - 1], coords[i], coords[i + 1]
        v1 = p_cur - p_prev
        v2 = p_next - p_cur
        l1 = float(np.linalg.norm(v1))
        l2 = float(np.linalg.norm(v2))
        if l1 < 1e-12 or l2 < 1e-12:
            continue
        u1 = v1 / l1
        u2 = v2 / l2
        cross = u1[0] * u2[1] - u1[1] * u2[0]
        dot = float(np.clip(np.dot(u1, u2), -1.0, 1.0))
        if abs(cross) < 1e-9:
            out.append(p_cur)
            continue
        turn = math.acos(dot)
        trim = radius * math.tan(turn / 2.0)
        q1 = p_cur - trim * u1
        q2 = p_cur + trim * u2
        sgn = 1.0 if cross > 0 else -1.0
        n_in = np.array([-sgn * u1[1], sgn * u1[0]])
        center = q1 + radius * n_in
        a1 = math.atan2(q1[1] - center[1], q1[0] - center[0])
        a2 = math.atan2(q2[1] - center[1], q2[0] - center[0])
        da = a2 - a1
        if sgn > 0 and da < 0:
            da += 2.0 * math.pi
        elif sgn < 0 and da > 0:
            da -= 2.0 * math.pi
        angles = np.linspace(a1, a1 + da, n_arc + 1)
        for ang in angles:
            out.append(center + radius * np.array([math.cos(ang), math.sin(ang)]))
    out.append(coords[-1])
    return LineString(out)


def design_to_metal_um(design):
    """Build exact chip positive-metal Shapely geometry in um directly from QGeometry."""
    sx = float(design.parse_value(design.chips.main.size.size_x)) * 1e3
    sy = float(design.parse_value(design.chips.main.size.size_y)) * 1e3
    cx = float(design.parse_value(design.chips.main.size.center_x)) * 1e3
    cy = float(design.parse_value(design.chips.main.size.center_y)) * 1e3
    chip_rect = box(cx - sx / 2.0, cy - sy / 2.0, cx + sx / 2.0, cy + sy / 2.0)

    pos_geoms = []
    etch_geoms = []

    def _scale_um(geom):
        from shapely.affinity import scale
        return scale(geom, xfact=1e3, yfact=1e3, origin=(0.0, 0.0))

    poly_tbl = design.qgeometry.tables['poly']
    for _, row in poly_tbl.iterrows():
        g = _scale_um(row['geometry'])
        if bool(row['subtract']):
            etch_geoms.append(g)
        else:
            pos_geoms.append(g)

    path_tbl = design.qgeometry.tables['path']
    for _, row in path_tbl.iterrows():
        line = _scale_um(row['geometry'])
        fillet_um = float(row['fillet']) * 1e3 if row['fillet'] else 0.0
        width_um = float(row['width']) * 1e3
        fline = _fillet_linestring(line, fillet_um)
        buf = fline.buffer(width_um / 2.0, cap_style=2, join_style=2, quad_segs=16)
        if bool(row['subtract']):
            etch_geoms.append(buf)
        else:
            pos_geoms.append(buf)

    etch_union = unary_union(etch_geoms) if etch_geoms else Polygon()
    pos_union = unary_union(pos_geoms) if pos_geoms else Polygon()
    metal = chip_rect.difference(etch_union).union(pos_union)
    return metal, pos_union, etch_union, [cx - sx / 2.0, cx + sx / 2.0, cy - sy / 2.0, cy + sy / 2.0]


def _polygon_to_mpl_path(poly):
    """Convert a Shapely Polygon (with optional interior holes) to a Matplotlib Path."""
    from matplotlib.path import Path as MplPath

    def _ring_codes(coords, ccw=True):
        pts = np.asarray(coords, dtype=float)
        # Ensure ring orientation: exterior CCW, interior CW for non-zero winding fill
        area2 = float(np.sum(pts[:-1, 0] * pts[1:, 1] - pts[1:, 0] * pts[:-1, 1]))
        if (ccw and area2 < 0) or ((not ccw) and area2 > 0):
            pts = pts[::-1]
        codes = np.full(len(pts), MplPath.LINETO, dtype=np.uint8)
        codes[0] = MplPath.MOVETO
        codes[-1] = MplPath.CLOSEPOLY
        return pts, codes

    verts_list, codes_list = [], []
    v, c = _ring_codes(poly.exterior.coords, ccw=True)
    verts_list.append(v)
    codes_list.append(c)
    for hole in poly.interiors:
        v, c = _ring_codes(hole.coords, ccw=False)
        verts_list.append(v)
        codes_list.append(c)
    return MplPath(np.vstack(verts_list), np.concatenate(codes_list))


def _draw_metal(ax, metal, window):
    from matplotlib.patches import PathPatch

    x0, x1, y0, y1 = window
    clipped = metal.intersection(box(x0, y0, x1, y1))
    if clipped.is_empty:
        polys = []
    elif isinstance(clipped, Polygon):
        polys = [clipped]
    else:
        polys = [g for g in getattr(clipped, 'geoms', []) if isinstance(g, Polygon) and g.area > 0]
    ax.set_facecolor('white')
    for poly in polys:
        if poly.area <= 0:
            continue
        ax.add_patch(
            PathPatch(
                _polygon_to_mpl_path(poly),
                facecolor=METAL_COLOR,
                edgecolor=METAL_COLOR,
                linewidth=0.25,
                zorder=2,
            )
        )
    ax.set(xlim=(x0, x1), ylim=(y0, y1), aspect='equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_color(EDGE_COLOR)
        side.set_linewidth(0.8)


def _scale_bar(ax, window):
    from matplotlib.patches import Rectangle

    x0, x1, y0, y1 = window
    span = x1 - x0
    nice = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
    length = max(v for v in nice if v <= span / 4.0)
    xs = x0 + 0.05 * span
    ys = y0 + 0.05 * (y1 - y0)
    ax.add_patch(Rectangle((xs, ys), length, 0.014 * (y1 - y0), fc='black', ec='none', zorder=5))
    label = f'{length / 1000:g} mm' if length >= 1000 else f'{length:g} µm'
    ax.text(
        xs + length / 2.0,
        ys + 0.032 * (y1 - y0),
        label,
        ha='center',
        va='bottom',
        fontsize=9,
        zorder=5,
        bbox=dict(fc='white', ec='none', alpha=0.85, pad=1.0),
    )


def _fit_window(window, low=0.75, high=2.5):
    x0, x1, y0, y1 = window
    w, h = x1 - x0, y1 - y0
    cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
    if w / h < low:
        w = h * low
    elif w / h > high:
        h = w / high
    return [round(cx - w / 2.0, 2), round(cx + w / 2.0, 2), round(cy - h / 2.0, 2), round(cy + h / 2.0, 2)]


def render_layout_figures(design, metal, chip_window):
    """Render standard repo-style chip.png and part close-ups directly from QGeometry (no GDS)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    rev = design.metadata['geometry_revision']
    fig_dir = DESIGN / 'fig'
    fig_dir.mkdir(parents=True, exist_ok=True)

    parts = {
        'launch_pad': {
            'label': 'LP_IN + feed_in + T1',
            'window': _fit_window([-4980, -2800, -650, 650]),
        },
        't_junction_stub': {
            'label': 'T3 + folded λ/4 stub S3',
            'window': _fit_window([750, 1950, -250, 2600]),
        },
        'connecting_loop': {
            'label': 'λ/4 connecting U-loop L23 (T2–T3)',
            'window': _fit_window([-1250, 1250, -1800, 320]),
        },
        'tap_closeup': {
            'label': 'T1 & S1_short (10/6 µm CPW)',
            'window': _fit_window([-3180, -2200, -180, 520]),
        },
    }

    # Whole chip figure
    w_c, h_c = chip_window[1] - chip_window[0], chip_window[3] - chip_window[2]
    fig, ax = plt.subplots(figsize=(9.2, 9.2 * h_c / w_c))
    _draw_metal(ax, metal, chip_window)
    for name, info in parts.items():
        x0, x1, y0, y1 = info['window']
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=BOX_COLOR, lw=1.4, zorder=4))
        ha = 'right' if name == 't_junction_stub' else 'left'
        tx = x1 if name == 't_junction_stub' else x0
        ax.text(
            tx,
            y1 + 30,
            f' {name} ',
            ha=ha,
            va='bottom',
            fontsize=8,
            color=BOX_COLOR,
            weight='bold',
            zorder=4,
            bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.5),
        )
    # Component annotations on whole chip for quick visual review
    for i, x_mm in enumerate(design.metadata['config']['stubs']['x_mm'], start=1):
        ax.text(x_mm * 1e3, -240, f'T{i}', ha='center', va='top', fontsize=8, color='#1b4f72',
                weight='bold', zorder=5, bbox=dict(fc='white', ec='none', alpha=0.75, pad=0.5))
    for rep in design.metadata['stubs']:
        xc = 0.5 * (rep['x_tap_mm'] + rep['x_short_mm']) * 1e3
        ax.text(xc, rep['y_top_mm'] * 1e3 + 120, rep['name'], ha='center', va='bottom', fontsize=8.5,
                color='#1b4f72', weight='bold', zorder=5, bbox=dict(fc='white', ec='none', alpha=0.75, pad=0.5))
    for rep, xm in zip(design.metadata['connecting_lines'], [-2.0, 0.0, 2.0]):
        ax.text(xm * 1e3, rep['y_loop_mm'] * 1e3 - 150, rep['name'], ha='center', va='top', fontsize=8.5,
                color='#1b4f72', weight='bold', zorder=5, bbox=dict(fc='white', ec='none', alpha=0.75, pad=0.5))
    ax.text(-4550, -560, 'LP_IN', ha='center', va='top', fontsize=8.5, color='#1b4f72', weight='bold',
            zorder=5, bbox=dict(fc='white', ec='none', alpha=0.75, pad=0.5))
    ax.text(+4550, -560, 'LP_OUT', ha='center', va='top', fontsize=8.5, color='#1b4f72', weight='bold',
            zorder=5, bbox=dict(fc='white', ec='none', alpha=0.75, pad=0.5))

    _scale_bar(ax, chip_window)
    fig.suptitle(f'{DESIGN.name} · {rev} · whole chip (10 × 6 mm, 10/6 µm CPW, draft preview)',
                 x=0.01, ha='left', fontsize=11, weight='bold')
    fig.tight_layout()
    chip_png = fig_dir / 'chip.png'
    fig.savefig(chip_png, dpi=300, facecolor='white')
    plt.close(fig)

    # Part close-ups
    rows = [('chip.png', 'Whole chip (10 × 6 mm)', chip_window)]
    for part_name, info in parts.items():
        win = info['window']
        wp, hp = win[1] - win[0], win[3] - win[2]
        fig, ax = plt.subplots(figsize=(6.5, max(2.8, min(8.0, 6.5 * hp / wp))))
        _draw_metal(ax, metal, win)
        _scale_bar(ax, win)
        fig.suptitle(f"{DESIGN.name} · {rev} · {info['label']}", x=0.02, ha='left', fontsize=9.5, weight='bold')
        fig.tight_layout()
        png = fig_dir / f'{part_name}.png'
        fig.savefig(png, dpi=220, facecolor='white')
        plt.close(fig)
        rows.append((png.name, info['label'], win))

    lines = [
        f'# {DESIGN.name} figures',
        '',
        'Stage 1 draft layout figures rendered directly from Qiskit Metal QGeometry (with exact circular-arc',
        'fillets) by [`layout/build.py`](../layout/build.py). **No GDS file has been generated yet**',
        '(per user instruction; awaiting layout review).',
        '',
        '- Grey is metal and white is bare silicon substrate.',
        '- Revision: `' + rev + '`',
        '',
        '| Figure | Shows | Window x0, x1, y0, y1 (µm) |',
        '|---|---|---|',
    ]
    for fname, label, win in rows:
        lines.append(f"| [{fname}]({fname}) | {label} | {', '.join(f'{v:g}' for v in win)} |")
    (fig_dir / 'README.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return [fig_dir / r[0] for r in rows]


def render_analytic_s_params(design):
    """Plot ideal lossless ABCD S21/S11 model for the 4-stub quarter-wave CPW BPF."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    geo = design.metadata['config']
    z0 = float(geo['cpw']['assumed_z0_ohm'])
    f0 = float(geo['filter_target']['f0_GHz'])
    stub_scales = [float(v) for v in geo['stubs']['electrical_scale']]
    conn_scales = [float(v) for v in geo['connecting_lines']['electrical_scale']]

    freqs = np.linspace(0.2, 13.8, 2721)
    theta = (np.pi / 2.0) * (freqs / f0)

    def abcd_tline(z, th):
        c, s = np.cos(th), np.sin(th)
        return np.array([[c, 1j * z * s], [1j * s / z, c]])

    def abcd_sc(z, th):
        t = np.tan(th)
        y = np.where(np.abs(t) < 1e-12, -1j * 1e12, -1j / (z * t))
        return np.array([[np.ones_like(th), np.zeros_like(th)], [y, np.ones_like(th)]])

    def mm(a, b):
        return np.array([
            [a[0, 0] * b[0, 0] + a[0, 1] * b[1, 0], a[0, 0] * b[0, 1] + a[0, 1] * b[1, 1]],
            [a[1, 0] * b[0, 0] + a[1, 1] * b[1, 0], a[1, 0] * b[0, 1] + a[1, 1] * b[1, 1]],
        ])

    m = abcd_sc(z0, theta * stub_scales[0])
    for cs, ss in zip(conn_scales, stub_scales[1:]):
        m = mm(m, abcd_tline(z0, theta * cs))
        m = mm(m, abcd_sc(z0, theta * ss))

    a, b, c, d = m[0, 0], m[0, 1], m[1, 0], m[1, 1]
    den = a + b / z0 + c * z0 + d
    s21_db = 20.0 * np.log10(np.maximum(np.abs(2.0 / den), 1e-12))
    s11_db = 20.0 * np.log10(np.maximum(np.abs((a + b / z0 - c * z0 - d) / den), 1e-12))

    pb = (freqs >= 4.0) & (freqs <= 8.0)
    idx3 = np.where(s21_db >= -3.0)[0]
    idx1 = np.where(s21_db >= -1.0)[0]
    summary = dict(
        f3_lo_GHz=float(freqs[idx3[0]]),
        f3_hi_GHz=float(freqs[idx3[-1]]),
        f1_lo_GHz=float(freqs[idx1[0]]),
        f1_hi_GHz=float(freqs[idx1[-1]]),
        min_s21_pb_dB=float(s21_db[pb].min()),
        max_s11_pb_dB=float(s11_db[pb].max()),
        s21_2GHz_dB=float(s21_db[np.argmin(np.abs(freqs - 2.0))]),
        s21_3GHz_dB=float(s21_db[np.argmin(np.abs(freqs - 3.0))]),
        s21_9GHz_dB=float(s21_db[np.argmin(np.abs(freqs - 9.0))]),
        s21_10GHz_dB=float(s21_db[np.argmin(np.abs(freqs - 10.0))]),
    )

    out_dir = DESIGN / 'notes' / 'fig'
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / 'analytic_s_params.png'

    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.axvspan(4.0, 8.0, color='#d6eaf8', alpha=0.55, label='Target passband (4–8 GHz)')
    ax.plot(freqs, s21_db, color='#1f618d', lw=2.0, label=r'$|S_{21}|$ (ideal ABCD)')
    ax.plot(freqs, s11_db, color='#c0392b', lw=1.6, ls='--', label=r'$|S_{11}|$ (ideal ABCD)')
    ax.axhline(-3.0, color='#555555', lw=0.8, ls=':')
    ax.axhline(-15.0, color='#999999', lw=0.8, ls=':')
    ax.set_xlim(0.2, 13.8)
    ax.set_ylim(-50, 2)
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('Magnitude (dB)')
    ax.set_title(
        'Design9 Analytic Pre-EM Model (Lossless TEM ABCD, Unsimulated)\n'
        f"4-stub λ/4 SC shunt-stub CPW (10/6 µm, Z0={z0:.1f} Ω): "
        f"4–8 GHz |S21| ≥ {summary['min_s21_pb_dB']:.2f} dB, |S11| ≤ {summary['max_s11_pb_dB']:.1f} dB",
        fontsize=10,
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc='lower center', frameon=True, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, facecolor='white')
    plt.close(fig)
    return out_png, summary


def verify_stage1_drc(design, metal):
    """Run Stage 1 layout DRC checks directly on QGeometry & Shapely metal."""
    geo = design.metadata['config']
    sx_mm = geo['chip']['size_x_mm']
    sy_mm = geo['chip']['size_y_mm']
    tw_um = geo['cpw']['trace_width_um']
    tg_um = geo['cpw']['trace_gap_um']
    cpw_outer_um = tw_um + 2.0 * tg_um

    print('\n=== Stage 1 DRC Verification ===')
    # 1. Filleted lengths (DRC B8)
    print('  [B8] Quarter-wave connecting lines (node-to-node filleted):')
    for rep in design.metadata['connecting_lines']:
        err_um = (rep['actual_filleted_mm'] - rep['target_len_mm']) * 1e3
        print(
            f"    {rep['name']}: target={rep['target_len_mm']:.4f} mm, "
            f"actual={rep['actual_filleted_mm']:.4f} mm, err={err_um:+.3f} um, "
            f"loop={rep['loop_width_mm']*1e3:.0f} x {rep['loop_depth_mm']*1e3:.1f} um"
        )
        assert abs(err_um) < 0.5, f"{rep['name']} length error {err_um:.3f} um"

    print('  [B8] Quarter-wave short-circuited stubs (node-to-short filleted):')
    for rep in design.metadata['stubs']:
        err_um = (rep['actual_filleted_mm'] - rep['target_len_mm']) * 1e3
        print(
            f"    {rep['name']}: target={rep['target_len_mm']:.4f} mm, "
            f"actual={rep['actual_filleted_mm']:.4f} mm, err={err_um:+.3f} um, "
            f"top={rep['y_top_mm']:.4f} mm, short=({rep['x_short_mm']:+.2f}, {rep['y_short_mm']:+.2f}) mm"
        )
        assert abs(err_um) < 0.5, f"{rep['name']} length error {err_um:.3f} um"

    # 2. Pin normals check (DRC B13)
    print('  [B13] Pin positions & outward normals:')
    for cname in ('LP_IN', 'LP_OUT', 'T1', 'T2', 'T3', 'T4', 'S1_short', 'S2_short', 'S3_short', 'S4_short'):
        comp = design.components[cname]
        for pname, pinfo in comp.pins.items():
            mid = np.asarray(pinfo['middle'], float)
            norm = np.asarray(pinfo['normal'], float)
            assert np.linalg.norm(norm) > 0.99, f'{cname}.{pname} zero normal'
            print(f'    {cname:8s}.{pname:11s} mid=({mid[0]:+7.3f}, {mid[1]:+7.3f}) mm  n=({norm[0]:+.1f}, {norm[1]:+.1f})')

    # 3. Parallel trench-to-trench clearance (DRC B7 >= 300 um)
    fold_trench_um = geo['stubs']['fold_width_um'] - cpw_outer_um
    short_to_bus_um = geo['stubs']['short_y_um'] - cpw_outer_um
    s1_s2_gap_um = (abs(geo['stubs']['x_mm'][1] - geo['stubs']['x_mm'][0]) * 1e3
                    - 2.0 * geo['stubs']['fold_width_um'] - cpw_outer_um)
    loop_leg_um = (abs(geo['stubs']['x_mm'][1] - geo['stubs']['x_mm'][0]) * 1e3
                   - 2.0 * geo['connecting_lines']['lead_offset_x_um'] - cpw_outer_um)
    loop_to_loop_um = 2.0 * geo['connecting_lines']['lead_offset_x_um'] - cpw_outer_um
    min_trench_um = min(fold_trench_um, short_to_bus_um, s1_s2_gap_um, loop_leg_um, loop_to_loop_um)
    print(f'  [B7] Trench-to-trench spacings: stub_legs={fold_trench_um:.0f} um, '
          f'short_to_bus={short_to_bus_um:.0f} um, S1-S2={s1_s2_gap_um:.0f} um, '
          f'loop_legs={loop_leg_um:.0f} um, loop-loop={loop_to_loop_um:.0f} um (min={min_trench_um:.0f} um >= 300 um)')
    assert min_trench_um >= 300.0, f'DRC B7 violated: min trench spacing {min_trench_um:.1f} um < 300 um'

    # 4. Top / bottom dicing edge keep-out (DRC B11 >= 500 um)
    top_margin_um = (sy_mm / 2.0 - max(r['y_top_mm'] for r in design.metadata['stubs'])) * 1e3 - cpw_outer_um / 2.0
    bot_margin_um = (sy_mm / 2.0 + min(r['y_loop_mm'] for r in design.metadata['connecting_lines'])) * 1e3 - cpw_outer_um / 2.0
    print(f'  [B11] Edge keep-out: top={top_margin_um:.1f} um, bottom={bot_margin_um:.1f} um (>= 500 um)')
    assert top_margin_um >= 500.0 and bot_margin_um >= 500.0, 'DRC B11 violated'

    # 5. Route crossings (DRC B5/B6): check pairwise disjointness of positive path traces
    path_tbl = design.qgeometry.tables['path']
    pos_paths = path_tbl[~path_tbl['subtract']]
    id_to_name = {c.id: c.name for c in design.components.values()}
    lines = [(id_to_name.get(int(row['component']), str(row['component'])),
              _fillet_linestring(row['geometry'], float(row['fillet']) if row['fillet'] else 0.0))
             for _, row in pos_paths.iterrows()]
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            inter = lines[i][1].intersection(lines[j][1])
            if not inter.is_empty:
                raise AssertionError(f'DRC B6 route crossing between {lines[i][0]} and {lines[j][0]}: {inter}')
    print(f'  [B5/B6] Route crossing check across {len(lines)} CPW routes: 0 crossings')

    # 6. Minimum metal width >= 4 um check (DRC C12a) via morphological opening (radius = 2 um)
    eroded = metal.buffer(-2.0, join_style=2)
    opened = eroded.buffer(+2.0, join_style=2)
    thin = metal.difference(opened)
    # Ignore tiny corner-rounding artifacts (< 5 um^2) from mitered vs buffered opening
    thin_parts = [g for g in getattr(thin, 'geoms', [thin]) if g.area > 5.0]
    print(f'  [C12a] Minimum metal width (>= 4 um) check on Shapely metal: {len(thin_parts)} violation(s)')
    assert len(thin_parts) == 0, f'Thin metal violations found: {[g.bounds for g in thin_parts[:5]]}'

    # 7. Ground-plane continuity check (DRC B9): ensure no isolated ground islands
    # Metal should consist of exactly 1 connected polygon because every stub is shorted to ground
    # AND wait: the center conductor is shorted to ground via S1..S4, so the entire metal is 1 polygon!
    n_metal_polys = len(metal.geoms) if isinstance(metal, MultiPolygon) else 1
    print(f'  [B9] Connected metal polygon count: {n_metal_polys} (1 = continuous ground + DC-shorted signal net, 0 isolated islands)')
    assert n_metal_polys == 1, f'Expected 1 connected metal polygon, got {n_metal_polys}'


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        '--export-sim-gds',
        action='store_true',
        help='Export simulation companion GDS (disabled by default during Stage 1 preview).',
    )
    args = parser.parse_args()

    design = build_design()
    metal, _, _, chip_window = design_to_metal_um(design)
    verify_stage1_drc(design, metal)

    written = render_layout_figures(design, metal, chip_window)
    print('\nRendered layout figures (no GDS written):')
    for p in written:
        print(f'  {p}')

    s_png, s_sum = render_analytic_s_params(design)
    print(f'\nRendered analytic pre-EM S-parameter plot: {s_png}')
    print(
        f"  Ideal ABCD summary: 3 dB = {s_sum['f3_lo_GHz']:.2f}-{s_sum['f3_hi_GHz']:.2f} GHz, "
        f"1 dB = {s_sum['f1_lo_GHz']:.2f}-{s_sum['f1_hi_GHz']:.2f} GHz, "
        f"4-8 GHz min |S21| = {s_sum['min_s21_pb_dB']:.2f} dB, max |S11| = {s_sum['max_s11_pb_dB']:.2f} dB"
    )

    if args.export_sim_gds:
        from gds_tools import export_gds
        import check_metal_width as cmw

        rev = design.metadata['geometry_revision']
        gds_dir = DESIGN / 'gds'
        gds_dir.mkdir(parents=True, exist_ok=True)
        sim_gds = gds_dir / f'{DESIGN.name}_{rev}_sim.gds'
        export_gds(design, str(sim_gds), fab=False)
        print(f'\nWrote simulation companion GDS: {sim_gds}')
        exempt = cmw.design_exemptions(design)
        thin = cmw.check(sim_gds, exempt=exempt)
        print(f'Metal-width check (sim GDS, >= 4 um): {len(thin)} violation(s)')


if __name__ == '__main__':
    main()
