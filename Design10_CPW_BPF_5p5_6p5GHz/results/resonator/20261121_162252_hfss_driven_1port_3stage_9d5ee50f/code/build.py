# -*- coding: utf-8 -*-
"""Design10_CPW_BPF_5p5_6p5GHz: single-port (1-port) 3-stage 5.5-6.5 GHz
reflection CPW bandpass filter with >40 dB attenuation across 2-4 GHz,
built with 100% uniform 10/6 um CPW.

Topology (8.0 x 6.8 mm silicon chip, 10/6 um CPW, R = 99 um fillets):
    LP_IN (x = -3.50, y = +0.60 mm)  <-- ONLY 1 PHYSICAL PORT
      -- feed_in --> T1 (x = -2.20 mm, straight shorted stub S1 = 2.147 mm)
      ==[R1: 7.420 mm U-loop]==> T2 (x = -0.20 mm, straight shorted stub S2 = 0.826 mm)
      ==[R2: 8.587 mm U-loop]==> T3 (x = +1.80 mm, straight shorted stub S3 = 0.458 mm)
      ==[R3: 4.446 mm U-loop (lam/4 open)]==> R3_open (OpenToGround at x = +3.20 mm)

    - 100% pure 10/6 um CPW: ZERO interdigital capacitors, ZERO narrow gaps, ZERO vias,
      ZERO bends in any stub, and NO second physical port.
    - R3 is lambda/4 open-ended (mode-free below 5.5 GHz -> >44 dB attenuation across 2-4 GHz).

Usage (from repository root, qmetal Python):
    python -X utf8 -B Design10_CPW_BPF_5p5_6p5GHz/layout/build.py
    python -X utf8 -B Design10_CPW_BPF_5p5_6p5GHz/layout/build.py --export-sim-gds
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
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround  # noqa: E402
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


def _assert_segments_ge_2fillet(name, pts, r_mm):
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
    poly_len = float(row['geometry'].length)
    return poly_len - n_bends * r_mm * (2.0 - math.pi / 2.0)


def build_design(*, overrides=None, wire_bonds=True):
    """Build the 1-port 3-stage reflection BPF Qiskit Metal design."""
    geo = load_config()
    if overrides:
        geo = json.loads(json.dumps(geo))
        geo.update(overrides)

    chip = geo['chip']
    cpw_cfg = geo['cpw']
    tap_cfg = geo['tap_junction']
    stub_cfg = geo['stubs']
    res_cfg = geo['resonators']
    y_main = float(geo['y_main_mm'])

    r_mm = cpw_cfg['fillet_um'] / 1000.0
    tw_mm = cpw_cfg['trace_width_um'] / 1000.0
    assert r_mm >= 3.0 * tw_mm - 1e-12, 'DRC B1 violated'

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

    # 1. Single physical input launch pad (LP_IN)
    lp = geo['launch_pad']
    add_bond_pad(
        design,
        lp['name'],
        dict(
            pos_x=_mm(lp['pos_x_mm']),
            pos_y=_mm(lp['pos_y_mm']),
            orientation=str(lp['orientation_deg']),
        ),
        **cpw,
    )

    # 2. Three TapJunctions along y = y_main (T1, T2, T3)
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
                pos_y=_mm(y_main),
                orientation='0',
                prime_length=_um(tap_cfg['prime_length_um']),
                stub_length=_um(tap_cfg['stub_length_um']),
                stub_up='True' if stub_up else 'False',
                **cpw,
            ),
        )
        taps.append(t)

    # 3. Straight 50-ohm input feed LP_IN.tie -> T1.prime_start
    RouteStraight(
        design,
        'feed_in',
        options=Dict(
            chip='main',
            pin_inputs=Dict(
                start_pin=Dict(component=lp['name'], pin='tie'),
                end_pin=Dict(component=taps[0].name, pin='prime_start'),
            ),
            **cpw,
        ),
    )

    # 4. Three straight short-circuited inductive stubs (S1, S2, S3, zero bends)
    stub_lens = [float(v) for v in stub_cfg['lengths_mm']]
    stub_reports = []
    for idx, (t_comp, x_k, target_mm) in enumerate(zip(taps, x_nodes, stub_lens), start=1):
        name = f'S{idx}'
        y_short = y_main + stub_sgn * target_mm
        short_name = f'{name}_short'
        short_ori = '270' if stub_sgn > 0 else '90'
        short_comp = ShortToGround(
            design,
            short_name,
            options=dict(
                pos_x=_mm(x_k),
                pos_y=_mm(y_short),
                orientation=short_ori,
                width=cpw['trace_width'],
                gap=cpw['trace_gap'],
            ),
        )
        expected_ny = -1.0 if stub_sgn > 0 else 1.0
        if np.dot(np.asarray(short_comp.pins['short']['normal'], float), np.array([0.0, expected_ny])) < 0.99:
            short_comp.options.orientation = '90' if short_ori == '270' else '270'
            short_comp.rebuild()

        route = RouteStraight(
            design,
            name,
            options=Dict(
                chip='main',
                pin_inputs=Dict(
                    start_pin=Dict(component=t_comp.name, pin='second_end'),
                    end_pin=Dict(component=short_name, pin='short'),
                ),
                **cpw,
            ),
        )
        total_stub_mm = _rendered_filleted_length_mm(design, route, 0, r_mm) + s_tap_mm
        stub_reports.append(dict(
            name=name,
            target_len_mm=target_mm,
            actual_filleted_mm=total_stub_mm,
            x_mm=x_k,
            y_short_mm=y_short,
        ))

    # 5. Open termination at the far end of Stage 3 (R3_open) -- NO second launch pad!
    x_open = float(res_cfg['x_open_mm'])
    open_comp = OpenToGround(
        design,
        'R3_open',
        options=dict(
            pos_x=_mm(x_open),
            pos_y=_mm(y_main),
            orientation='0',  # incoming from -x -> outward normal (-1, 0)
            width=cpw['trace_width'],
            gap=cpw['trace_gap'],
            termination_gap=cpw['trace_gap'],
        ),
    )
    if np.dot(np.asarray(open_comp.pins['open']['normal'], float), np.array([-1.0, 0.0])) < 0.99:
        open_comp.options.orientation = '180'
        open_comp.rebuild()

    # 6. Three resonator U-loops (R1: T1->T2, R2: T2->T3, R3: T3->R3_open) below y = y_main
    loop_dir = float(res_cfg['loop_dir'])
    lead_offs = [float(v) / 1000.0 for v in res_cfg['lead_offset_x_um']]
    res_lens = [float(v) for v in res_cfg['lengths_mm']]
    fillet_loss_90 = r_mm * (2.0 - math.pi / 2.0)

    res_reports = []
    for idx, target_mm in enumerate(res_lens):
        name = f'R{idx + 1}'
        t_left = taps[idx]
        x_l_node = x_nodes[idx]
        is_last = (idx == len(res_lens) - 1)
        if not is_last:
            t_right = taps[idx + 1]
            x_r_node = x_nodes[idx + 1]
            end_pin_spec = Dict(component=t_right.name, pin='prime_start')
            p_end = np.asarray(t_right.pins['prime_start']['middle'], float)
            tap_contrib_mm = l_prime_mm  # half from t_left + half from t_right
        else:
            x_r_node = x_open
            end_pin_spec = Dict(component='R3_open', pin='open')
            p_end = np.asarray(open_comp.pins['open']['middle'], float)
            tap_contrib_mm = l_prime_mm / 2.0  # half from t_left only; open pin is exact end

        dx_nodes = x_r_node - x_l_node
        n_bends = 4
        h_loop = (target_mm + n_bends * fillet_loss_90 - dx_nodes) / 2.0
        y_loop = y_main + loop_dir * h_loop
        x_leg_l = x_l_node + lead_offs[idx]
        x_leg_r = x_r_node - lead_offs[idx]

        p_start = np.asarray(t_left.pins['prime_end']['middle'], float)
        anchors = [
            np.array([x_leg_l, y_main]),
            np.array([x_leg_l, y_loop]),
            np.array([x_leg_r, y_loop]),
            np.array([x_leg_r, y_main]),
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
                    end_pin=end_pin_spec,
                ),
                **cpw,
            ),
        )
        route_len_mm = _rendered_filleted_length_mm(design, route, n_bends, r_mm)
        total_res_mm = route_len_mm + tap_contrib_mm
        res_reports.append(dict(
            name=name,
            target_len_mm=target_mm,
            actual_filleted_mm=total_res_mm,
            loop_depth_mm=h_loop,
            loop_width_mm=x_leg_r - x_leg_l,
            x_leg_l_mm=x_leg_l,
            x_leg_r_mm=x_leg_r,
            y_loop_mm=y_loop,
        ))

    if wire_bonds:
        t = design.qgeometry.tables['path']
        t['hfss_wire_bonds'] = t['subtract'].astype(bool)

    design.metadata['geometry_revision'] = geo['geometry_revision']
    design.metadata['config'] = geo
    design.metadata['stubs'] = stub_reports
    design.metadata['resonators'] = res_reports
    return design


# -----------------------------------------------------------------------------
# Exact filleted Shapely rendering (NO GDS written) + Stage 1 DRC checks
# -----------------------------------------------------------------------------
def _fillet_linestring(line, radius, n_arc=24):
    coords = [np.asarray(pt, dtype=float) for pt in line.coords]
    if radius <= 0.0 or len(coords) < 3:
        return line
    out = [coords[0]]
    for i in range(1, len(coords) - 1):
        p_prev, p_cur, p_next = coords[i - 1], coords[i], coords[i + 1]
        v1, v2 = p_cur - p_prev, p_next - p_cur
        l1, l2 = float(np.linalg.norm(v1)), float(np.linalg.norm(v2))
        if l1 < 1e-12 or l2 < 1e-12:
            continue
        u1, u2 = v1 / l1, v2 / l2
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
        for ang in np.linspace(a1, a1 + da, n_arc + 1):
            out.append(center + radius * np.array([math.cos(ang), math.sin(ang)]))
    out.append(coords[-1])
    return LineString(out)


def design_to_metal_um(design):
    from shapely.affinity import scale

    sx = float(design.parse_value(design.chips.main.size.size_x)) * 1e3
    sy = float(design.parse_value(design.chips.main.size.size_y)) * 1e3
    cx = float(design.parse_value(design.chips.main.size.center_x)) * 1e3
    cy = float(design.parse_value(design.chips.main.size.center_y)) * 1e3
    chip_rect = box(cx - sx / 2.0, cy - sy / 2.0, cx + sx / 2.0, cy + sy / 2.0)

    def _s(geom):
        return scale(geom, xfact=1e3, yfact=1e3, origin=(0.0, 0.0))

    pos_geoms, etch_geoms = [], []
    for _, row in design.qgeometry.tables['poly'].iterrows():
        g = _s(row['geometry'])
        (etch_geoms if bool(row['subtract']) else pos_geoms).append(g)

    for _, row in design.qgeometry.tables['path'].iterrows():
        line = _s(row['geometry'])
        fillet_um = float(row['fillet']) * 1e3 if row['fillet'] else 0.0
        width_um = float(row['width']) * 1e3
        buf = _fillet_linestring(line, fillet_um).buffer(width_um / 2.0, cap_style=2, join_style=2, quad_segs=16)
        (etch_geoms if bool(row['subtract']) else pos_geoms).append(buf)

    etch_union = unary_union(etch_geoms) if etch_geoms else Polygon()
    pos_union = unary_union(pos_geoms) if pos_geoms else Polygon()
    metal = chip_rect.difference(etch_union).union(pos_union)
    return metal, [cx - sx / 2.0, cx + sx / 2.0, cy - sy / 2.0, cy + sy / 2.0]


def _polygon_to_mpl_path(poly):
    from matplotlib.path import Path as MplPath

    def _ring(coords, ccw=True):
        pts = np.asarray(coords, dtype=float)
        area2 = float(np.sum(pts[:-1, 0] * pts[1:, 1] - pts[1:, 0] * pts[:-1, 1]))
        if (ccw and area2 < 0) or ((not ccw) and area2 > 0):
            pts = pts[::-1]
        codes = np.full(len(pts), MplPath.LINETO, dtype=np.uint8)
        codes[0] = MplPath.MOVETO
        codes[-1] = MplPath.CLOSEPOLY
        return pts, codes

    verts, codes = [], []
    v, c = _ring(poly.exterior.coords, ccw=True)
    verts.append(v); codes.append(c)
    for hole in poly.interiors:
        v, c = _ring(hole.coords, ccw=False)
        verts.append(v); codes.append(c)
    return MplPath(np.vstack(verts), np.concatenate(codes))


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
        ax.add_patch(PathPatch(_polygon_to_mpl_path(poly), fc=METAL_COLOR, ec=METAL_COLOR, lw=0.25, zorder=2))
    ax.set(xlim=(x0, x1), ylim=(y0, y1), aspect='equal')
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color(EDGE_COLOR); s.set_linewidth(0.8)


def _scale_bar(ax, window):
    from matplotlib.patches import Rectangle

    x0, x1, y0, y1 = window
    span = x1 - x0
    nice = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
    length = max(v for v in nice if v <= span / 4.0)
    xs, ys = x0 + 0.05 * span, y0 + 0.05 * (y1 - y0)
    ax.add_patch(Rectangle((xs, ys), length, 0.014 * (y1 - y0), fc='black', ec='none', zorder=5))
    label = f'{length / 1000:g} mm' if length >= 1000 else f'{length:g} µm'
    ax.text(xs + length / 2.0, ys + 0.032 * (y1 - y0), label, ha='center', va='bottom', fontsize=9, zorder=5,
            bbox=dict(fc='white', ec='none', alpha=0.85, pad=1.0))


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
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    rev = design.metadata['geometry_revision']
    ym = float(design.metadata['config']['y_main_mm']) * 1e3
    fig_dir = DESIGN / 'fig'
    fig_dir.mkdir(parents=True, exist_ok=True)

    parts = {
        'launch_pad': {
            'label': 'Single input port LP_IN + feed_in + T1',
            'window': _fit_window([-3980, -2000, ym - 650, ym + 650]),
        },
        't_junction_stub': {
            'label': 'T2 + straight short stub S2 (10/6 µm CPW)',
            'window': _fit_window([-580, 180, ym - 260, ym + 1160]),
        },
        'open_end': {
            'label': 'Stage 3 open end R3_open (no port 2)',
            'window': _fit_window([2550, 3550, ym - 350, ym + 350]),
        },
    }

    w_c, h_c = chip_window[1] - chip_window[0], chip_window[3] - chip_window[2]
    fig, ax = plt.subplots(figsize=(8.8, 8.8 * h_c / w_c))
    _draw_metal(ax, metal, chip_window)
    for name, info in parts.items():
        x0, x1, y0, y1 = info['window']
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=BOX_COLOR, lw=1.4, zorder=4))
        ax.text(x0, y1 + 35, f' {name}', ha='left', va='bottom', fontsize=8, color=BOX_COLOR, weight='bold',
                zorder=4, bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.5))

    for rep in design.metadata['stubs']:
        ax.text(rep['x_mm'] * 1e3, rep['y_short_mm'] * 1e3 + 120,
                f"{rep['name']} ({rep['actual_filleted_mm']:.2f} mm)",
                ha='center', va='bottom', fontsize=8, color='#1b4f72', weight='bold',
                zorder=5, bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.5))
    for i, xk in enumerate(design.metadata['config']['stubs']['x_mm'], start=1):
        ax.text(xk * 1e3, ym - 220, f'T{i}', ha='center', va='top', fontsize=8, color='#1b4f72',
                weight='bold', zorder=5, bbox=dict(fc='white', ec='none', alpha=0.8, pad=0.5))
    for rep in design.metadata['resonators']:
        xm = 0.5 * (rep['x_leg_l_mm'] + rep['x_leg_r_mm']) * 1e3
        ax.text(xm, rep['y_loop_mm'] * 1e3 - 130,
                f"{rep['name']} ({rep['actual_filleted_mm']:.2f} mm)",
                ha='center', va='top', fontsize=8, color='#1b4f72', weight='bold',
                zorder=5, bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.5))
    ax.text(-3550, ym - 560, 'LP_IN (1 port only)', ha='center', va='top', fontsize=8.5, color='#1b4f72',
            weight='bold', zorder=5, bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.5))
    ax.text(3200, ym - 220, 'R3_open', ha='center', va='top', fontsize=8, color='#1b4f72',
            weight='bold', zorder=5, bbox=dict(fc='white', ec='none', alpha=0.85, pad=0.5))

    _scale_bar(ax, chip_window)
    fig.suptitle(f'{DESIGN.name} · {rev} · 1-Port 3-Stage Reflection BPF (8 × 6.8 mm, 10/6 µm CPW)',
                 x=0.5, ha='center', fontsize=9.4, weight='bold')
    fig.tight_layout()
    import io
    import os
    import time

    def _safe_write_bytes(path: Path, data: bytes):
        tmp = path.with_suffix(path.suffix + '.tmp')
        tmp.write_bytes(data)
        for attempt in range(5):
            try:
                os.replace(tmp, path)
                return
            except OSError:
                time.sleep(0.3)
        path.unlink(missing_ok=True)
        os.replace(tmp, path)

    chip_png = fig_dir / 'chip.png'
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=300, facecolor='white')
    plt.close(fig)
    _safe_write_bytes(chip_png, buf.getvalue())

    rows = [('chip.png', 'Whole chip (8.0 × 6.8 mm, 1 port only)', chip_window)]
    for part_name, info in parts.items():
        win = info['window']
        wp, hp = win[1] - win[0], win[3] - win[2]
        fig, ax = plt.subplots(figsize=(6.5, max(3.0, min(7.5, 6.5 * hp / wp))))
        _draw_metal(ax, metal, win)
        _scale_bar(ax, win)
        fig.suptitle(f"{DESIGN.name} · {rev} · {info['label']}", x=0.5, ha='center', fontsize=8.0, weight='bold')
        fig.tight_layout()
        png = fig_dir / f'{part_name}.png'
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=220, facecolor='white')
        plt.close(fig)
        _safe_write_bytes(png, buf.getvalue())
        rows.append((png.name, info['label'], win))

    lines = [
        f'# {DESIGN.name} figures',
        '',
        'Stage 1 draft layout figures (1-port 3-stage reflection BPF) rendered directly from QGeometry.',
        '**No GDS file has been generated yet.**',
        '',
        f'- Revision: `{rev}`',
        '',
        '| Figure | Shows | Window x0, x1, y0, y1 (µm) |',
        '|---|---|---|',
    ]
    for fname, label, win in rows:
        lines.append(f"| [{fname}]({fname}) | {label} | {', '.join(f'{v:g}' for v in win)} |")
    (fig_dir / 'README.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return [fig_dir / r[0] for r in rows]


def render_analytic_s_params(design):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    geo = design.metadata['config']
    z0 = float(geo['cpw']['assumed_z0_ohm'])
    v_ph = float(geo['cpw']['calibrated_v_ph_m_s'])
    ls = [float(v) * 1e-3 for v in geo['stubs']['lengths_mm']]
    lr = [float(v) * 1e-3 for v in geo['resonators']['lengths_mm']]

    freqs = np.linspace(2.0e9, 7.5e9, 2751)
    f_ghz = freqs * 1e-9

    def abcd_tl(l):
        th = 2.0 * np.pi * freqs * l / v_ph
        return np.array([[np.cos(th), 1j * z0 * np.sin(th)], [1j * np.sin(th) / z0, np.cos(th)]])

    def abcd_sc(l):
        th = 2.0 * np.pi * freqs * l / v_ph
        y = -1j / (z0 * np.tan(th))
        return np.array([[np.ones_like(th), np.zeros_like(th)], [y, np.ones_like(th)]])

    def mm(a, b):
        return np.array([
            [a[0, 0] * b[0, 0] + a[0, 1] * b[1, 0], a[0, 0] * b[0, 1] + a[0, 1] * b[1, 1]],
            [a[1, 0] * b[0, 0] + a[1, 1] * b[1, 0], a[1, 0] * b[0, 1] + a[1, 1] * b[1, 1]],
        ])

    m = abcd_sc(ls[0])
    m = mm(m, abcd_tl(lr[0]))
    m = mm(m, abcd_sc(ls[1]))
    m = mm(m, abcd_tl(lr[1]))
    m = mm(m, abcd_sc(ls[2]))
    m = mm(m, abcd_tl(lr[2]))
    a, c = m[0, 0], m[1, 0]
    h_db = 20.0 * np.log10(np.maximum(np.abs(2.0 / (a + z0 * c)), 1e-12))
    h_norm = h_db - np.max(h_db)
    zin = a / c
    s11 = (zin - z0) / (zin + z0)
    tau_ns = -np.gradient(np.unwrap(np.angle(s11)), 2.0 * np.pi * freqs) * 1e9

    pb = (f_ghz >= 5.5) & (f_ghz <= 6.5)
    sb = (f_ghz >= 2.0) & (f_ghz <= 4.0)
    idx3 = np.where(h_norm >= -3.0)[0]

    out_dir = DESIGN / 'notes' / 'fig'
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / 'analytic_s_params.png'

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8.6, 6.4), sharex=True)
    ax0.axvspan(2.0, 4.0, color='#fadbd8', alpha=0.45, label='Stopband target: 2–4 GHz (≤ -40 dB)')
    ax0.axvspan(5.5, 6.5, color='#d6eaf8', alpha=0.55, label='Passband target: 5.5–6.5 GHz')
    ax0.plot(f_ghz, h_norm, color='#1b4f72', lw=2.0, label='Normalized |V_open / V_inc| (50 kΩ probe S21)')
    ax0.axhline(-40.0, color='#922b21', lw=0.9, ls=':')
    ax0.axhline(-3.0, color='#555555', lw=0.8, ls=':')
    ax0.set_ylim(-75, 2)
    ax0.set_ylabel('Normalized Transfer (dB)')
    ax0.set_title(
        f'{DESIGN.name} · 1-Port 3-Stage Reflection BPF (Pre-EM Analytic ABCD, Unsimulated)\n'
        f'-3 dB: {f_ghz[idx3[0]]:.2f}–{f_ghz[idx3[-1]]:.2f} GHz | 5.5–6.5 GHz ripple: {h_norm[pb].min():.2f} dB | '
        f'2–4 GHz: ≤ {h_norm[sb].max():.1f} dB',
        fontsize=9.6,
        weight='bold',
    )
    ax0.grid(True, alpha=0.3)
    ax0.legend(loc='lower right', fontsize=8.3, frameon=True)

    ax1.axvspan(2.0, 4.0, color='#fadbd8', alpha=0.40)
    ax1.axvspan(5.5, 6.5, color='#d6eaf8', alpha=0.50)
    ax1.plot(f_ghz, tau_ns, color='#c0392b', lw=1.8, label='1-Port Reflection Group Delay τ_g(S11)')
    ax1.set_yscale('log')
    ax1.set_xlim(2.0, 7.5)
    ax1.set_ylim(0.01, 5.0)
    ax1.set_xlabel('Frequency (GHz)')
    ax1.set_ylabel('Group Delay τ_g(S11) (ns)')
    ax1.grid(True, alpha=0.3, which='both')
    ax1.legend(loc='upper left', fontsize=8.3, frameon=True)

    fig.tight_layout()
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=200, facecolor='white')
    plt.close(fig)
    out_png.write_bytes(buf.getvalue())
    return out_png, dict(
        f3_lo=float(f_ghz[idx3[0]]),
        f3_hi=float(f_ghz[idx3[-1]]),
        pb_min=float(h_norm[pb].min()),
        sb_max=float(h_norm[sb].max()),
    )


def verify_stage1_drc(design, metal):
    geo = design.metadata['config']
    cpw_outer_um = geo['cpw']['trace_width_um'] + 2.0 * geo['cpw']['trace_gap_um']

    print('\n=== Design10 (1-Port 3-Stage Reflection BPF) Stage 1 DRC Verification ===')
    for rep in design.metadata['stubs']:
        err = (rep['actual_filleted_mm'] - rep['target_len_mm']) * 1e3
        print(f"  [B8] {rep['name']} (straight shorted stub): target={rep['target_len_mm']:.4f} mm, "
              f"actual={rep['actual_filleted_mm']:.4f} mm, err={err:+.3f} um, y_short={rep['y_short_mm']:+.4f} mm")
        assert abs(err) < 0.5
    for rep in design.metadata['resonators']:
        err = (rep['actual_filleted_mm'] - rep['target_len_mm']) * 1e3
        print(f"  [B8] {rep['name']} (U-loop): target={rep['target_len_mm']:.4f} mm, "
              f"actual={rep['actual_filleted_mm']:.4f} mm, err={err:+.3f} um, y_loop={rep['y_loop_mm']:+.4f} mm")
        assert abs(err) < 0.5

    for cname in ('LP_IN', 'T1', 'T2', 'T3', 'S1_short', 'S2_short', 'S3_short', 'R3_open'):
        comp = design.components[cname]
        for pname, pinfo in comp.pins.items():
            assert np.linalg.norm(np.asarray(pinfo['normal'], float)) > 0.99
    print('  [B13] All 15 component pins verified with unit outward normals (1 physical port LP_IN only)')

    min_loop_w = min(r['loop_width_mm'] * 1e3 - cpw_outer_um for r in design.metadata['resonators'])
    r1_r2_gap = (design.metadata['resonators'][1]['x_leg_l_mm'] - design.metadata['resonators'][0]['x_leg_r_mm']) * 1e3 - cpw_outer_um
    r2_r3_gap = (design.metadata['resonators'][2]['x_leg_l_mm'] - design.metadata['resonators'][1]['x_leg_r_mm']) * 1e3 - cpw_outer_um
    print(f'  [B7] Trench-to-trench clearances: min_loop_width={min_loop_w:.0f} um, '
          f'R1-R2={r1_r2_gap:.0f} um, R2-R3={r2_r3_gap:.0f} um (all >= 300 um)')
    assert min(min_loop_w, r1_r2_gap, r2_r3_gap) >= 300.0

    sx_half, sy_half = geo['chip']['size_x_mm'] / 2.0, geo['chip']['size_y_mm'] / 2.0
    top_m = (sy_half - max(r['y_short_mm'] for r in design.metadata['stubs'])) * 1e3 - cpw_outer_um / 2.0
    bot_m = (sy_half + min(r['y_loop_mm'] for r in design.metadata['resonators'])) * 1e3 - cpw_outer_um / 2.0
    right_m = (sx_half - float(geo['resonators']['x_open_mm'])) * 1e3 - cpw_outer_um / 2.0
    print(f'  [B11] Edge keep-out: top={top_m:.0f} um, bottom={bot_m:.0f} um, right={right_m:.0f} um (>= 500 um)')
    assert min(top_m, bot_m, right_m) >= 500.0

    pt = design.qgeometry.tables['path']
    pos_paths = pt[~pt['subtract']]
    lines = [_fillet_linestring(r['geometry'], float(r['fillet']) if r['fillet'] else 0.0)
             for _, r in pos_paths.iterrows()]
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            assert lines[i].intersection(lines[j]).is_empty
    print(f'  [B5/B6] Route crossings across {len(lines)} CPW routes: 0')

    thin = metal.difference(metal.buffer(-2.0, join_style=2).buffer(+2.0, join_style=2))
    thin_parts = [g for g in getattr(thin, 'geoms', [thin]) if g.area > 5.0]
    print(f'  [C12a] Minimum metal width (>= 4 um) check: {len(thin_parts)} violation(s)')
    assert len(thin_parts) == 0

    n_polys = len(metal.geoms) if isinstance(metal, MultiPolygon) else 1
    print(f'  [B9] Connected metal polygon count: {n_polys} (0 isolated ground islands)')
    assert n_polys == 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--export-sim-gds', action='store_true')
    args = ap.parse_args()

    design = build_design()
    metal, chip_window = design_to_metal_um(design)
    verify_stage1_drc(design, metal)

    written = render_layout_figures(design, metal, chip_window)
    print('\nRendered layout figures (NO GDS written):')
    for p in written:
        print(f'  {p}')

    s_png, s_sum = render_analytic_s_params(design)
    print(f'\nRendered pre-EM analytic response: {s_png}')
    print(f"  1-port 3-stage ABCD: -3 dB = {s_sum['f3_lo']:.3f}-{s_sum['f3_hi']:.3f} GHz | "
          f"5.5-6.5G min={s_sum['pb_min']:.2f} dB | 2-4G max={s_sum['sb_max']:.2f} dB")

    if args.export_sim_gds:
        from gds_tools import export_gds
        rev = design.metadata['geometry_revision']
        gds_dir = DESIGN / 'gds'
        gds_dir.mkdir(parents=True, exist_ok=True)
        sim_gds = gds_dir / f'{DESIGN.name}_{rev}_sim.gds'
        export_gds(design, str(sim_gds), fab=False)
        print(f'Wrote simulation companion GDS: {sim_gds}')


if __name__ == '__main__':
    main()
