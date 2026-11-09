"""Design5 chip: 6 readout channels (lambda/4 resonator + individual lambda/4 Purcell filter, W1 / Wallraff-type)
on one feed line, 10 mm x 10 mm chip, input and output pads in the middle of the left and right edges, facing
each other (y = 0): one straight feed line, no bends.

    Pin (-4.5, 0) -[C_in]- n0 ch2 - n0 ch4 - n0 ch6 - n0 ch1 - n0 ch3 - n0 ch5 - Pout (4.5, 0)
         up (mirrored):        ch4 6.9          ch1 6.3          ch5 7.1
         down:         ch2 6.5          ch6 7.3          ch3 6.7

  rows staggered by half the 2.6 mm pitch (T junctions 1.3 mm apart); frequency order chosen by brute force:
  neighbours in the same row differ >= 600 MHz, neighbouring T junctions >= 400 MHz.
  Serpentines: short run next to the junction, 2 equal legs running away from the feed line -> 1.5 mm wide channels.
  NOTE: the pads are not at the 8-port PCB snapshot positions (side ports there are at y = +-1.6 mm).

Each channel = the Design5 W1 channel (design5_layout.py): C_p hugger (plate 165 um, frame on the filter side),
C_j hammer (gap 19 um), taps L_pf 1.2 mm / L_rr 1.3 mm; lengths per channel from simulation/solve_chip6.py
(results/chip6_lengths.json). One input capacitor C_in (enclosing cap, as Wallraff: ~50 fF at the feed-line input)
on the vertical run after P1, no per-channel C_in. Filter / resonator are drawn as deterministic serpentines
(exact length incl. fillets) because RouteMixed meanders fail or miss the length for the short high-f lines.

Usage: python chip6_layout.py   -> gds/chip6.gds, fig/chip6_layout.png
"""
import os, sys, json, warnings
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import design5_layout as d5                                       # constants + component classes (W1)
from design5_layout import (designs, Dict, LaunchpadWirebond, RouteStraight, ResEnd,
                            TapJunction, EnclosingCap, HammerCap, pcb_ports, f3, pl)

LEN_JSON = os.path.join(ROOT, 'simulation', 'results', 'chip6_lengths.json')
GDS_FILE = os.path.join(ROOT, 'gds', 'chip6.gds')
# ports: middle of the left and right chip edges, facing each other (NOT the 8-port PCB snapshot positions,
# whose side ports sit at y = +-1.6 mm) -> one straight feed line at y = 0, no bends
PORTS = {'Pin': dict(pos_x='-4.5mm', pos_y='0mm', orientation='0'), 'Pout': dict(pos_x='4.5mm', pos_y='0mm', orientation='180')}
Y_FL = 0.0                               # feed line
X_CIN = -3.9                             # C_in centre, in line right after the input pad
# serpentines: short horizontal run next to the junction, legs running AWAY from the feed line, leg height chosen
# close to a target (same look for all channels) within [HMIN, HMAX]
XRUN_F, SP_F, HT_F, HMAX_F = 0.30, 0.25, 0.90, 1.40     # filter: run, leg pitch, target / max leg height (mm)
XRUN_R, SP_R, HT_R, HMAX_R = 0.40, 0.35, 1.20, 1.80     # resonator (pitch = MSP_R 350 um)
HMIN, N_LEGS = 0.40, 2                  # min leg height (mm); preferred number of legs (all channels alike)
# (channel, n0 x position, up?) in feed-line order (left -> right); rows staggered by half the 2.6 mm pitch
ORDER = [(2, -3.20, False), (4, -1.90, True), (6, -0.60, False), (1, 0.70, True), (3, 2.00, False), (5, 3.30, True)]


def add_channel(design, ch, x0, y0, Lf, Lr, cpw, up=False):
    """One W1 channel hanging below (up=False) or mirrored above (up=True) the feed line at n0 = (x0, y0)."""
    s = f'_ch{ch}'; sg = 1 if up else -1                           # sg: direction away from the feed line
    straight = lambda name, a, b: RouteStraight(design, name, Dict(chip='main', pin_inputs=Dict(
        start_pin=Dict(component=a[0], pin=a[1]), end_pin=Dict(component=b[0], pin=b[1]))))
    E, EP = d5.ENC, d5.ENC_P; e = d5.enc_xpin() / 1000
    TapJunction(design, 'n0' + s, options=dict(pos_x=f'{x0}mm', pos_y=f'{y0}mm', orientation='0', prime_length='0.3mm',
                                               stub_length='0.1mm', stub_up=str(up), **cpw))
    y = y0 + sg * (0.15 + e)
    o_cp = (-90 if d5.CP_TYPE == 'T' else 90) * (-1 if up else 1)
    EnclosingCap(design, 'cp' + s, options=dict(
        chip='main', pos_x=f'{x0}mm', pos_y=f'{y:.5f}mm', orientation=str(o_cp),
        plate_width=f"{EP['plate_width']}um", plate_length=f"{E['plate_length']}um", gap=f"{E['gap']}um",
        wall=f"{E['wall']}um", wall_side=f"{E['wall_side']}um", ground_gap=f"{E['ground_gap']}um",
        leg_ext=f"{E['leg_ext']}um", lead_length=f"{E['lead_length']}um", **cpw))
    fpin = 'b' if d5.CP_TYPE == 'T' else 'a'
    straight('stub_p' + s, ('n0' + s, 'second_end'), ('cp' + s, 'a' if d5.CP_TYPE == 'T' else 'b'))
    y_c = y + sg * e; y_j = y_c + sg * d5.D_CJ
    DX, prime, half = 0.20, 0.3, 0.15
    stub_len = (DX - 2 * d5.CJ_TAPER_X / 1000 - 2 * d5.CJ_GAP) / 2
    o_j = '-90' if up else '90'                                    # prime_end always toward C_p / the feed line
    TapJunction(design, 'jf' + s, options=dict(pos_x=f'{x0}mm', pos_y=f'{y_j:.5f}mm', orientation=o_j, prime_length=f'{prime}mm',
                                               stub_length=f'{stub_len:.5f}mm', stub_up=str(up), **cpw))
    straight('feed_jf' + s, ('cp' + s, fpin), ('jf' + s, 'prime_end'))
    TapJunction(design, 'jr' + s, options=dict(pos_x=f'{x0 + DX}mm', pos_y=f'{y_j:.5f}mm', orientation=o_j, prime_length=f'{prime}mm',
                                               stub_length=f'{stub_len:.5f}mm', stub_up=str(not up), **cpw))
    HammerCap(design, 'cj' + s, options=dict(chip='main', pos_x=f'{x0 + DX / 2:.5f}mm', pos_y=f'{y_j:.5f}mm', orientation='0',
                                             head_width=f'{d5.CJ_HW}um', head_height=f'{d5.CJ_HH}um', head_gap=f'{d5.CJ_HG}um',
                                             lead_length=f'{d5.CJ_LEAD}um', ground_gap=f'{d5.CJ_GG}um', **cpw))
    straight('stub_f' + s, ('jf' + s, 'second_end'), ('cj' + s, 'a'))
    straight('stub_r' + s, ('jr' + s, 'second_end'), ('cj' + s, 'b'))
    # filter: from jf away from the feed line to VERT, short run left, then an equal-leg serpentine (away first), short
    snake_route(design, 'filt' + s, 'end_filter' + s, ('jf' + s, 'prime_start'), (x0, y_j + sg * half), sg,
                d5.VERT - d5.D_CJ - half, -1, XRUN_F, sg, SP_F, Lf - (d5.D_CJ - half) - prime, HT_F, HMAX_F)
    # resonator: open end beyond the junction, straight res_a to jr
    y_ro = y_c + sg * (d5.VERT + d5.OPEN_DY)
    ResEnd(design, 'open_res' + s, options=dict(pos_x=f'{x0 + DX}mm', pos_y=f'{y_ro:.5f}mm', orientation='270' if up else '90'))
    straight('res_a' + s, ('open_res' + s, 'end'), ('jr' + s, 'prime_start'))
    # resonator: from jr toward the feed line to 0.44 mm from C_p pin, right, then an equal-leg serpentine away, short
    snake_route(design, 'res' + s, 'end_res' + s, ('jr' + s, 'prime_end'), (x0 + DX, y_j - sg * half), -sg,
                d5.D_CJ - 0.44 - half, +1, XRUN_R, sg, SP_R, Lr - (d5.VERT + d5.OPEN_DY - d5.D_CJ - half) - prime, HT_R, HMAX_R)


def snake_points(A0, sv, s1, sh, xrun, sl, sp, L, htarget, hmax, r):
    """Deterministic serpentine from point A0: vertical s1 (direction sv), horizontal xrun (direction sh), then n + 1
    equal vertical legs of height h (first direction sl, alternating) separated by horizontal steps sp.
    Exact length with fillets of radius r: every 90 deg corner shortens the path by c = r (2 - pi / 2):
        h(n) = (L - s1 - xrun - n sp + (2n + 2) c) / (n + 1)
    n is chosen so that h is closest to htarget within [HMIN, hmax] (uniform look, narrow channel)."""
    c = r * (2 - np.pi / 2)
    hs = [(n, (L - s1 - xrun - n * sp + (2 * n + 2) * c) / (n + 1)) for n in range(0, 40)]
    ok = [(n, h) for n, h in hs if HMIN <= h <= hmax]
    pref = [(n, h) for n, h in ok if n + 1 == N_LEGS]
    if pref:                                         # same number of legs on every channel when possible
        ok = pref
    if not ok:
        ok = [min(hs, key=lambda t: abs(t[1] - htarget))]
    n, h = min(ok, key=lambda t: abs(t[1] - htarget))
    if h < 2.2 * r:
        raise ValueError(f'serpentine leg {h:.3f} mm too short for the fillet (L = {L:.3f} mm)')
    x, y = A0
    pts = [(x, y), (x, y + sv * s1), (x + sh * xrun, y + sv * s1)]
    x, y = pts[-1]; d = sl
    for i in range(n + 1):
        y += d * h; pts.append((x, y))
        if i < n:
            x += sh * sp; pts.append((x, y)); d = -d
    return pts, n, h, ('90' if d < 0 else '270')        # short end: the route leaves it opposite to the last leg


def snake_route(design, name, end_name, pin, A0, sv, s1, sh, xrun, sl, sp, L, htarget, hmax):
    """Shorted end (ResEnd) at the serpentine end, RouteAnchors from it back to `pin` (jf / jr)."""
    r = float(str(f3.FILLET).replace('um', '')) / 1000
    pts, n, h, orient = snake_points(A0, sv, s1, sh, xrun, sl, sp, L, htarget, hmax, r)
    SNAKES[name] = dict(legs=n + 1, h=round(h, 3))
    E = pts[-1]
    ResEnd(design, end_name, options=dict(pos_x=f'{E[0]:.5f}mm', pos_y=f'{E[1]:.5f}mm', orientation=orient))
    f3.route(design, name, (end_name, 'end'), pin, list(reversed(pts[1:-1])))
    return n, h


def build_chip():
    warnings.filterwarnings('ignore')
    L = {c['ch']: c for c in json.load(open(LEN_JSON))['channels']}
    design = designs.DesignPlanar(); design.overwrite_enabled = True
    design.chips.main.material = 'silicon'
    design.variables['cpw_width'] = f3.CPW_WIDTH; design.variables['cpw_gap'] = f3.CPW_GAP
    design.chips.main.size.size_x = '10mm'; design.chips.main.size.size_y = '10mm'; design.chips.main.size.size_z = '-525um'
    cpw = dict(trace_width=f3.CPW_WIDTH, trace_gap=f3.CPW_GAP)
    for name, opt in PORTS.items():
        LaunchpadWirebond(design, 'LP' + name, options=dict(f3.LP, **opt))
    # input capacitor in line right after the input pad (frame / pin a toward the input pad)
    E = d5.ENC
    EnclosingCap(design, 'cin', options=dict(
        chip='main', pos_x=f'{X_CIN}mm', pos_y=f'{Y_FL}mm', orientation='0',
        plate_width=f"{E['plate_width']}um", plate_length=f"{E['plate_length']}um", gap=f"{E['gap']}um",
        wall=f"{E['wall']}um", wall_side=f"{E['wall_side']}um", ground_gap=f"{E['ground_gap']}um", leg_ext=f"{E['leg_ext']}um",
        lead_length=f"{E['lead_length']}um", **cpw))
    SNAKES.clear()
    for ch, x0, up in ORDER:
        add_channel(design, ch, x0, Y_FL, L[ch]['L_filter'], L[ch]['L_res'], cpw, up)
    n = lambda k: f'n0_ch{ORDER[k][0]}'
    st = lambda name, a, b: RouteStraight(design, name, Dict(chip='main', pin_inputs=Dict(
        start_pin=Dict(component=a[0], pin=a[1]), end_pin=Dict(component=b[0], pin=b[1]))))
    st('feed_in', ('LPPin', 'tie'), ('cin', 'a'))
    st('fl_0', ('cin', 'b'), (n(0), 'prime_start'))
    for k in range(1, len(ORDER)):
        st(f'fl_{k}', (n(k - 1), 'prime_end'), (n(k), 'prime_start'))
    st(f'fl_{len(ORDER)}', (n(len(ORDER) - 1), 'prime_end'), ('LPPout', 'tie'))
    return design


SNAKES = {}


def channel_geoms(design):
    """Union of every channel's metal (n0 excluded), for spacing checks."""
    from shapely.ops import unary_union
    id2n = {c.id: nm for nm, c in design.components.items()}
    out = {}
    for ch, _, _ in ORDER:
        geoms = []
        for kind in ('path', 'poly'):
            t = design.qgeometry.tables[kind]
            for r in t.itertuples():
                nm = id2n[int(r.component)]
                if nm.endswith(f'_ch{ch}') and not nm.startswith(('n0', 'stub_p')) and not bool(r.subtract):
                    geoms.append(r.geometry.buffer(float(r.width) / 2) if kind == 'path' else r.geometry)
        out[ch] = unary_union(geoms)
    return out


if __name__ == '__main__':
    import geopandas as gpd
    from shapely.ops import unary_union
    if not hasattr(gpd.GeoDataFrame, 'append'):
        gpd.GeoDataFrame.append = gpd.GeoDataFrame._append
    d = build_chip()
    G = channel_geoms(d)
    for ch, g in G.items():
        print(f'ch{ch}: x {g.bounds[0]:.3f} .. {g.bounds[2]:.3f}, y {g.bounds[1]:.3f} .. {g.bounds[3]:.3f}')
    chs = list(G)
    gaps = {(a, b): G[a].distance(G[b]) for i, a in enumerate(chs) for b in chs[i + 1:]}
    print('min channel-channel metal distance (mm):', round(min(gaps.values()), 3), 'between ch', min(gaps, key=gaps.get))
    print('serpentines (legs, leg height mm):', SNAKES)
    id2n = {c.id: nm for nm, c in d.components.items()}
    pads = unary_union([r.geometry for r in d.qgeometry.tables['poly'].itertuples()
                        if id2n[int(r.component)] in ('LPPin', 'LPPout', 'cin') and not bool(r.subtract)])
    print('min pad / C_in -> channel distance (mm):', round(min(pads.distance(g) for g in G.values()), 3))
    allg = unary_union(list(G.values()))
    b = allg.bounds
    print(f'chip extent of channels: x {b[0]:.3f} .. {b[2]:.3f}, y {b[1]:.3f} .. {b[3]:.3f} (chip +-5)')
    os.makedirs(os.path.dirname(GDS_FILE), exist_ok=True)
    pl.export_gds(d, GDS_FILE)
    png = os.path.join(ROOT, 'fig', 'chip6_layout.png')
    pl.plot_gds(GDS_FILE, png, window=(-5, 5, -5, 5))
    print('gds', GDS_FILE, 'png', png)
