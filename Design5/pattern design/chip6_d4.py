"""Design5 chip, Design4 style (Jay 2026-11-12): 6 readout channels, 6.0-7.0 GHz, each channel = the Design4 Cai 2024
channel (Design4/pattern design/purcell_6500.py): CPW 16/9 um, Si eps 11.5 / 500 um, T-type enclosing C_p (plate 182 um,
gap 12 um, ~40 fF), hammer C_j 25 x 65 um plates, gap 42 um (~2 fF), f_p = f_r, wide serpentines (every CPW >= 400 um apart),
qiskit auto wire bonds every 350 um + hammer / tee straps.

Same chip arrangement as chip6_layout.py: 10 x 10 mm, pads in the middle of the left / right edges, one straight feed
line at y = 0, one C_in (Cai enclosing cap) after the input pad, 3 channels up / 3 mirrored down, staggered:

    Pin -[C_in]- n0 ch2 - n0 ch4 - n0 ch6 - n0 ch1 - n0 ch3 - n0 ch5 - Pout
          up:            ch4 6.6          ch1 6.0          ch5 6.8
          down:  ch2 6.2          ch6 7.0          ch3 6.4

Lengths: purcell_6500.tune_lengths at f_model = f * HFSS_CORR (Design4 eigenmode: model 6.08 GHz -> HFSS 5.975 GHz,
both modes, 2026-11-12), saved in parameters/Design5_<date>.json.

Usage: python chip6_d4.py [--retune]   -> gds/Design5_<date>_nocheese.gds, fig/Design5_<date>_layout.png
(fab GDS with cheese holes: chip6_d4_fab.py -> gds/Design5_<date>.gds)
"""
import json, os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(os.path.dirname(ROOT), 'Design4', 'pattern design')]
import purcell_6500 as P                     # sets SI_EPS 11.5 before design5_layout is imported
import design5_layout as d5
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround

F_CHIP = [6.0, 6.2, 6.4, 6.6, 6.8, 7.0]      # GHz, ch1..ch6
HFSS_CORR = 6.08 / 5.975                     # Design4 eigenmode vs circuit model
ORDER = [(2, -3.20, False), (4, -1.90, True), (6, -0.60, False), (1, 0.70, True), (3, 2.00, False), (5, 3.30, True)]
PORT_X, X_CIN = 4.5, -3.9
# file naming (Jay 2026-11-13): <folder>_<date>, e.g. gds/Design5_2026-11-13.gds, parameters/Design5_2026-11-13.json
import glob, time
TAG = f"{os.path.basename(ROOT)}_{time.strftime('%Y-%m-%d')}"
PAR_DIR = os.path.join(ROOT, 'parameters')
PAR_FILE = os.path.join(PAR_DIR, TAG + '.json')                          # written today
_par = sorted(glob.glob(os.path.join(PAR_DIR, os.path.basename(ROOT) + '_*.json')))
LEN_JSON = _par[-1] if _par else PAR_FILE                                 # newest parameter file
GDS_FILE = os.path.join(ROOT, 'gds', TAG + '_nocheese.gds')
CAP_COMPS = ('cp', 'cj', 'cin')


def tune_all():
    out = []
    for k, f in enumerate(F_CHIP, 1):
        P.F0 = f * HFSS_CORR
        p, r = P.tune_lengths()
        out.append(dict(ch=k, f=f, f_model=P.F0, L_filter=float(p['L_filt']), L_res=float(p['L_res']),
                        J_model=float(r['J']), kp_model=float(r['kp'])))
        print(f"ch{k} {f:.2f} GHz: L_filter {p['L_filt']:.4f} L_res {p['L_res']:.4f} mm, J {r['J']:.1f} kp {r['kp']:.1f} MHz")
    save_params(out)
    return out


def save_params(channels, extra=None):
    """Update parameters/<folder>_<date>.json (sections: substrate, chip, cpw, C_p, C_j_hammer, channels, fab ...).
    Only `channels` (and `extra`) are rewritten; the other sections are kept from the newest parameter file."""
    os.makedirs(PAR_DIR, exist_ok=True)
    src = PAR_FILE if os.path.exists(PAR_FILE) else LEN_JSON
    par = json.load(open(src)) if os.path.exists(src) else {'design': os.path.basename(ROOT)}
    old = {c['ch']: c for c in par.get('channels', [])}
    par['date'] = time.strftime('%Y-%m-%d')
    par['channels'] = [dict(old.get(c['ch'], {}), ch=c['ch'], f_target_GHz=c['f'], f_model_GHz=round(c['f_model'], 4),
                            L_filter_mm=round(c['L_filter'], 4), L_res_mm=round(c['L_res'], 4),
                            J_model_MHz=round(c['J_model'], 2), kappa_p_model_MHz=round(c['kp_model'], 2)) for c in channels]
    par.update(extra or {})
    json.dump(par, open(PAR_FILE, 'w'), indent=1, default=str)
    return PAR_FILE


def add_channel(design, ch, x0, up, Lf, Lr):
    """Design4 channel (purcell_6500.build, wide) at n0 = (x0, 0); up=False mirrors it below the feed line."""
    s, sg = f'_ch{ch}', 1 if up else -1
    pos = lambda x, y: P.position(x0 + x, sg * y)
    d5.TapJunction(design, 'n0' + s, options=dict(**P.CPW, **pos(0, 0), orientation='0', stub_length='100um', stub_up=str(up)))
    opts = {k: f'{v}um' for k, v in dict(d5.ENC, **P.ENC_CAI).items()}
    d5.EnclosingCap(design, 'cp' + s, options=dict(opts, **P.CPW, **pos(0, P.CP_Y), orientation='90' if up else '-90'))
    P.straight(design, 'stub_p' + s, ('n0' + s, 'second_end'), ('cp' + s, 'a'))
    H = P.PAPER_HAMMER
    spacing = P.MIN_CPW_SPACING          # Design4 final: filter-resonator 400 um (W400g42)
    stub = spacing / 2 - (H['gap_um'] / 2 + H['width_um'] + d5.CJ_LEAD) / 1000 - 0.02
    for name, x, su in [('jf', 0, up), ('jr', spacing, not up)]:
        d5.TapJunction(design, name + s, options=dict(**P.CPW, **pos(x, P.TAP_Y), orientation='-90' if up else '90',
                                                      stub_up=str(su), stub_length=f'{stub:.9f}mm'))
    d5.HammerCap(design, 'cj' + s, options=dict(**P.CPW, **pos(spacing / 2, P.TAP_Y), head_width=f"{H['width_um']}um",
                                                head_height=f"{H['height_um']}um", head_gap=f"{H['gap_um']}um",
                                                lead_length=f'{d5.CJ_LEAD}um', ground_gap=f"{H['ground_um']}um"))
    P.straight(design, 'stub_f' + s, ('jf' + s, 'second_end'), ('cj' + s, 'a'))
    P.straight(design, 'stub_r' + s, ('cj' + s, 'b'), ('jr' + s, 'second_end'))
    P.straight(design, 'feed_jf' + s, ('cp' + s, 'b'), ('jf' + s, 'prime_end'))
    OpenToGround(design, 'readout_open' + s, options=dict(**pos(spacing, P.TAP_Y + P.RES_TAP), orientation='90' if up else '270',
                                                          width=P.CPW['trace_width'], gap=P.CPW['trace_gap']))
    P.straight(design, 'res_open_run' + s, ('jr' + s, 'prime_start'), ('readout_open' + s, 'open'))
    W = P.WIDE
    f_lead = W['filter']['lead']
    P.serpentine_tail(design, 'filter_tail' + s, ('jf' + s, 'prime_start'), Lf - P.FILTER_TAP - P.TEE_HALF, f_lead, -1, sg,
                      height=P.TAP_Y + P.TEE_HALF + f_lead - W['y_bottom'])
    P.serpentine_tail(design, 'res_tail' + s, ('jr' + s, 'prime_end'), Lr - P.RES_TAP - P.TEE_HALF,
                      (P.TAP_Y - P.TEE_HALF) - W['y_bottom'], 1, -sg, height=P.TAP_Y + P.RES_TAP - W['y_bottom'])
    for name, target, fixed in [('filter_tail' + s, Lf, P.FILTER_TAP), ('res_tail' + s, Lr, P.RES_TAP)]:
        actual = d5.f3.path_length(design, name) + fixed + P.TEE_HALF
        if abs(actual - target) > 0.001:
            raise ValueError(f'{name}: drawn {actual:.5f} mm vs {target:.5f} mm')


def build_chip(lengths=None):
    if lengths is None:
        lengths = json.load(open(LEN_JSON))['channels']
    L = {c['ch']: dict(L_filter=c.get('L_filter', c.get('L_filter_mm')), L_res=c.get('L_res', c.get('L_res_mm')))
         for c in lengths}
    design = d5.designs.DesignPlanar(); design.overwrite_enabled = True
    design.chips.main.material = d5.SI_MATERIAL
    design.chips.main.size.update(dict(size_x='10mm', size_y='10mm', size_z='-500um', center_x='0mm', center_y='0mm'))
    design.variables.update(cpw_width=P.CPW['trace_width'], cpw_gap=P.CPW['trace_gap'])
    for name, x, o in [('LPPin', -PORT_X, '0'), ('LPPout', PORT_X, '180')]:
        d5.LaunchpadWirebond(design, name, options=dict(d5.f3.LP, **P.position(x, 0), orientation=o))
    P.enclosing(design, 'cin', X_CIN, 0, 0)
    for ch, x0, up in ORDER:
        add_channel(design, ch, x0, up, L[ch]['L_filter'], L[ch]['L_res'])
    n = lambda k: f'n0_ch{ORDER[k][0]}'
    P.straight(design, 'feed_in', ('LPPin', 'tie'), ('cin', 'a'))
    P.straight(design, 'fl_0', ('cin', 'b'), (n(0), 'prime_start'))
    for k in range(1, len(ORDER)):
        P.straight(design, f'fl_{k}', (n(k - 1), 'prime_end'), (n(k), 'prime_start'))
    P.straight(design, f'fl_{len(ORDER)}', (n(len(ORDER) - 1), 'prime_end'), ('LPPout', 'tie'))
    return design


def channel_geoms(design):
    from shapely.ops import unary_union
    id2n = {c.id: nm for nm, c in design.components.items()}
    out = {}
    for ch, _, _ in ORDER:
        g = [r.geometry.buffer(float(r.width) / 2) if k == 'path' else r.geometry
             for k in ('path', 'poly') for r in design.qgeometry.tables[k].itertuples()
             if id2n[int(r.component)].endswith(f'_ch{ch}') and not id2n[int(r.component)].startswith(('n0', 'stub_p'))
             and not bool(r.subtract)]
        out[ch] = unary_union(g)
    return out


if __name__ == '__main__':
    import geopandas as gpd
    if not hasattr(gpd.GeoDataFrame, 'append'):
        gpd.GeoDataFrame.append = gpd.GeoDataFrame._append
    lengths = tune_all() if ('--retune' in sys.argv or not os.path.exists(LEN_JSON)) else None
    d = build_chip(lengths)
    G = channel_geoms(d)
    for ch, g in G.items():
        b = g.bounds; print(f'ch{ch}: x {b[0]:.3f} .. {b[2]:.3f}, y {b[1]:.3f} .. {b[3]:.3f}')
    chs = list(G)
    gaps = {(a, b): G[a].distance(G[b]) for i, a in enumerate(chs) for b in chs[i + 1:]}
    print('min channel-channel metal distance (mm):', round(min(gaps.values()), 3), 'between ch', min(gaps, key=gaps.get))
    try:
        import cpw_spacing_check as S
        S.report(d, 400)
    except Exception as e:
        print('spacing check skipped:', repr(e)[:120])
    os.makedirs(os.path.dirname(GDS_FILE), exist_ok=True)
    d5.pl.export_gds(d, GDS_FILE)
    png = os.path.join(ROOT, 'fig', TAG + '_layout.png'); os.makedirs(os.path.dirname(png), exist_ok=True)
    d5.pl.plot_gds(GDS_FILE, png, window=(-5, 5, -5, 5))
    print('gds', GDS_FILE, 'png', png)
