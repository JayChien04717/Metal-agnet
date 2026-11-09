"""Design6: verification layouts - do the capacitors predict the couplings that HFSS S21 shows?

One W1 channel (Design5 / Wallraff type, 6.7 GHz lengths of chip6 ch3), hairpin routing, NO C_in: the feed line runs
straight from the left pad to the right pad (matched on both sides), so S21 = 2 / (2 + Z0 Y_branch) up to a phase.

  VARIANT=filter : feed line - C_p - Purcell filter (shorted). Resonator lines removed, the C_j hammer kept (its
                   resonator plate floating) so the filter is loaded exactly as in the full channel.
                   -> ONE notch: omega_p and kappa_p (Purcell <-> feed line, set by C_p)
  VARIANT=full   : the whole channel: + C_j + resonator (open end, shorted end)
                   -> two hybridised notches: J (resonator <-> Purcell, set by C_j)

Usage: python verify_layout.py [filter|full|both]   -> gds/verify_<variant>.gds, fig/verify_<variant>.png (+ zoom)
Parameters for the C sweeps: env CP_PLATE (um, hugger plate width, default 165) and CJ_GAP (um, hammer gap, default 19).
"""
import os, sys, json, warnings
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); HFSS = os.path.dirname(ROOT)
os.environ.setdefault('CHIP6_STYLE', 'hairpin')
sys.path.insert(0, os.path.join(HFSS, 'Design5', 'pattern design'))
import design5_layout as d5                                       # W1 constants (caps, taps), shared components
CP_PLATE = float(os.environ.get('CP_PLATE', d5.ENC_P['plate_width']))
CJ_GAP_UM = float(os.environ.get('CJ_GAP', d5.CJ_HG))
d5.ENC_P = dict(d5.ENC_P, plate_width=CP_PLATE)
d5.CJ_HG = CJ_GAP_UM; d5.CJ_TAPER_X = d5.CJ_HG / 2 + d5.CJ_HW + d5.CJ_LEAD
import chip6_layout as C6                                         # add_channel (hairpin), PORTS
from design5_layout import designs, Dict, LaunchpadWirebond, RouteStraight, f3, pl

LEN_JSON = os.path.join(HFSS, 'Design5', 'simulation', 'results', 'chip6_lengths.json')
CH = 3                                                            # 6.7 GHz lengths
X0, Y_FL = 0.0, 0.0                                               # n0 at the chip centre, channel below the line
# short feed line (Jay 2026-11-09): pads 1.2 mm from n0 instead of at the 10 mm chip edges; small chip around the
# channel (x -2..2, y +1.0..-3.6 mm) -> smaller HFSS box
PORTS = {'Pin': dict(pos_x='-1.2mm', pos_y='0mm', orientation='0'), 'Pout': dict(pos_x='1.2mm', pos_y='0mm', orientation='180')}
CHIP = dict(size_x='4.0mm', size_y='4.6mm', center_x='0mm', center_y='-1.3mm')
WIREBONDS = os.environ.get('WIREBONDS', '1') == '1'               # HFSS wire bonds across every CPW (ground straps)


def add_wirebonds(design):
    """Flag the CPW gap paths for HFSS wire bonds. This qiskit-metal version ignores the route option
    'hfss_wire_bonds'; the HFSS renderer reads the column of the path table instead (ansys_renderer.auto_wirebonds:
    one bond every wb_threshold = 400 um on straight segments longer than that). Call after the last rebuild."""
    t = design.qgeometry.tables['path']
    m = t['subtract'].astype(bool)
    t.loc[m, 'hfss_wire_bonds'] = True
    return int(m.sum())
RES_PARTS = ['stub_r', 'res', 'end_res', 'res_a', 'open_res', 'jr']   # removed in VARIANT=filter (routes first)


def lengths():
    c = {c['ch']: c for c in json.load(open(LEN_JSON))['channels']}[CH]
    return c['L_filter'], c['L_res']


def build(variant='full'):
    warnings.filterwarnings('ignore')
    design = designs.DesignPlanar(); design.overwrite_enabled = True
    design.chips.main.material = 'silicon'
    design.variables['cpw_width'] = f3.CPW_WIDTH; design.variables['cpw_gap'] = f3.CPW_GAP
    for k, v in CHIP.items():
        design.chips.main.size[k] = v
    design.chips.main.size.size_z = '-525um'
    cpw = dict(trace_width=f3.CPW_WIDTH, trace_gap=f3.CPW_GAP)
    for name, opt in PORTS.items():
        LaunchpadWirebond(design, 'LP' + name, options=dict(f3.LP, **opt))
    Lf, Lr = lengths()
    C6.SNAKES.clear()
    C6.add_channel(design, CH, X0, Y_FL, Lf, Lr, cpw, up=False)
    st = lambda name, a, b: RouteStraight(design, name, Dict(chip='main', pin_inputs=Dict(
        start_pin=Dict(component=a[0], pin=a[1]), end_pin=Dict(component=b[0], pin=b[1]))))
    st('feed_in', ('LPPin', 'tie'), (f'n0_ch{CH}', 'prime_start'))      # no C_in: straight, matched feed line
    st('feed_out', (f'n0_ch{CH}', 'prime_end'), ('LPPout', 'tie'))
    if variant == 'filter':
        for p in RES_PARTS:
            design.delete_component(f'{p}_ch{CH}', force=True)
    if WIREBONDS:
        add_wirebonds(design)
    return design


def info():
    Lf, Lr = lengths()
    return dict(ch=CH, L_filter=Lf, L_res=Lr, cp_plate_um=CP_PLATE, cj_gap_um=CJ_GAP_UM, D_cj=d5.D_CJ,
                L_rr=d5.VERT + d5.OPEN_DY - d5.D_CJ, style=C6.STYLE, pads_x_mm=1.2, chip=CHIP, hfss_wirebonds=WIREBONDS)


if __name__ == '__main__':
    import geopandas as gpd
    if not hasattr(gpd.GeoDataFrame, 'append'):
        gpd.GeoDataFrame.append = gpd.GeoDataFrame._append
    which = sys.argv[1] if len(sys.argv) > 1 else 'both'
    tag = '' if (CP_PLATE == 165.0 and CJ_GAP_UM == 19.0) else f'_cp{CP_PLATE:g}_cj{CJ_GAP_UM:g}'
    for v in (['filter', 'full'] if which == 'both' else [which]):
        d = build(v)
        gds = os.path.join(ROOT, 'gds', f'verify_{v}{tag}.gds'); os.makedirs(os.path.dirname(gds), exist_ok=True)
        pl.export_gds(d, gds)
        png = os.path.join(ROOT, 'fig', f'verify_{v}{tag}.png'); os.makedirs(os.path.dirname(png), exist_ok=True)
        pl.plot_gds(gds, png, window=(-2.0, 2.0, -3.6, 1.0))
        pl.plot_gds(gds, png.replace('.png', '_zoom.png'), window=(-0.9, 1.1, -3.1, 0.3))
        print(v, info(), 'components', len(d.components), '->', gds)
