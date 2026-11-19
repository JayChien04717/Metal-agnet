"""Design5 six-channel chip with the shared PurcellReadout module (2026-11-19; replaces chip6_d4.py channels).

Same chip as chip6_d4.py: 10 x 10 mm, 16/9 um CPW, Si eps_r 11.5, launches at x = +-4.5 mm on one straight feed
line (y = 0), C_in (Cai T-type enclosing cap) at x = -3.9 mm, channels 1-6 at 6.0-7.0 GHz, three above and three
mirrored below the feed, staggered (positions as chip6_d4.ORDER). Each channel: n0 tee - C_p (same T-type cap) -
PurcellReadout (400 um line spacing, taps 1.2 / 1.3 mm, hammer from parameters/readout_prr.json) - open resonator
end (OpenToGround; no qubit). Lengths: purcell_readout_model.length_for_f with the loads of readout_prr.json or of a
Design5 q3d_readout run (--loads).

    python -I -B Design5/layout/chip6_prr.py [--loads <q3d run>] [--rev 2026-11-19_prr_r1]
        -> gds/Design5_<rev>_sim.gds, fig/Design5_<rev>_{layout,channel}.png, parameters/<rev>_build_checks.json
"""
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'Design5' / 'layout'), str(ROOT / 'Design4' / 'layout'),
                str(ROOT / 'Design8' / 'layout')]
import project_paths  # noqa: E402
project_paths.configure(__file__)
import numpy as np  # noqa: E402
import purcell_6500 as P  # noqa: E402  (sets SI_EPS 11.5 before design5_layout)
import design5_layout as d5  # noqa: E402
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround  # noqa: E402
from purcell_readout import PurcellReadout  # noqa: E402
import purcell_readout_model as PRM  # noqa: E402

HOME = ROOT / 'Design5'
CFG_FILE = HOME / 'parameters' / 'readout_prr.json'
ORDER = [(2, -3.20, False), (4, -1.90, True), (6, -0.60, False), (1, 0.70, True), (3, 2.00, False), (5, 3.30, True)]
PORT_X, X_CIN = 4.5, -3.9


def config():
    return json.loads(CFG_FILE.read_text(encoding='utf-8'))


def loads(cfg, q3d=None):
    h = cfg['hammer']
    out = dict(cfg['loads'])
    out.update(q3d or {})
    if out.get('cj') is None:
        out['cj'] = PRM.hammer_cj(h['head_gap_um'], h['head_height_um'], h['head_width_um'])
    if out.get('head_excess') is None:
        out['head_excess'] = PRM.hammer_excess_ground(h['head_gap_um'], h['head_height_um'])
    if out.get('open_end') is None:
        out['open_end'] = PRM.open_end_c(cfg['cpw']['model_preset'])
    return out


def lengths(cfg=None, q3d=None):
    cfg = cfg or config()
    L = loads(cfg, q3d)
    m, cpw = cfg['module'], cfg['cpw']['model_preset']
    stub, c_stub, c_end_f = m['line_spacing_um'] / 2000, L['cj'] + L['head_excess'], L['cp_mutual'] + L['cp_plate_ground']
    out = {}
    for ch, f in cfg['channels'].items():
        lf = PRM.length_for_f(f, m['filter_tap_mm'], c_end_f, stub, c_stub, cpw)
        lr = PRM.length_for_f(f, m['resonator_tap_mm'], L['open_end'], stub, c_stub, cpw)
        out[int(ch)] = dict(f_target_GHz=f, L_filter=lf, L_res=lr,
                            J_model_MHz=PRM.j_hammer(f, L['cj'], lf, m['filter_tap_mm'], c_end_f, lr,
                                                     m['resonator_tap_mm'], L['open_end'], cpw),
                            kappa_p_model_MHz=PRM.kappa_hanger(f, L['cp_mutual'], cpw), c_end_filter_fF=c_end_f,
                            c_end_resonator_fF=L['open_end'], c_stub_fF=c_stub)
    return out, L


def module_options(cfg, Lf, Lr, **extra):
    m, h = cfg['module'], cfg['hammer']

    def um(v):
        return f'{v:.6f}um'
    return dict(filter_length=f'{Lf:.9f}mm', resonator_length=f'{Lr:.9f}mm', filter_tap=f"{m['filter_tap_mm']}mm",
                resonator_tap=f"{m['resonator_tap_mm']}mm", line_spacing=um(m['line_spacing_um']),
                tee_half=um(m['tee_half_um']), lead=um(m['lead_um']), jog=um(m['jog_um']),
                meander_pitch=um(m['meander_pitch_um']), fillet=um(m['fillet_um']),
                filter_meander_height=f"{m['filter_meander_height_mm']}mm",
                resonator_meander_height=f"{m['resonator_meander_height_mm']}mm",
                filter_first_run=m['filter_first_run'], resonator_first_run=m['resonator_first_run'],
                head_width=um(h['head_width_um']), head_height=um(h['head_height_um']), head_gap=um(h['head_gap_um']),
                hammer_ground_gap=um(h['ground_gap_um']), corner_radius=um(h['corner_radius_um']), **P.CPW, **extra)


def add_channel(design, cfg, ch, x0, up, Lf, Lr):
    """n0 tee at (x0, 0); C_p; PurcellReadout with its cp pin on C_p pin b (orientation 180 above the feed)."""
    s, sg = f'_ch{ch}', 1 if up else -1
    d5.TapJunction(design, 'n0' + s, options=dict(**P.CPW, **P.position(x0, 0), orientation='0', stub_length='100um',
                                                  stub_up=str(up)))
    enc_half = (d5.ENC['plate_length'] / 2 + d5.ENC['leg_ext'] + d5.ENC['lead_length']) / 1000
    opts = {k: f'{v}um' for k, v in dict(d5.ENC, **P.ENC_CAI).items()}
    d5.EnclosingCap(design, 'cp' + s, options=dict(opts, **P.CPW, **P.position(x0, sg * (0.100 + enc_half)),
                                                   orientation='90' if up else '-90'))
    P.straight(design, 'stub_p' + s, ('n0' + s, 'second_end'), ('cp' + s, 'a'))
    pb = np.asarray(design.components['cp' + s].pins['b'].middle, float)
    m = cfg['module']
    # local cp pin at (-s/2, filter_tap); orientation 180 maps local (x, y) -> (-x, -y)
    lx, ly = -m['line_spacing_um'] / 2000, m['filter_tap_mm']
    x, y = (pb[0] + lx, pb[1] + ly) if up else (pb[0] - lx, pb[1] - ly)
    prr = PurcellReadout(design, 'prr' + s, options=module_options(cfg, Lf, Lr, **P.position(x, y),
                                                                   orientation='180' if up else '0'))
    q = prr.pins['qubit']
    OpenToGround(design, 'readout_open' + s, options=dict(**P.position(*q.middle), orientation='90' if up else '270',
                                                          width=P.CPW['trace_width'], gap=P.CPW['trace_gap']))
    for key, tgt in (('filter_mm', Lf), ('resonator_mm', Lr)):
        if abs(prr.length_report[key] - tgt) > 1e-6:
            raise ValueError(f'ch{ch}: drawn {key} {prr.length_report[key]:.6f} vs {tgt:.6f}')
    return prr


def build_chip(lens=None, cfg=None, q3d=None, channels=None, wire_bonds=True):
    """channels: subset of channel numbers to draw (simulation fixtures); the feed line, C_in and all n0 tees stay."""
    cfg = cfg or config()
    if lens is None:
        lens, _ = lengths(cfg, q3d)
    design = d5.designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.material = d5.SI_MATERIAL
    design.chips.main.size.update(dict(size_x='10mm', size_y='10mm', size_z='-500um', center_x='0mm', center_y='0mm'))
    design.variables.update(cpw_width=P.CPW['trace_width'], cpw_gap=P.CPW['trace_gap'])
    for name, x, o in [('LPPin', -PORT_X, '0'), ('LPPout', PORT_X, '180')]:
        d5.add_bond_pad(design, name, dict(P.position(x, 0), orientation=o))
    P.enclosing(design, 'cin', X_CIN, 0, 0)
    for ch, x0, up in ORDER:
        if channels is None or ch in channels:
            add_channel(design, cfg, ch, x0, up, lens[ch]['L_filter'], lens[ch]['L_res'])
        else:      # tee only, its stub ends open (OpenToGround)
            d5.TapJunction(design, f'n0_ch{ch}', options=dict(**P.CPW, **P.position(x0, 0), orientation='0',
                                                              stub_length='100um', stub_up=str(up)))
            e = design.components[f'n0_ch{ch}'].pins['second_end']
            OpenToGround(design, f'n0_open_ch{ch}', options=dict(**P.position(*e.middle), orientation='90' if up else '270',
                                                                width=P.CPW['trace_width'], gap=P.CPW['trace_gap']))
    n = lambda k: f'n0_ch{ORDER[k][0]}'
    P.straight(design, 'feed_in', ('LPPin', 'tie'), ('cin', 'a'))
    P.straight(design, 'fl_0', ('cin', 'b'), (n(0), 'prime_start'))
    for k in range(1, len(ORDER)):
        P.straight(design, f'fl_{k}', (n(k - 1), 'prime_end'), (n(k), 'prime_start'))
    P.straight(design, f'fl_{len(ORDER)}', (n(len(ORDER) - 1), 'prime_end'), ('LPPout', 'tie'))
    if wire_bonds:
        t = design.qgeometry.tables['path']
        t['hfss_wire_bonds'] = t['subtract'].astype(bool)
    design.metadata['lengths'] = lens
    return design


def channel_distance(design):
    from shapely.ops import unary_union
    id2n = {c.id: nm for nm, c in design.components.items()}
    geo = {}
    for ch, _, _ in ORDER:
        if f'prr_ch{ch}' not in design.components:
            continue
        geo[ch] = unary_union([r.geometry.buffer(float(r.width) / 2) if k == 'path' else r.geometry
                               for k in ('path', 'poly') for r in design.qgeometry.tables[k].itertuples()
                               if id2n[int(r.component)] == f'prr_ch{ch}' and bool(r.subtract)])
    keys = list(geo)
    d = {f'{a}-{b}': float(geo[a].distance(geo[b])) for i, a in enumerate(keys) for b in keys[i + 1:]}
    return min(d.values()) if d else None, d


def main():
    import gds_tools
    import check_metal_width as cmw
    from build import spacing_check
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--loads', type=Path)
    ap.add_argument('--rev', default='2026-11-19_prr_r1')
    a = ap.parse_args()
    q3d = None
    if a.loads:
        dd = json.loads((a.loads / 'summary.json').read_text(encoding='utf-8'))['derived']
        q3d = {k: dd[k] for k in ('cj', 'head_excess', 'cp_mutual', 'cp_plate_ground') if k in dd}
    cfg = config()
    lens, L = lengths(cfg, q3d)
    d = build_chip(lens, cfg)
    exempt = cmw.design_exemptions(d)
    spacing = spacing_check(d)
    dmin, dists = channel_distance(d)
    sim = HOME / 'gds' / f'Design5_{a.rev}_sim.gds'
    gds_tools.export_gds(d, str(sim))
    found = cmw.check(sim, exempt=exempt)
    fig = HOME / 'fig' / f'Design5_{a.rev}'
    gds_tools.plot_gds(str(sim), f'{fig}_layout.png', window=(-5, 5, -5, 5))
    gds_tools.plot_gds(str(sim), f'{fig}_channel.png', window=(-1.0, 2.4, -0.2, 3.3), figsize=(8, 8), dpi=140)
    out = dict(revision=a.rev, readout='prr', loads=L, loads_source=str(a.loads) if a.loads else
               'parameters/readout_prr.json estimates', channels={str(k): v for k, v in lens.items()},
               length_report={ch: d.components[f'prr_ch{ch}'].length_report for ch, _, _ in ORDER},
               spacing_violations=spacing, min_channel_etch_distance_mm=dmin,
               metal_width_regions_below_4um=[dict(width_um=w, bounds=g.bounds) for g, w in found],
               gds=str(sim.relative_to(ROOT)), figures=[str(Path(f'{fig}_{v}.png').relative_to(ROOT)) for v in ('layout', 'channel')])
    (HOME / 'parameters' / f'{a.rev}_build_checks.json').write_text(json.dumps(out, indent=1, default=float), encoding='utf-8')
    for ch, v in sorted(lens.items()):
        print(f"ch{ch} {v['f_target_GHz']} GHz: L_f {v['L_filter']:.4f} L_r {v['L_res']:.4f} mm, J {v['J_model_MHz']:.2f}, "
              f"kappa_p {v['kappa_p_model_MHz']:.1f} MHz")
    print('spacing violations', len(spacing), '| min channel etch distance', dmin, '| <4um', len(found))


if __name__ == '__main__':
    main()
