"""Design6_TXT PurcellReadout candidate: simulation GDS, junction-free mask preview, layout figures and checks.

Builds design6_layout.build_design(readout='purcell_resonator_pair') (PurcellReadout channels + C_in), exports
  gds/Design6_<rev>_sim.gds        simulation geometry (junction markers on layer 50, as export_simulation_gds)
  gds/Design6_<rev>_mask.gds       junction-free, uncheesed mask preview (not the wafer fab mask)
  fig/Design6_<rev>_{layout,channel_q1,channel_q2,cin}.png
and records pins, CPW spacing (>= 400 um), 4 um metal width and drawn lengths in parameters/<rev>_build_checks.json.
The wafer fab mask (cheese, die ID, rounded corners) is produced by code/general/export_wafer_test_junctions.py
after the candidate is confirmed.

    python -I -B Design6_TXT/layout/purcell_resonator_pair_candidate.py [--loads <q3d_readout run folder>] [--rev 2026-11-19_purcell_resonator_pair_r1]
"""
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'Design6_TXT' / 'layout'), str(ROOT / 'Design8' / 'layout')]
import numpy as np  # noqa: E402

import design6_layout as L  # noqa: E402


def loads_from(run):
    if run is None:
        return None
    d = json.loads((Path(run) / 'summary.json').read_text(encoding='utf-8'))['derived']
    return {k: d[k] for k in ('cj', 'head_excess', 'cp_mutual', 'cp_plate_ground', 'claw_end') if k in d}


def pin_pairs(design):
    out = []
    for name in ('Q1', 'Q2'):
        s = '_' + name
        for a, b in ((('cp' + s, 'b'), ('purcell_resonator_pair' + s, 'cp')),):
            pa, pb = design.components[a[0]].pins[a[1]], design.components[b[0]].pins[b[1]]
            out.append(dict(pair=f'{a[0]}.{a[1]} - {b[0]}.{b[1]}',
                            offset_um=float(np.linalg.norm(np.asarray(pa.middle) - np.asarray(pb.middle)) * 1e3),
                            normal_dot=float(np.dot(pa.normal, pb.normal))))
    return out


def main():
    import gds_tools
    import check_metal_width as cmw
    from build import spacing_check
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--loads', type=Path)
    ap.add_argument('--rev', default='2026-11-19_purcell_resonator_pair_r1')
    a = ap.parse_args()
    loads = loads_from(a.loads)
    home = ROOT / 'Design6_TXT'
    sim = home / 'gds' / f'Design6_{a.rev}_sim.gds'
    mask = home / 'gds' / f'Design6_{a.rev}_mask.gds'
    d = L.build_design(purcell_resonator_pair_loads_fF=loads)
    L.export_simulation_gds(d, sim)
    m = L.build_design(jj_metal=False, purcell_resonator_pair_loads_fF=loads)
    exempt = cmw.design_exemptions(m)
    spacing = spacing_check(m)
    pins = pin_pairs(m)
    channels = m.metadata['channels']
    lengths = m.metadata['lengths']
    gds_tools.export_gds(m, str(mask))
    found = cmw.check(mask, exempt=exempt)
    fig = home / 'fig' / f'Design6_{a.rev}'
    gds_tools.plot_gds(str(mask), f'{fig}_layout.png', window=(-5, 5, -5, 5))
    for name in ('Q1', 'Q2'):
        x = channels[name]['jr_x']
        gds_tools.plot_gds(str(mask), f'{fig}_channel_{name.lower()}.png',
                           window=(x - 1.6, x + 1.6, -0.6, 2.6), figsize=(9, 7), dpi=140)
    cin = m.components['cin'].pins['b'].middle
    gds_tools.plot_gds(str(mask), f'{fig}_cin.png', window=(cin[0] - .5, cin[0] + .5, cin[1] - .3, cin[1] + .9),
                       figsize=(6, 7), dpi=140)
    out = dict(revision=a.rev, readout='purcell_resonator_pair', loads=L.purcell_resonator_pair_loads(None, loads), loads_source=str(a.loads) if a.loads else
               'parameters/readout_purcell_resonator_pair.json estimates', lengths_mm=lengths,
               channels={k: dict({kk: vv for kk, vv in v.items() if kk != 'length_report'}, length_report=v['length_report'])
                         for k, v in channels.items()},
               feed_y_mm=m.metadata['feed_y_mm'], pins=pins, spacing_violations=spacing,
               metal_width_regions_below_4um=[dict(width_um=w, bounds=g.bounds) for g, w in found],
               gds=[str(p.relative_to(ROOT)) for p in (sim, mask)],
               figures=[str(Path(f'{fig}_{v}.png').relative_to(ROOT)) for v in ('layout', 'channel_q1', 'channel_q2', 'cin')])
    (home / 'parameters' / f'{a.rev}_build_checks.json').write_text(json.dumps(out, indent=1, default=float), encoding='utf-8')
    print(json.dumps(dict(lengths=lengths, pins=pins, spacing=spacing, below_4um=len(found)), indent=1, default=float))


if __name__ == '__main__':
    main()
