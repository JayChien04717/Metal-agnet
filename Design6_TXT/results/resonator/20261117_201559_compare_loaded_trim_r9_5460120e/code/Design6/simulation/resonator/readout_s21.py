"""Design6 readout-only driven S21: feed line + both Purcell filters + both readout resonators, no qubits.

Geometry: design6_layout.build_design(wafer_ports=True, readout_hammer=design.json wafer_candidate.readout_hammer_um,
lengths=lengths_from_params()), i.e. the same call as the r9 wafer export. Rendered: feed_in / feed_mid / feed_out and
every channel component of Q1 and Q2 (no Xmon, claw, coupler or control lines). The readout open run ends in an
explicit open termination where the claw would be, so resonator frequencies are NOT the loaded (claw) values.
Ground: Design4 'all'-style PEC straps per channel (hammer stubs, tee arms, CPW pieces <= 350 um), no auto bonds.
Ports: 50 ohm lumped at the two feed ends (launch pads not rendered).

Q3D (--q3d): local hammer coupon with the selected head and 55 um arms, and the enclosing Cp coupon, 20 passes / 0.1 %.
The straight-CPW 55 um reference and the 200 um parallel-pair C'm are geometrically identical to Design4 coupons and
are reused by explicit path after geometry equality checks.
"""
from argparse import ArgumentParser
import json
import shutil
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / p) for p in ('code/general', 'Design6/layout', 'Design6/simulation/resonator',
                                         'Design5/simulation/resonator', 'Design4/simulation/resonator')]
import design6_layout as L
from run_record import create_run, finish_run

D4 = ROOT / 'Design4/results/resonator'
D4_REF = D4 / '20261116_195508_h30_g40_w25_a55_q3d_cae2f56a/raw/ref_55um'
D4_PAIR = D4 / '20261116_220124_hammer_coupled_model_e1d28295'
MESH_LINE, MESH_CAP = '10um', '10um'


def build(res_extra_mm=0.0):
    home = ROOT / 'Design6'
    hammer = json.loads((home / 'design.json').read_text(encoding='utf-8'))['wafer_candidate']['readout_hammer_um']
    lengths = L.lengths_from_params()
    lengths = {q: dict(v, L_res=v['L_res'] + res_extra_mm) for q, v in lengths.items()}
    d = L.build_design(wafer_ports=True, readout_hammer=hammer, lengths=lengths)
    return d, hammer, lengths


def components(d):
    names = ['feed_in', 'feed_mid', 'feed_out']
    for q in ('Q1', 'Q2'):
        names += [n for n in d.components if n.endswith('_' + q) and not n.startswith(('xy', 'cl_', 'fl_'))]
    return names


def run_q3d(folder, hammer):
    import hammer_deembed as H
    P = H.P
    arm = hammer['cpw_centres_um'] / 2 - hammer['gap_um'] / 2 - hammer['width_um']
    P.PAPER_HAMMER.update(gap_um=hammer['gap_um'], height_um=hammer['height_um'], width_um=hammer['width_um'],
                          ground_um=hammer['ground_um'])
    out = {}
    local, comp = H.coupon(arm)
    out['hammer'] = H.solve_cap(local, comp, folder / f'raw/cj_{arm:g}um', 'D6CJ_' + folder.name[-8:],
                                max_passes=20, percent_error=.1)
    ref_local, _ = H.coupon(arm, reference=True)
    old = json.loads((D4_REF / 'geometry.json').read_text())
    new = H.geometry_dump(ref_local)
    assert all(old[k] == new[k] for k in ('poly', 'path')), 'Straight reference differs from the Design4 coupon'
    shutil.copytree(D4_REF, folder / f'raw/ref_{arm:g}um')
    out['reference'] = json.loads((D4_REF / 'caps.json').read_text())
    out['reference_source'] = str(D4_REF)
    cp = H.designs.DesignPlanar()
    cp.chips.main.material = P.d5.SI_MATERIAL
    cp.chips.main.size.update(size_x='5.2mm', size_y='4.3mm', size_z='-500um', center_x='0mm', center_y='0mm')
    P.enclosing(cp, 'cp', 0, 0, 0)
    out['cp'] = H.solve_cap(cp, 'cp', folder / 'raw/cp_coupon', 'D6CP_' + folder.name[-8:], max_passes=20, percent_error=.1)
    pair = json.loads((D4_PAIR / 'comparison.json').read_text())
    assert np.isclose(pair['pair_coupons'][0]['spacing_mm'] * 1000, hammer['cpw_centres_um'])
    assert L.CPW_RO == {'trace_width': '16um', 'trace_gap': '9um'}
    out['parallel_pair'] = dict(mutual_per_mm_fF=pair['mutual_per_mm_fF'], ground_per_mm_fF=pair['ground_per_mm_fF'],
                                source=str(D4_PAIR), reason='Identical 16/9 um CPW pair at the same 200 um centre spacing')
    half = out['reference']['ground_fF'][0] / 2
    out['hammer_excess_ground_fF'] = [g - half for g in out['hammer']['ground_fF']]
    out['arm_um'] = arm
    return out


def run_hfss(folder, d, comps, args):
    import purcell_6500 as P
    import hfss_readout as R
    import save_aedt, hfss_dirs
    from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
    from pyEPR.ansys import parse_units
    from scipy.signal import find_peaks
    from qiskit_metal import Dict
    hfss_dirs.use_d_drive()
    save_aedt.use_project(save_aedt.project_file(ROOT / 'Design6/simulation'))
    P.set_wire_bonds(d, False)
    straps, seen = [], set()
    for q in ('Q1', 'Q2'):
        for s in R.strap_list(d, [c for c in comps if c.endswith('_' + q) or c.startswith('feed')], q):
            key = tuple(np.round(s['pos'], 6))
            if key not in seen:
                seen.add(key); straps.append(s)
    straps = P._drop_crossing(straps)
    (folder / 'raw/ground_straps.json').write_text(json.dumps(straps, indent=1, default=float))
    id2n = {c.id: n for n, c in d.components.items()}
    lines, caps = [], []
    for kind in ('path', 'poly'):
        for r in d.qgeometry.tables[kind].itertuples():
            n = id2n[int(r.component)]
            if not bool(r.subtract) and n in comps:
                (caps if n.startswith(('cj_', 'cp_')) else lines).append(f'{r.name}_{n}')
    opens = getattr(args, 'open_terminations', [(f'res_open_{q}', 'end') for q in ('Q1', 'Q2')])

    class Driven(ScatteringImpedanceSim):
        def _analyze(self):
            for s in straps:
                self.renderer.modeler.draw_wirebond(
                    pos=parse_units(np.asarray(s['pos'])), ori=np.asarray(s['ori']), width=parse_units(s['span']),
                    height=parse_units(s['span']), z=0, wire_diameter=f"{s['diameter']}mm", NumSides=6,
                    name='g_wb', material='pec', solve_inside=False)
            self.renderer.modeler.mesh_length('cpw_line_mesh', lines, MaxLength=MESH_LINE)
            self.renderer.modeler.mesh_length('cap_mesh', caps, MaxLength=MESH_CAP)
            print(f'{len(straps)} straps; mesh {MESH_LINE} on {len(lines)} lines, {MESH_CAP} on {len(caps)} caps', flush=True)
            super()._analyze()

    sim = Driven(d, 'hfss')
    sim.setup.update(freq_ghz=args.freq, max_delta_s=.01, max_passes=args.max_passes, min_passes=3, min_converged=1,
                     pct_refinement=20, basis_order=1, vars=Dict())
    lo, hi = args.wide
    sim.setup.sweep_setup.update(name='Wide', start_ghz=lo, stop_ghz=hi, count=int(round((hi - lo) * 1000)) + 1,
                                 type='Fast', save_fields=False)
    sim.run(name='D6S21_' + folder.name[-8:], components=comps, open_terminations=opens,
            port_list=[('feed_in', 'start', 50), ('feed_out', 'start', 50)])
    _, _, S = sim.renderer.get_params(['S11', 'S21'])
    S.to_csv(folder / 'raw/S_wide.csv')
    conv = sim.renderer.pinfo.setup.get_convergence()[0]
    conv.to_csv(folder / 'raw/hfss_convergence.csv')
    f = S.index.to_numpy(float); s21 = S['S21'].to_numpy(complex)
    idx, _ = find_peaks(-20*np.log10(np.maximum(abs(s21), 1e-12)), prominence=1)
    bands = []
    if len(idx) and getattr(args, 'fine_step_ghz', .000025):
        lo2, hi2 = max(lo, float(f[idx].min()) - .05), min(hi, float(f[idx].max()) + .05)
        fine_step = getattr(args, 'fine_step_ghz', .000025)
        sim.renderer.add_sweep(setup_name='Setup', name='Fine', start_ghz=lo2, stop_ghz=hi2,
                               count=min(25000, int(np.ceil((hi2 - lo2) / fine_step)) + 1), type='Fast', save_fields=False)
        sim.renderer.analyze_sweep('Fine', 'Setup')
        _, _, S2 = sim.renderer.get_params(['S11', 'S21'])
        S2.to_csv(folder / 'raw/S_fine.csv')
        bands = [lo2, hi2]
    after = sim.renderer.pinfo.setup.get_convergence()[0]
    after.to_csv(folder / 'raw/hfss_convergence_after_fine.csv')
    sim.renderer.pinfo.project._project.Save()
    convergence = conv.astype(object).where(conv.notna(), None).to_dict()
    info = dict(design=sim.renderer.pinfo.design.name, solver_version=save_aedt._desktop().get_version(),
                wide_GHz=[lo, hi], fine_GHz=bands, wide_dips_GHz=f[idx].tolist(), straps=len(straps),
                adaptive_history_unchanged_by_fine=bool(conv.equals(after)),
                convergence=convergence, nominal_pass_count_is_not_convergence=True)
    (folder / 'raw/hfss.json').write_text(json.dumps(info, indent=1, default=float))
    return info


def plot(folder, title='Design6 r9 readout only: feed + Purcell + resonator, both channels, no qubits (open where claws would be)'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(10, 7))
    for a, name in zip(ax, ('S_wide.csv', 'S_fine.csv')):
        p = folder / 'raw' / name
        if not p.exists():
            continue
        S = pd.read_csv(p, index_col=0)
        for k, c in (('S21', 'black'), ('S11', 'tab:gray')):
            z = S[k].map(lambda x: complex(str(x).replace('i', 'j'))).to_numpy()
            a.plot(S.index, 20*np.log10(np.maximum(abs(z), 1e-12)), color=c, lw=1.4 if k == 'S21' else 1, label=k)
        a.set_ylabel('dB'); a.grid(alpha=.3); a.legend(); a.set_title(name.replace('.csv', '') + ' (HFSS Driven)')
    ax[-1].set_xlabel('Frequency (GHz)')
    fig.suptitle(title)
    fig.tight_layout(); fig.savefig(folder / 's21.png', dpi=160); plt.close(fig)


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--q3d', action='store_true', help='Run the local Q3D coupons only')
    ap.add_argument('--max-passes', type=int, default=12)
    ap.add_argument('--freq', type=float, default=7.0, help='Adaptive frequency, GHz')
    ap.add_argument('--wide', type=float, nargs=2, default=[6.3, 7.8])
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--res-extra-mm', type=float, default=0.0,
                    help='Add to every L_res (e.g. 0.1 = L.CLAW_DL, undoing the claw shortening in the no-claw model)')
    args = ap.parse_args()
    d, hammer, lengths = build(args.res_extra_mm)
    comps = components(d)
    tag = 'readout_caps' if args.q3d else 'readout_s21' + ('_clawdl' if args.res_extra_mm else '')
    folder = create_run('Design6', 'resonator', tag, entrypoint=__file__,
                        inputs=[ROOT / 'Design6/design.json', Path(L._length_file())],
                        upstream=[str(D4_REF), str(D4_PAIR)] if args.q3d else [],
                        parameters=dict(geometry='r9 wafer readout geometry, simulation only (no qubits)',
                                        readout_hammer_um=hammer, lengths_mm=lengths, components=comps,
                                        claw='omitted; open termination at res_open end', HFSS_max_passes=args.max_passes,
                                        HFSS_max_delta_s=.01, adaptive_GHz=args.freq, wide_GHz=args.wide,
                                        mesh=dict(line=MESH_LINE, cap=MESH_CAP), ports='50 ohm lumped at feed ends',
                                        Q3D=dict(max_passes=20, percent_error=.1) if args.q3d else None,
                                        res_extra_mm=args.res_extra_mm, claw_dl_mm=L.CLAW_DL,
                                        res_extra_reason='undo CLAW_DL so f_p ~ f_r without the claw load' if args.res_extra_mm else None,
                                        dry_run=args.dry_run))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        import hammer_deembed as H
        H.write(folder / 'raw/geometry.json', H.geometry_dump(d))
        if args.dry_run:
            finish_run(folder, status='completed', summary=dict(dry_run=True, components=comps)); return
        if args.q3d:
            summary = run_q3d(folder, hammer)
        else:
            summary = run_hfss(folder, d, comps, args)
            plot(folder)
        (folder / 'summary_out.json').write_text(json.dumps(summary, indent=1, default=float))
        finish_run(folder, status='completed', summary=summary)
        print(json.dumps(summary, default=float)[:3000], flush=True)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
