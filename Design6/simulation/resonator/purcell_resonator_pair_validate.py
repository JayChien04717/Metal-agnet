"""Design6 PurcellReadout readout: Q3D loads, then HFSS eigenmode and driven S21 of the loaded readout fixture.

    python -I -B Design6/simulation/resonator/purcell_resonator_pair_validate.py q3d
    python -I -B Design6/simulation/resonator/purcell_resonator_pair_validate.py eigen  --q3d <q3d run> [--passes 6]
    python -I -B Design6/simulation/resonator/purcell_resonator_pair_validate.py driven --q3d <q3d run> [--passes N]

q3d (8 passes; local solves): cj (module hammer coupon), ref (two straight 16/9 CPWs, same outer ends), cp (C_p / C_in
T-type cap), claw (Q1 Xmon alone, junction-free, flux line as ground: claw self-C with the cross island floating =
resonator end load), module (purcell_resonator_pair_Q1 alone, shorted ends rendered open: total filter-resonator static mutual).
eigen / driven: build_design(readout='purcell_resonator_pair', control_lines=False, coupler flux as ground) with the lengths of
purcell_resonator_pair_lengths(Q3D loads); rendered: feed (C_in, both tees, both C_p), both PurcellReadout modules and open
extensions, Q1 and Q2 (junctions as lumped Lj; TC not rendered, as in the r9 / r10 loaded runs). Ports: 50 Ohm at
feed_launch start (outside C_in) and feed_out start. Predictions: f_p = f_r = f_target from the model; J closed form.
"""
from pathlib import Path
import argparse
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator'), str(ROOT / 'Design6' / 'layout')]
import numpy as np  # noqa: E402
import design6_layout as L  # noqa: E402  (imports purcell_6500: SI_EPS 11.5)
import purcell_readout_model as PRM  # noqa: E402
import purcell_resonator_pair_sim as S  # noqa: E402
from purcell_readout import PurcellReadout  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

DESIGN = 'Design6'
CHANNELS = ('Q1', 'Q2')


def small_design(sx, sy):
    d = L.d5.designs.DesignPlanar()
    d.overwrite_enabled = True
    d.chips.main.material = L.d5.SI_MATERIAL
    d.chips.main.size.update(size_x=f'{sx}mm', size_y=f'{sy}mm', size_z='-500um', center_x='0mm', center_y='0mm')
    return d


def q3d_cases(cfg):
    lens = L.purcell_resonator_pair_lengths(cfg)
    opts = L._purcell_resonator_pair_options(cfg, lens['Q1']['L_filter'], lens['Q1']['L_res'], 0.0)
    cases = {}
    d = small_design(2.0, 2.0)
    PurcellReadout(d, 'cj', options=dict(opts, draw_mode='coupon', coupon_lead='400um'))
    cases['cj'] = (d, ['cj'], [], ('filter', 'resonator'))
    d = small_design(2.0, 2.0)
    S.TwoLines(d, 'ref', options=dict(spacing=opts['line_spacing'], half_length='400um', **L.CPW_RO))
    cases['ref'] = (d, ['ref'], [], ('line_f', 'line_r'))
    d = small_design(1.5, 1.0)
    o = {k: f'{v}um' for k, v in dict(L.d5.ENC, **L.P.ENC_CAI).items()}
    L.d5.EnclosingCap(d, 'cp', options=dict(o, **L.CPW_RO))
    cases['cp'] = (d, ['cp'], [('cp', 'a'), ('cp', 'b')], ('frame|lead_a', 'plate|lead_b'))
    d = L.build_design(jj_metal=False, landing_pads=False, channels=False, control_lines=False, marks=False,
                       coupler=dict(flux_as_ground='True'))
    cases['claw'] = (d, ['Q1'], [('Q1', 'readout')], ('readout', 'cross'))
    d = L.build_design(jj_metal=False, landing_pads=False, control_lines=False, marks=False,
                       coupler=dict(flux_as_ground='True'))
    d.components['purcell_resonator_pair_Q1'].options.short_ends = 'open'
    d.components['purcell_resonator_pair_Q1'].rebuild()
    cases['module'] = (d, ['purcell_resonator_pair_Q1'], [('purcell_resonator_pair_Q1', 'cp'), ('purcell_resonator_pair_Q1', 'qubit')], ('filter', 'resonator'))
    return cases


def derive(r, cfg):
    out = {}
    cl = PRM.line(cfg['cpw']['model_preset'])['c_line']
    if 'cj' in r and 'ref' in r:
        out.update(S.hammer_loads(r['cj'], r['ref'], cl, cfg['module']['line_spacing_um'] / 1000))
        h = cfg['hammer']
        out['cj_fit_estimate'] = PRM.hammer_cj(h['head_gap_um'], h['head_height_um'], h['head_width_um'])
    if 'cp' in r:
        out.update(cp_mutual=r['cp']['mutual_fF'], cp_plate_ground=r['cp']['ground_fF']['plate'],
                   cp_frame_ground=r['cp']['ground_fF']['frame'])
    if 'claw' in r:
        c = r['claw']
        names = c['maxwell_order']
        i_claw = names.index(c['signals']['readout'])
        floating = [names.index(n) for n in names if n != c['signals']['readout']]
        K = S.kron_floating(c['maxwell_fF'], [i_claw], floating)
        out.update(claw_end=float(K[0, 0]), claw_self_maxwell=float(np.asarray(c['maxwell_fF'])[i_claw, i_claw]),
                   claw_to_cross=float(-np.asarray(c['maxwell_fF'])[i_claw, names.index(c['signals']['cross'])]))
    if 'module' in r and 'cj' in out:
        out.update(module_mutual=r['module']['mutual_fF'], parasitic_C=r['module']['mutual_fF'] - out['cj'],
                   parasitic_fraction=(r['module']['mutual_fF'] - out['cj']) / out['cj'])
    return out


def q3d_loads(run):
    d = json.loads((Path(run) / 'summary.json').read_text(encoding='utf-8'))['derived']
    return {k: float(d[k]) for k in ('cj', 'head_excess', 'cp_mutual', 'cp_plate_ground', 'claw_end')}


def fixture(loads):
    d = L.build_design(control_lines=False, coupler=dict(flux_as_ground='True'), marks=False, purcell_resonator_pair_loads_fF=loads)
    t = d.qgeometry.tables['path']
    t['hfss_wire_bonds'] = t['subtract'].astype(bool)
    comps = ['feed_launch', 'cin', 'feed_in', 'feed_mid', 'feed_out', 'Q1', 'Q2']
    for q in CHANNELS:
        comps += [f'{p}_{q}' for p in ('n0', 'stub_p', 'cp', 'purcell_resonator_pair', 'res_open')]
    ports = [('feed_launch', 'start', 50), ('feed_out', 'start', 50)]
    opens = [(q, 'flux_bias_line') for q in CHANNELS]
    return d, comps, ports, opens


def predictions(loads):
    lens = L.purcell_resonator_pair_lengths(None, loads)
    return {q: dict(f_p_GHz=v['f_target_GHz'], f_r_GHz=v['f_target_GHz'], J_hammer_MHz=v['J_model_MHz'],
                    kappa_p_MHz=v['kappa_p_model_MHz'], L_filter_mm=v['L_filter'], L_res_mm=v['L_res'])
            for q, v in lens.items()}


def compare(pred, got):
    out = {}
    for q, g in got.items():
        p = pred[q]
        if 'f_p_GHz' not in g:
            out[q] = g
            continue
        out[q] = dict(df_p_MHz=(g['f_p_GHz'] - p['f_p_GHz']) * 1e3, df_r_MHz=(g['f_r_GHz'] - p['f_r_GHz']) * 1e3,
                      J_ratio=g['J_MHz'] / p['J_hammer_MHz'], kappa_p_ratio=g['kappa_p_MHz'] / p['kappa_p_MHz'])
    return out


def driven_channels(folder, raw, pred):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    f, s = S.load_s(Path(folder) / 'raw' / raw.get('s_file', 'S_fine.csv'))
    dips = sorted(raw['wide_dips_GHz'])
    if len(dips) != 4:
        raise RuntimeError(f'expected four dips (two channels), found {dips}')
    mid = (dips[1] + dips[2]) / 2
    windows = {'Q2': (f.min(), mid), 'Q1': (mid, f.max())} if pred['Q1']['f_r_GHz'] > pred['Q2']['f_r_GHz'] else \
        {'Q1': (f.min(), mid), 'Q2': (mid, f.max())}
    out = {}
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(f, 20 * np.log10(abs(s)), 'k-', lw=2, label='HFSS Driven S21')
    for q, w in windows.items():
        r = S.s21_channel(f, s, w)
        ax.plot(r['f'], 20 * np.log10(abs(r['fit'])), '--', lw=1, label=f'{q} two-pole fit (zeros fixed)')
        out[q] = {k: v for k, v in r.items() if k not in ('f', 's', 'fit')}
    ax.set_xlabel('Frequency (GHz)'); ax.set_ylabel('|S21| (dB)'); ax.grid(alpha=.3); ax.legend(fontsize=8)
    ax.set_title(f"Design6 PurcellReadout loaded S21 (HFSS {raw['passes']} passes, dS {raw['last_delta_s']:.3g})", fontsize=9)
    fig.tight_layout(); fig.savefig(Path(folder) / 's21_channels.png', dpi=150); plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['q3d', 'eigen', 'driven'])
    ap.add_argument('--q3d', type=Path)
    ap.add_argument('--passes', type=int)
    ap.add_argument('--cases', default='cj,ref,cp,claw,module')
    a = ap.parse_args()
    cfg = L.purcell_resonator_pair_config()
    inputs = [ROOT / 'Design6/parameters/readout_purcell_resonator_pair.json']
    if a.mode == 'q3d':
        params = dict(mode='q3d', cases=a.cases.split(','), max_passes=8, config=cfg)
        upstream = []
    else:
        if a.q3d is None:
            ap.error('--q3d is required for eigen / driven')
        loads = q3d_loads(a.q3d.resolve())
        pred = predictions(loads)
        passes = a.passes or (6 if a.mode == 'eigen' else 3)
        params = dict(mode=a.mode, q3d_loads=loads, prediction=pred, passes=passes, mesh=S.MESH, wire_bonds=S.WB,
                      straps='PurcellReadout.strap_points, PEC 10 um', config=cfg,
                      pass_budget_note='skill.md section 12 default' if a.passes is None else
                      'explicit pass budget (deliberate deviation from skill.md section 12; 2 passes for eigenmode and S21 approved by the user on 2026-11-19 because the full-feed chip fixture does not fit in 32 GB at more passes)')
        inputs.append(a.q3d.resolve() / 'summary.json')
        upstream = [str(a.q3d.resolve())]
    folder = create_run(DESIGN, 'resonator', f'purcell_resonator_pair_{a.mode}', entrypoint=__file__, inputs=inputs, upstream=upstream,
                        parameters=params)
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        import save_aedt
        import hfss_dirs
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(ROOT / DESIGN / 'simulation'))
        rec = json.loads((folder / 'run.json').read_text(encoding='utf-8'))
        rec['environment']['solver_version'] = save_aedt._desktop().get_version()
        rec['geometry_revision'] = cfg['revision']
        S.write(folder / 'run.json', rec)
        if a.mode == 'q3d':
            cases = q3d_cases(cfg)
            res = {c: S.q3d_case(c, *cases[c], folder, folder.name[-8:]) for c in a.cases.split(',')}
            summary = dict(results=res, derived=derive(res, cfg),
                           validation='Q3D 8-pass budget; solver_converged / last-pass change per case')
            print(json.dumps(summary['derived'], indent=1), flush=True)
        else:
            d, comps, ports, opens = fixture(loads)
            S.write(folder / 'raw' / 'geometry.json', S.geometry_dump(d))
            if a.mode == 'eigen':
                raw = S.run_eigen(folder, d, comps, ports, opens, minf=6.55, n_modes=6, passes=passes, tag='D6E')
                got = S.eigen_channels(raw, {q: pred[q]['f_r_GHz'] for q in CHANNELS})
            else:
                raw = S.run_driven(folder, d, comps, ports, opens, f0=6.9, band=(6.55, 7.25), passes=passes, tag='D6S',
                                   fine=False, wide_step_ghz=0.0002)
                got = driven_channels(folder, raw, pred)
                got = {q: dict(v['poles'], zeros=v['zeros']) for q, v in got.items()}
            summary = dict(prediction=pred, hfss_raw=raw, hfss=got, comparison=compare(pred, got),
                           validation='fixed pass budget; convergence recorded, not claimed converged')
            print(json.dumps(dict(hfss=got, comparison=summary['comparison']), indent=1, default=float), flush=True)
        finish_run(folder, status='completed', summary=summary)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
