"""Design5_6_purcell-res PurcellReadout chip: Q3D loads, then HFSS eigenmode and driven S21 per channel pair.

    python -I -B Design5_6_purcell-res/simulation/resonator/purcell_resonator_pair_validate.py q3d
    python -I -B Design5_6_purcell-res/simulation/resonator/purcell_resonator_pair_validate.py eigen  --q3d <q3d run> --pair 2,4 [--passes 6]
    python -I -B Design5_6_purcell-res/simulation/resonator/purcell_resonator_pair_validate.py driven --q3d <q3d run> --pair 2,4 [--passes N]

q3d (8 passes): cj (hammer coupon), ref (two straight 16/9 CPWs), cp (C_p / C_in T-type cap), module (channel 4 alone,
shorted ends rendered open). The open resonator end uses the CPW open-end estimate (OpenToGround, as drawn).
eigen / driven fixture: chip6_purcell_resonator_pair.build_chip(channels=pair) with the Q3D-load lengths: the whole feed line (C_in, all
six n0 tees; tees of channels not built end in an open stub) and the channels in --pair. Since 2026-11-19 the check
is one run with all six channels (--pair 1,2,3,4,5,6), 1 pass for eigen and driven (user quick check; 10 um seed,
out-of-core on 32 GB). Driven uses one Fast wide sweep at 0.2 MHz, no fine sweep (user choice). Ports: 50 Ohm at feed_in start and at the last feed segment end.
"""
from pathlib import Path
import argparse
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator'), str(ROOT / 'Design5_6_purcell-res' / 'layout')]
import numpy as np  # noqa: E402
import chip6_purcell_resonator_pair as C  # noqa: E402  (purcell_6500: SI_EPS 11.5)
import purcell_readout_model as PRM  # noqa: E402
import purcell_resonator_pair_sim as S  # noqa: E402
from purcell_readout import PurcellReadout  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

DESIGN = 'Design5_6_purcell-res'


def small_design(sx, sy):
    d = C.d5.designs.DesignPlanar()
    d.overwrite_enabled = True
    d.chips.main.material = C.d5.SI_MATERIAL
    d.chips.main.size.update(size_x=f'{sx}mm', size_y=f'{sy}mm', size_z='-500um', center_x='0mm', center_y='0mm')
    return d


def q3d_cases(cfg):
    lens, _ = C.lengths(cfg)
    opts = C.module_options(cfg, lens[4]['L_filter'], lens[4]['L_res'])
    cases = {}
    d = small_design(2.0, 2.0)
    PurcellReadout(d, 'cj', options=dict(opts, draw_mode='coupon', coupon_lead='400um'))
    cases['cj'] = (d, ['cj'], [], ('filter', 'resonator'))
    d = small_design(2.0, 2.0)
    S.TwoLines(d, 'ref', options=dict(spacing=opts['line_spacing'], half_length='400um', **C.P.CPW))
    cases['ref'] = (d, ['ref'], [], ('line_f', 'line_r'))
    d = small_design(1.5, 1.0)
    o = {k: f'{v}um' for k, v in dict(C.d5.ENC, **C.P.ENC_CAI).items()}
    C.d5.EnclosingCap(d, 'cp', options=dict(o, **C.P.CPW))
    cases['cp'] = (d, ['cp'], [('cp', 'a'), ('cp', 'b')], ('frame|lead_a', 'plate|lead_b'))
    d = small_design(4.0, 4.0)
    PurcellReadout(d, 'purcell_resonator_pair', options=dict(opts, short_ends='open'))
    cases['module'] = (d, ['purcell_resonator_pair'], [('purcell_resonator_pair', 'cp'), ('purcell_resonator_pair', 'qubit')], ('filter', 'resonator'))
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
    if 'module' in r and 'cj' in out:
        out.update(module_mutual=r['module']['mutual_fF'], parasitic_C=r['module']['mutual_fF'] - out['cj'],
                   parasitic_fraction=(r['module']['mutual_fF'] - out['cj']) / out['cj'])
    return out


def q3d_loads(run):
    d = json.loads((Path(run) / 'summary.json').read_text(encoding='utf-8'))['derived']
    return {k: float(d[k]) for k in ('cj', 'head_excess', 'cp_mutual', 'cp_plate_ground')}


def fixture(loads, pair):
    lens, _ = C.lengths(C.config(), loads)
    d = C.build_chip(lens, channels=pair)
    comps = [n for n, c in d.components.items() if type(c).__name__ not in ('BondPad', 'LaunchpadWirebond')]
    last = f'fl_{len(C.order(C.config()))}'
    ports = [('feed_in', 'start', 50), (last, 'end', 50)]
    return d, comps, ports, [], lens


def predictions(lens, pair):
    return {f'ch{c}': dict(f_p_GHz=lens[c]['f_target_GHz'], f_r_GHz=lens[c]['f_target_GHz'],
                           J_hammer_MHz=lens[c]['J_model_MHz'], kappa_p_MHz=lens[c]['kappa_p_model_MHz'],
                           L_filter_mm=lens[c]['L_filter'], L_res_mm=lens[c]['L_res']) for c in pair}


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
    names = sorted(pred, key=lambda k: pred[k]['f_r_GHz'])
    windows = S.channel_windows(f, s, {q: pred[q]['f_r_GHz'] for q in names})
    out = {}
    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.plot(f, 20 * np.log10(abs(s)), 'k-', lw=2, label='HFSS Driven S21')
    for q, w in windows.items():
        r = S.s21_channel(f, s, w)
        ax.plot(r['f'], 20 * np.log10(abs(r['fit'])), '--', lw=1, label=f'{q} two-pole fit (zeros fixed)')
        out[q] = dict(r['poles'], zeros=r['zeros'])
    ax.set_xlabel('Frequency (GHz)'); ax.set_ylabel('|S21| (dB)'); ax.grid(alpha=.3); ax.legend(fontsize=8)
    ax.set_title(f"Design5_6_purcell-res PurcellReadout {'+'.join(names)} (HFSS {raw['passes']} passes, dS {raw['last_delta_s']:.3g})",
                 fontsize=9)
    fig.tight_layout(); fig.savefig(Path(folder) / 's21_channels.png', dpi=150); plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['q3d', 'eigen', 'driven'])
    ap.add_argument('--q3d', type=Path)
    ap.add_argument('--pair', default='2,4')
    ap.add_argument('--passes', type=int)
    ap.add_argument('--cases', default='cj,ref,cp,module')
    a = ap.parse_args()
    cfg = C.config()
    inputs = [ROOT / 'Design5_6_purcell-res/parameters/readout_purcell_resonator_pair.json']
    pair = [int(x) for x in a.pair.split(',')]
    if a.mode == 'q3d':
        params = dict(mode='q3d', cases=a.cases.split(','), max_passes=8, config=cfg)
        upstream, tag = [], 'purcell_resonator_pair_q3d'
    else:
        if a.q3d is None:
            ap.error('--q3d is required for eigen / driven')
        loads = q3d_loads(a.q3d.resolve())
        passes = a.passes or (6 if a.mode == 'eigen' else 3)
        lens, _ = C.lengths(cfg, loads)
        pred = predictions(lens, pair)
        params = dict(mode=a.mode, pair=pair, q3d_loads=loads, prediction=pred, passes=passes, mesh=S.MESH,
                      wire_bonds=S.WB, straps='PurcellReadout.strap_points, PEC 10 um', config=cfg,
                      pass_budget_note='skill.md section 12 default' if a.passes is None else
                      'explicit pass budget (deliberate deviation from skill.md section 12; 1 pass for eigenmode and S21 as a quick check, approved by the user on 2026-11-19 (full-feed pair fixture too large for more passes; 6-pass eigen ran out of core))')
        inputs.append(a.q3d.resolve() / 'summary.json')
        upstream, tag = [str(a.q3d.resolve())], f"purcell_resonator_pair_{a.mode}_ch{'_'.join(map(str, pair))}"
    folder = create_run(DESIGN, 'resonator', tag, entrypoint=__file__, inputs=inputs, upstream=upstream, parameters=params)
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
            d, comps, ports, opens, _ = fixture(loads, pair)
            S.write(folder / 'raw' / 'geometry.json', S.geometry_dump(d))
            fs = [pred[k]['f_r_GHz'] for k in pred]
            if a.mode == 'eigen':
                raw = S.run_eigen(folder, d, comps, ports, opens, minf=min(fs) - 0.15, n_modes=2 * len(pair), passes=passes,
                                  tag='D5E')
                got = S.eigen_channels(raw, {k: pred[k]['f_r_GHz'] for k in pred})
            else:
                raw = S.run_driven(folder, d, comps, ports, opens, f0=float(np.mean(fs)),
                                   band=(min(fs) - 0.25, max(fs) + 0.25), passes=passes, tag='D5S',
                                   fine=False, wide_step_ghz=0.0002)
                got = driven_channels(folder, raw, pred)
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
