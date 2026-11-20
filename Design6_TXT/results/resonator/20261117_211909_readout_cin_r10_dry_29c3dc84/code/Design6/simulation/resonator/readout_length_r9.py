"""Propose Design6 r9 tail length corrections with the saved variant G circuit."""
from pathlib import Path
import json
import sys
import traceback

import numpy as np
from scipy.optimize import root
from shapely import wkt

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code/general'))
import readout_model_q3d as G
import readout_s21 as R
from run_record import create_run, finish_run

DESIGN = ROOT / 'Design6'
RESULTS = DESIGN / 'results/resonator'
PARENT = RESULTS / '20261117_133256_readout_model_q3d_07c856b2'
MODEL = RESULTS / '20261117_073940_readout_model_a50c32aa'
GEOMETRY = RESULTS / '20261117_120307_readout_q3d_full_40c53be2/raw/geometry.json'
NAMES = ('Q1', 'Q2')
TARGET = {q: R.L.READOUT[q]['f'] for q in NAMES}
FREQ_KEYS = ('bare_fp_GHz', 'bare_fr_GHz')
METRICS = FREQ_KEYS + ('J_MHz', 'kappa_p_MHz')


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def coupled_sections(base, geom, q, delta_filter):
    ids = {name: i for i, name in enumerate(geom['options'], 1)}
    lines = {r['component']: wkt.loads(r['wkt']) for r in geom['path'] if not r['subtract']}
    res_open = lines[ids['res_open_' + q]]
    filter_tail = lines[ids['filter_tail_' + q]]
    below = min(res_open.coords[0][1] - res_open.coords[1][1],
                filter_tail.coords[0][1] - filter_tail.coords[1][1] + delta_filter / 2)
    return dict(base, below_mm=below, res_open_free_mm=G.P.RES_TAP - G.P.TEE_HALF - below,
                overlap_mm=base['above_mm'] + below + 2 * G.P.TEE_HALF)


def complex_pole(seed, model_args):
    def residual(offset):
        freq = seed + complex(*offset) / 1000
        det = np.linalg.det(np.eye(2) + 50 * G.port_Y(freq, *model_args[0], **model_args[1]))
        return (det.real, det.imag)

    solution = root(residual, [0, 0], tol=1e-10)
    if np.linalg.norm(residual(solution.x)) > 1e-6:
        raise RuntimeError(f'Could not recover complex pole near {seed}: {solution.message}')
    return seed + complex(*solution.x) / 1000


def main():
    ref = load(PARENT / 'compare.json')
    pm = load(MODEL / 'run.json')['parameters']
    geom = load(GEOMETRY)
    wafer = load(DESIGN / 'design.json')['wafer_candidate']
    base = {q: np.asarray(pm['lengths_mm'][q]) for q in NAMES}
    cw = R.L.CLAW_DL
    for q in NAMES:
        lengths = wafer['readout_lengths_mm'][q]
        np.testing.assert_allclose(base[q], [lengths['L_filter'], lengths['L_res'] + cw], atol=1e-6, rtol=0)
    c = ref['corrections']
    geo = dict(pm['geometry'])
    for k in ('z0_ohm', 'velocity_m_s'):
        geo[k] *= c['line_scale_F']
    _, cpa, cpb = pm['cp_fF']
    cp = [c['cp_mutual_full_q3d_fF'], cpa, cpb - c['delta_used_fF']]
    kw = dict(cm=G.CM_PAIR, feed_res=c['feed_res_fF'], extra_cj=c['extra_lumped_cj_fF'])
    baseline = {q: {k: ref['variants']['G'][q][k]['model'] for k in METRICS} for q in NAMES}
    folder = create_run('Design6', 'resonator', 'readout_length_r9', entrypoint=__file__,
                        inputs=[PARENT / 'compare.json', MODEL / 'run.json', GEOMETRY, DESIGN / 'design.json'],
                        upstream=[str(PARENT), str(MODEL), str(GEOMETRY.parents[1])],
                        parameters=dict(parent_geometry=wafer['revision'], f_target_GHz=TARGET, model='variant G',
                                        claw_dl_mm=cw, EM_rerun=False, hfss_fit=False,
                                        sweep_step_um=20, tuning='both shorted tail lengths; fixed taps and Cp',
                                        candidate_only=True))
    print('RUN_FOLDER=' + str(folder), flush=True)
    trials = []
    try:
        def evaluate(deltas_um):
            offsets = {q: np.asarray(deltas_um[2 * i:2 * i + 2]) / 1000 for i, q in enumerate(NAMES)}
            lens = {q: base[q] + offsets[q] for q in NAMES}
            sec = {q: coupled_sections(ref['sections_mm'][q], geom, q, offsets[q][0]) for q in NAMES}
            args = (lens, geo, cp, pm['hammer_deembedded_fF'], sec)
            out = {}
            for q in NAMES:
                b = baseline[q]
                fp, fr = np.asarray([b[k] for k in FREQ_KEYS]) * base[q] / lens[q]
                J, kappa = b['J_MHz'] / 1000, b['kappa_p_MHz'] / 1000
                seeds = np.linalg.eigvals([[fp + .5j * kappa, J], [J, fr]])
                poles = [complex_pole(p, (args, kw)) for p in seeds]
                if abs(poles[0] - poles[1]) < .003 or min(p.imag for p in poles) <= 0:
                    raise RuntimeError(f'Invalid {q} pole pair: {poles}')
                out[q] = G.invert(*poles)
                out[q]['poles_GHz'] = [[p.real, p.imag] for p in sorted(poles, key=lambda p: p.real)]
            trials.append(dict(deltas_um=np.asarray(deltas_um).tolist(), channels=out))
            return out, args

        start, _ = evaluate(np.zeros(4))
        for q in NAMES:
            for k in METRICS:
                factor = 1000 if k.endswith('GHz') else 1
                if abs(start[q][k] - baseline[q][k]) * factor > 1:
                    raise RuntimeError(f'{q} {k}: complex-pole baseline disagrees with the saved variant G fit')
        for direction in np.r_[20 * np.eye(4), -20 * np.eye(4)]:
            evaluate(direction)

        def targets(delta):
            channels, _ = evaluate(delta)
            return [(channels[q][k] - TARGET[q]) * 1000 for q in NAMES for k in FREQ_KEYS]

        delta0 = [(start[q][k] / TARGET[q] - 1) * base[q][i] * 1000
                  for q in NAMES for i, k in enumerate(FREQ_KEYS)]
        opt = root(targets, delta0)
        end, args = evaluate(opt.x)
        if max(abs(r) for r in targets(opt.x)) > .01:
            raise RuntimeError(f'Length solve failed: {opt.message}')
        lengths = {q: dict(L_filter=float(args[0][q][0]), L_res=float(args[0][q][1] - cw)) for q in NAMES}
        fitted, spectrum = G.bare_of(*args, **kw)
        for q in NAMES:
            assert all(abs(fitted[q][k] - end[q][k]) * (1000 if k.endswith('GHz') else 1) < 1
                       for k in METRICS), (q, fitted[q], end[q])
        R.L.build_design(wafer_ports=True, readout_hammer=wafer['readout_hammer_um'], lengths=lengths)
        np.savez(folder / 'raw/candidate_s21.npz', f_GHz=G.M.GRID, S21=spectrum)
        summary = dict(baseline=start, tuned=end, tuned_s21_fit=fitted, lengths_fab_mm=lengths,
                       lengths_no_claw_mm={q: args[0][q].tolist() for q in NAMES},
                       delta_um={q: dict(filter=float(opt.x[2*i]), res=float(opt.x[2*i+1]))
                                 for i, q in enumerate(NAMES)},
                       sections_mm=args[-1], load='No claw or qubit, 0.100 mm restoration carried over from r9',
                       status='circuit-model candidate; loaded HFSS verification pending')
        (folder / 'raw/length_trials.json').write_text(json.dumps(trials, indent=2), encoding='utf-8')
        finish_run(folder, status='completed', summary=summary)
        print(json.dumps(summary, indent=2))
    except Exception:
        (folder / 'raw/length_trials.json').write_text(json.dumps(trials, indent=2), encoding='utf-8')
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
