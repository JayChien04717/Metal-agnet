# -*- coding: utf-8 -*-
r"""HFSS DrivenModal flux-line mutual inductance + coupler-to-coupler crosstalk for Design13 (rail-SQUID coupler).

Method (skill/mutual_inductance_analysis.md):
- Fixture: qubit Q0 at (0, 0) and its two nearest couplers TC_E (orientation 0, at (+1, 0) mm) and TC_N
  (orientation 90, at (0, +1) mm); the two SQUIDs are 1.41 mm apart (nearest foreign SQUID in the lattice).
- Ports: both flux-bias-line pins (50 Ohm lumped) + jj1 of each coupler as PORT ONLY (open-loop EMF);
  jj2 keeps a lumped Lj (TC_E 6.5 nH, TC_N 7.5 nH) so the two JJ ports are identified by Im(Z_ii)/w.
- M_ij = Im(Z_ji) / (2 pi f) at f_min = 10 MHz (discrete 10-100 MHz); reciprocity, Re Z, plateau checks.
- Crosstalk = M(loop_i <- line_j) / M(loop_j <- line_j).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parents[3]
DESIGN_DIR = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'aedt'), str(DESIGN_DIR / 'layout')]
import project_paths
project_paths.configure(__file__)

import numpy as np
from qiskit_metal import Dict, designs

import build as D13
from run_record import create_run, finish_run

PHI0 = 2.067833848e-15
DEVS = [('TC_E', 1.0, 0.0, '0', '6.5nH'), ('TC_N', 0.0, 1.0, '90', '7.5nH')]


def build_fixture(q_over, c_over, fbl_over):
    d = designs.DesignPlanar()
    d.overwrite_enabled = True
    d.chips.main.material = 'silicon_eps11p45'
    d.chips.main.size.update(dict(size_x='5mm', size_y='5mm', size_z='-500um', center_x='0mm', center_y='0mm'))
    q_opt = dict(D13.QUBIT_OPT, jj_metal='False', make_readout_pad='False', **q_over)
    D13.ScalablePadQubit(d, 'Q0', options=dict(**D13._pos(0.0, 0.0), ro_corner='NE', jj_corner='SW', **q_opt))
    for nm, x, y, ori, lj in DEVS:
        c_opt = dict(D13.COUPLER_OPT, **c_over)
        c_opt['Lj'] = lj
        c_opt['jj_metal'] = 'True'
        c_opt['flux_as_ground'] = 'False'
        from floating_tunable_coupler import FloatingTunableCoupler
        base_fbl = dict(FloatingTunableCoupler.default_options.flux_line_options)
        base_fbl.update(fbl_over)
        c_opt['flux_line_options'] = Dict(**base_fbl)
        D13.FloatingTunableCoupler(d, nm, options=dict(**D13._pos(x, y), orientation=ori, **c_opt))
    return d


def analyse(Z, freqs):
    """Port order: [fl_E, fl_N, jj1_E, jj1_N] as rendered; identify by Im Zii (jj ports have nH)."""
    w = 2 * np.pi * freqs
    n = Z.shape[1]
    L0 = np.array([Z[0, i, i].imag / w[0] for i in range(n)])
    jj = sorted(np.argsort(L0)[-2:], key=lambda i: L0[i])        # smaller Lj -> TC_E (6.5 nH)
    jj_map = {'TC_E': int(jj[0]), 'TC_N': int(jj[1])}
    fl_pool = [i for i in range(n) if i not in jj_map.values()]
    fl_map = {}
    for dev in ('TC_E', 'TC_N'):
        cand = [i for i in fl_pool if i not in fl_map.values()]
        fl_map[dev] = int(max(cand, key=lambda i: abs(Z[0, jj_map[dev], i].imag)))
    M, M_curve = {}, {}
    recip, reim = {}, {}
    for li in ('TC_E', 'TC_N'):
        for lj in ('TC_E', 'TC_N'):
            a, b = jj_map[li], fl_map[lj]
            zm = 0.5 * (Z[:, a, b] + Z[:, b, a])
            M_curve[f'{li}<-{lj}'] = (zm.imag / w * 1e12).tolist()
            M[f'{li}<-{lj}'] = float(zm[0].imag / w[0] * 1e12)
            recip[f'{li}<-{lj}'] = float(np.max(np.abs(Z[:, a, b] - Z[:, b, a]) / np.maximum(np.abs(zm), 1e-18)))
            reim[f'{li}<-{lj}'] = float(abs(zm[0].real) / max(abs(zm[0].imag), 1e-18))
    own = {'TC_E': M['TC_E<-TC_E'], 'TC_N': M['TC_N<-TC_N']}
    xt = {'TC_N<-fl_E': 100 * M['TC_N<-TC_E'] / own['TC_E'], 'TC_E<-fl_N': 100 * M['TC_E<-TC_N'] / own['TC_N']}
    spread = {d: float(np.ptp(M_curve[f'{d}<-{d}']) / abs(own[d]) * 100) for d in own}
    return dict(port_map=dict(jj=jj_map, flux=fl_map), L0_by_port_nH=(L0 * 1e9).tolist(),
                M_pH=M, M_self_pH=own, crosstalk_pct=xt, M_spread_pct=spread,
                I_phi0_mA={d: float(PHI0 / abs(own[d] * 1e-12) * 1e3) for d in own},
                L_loop_nH={d: float(Z[0, jj_map[d], jj_map[d]].imag / w[0] * 1e9) for d in own},
                reciprocity_rel_err=recip, re_over_im=reim, f_MHz=(freqs / 1e6).tolist(), M_curves_pH=M_curve)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tag', required=True)
    ap.add_argument('--qubit-overrides', default='{}')
    ap.add_argument('--coupler-overrides', default='{}')
    ap.add_argument('--fbl-overrides', default='{}', help='e.g. {"fbl_type": "T", "fbl_length": "100um", "fbl_x": "30um", "fbl_height": "80um"}')
    ap.add_argument('--min-passes', type=int, default=12)
    ap.add_argument('--max-passes', type=int, default=12)
    ap.add_argument('--max-delta-s', type=float, default=0.0002)
    ap.add_argument('--npts', type=int, default=10)
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()
    q_over, c_over, f_over = (json.loads(x) for x in (a.qubit_overrides, a.coupler_overrides, a.fbl_overrides))
    d = build_fixture(q_over, c_over, f_over)
    comps = ['Q0'] + [n for n, *_ in DEVS]
    ports = [(n, 'flux_bias_line', 50) for n, *_ in DEVS]
    jj_ports = [(n, 'jj1', 50, False) for n, *_ in DEVS]
    setup = dict(freq_ghz=0.1, max_delta_s=a.max_delta_s, max_passes=a.max_passes, min_passes=a.min_passes,
                 min_converged=3, pct_refinement=30, basis_order=1,
                 sweep=dict(start_ghz=0.01, stop_ghz=0.10, count=a.npts, type='Discrete'))
    run_dir = create_run('Design13_ScalableFloatingCoupler2Q', 'flux', f'hfss_mutual_{a.tag}', entrypoint=__file__,
                         inputs=[DESIGN_DIR / 'design.json', DESIGN_DIR / 'layout' / 'build.py'],
                         parameters=dict(revision=D13.REVISION, components=comps, ports=ports, jj_to_port=jj_ports,
                                         qubit_overrides=q_over, coupler_overrides=c_over, fbl_overrides=f_over,
                                         setup=setup, note='Stage 2 check requested by user (flux line M / crosstalk)'))
    print('RUN_FOLDER=' + str(run_dir), flush=True)
    if a.dry_run:
        finish_run(run_dir, status='completed', summary=dict(dry_run=True))
        return

    from qiskit_metal.analyses.simulation.scattering_impedance import ScatteringImpedanceSim
    import hfss_dirs
    import save_aedt
    hfss_dirs.use_d_drive()
    project = DESIGN_DIR / 'simulation' / 'aedt' / 'D13F.aedt'
    project.parent.mkdir(parents=True, exist_ok=True)
    save_aedt.use_project(project)

    sim = ScatteringImpedanceSim(d, 'hfss')
    for k, v in dict(freq_ghz=0.1, max_delta_s=a.max_delta_s, max_passes=a.max_passes, min_passes=a.min_passes,
                     min_converged=3, pct_refinement=30, basis_order=1,
                     sweep_setup=Dict(name='Sweep', start_ghz=0.01, stop_ghz=0.10, count=a.npts, step_ghz=None,
                                      type='Discrete', save_fields=False)).items():
        sim.setup[k] = v
    t0 = time.time()
    try:
        sim.run_sim(name=f'F_{run_dir.name[-4:]}', components=comps, open_terminations=[], port_list=ports,
                    jj_to_port=jj_ports, box_plus_buffer=True)
    except Exception:
        try:
            for msg in sim.renderer._pinfo.desktop._desktop.GetMessages('', '', 0):
                print('  AEDT:', msg, flush=True)
        except Exception as exc:
            print('Could not read AEDT messages:', exc, flush=True)
        finish_run(run_dir, status='failed', summary=dict(error='HFSS Analyze failed'))
        raise
    n = 4
    pnames = [f'Z{i}{j}' for i in range(1, n + 1) for j in range(1, n + 1)]
    freqs, curves, _ = sim.renderer.get_params(pnames)
    freqs = np.asarray(freqs, float)
    Z = np.zeros((len(freqs), n, n), complex)
    for nm, cv in zip(pnames, curves):
        Z[:, int(nm[1]) - 1, int(nm[2]) - 1] = np.asarray(cv, complex)
    (run_dir / 'raw').mkdir(parents=True, exist_ok=True)
    np.save(run_dir / 'raw' / 'Z.npy', Z)
    np.save(run_dir / 'raw' / 'f.npy', freqs)
    conv_info = None
    try:
        conv_df, _ = sim.renderer.pinfo.setup.get_convergence()
        conv_df.to_csv(run_dir / 'raw' / 'convergence.csv')
        conv_info = dict(passes=int(len(conv_df)), last=str(conv_df.iloc[-1].to_dict()))
    except Exception as exc:
        print('convergence not saved:', repr(exc), flush=True)
    res = analyse(Z, freqs)
    res.update(elapsed_s=round(time.time() - t0, 1), setup=setup, convergence=conv_info, tag=a.tag)
    (run_dir / 'summaries').mkdir(parents=True, exist_ok=True)
    (run_dir / 'summaries' / 'summary.json').write_text(json.dumps(res, indent=2), encoding='utf-8')
    finish_run(run_dir, status='completed', validation='single pass count; 12-vs-15 pass check not run', summary=res)
    print(json.dumps({k: res[k] for k in ('tag', 'M_self_pH', 'crosstalk_pct', 'M_spread_pct', 'I_phi0_mA',
                                          'L_loop_nH', 'reciprocity_rel_err', 're_over_im', 'port_map',
                                          'elapsed_s', 'convergence')}, indent=2), flush=True)


if __name__ == '__main__':
    main()
