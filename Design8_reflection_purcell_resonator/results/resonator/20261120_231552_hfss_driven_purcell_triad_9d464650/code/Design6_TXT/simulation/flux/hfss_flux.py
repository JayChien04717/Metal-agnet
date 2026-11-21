"""Flux lines of Design6_TXT: full 3 x 3 mutual-inductance (crosstalk) matrix and its frequency dependence, HFSS driven.

Method (Design1_FTF_v1/notes/mutual_inductance_method.md, extended to 6 ports):
  ports 1-3: the flux-line pins of Q1, TC, Q2 (local model, default) or, with --routed, the launch-pad ends of the
             routed lines fl_Q1, fl_TC, fl_Q2 (50 Ohm lumped ports);
  ports 4-6: one junction of each SQUID (rect_jj1), port only (I = 0): no loop current flows, the port voltage
             is the open-loop EMF, V_j = j w M_ij I_i  ->  M_ij = Im Z_ji / w.
  The other SQUID junction keeps its lumped Lj; it carries no current with the port open, but it makes the JJ ports
  distinguishable: Im Z_jj / w ~ Lj(other JJ) = Q1 31.8, Q2 37.6, TC 12.9 nH.
Rendered: Q1, TC, Q2 (claws, SQUIDs, local L flux lines); ground = bounding box + buffer.
Setup: adaptive at 100 MHz (max dS 0.0005, 8 passes), discrete sweep 10 MHz - 1 GHz.

Figures of merit
  * M_ii (pH) and I_Phi0 = Phi0 / M_ii of each element's own line.
  * crosstalk C_ij = M_ji / M_jj (% of the flux line j's own-loop coupling seen by loop i).
  * distortion: (M_ii(f) - M_ii(10 MHz)) / M_ii(10 MHz) and the phase of Z_ji over 10 MHz - 1 GHz (a flux pulse
    sees the on-chip transfer function M(f); flat = no on-chip distortion), plus Re Z11 / Im Z11 of each line.

Usage:  python driven/hfss_flux.py --tag v8 [--sep 10] [--overrides '{...}']
Output: results/<date>/flux_<tag>/{Z.csv, M_matrix.json, convergence.csv}
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd

SIM = _S('Design6_TXT/simulation')
DESIGN = _S('Design6_TXT'); REPO = _S('.')
sys.path[:0] = [_I('Design6_TXT/layout'), _I('Design5_6_purcell-res/simulation')]
PHI0 = 2.067833848e-15
NAMES = ['Q1', 'TC', 'Q2']


def run(a, out):
    import design6_layout as L
    from qiskit_metal import Dict
    from qiskit_metal.analyses.simulation.scattering_impedance import ScatteringImpedanceSim
    import save_aedt
    ov = json.loads(a.overrides)
    d = L.build_design(channels=False, **ov)
    # local model (Jay 2026-11-14: the long routed CPWs are not needed): Q1, TC, Q2 with the ports on their own
    # flux-line pins; --routed renders the routed lines to the launch pads as well (slow: > 1 h for the sweep)
    comps = ['Q1', 'TC', 'Q2'] + (['fl_Q1', 'fl_TC', 'fl_Q2'] if a.routed else [])
    ports = ([('fl_Q1', 'end', 50), ('fl_TC', 'end', 50), ('fl_Q2', 'end', 50)] if a.routed else
             [('Q1', 'flux_bias_line', 50), ('TC', 'flux_bias_line', 50), ('Q2', 'flux_bias_line', 50)])
    save_aedt.use_project(os.environ.get('AEDT_PROJECT', _S('Design6_TXT/simulation/aedt/Design6_TXT.aedt')))
    fine = [f'{r}_{c}' for c in ('Q1', 'Q2', 'TC') for r in ('flux_bias_line', 'squid_ground')
            if any((t['name'] == r) & (t['component'] == d.components[c].id) for _, t in d.qgeometry.tables['poly'].iterrows())]

    class Sim(ScatteringImpedanceSim):
        def _analyze(self):           # seed the mesh on the SQUIDs and the local flux lines (pH-level M)
            self.renderer.modeler.mesh_length('squid_flux_mesh', fine, MaxLength=a.mesh)
            print(f'mesh {a.mesh} on {fine}', flush=True)
            super()._analyze()

    sim = Sim(d, 'hfss')
    for k, v in dict(freq_ghz=0.1, max_delta_s=0.0002, max_passes=a.maxp, min_passes=a.maxp, min_converged=1,
                     pct_refinement=30, basis_order=1,
                     sweep_setup=Dict(name='Sweep', start_ghz=0.01, stop_ghz=1.0, count=a.npts, step_ghz=None,
                                      type='Discrete', save_fields=False)).items():
        sim.setup[k] = v
    t0 = time.time()
    sim.run_sim(name=save_aedt.design_name(f'flux_{a.tag}'), components=comps,
                open_terminations=[('Q1', 'readout'), ('Q2', 'readout')],
                port_list=ports,
                jj_to_port=[('Q1', 'rect_jj1', 50, False), ('TC', 'rect_jj1', 50, False), ('Q2', 'rect_jj1', 50, False)],
                box_plus_buffer=True)
    names = [f'Z{i}{j}' for i in range(1, 7) for j in range(1, 7)]
    freqs, curves, _ = sim.renderer.get_params(names)
    Z = np.zeros((len(freqs), 6, 6), complex)
    for n, c in zip(names, curves):
        Z[:, int(n[1]) - 1, int(n[2]) - 1] = np.asarray(c)
    np.save(_relocate(os.path.join(out, 'Z.npy')), Z); np.save(_relocate(os.path.join(out, 'f.npy')), np.asarray(freqs))
    try:
        sim.renderer.pinfo.setup.get_convergence()[0].to_csv(_relocate(os.path.join(out, 'convergence.csv')))
    except Exception as ex:
        print('convergence not saved:', repr(ex)[:100])
    json.dump(dict(overrides=ov, minutes=round((time.time() - t0) / 60, 1)), open(_relocate(os.path.join(out, 'run.json')), 'w'))
    try:
        save_aedt.snapshot(out, f'flux_{a.tag}')
    except Exception as ex:
        print('snapshot failed:', repr(ex)[:100])


def analyse(out):
    Z, f = np.load(_relocate(os.path.join(out, 'Z.npy'))), np.load(_relocate(os.path.join(out, 'f.npy')))
    w = 2 * np.pi * f
    L0 = np.array([Z[0, i, i].imag / w[0] for i in range(6)])        # H
    # JJ ports: the three largest self inductances; identify by the other junction's Lj (TC smallest, Q2 largest).
    # (Checked 2026-11-14: equals the creation order, flux lines Q1, Q2, TC = ports 1-3, junctions Q1, Q2, TC = 4-6.)
    jj = list(np.argsort(L0)[-3:]); jj.sort(key=lambda i: L0[i])        # TC, Q1, Q2
    jj = {'TC': jj[0], 'Q1': jj[1], 'Q2': jj[2]}
    fl_ports = [i for i in range(6) if i not in jj.values()]
    # flux port of element n = the one with the largest |M| to n's loop
    fl = {}
    for n in NAMES:
        cand = [i for i in fl_ports if i not in fl.values()]
        fl[n] = max(cand, key=lambda i: abs(Z[0, jj[n], i].imag))
    M = {n: {m: (Z[:, jj[n], fl[m]].imag / w) for m in NAMES} for n in NAMES}    # M[loop][line] vs f
    rows = []
    for k, fk in enumerate(f):
        row = dict(f_MHz=fk / 1e6)
        for n in NAMES:
            for m in NAMES:
                row[f'M_{n}<-{m}_pH'] = M[n][m][k] * 1e12
            row[f'phase_Z_{n}<-{n}_deg'] = float(np.degrees(np.angle(Z[k, jj[n], fl[n]])))
            row[f'L_line_{n}_pH'] = Z[k, fl[n], fl[n]].imag / w[k] * 1e12
            row[f'R_line_{n}_ohm'] = Z[k, fl[n], fl[n]].real
        rows.append(row)
    df = pd.DataFrame(rows); df.to_csv(_relocate(os.path.join(out, 'Z.csv')), index=False)
    ref = {n: M[n][n][0] for n in NAMES}
    res = dict(port_map=dict(jj={k: int(v) + 1 for k, v in jj.items()}, flux={k: int(v) + 1 for k, v in fl.items()}),
               M_self_pH={n: ref[n] * 1e12 for n in NAMES},
               I_phi0_mA={n: PHI0 / abs(ref[n]) * 1e3 for n in NAMES},
               crosstalk_pct={f'{n}<-{m}': 100 * M[n][m][0] / ref[m] for n in NAMES for m in NAMES if n != m},
               M_flatness_pct={n: float(100 * np.max(np.abs(M[n][n] / ref[n] - 1))) for n in NAMES},
               f_range_MHz=[float(f[0] / 1e6), float(f[-1] / 1e6)])
    json.dump(res, open(_relocate(os.path.join(out, 'M_matrix.json')), 'w'), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for n in NAMES:
        ax[0].semilogx(f / 1e6, M[n][n] / ref[n], label=f'{n}: M_self(f)/M(10 MHz)')
        for m in NAMES:
            if m != n:
                ax[1].semilogx(f / 1e6, 100 * M[n][m] / ref[m], label=f'loop {n} <- line {m}')
    ax[0].set_xlabel('f (MHz)'); ax[0].set_ylabel('normalised M'); ax[0].legend(); ax[0].set_title('on-chip flux transfer (distortion)')
    ax[1].set_xlabel('f (MHz)'); ax[1].set_ylabel('crosstalk (% of own M)'); ax[1].legend(fontsize=7); ax[1].set_title('flux crosstalk')
    fig.tight_layout(); fig.savefig(_relocate(os.path.join(out, 'flux_M.png')), dpi=140); plt.close(fig)
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', default='v8')
    ap.add_argument('--date', default=time.strftime('%Y-%m-%d'))
    ap.add_argument('--maxp', type=int, default=6)   # 4 um seed mesh on SQUIDs + local lines from pass 1 (skill.md 12); 12 passes reached 13 GB
    ap.add_argument('--mesh', default='4um')
    ap.add_argument('--npts', type=int, default=4, help='discrete sweep points 10 MHz - 1 GHz')
    ap.add_argument('--routed', action='store_true', help='include the routed flux lines to the launch pads')
    ap.add_argument('--overrides', default='{}')
    ap.add_argument('--analyze-only', action='store_true')
    a = ap.parse_args()
    out = _relocate(os.path.join(SIM, 'results', a.date, f'flux_{a.tag}'))
    os.makedirs(out, exist_ok=True)
    if not a.analyze_only:
        run(a, out)
    analyse(out)
