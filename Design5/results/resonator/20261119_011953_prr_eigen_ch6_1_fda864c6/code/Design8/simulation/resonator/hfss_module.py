"""Design8 HFSS: eigenmode and driven S21 of the PurcellReadout channel (open fixture, no qubit).

Geometry: Design8/layout/build.py build('open') with the filter / resonator lengths that the frequency model
(purcell_readout_model.length_for_f) gives for f_target when it is fed the Q3D loads of an upstream q3d_module run
(C_j, head excess, C_p mutual + plate ground). So the HFSS result tests "frequency from length" and
"J from the Q3D hammer" without any fitting.

    python -I -B Design8/simulation/resonator/hfss_module.py eigen  --q3d <q3d_module run folder>
    python -I -B Design8/simulation/resonator/hfss_module.py driven --q3d <q3d_module run folder>

Solver rules (skill.md sections 12-13): eigenmode 6 passes, driven 3 passes, first order, 10 um length mesh on every
metal object; qiskit auto wire bonds every 350 um on every CPW path plus PEC straps at PurcellReadout.strap_points()
(both hammer stubs, both main lines on each side of each T). Two 50 Ohm lumped feed ports (eigenmode: as port sheets,
so the filter mode is lossy). Bare f_p, f_r, J, kappa_p from the two-mode inversion
[[f_p + i kappa_p/2, J], [J, f_r]] of the two complex eigenvalues / S21 poles (one lossy mode, real J, lossless
resonator; note/purcell_filter_method.md section 4, skill/purcell_resonator_design.md step 5).
"""
from pathlib import Path
import argparse
import json
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator'), str(ROOT / 'Design8' / 'layout')]
import project_paths  # noqa: E402
project_paths.configure(__file__)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from qiskit_metal import Dict  # noqa: E402

import build as B  # noqa: E402
import purcell_readout_model as M  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

WB = dict(wb_threshold='350um', wb_offset='0um', wb_size=2)
MESH = '10um'
STRAP_DIAMETER_MM = 0.010


def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, default=float) + '\n', encoding='utf-8')


def q3d_loads(q3d_run):
    s = json.loads((Path(q3d_run) / 'summary.json').read_text(encoding='utf-8'))['derived']
    return {k: float(s[k]) for k in ('cj', 'head_excess', 'cp_mutual', 'cp_plate_ground')}


def prediction(lengths):
    g = B.GEOMETRY
    cpw = lengths['cpw']
    m = g['module']
    loads = lengths['loads']
    fp = M.tapped_lambda4_f(lengths['filter_mm'], m['filter_tap_mm'], lengths['c_end_filter_fF'],
                            lengths['stub_mm'], lengths['c_stub_fF'], cpw)
    fr = M.tapped_lambda4_f(lengths['resonator_mm'], m['resonator_tap_mm'], lengths['c_end_resonator_fF'],
                            lengths['stub_mm'], lengths['c_stub_fF'], cpw)
    j = M.j_hammer((fp + fr) / 2, loads['cj'], lengths['filter_mm'], m['filter_tap_mm'], lengths['c_end_filter_fF'],
                   lengths['resonator_mm'], m['resonator_tap_mm'], lengths['c_end_resonator_fF'], cpw)
    kp = M.kappa_hanger(fp, loads['cp_mutual'], cpw)
    return dict(f_p_GHz=fp, f_r_GHz=fr, J_hammer_MHz=j, kappa_p_MHz=kp,
                method='purcell_readout_model: tapped lambda/4 TL root, closed-form J (hammer only), hanger kappa_p')


def invert(z1, z2):
    """Two complex modes (GHz, Im = kappa/2 > 0) -> bare f_p, f_r, J, kappa_p (one lossy mode, real J)."""
    s, p = z1 + z2, z1 * z2
    fr = p.imag / s.imag
    fp = s.real - fr
    j2 = fp * fr - p.real
    return dict(f_p_GHz=fp, f_r_GHz=fr, J_MHz=float(np.sign(j2) * np.sqrt(abs(j2)) * 1e3),
                kappa_p_MHz=2 * s.imag * 1e3, J2_sign=int(np.sign(j2)))


def draw_straps(renderer, design):
    from pyEPR.ansys import parse_units
    names = []
    for st in design.components['prr'].strap_points():
        a = np.radians(st['angle_deg'])
        names.append(renderer.modeler.draw_wirebond(
            pos=parse_units(np.array([st['x'], st['y']])), ori=np.array([np.cos(a), np.sin(a)]),
            width=parse_units(st['span']), height=parse_units(st['span']), z=0,
            wire_diameter=f'{STRAP_DIAMETER_MM}mm', NumSides=6, name='g_wb', material='pec', solve_inside=False))
    return names


def mesh_all(renderer, design, comps):
    id2n = {c.id: n for n, c in design.components.items()}
    objs = []
    for kind in ('path', 'poly'):
        for r in design.qgeometry.tables[kind].itertuples():
            comp = id2n[int(r.component)]
            if not bool(r.subtract) and comp in comps:
                objs.append(f'{r.name}_{comp}')
    renderer.modeler.mesh_length('metal_mesh', objs, MaxLength=MESH)
    return objs


def fit_two_poles(f, s, guesses, window, degree=2):
    """S21 = polynomial(degree) + sum r_k / (f - p_k), complex poles p_k (Im > 0), linear least squares in the
    residues, nonlinear in the poles (as Design6 readout_poles.fit)."""
    from scipy.optimize import least_squares
    m = (f >= window[0]) & (f <= window[1])
    f, s = f[m], s[m]
    c = np.mean(window)
    x = (f - c) * 1e3

    def basis(q):
        poles = q[:2] + 1j * np.exp(q[2:])
        return np.column_stack([*((x / 100) ** k for k in range(degree + 1)), *(1 / (x - p) for p in poles)]), poles

    def res(q):
        b, _ = basis(q)
        r = b @ np.linalg.lstsq(b, s, rcond=None)[0] - s
        return np.r_[r.real, r.imag]
    p0 = (np.array([g[0] for g in guesses]) - c) * 1e3
    w0 = np.log([g[1] for g in guesses])
    q = least_squares(res, np.r_[p0, w0], max_nfev=5000).x
    b, poles = basis(q)
    rms = float(np.sqrt(np.mean(np.abs(b @ np.linalg.lstsq(b, s, rcond=None)[0] - s) ** 2)))
    return c + poles / 1e3, rms


def build_case(lengths):
    d, info = B.build('open', lengths=lengths)
    comps = [n for n, c in d.components.items() if type(c).__name__ != 'ResEnd']
    return d, comps, info


def run_eigen(folder, lengths, args):
    from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim
    d, comps, info = build_case(lengths)

    class Eig(EigenmodeSim):
        def _analyze(self):
            info['straps'] = draw_straps(self.renderer, d)
            info['mesh_objects'] = len(mesh_all(self.renderer, d, comps))
            super()._analyze()
    eig = Eig(d, 'hfss')
    eig.renderer._options.update(WB)
    for k, v in dict(min_freq_ghz=args.minf, n_modes=2, max_delta_f=0.1, max_passes=args.passes or 6, min_passes=2,
                     min_converged=1, pct_refinement=20, basis_order=1, vars=Dict()).items():
        eig.setup[k] = v
    t0 = time.time()
    eig.run(name=f'D8E_{folder.name[-8:]}', components=comps, open_terminations=[], port_list=B.PORTS)
    setup = eig.renderer.pinfo.setup
    conv = setup.get_convergence()[0]
    conv.to_csv(folder / 'raw' / 'eigen_convergence.csv')
    f, im = setup.get_solutions().eigenmodes()
    z = [complex(a, b) for a, b in zip(f, im)]
    eig.renderer.pinfo.project._project.Save()
    raw = dict(frequencies_GHz=list(map(float, f)), imag_GHz=list(map(float, im)), passes=int(len(conv)),
               last_delta_f_pct=float(conv.iloc[-1, -1]), minutes=(time.time() - t0) / 60,
               aedt_design=eig.renderer.pinfo.design.name, straps=info['straps'], mesh_objects=info['mesh_objects'])
    write(folder / 'raw' / 'eigen_modes.json', raw)
    return dict(raw=raw, inversion=invert(z[0], z[1]) if len(z) == 2 else None)


def run_driven(folder, lengths, args):
    from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
    from scipy.signal import find_peaks
    d, comps, info = build_case(lengths)

    class Drv(ScatteringImpedanceSim):
        def _analyze(self):
            info['straps'] = draw_straps(self.renderer, d)
            info['mesh_objects'] = len(mesh_all(self.renderer, d, comps))
            super()._analyze()
    sim = Drv(d, 'hfss')
    sim.renderer._options.update(WB)
    f0 = B.GEOMETRY['f_target_GHz']
    lo, hi = f0 - 0.5, f0 + 0.5
    sim.setup.update(freq_ghz=f0, max_delta_s=0.02, max_passes=args.passes or 3, min_passes=1, min_converged=1,
                     pct_refinement=20, basis_order=1, vars=Dict())
    sim.setup.sweep_setup.update(name='Wide', start_ghz=lo, stop_ghz=hi, count=1001, type='Fast', save_fields=False)
    t0 = time.time()
    sim.run(name=f'D8S_{folder.name[-8:]}', components=comps, open_terminations=[], port_list=B.PORTS)
    _, _, S = sim.renderer.get_params(['S11', 'S21'])
    S.to_csv(folder / 'raw' / 'S_wide.csv')
    conv = sim.renderer.pinfo.setup.get_convergence()[0]
    conv.to_csv(folder / 'raw' / 'hfss_convergence.csv')
    fw = S.index.to_numpy(float)
    sw = S['S21'].to_numpy(complex)
    idx, _ = find_peaks(-20 * np.log10(np.maximum(abs(sw), 1e-12)), prominence=1)
    if len(idx):
        flo, fhi = max(lo, fw[idx].min() - 0.06), min(hi, fw[idx].max() + 0.06)
    else:
        flo, fhi = f0 - 0.1, f0 + 0.1
    count = int(np.ceil((fhi - flo) / 0.00005)) + 1
    sim.renderer.add_sweep(setup_name='Setup', name='Fine', start_ghz=flo, stop_ghz=fhi, count=min(count, 10001),
                           type='Fast', save_fields=False)
    sim.renderer.analyze_sweep('Fine', 'Setup')
    _, _, S = sim.renderer.get_params(['S11', 'S21'])
    S.to_csv(folder / 'raw' / 'S_fine.csv')
    sim.renderer.pinfo.project._project.Save()
    raw = dict(passes=int(len(conv)), last_delta_s=float(conv.iloc[-1, -1]), minutes=(time.time() - t0) / 60,
               wide_GHz=[lo, hi], fine_GHz=[flo, fhi], wide_dips_GHz=fw[idx].tolist(),
               aedt_design=sim.renderer.pinfo.design.name, straps=info['straps'], mesh_objects=info['mesh_objects'])
    write(folder / 'raw' / 'driven.json', raw)
    import prr_sim
    f, s = prr_sim.load_s(folder / 'raw' / 'S_fine.csv')
    r = prr_sim.s21_channel(f, s, (f.min(), f.max()))
    return dict(raw=raw, fit=dict(inversion=r['poles'], zeros=r['zeros']))


def analyse_s21(path, folder, dips):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    S = pd.read_csv(path, index_col=0)
    f = S.index.to_numpy(float)
    s = S['S21'].map(lambda v: complex(str(v).replace('i', 'j'))).to_numpy()
    s = np.conj(s)          # HFSS e^{+jwt} -> poles with Im > 0 in the fitted form (method note section 3)
    if len(dips) >= 2:
        guesses = [(dips[0], 5.0), (dips[-1], 5.0)]
    else:
        f0 = f[np.argmin(abs(s))]
        guesses = [(f0 - 0.01, 5.0), (f0 + 0.01, 5.0)]
    window = (f.min(), f.max())
    poles, rms = fit_two_poles(f, s, guesses, window)
    inv = invert(poles[0], poles[1])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(f, 20 * np.log10(abs(s)), 'k-', label='HFSS Driven S21')
    for p in poles:
        ax.axvline(p.real, color='tab:red', ls=':', lw=1)
    ax.set_xlabel('Frequency (GHz)'); ax.set_ylabel('|S21| (dB)'); ax.grid(alpha=.3)
    ax.set_title(f"fit: f_p {inv['f_p_GHz']:.4f}, f_r {inv['f_r_GHz']:.4f} GHz, J {inv['J_MHz']:.2f} MHz, "
                 f"kappa_p {inv['kappa_p_MHz']:.1f} MHz", fontsize=9)
    ax.legend(); fig.tight_layout(); fig.savefig(folder / 's21_fit.png', dpi=150); plt.close(fig)
    return dict(poles=[[p.real, p.imag] for p in poles], fit_rms=rms, inversion=inv,
                assumptions='two poles + quadratic background over the fine band, conj(S21); one lossy mode, real J')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['eigen', 'driven'])
    ap.add_argument('--q3d', type=Path, required=True)
    ap.add_argument('--minf', type=float, default=6.0)
    ap.add_argument('--passes', type=int)
    args = ap.parse_args()
    q3d = args.q3d.resolve()
    loads = q3d_loads(q3d)
    lengths = B.model_lengths('open', caps=loads)
    pred = prediction(lengths)
    params = dict(mode=args.mode, lengths=lengths, prediction=pred, q3d_loads=loads, mesh=MESH, wire_bonds=WB,
                  straps='PurcellReadout.strap_points, PEC, 10 um wire', ports=B.PORTS,
                  passes=args.passes or (6 if args.mode == 'eigen' else 3), geometry='Design8 build open (no qubit)',
                  geometry_revision=B.GEOMETRY['revision'])
    folder = create_run('Design8', 'resonator', f'hfss_{args.mode}', entrypoint=__file__,
                        inputs=[ROOT / 'Design8/parameters/geometry.json', q3d / 'summary.json'], upstream=[str(q3d)],
                        parameters=params)
    print('RUN_FOLDER=' + str(folder), flush=True)
    print('prediction', json.dumps(pred), flush=True)
    try:
        import save_aedt
        import hfss_dirs
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(ROOT / 'Design8' / 'simulation'))
        rec = json.loads((folder / 'run.json').read_text(encoding='utf-8'))
        rec['environment']['solver_version'] = save_aedt._desktop().get_version()
        rec['geometry_revision'] = B.GEOMETRY['revision']
        write(folder / 'run.json', rec)
        result = run_eigen(folder, lengths, args) if args.mode == 'eigen' else run_driven(folder, lengths, args)
        inv = result['inversion'] if args.mode == 'eigen' else result['fit']['inversion']
        compare = dict(df_p_MHz=(inv['f_p_GHz'] - pred['f_p_GHz']) * 1e3, df_r_MHz=(inv['f_r_GHz'] - pred['f_r_GHz']) * 1e3,
                       J_ratio_hfss_over_hammer=inv['J_MHz'] / pred['J_hammer_MHz'],
                       kappa_p_ratio=inv['kappa_p_MHz'] / pred['kappa_p_MHz'])
        summary = dict(prediction=pred, hfss=result, comparison=compare,
                       validation='fixed pass budget; convergence values recorded in raw/, not claimed converged')
        print(json.dumps(dict(hfss=inv, comparison=compare), indent=1), flush=True)
        finish_run(folder, status='completed', summary=summary)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
