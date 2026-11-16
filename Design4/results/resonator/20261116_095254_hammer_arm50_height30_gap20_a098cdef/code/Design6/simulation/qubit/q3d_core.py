"""Q3D capacitance matrix of the Design6 core: Q1, TC, Q2 with their readout claws, XY lines and the TC charge stub.

Rendered: Q1, Q2 (Xmon + readout claw), TC (claw coupler, flux line drawn as ground), XY routes + their open ends.
Open terminations: the readout claw pins (resonator side), the XY route starts (launch pad side), TC charge pin.
Solver (skill.md section 12): 8 passes, Medium order.

From the Maxwell matrix (all nets other than the three islands grounded):
    Ec_i  = e^2 / 2 * Cinv[i,i] / h                      g_ij(n) = 4 e^2 Cinv[i,j] / h   (H_int = g n_i n_j)
    readout g (lambda/4, f_r of the channel)            cap_matrix.readout_g (lam=4, transmon convention)
    XY / TC charge-line T1 limit and drive power        cap_matrix.charge_line

Usage:  python q3d_core.py <tag> [json overrides for build_design, e.g. '{"xmon": {"cross_length": "170um"}}']
Output: simulation/results/<date>/q3d_<tag>/{cmat_fF.csv, summary.json}
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np

SIM = _S('Design6/simulation')
DESIGN = _S('Design6'); REPO = _S('.')
sys.path[:0] = [_I('Design6/layout'), _I('note'),
                _I('Design5/simulation')]
import design6_layout as L   # noqa: E402
import cap_matrix as cm      # noqa: E402

e, h = 1.602176634e-19, 6.62607015e-34
F_RO = {k: v['f'] * 1e9 for k, v in L.READOUT.items()}
F_Q = dict(Q1=4.80e9, Q2=4.40e9, TC=6.3e9)


def nodes_of(cmat):
    f = lambda *k: cm.find_node(cmat, *k)
    return dict(Q1=f('cross_Q1'), Q2=f('cross_Q2'), TC=f('tmon_TC'),
                R1=f('readout_connector_arm_Q1'), R2=f('readout_connector_arm_Q2'),
                XY1=f('xy_Q1', 'xy_end_Q1'), XY2=f('xy_Q2', 'xy_end_Q2'))


def analyse(cmat, verbose=True):
    n = nodes_of(cmat)
    names = [k for k in n]
    C = cm.reduced_matrix(cmat, [n[k] for k in names])                  # F (cap_matrix converts), others grounded
    Ci = np.linalg.inv(C)
    idx = {k: i for i, k in enumerate(names)}
    # qubit subsystem: Kron-free (claws / lines are AC ground for the qubit Hamiltonian)
    q = ['Q1', 'TC', 'Q2']
    Cq = C[np.ix_([idx[k] for k in q], [idx[k] for k in q])]
    Cqi = np.linalg.inv(Cq)
    out = dict(nodes=n, C_maxwell_fF={a: {b: float(C[idx[a], idx[b]] * 1e15) for b in names} for a in names})
    out['Ec_GHz'] = {k: float(e ** 2 / 2 * Cqi[i, i] / h / 1e9) for i, k in enumerate(q)}
    out['C_sigma_fF'] = {k: float(1e15 / Cqi[i, i]) for i, k in enumerate(q)}
    out['g_n_MHz'] = {f'{q[i]}-{q[j]}': float(4 * e ** 2 * Cqi[i, j] / h / 1e6)
                      for i in range(3) for j in range(i + 1, 3)}
    out['C_mutual_fF'] = {'Q1-TC': float(-C[idx['Q1'], idx['TC']] * 1e15), 'Q2-TC': float(-C[idx['Q2'], idx['TC']] * 1e15),
                          'Q1-Q2': float(-C[idx['Q1'], idx['Q2']] * 1e15)}
    ro = {}
    for qn, rn in (('Q1', 'R1'), ('Q2', 'R2')):
        ec = out['Ec_GHz'][qn]; ej = (F_Q[qn] / 1e9 + ec) ** 2 / (8 * ec)
        ro[qn] = cm.readout_g(cmat, (n[qn],), n[rn], Z0=50.0, f_r=F_RO[qn], lam=4,
                              transmon_EJ_over_Ec=ej / ec, verbose=verbose)
    out['readout'] = ro
    ch = {}
    for qn, ln in (('Q1', 'XY1'), ('Q2', 'XY2')):
        ec = out['Ec_GHz'][qn]; ej = (F_Q[qn] / 1e9 + ec) ** 2 / (8 * ec)
        n01 = (ej / (8 * ec)) ** 0.25 / np.sqrt(2)
        ch[qn] = cm.charge_line(cmat, (n[qn],), n[ln], F_Q[qn], n01)
    out['charge_lines'] = ch
    full = [x for x in cmat.index if 'ground' not in x.lower()]
    Cf = cm.reduced_matrix(cmat, full); Cfi = np.linalg.inv(Cf); ix = {k: i for i, k in enumerate(full)}
    xt = {}
    for ln, own in (('XY1', 'Q1'), ('XY2', 'Q2')):
        b = {qq: float(Cfi[ix[n[qq]], ix[n[ln]]] / Cfi[ix[n[ln]], ix[n[ln]]]) for qq in q}
        cd = {qq: float(-cmat.loc[n[ln], n[qq]]) for qq in q}
        xt[ln] = dict(beta=b, crosstalk_pct={qq: 100 * abs(b[qq] / b[own]) for qq in q if qq != own},
                      C_line_fF=cd, direct_C_ratio_pct={qq: 100 * cd[qq] / cd[own] for qq in q if qq != own})
    out['xy_crosstalk'] = xt
    if verbose:
        print('Ec (GHz):', {k: round(v, 4) for k, v in out['Ec_GHz'].items()})
        print('C_sigma (fF):', {k: round(v, 2) for k, v in out['C_sigma_fF'].items()})
        print('g_n (MHz, H = g n_i n_j):', {k: round(v, 2) for k, v in out['g_n_MHz'].items()})
        print('C mutual (fF):', {k: round(v, 3) for k, v in out['C_mutual_fF'].items()})
        for k, v in ch.items():
            print(f'charge line {k}:', {a: round(b, 4) for a, b in v.items()})
        for k, v in out['xy_crosstalk'].items():
            print(f'{k} crosstalk: drive ratio % of own', {a: round(b, 2) for a, b in v['crosstalk_pct'].items()},
                  '| direct C ratio %', {a: round(b, 2) for a, b in v['direct_C_ratio_pct'].items()})
    return out


def run(tag, overrides=None):
    from qiskit_metal.analyses.quantization import LOManalysis
    import save_aedt
    overrides = overrides or {}
    # simulation geometry: shorted flux line drawn as ground (skill 8.2); SQUID electrodes end flush with the ground
    # (the 3 um mask overlap makes Q3D fail: net squid_ground overlaps ground_main_plane)
    overrides.setdefault('coupler', {}).update(flux_as_ground='True', squid_contact_overlap='0um', squid_shift_toward_stem='0um')
    d = L.build_design(channels=False, **overrides)
    save_aedt.use_project(_S('Design6/simulation/aedt/Design6.aedt'))
    lom = LOManalysis(d, 'q3d')
    for k, v in dict(max_passes=8, min_passes=2, percent_error=0.5, freq_ghz=5.0,
                     solution_order='Medium', auto_increase_solution_order=False).items():
        lom.sim.setup[k] = v
    comps = ['Q1', 'Q2', 'TC', 'xy_Q1', 'xy_Q2', 'xy_end_Q1', 'xy_end_Q2']
    opens = [('Q1', 'readout'), ('Q2', 'readout'), ('xy_Q1', 'start'), ('xy_Q2', 'start')]
    t0 = time.time()
    lom.sim.run(name=save_aedt.design_name(f'core_{tag}'), components=comps, open_terminations=opens)
    cmat = lom.sim.capacitance_matrix
    out_dir = _relocate(os.path.join(SIM, 'results', time.strftime('%Y-%m-%d'), f'q3d_{tag}'))
    os.makedirs(out_dir, exist_ok=True)
    cmat.to_csv(_relocate(os.path.join(out_dir, 'cmat_fF.csv')))
    print(cmat.round(3).to_string())
    res = analyse(cmat)
    res.update(tag=tag, overrides=overrides, minutes=round((time.time() - t0) / 60, 2))
    try:
        conv = lom.sim.renderer.pinfo.setup.get_convergence()
        (conv[0] if isinstance(conv, tuple) else conv).to_csv(_relocate(os.path.join(out_dir, 'convergence.csv')))
    except Exception as ex:
        print('convergence not saved:', repr(ex)[:100])
    json.dump(res, open(_relocate(os.path.join(out_dir, 'summary.json')), 'w'), indent=1, default=float)
    print('saved', out_dir, f'({res["minutes"]} min)', flush=True)
    return res


if __name__ == '__main__':
    if sys.argv[1:2] in (['--help'], ['-h']):
        print('Usage: q3d_core.py [tag] [overrides_json]\n'
              '       q3d_core.py --analyze <run folder>\n'
              'A normal run connects to AEDT; --analyze reads a saved matrix.')
        sys.exit(0)
    if sys.argv[1:2] == ['--analyze']:          # re-analyse a saved matrix: --analyze <run folder>
        import pandas as pd
        print(json.dumps(analyse(pd.read_csv(_relocate(os.path.join(sys.argv[2], 'cmat_fF.csv')), index_col=0), verbose=False)['Ec_GHz']))
        sys.exit()
    tag = sys.argv[1] if len(sys.argv) > 1 else 'base'
    ov = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
    run(tag, ov)
