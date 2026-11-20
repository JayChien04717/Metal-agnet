"""Shared Q3D / HFSS plumbing for PurcellReadout channels (Design4, Design5, Design6 drivers).

Solver rules (skill.md sections 12-13): Q3D 8 passes; HFSS eigenmode 6 passes, driven 3 passes (a driver may pass a
larger budget deliberately and records it), first order, 10 um length mesh on every rendered metal object, qiskit
auto bonds every 350 um on every etched path, plus PEC straps at every PurcellReadout.strap_points() (hammer stubs
and both sides of each T).

Analysis (no fitted lengths or capacitances):
  invert(z1, z2)    two complex modes (GHz, Im = kappa/2 > 0) -> bare f_p, f_r, J, kappa_p (one lossy mode, real J)
  s21_channel(...)  one channel's window of a driven S21: zeros + S21 = 1 point, and a pole fit with fixed zeros
Q3D:
  q3d_case(...)     one local solve; Maxwell sub-matrix of the named nets, pass history, convergence flags
  hammer_loads(...) C_j = mutual(coupon) - mutual(straight reference); head excess per side
"""
from contextlib import chdir
from pathlib import Path
import json
import time

import numpy as np
import pandas as pd
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import QComponent

WB = dict(wb_threshold='350um', wb_offset='0um', wb_size=2)
MESH = '10um'
STRAP_DIAMETER_MM = 0.010


def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, default=float) + '\n', encoding='utf-8')


class TwoLines(QComponent):
    """Straight-line reference for the hammer coupon: two CPWs at x = +-spacing/2, y in [-half_length, half_length],
    etched open ends; nets 'line_f' and 'line_r'."""
    default_options = Dict(spacing='400um', half_length='400um', trace_width='10um', trace_gap='6um')

    def make(self):
        p = self.p
        for nm, x in (('line_f', -p.spacing / 2), ('line_r', p.spacing / 2)):
            ln = draw.LineString([(x, -p.half_length), (x, p.half_length)])
            cut = draw.LineString([(x, -p.half_length - p.trace_gap), (x, p.half_length + p.trace_gap)])
            self.add_qgeometry('poly', {nm: ln.buffer(p.trace_width / 2, cap_style=2)}, layer=1)
            self.add_qgeometry('poly', {nm + '_cut': cut.buffer(p.trace_width / 2 + p.trace_gap, cap_style=2)},
                               layer=1, subtract=True)


def geometry_dump(d):
    return {kind: [dict(name=r.name, subtract=bool(r.subtract), wkt=r.geometry.wkt,
                        width=float(getattr(r, 'width', 0) or 0)) for r in d.qgeometry.tables[kind].itertuples()]
            for kind in ('poly', 'path')}


def q3d_case(case, d, comps, opens, keys, folder, tag, max_passes=8, pct_error=0.3, freq_ghz=6.5):
    """keys: net selectors, each 'a|b' matches a net whose name contains a or b. Returns the result dict."""
    from qiskit_metal.analyses.quantization import LOManalysis
    out = Path(folder) / 'raw' / case
    out.mkdir(parents=True, exist_ok=False)
    write(out / 'geometry.json', geometry_dump(d))
    lom = LOManalysis(d, 'q3d')
    lom.sim.renderer.options['x_buffer_width_mm'] = 0
    lom.sim.renderer.options['y_buffer_width_mm'] = 0
    lom.sim.setup.update(max_passes=max_passes, min_passes=2, percent_error=pct_error, freq_ghz=freq_ghz,
                         solution_order='Medium', auto_increase_solution_order=False)
    t0 = time.time()
    lom.sim.run(name=f'Q_{case}_{tag}', components=comps, open_terminations=opens, box_plus_buffer=False)
    cm = lom.sim.capacitance_matrix
    cm.to_csv(out / 'matrix.csv')
    if str(lom.sim.units).lower() not in ('ff', 'femtofarad', 'femtofarads'):
        raise RuntimeError(f'units {lom.sim.units}')
    names = [n for n in cm.index if 'ground' not in n.lower() and not n.startswith('g_wb')]
    signals = []
    for key in keys:
        hit = [n for n in names if any(k in n.lower() for k in key.split('|'))]
        if len(hit) != 1:
            raise RuntimeError(f'{case}: net for {key!r} not unique in {names}')
        signals.append(hit[0])
    others = [n for n in names if n not in signals]
    sub = cm.loc[signals, signals].to_numpy(float)
    passes = []
    for k, m in lom.sim.capacitance_all_passes.items():
        frame = m if isinstance(m, pd.DataFrame) else pd.DataFrame(m, index=cm.index, columns=cm.columns)
        frame.to_csv(out / f'pass_{k}.csv')
        passes.append(frame.loc[signals, signals].to_numpy(float))
    delta = None
    if len(passes) > 1:
        delta = float(np.max(np.abs(passes[-1] - passes[-2]) / np.maximum(np.abs(passes[-1]), 1e-30)))
    labels = [k.split('|')[0] for k in keys]
    res = dict(case=case, signals=dict(zip(labels, signals)), other_nets=others,
               maxwell_fF=cm.loc[signals + others, signals + others].to_numpy(float).tolist(),
               maxwell_order=signals + others,
               mutual_fF=float(-sub[0, 1]) if len(signals) > 1 else None,
               ground_fF=dict(zip(labels, sub.sum(axis=1).tolist())), max_passes=max_passes, percent_error=pct_error,
               passes_saved=len(passes), last_pass_max_relative_change=delta,
               solver_converged=bool(lom.sim.is_converged), minutes=(time.time() - t0) / 60,
               aedt_design=lom.sim.renderer.pinfo.design.name)
    write(out / 'caps.json', res)
    lom.sim.renderer.pinfo.project._project.Save()
    print(json.dumps({k: v for k, v in res.items() if k != 'maxwell_fF'}), flush=True)
    return res


def hammer_loads(cj_case, ref_case, c_line_fF_per_mm, spacing_mm):
    cj = cj_case['mutual_fF'] - ref_case['mutual_fF']
    ex = [cj_case['ground_fF'][a] - ref_case['ground_fF'][b] - c_line_fF_per_mm * spacing_mm / 2
          for a, b in (('filter', 'line_f'), ('resonator', 'line_r'))]
    return dict(cj=cj, head_excess=float(np.mean(ex)), head_excess_sides=ex, ref_mutual_fF=ref_case['mutual_fF'])


def kron_floating(maxwell, keep, floating):
    """Maxwell matrix (fF) reduced to `keep` with `floating` nets eliminated (net charge 0)."""
    M = np.asarray(maxwell, float)
    K = M[np.ix_(keep, keep)] - M[np.ix_(keep, floating)] @ np.linalg.solve(M[np.ix_(floating, floating)],
                                                                             M[np.ix_(floating, keep)])
    return K


# ---- HFSS ----------------------------------------------------------------------------------------------------------------
def purcell_resonator_pair_components(design, comps):
    return [n for n in comps if type(design.components[n]).__name__ == 'PurcellReadout']


def draw_straps(renderer, design, comps):
    from pyEPR.ansys import parse_units
    names = []
    for c in purcell_resonator_pair_components(design, comps):
        for st in design.components[c].strap_points():
            a = np.radians(st['angle_deg'])
            names.append(renderer.modeler.draw_wirebond(
                pos=parse_units(np.array([st['x'], st['y']])), ori=np.array([np.cos(a), np.sin(a)]),
                width=parse_units(st['span']), height=parse_units(st['span']), z=0,
                wire_diameter=f'{STRAP_DIAMETER_MM}mm', NumSides=6, name='g_wb', material='pec', solve_inside=False))
    return names


def mesh_all(renderer, design, comps):
    id2n = {c.id: n for n, c in design.components.items()}
    objs = [f'{r.name}_{id2n[int(r.component)]}' for kind in ('path', 'poly')
            for r in design.qgeometry.tables[kind].itertuples()
            if not bool(r.subtract) and id2n[int(r.component)] in comps]
    renderer.modeler.mesh_length('metal_mesh', objs, MaxLength=MESH)
    return objs


def run_eigen(folder, d, comps, ports, opens, minf, n_modes, passes=6, tag='E'):
    from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim
    info = {}

    class Eig(EigenmodeSim):
        def _analyze(self):
            info['straps'] = draw_straps(self.renderer, d, comps)
            info['mesh_objects'] = len(mesh_all(self.renderer, d, comps))
            super()._analyze()
    eig = Eig(d, 'hfss')
    eig.renderer._options.update(WB)
    for k, v in dict(min_freq_ghz=minf, n_modes=n_modes, max_delta_f=0.1, max_passes=passes, min_passes=min(2, passes),
                     min_converged=1, pct_refinement=20, basis_order=1, vars=Dict()).items():
        eig.setup[k] = v
    t0 = time.time()
    # The qiskit-metal renderer exports hfss_eig_f_convergence.csv to the current directory: keep it in this run.
    with chdir(Path(folder) / 'raw'):
        eig.run(name=f'{tag}_{Path(folder).name[-8:]}', components=comps, open_terminations=opens, port_list=ports)
    setup = eig.renderer.pinfo.setup
    conv = setup.get_convergence()[0]
    conv.to_csv(Path(folder) / 'raw' / 'eigen_convergence.csv')
    f, im = setup.get_solutions().eigenmodes()
    eig.renderer.pinfo.project._project.Save()
    raw = dict(frequencies_GHz=list(map(float, f)), imag_GHz=list(map(float, im)), passes=int(len(conv)),
               last_delta_f_pct=float(conv.iloc[-1, -1]), minutes=(time.time() - t0) / 60,
               aedt_design=eig.renderer.pinfo.design.name, n_straps=len(info['straps']),
               mesh_objects=info['mesh_objects'], max_passes=passes, min_freq_GHz=minf)
    write(Path(folder) / 'raw' / 'eigen_modes.json', raw)
    return raw


def run_driven(folder, d, comps, ports, opens, f0, band, passes=3, fine_step_ghz=0.00005, tag='S', max_delta_s=0.02,
               fine=True, wide_step_ghz=0.001):
    """Driven S21. fine=False: one Fast wide sweep at wide_step_ghz only; raw['s_file'] names the sweep to analyse."""
    from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
    from scipy.signal import find_peaks
    info = {}

    class Drv(ScatteringImpedanceSim):
        def _analyze(self):
            info['straps'] = draw_straps(self.renderer, d, comps)
            info['mesh_objects'] = len(mesh_all(self.renderer, d, comps))
            super()._analyze()
    sim = Drv(d, 'hfss')
    sim.renderer._options.update(WB)
    lo, hi = band
    sim.setup.update(freq_ghz=f0, max_delta_s=max_delta_s, max_passes=passes, min_passes=1, min_converged=1,
                     pct_refinement=20, basis_order=1, vars=Dict())
    sim.setup.sweep_setup.update(name='Wide', start_ghz=lo, stop_ghz=hi, count=int(round((hi - lo) / wide_step_ghz)) + 1,
                                 type='Fast', save_fields=False)
    t0 = time.time()
    sim.run(name=f'{tag}_{Path(folder).name[-8:]}', components=comps, open_terminations=opens, port_list=ports)
    _, _, S = sim.renderer.get_params(['S11', 'S21'])
    S.to_csv(Path(folder) / 'raw' / 'S_wide.csv')
    conv = sim.renderer.pinfo.setup.get_convergence()[0]
    conv.to_csv(Path(folder) / 'raw' / 'hfss_convergence.csv')
    fw, sw = S.index.to_numpy(float), S['S21'].to_numpy(complex)
    idx, _ = find_peaks(-20 * np.log10(np.maximum(abs(sw), 1e-12)), prominence=1)
    if fine:
        flo, fhi = (max(lo, fw[idx].min() - 0.06), min(hi, fw[idx].max() + 0.06)) if len(idx) else (lo, hi)
        count = min(int(np.ceil((fhi - flo) / fine_step_ghz)) + 1, 25000)
        sim.renderer.add_sweep(setup_name='Setup', name='Fine', start_ghz=flo, stop_ghz=fhi, count=count, type='Fast',
                               save_fields=False)
        sim.renderer.analyze_sweep('Fine', 'Setup')
        _, _, S = sim.renderer.get_params(['S11', 'S21'])
        S.to_csv(Path(folder) / 'raw' / 'S_fine.csv')
    else:
        flo, fhi, count = None, None, 0
    sim.renderer.pinfo.project._project.Save()
    raw = dict(passes=int(len(conv)), max_passes=passes, last_delta_s=float(conv.iloc[-1, -1]),
               minutes=(time.time() - t0) / 60, wide_GHz=[lo, hi], fine_GHz=[flo, fhi], fine_points=count,
               wide_step_GHz=wide_step_ghz, s_file='S_fine.csv' if fine else 'S_wide.csv',
               wide_dips_GHz=fw[idx].tolist(), aedt_design=sim.renderer.pinfo.design.name,
               n_straps=len(info['straps']), mesh_objects=info['mesh_objects'])
    write(Path(folder) / 'raw' / 'driven.json', raw)
    return raw


# ---- analysis ------------------------------------------------------------------------------------------------------
def invert(z1, z2):
    s, p = z1 + z2, z1 * z2
    fr = p.imag / s.imag
    fp = s.real - fr
    j2 = fp * fr - p.real
    return dict(f_p_GHz=float(fp), f_r_GHz=float(fr), J_MHz=float(np.sign(j2) * np.sqrt(abs(j2)) * 1e3),
                kappa_p_MHz=float(2 * abs(s.imag) * 1e3))


def load_s(path):
    S = pd.read_csv(path, index_col=0)
    return S.index.to_numpy(float), S['S21'].map(lambda v: complex(str(v).replace('i', 'j'))).to_numpy()


def channel_windows(f, s, targets):
    """{name: (lo, hi)} S21 windows. targets: {name: f_target_GHz}. If the trace has exactly two dips (>= 3 dB
    prominence) per channel, split halfway between neighbouring dip pairs (robust to common frequency offsets);
    otherwise split halfway between the targets."""
    from scipy.signal import find_peaks
    names = sorted(targets, key=targets.get)
    idx, _ = find_peaks(-20 * np.log10(np.maximum(abs(s), 1e-12)), prominence=3)
    dips = np.sort(f[idx])
    if len(dips) == 2 * len(names):
        cuts = [(dips[2 * i + 1] + dips[2 * i + 2]) / 2 for i in range(len(names) - 1)]
    else:
        fr = [targets[q] for q in names]
        cuts = [(x + y) / 2 for x, y in zip(fr, fr[1:])]
    edges = [f.min()] + cuts + [f.max()]
    return {q: (float(edges[i]), float(edges[i + 1])) for i, q in enumerate(names)}


def s21_channel(f, s, window):
    """Bare f_p, f_r, J (and kappa_p from the pole fit) of the one channel inside `window` (GHz)."""
    from scipy.optimize import least_squares
    from scipy.signal import find_peaks
    m = (f >= window[0]) & (f <= window[1])
    f, s = f[m], s[m]
    # Dips in dB (a low off-resonance baseline at coarse meshes hides them in linear |S21|); keep the two most
    # prominent ones (the channel's zeros) and report any weaker extra dip in the window.
    idx, props = find_peaks(-20 * np.log10(np.maximum(abs(s), 1e-12)), prominence=3)
    if len(idx) < 2:
        raise RuntimeError(f'window {window}: expected two dips, found {len(idx)} at {f[idx]}')
    order = np.argsort(props['prominences'])[::-1]
    extra = [dict(f_GHz=float(f[idx[i]]), prominence_dB=float(props['prominences'][i])) for i in order[2:]]
    idx = np.sort(idx[order[:2]])
    z = []
    for i in idx:
        sl = slice(max(i - 3, 0), i + 4)
        b, a = np.polyfit(f[sl] - f[i], s[sl], 1)
        z.append(f[i] - a / b)
    k = idx[0] + int(np.argmax(abs(s[idx[0]:idx[1]])))
    sl = slice(k - 3, k + 4)
    c = np.polyfit(f[sl] - f[k], abs(s[sl]), 2)
    fr0 = f[k] - c[1] / (2 * c[0])
    fp0 = z[0].real + z[1].real - fr0
    j2 = fp0 * fr0 - (z[0] * z[1]).real
    zeros = dict(zeros_GHz=[[v.real, v.imag] for v in z], f_r_GHz=float(fr0), f_p_GHz=float(fp0),
                 J_MHz=float(np.sqrt(abs(j2)) * 1e3), J2_sign=int(np.sign(j2)))
    cen = f.mean()
    x = (f - cen) * 1e3
    zz = [(v - cen) * 1e3 for v in z]
    num = (x - zz[0]) * (x - zz[1])

    def model(q):
        p1, p2 = q[0] + 1j * q[1], q[2] + 1j * q[3]
        shape = num / ((x - p1) * (x - p2))
        B = np.column_stack([shape, shape * x / 100])
        return B @ np.linalg.lstsq(B, s, rcond=None)[0]

    def res(q):
        r = model(q) - s
        return np.r_[r.real, r.imag]
    best = min((least_squares(res, [zz[0].real, sg * 10, zz[1].real, sg * 10], max_nfev=20000) for sg in (1, -1)),
               key=lambda r: r.cost)
    q = best.x
    poles = [cen + complex(q[0], q[1]) / 1e3, cen + complex(q[2], q[3]) / 1e3]
    if np.mean([p.imag for p in poles]) < 0:
        poles = [p.conjugate() for p in poles]
    rms = float(np.sqrt(np.mean(abs(model(q) - s) ** 2)))
    return dict(window_GHz=list(window), zeros=zeros, extra_dips=extra,
                poles=dict(poles_GHz=[[p.real, p.imag] for p in poles], fit_rms=rms, **invert(*poles)),
                f=f, s=s, fit=model(q))


def eigen_channels(raw, targets, half_window=0.12):
    """Group eigenmodes by channel (two modes nearest each target, within half_window GHz) and invert each pair."""
    z = [complex(a, b) for a, b in zip(raw['frequencies_GHz'], raw['imag_GHz'])]
    out = {}
    for name, f0 in targets.items():
        near = sorted([v for v in z if abs(v.real - f0) < half_window], key=lambda v: abs(v.real - f0))[:2]
        out[name] = invert(*near) if len(near) == 2 else dict(error=f'{len(near)} modes within {half_window} GHz')
        out[name]['modes_GHz'] = [[v.real, v.imag] for v in near]
    return out
