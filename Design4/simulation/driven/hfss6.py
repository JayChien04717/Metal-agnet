"""HFSS driven modal of a Design6 verification layout (no C_in, matched feed line), then analyze6.

Usage: python driven/hfss6.py <filter|full> <tag> [f_start f_stop]
       env SOL_F, MAXP, MINP, PCT_REF, BASIS, MESH_LEN, MESH_CAP, RUN_NOTE
Defaults: Fast sweep, 0.25 MHz step, basis 2, max dS 0.003, 3-15 passes, 20 % refinement, adapt at 6.71 GHz,
initial mesh 10 um on the lines / 5 um on the capacitor electrodes, HFSS wire bonds on every CPW.
The full channel alone is enough: with no C_in its S21 fixes dL_filter, dL_res, C_p and C_j together (synthetic test).
"""
import os, sys, json, time, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES = os.path.join(HERE, 'results')   # simulation/ (this script: simulation/<kind>/)
import sys as _sys; _sys.path[:0] = [os.path.join(HERE, _k) for _k in ('driven', 'eigenmode', 'q3d')] + [HERE]
HFSS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, HERE)
sys.path.insert(2, os.path.join(HFSS, 'Design5', 'simulation'))

variant, tag = sys.argv[1], sys.argv[2]
f0, f1 = (float(sys.argv[3]), float(sys.argv[4])) if len(sys.argv) > 4 else (6.50, 6.95)
count = int(round((f1 - f0) * 1000 / 0.25)) + 1
setup = dict(freq_ghz=float(os.environ.get('SOL_F', 6.71)), max_delta_s=0.003, max_passes=int(os.environ.get('MAXP', 15)),
             min_passes=int(os.environ.get('MINP', 3)), min_converged=2, pct_refinement=int(os.environ.get('PCT_REF', 20)),
             basis_order=int(os.environ.get('BASIS', 2)))     # fine initial mesh -> fewer passes, gentler refinement
import verify_layout as V, hfss_dirs
from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
os.makedirs(RES, exist_ok=True)
# initial mesh (length based) so the adaptive passes start from a mesh that already resolves the lines:
#   CPW centre conductors MESH_LEN (10 um = trace width), capacitor electrodes (C_p hugger, C_j hammer) MESH_CAP (5 um).
# 5 um on every line (2026-11-09 v6_filter) started at 428k elements and filled 17.6 of 32 GB by pass 3 -> too heavy.
MESH_LEN = os.environ.get('MESH_LEN', '10um')
MESH_CAP = os.environ.get('MESH_CAP', '5um')
CAP_COMPS = (f'cp_ch{V.CH}', f'cj_ch{V.CH}')


class SimMesh(ScatteringImpedanceSim):
    def _analyze(self):
        lines, caps = [], []
        id2n = {c.id: n for n, c in d.components.items()}
        for kind in ('path', 'poly'):
            for r in d.qgeometry.tables[kind].itertuples():
                comp = id2n[int(r.component)]
                if not bool(r.subtract) and comp in comps:
                    (caps if comp in CAP_COMPS else lines).append(f'{r.name}_{comp}')
        try:
            existing = set(self.renderer.modeler.get_objects_in_group('Sheets'))
        except Exception:
            existing = set(lines + caps)
        lines = sorted(set(lines) & existing) or lines; caps = sorted(set(caps) & existing) or caps
        print(f'length-based mesh: {MESH_LEN} on {len(lines)} line objects, {MESH_CAP} on {len(caps)} capacitor objects', flush=True)
        if lines: self.renderer.modeler.mesh_length('cpw_line_mesh', lines, MaxLength=MESH_LEN)
        if caps: self.renderer.modeler.mesh_length('cap_mesh', caps, MaxLength=MESH_CAP)
        MESH_INFO.update(mesh_lines=MESH_LEN, mesh_caps=MESH_CAP, mesh_objects=len(lines) + len(caps))
        super()._analyze()


MESH_INFO = {}
temp_dirs = hfss_dirs.use_d_drive()
import save_aedt; save_aedt.use_project(save_aedt.project_file(HERE))   # one AEDT project per Design, never closed
d = V.build(variant)
sim = SimMesh(d, 'hfss')
for k, v in setup.items():
    sim.setup[k] = v
sim.setup.sweep_setup.update(name='Sweep', start_ghz=f0, stop_ghz=f1, count=count, type='Fast', save_fields=False)
comps = [n for n, c in d.components.items() if type(c).__name__ not in ('ResEnd', 'LaunchpadWirebond')]
opens = [(f'res_a_ch{V.CH}', 'start')] if variant == 'full' else []          # resonator open end; other ends shorted
t0 = time.time()
sim.run(name=save_aedt.design_name(tag), components=comps, open_terminations=opens, port_list=[('feed_in', 'start', 50), ('feed_out', 'end', 50)])
conv_info = {}
try:
    conv = sim.renderer.pinfo.setup.get_convergence()[0]
    conv.to_csv(os.path.join(RES, f'{tag}_convergence.csv')); print(conv.tail(4).to_string())
    conv_info = dict(passes=int(len(conv)), last_dS=float(conv.iloc[-1, -1]))
except Exception as e:
    print('convergence not read:', repr(e)[:150])
_, _, S = sim.renderer.get_params(['S11', 'S21'])
S.index.name = 'f_GHz'; S.to_csv(os.path.join(RES, f'{tag}_S.csv'))
minutes = round((time.time() - t0) / 60, 1)
json.dump(dict(tag=tag, variant=variant, note=os.environ.get('RUN_NOTE', ''),
               hfss=dict(sweep=f'{f0}-{f1} GHz, {count} pts, Fast', **setup, **conv_info, **MESH_INFO, minutes=minutes, temp_dir=temp_dirs['after']),
               layout=V.info()), open(os.path.join(RES, f'{tag}_meta.json'), 'w'), indent=1)
print(f'=== {tag} HFSS done ({minutes} min) ===')
try:                                                   # keep the AEDT project + solved results with the run
    import save_aedt; save_aedt.snapshot(os.path.join(RES, f'{tag}_aedt'), tag)
except Exception as e:
    print('AEDT project not saved:', repr(e)[:150])
try: sim.close()
except Exception: pass
import analyze6
o, run_dir = analyze6.analyze(tag)
print(json.dumps({k: o[k] for k in ('caps_q3d', 'caps_fit', 'couplings_q3d', 'couplings_fit', 'rms_S21')}, indent=1, default=float))
print('run folder', run_dir)
