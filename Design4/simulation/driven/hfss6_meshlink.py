"""Design6 full channel: adapt the mesh with an EIGENMODE solve, import that mesh into the DRIVEN modal design.

1. eigenmode design <tag>_eig_hfss: same geometry (the two ports drawn as 50 Ohm RLC sheets = the same objects as the
   driven lumped ports, same wire bonds), initial length-based mesh (10 um lines / 5 um capacitor electrodes),
   adaptive passes converged on the two mode frequencies (~6.69 / 6.73 GHz): max delta f MAX_DF % (default 0.1 %,
   Jay 2026-11-09), 3-15 passes.
2. driven design <tag>_hfss: setup with Mesh Link -> <tag>_eig_hfss : Setup : LastAdaptive (+ its own mesh ops),
   a few more adaptive passes on S (DRV_MAXP), then the Fast sweep.
3. analyze6 (capacitances fitted to S21 vs Q3D, kappa_p, J) and the AEDT project saved with the run.

Usage: python driven/hfss6_meshlink.py [tag] [f_start f_stop]      (env MAX_DF, EIG_MAXP, DRV_MAXP, MESH_LEN, MESH_CAP, RUN_NOTE)
       SKIP_EIG=1: reuse <tag>_eig_hfss already solved in the open AEDT project (EIG_NOTE describes how it was solved)
"""
import os, sys, json, time, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # simulation/ (this script: simulation/driven/)
RES = os.path.join(HERE, 'results'); HFSS = os.path.dirname(os.path.dirname(HERE))
sys.path[:0] = [os.path.join(HERE, k) for k in ('driven', 'eigenmode', 'q3d')] + [HERE]
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, os.path.join(HFSS, 'Design5', 'simulation'))
from qiskit_metal import Dict
from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim
import verify_layout as V, hfss_dirs

tag = sys.argv[1] if len(sys.argv) > 1 else 'v6_full'
f0, f1 = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (6.50, 6.95)
count = int(round((f1 - f0) * 1000 / 0.25)) + 1
MESH_LEN, MESH_CAP = os.environ.get('MESH_LEN', '10um'), os.environ.get('MESH_CAP', '5um')
EIG = dict(min_freq_ghz=float(os.environ.get('EIG_MINF', 6.5)), n_modes=int(os.environ.get('EIG_NMODES', 2)), max_delta_f=float(os.environ.get('MAX_DF', 0.1)), max_passes=int(os.environ.get('EIG_MAXP', 15)),
           min_passes=3, min_converged=2, pct_refinement=20, basis_order=int(os.environ.get('EIG_BASIS', 1)), vars=Dict())
# eigenmode only builds the mesh -> first-order basis (2nd order ran out of the 32 GB at pass 1, 2026-11-09);
# the driven solve on the imported mesh stays second order
DRV = dict(freq_ghz=6.71, max_delta_s=0.003, max_passes=int(os.environ.get('DRV_MAXP', 1)), min_passes=1, min_converged=1,
           # DRV_MAXP=1 (default, Jay 2026-11-09): solve once on the imported eigenmode mesh, no further refinement
           pct_refinement=10, basis_order=2, vars=Dict())
CAP_COMPS = (f'cp_ch{V.CH}', f'cj_ch{V.CH}')
PORTS = [('feed_in', 'start', 50), ('feed_out', 'end', 50)]
INFO = {}
os.makedirs(RES, exist_ok=True)
temp_dirs = hfss_dirs.use_d_drive()
import save_aedt; save_aedt.use_project(save_aedt.project_file(HERE))   # one AEDT project per Design, never closed
d = V.build('full')
comps = [n for n, c in d.components.items() if type(c).__name__ not in ('ResEnd', 'LaunchpadWirebond')]
opens = [(f'res_a_ch{V.CH}', 'start')]


def mesh_ops(renderer):
    lines, caps = [], []
    id2n = {c.id: n for n, c in d.components.items()}
    for kind in ('path', 'poly'):
        for r in d.qgeometry.tables[kind].itertuples():
            comp = id2n[int(r.component)]
            if not bool(r.subtract) and comp in comps:
                (caps if comp in CAP_COMPS else lines).append(f'{r.name}_{comp}')
    renderer.modeler.mesh_length('cpw_line_mesh', lines, MaxLength=MESH_LEN)
    renderer.modeler.mesh_length('cap_mesh', caps, MaxLength=MESH_CAP)
    print(f'length-based mesh: {MESH_LEN} on {len(lines)} lines, {MESH_CAP} on {len(caps)} capacitor objects', flush=True)


class Eig(EigenmodeSim):
    def _analyze(self):
        mesh_ops(self.renderer)
        super()._analyze()


class DrvLinked(ScatteringImpedanceSim):
    def _analyze(self):
        self.sim_setup_name, self.sweep_name = self.renderer.initialize_drivenmodal(**self.setup)
        mod = self.renderer.pinfo.design._design.GetModule('AnalysisSetup')
        s = self.setup
        mod.EditSetup(self.sim_setup_name, [
            f'NAME:{self.sim_setup_name}', 'SolveType:=', 'Single', 'Frequency:=', f'{s.freq_ghz}GHz',
            'MaxDeltaS:=', s.max_delta_s, 'UseMatrixConv:=', False, 'MaximumPasses:=', s.max_passes,
            'MinimumPasses:=', s.min_passes, 'MinimumConvergedPasses:=', s.min_converged, 'PercentRefinement:=', s.pct_refinement,
            'IsEnabled:=', True,
            'MeshLink:=', ['NAME:MeshLink', 'ImportMesh:=', True, 'Project:=', 'This Project*', 'Product:=', 'HFSS', 'Design:=', EIG_DESIGN,
                           'Soln:=', 'Setup : LastAdaptive', ['NAME:Params'], 'ForceSourceToSolve:=', True,
                           'PreservePartnerSoln:=', True, 'PathRelativeTo:=', 'TargetProject', 'ApplyMeshOp:=', True],
            'BasisOrder:=', s.basis_order, 'DoLambdaRefine:=', True, 'DoMaterialLambda:=', True, 'SetLambdaTarget:=', False,
            'Target:=', 0.3333, 'UseMaxTetIncrease:=', False, 'PortAccuracy:=', 2, 'UseABCOnPort:=', False,
            'SetPortMinMaxTri:=', False, 'UseDomains:=', False, 'UseIterativeSolver:=', False, 'SaveRadFieldsOnly:=', False,
            'SaveAnyFields:=', True, 'IESolverType:=', 'Auto', 'LambdaTargetForIESolver:=', 0.15, 'UseDefaultLambdaTgtForIESolver:=', True])
        # the nested list is silently ignored unless it starts with 'NAME:MeshLink' (run 1601: ImportMesh stayed false, 10,082 elements)
        pr = self.renderer.pinfo.project._project; pr.Save()
        txt = open(os.path.join(pr.GetPath(), pr.GetName() + '.aedt'), encoding='utf-8', errors='ignore').read()
        blk = txt[txt.find(f"Name='{self.renderer.pinfo.design.name}'"):]
        i = blk.find("$begin 'MeshLink'")
        assert 'ImportMesh=true' in blk[i:i + 200], 'Mesh Link not registered in the driven setup'
        print(f'driven setup: mesh imported from {EIG_DESIGN} : Setup : LastAdaptive (ImportMesh=true confirmed in project)', flush=True)
        self.renderer.analyze_sweep(self.sweep_name, self.sim_setup_name)


t0 = time.time()
EIG_DESIGN = os.environ.get('EIG_DESIGN') or (save_aedt.design_name(f'{tag}_eig') + '_hfss')
if os.environ.get('SKIP_EIG') == '1':
    # ---- 1'. reuse an eigenmode design already solved in the open project (e.g. finished by hand in AEDT) ----
    from pyEPR.ansys import HfssApp
    proj = HfssApp().get_app_desktop().get_active_project()
    ed = [x for x in proj.get_designs() if x.name == EIG_DESIGN][0]
    if os.environ.get('CONTINUE_EIG') == '1':       # add adaptive passes from the last one until its own criterion
        te = time.time(); print(f'continuing {EIG_DESIGN} : Setup (Delta F', ed._design.GetPropertyValue('HfssTab', 'AnalysisSetup:Setup', 'Delta F'), '%)', flush=True)
        ed._design.Analyze('Setup'); INFO['eig_continue_min'] = round((time.time() - te) / 60, 1)
    conv = ed.get_setup('Setup').get_convergence()[0]; conv.to_csv(os.path.join(RES, f'{tag}_eig_convergence.csv'))
    print(f'reusing solved {EIG_DESIGN}:', conv.tail(4).to_string(), sep='\n', flush=True)
    INFO['eig'] = dict(passes=int(len(conv)), elements=int(conv.iloc[-1, 0]), last_df_pct=float(conv.iloc[-1, -1]),
                       note=os.environ.get('EIG_NOTE', 'solved before this script (reused)'))
else:
    # ---- 1. eigenmode: adapt the mesh ----
    eig = Eig(d, 'hfss')
    for k, v in EIG.items():
        eig.setup[k] = v
    eig.run(name=EIG_DESIGN[:-len('_hfss')], components=comps, open_terminations=opens, port_list=PORTS)
    try:
        conv = eig.renderer.pinfo.setup.get_convergence()[0]; conv.to_csv(os.path.join(RES, f'{tag}_eig_convergence.csv'))
        print(conv.tail(4).to_string(), flush=True)
        INFO['eig'] = dict(passes=int(len(conv)), elements=int(conv.iloc[-1, 0]), minutes=round((time.time() - t0) / 60, 1))
    except Exception as e:
        print('eigen convergence not read:', repr(e)[:150])
if os.environ.get('EIG_ONLY') == '1':          # eigenmode diagnostic only: print the mode frequencies / Q and stop
    try:
        fr, qq = eig.renderer.pinfo.setup.get_solutions().eigenmodes()
        for i, (f_, g_) in enumerate(zip(fr, qq)): print(f'MODE {i + 1}: f = {f_:.4f} GHz, imaginary part = {g_:.5f} GHz (Q = f/2/imag = {f_ / (2 * g_):.0f})', flush=True)
    except Exception as e: print('modes not read:', repr(e)[:150])
    print('=== eig only done', flush=True); sys.exit(0)
t1 = time.time()
# ---- 2. driven with the imported mesh ----
sim = DrvLinked(d, 'hfss')
for k, v in DRV.items():
    sim.setup[k] = v
sim.setup.sweep_setup.update(name='Sweep', start_ghz=f0, stop_ghz=f1, count=count, type='Fast', save_fields=False)
sim.run(name=save_aedt.design_name(tag), components=comps, open_terminations=opens, port_list=PORTS)
try:
    conv = sim.renderer.pinfo.setup.get_convergence()[0]; conv.to_csv(os.path.join(RES, f'{tag}_convergence.csv'))
    print(conv.tail(4).to_string(), flush=True)
    INFO['driven'] = dict(passes=int(len(conv)), elements=int(conv.iloc[-1, 0]), last_dS=float(conv.iloc[-1, -1]),
                          minutes=round((time.time() - t1) / 60, 1))
except Exception as e:
    print('driven convergence not read:', repr(e)[:150])
_, _, S = sim.renderer.get_params(['S11', 'S21'])
S.index.name = 'f_GHz'; S.to_csv(os.path.join(RES, f'{tag}_S.csv'))
minutes = round((time.time() - t0) / 60, 1)
json.dump(dict(tag=tag, variant='full', note=os.environ.get('RUN_NOTE', ''),
               hfss=dict(method='eigenmode mesh -> driven (mesh link)', sweep=f'{f0}-{f1} GHz, {count} pts, Fast',
                         eig_setup={k: v for k, v in EIG.items() if k != 'vars'}, driven_setup={k: v for k, v in DRV.items() if k != 'vars'},
                         mesh_lines=MESH_LEN, mesh_caps=MESH_CAP, **INFO, minutes=minutes, temp_dir=temp_dirs['after']),
               layout=V.info()), open(os.path.join(RES, f'{tag}_meta.json'), 'w'), indent=1)
print(f'=== {tag} done ({minutes} min) ===', flush=True)
try:
    import save_aedt; save_aedt.snapshot(os.path.join(RES, f'{tag}_aedt'), tag)
except Exception as e:
    print('AEDT project not saved:', repr(e)[:150])
try: sim.close()
except Exception: pass
import analyze6
o, run_dir = analyze6.analyze(tag)
print(json.dumps({k: o[k] for k in ('caps_q3d', 'caps_fit', 'couplings_q3d', 'couplings_fit', 'rms_S21')}, indent=1, default=float))
print('run folder', run_dir)
