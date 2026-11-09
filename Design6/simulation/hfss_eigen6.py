"""Design6 eigenmode check: resonance frequencies of the coupled filter + resonator (no ports, no sweep).
The feed line is cut at C_p (only the short stub_p/n0 junction remains, open), so the modes are the bare
filter / resonator modes hybridised by C_j. Usage: python hfss_eigen6.py [tag]   (env CAP_MESH=12um, BASIS=-1, MINF=5.8)"""
import os, sys, time, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, HERE)
from qiskit_metal.analyses.simulation import EigenmodeSim
import design6_layout as d6
tag = sys.argv[1] if len(sys.argv) > 1 else 'e6_1'
d = d6.build_design()
sim = EigenmodeSim(d, 'hfss')
for k, v in dict(min_freq_ghz=float(os.environ.get('MINF', 5.8)), n_modes=4, max_delta_f=0.1, max_passes=15, min_passes=4,
                 min_converged=2, pct_refinement=30, basis_order=int(os.environ.get('BASIS', 2))).items():
    sim.setup[k] = v
CAP_MESH = os.environ.get('CAP_MESH', '')
if CAP_MESH:
    _orig = sim.renderer.add_mesh
    def _add_mesh_caps():
        _orig()
        objs = ['head_a_cj', 'head_b_cj', 'plate_cin', 'frame_cin', 'plate_cp', 'frame_cp']
        try:
            sim.renderer.modeler.mesh_length('cap_mesh', objs, MaxLength=CAP_MESH); print('MESH op cap_mesh', CAP_MESH)
        except Exception as e:
            print('MESH op failed:', repr(e)[:200])
    sim.renderer.add_mesh = _add_mesh_caps
comps = ['n0', 'stub_p', 'cp', 'feed_jf', 'jf', 'jr', 'cj', 'stub_f', 'stub_r', 'filt', 'res', 'res_a']
t0 = time.time()
sim.run(name=tag, components=comps, open_terminations=[('n0', 'prime_start'), ('n0', 'prime_end'), ('res_a', 'start')])
print(f'=== {tag} ({(time.time() - t0) / 60:.1f} min) ===')
try:
    cf = sim.convergence_f; ct = sim.convergence_t
    print('CONVERGENCE f (GHz per pass):'); print(cf.to_string()); print(ct.to_string())
    cf.to_csv(os.path.join(HERE, 'results', f'{tag}_eigen_f.csv')); ct.to_csv(os.path.join(HERE, 'results', f'{tag}_eigen_t.csv'))
except Exception as e:
    print('convergence not read:', repr(e)[:200])
try: sim.close()
except Exception: pass
