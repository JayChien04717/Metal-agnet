"""Design6 single channel: driven S21 P1 -> P2 (Interpolating, fine mesh). The Purcell filter peak/dips show as
the dip + double-dip structure of the filter-resonator pair on the feedline.
Usage: python hfss_chan6.py [tag] [f_start f_stop]"""
import os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, HERE)
from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
import design6_layout as d6
tag = sys.argv[1] if len(sys.argv) > 1 else 'c6_v1'
f0, f1 = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (6.3, 7.1)
FO = os.environ.get('FILTER_ONLY') == '1'
d = d6.build_design(filter_only=FO)
sim = ScatteringImpedanceSim(d, 'hfss')
for k, v in dict(freq_ghz=float(os.environ.get('SOL_F', 6.7)),max_delta_s=0.003, max_passes=int(os.environ.get('MAXP', 15)), min_passes=6, min_converged=2,
                 pct_refinement=30, basis_order=int(os.environ.get('BASIS', 2))).items():   # BASIS=-1 -> Mixed Order
    sim.setup[k] = v
sim.setup.sweep_setup.update(name='Sweep', start_ghz=f0, stop_ghz=f1, count=(int(sys.argv[4]) if len(sys.argv) > 4 else int(round((f1 - f0) * 1000 / 0.5)) + 1),
                             type=(sys.argv[5] if len(sys.argv) > 5 else 'Interpolating'), save_fields=False)
comps = ['feed_in', 'cin', 'feed_c', 'n0', 'out_line', 'stub_p', 'cp', 'feed_jf', 'jf', 'jr', 'cj', 'stub_f', 'stub_r', 'filt', 'res', 'res_a']
t0 = time.time()
CAP_MESH = os.environ.get('CAP_MESH', '')          # e.g. '12um': length-based mesh operation on the three capacitor electrodes
if CAP_MESH:
    _orig_add_mesh = sim.renderer.add_mesh
    def _add_mesh_caps():
        _orig_add_mesh()
        objs = ['head_a_cj', 'head_b_cj', 'plate_cin', 'frame_cin', 'plate_cp', 'frame_cp']
        try:
            sim.renderer.modeler.mesh_length('cap_mesh', objs, MaxLength=CAP_MESH)
            print('MESH op cap_mesh', CAP_MESH, objs)
        except Exception as e:
            print('MESH op failed:', repr(e)[:200])
    sim.renderer.add_mesh = _add_mesh_caps
if FO: comps = ['feed_in', 'cin', 'feed_c', 'n0', 'out_line', 'stub_p', 'cp', 'filt']
sim.run(name=tag, components=comps, open_terminations=([] if FO else [('res_a', 'start')]), port_list=[('feed_in', 'start', 50), ('out_line', 'end', 50)])
try:      # how many adaptive passes were really used, and the last Delta S (to decide whether max_passes can be lowered)
    conv = sim.renderer.pinfo.setup.get_convergence()[0]
    conv.to_csv(os.path.join(HERE, 'results', f'{tag}_convergence.csv')); print('CONVERGENCE (passes used, last rows):'); print(conv.tail(4).to_string())
except Exception as e:
    print('convergence not read:', repr(e)[:150])
_, _, S = sim.renderer.get_params(['S11', 'S21'])
S.index.name = 'f_GHz'; S.to_csv(os.path.join(HERE, 'results', f'{tag}_S.csv'))
import run_record; run_record.save(tag, d6, FO, os.environ.get('RUN_NOTE', ''))
f = np.asarray(S.index, float); s21 = 20 * np.log10(np.abs(S['S21'].values)); s11 = 20 * np.log10(np.abs(S['S11'].values))
from scipy.signal import find_peaks
print(f'=== {tag} ({(time.time() - t0) / 60:.1f} min) ===')
idx, pr = find_peaks(-s21, prominence=0.5)
for i, p in zip(idx, pr['prominences']):
    print(f'  dip at {f[i]:.4f} GHz: S21 {s21[i]:.1f} dB (prominence {p:.1f})')
i = int(np.argmax(s21)); print(f'max S21 {s21[i]:.2f} dB at {f[i]:.4f} GHz; min S21 {s21.min():.1f} dB')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 4.5)); ax.plot(f, s21, label='S21'); ax.plot(f, s11, label='S11')
ax.axvline(d6.F0, color='gray', ls=':'); ax.set_xlabel('f (GHz)'); ax.set_ylabel('dB'); ax.legend(); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'results', f'{tag}_S.png'), dpi=120)
import theory_abcd; print('run folder', theory_abcd.archive(tag))   # move into results/<time>_<tag>/ + theory curve + overlay
try: sim.close()
except Exception: pass
