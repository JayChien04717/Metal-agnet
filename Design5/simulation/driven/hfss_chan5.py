"""Design5 full channel: HFSS driven modal S11 / S21 (P1 -> P2), then analyze5 (capacitances, S21 / S11, plot raw / Q3D model / fitted model).
Defaults from SQL HFSS/skill/skill.md sec. 10: Fast sweep (Interpolating can miss the narrow dips), basis 2,
max_delta_s 0.003, min_passes 6, max passes 16.

Usage: python hfss_chan5.py <tag> [f_start f_stop count type]     (env SOL_F, MAXP, BASIS, RUN_NOTE)
"""
import os, sys, json, time, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # simulation/ (this script: simulation/<kind>/)
import sys as _sys; _sys.path[:0] = [os.path.join(HERE, _k) for _k in ('driven', 'eigenmode', 'q3d')] + [HERE]
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, HERE)
from qiskit_metal.analyses.simulation import ScatteringImpedanceSim
import design5_layout as d5

tag = sys.argv[1] if len(sys.argv) > 1 else 'c5_W1'
f0, f1 = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (6.2, 7.2)
count = int(sys.argv[4]) if len(sys.argv) > 4 else int(round((f1 - f0) * 1000 / 0.25)) + 1      # 0.25 MHz step
stype = sys.argv[5] if len(sys.argv) > 5 else 'Fast'
setup = dict(freq_ghz=float(os.environ.get('SOL_F', 6.71)), max_delta_s=0.003, max_passes=int(os.environ.get('MAXP', 16)),
             min_passes=6, min_converged=2, pct_refinement=30, basis_order=int(os.environ.get('BASIS', 2)))
import hfss_dirs
temp_dirs = hfss_dirs.use_d_drive()          # HFSS solver scratch on D: (C: fills up -> COM error at Analyze)
d = d5.build_design()
sim = ScatteringImpedanceSim(d, 'hfss')
for k, v in setup.items():
    sim.setup[k] = v
sim.setup.sweep_setup.update(name='Sweep', start_ghz=f0, stop_ghz=f1, count=count, type=stype, save_fields=False)
comps = ['feed_in', 'cin', 'feed_c', 'n0', 'out_line', 'stub_p', 'cp', 'feed_jf', 'jf', 'jr', 'cj', 'stub_f', 'stub_r', 'filt', 'res', 'res_a']
t0 = time.time()
# the resonator open end must be listed (ResEnd is pin-only; unlisted ends are shorted = the two lambda/4 shorts)
sim.run(name=tag, components=comps, open_terminations=[('res_a', 'start')], port_list=[('feed_in', 'start', 50), ('out_line', 'end', 50)])
conv_info = {}
try:
    conv = sim.renderer.pinfo.setup.get_convergence()[0]
    conv.to_csv(os.path.join(HERE, 'results', f'{tag}_convergence.csv')); print(conv.tail(4).to_string())
    conv_info = dict(passes=int(len(conv)), last_dS=float(conv.iloc[-1, -1]))
except Exception as e:
    print('convergence not read:', repr(e)[:150])
_, _, S = sim.renderer.get_params(['S11', 'S21'])
S.index.name = 'f_GHz'; S.to_csv(os.path.join(HERE, 'results', f'{tag}_S.csv'))
minutes = round((time.time() - t0) / 60, 1)
json.dump(dict(tag=tag, note=os.environ.get('RUN_NOTE', ''), hfss=dict(sweep=f'{f0}-{f1} GHz, {count} pts, {stype}', **setup, **conv_info, minutes=minutes,
                         temp_dir=temp_dirs['after']),
               layout=dict(variant=d5.VARIANT, L_FILTER=d5.L_FILTER, L_RES_TOTAL=d5.L_RES_TOTAL, D_CJ=d5.D_CJ, CP_TYPE=d5.CP_TYPE,
                           Cp_plate_um=d5.ENC_P['plate_width'], Cj_gap_um=d5.CJ_HG)),
          open(os.path.join(HERE, 'results', f'{tag}_meta.json'), 'w'), indent=1)
print(f'=== {tag} HFSS done ({minutes} min) ===')
try:                                                   # keep the AEDT project + solved results with the run
    import save_aedt; save_aedt.save_active(os.path.join(HERE, 'results', f'{tag}_aedt'), tag, close=True)
except Exception as e:
    print('AEDT project not saved:', repr(e)[:150])
try: sim.close()
except Exception: pass
import analyze5
o, run_dir = analyze5.analyze(tag)
print(json.dumps(dict(caps_fit=o['caps_fit'], dips=o['dips'], rms=o['rms_S21']), indent=1, default=float)); print('run folder', run_dir)
