"""Solve + export a driven design whose Mesh Link (eigenmode mesh) was set by hand in the AEDT GUI, then fit.
Scripted EditSetup could not set ImportMesh (the GUI writes PathRelativeTo='SourceProduct'; see run_record.md).

Usage: python driven/run_linked_gui.py <driven design> <eigenmode design> [tag]      (AEDT with Project1 already open)
Checks ImportMesh=true in the project, runs Analyze (setup + sweep), exports S11/S21 to results/<tag>_S.csv,
convergence (solved elements must be ~ the eigenmode mesh), meta.json, project snapshot, then analyze6.
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import os, sys, re, json, time, warnings
warnings.filterwarnings('ignore')
HERE = _S('Design4/simulation')
RES = _S('Design4/results/resonator/legacy'); HFSS = _S('.')
sys.path[:0] = [HERE, _I('Design5/simulation'), _I('Design4/layout')]
import pandas as pd
from pyEPR.ansys import HfssApp
import save_aedt

drv_name, eig_name = sys.argv[1], sys.argv[2]
tag = sys.argv[3] if len(sys.argv) > 3 else 'v6_full'
proj = HfssApp().get_app_desktop().get_active_project()
drv = [d for d in proj.get_designs() if d.name == drv_name][0]
eig = [d for d in proj.get_designs() if d.name == eig_name][0]
proj._project.Save()
txt = open(_relocate(os.path.join(proj._project.GetPath(), proj._project.GetName() + '.aedt')), encoding='utf-8', errors='ignore').read()
blk = txt[txt.find(f"Name='{drv_name}'"):]; i = blk.find("$begin 'MeshLink'")
assert f"ImportMesh=true" in blk[i:i + 120] and f"Design='{eig_name}'" in blk[i:i + 300], 'Mesh Link not set in the driven setup'
econv = eig.get_setup('Setup').get_convergence()[0]; n_eig = int(econv.iloc[-1, 0])
print('Mesh Link confirmed ->', eig_name, '| eigenmode last pass elements', n_eig, flush=True)
t0 = time.time(); drv._design.Analyze('Setup')
minutes = round((time.time() - t0) / 60, 1)
setup = drv.get_setup('Setup'); conv = setup.get_convergence()[0]; conv.to_csv(_relocate(os.path.join(RES, f'{tag}_convergence.csv')))
n_drv = int(conv.iloc[-1, 0]); print(conv.tail(3).to_string(), f'\ndriven solved elements {n_drv} (eigenmode {n_eig}), {minutes} min', flush=True)
freqs, curves = setup.get_sweep().get_network_data(['S11', 'S21'])
S = pd.DataFrame(curves, columns=freqs / 1e9, index=['S11', 'S21']).transpose(); S.index.name = 'f_GHz'
S.to_csv(_relocate(os.path.join(RES, f'{tag}_S.csv'))); econv.to_csv(_relocate(os.path.join(RES, f'{tag}_eig_convergence.csv')))
json.dump(dict(tag=tag, variant='full', note=f'Design6 full channel, eigenmode mesh ({eig_name}, {n_eig} el) imported into {drv_name} (Mesh Link set in AEDT GUI), 1 pass 2nd order, Fast sweep',
               hfss=dict(method='eigenmode mesh -> driven (mesh link, set in GUI)', sweep='6.5-6.95 GHz, 1801 pts, Fast',
                         driven_design=drv_name, eig_design=eig_name, eig=dict(passes=int(len(econv)), elements=n_eig),
                         driven=dict(passes=int(len(conv)), elements=n_drv, minutes=minutes), minutes=minutes)),
          open(_relocate(os.path.join(RES, f'{tag}_meta.json')), 'w'), indent=1)
try: save_aedt.snapshot(_relocate(os.path.join(RES, f'{tag}_aedt')), tag)
except Exception as e: print('AEDT project not saved:', repr(e)[:150])
import analyze6
o, run_dir = analyze6.analyze(tag)
print(json.dumps({k: o[k] for k in ('caps_q3d', 'caps_fit', 'couplings_q3d', 'couplings_fit', 'rms_S21')}, indent=1, default=float))
print('run folder', run_dir, '=== done')
