
# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
from pathlib import Path
import os
import sys

design_dir = _P('Design2')
root = _P('.')
out = _P('Design2/results/general/legacy/layout3_full')
out.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, _I('Design2/layout'))

import pandas as pd
sys.path.insert(0,_I('Design5/simulation'))
import save_aedt
from layout3_full import build_design
from qiskit_metal.analyses.quantization import LOManalysis
src=_P('Design2/gds/layout3 ftf.gds')
d=build_design(src)
save_aedt.use_project(Path(os.environ.get('AEDT_PROJECT', _S('Design2/simulation/aedt/D2_Layout3.aedt'))))
sim=LOManalysis(d,'q3d').sim
sim.setup.update(dict(max_passes=8,min_passes=2,percent_error=.5,freq_ghz=6.,
   solution_order='Medium',auto_increase_solution_order=False))
sim.run(name='D2_Layout3',open_terminations=[],box_plus_buffer=False)
sim.capacitance_matrix.to_csv(_P('Design2/results/general/legacy/layout3_full/capacitance_fF.csv'))
conv,_=sim.renderer.pinfo.setup.get_convergence()
conv.to_csv(_P('Design2/results/general/legacy/layout3_full/q3d_convergence.csv'))
save_aedt._desktop().get_active_project()._project.Save()
print(sim.capacitance_matrix.to_string(),flush=True)
print(conv.to_string(),flush=True)
