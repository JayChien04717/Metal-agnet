from pathlib import Path
import os
import sys

design_dir = Path(__file__).resolve().parents[1]
root = design_dir.parent
out = design_dir / 'simulation/results/layout3_full'
out.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(design_dir / 'pattern design'))

import pandas as pd
sys.path.insert(0,str(root/'Design5/simulation'))
import save_aedt
from layout3_full import build_design
from qiskit_metal.analyses.quantization import LOManalysis
src=design_dir/'gds/layout3 ftf.gds'
d=build_design(src)
save_aedt.use_project(Path(os.environ.get('AEDT_PROJECT', str(design_dir / 'simulation/aedt/D2_Layout3.aedt'))))
sim=LOManalysis(d,'q3d').sim
sim.setup.update(dict(max_passes=8,min_passes=2,percent_error=.5,freq_ghz=6.,
   solution_order='Medium',auto_increase_solution_order=False))
sim.run(name='D2_Layout3',open_terminations=[],box_plus_buffer=False)
sim.capacitance_matrix.to_csv(out/'capacitance_fF.csv')
conv,_=sim.renderer.pinfo.setup.get_convergence()
conv.to_csv(out/'q3d_convergence.csv')
save_aedt._desktop().get_active_project()._project.Save()
print(sim.capacitance_matrix.to_string(),flush=True)
print(conv.to_string(),flush=True)
