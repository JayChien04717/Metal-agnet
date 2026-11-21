
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

design_dir = _P('Design2_FTF_Andersen')
root = _P('.')
out = _P('Design2_FTF_Andersen/results/general/legacy/layout3_full')
out.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, _I('Design2_FTF_Andersen/layout'))

from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim
sys.path.insert(0,_I('Design5_6_purcell-res/simulation'))
import save_aedt
from layout3_full import build_linear_design

src=_P('Design2_FTF_Andersen/gds/layout3 ftf.gds')
d = build_linear_design(src, readout_ports=True)
save_aedt.use_project(Path(os.environ.get('AEDT_PROJECT', _S('Design2_FTF_Andersen/simulation/aedt/D2_Layout3.aedt'))))
os.chdir(out)
for side in ('left','right'):
    # Remove the feed conductor and opposite resonator to identify each resonance.
    keep = [name for name in d.components if name not in ('feed','res_'+('right' if side=='left' else 'left'))]
    sim=EigenmodeSim(d,'hfss')
    sim.setup.update(dict(min_freq_ghz=5.5,n_modes=4,max_delta_f=.3,max_passes=6,
                          min_passes=3,min_converged=2,pct_refinement=20,basis_order=1))
    name='D2_Layout3_Res_'+side
    sim.run_sim(name=name,components=keep,open_terminations=[],box_plus_buffer=False)
    sim.convergence_f.to_csv(_relocate(out/f'eigen_readout_{side}_freqs.csv'))
    sim.convergence_t.to_csv(_relocate(out/f'eigen_readout_{side}_mesh.csv'))
    save_aedt._desktop().get_active_project()._project.Save()
