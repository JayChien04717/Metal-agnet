
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

import numpy as np
import pandas as pd
from qiskit_metal import Dict
from qiskit_metal.analyses.simulation.scattering_impedance import ScatteringImpedanceSim
sys.path.insert(0,_I('Design5/simulation'))
import save_aedt
from layout3_full import build_design

src=_P('Design2/gds/layout3 ftf.gds')
d=build_design(src,readout_ports=True)
save_aedt.use_project(Path(os.environ.get('AEDT_PROJECT', _S('Design2/simulation/aedt/D2_Layout3.aedt'))))
sim=ScatteringImpedanceSim(d,'hfss')
sim.setup.update(dict(freq_ghz=6., max_delta_s=.02,max_passes=12,min_passes=3,
    min_converged=2,pct_refinement=20,basis_order=1,
    sweep_setup=Dict(name='Sweep',start_ghz=4.8,stop_ghz=7.5,count=1081,
                     step_ghz=None,type='Interpolating',save_fields=False)))
sim.run_sim(name='D2_Layout3_Readout',components=list(d.components),open_terminations=[],
            port_list=[('feed','left',50),('feed','right',50)],box_plus_buffer=False)
f,params,_=sim.renderer.get_params(['S11','S21','S12','S22'])
table={'frequency_Hz':f}
for name,z in zip(('S11','S21','S12','S22'),params):
    table[name+'_re']=np.real(z)
    table[name+'_im']=np.imag(z)
pd.DataFrame(table).to_csv(_P('Design2/results/general/legacy/layout3_full/readout_sparams.csv'),index=False)
conv,_=sim.renderer.pinfo.setup.get_convergence()
conv.to_csv(_P('Design2/results/general/legacy/layout3_full/readout_hfss_convergence.csv'))
desk=save_aedt._desktop()
desk.get_active_project()._project.Save()
(_P('Design2/results/general/legacy/layout3_full/readout_hfss_messages.log')).write_text('\n'.join(desk._desktop.GetMessages(
    'D2_Layout3','D2_Layout3_Readout_hfss',0)),encoding='utf8')
