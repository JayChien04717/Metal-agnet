from pathlib import Path
import os
import sys

design_dir = Path(__file__).resolve().parents[1]
root = design_dir.parent
out = design_dir / 'simulation/results/layout3_full'
out.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(design_dir / 'pattern design'))

import numpy as np
import pandas as pd
from qiskit_metal import Dict
from qiskit_metal.analyses.simulation.scattering_impedance import ScatteringImpedanceSim
sys.path.insert(0,str(root/'Design5/simulation'))
import save_aedt
from layout3_full import build_linear_design

src=design_dir/'gds/layout3 ftf.gds'
d = build_linear_design(src, readout_ports=True)
save_aedt.use_project(Path(os.environ.get('AEDT_PROJECT', str(design_dir / 'simulation/aedt/D2_Layout3.aedt'))))
sim=ScatteringImpedanceSim(d,'hfss')
sim.setup.update(dict(freq_ghz=6.6, max_delta_s=.02,max_passes=12,min_passes=8,
    min_converged=2,pct_refinement=20,basis_order=1,
    sweep_setup=Dict(name='Sweep',start_ghz=6.25,stop_ghz=6.85,count=2401,
                     step_ghz=None,type='Interpolating',save_fields=False)))
sim.run_sim(name='D2_Layout3_Readout_Linear',components=list(d.components),open_terminations=[],
            port_list=[('feed','left',50),('feed','right',50)],box_plus_buffer=False)
f,params,_=sim.renderer.get_params(['S11','S21','S12','S22'])
table={'frequency_Hz':f}
for name,z in zip(('S11','S21','S12','S22'),params):
    table[name+'_re']=np.real(z)
    table[name+'_im']=np.imag(z)
pd.DataFrame(table).to_csv(out/'linear_readout_sparams.csv',index=False)
conv,_=sim.renderer.pinfo.setup.get_convergence()
conv.to_csv(out/'linear_readout_hfss_convergence.csv')
desk=save_aedt._desktop()
desk.get_active_project()._project.Save()
(out/'linear_readout_hfss_messages.log').write_text('\n'.join(desk._desktop.GetMessages(
    'D2_Layout3','D2_Layout3_Readout_Linear_hfss',0)),encoding='utf8')
