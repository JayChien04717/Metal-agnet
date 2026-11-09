"""Eigenmodes of one tee + resonator from 1 GHz up: is the fundamental ~3.4 GHz (resonator shorted)?"""
import os, sys, warnings
warnings.filterwarnings('ignore')
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim
import design3_purcell as dp
name = 'FX1'
d = dp.build_design()
sim = EigenmodeSim(d, 'hfss')
for k, v in dict(min_freq_ghz=5.0, n_modes=2, max_delta_f=0.01, max_passes=17, min_passes=8,
                 min_converged=3, pct_refinement=30, basis_order=1).items():
    sim.setup[k] = v
sim.run(name='eiglow', components=[f'res_{name}', f'tee_{name}'],
        open_terminations=[(f'tee_{name}', 'prime_start'), (f'tee_{name}', 'prime_end'), (f'res_{name}', 'end')])
conv = sim.convergence_f
print('CONV_PER_PASS'); print(conv.round(4).to_string())
last = conv.iloc[-1]
print('LOWMODES', sorted(round(float(v), 4) for v in last.values if pd.notna(v)))
print('resonator length mm', dp.RESONATORS[name]['length'])
try: sim.close()
except Exception: pass
