"""Q3D of one finger coupler of design3_purcell (COUPLER='finger'): C between the two combs vs target Cc.
Usage: python q3d_finger.py [FX1|TC|FX2]"""
import os, sys, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # simulation/ (this script: simulation/<kind>/)
import sys as _sys; _sys.path[:0] = [os.path.join(HERE, _k) for _k in ('driven', 'eigenmode', 'q3d')] + [HERE]
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from qiskit_metal.analyses.quantization import LOManalysis
import design3_purcell as dp
SETUP = dict(max_passes=12, min_passes=2, percent_error=0.2, freq_ghz=6.7,
             solution_order='Medium', auto_increase_solution_order=False)
name = sys.argv[1] if len(sys.argv) > 1 else 'FX1'
d = dp.build_design()
lom = LOManalysis(d, 'q3d')
for k, v in SETUP.items():
    lom.sim.setup[k] = v
lom.sim.run(name=f'q3dfingerB_{name}', components=[f'cc_{name}'],
            open_terminations=[(f'cc_{name}', 'a'), (f'cc_{name}', 'b')])
cm = lom.sim.capacitance_matrix
cm.to_csv(os.path.join(HERE, 'results', f'q3dfinger_{name}_cmat_fF.csv'))
print(cm.round(3).to_string())
sig = [n for n in cm.index if 'ground' not in n.lower()]
r = dp.RESONATORS[name]
print(f'FINGER {name}: Cc = {-cm.loc[sig[0], sig[1]]:.3f} fF  (target {r["Cc"]:.2f} fF, overlap {r["overlap"]} um)')
try: lom.sim.close()
except Exception: pass
