"""Q3D of the three Design6 capacitors (cin, cp, cj): coupling C between the two electrodes.
Usage: python q3d_caps6.py [cin cp cj]"""
import os, sys, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
from qiskit_metal.analyses.quantization import LOManalysis
import design6_layout as d6
TARGET = {'cin': d6.C_IN_FF, 'cp': d6.C_P_FF, 'cj': d6.C_J_FF}
for name in (sys.argv[1:] or ['cin', 'cp', 'cj']):
    d = d6.build_design()
    if os.environ.get('CJ_OPTS') and name == 'cj':
        import json; d.components['cj'].options.update(json.loads(os.environ['CJ_OPTS'])); d.components['cj'].rebuild()
    if os.environ.get('CJ_LEAD') and name == 'cj':
        d.components['cj'].options.lead_length = os.environ['CJ_LEAD'] + 'um'; d.components['cj'].rebuild()
    lom = LOManalysis(d, 'q3d')
    for k, v in dict(max_passes=12, min_passes=2, percent_error=0.2, freq_ghz=6.7, solution_order='Medium',
                     auto_increase_solution_order=False).items():
        lom.sim.setup[k] = v
    lom.sim.run(name=f'q3d6r_{name}', components=[name], open_terminations=[(name, 'a'), (name, 'b')])
    cm = lom.sim.capacitance_matrix
    cm.to_csv(os.path.join(HERE, 'results', f'q3d6_{name}_cmat_fF.csv'))
    import time; cm.to_csv(os.path.join(HERE, 'results', f'q3d6_{name}_cmat_fF_{time.strftime("%m%d_%H%M")}.csv'))
    sig = [n for n in cm.index if 'ground' not in n.lower()]
    print(f'CAP6 {name}: C = {-cm.loc[sig[0], sig[1]]:.2f} fF (target {TARGET[name]}); '
          f'to ground {[round(cm.loc[s, s] + cm.loc[s, sig[1 - i]], 1) for i, s in enumerate(sig)]} fF', flush=True)
    try: lom.sim.close()
    except Exception: pass
