"""Q3D capacitance simulation of the Design2 FTF set (Q1 - Q3 - Q2, FXTFX24 layout, CPW 10 / 6 um).

Same procedure as Design1/simulation/q3d_ftf.py and as cell 33 of the lab notebook
FXTFX24-W7-0x-V.3.1-ALL_Sims.ipynb: render Q1, Q2, Q3, open-terminate every pin (the flux
lines of this lab design are not shorted to ground in the model), run Q3D through
LOManalysis, save the Maxwell matrix and print Ec / J / readout beta, g / charge line
(method: SQL HFSS/note/capacitance_method.md, SQL HFSS/skill/cap_matrix_analysis.md).

Usage:
    python q3d_ftf.py [tag]
or:
    from q3d_ftf import Q3DRunner
    r = Q3DRunner(); r.run('test', tmon_options=dict(pad_arm_length='400um')); r.close()
"""
import json
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # simulation/ (this script: simulation/<kind>/)
import sys as _sys; _sys.path[:0] = [os.path.join(HERE, _k) for _k in ('driven', 'eigenmode', 'q3d')] + [HERE]
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'pattern design'))
# the capacitance-matrix maths is shared with Design1 (SQL HFSS/note/cap_matrix.py)
for cand in (os.path.join(ROOT, 'note'), os.path.join(ROOT, '..', 'note')):
    if os.path.exists(os.path.join(cand, 'cap_matrix.py')):
        sys.path.insert(0, cand)
        break
RESULTS = os.path.join(HERE, 'results', 'runs')
os.makedirs(RESULTS, exist_ok=True)

from qiskit_metal.analyses.quantization import LOManalysis  # noqa: E402
from ftf_layout import build_design, CPW_Z0  # noqa: E402
import cap_matrix  # noqa: E402

COMPONENTS = ['Q1', 'Q2', 'Q3']
OPEN_PINS = [
    ('Q1', 'charge_line'), ('Q1', 'flux_bias_line'), ('Q1', 'fake_flux_bias_line'), ('Q1', 'readout_line'),
    ('Q2', 'charge_line'), ('Q2', 'flux_bias_line'), ('Q2', 'fake_flux_bias_line'), ('Q2', 'readout_line'),
    ('Q3', 'flux_bias_line'), ('Q3', 'fake_flux_bias_line'),
]
# 0-based: 0,1 Q1 pads; 2 coupler; 3,4 Q2 pads; 5,6 readout pads
NODES5 = ['pad_bot_Q1', 'pad_top_Q1', 'tmon_0_Q3', 'pad_top_Q2', 'pad_bot_Q2']
NODES7 = NODES5 + ['readout_padNline_Q1', 'readout_padNline_Q2']
# resonator frequencies of the lab design (notebook: cavity1 / cavity2)
F_RES = {'Q1': 6.109e9, 'Q2': 6.475e9}

SETUP = dict(max_passes=6, min_passes=2, percent_error=0.5, freq_ghz=5.0,
             solution_order='Medium', auto_increase_solution_order=False)


class Q3DRunner:
    def __init__(self, setup=None):
        self.setup = dict(SETUP, **(setup or {}))
        self.lom = None

    def run(self, tag, fx1_options=None, fx2_options=None, tmon_options=None):
        design = build_design(fx1_options=fx1_options, fx2_options=fx2_options, tmon_options=tmon_options)
        self.lom = LOManalysis(design, 'q3d')
        for k, v in self.setup.items():
            self.lom.sim.setup[k] = v
        t0 = time.time()
        self.lom.sim.run(components=COMPONENTS, open_terminations=OPEN_PINS)
        cmat = self.lom.sim.capacitance_matrix
        dt = time.time() - t0
        cmat.to_csv(os.path.join(RESULTS, f'{tag}_cmat_fF.csv'))
        print(f'\n=== {tag}  ({dt / 60:.1f} min) ===')
        print(cmat.round(3).to_string())

        res = {}
        # Ec / J: readout pads, flux and charge lines grounded (5 nodes), and with the readout pads (7 nodes)
        print('--- Ec / J, 5 nodes ---')
        res['5nodes'] = cap_matrix.analyse(cmat, nodes=NODES5)
        res['C_sigma [fF]'] = cap_matrix.total_capacitances(cmat, NODES5)
        print('C_sigma [fF]', {k: round(v, 2) for k, v in res['C_sigma [fF]'].items()})
        print('--- Ec / J, 7 nodes (with readout pads) ---')
        res['7nodes'] = cap_matrix.analyse(cmat, nodes=NODES7)
        print('--- readout beta / g (pad only, no resonator rendered; Z0 = %.1f Ohm) ---' % CPW_Z0)
        res['readout_Q1'] = cap_matrix.readout_g(cmat, ('pad_bot_Q1', 'pad_top_Q1'), 'readout_padNline_Q1',
                                                 Z0=CPW_Z0, f_r=F_RES['Q1'])
        res['readout_Q2'] = cap_matrix.readout_g(cmat, ('pad_top_Q2', 'pad_bot_Q2'), 'readout_padNline_Q2',
                                                 Z0=CPW_Z0, f_r=F_RES['Q2'])
        print('--- charge lines (f01 0.5 GHz, |n01| 0.1 assumed) ---')
        for q, pads in (('Q1', ('pad_bot_Q1', 'pad_top_Q1')), ('Q2', ('pad_top_Q2', 'pad_bot_Q2'))):
            node = cap_matrix.find_node(cmat, f'charge_line_{q}')
            res[f'charge_{q}'] = cap_matrix.charge_line(cmat, pads, node, 0.5e9, 0.1)
            print(f'  {q}:', ', '.join(f'{a} = {b:.4g}' for a, b in res[f'charge_{q}'].items()))
        summary = dict(tag=tag, fx1_options=fx1_options or {}, fx2_options=fx2_options or {},
                       tmon_options=tmon_options or {}, setup=self.setup, minutes=round(dt / 60, 2),
                       results=res)
        with open(os.path.join(RESULTS, f'{tag}_summary.json'), 'w') as f:
            json.dump(summary, f, indent=2, default=float)
        return res

    def close(self):
        if self.lom is not None:
            try:
                self.lom.sim.close()
            except Exception:
                pass


if __name__ == '__main__':
    tag = sys.argv[1] if len(sys.argv) > 1 else 'cpw10_6'
    r = Q3DRunner()
    r.run(tag)
    r.close()
