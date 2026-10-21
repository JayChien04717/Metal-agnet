"""Q3D capacitance simulation of one FTF set, rendered with qiskit-metal.

Renders FX1, TC and FX2 of build_design(for_process=True), so each fluxonium's JJ
and array are one lumped inductor, with every qubit pin open-terminated. It runs
Q3D through LOManalysis, saves the Maxwell matrix and prints Ec / J / g
(method: ../note/capacitance_method.md).

Usage:
    python q3d_ftf.py                       # current design
or from a notebook / another script:
    from q3d_ftf import Q3DRunner
    r = Q3DRunner(); res = r.run('test', fx_options=dict(pad_width='200um')); r.close()
"""
import json
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'pattern design'))
sys.path.insert(0, os.path.join(ROOT, 'note'))
# new runs are written to results/runs (latest/ and previous/ hold the runs behind the notes)
RESULTS = os.path.join(HERE, 'results', 'runs')
os.makedirs(RESULTS, exist_ok=True)

from qiskit_metal.analyses.quantization import LOManalysis  # noqa: E402
from ftf_layout import build_design, RESONATORS, CPW_Z0  # noqa: E402
import cap_matrix  # noqa: E402

ANSYS_PATH = r'C:\Program Files\AnsysEM\AnsysEM21.2\Win64'
COMPONENTS = ['FX1', 'TC', 'FX2']
# flux lines are shorted to ground (drawn as ground plane in the simulation
# geometry), so only the readout pads and charge lines are open, floating nets
OPEN_PINS = [
    ('FX1', 'readout_line'), ('FX1', 'charge_line'),
    ('FX2', 'readout_line'), ('FX2', 'charge_line'),
    ('TC', 'readout_line'), ('TC', 'charge_line'),
]
# Medium order, no auto increase (as the FXTFX24 reference); 6 passes is enough here (user's choice)
# with_resonators=True: also render the two readout resonators. Their far (tee) end
# is open; the other end merges with the FX readout pad into one net.
RES_COMPONENTS = ['res_FX1', 'res_TC', 'res_FX2']
RES_OPEN_PINS = [('res_FX1', 'start'), ('res_TC', 'start'), ('res_FX2', 'start')]

# resonator design frequencies (used for g)
F_RES = {k: r['f'] * 1e9 for k, r in RESONATORS.items()}

SETUP = dict(max_passes=6, min_passes=2, percent_error=0.5, freq_ghz=5.0,
             solution_order='Medium', auto_increase_solution_order=False)


class Q3DRunner:
    """Runs Q3D for one or more geometries in the same Ansys session."""

    def __init__(self, setup=None):
        self.setup = dict(SETUP, **(setup or {}))
        self.design = None
        self.lom = None

    def _prepare(self, design):
        # a new LOManalysis per geometry; it reconnects to the Ansys session already open
        self.lom = LOManalysis(design, 'q3d')
        for k, v in self.setup.items():
            self.lom.sim.setup[k] = v
        return self.lom

    def run(self, tag, fx_options=None, tmon_options=None, fx_x=None, nodes=None,
            with_resonators=False, Z0_res=CPW_Z0, fx1_options=None, fx2_options=None):
        design = build_design(for_process=True, fx_options=fx_options, tmon_options=tmon_options,
                              fx_x=fx_x, fx1_options=fx1_options, fx2_options=fx2_options)
        lom = self._prepare(design)
        t0 = time.time()
        comps, pins = list(COMPONENTS), list(OPEN_PINS)
        if with_resonators:
            comps += RES_COMPONENTS
            pins = [p for p in pins if p[1] != 'readout_line'] + RES_OPEN_PINS
        lom.sim.run(components=comps, open_terminations=pins)
        cmat = lom.sim.capacitance_matrix
        dt = time.time() - t0

        cmat.to_csv(os.path.join(RESULTS, f'{tag}_cmat_fF.csv'))
        print(f'\n=== {tag}  ({dt / 60:.1f} min) ===')
        print(cmat.round(3).to_string())
        if with_resonators:
            r1 = cap_matrix.find_node(cmat, 'res_FX1', 'readout_line_FX1')
            r2 = cap_matrix.find_node(cmat, 'res_FX2', 'readout_line_FX2')
            rt = cap_matrix.find_node(cmat, 'res_TC', 'readout_line_TC')
            nodes = nodes or cap_matrix.DEFAULT_NODES[:5] + [r1, r2]
        nodes = nodes or cap_matrix.DEFAULT_NODES
        res = cap_matrix.analyse(cmat, nodes=nodes)
        res['C_sigma [fF]'] = cap_matrix.total_capacitances(cmat, nodes)
        if with_resonators:
            print('readout coupling (beta * V_zpf, lambda/2, Z0 = %.0f Ohm, design f_r):' % Z0_res)
            res['readout_FX1'] = cap_matrix.readout_g(cmat, ('pad_bot_FX1', 'pad_top_FX1'), r1,
                                                      Z0=Z0_res, f_r=F_RES['FX1'])
            res['readout_FX2'] = cap_matrix.readout_g(cmat, ('pad_top_FX2', 'pad_bot_FX2'), r2,
                                                      Z0=Z0_res, f_r=F_RES['FX2'])
            res['readout_TC'] = cap_matrix.readout_g(cmat, ('tmon_TC',), rt, Z0=Z0_res,
                                                     f_r=F_RES['TC'], transmon_EJ_over_Ec=50)
            # charge lines: example operating points (FX f01 0.5 GHz, |n01| 0.1; TC 5 GHz, EJ/Ec 50)
            print('charge lines:')
            res['charge_FX1'] = cap_matrix.charge_line(cmat, ('pad_bot_FX1', 'pad_top_FX1'),
                                                       'charge_line_FX1', 0.5e9, 0.1)
            res['charge_TC'] = cap_matrix.charge_line(cmat, ('tmon_TC',), 'charge_line_TC',
                                                      5e9, (50 / 32) ** 0.25)
            for k in ('charge_FX1', 'charge_TC'):
                print(f'  {k}:', ', '.join(f'{a} = {b:.4g}' for a, b in res[k].items()))
        summary = dict(tag=tag, fx_options=fx_options or {}, tmon_options=tmon_options or {}, fx_x=fx_x,
                       setup=self.setup, minutes=round(dt / 60, 2), results=res)
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
    r = Q3DRunner()
    r.run('baseline')
    r.close()
