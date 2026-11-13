"""Design2: HFSS eigenmode frequency of each readout resonator (same method as Design1 / 3).

For every resonator one HFSS eigenmode design renders the resonator route, its OpenToGround
end, the feedline TL (open at both ends) and its qubit (simulation loop: lumped inductors), so
the readout pad and the parallel coupling section load the resonator as in the real device.
The open 8 mm feedline has its own lambda/2 mode near 7.5 GHz: three modes are solved and the
one nearest the target is taken.

Usage:  python hfss_eigen_res.py [Q1 Q2] [--res '{"Q1": {"length": 9.3}}']
"""
import json
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # simulation/
sys.path[:0] = [os.path.join(HERE, 'q3d'), HERE, os.path.join(os.path.dirname(HERE), 'pattern design')]
from qiskit_metal.analyses.simulation.eigenmode import EigenmodeSim  # noqa: E402
from ftf_layout import build_design, RESONATORS  # noqa: E402
from q3d_ftf import RESULTS  # noqa: E402

# as Design1 / 3 (17 passes, 0.05 %), three modes (resonator + feedline mode)
SETUP = dict(min_freq_ghz=4.0, n_modes=3, max_delta_f=0.05, max_passes=17, min_passes=3,
             min_converged=2, pct_refinement=30, basis_order=1)


def run(names, resonators=None, tag='eigen'):
    rows = []
    for name in names:
        d = build_design(for_process=True, resonators=resonators)
        sim = EigenmodeSim(d, 'hfss')
        for k, v in SETUP.items():
            sim.setup[k] = v
        comps = [f'res_{name}', f'otg_{name}', 'TL', name]
        opens = [('TL', 'start'), ('TL', 'end')]
        t0 = time.time()
        sim.run(components=comps, open_terminations=opens)
        conv = sim.convergence_f
        conv.to_csv(os.path.join(RESULTS, f'{tag}_{name}_convergence.csv'))
        last = conv.iloc[-1]
        freqs = sorted(float(v) for v in last.values if pd.notna(v))
        r = dict(RESONATORS[name], **(resonators or {}).get(name, {}))
        f_res = min(freqs, key=lambda f: abs(f - r['f']))
        rows.append(dict(resonator=name, f_target_GHz=r['f'], f_hfss_GHz=f_res,
                         error_MHz=(f_res - r['f']) * 1e3, all_modes_GHz=str([round(f, 4) for f in freqs]),
                         passes=len(conv), minutes=round((time.time() - t0) / 60, 1),
                         length_mm=r['length'], start_straight_mm=r['start_straight']))
        print(f'{name}: modes {freqs} GHz -> resonator {f_res:.4f} GHz (target {r["f"]}) '
              f'after {len(conv)} passes, {(time.time() - t0) / 60:.1f} min', flush=True)
        try:
            sim.close()
        except Exception:
            pass
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, f'{tag}_resonators.csv'), index=False)
    print('\n', df.to_string(index=False))
    return df


if __name__ == '__main__':
    args = sys.argv[1:]
    res = None
    if '--res' in args:
        i = args.index('--res')
        res = json.loads(args[i + 1])
        args = args[:i] + args[i + 2:]
    tag = 'eigen'
    if '--tag' in args:
        i = args.index('--tag')
        tag = args[i + 1]
        args = args[:i] + args[i + 2:]
    run(args or list(RESONATORS), res, tag)
