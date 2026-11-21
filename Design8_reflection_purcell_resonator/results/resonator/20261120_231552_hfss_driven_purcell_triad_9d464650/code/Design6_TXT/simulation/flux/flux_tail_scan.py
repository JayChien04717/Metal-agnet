"""Tolerance of the 80 ns CZ to an uncorrected off-chip flux-pulse tail (y = x + A x * exp(-t/tau)/tau), applied to the
crosstalk-compensated drives through the HFSS local network (flux_pulse.py). Output: results/<date>/flux_pulse/tail_scan.json"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
HERE = _S('Design6_TXT/simulation/flux'); SIM = _S('Design6_TXT/simulation')
sys.path[:0] = [_I('Design6_TXT/simulation/qubit'), HERE]
import cz_scqubits as C
import flux_pulse as FP
R = _S('Design6_TXT/results/general/legacy/2026-11-14')
f, H = FP.transfer(_S('Design6_TXT/results/flux/2026-11-14/flux_loc_fl10'))
H0 = H[0].real
opt = json.load(open(_S('Design6_TXT/results/flux/2026-11-14/flux_pulse/flux_pulse.json')))['cz']['ideal_pulses']
m = C.Model(C.load_params(_S('Design6_TXT/results/qubit/2026-11-14/hfss_core_v8/fit.json'), dict(Q1=4.80, Q2=4.40, TC=6.30)))
t, w = FP.pulse(80e-9, 8e-9, 1.0, pre=0.0, post=0.0)
env = np.concatenate([w, np.zeros(int(round(60e-9 / FP.DT)))])
want = np.vstack([opt['a_q'] * env, opt['a_c'] * env, 0 * env])
drive = np.linalg.inv(H0) @ want
out = []
for tau in (10e-9, 30e-9, 100e-9):
    for A in (0.001, 0.003, 0.01, 0.02):
        phi = FP.apply(f, H, FP.tail(drive, A, tau))
        r = C.analyse_U(FP.evolve_wave(m, 0.0, phi[0], phi[1], phi[2], FP.DT * 1e9))
        out.append(dict(A_pct=100 * A, tau_ns=tau * 1e9, **r))
        print(f"tau {tau*1e9:5.0f} ns, A {100*A:4.1f} %: phase {r['cond_phase_deg']:7.2f}, leakage {r['leakage']:.2e}, F {r['fidelity']:.6f}", flush=True)
json.dump(out, open(_S('Design6_TXT/results/flux/2026-11-14/flux_pulse/tail_scan.json'), 'w'), indent=1)
