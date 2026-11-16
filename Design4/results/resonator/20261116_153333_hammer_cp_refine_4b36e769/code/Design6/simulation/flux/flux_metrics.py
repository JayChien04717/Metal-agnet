"""Flux-line acceptance metrics of one Design6 flux variant (HFSS local model, dense sweep):

  1. own-loop step response, first 30 ns (ring-up): a step (2 ns raised-cosine edge) on each line through the
     on-chip transfer H(f) (flux_pulse.transfer: 50 Ohm source, matched line, 6-port Z) -> 10-90 % rise, overshoot,
     |error| at 2 / 5 / 10 / 20 / 30 ns after the edge, and the impulse response h(t) (derivative of the step).
  2. DC crosstalk matrix C[n, m] = Phi_n / Phi_m(own) for a DC change on line m, and what it does to the qubits:
     frequency shift of loop n when line m moves its own loop by 0.1 Phi0, at the idle sweet spot (2nd order) and
     at the CZ operating point of Q1 (Phi_Q1 = 0.134 Phi0, 1st order), from the HFSS-fitted Hamiltonian.
  3. own mutual inductance M (pH) and I_Phi0 (mA); requirement M > 1 pH.

Usage (scqenv): python driven/flux_metrics.py --run results/2026-11-14/flux_loc_fl10 --tag A
Output: <run>/metrics.json, <run>/step_30ns.png
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse, json, os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np

HERE = _S('Design6/simulation/flux'); SIM = _S('Design6/simulation')
sys.path[:0] = [_I('Design6/simulation/qubit'), HERE]
import flux_pulse as FP
NAMES = FP.NAMES
PHI0 = FP.PHI0


def step_metrics(t, y, t0, edge):
    """y: own flux normalised to the DC gain (ideal step = 1 after the edge)."""
    tt = t - t0
    i10 = np.argmax(y >= 0.1); i90 = np.argmax(y >= 0.9)
    after = tt >= edge
    out = dict(rise_10_90_ns=float((t[i90] - t[i10]) * 1e9), overshoot_pct=float(100 * max(0.0, y[after].max() - 1)))
    for k in (2, 5, 10, 20, 30):
        sel = (tt >= edge + k * 1e-9 - 0.5e-9) & (tt <= edge + k * 1e-9 + 0.5e-9)
        out[f'err_at_{k}ns_pct'] = float(100 * np.max(np.abs(y[sel] - 1)))
    win = (tt >= edge) & (tt <= edge + 30e-9)
    out['max_err_edge_to_30ns_pct'] = float(100 * np.max(np.abs(y[win] - 1)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', required=True)
    ap.add_argument('--tag', default='')
    ap.add_argument('--fit', default=_S('Design6/results/qubit/2026-11-14/hfss_core_v8/fit.json'))
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else _relocate(os.path.join(SIM, a.run))
    f, H = FP.transfer(run)
    H0 = H[0].real
    Z = np.load(_relocate(os.path.join(run, 'Z_dense.npy')))
    pm = json.load(open(_relocate(os.path.join(run, 'M_matrix.json'))))['port_map']
    w0 = 2 * np.pi * f[0]
    M = {n: abs(Z[0, pm['jj'][n] - 1, pm['flux'][n] - 1].imag / w0) * 1e12 for n in NAMES}
    res = dict(variant=a.tag, run=run,
               M_self_pH=M, I_phi0_mA={n: PHI0 / (M[n] * 1e-12) * 1e3 for n in NAMES},
               M_above_1pH={n: bool(M[n] > 1.0) for n in NAMES},
               crosstalk_dc_pct={f'loop {NAMES[n]} <- line {NAMES[m]}': 100 * H0[n, m] / H0[m, m]
                                 for n in range(3) for m in range(3) if n != m})

    # ---- 1. step response, first 30 ns ----
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    dt = 0.02e-9
    t = np.arange(0, 200e-9, dt); t0, edge = 10e-9, 2e-9
    s = np.clip((t - t0) / edge, 0, 1); s = 0.5 * (1 - np.cos(np.pi * s))   # 2 ns raised-cosine step (< ~0.5 GHz)
    s = np.where(t > 150e-9, 0.5 * (1 + np.cos(np.pi * np.clip((t - 150e-9) / edge, 0, 1))), s)  # back to 0 (no FFT wrap)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    res['step_first_30ns'] = {}
    for k, nk in enumerate(NAMES):
        src = np.zeros((3, len(t))); src[k] = s
        phi = FP.apply(f, H, src, dt=dt, t_win=4e-6)
        y = phi[k] / H0[k, k]
        res['step_first_30ns'][nk] = step_metrics(t, y, t0, edge)
        ax[0].plot((t - t0) * 1e9, 100 * (y - s), label=f'{nk}')
        h = np.gradient(y, dt) * 1e-9
        ax[1].plot((t - t0) * 1e9, h, label=nk)
    ax[0].set_xlim(-2, 30); ax[0].set_xlabel('t after step (ns)'); ax[0].set_ylabel('own flux - ideal step (% of step)')
    ax[0].set_title('step response error, first 30 ns'); ax[0].legend()
    ax[1].set_xlim(-1, 5); ax[1].set_xlabel('t (ns)'); ax[1].set_ylabel('impulse response (1/ns)'); ax[1].legend()
    ax[1].set_title('impulse response (2 ns edge)')
    fig.suptitle(f'Design6 flux line {a.tag}: on-chip (HFSS local model)'); fig.tight_layout()
    fig.savefig(_relocate(os.path.join(run, 'step_30ns.png')), dpi=130); plt.close(fig)

    # ---- 2. what a DC change on a neighbour does to the qubit frequencies ----
    import cz_scqubits as C
    m = C.Model(C.load_params(a.fit, dict(Q1=4.80, Q2=4.40, TC=6.30)))
    def freqs(f1=0.0, fc=0.0, f2=0.0):
        st = m.static(f1=f1, fc=fc, f2=f2)
        return dict(Q1=st['f1'], TC=st['fc'], Q2=st['f2'])
    base_idle = freqs()
    base_gate = freqs(f1=0.134)
    shifts = {}
    for k, nk in enumerate(NAMES):                       # line k moves its own loop by +0.1 Phi0
        dphi = 0.1 * H0[:, k] / H0[k, k]                 # flux in loops Q1, TC, Q2
        idle = freqs(f1=dphi[0], fc=dphi[1], f2=dphi[2])
        gate = freqs(f1=0.134 + dphi[0], fc=dphi[1], f2=dphi[2])
        shifts[f'line {nk} +0.1 Phi0'] = {
            'idle (all at sweet spot)': {n: 1e3 * (idle[n] - base_idle[n]) for n in NAMES if n != nk},
            'Q1 parked at 0.134 Phi0': {n: 1e3 * (gate[n] - base_gate[n]) for n in NAMES if n != nk}}
    res['neighbour_dc_shift_MHz'] = shifts
    json.dump(res, open(_relocate(os.path.join(run, 'metrics.json')), 'w'), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
