"""Simulated flux-pulse spectroscopy (qubit frequency vs excitation-pulse delay during a square flux pulse), ns scale.

Square pulse on one flux line (raised-cosine edges `edge` ns, flat `T` ns, amplitude = own-loop flux `amp` Phi0),
filtered by the on-chip transfer H(f) of the HFSS local model (flux_pulse.transfer: 50 Ohm source, matched line,
6-port Z), optionally followed by an assumed off-chip tail (y = x + A x * exp(-t/tau)/tau, not pre-distorted).
Loop flux -> qubit frequencies through the HFSS-fitted Hamiltonian (cz_scqubits.Model, retargeted junctions):
f_Q1(t), plus the neighbours (crosstalk). The "measured" frequency at delay t is the instantaneous transition frequency
(a short excitation pulse samples it, as in the measurement).

Usage (scqenv): python driven/flux_square_spectro.py [--line Q1] [--amp 0.24] [--T 50] [--edge 1]
Output: results/<date>/flux_pulse/square_spectro_<line>.png / .json
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np

HERE = _S('Design6/simulation/flux'); SIM = _S('Design6/simulation')
sys.path[:0] = [_I('Design6/simulation/qubit'), HERE]
import flux_pulse as FP
NAMES = FP.NAMES


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default=_S('Design6/results/flux/2026-11-14/flux_loc_fl10'))
    ap.add_argument('--fit', default=_S('Design6/results/qubit/2026-11-14/hfss_core_v8/fit.json'))
    ap.add_argument('--line', default='Q1', choices=NAMES)
    ap.add_argument('--amp', type=float, default=0.24, help='own-loop flux of the pulse (Phi0)')
    ap.add_argument('--T', type=float, default=50.0, help='flat-top length (ns)')
    ap.add_argument('--edge', type=float, default=1.0, help='rise / fall time (ns)')
    ap.add_argument('--tail', default='0.02,10', help='assumed off-chip tail A,tau_ns for the comparison curve')
    ap.add_argument('--tag', default='', help='suffix of the output files (flux-line variant)')
    ap.add_argument('--date', default=time.strftime('%Y-%m-%d'))
    a = ap.parse_args()
    out = _relocate(os.path.join(SIM, 'results', a.date, 'flux_pulse')); os.makedirs(out, exist_ok=True)
    import cz_scqubits as C
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    f, H = FP.transfer(a.run)
    H0 = H[0].real
    k = NAMES.index(a.line)

    dt = 0.02e-9
    t = np.arange(-20e-9, a.T * 1e-9 + 120e-9, dt)
    e, T = a.edge * 1e-9, a.T * 1e-9
    up = 0.5 * (1 - np.cos(np.pi * np.clip(t / e, 0, 1)))
    down = 0.5 * (1 - np.cos(np.pi * np.clip((t - T - e) / e, 0, 1)))
    ideal = up - down                                            # square pulse, own-loop units
    src = np.zeros((3, len(t))); src[k] = a.amp * ideal / H0[k, k]
    phi = FP.apply(f, H, src, dt=dt, t_win=4e-6)                  # on-chip only
    A, tau = (float(x) for x in a.tail.split(','))
    phi_tail = FP.apply(f, H, FP.tail(src, A, tau * 1e-9, dt=dt), dt=dt, t_win=4e-6)

    m = C.Model(C.load_params(a.fit, dict(Q1=4.80, Q2=4.40, TC=6.30)))
    grid = np.linspace(-0.05, 0.45, 201)
    fmap = {n: np.array([m.static(**{dict(Q1='f1', TC='fc', Q2='f2')[n]: x})[dict(Q1='f1', TC='fc', Q2='f2')[n]]
                         for x in grid]) for n in NAMES}
    tofreq = lambda n, x: np.interp(x, grid, fmap[n])
    fr = {n: tofreq(n, phi[i]) for i, n in enumerate(NAMES)}
    fr_tail = {n: tofreq(n, phi_tail[i]) for i, n in enumerate(NAMES)}
    f_ideal = tofreq(a.line, a.amp * ideal)
    f_on = tofreq(a.line, a.amp)

    def settle(fq, t_edge, sign):
        """|f - f_target| (MHz) at 2, 5, 10, 20, 30 ns after an edge."""
        tgt = f_on if sign > 0 else tofreq(a.line, 0.0)
        return {f'{x} ns': float(1e3 * abs(np.interp(t_edge + x * 1e-9, t, fq) - tgt)) for x in (2, 5, 10, 20, 30)}

    t_rise_end, t_fall_end = e, T + 2 * e
    res = dict(line=a.line, amp_phi0=a.amp, T_ns=a.T, edge_ns=a.edge, f_idle_GHz=float(tofreq(a.line, 0)),
               f_on_GHz=float(f_on),
               on_chip=dict(after_rise_MHz=settle(fr[a.line], t_rise_end, 1), after_fall_MHz=settle(fr[a.line], t_fall_end, -1)),
               with_tail=dict(A=A, tau_ns=tau, after_rise_MHz=settle(fr_tail[a.line], t_rise_end, 1),
                              after_fall_MHz=settle(fr_tail[a.line], t_fall_end, -1)),
               neighbours_max_shift_MHz={n: float(1e3 * np.max(np.abs(fr[n] - tofreq(n, 0.0)))) for n in NAMES if n != a.line})

    fig = plt.figure(figsize=(14, 8))
    g = fig.add_gridspec(2, 3)
    ax0 = fig.add_subplot(g[0, :])
    ax0.plot(t * 1e9, f_ideal, 'k--', lw=1, label='ideal square')
    ax0.plot(t * 1e9, fr[a.line], 'r', lw=1.5, label='on-chip (HFSS H)')
    ax0.plot(t * 1e9, fr_tail[a.line], 'b', lw=1.5, label=f'+ off-chip tail (assumed {100*A:.0f} %, {tau:.0f} ns)')
    ax0.set_xlabel('excitation-pulse delay (ns)'); ax0.set_ylabel(f'{a.line} frequency (GHz)'); ax0.legend()
    ax0.set_title(f'{a.line}: square flux pulse {a.T:.0f} ns, {a.edge:.0f} ns edges, {a.amp} $\\Phi_0$')
    for j, (t1, t2, ttl) in enumerate(((-2, 32, 'rising edge (first 30 ns)'), (a.T + a.edge - 2, a.T + a.edge + 32, 'falling edge (next 30 ns)'))):
        ax = fig.add_subplot(g[1, j])
        ax.plot(t * 1e9, 1e3 * (fr[a.line] - f_ideal), 'r', label='on-chip')
        ax.plot(t * 1e9, 1e3 * (fr_tail[a.line] - f_ideal), 'b', label='+ off-chip tail')
        ax.set_xlim(t1, t2); ax.set_xlabel('delay (ns)'); ax.set_ylabel('f - f_ideal (MHz)'); ax.set_title(ttl)
        lim = 1e3 * np.max(np.abs((fr_tail[a.line] - f_ideal)[(t * 1e9 > t1 + 4) & (t * 1e9 < t2)]))
        ax.set_ylim(-1.3 * lim, 1.3 * lim); ax.axhline(0, c='k', lw=0.5); ax.legend(fontsize=8)
    ax = fig.add_subplot(g[1, 2])
    for n in NAMES:
        if n != a.line:
            ax.plot(t * 1e9, 1e3 * (fr[n] - tofreq(n, 0.0)), label=f'{n} (crosstalk)')
    ax.set_xlabel('delay (ns)'); ax.set_ylabel('frequency shift (MHz)'); ax.set_title('neighbours during the pulse'); ax.legend()
    fig.tight_layout()
    png = _relocate(os.path.join(out, f'square_spectro_{a.line}{"_" + a.tag if a.tag else ""}.png')); fig.savefig(png, dpi=120); plt.close(fig)
    json.dump(res, open(_relocate(os.path.join(out, f'square_spectro_{a.line}{"_" + a.tag if a.tag else ""}.json')), 'w'), indent=1)
    print(json.dumps(res, indent=1)); print('saved', png)


if __name__ == '__main__':
    main()
