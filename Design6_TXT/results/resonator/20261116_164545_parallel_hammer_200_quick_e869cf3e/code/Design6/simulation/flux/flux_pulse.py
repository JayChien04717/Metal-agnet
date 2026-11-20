"""Flux-pulse distortion and crosstalk of the Design6 local flux lines, from the HFSS 6-port Z matrix, and their effect
on the CZ gate (simulated "Cryoscope", cf. Rol et al., APL 116, 054001 (2020)).

Network (per frequency, HFSS dense sweep of the local model, flux_resweep.py --dense):
  flux ports F (Q1, Q2, TC lines): driven by a 50 Ohm source through a matched 50 Ohm line (the routed CPW is then a
      pure delay and is left out): Thevenin V_s, 50 Ohm  ->  I_F = (Z_FF + 50)^-1 V_s
  junction ports J (one per SQUID, open, I = 0): V_J = Z_JF I_F = j w Phi  ->  Phi_J = Z_JF (Z_FF + 50)^-1 V_s / (j w)
  H[n, m](f) = Phi_n / V_s,m.  DC: H(f_min) (quasi static, f_min = 1 MHz).
Sign convention: each line's drive polarity is chosen so that its own loop sees +Phi (H_hat = H s_m); each loop's
orientation is then the one in which its own line is positive, and crosstalk has a physical sign.

Pulses: flat top with cosine ramps (t_ramp), sampled at 0.05 ns, filtered through H by FFT (window 1 us).
  1. step / 100 ns pulse on each line alone: own-loop waveform error, settling, crosstalk waveforms;
  2. the CZ pulses (Q1 + TC together, cz_scqubits.optimise at T): distorted own waveforms + crosstalk into all three
     loops -> CZ evolution (cz_scqubits model) with (a) ideal pulses, (b) on-chip H without crosstalk,
     (c) on-chip H with crosstalk, (d) (c) after crosstalk compensation (inverse DC matrix applied to the drive),
     (e) (d) + a typical off-chip tail (2 % overshoot, tau 30 ns; Rol 2020 order of magnitude), not pre-distorted.

Usage (scqenv; needs results/<date>/flux_loc_fl10/{Z_dense,f_dense}.npy and the HFSS core fit):
  python driven/flux_pulse.py [--T 100 test pulse] [--cz-T 80] [--ramp 8]
Output: results/<date>/flux_pulse/ (flux_pulse.json, step_response.png, cz_pulses.png)
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

HERE = _S('Design6/simulation/flux')
SIM = _S('Design6/simulation')
sys.path[:0] = [_I('Design6/simulation/qubit'), HERE]
PHI0 = 2.067833848e-15
NAMES = ['Q1', 'TC', 'Q2']
DT = 0.05e-9


F_MAX = 2.0e9   # Hz. Above ~2.5 GHz H has a sharp resonance that belongs to the fixture, not the flux line: with one
                # junction of each SQUID replaced by the open EMF port, the other junction (~29 nH) and the qubit capacitance
                # (~100 fF) resonate at ~2.6-2.9 GHz (the qubit itself, linearised). Flux pulses live below ~1 GHz, so H is
                # truncated at F_MAX and extended flat above it.


def transfer(run_dir, f_max=F_MAX):
    Z = np.load(_relocate(os.path.join(run_dir, 'Z_dense.npy'))); f = np.load(_relocate(os.path.join(run_dir, 'f_dense.npy')))
    keep = f <= f_max
    Z, f = Z[keep], f[keep]
    pm = json.load(open(_relocate(os.path.join(run_dir, 'M_matrix.json'))))['port_map']
    F = [pm['flux'][n] - 1 for n in NAMES]; J = [pm['jj'][n] - 1 for n in NAMES]
    H = np.zeros((len(f), 3, 3), complex)                       # Phi0 per volt of source
    for k, fk in enumerate(f):
        IF = np.linalg.solve(Z[k][np.ix_(F, F)] + 50 * np.eye(3), np.eye(3))      # A per V
        H[k] = Z[k][np.ix_(J, F)] @ IF / (2j * np.pi * fk) / PHI0
    s = np.sign(np.real(np.diag(H[0])))
    H = H * s[None, :]          # line polarity chosen so each line gives +Phi in its own loop (columns only: the loop
                                # orientation is then the one in which its own line is positive; flipping rows as well
                                # flipped the sign of the crosstalk in the rows with s = -1, fixed 2026-11-14 evening)
    return f, H


def apply(f, H, waves, dt=DT, t_win=1e-6):
    """waves: (3, N) source voltages (V) for lines Q1, TC, Q2 -> (3, N) flux (Phi0) in loops Q1, TC, Q2."""
    n = waves.shape[1]
    N = int(2 ** np.ceil(np.log2(max(n * 2, t_win / dt))))
    fr = np.fft.rfftfreq(N, dt)
    Hi = np.zeros((len(fr), 3, 3), complex)
    for a in range(3):
        for b in range(3):
            re = np.interp(fr, f, H[:, a, b].real, left=H[0, a, b].real, right=H[-1, a, b].real)
            im = np.interp(fr, f, H[:, a, b].imag, left=0.0, right=H[-1, a, b].imag)
            Hi[:, a, b] = re + 1j * im
    Hi[0] = H[0].real                                            # DC: quasi-static, real
    W = np.fft.rfft(np.pad(waves, ((0, 0), (0, N - n))), axis=1)
    out = np.einsum('fab,bf->af', Hi, W)
    return np.fft.irfft(out, n=N, axis=1)[:, :n]


def pulse(T, ramp, amp, pre=20e-9, post=200e-9, dt=DT):
    import cz_scqubits as C
    t = np.arange(0, pre + T + post, dt)
    w = np.where((t >= pre) & (t <= pre + T), C.envelope((t - pre) * 1e9, T * 1e9, ramp * 1e9), 0.0)
    return t, amp * w


def tail(w, A=0.02, tau=30e-9, dt=DT):
    """Typical uncorrected off-chip response with the DC amplitude calibrated (unity DC gain): step response
    s(t) = 1 - A exp(-t/tau), i.e. y = x - A x + A (x convolved with exp(-t/tau)/tau): the flux reaches only (1 - A) at
    the edge and settles to the target with time constant tau (the slow approach of a flux-pulse spectroscopy trace).
    (Before 2026-11-14 evening: y = x + A (x * k), DC gain 1 + A; tail_scan / case e re-run with this version.)"""
    k = np.exp(-np.arange(0, 10 * tau, dt) / tau) * dt / tau
    conv = np.array([np.convolve(x, k)[:len(x)] for x in np.atleast_2d(w)]).reshape(w.shape)
    return w - A * w + A * conv


def evolve_wave(m, idle, f1, fc, f2, dt_ns):
    """cz_scqubits.evolve with arbitrary flux waveforms (Phi0) for Q1, TC (offset idle), Q2."""
    from scipy.linalg import eigh
    import cz_scqubits as C
    w0, v0, lab = m.dressed(fc=idle)
    psi = np.stack([v0[:, lab[s]] for s in C.COMP], axis=1).astype(complex)
    for a, c, b in zip(f1, fc, f2):
        ev, V = eigh(m.H(f1=a, fc=idle + c, f2=b))
        psi = V @ (np.exp(-2j * np.pi * ev * dt_ns)[:, None] * (V.conj().T @ psi))
    T = len(f1) * dt_ns
    cc = (v0.conj().T @ psi) * np.exp(2j * np.pi * w0 * T)[:, None]
    return np.array([[cc[lab[s], k] for k in range(4)] for s in C.COMP])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default=_S('Design6/results/flux/2026-11-14/flux_loc_fl10'))
    ap.add_argument('--fit', default=_S('Design6/results/qubit/2026-11-14/hfss_core_v8/fit.json'))
    ap.add_argument('--T', type=float, default=100.0, help='CZ / test pulse length (ns)')
    ap.add_argument('--ramp', type=float, default=8.0)
    ap.add_argument('--cz-T', type=float, default=80.0, help='CZ gate time (ns); 80 ns = best of cz_hfss_v8')
    ap.add_argument('--date', default=time.strftime('%Y-%m-%d'))
    a = ap.parse_args()
    out = _relocate(os.path.join(SIM, 'results', a.date, 'flux_pulse')); os.makedirs(out, exist_ok=True)
    import cz_scqubits as C
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    f, H = transfer(a.run)
    H0 = H[0].real
    res = dict(H_dc_phi0_per_V=H0.tolist(),
               crosstalk_dc_pct={f'{NAMES[n]}<-{NAMES[k]}': 100 * H0[n, k] / H0[k, k] for n in range(3) for k in range(3) if n != k},
               own_gain_vs_f={NAMES[n]: dict(f_MHz=[float(x) for x in f[::max(1, len(f) // 30)] / 1e6],
                                             mag=[float(x) for x in np.abs(H[::max(1, len(f) // 30), n, n] / H0[n, n])],
                                             phase_deg=[float(x) for x in np.degrees(np.angle(H[::max(1, len(f) // 30), n, n]))])
                              for n in range(3)})

    # ---- 1. a T-ns pulse on each line alone (amplitude for 0.1 Phi0 in the own loop) ----
    T, ramp = a.T * 1e-9, a.ramp * 1e-9
    fig, ax = plt.subplots(3, 2, figsize=(12, 9))
    step = {}
    for k, nk in enumerate(NAMES):
        t, w = pulse(T, ramp, 0.1 / H0[k, k])
        src = np.zeros((3, len(t))); src[k] = w
        phi = apply(f, H, src)
        ideal = w * H0[k, k]
        err = phi[k] - ideal
        flat = (t > 20e-9 + ramp) & (t < 20e-9 + T - ramp)
        after = t > 20e-9 + T + 5e-9
        step[nk] = dict(max_err_pct_of_amp=float(100 * np.max(np.abs(err)) / 0.1),
                        flat_top_err_pct=float(100 * np.max(np.abs(err[flat])) / 0.1),
                        tail_after_pulse_pct=float(100 * np.max(np.abs(phi[k][after])) / 0.1),
                        crosstalk_peak_phi0={NAMES[n]: float(np.max(np.abs(phi[n]))) for n in range(3) if n != k})
        ax[k, 0].plot(t * 1e9, ideal, 'k--', lw=1, label='ideal (DC gain)')
        ax[k, 0].plot(t * 1e9, phi[k], label=f'loop {nk} (HFSS H)')
        for n in range(3):
            if n != k:
                ax[k, 0].plot(t * 1e9, phi[n], label=f'crosstalk -> loop {NAMES[n]}')
        ax[k, 0].set_ylabel('flux ($\\Phi_0$)'); ax[k, 0].legend(fontsize=7); ax[k, 0].set_title(f'{a.T:.0f} ns pulse on line {nk}')
        ax[k, 1].plot(t * 1e9, 1e4 * err / 0.1); ax[k, 1].set_ylabel('own-loop error (1e-4 of amplitude)')
        ax[k, 1].set_title(f'{nk}: max {step[nk]["max_err_pct_of_amp"]:.3f} %, after pulse {step[nk]["tail_after_pulse_pct"]:.4f} %')
    for x in ax[-1]:
        x.set_xlabel('t (ns)')
    fig.tight_layout(); fig.savefig(_relocate(os.path.join(out, 'step_response.png')), dpi=130); plt.close(fig)
    res['pulse_on_one_line'] = step
    print(json.dumps(step, indent=1), flush=True)

    # ---- 2. CZ pulses through the on-chip network ----
    p = C.load_params(a.fit, dict(Q1=4.80, Q2=4.40, TC=6.30))
    m = C.Model(p)
    idle = 0.0
    si = m.static(fc=idle)
    target = si['f2'] - si['a1']
    fq = np.linspace(0, 0.3, 61)
    a_q0 = float(np.interp(target, np.array([m.static(f1=x, fc=idle)['f1'] for x in fq])[::-1], fq[::-1]))
    fl = np.linspace(0, 0.4, 81); fcs = np.array([m.static(fc=x)['fc'] for x in fl])
    acq = max(x for x, y in zip(fl, fcs) if y >= target + C.FC_MARGIN)
    opt = C.optimise(m, idle, a.cz_T, a.ramp, a_q0, acq)
    print('CZ ideal pulses:', {k: (round(v, 6) if isinstance(v, float) else v) for k, v in opt.items()}, flush=True)
    t, w = pulse(a.cz_T * 1e-9, ramp, 1.0, pre=0.0, post=0.0)
    n_g = len(t)
    pad = int(round(60e-9 / DT))                                  # also evolve 60 ns after the pulse (tails)
    t_all = np.arange(n_g + pad) * DT
    env = np.concatenate([w, np.zeros(pad)])
    want = np.vstack([opt['a_q'] * env, opt['a_c'] * env, 0 * env])     # desired flux in loops Q1, TC, Q2
    Hdc_inv = np.linalg.inv(H0)
    cases = {}

    def run(name, phi):
        U = evolve_wave(m, idle, phi[0], phi[1], phi[2], DT * 1e9)
        r = C.analyse_U(U)
        cases[name] = dict(**r, max_flux_dev_phi0=[float(np.max(np.abs(phi[i] - want[i]))) for i in range(3)])
        print(f"{name:42s}: phase {r['cond_phase_deg']:7.2f} deg, leakage {r['leakage']:.2e}, F {r['fidelity']:.6f}", flush=True)
        return phi

    run('a ideal flux', want)
    drive_own = np.diag(1 / np.diag(H0)) @ want                   # drive each line for its own loop only
    H_diag = H.copy()
    for i in range(3):
        for j in range(3):
            if i != j:
                H_diag[:, i, j] = 0
    run('b on-chip H, no crosstalk', apply(f, H_diag, drive_own))
    phi_c = run('c on-chip H with crosstalk', apply(f, H, drive_own))
    drive_comp = Hdc_inv @ want                                   # DC crosstalk matrix compensation
    phi_d = run('d crosstalk-compensated drive', apply(f, H, drive_comp))
    phi_e = run('e (d) + off-chip tail 2 %, 30 ns (no predist.)', apply(f, H, tail(drive_comp)))
    res['cz'] = dict(T_ns=a.cz_T, ramp_ns=a.ramp, ideal_pulses=opt, cases=cases)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    for i, nm in enumerate(NAMES):
        ax[i].plot(t_all * 1e9, want[i], 'k--', label='target')
        ax[i].plot(t_all * 1e9, phi_c[i], label='(c) crosstalk, uncompensated')
        ax[i].plot(t_all * 1e9, phi_d[i], label='(d) compensated')
        ax[i].plot(t_all * 1e9, phi_e[i], label='(e) + off-chip tail')
        ax[i].set_title(f'flux in loop {nm}'); ax[i].set_xlabel('t (ns)'); ax[i].set_ylabel('$\\Phi$ ($\\Phi_0$)')
    ax[0].legend(fontsize=7)
    fig.suptitle(f'CZ flux pulses through the HFSS local flux network (T = {a.cz_T:.0f} ns)')
    fig.tight_layout(); fig.savefig(_relocate(os.path.join(out, 'cz_pulses.png')), dpi=130); plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for n in range(3):
        ax[0].semilogx(f / 1e6, np.abs(H[:, n, n] / H0[n, n]), label=NAMES[n])
        ax[1].semilogx(f / 1e6, np.degrees(np.angle(H[:, n, n])), label=NAMES[n])
    ax[0].set_ylabel('|H(f)| / H(DC)'); ax[1].set_ylabel('phase (deg)')
    for x in ax:
        x.set_xlabel('f (MHz)'); x.legend()
    fig.suptitle('on-chip flux transfer (line -> own SQUID), HFSS local model'); fig.tight_layout()
    fig.savefig(_relocate(os.path.join(out, 'transfer.png')), dpi=130); plt.close(fig)
    json.dump(res, open(_relocate(os.path.join(out, 'flux_pulse.json')), 'w'), indent=1, default=float)
    print('saved', out)


if __name__ == '__main__':
    main()
