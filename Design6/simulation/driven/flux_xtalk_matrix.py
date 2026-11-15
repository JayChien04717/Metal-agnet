"""Flux crosstalk matrix C (Phi0 / V) of Design6 "measured" with the transmon model, as in the experiment.

Physics input: the DC transfer H0[n, m] = Phi_n / V_m of the HFSS local flux model (flux_pulse.transfer; V = voltage
of the wave incident at the chip end of line m, 50 Ohm line shorted on chip -> I = 2 V / 50 Ohm).
Procedure (simulated spectroscopy, cz_scqubits.Model = HFSS-fitted Hamiltonian, retargeted junctions):
  * diagonal C[i, i]: sweep line i over a full period, transmon arc f_i(V) -> fit f_i(Phi = C_ii V + b);
  * off-diagonal C[i, j]: park qubit i at Phi_i = 0.25 Phi0 with its own line (steep slope, no sweet-spot sign
    ambiguity), sweep line j by +-V, fit f_i(V_j) = F_i(0.25 + C_ij V_j)  (F_i = single-element transmon map).
The fitted C contains the dressing by the other elements, as a measured matrix does.
Also given: C scaled to the room-temperature AWG through an assumed attenuation (--atten-db, flux lines typically
~20 dB) and the normalised matrix C_ij / C_jj (%).

Usage (scqenv): python driven/flux_xtalk_matrix.py [--run results/2026-11-14/flux_loc_fl10] [--tag A] [--atten-db 20]
Output: <run>/xtalk_matrix.json, <run>/xtalk_matrix.png
"""
import argparse, json, os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__)); SIM = os.path.dirname(HERE)
sys.path[:0] = [os.path.join(SIM, 'cz'), HERE]
import flux_pulse as FP
NAMES = FP.NAMES
KEY = dict(Q1='f1', TC='fc', Q2='f2')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default=os.path.join(SIM, 'results', '2026-11-14', 'flux_loc_fl10'))
    ap.add_argument('--fit', default=os.path.join(SIM, 'results', '2026-11-14', 'hfss_core_v8', 'fit.json'))
    ap.add_argument('--tag', default='A')
    ap.add_argument('--atten-db', type=float, default=20.0)
    ap.add_argument('--park', type=float, default=0.25)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(SIM, a.run)
    import cz_scqubits as C
    f, H = FP.transfer(run)
    Hs = H[0].real / 2.0          # transfer() is per Thevenin V_s (= 2 V_incident): convert to per incident volt
    m = C.Model(C.load_params(a.fit, dict(Q1=4.80, Q2=4.40, TC=6.30)))

    def freq(i, V):                                # frequency of element i for line voltages V (3,)
        phi = Hs @ V
        st = m.static(f1=phi[0], fc=phi[1], f2=phi[2])
        return st[KEY[NAMES[i]]]

    grid = np.linspace(-0.6, 0.6, 241)
    F = {n: np.interp for n in NAMES}
    fmap = {n: np.array([m.static(**{KEY[n]: x})[KEY[n]] for x in grid]) for n in NAMES}
    Fi = lambda n, x: np.interp(x, grid, fmap[n])

    Cm = np.zeros((3, 3)); fits = {}
    for i, ni in enumerate(NAMES):
        # diagonal: full arc
        Vp = 0.5 / Hs[i, i]
        V = np.linspace(-Vp, Vp, 41)
        fi = np.array([freq(i, np.eye(3)[i] * v) for v in V])
        sol = least_squares(lambda x: Fi(ni, x[0] * V + x[1]) - fi, [Hs[i, i] * 0.9, 0.0])
        Cm[i, i] = sol.x[0]
        fits[f'{ni}<-{ni}'] = dict(rms_MHz=float(1e3 * np.sqrt(np.mean(sol.fun ** 2))))
        V_park = a.park / Hs[i, i]
        for j, nj in enumerate(NAMES):
            if j == i:
                continue
            Vj = np.linspace(-0.3, 0.3, 21) / Hs[j, j]            # line j moves its own loop by +-0.3 Phi0
            fi = np.array([freq(i, np.eye(3)[i] * V_park + np.eye(3)[j] * v) for v in Vj])
            sol = least_squares(lambda x: Fi(ni, a.park + x[0] * Vj) - fi, [0.0])
            Cm[i, j] = sol.x[0]
            fits[f'{ni}<-{nj}'] = dict(rms_MHz=float(1e3 * np.sqrt(np.mean(sol.fun ** 2))),
                                       df_MHz_over_sweep=float(1e3 * (fi.max() - fi.min())))
    att = 10 ** (-a.atten_db / 20)
    res = dict(tag=a.tag, units='Phi0 per volt of the wave incident at the chip end of the line', rows='loop', cols='line',
               names=NAMES, C_chip=Cm.tolist(), C_hfss_direct=Hs.tolist(),
               C_room_temp=(Cm * att).tolist(), atten_db=a.atten_db,
               C_normalised_pct=(100 * Cm / np.diag(Cm)[None, :]).tolist(), fits=fits, park_phi0=a.park)
    json.dump(res, open(os.path.join(run, 'xtalk_matrix.json'), 'w'), indent=1)

    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for k, (mat, title, fmt) in enumerate(((Cm * att, f'C (Phi0/V) at the AWG, {a.atten_db:.0f} dB attenuation', '{:.3f}'),
                                           (100 * Cm / np.diag(Cm)[None, :], 'C_ij / C_jj (%)', '{:.1f}'))):
        lim = np.max(np.abs(mat))
        im = ax[k].imshow(mat, cmap='RdBu_r', vmin=-lim, vmax=lim)
        for r in range(3):
            for c in range(3):
                ax[k].text(c, r, fmt.format(mat[r, c]), ha='center', va='center',
                           color='white' if abs(mat[r, c]) > 0.6 * lim else 'black', fontsize=12)
        ax[k].set_xticks(range(3)); ax[k].set_xticklabels([f'line {n}' for n in NAMES])
        ax[k].set_yticks(range(3)); ax[k].set_yticklabels([f'loop {n}' for n in NAMES])
        ax[k].set_title(title); fig.colorbar(im, ax=ax[k], fraction=0.046)
    fig.suptitle(f'Design6 flux crosstalk matrix ({a.tag}), from transmon-model spectroscopy on the HFSS flux network')
    fig.tight_layout(); fig.savefig(os.path.join(run, 'xtalk_matrix.png'), dpi=130); plt.close(fig)
    np.set_printoptions(precision=4, suppress=True)
    print('C at chip (Phi0/V):\n', Cm, '\nHFSS direct:\n', Hs, '\nnormalised %:\n', 100 * Cm / np.diag(Cm)[None, :])
    print('fit rms (MHz):', {k: round(v['rms_MHz'], 3) for k, v in fits.items()})


if __name__ == '__main__':
    main()
