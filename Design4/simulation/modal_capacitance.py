"""Convert Design6 Q3D capacitance to modal coupling at the hammer tap.

The Q3D series mutual is an equipotential-net value. Using it as a capacitor
at the hammer tap is a localized approximation to the electrostatic coupling.
The HFSS inverse conversion uses the same normalization to report an equivalent
capacitance. It cannot separate other coherent coupling paths.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq, root

HERE = Path(__file__).resolve().parent
DESIGN_ROOT = HERE.parent
sys.path.insert(0, str(DESIGN_ROOT / 'pattern design'))
from purcell_6500 import tune_lengths, theory5


def tap_admittances(f_ghz, p):
    """Diagonal admittances with the Cj diagonal loading retained, in siemens."""
    omega = 2 * np.pi * f_ghz * 1e9
    beta = omega / p['v'] * 1e-3
    y0 = 1 / p['Z0']

    def line(length, load):
        tangent = np.tan(beta * length)
        return y0 * (load + 1j * y0 * tangent) / (y0 + 1j * load * tangent)

    def cap(value):
        return 1j * omega * value * 1e-15

    cp, cp_feed, cp_filter = p['cp']
    cj, cj_filter, cj_readout = p['cj']
    feed = line(p['l_stub'], theory5.n0_env_Y(f_ghz, p)) + cap(cp_feed)
    filter_input = cap(cp_filter) + cap(cp) * feed / (cap(cp) + feed)
    filter_tap = (line(p['D_cj'], filter_input)
                  - 1j * y0 / np.tan(beta * (p['L_filt'] - p['D_cj']))
                  + cap(cj + cj_filter))
    resonator_tap = (1j * y0 * np.tan(beta * p['d_r'])
                     - 1j * y0 / np.tan(beta * (p['L_res'] - p['d_r']))
                     + cap(cj + cj_readout))
    return np.array([filter_tap, resonator_tap])


def modal_normalization(p):
    """Find positive-slope susceptance zeros and local equivalent capacitances."""
    grid = np.linspace(5.8, 7.2, 401)
    frequencies, capacitances, linewidths = [], [], []
    for index in range(2):
        def susceptance(f):
            return tap_admittances(f, p)[index].imag

        samples = np.array([susceptance(f) for f in grid])
        crossings = np.flatnonzero((samples[:-1] <= 0) & (samples[1:] >= 0))
        if len(crossings) == 0:
            raise RuntimeError(f'No mode found for tap {index}')
        roots = [brentq(susceptance, grid[k], grid[k+1]) for k in crossings]
        f = min(roots, key=lambda value: abs(value - 6.5))
        step = 1e-5  # GHz
        slope = (susceptance(f+step) - susceptance(f-step)) / (2*step*2*np.pi*1e9)
        c_mode = slope / 2
        if not np.isfinite(c_mode) or c_mode <= 0:
            raise RuntimeError(f'Invalid modal capacitance for tap {index}')
        frequencies.append(f)
        capacitances.append(c_mode * 1e15)
        linewidths.append(tap_admittances(f, p)[index].real / c_mode / (2*np.pi*1e6))
    return np.array(frequencies), np.array(capacitances), linewidths


def mhz_per_ff(frequencies, capacitances):
    return float(500 * np.sqrt(np.prod(frequencies) / np.prod(capacitances)))


def coupled_poles(p, frequencies, widths, j_mhz):
    """Solve det(Ytap) = 0 with the full distributed self-admittances."""
    angular_cap = 2 * np.pi * 1e9 * p['cj'][0] * 1e-15

    def det(f):
        diagonal = tap_admittances(f, p)
        return diagonal.prod() + (f * angular_cap)**2

    diagonal = frequencies + .5j*np.asarray(widths)/1000
    offdiag = j_mhz/1000
    guesses = np.linalg.eigvals([[diagonal[0], offdiag], [offdiag, diagonal[1]]])
    poles = []
    for guess in guesses:
        def residual(x):
            value = 1e6 * det(complex(*x))
            return [value.real, value.imag]

        solution = root(residual, [guess.real, guess.imag], tol=1e-10)
        if np.linalg.norm(residual(solution.x)) > 1e-8:
            raise RuntimeError(f'Coupled-pole solve failed: {solution.message}')
        poles.append(complex(*solution.x))
    poles.sort(key=lambda pole: pole.real)
    if abs(poles[1] - poles[0]) < 1e-6:
        raise RuntimeError('Coupled-pole solve returned the same root twice')
    return [[z.real, 2000*z.imag] for z in poles]


def evaluate(q3d, lengths, hfss=None, shunts=None):
    p, _ = tune_lengths()
    p['L_filt'], p['L_res'] = lengths
    for key, cap in q3d['series_fF'].items():
        p[key] = (cap, *p[key][1:])
    if shunts is not None:
        p['cj'] = (p['cj'][0], *shunts)
    frequencies, capacitances, widths = modal_normalization(p)
    gain = mhz_per_ff(frequencies, capacitances)
    result = dict(
        lengths_mm=lengths, q3d_cj_fF=p['cj'][0], cj_local_shunts_fF=list(p['cj'][1:]),
        cj_shunt_source='supplied' if shunts is not None else 'Design5 localized Cj extraction',
        model_freq_GHz=frequencies.tolist(), tap_modal_C_fF=capacitances.tolist(),
        model_J_MHz=p['cj'][0]*gain, model_kappa_p_MHz=widths[0],
        mhz_per_fF=gain,
        model_hybrid_poles_GHz_and_FWHM_MHz=coupled_poles(p, frequencies, widths, p['cj'][0]*gain),
    )
    if hfss is not None:
        hfss_freqs = [hfss['filter_GHz'], hfss['resonator_GHz']]
        inverse_gain = mhz_per_ff(hfss_freqs, capacitances)
        result.update(hfss_freq_GHz=hfss_freqs, hfss_J_MHz=hfss['J_MHz'],
                      hfss_kappa_p_MHz=hfss['kappa_p_MHz'],
                      hfss_equivalent_Cj_fF=hfss['J_MHz']/inverse_gain)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q3d', type=Path)
    parser.add_argument('--lengths', nargs=2, type=float, required=True, metavar=('FILTER_MM','RES_MM'))
    parser.add_argument('--hfss', type=Path, help='S11/S21 pole fit JSON, optional')
    parser.add_argument('--cj-shunts', nargs=2, type=float, metavar=('FILTER_FF','RES_FF'))
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    data = evaluate(json.loads(args.q3d.read_text()), args.lengths,
                    json.loads(args.hfss.read_text()) if args.hfss else None, args.cj_shunts)
    formatted = json.dumps(data, indent=2)
    if args.out:
        args.out.write_text(formatted+'\n', encoding='utf-8')
    print(formatted)


if __name__ == '__main__':
    main()


