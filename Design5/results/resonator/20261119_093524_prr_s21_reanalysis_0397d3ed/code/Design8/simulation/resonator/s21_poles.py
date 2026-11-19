"""Re-analysis of a Design8 hfss_module driven run: bare f_p, f_r, J, kappa_p of the Purcell channel from S21.

Two independent extractions (note/purcell_filter_method.md section 2), no EM rerun, written to a new run that names
the driven run as upstream:
  zeros  S21 zeros z1, z2 (complex linear interpolation at each dip) and f_r = the S21 = 1 point between them
         (Wallraff: real zero of the branch admittance); f_p = Re(z1 + z2) - f_r, J^2 = f_p f_r - Re(z1 z2).
  poles  S21 = (c0 + c1 x) (f - z1)(f - z2) / ((f - p1)(f - p2)) over the fine band with the zeros fixed; the two
         complex poles (Im > 0) inverted with [[f_p + i kappa_p/2, J], [J, f_r]] (one lossy mode, real J).
The zeros method gives no kappa_p; the poles method does.

    python -I -B Design8/simulation/resonator/s21_poles.py <driven run folder>
"""
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(Path(__file__).resolve().parent)]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402
from scipy.signal import find_peaks  # noqa: E402

from run_record import create_run, finish_run  # noqa: E402


def load(path):
    S = pd.read_csv(path, index_col=0)
    f = S.index.to_numpy(float)
    return f, S['S21'].map(lambda v: complex(str(v).replace('i', 'j'))).to_numpy()


def zeros(f, s):
    idx, _ = find_peaks(-abs(s), prominence=0.3)
    if len(idx) != 2:
        raise RuntimeError(f'expected two dips, found {len(idx)} at {f[idx]}')
    out = []
    for i in idx:
        sl = slice(max(i - 3, 0), i + 4)
        b, a = np.polyfit(f[sl] - f[i], s[sl], 1)            # s ~ a + b (f - f_i)
        out.append(f[i] - a / b)
    i0, i1 = idx
    k = i0 + int(np.argmax(abs(s[i0:i1])))
    sl = slice(k - 3, k + 4)
    c = np.polyfit(f[sl] - f[k], abs(s[sl]), 2)
    fr = f[k] - c[1] / (2 * c[0])
    return out, fr, idx


def pole_fit(f, s, z):
    c = f.mean()
    x = (f - c) * 1e3
    zz = [(zi - c) * 1e3 for zi in z]
    num = (x - zz[0]) * (x - zz[1])

    def model(q):
        p1, p2 = q[0] + 1j * q[1], q[2] + 1j * q[3]
        shape = num / ((x - p1) * (x - p2))
        B = np.column_stack([shape, shape * x / 100])
        coef = np.linalg.lstsq(B, s, rcond=None)[0]
        return B @ coef

    def res(q):
        r = model(q) - s
        return np.r_[r.real, r.imag]
    best = None
    for sgn in (1, -1):
        mid = np.mean([zi.real for zi in zz])
        q0 = [zz[0].real, sgn * 10, zz[1].real, sgn * 10]
        r = least_squares(res, q0, max_nfev=20000)
        if best is None or r.cost < best.cost:
            best = r
    q = best.x
    poles = [c + complex(q[0], q[1]) / 1e3, c + complex(q[2], q[3]) / 1e3]
    rms = float(np.sqrt(np.mean(abs(model(q) - s) ** 2)))
    return poles, rms, model(q)


def invert(z1, z2):
    s, p = z1 + z2, z1 * z2
    fr = p.imag / s.imag
    fp = s.real - fr
    j2 = fp * fr - p.real
    return dict(f_p_GHz=fp, f_r_GHz=fr, J_MHz=float(np.sign(j2) * np.sqrt(abs(j2)) * 1e3),
                kappa_p_MHz=2 * abs(s.imag) * 1e3)


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('driven', type=Path)
    a = ap.parse_args()
    up = a.driven.resolve()
    folder = create_run('Design8', 'resonator', 's21_poles', entrypoint=__file__,
                        inputs=[up / 'raw' / 'S_fine.csv', up / 'run.json'], upstream=[str(up)],
                        parameters=dict(method='zeros + S21=1 peak; pole fit with fixed zeros', fitted_lengths=False))
    try:
        f, s = load(up / 'raw' / 'S_fine.csv')
        z, fr0, idx = zeros(f, s)
        fp0 = sum(zi.real for zi in z) - fr0
        j2 = fp0 * fr0 - (z[0] * z[1]).real
        zres = dict(zeros_GHz=[[zi.real, zi.imag] for zi in z], f_r_GHz=fr0, f_p_GHz=fp0,
                    J_MHz=float(np.sqrt(abs(j2)) * 1e3), J2_sign=int(np.sign(j2)))
        poles, rms, fit = pole_fit(f, s, z)
        if np.mean([p.imag for p in poles]) < 0:
            poles = [p.conjugate() for p in poles]
        pres = dict(poles_GHz=[[p.real, p.imag] for p in poles], fit_rms=rms, **invert(*poles))
        pred = json.loads((up / 'run.json').read_text(encoding='utf-8'))['parameters']['prediction']
        drv = json.loads((up / 'raw' / 'driven.json').read_text(encoding='utf-8'))
        cmp = {m: dict(df_p_MHz=(r['f_p_GHz'] - pred['f_p_GHz']) * 1e3, df_r_MHz=(r['f_r_GHz'] - pred['f_r_GHz']) * 1e3,
                       J_ratio_hfss_over_hammer=r['J_MHz'] / pred['J_hammer_MHz'])
               for m, r in (('zeros', zres), ('poles', pres))}
        cmp['poles']['kappa_p_ratio'] = pres['kappa_p_MHz'] / pred['kappa_p_MHz']
        fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
        ax[0].plot(f, 20 * np.log10(abs(s)), 'k-', lw=2, label='HFSS Driven S21')
        ax[0].plot(f, 20 * np.log10(abs(fit)), 'r--', lw=1, label='two-pole fit (zeros fixed)')
        for zi in z:
            ax[0].axvline(zi.real, color='tab:blue', ls=':', lw=1)
        ax[0].axvline(fr0, color='tab:green', ls=':', lw=1, label='S21 = 1 (f_r)')
        ax[0].set_ylabel('|S21| (dB)'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
        ax[1].plot(f, abs(fit - s), 'r-'); ax[1].set_ylabel('|fit - HFSS|'); ax[1].set_xlabel('Frequency (GHz)')
        ax[1].grid(alpha=.3)
        ax[0].set_title(f"poles: f_p {pres['f_p_GHz']:.4f}, f_r {pres['f_r_GHz']:.4f} GHz, J {pres['J_MHz']:.2f} MHz, "
                        f"kappa_p {pres['kappa_p_MHz']:.1f} MHz (HFSS {drv['passes']} passes, dS {drv['last_delta_s']:.3g})",
                        fontsize=8)
        fig.tight_layout(); fig.savefig(folder / 's21_poles.png', dpi=150); plt.close(fig)
        summary = dict(upstream=str(up), driven=dict(passes=drv['passes'], last_delta_s=drv['last_delta_s']),
                       prediction=pred, zeros=zres, poles=pres, comparison=cmp)
        print(json.dumps(summary, indent=1, default=float))
        finish_run(folder, status='completed', summary=summary)
    except Exception:
        import traceback
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
