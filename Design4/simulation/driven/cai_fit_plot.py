"""Ideal S21 from the Q3D capacitances vs HFSS driven S21 / S11, and the capacitor-model fit to S21 + S11.

Curves (one figure, |S21| and |S11|, full span + zoom):
  HFSS         driven S.csv (S11, S21)
  Q3D ideal    circuit model, drawn lengths (P.LENGTHS), C_in / C_p / C_j = Q3D series C   (no fitting)
  fit          same circuit, L_filt, L_res, C_in, C_p, C_j free, least squares on |S21| and |S11| (linear)
Derived per curve: poles -> f_p, f_r, kappa_p, J (2-mode inversion, point_J).
The pi-network ground capacitances are the Design5 Q3D values used by purcell_6500.tune_lengths (base).

Usage: python driven/cai_fit_plot.py <S21 run folder> <Q3D caps json>
Output: <run>/cai_fit.png, <run>/cai_fit.json
"""
import json, os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np
from scipy.optimize import least_squares

SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DESIGN = os.path.dirname(SIM); ROOT = os.path.dirname(DESIGN)
sys.path[:0] = [os.path.join(SIM, 'eigenmode'), os.path.join(DESIGN, 'pattern design'), os.path.join(ROOT, 'Design5', 'simulation')]
import purcell_6500 as P
import theory5
import s21_model_vs_eigen as M
from hfss_hammer_straps import point_J


def load_S(run):
    import pandas as pd
    d = pd.read_csv(os.path.join(run, 'S.csv'), index_col=0)
    f = d.index.values.astype(float)
    col = lambda name: d[[c for c in d.columns if name in c][0]].astype(complex).values
    return f, col('S11'), col('S21')


def params(base, x):
    """x = [L_filt, L_res, C_in, C_p, C_j]."""
    p = M.circuit(base, x[3], x[4], x[0], x[1])
    p['cin'] = (x[2],) + tuple(base['cin'][1:])
    return p


def derived(p, f, s21):
    i = int(np.argmin(np.abs(s21)))
    pz = M.poles(p, [(f[i] + 0.008, 40.0), (f[i] - 0.008, 5.0)])
    inv = point_J([z[0] for z in pz], [z[1] / 2e3 for z in pz])
    return dict(poles=pz, f_p_GHz=inv['f_p_GHz'], f_r_GHz=inv['f_r_GHz'], kappa_p_MHz=inv['kappa_p_MHz'], J_MHz=inv['J_MHz'])


def main(run, q3d_json):
    f, s11, s21 = load_S(run)
    q = json.load(open(q3d_json))['series_fF']
    base, _ = P.tune_lengths()
    xq = [P.LENGTHS[0], P.LENGTHS[1], q['cin'], q['cp'], q['cj']]
    model = lambda x, ff: theory5.s21_theory(ff, params(base, x))

    def r(x):
        m11, m21 = model(x, f)
        return np.concatenate([np.abs(m21) - np.abs(s21), np.abs(m11) - np.abs(s11)])
    best = None
    for dl in (-0.06, -0.03, 0.0, 0.03, 0.06):           # start lengths around the drawn ones (HFSS shift)
        x0 = [xq[0] + dl, xq[1] + dl, *xq[2:]]
        sl = least_squares(r, x0, bounds=([x0[0] - .15, x0[1] - .15, 20, 20, .2], [x0[0] + .15, x0[1] + .15, 80, 80, 8]),
                           x_scale=[1e-3, 1e-3, 1, 1, .05], diff_step=1e-6)
        if best is None or sl.cost < best.cost:
            best = sl
    xf = list(best.x)
    mq11, mq21 = model(xq, f); mf11, mf21 = model(xf, f)
    dq, df = derived(params(base, xq), f, mq21), derived(params(base, xf), f, mf21)
    res = dict(run=run, q3d_json=q3d_json,
               q3d=dict(L_mm=xq[:2], C_fF=dict(cin=xq[2], cp=xq[3], cj=xq[4]), **dq),
               fit=dict(L_mm=xf[:2], C_fF=dict(cin=xf[2], cp=xf[3], cj=xf[4]), **df,
                        rms=float(np.sqrt(np.mean(best.fun ** 2)))),
               cai2024=dict(f_GHz=6.08, C_fF=dict(cin=40, cp=40, cj=2), J_MHz=8.8, kappa_p_MHz=45.5))
    json.dump(res, open(os.path.join(run, 'cai_fit.json'), 'w'), indent=1, default=float)

    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    db = lambda z: 20 * np.log10(np.abs(z))
    c0 = f[int(np.argmin(np.abs(s21)))]
    fig, axs = plt.subplots(2, 2, figsize=(13, 8.5))
    for col, span in enumerate((None, 0.06)):
        for row, (h, mq, mf, name) in enumerate(((s21, mq21, mf21, 'S21'), (s11, mq11, mf11, 'S11'))):
            ax = axs[row, col]
            ax.plot(f, db(h), 'k', lw=2.6, alpha=.35, label=f'HFSS driven {name}')
            ax.plot(f, db(mq), 'C0--', lw=1.2, label=f'ideal: Q3D C, drawn lengths')
            ax.plot(f, db(mf), 'C3', lw=1.2, label=f'capacitor model fitted to |S21|,|S11|')
            if span:
                ax.set_xlim(c0 - span, c0 + span)
            ax.set_xlabel('f (GHz)'); ax.set_ylabel(f'|{name}| (dB)'); ax.grid(alpha=.3)
    axs[0, 0].legend(fontsize=8, loc='lower left')
    a, b = res['q3d'], res['fit']
    rows = [('C_in (fF)', 40, a['C_fF']['cin'], b['C_fF']['cin']), ('C_p (fF)', 40, a['C_fF']['cp'], b['C_fF']['cp']),
            ('C_j (fF)', 2, a['C_fF']['cj'], b['C_fF']['cj']), ('L_filt (mm)', None, a['L_mm'][0], b['L_mm'][0]),
            ('L_res (mm)', None, a['L_mm'][1], b['L_mm'][1]), ('f_p (GHz)', 6.08, a['f_p_GHz'], b['f_p_GHz']),
            ('f_r (GHz)', 6.08, a['f_r_GHz'], b['f_r_GHz']), ('κ_p/2π (MHz)', 45.5, a['kappa_p_MHz'], b['kappa_p_MHz']),
            ('J/2π (MHz)', 8.8, a['J_MHz'], b['J_MHz'])]
    txt = f"{'':13s}{'Cai':>7s}{'Q3D':>9s}{'fit':>9s}\n" + '\n'.join(
        f"{n:13s}{('' if c is None else f'{c:g}'):>7s}{u:9.3f}{v:9.3f}" for n, c, u, v in rows)
    axs[1, 1].text(1.02, 0.5, txt, transform=axs[1, 1].transAxes, family='monospace', fontsize=8.5, va='center')
    fig.suptitle('Design4 (Cai 2024 parameters, CPW 16/9, 400 µm spacing, bonds): '
                 'ideal S21 from Q3D C vs HFSS driven, and C fitted to S21 + S11', fontsize=10)
    fig.tight_layout(rect=(0, 0, 0.84, 0.96)); out = os.path.join(run, 'cai_fit.png'); fig.savefig(out, dpi=140)
    print(txt); print('figure:', out)


if __name__ == '__main__':
    main(*sys.argv[1:3])
