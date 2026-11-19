"""Plot HFSS S21 against the Q3D-based results for Design5 PurcellReadout revisions (no solver; new analysis run).

    python -I -B Design5/simulation/resonator/purcell_resonator_pair_compare_plot.py --model <circuit run> [--model <circuit run> ...]

Each --model is a purcell_resonator_pair_circuit_model run; its upstreams give the HFSS S21 fit (purcell_resonator_pair_s21_reanalysis) and the raw S21.
Per revision (one column): |S21| HFSS vs the circuit model with C_in; kappa_p (matched-feed hanger, circuit with
C_in, HFSS); J (closed-form hammer, circuit, HFSS); f_p and f_r offsets from target (circuit, HFSS).
Q3D-based = lengths, C_j, head excess, C_p / C_in pi networks from Q3D; nothing fitted to HFSS.
"""
from pathlib import Path
import argparse
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator')]
import numpy as np  # noqa: E402
import purcell_resonator_pair_sim as S  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

DESIGN = 'Design5'


def load_case(model_run):
    m = model_run.resolve()
    rec = json.loads((m / 'run.json').read_text(encoding='utf-8'))
    fit_run = Path([u for u in rec['upstream'] if 's21_reanalysis' in u][0])
    fit = json.loads((fit_run / 'summary.json').read_text(encoding='utf-8'))
    raw_run = Path(json.loads((fit_run / 'run.json').read_text(encoding='utf-8'))['upstream'][0])
    raw = json.loads((raw_run / 'raw' / 'driven.json').read_text(encoding='utf-8'))
    rev = json.loads((raw_run / 'run.json').read_text(encoding='utf-8'))['parameters']['config']['revision']
    fh, sh = S.load_s(raw_run / 'raw' / raw.get('s_file', 'S_wide.csv'))
    d = np.loadtxt(m / 'model_S21_with_cin.csv', delimiter=',', skiprows=1)
    return dict(rev=rev, model=json.loads((m / 'summary.json').read_text(encoding='utf-8'))['comparison'],
                fit=fit, fh=fh, sh=sh, fm=d[:, 0], sm=d[:, 1] + 1j * d[:, 2], passes=raw['passes'],
                runs=[str(m), str(fit_run), str(raw_run)])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--model', type=Path, action='append', required=True)
    a = ap.parse_args()
    cases = [load_case(p) for p in a.model]
    folder = create_run(DESIGN, 'resonator', 'purcell_resonator_pair_compare_plot', entrypoint=__file__,
                        inputs=[Path(c['runs'][0]) / 'summary.json' for c in cases],
                        upstream=[r for c in cases for r in c['runs']],
                        parameters=dict(revisions=[c['rev'] for c in cases]))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        n = len(cases)
        fig, axs = plt.subplots(4, n, figsize=(7.2 * n, 15), squeeze=False,
                                gridspec_kw=dict(height_ratios=[1.3, 1, 1, 1]))
        table = {}
        for j, c in enumerate(cases):
            mod, fit = c['model'], c['fit']
            names = sorted(mod, key=lambda q: mod[q]['target_GHz'])
            x = np.arange(len(names))
            tgt = np.array([mod[q]['target_GHz'] for q in names])
            ax = axs[0, j]
            ax.plot(c['fh'], 20 * np.log10(abs(c['sh'])), 'k-', lw=1.3, label=f"HFSS S21 ({c['passes']} pass)")
            ax.plot(c['fm'], 20 * np.log10(abs(c['sm'])), 'C3-', lw=.9, alpha=.85, label='Q3D circuit model (with C_in)')
            for t in tgt:
                ax.axvline(t, color='0.75', ls=':', lw=.7)
            ax.set(xlabel='Frequency (GHz)', ylabel='|S21| (dB)', title=f"{c['rev']}: |S21| (dotted = targets)")
            ax.grid(alpha=.3); ax.legend(fontsize=8, loc='lower left')

            k_h = np.array([mod[q]['kappa_hanger_MHz'] for q in names])
            k_m = np.array([mod[q]['kappa_model_cin_MHz'] for q in names])
            k_s = np.array([mod[q]['kappa_hfss_MHz'] for q in names])
            ax = axs[1, j]
            w = .26
            ax.bar(x - w, k_h, w, color='0.75', label='matched-feed hanger (Q3D C_p)')
            ax.bar(x, k_m, w, color='C3', label='Q3D circuit, with C_in')
            ax.bar(x + w, k_s, w, color='k', label='HFSS fit')
            ax.set_ylim(0, 1.3 * max(k_h.max(), k_m.max(), k_s.max()))
            ax.set(ylabel='kappa_p / 2pi (MHz)', title='kappa_p'); ax.legend(fontsize=8, loc='upper left', ncol=3)
            ax.grid(alpha=.3, axis='y')

            j_c = np.array([fit['prediction'][q]['J_hammer_MHz'] for q in names])
            j_m = np.array([mod[q]['J_model_cin_MHz'] for q in names])
            j_s = np.array([mod[q]['J_hfss_MHz'] for q in names])
            bad = [i for i, q in enumerate(names) if j_s[i] <= 0 or fit['comparison'][q]['fit_rms'] > 0.01]
            ax = axs[2, j]
            ax.bar(x - w, j_c, w, color='0.75', label='closed-form hammer (Q3D C_j)')
            ax.bar(x, j_m, w, color='C3', label='Q3D circuit')
            ax.bar(x + w, np.where(j_s > 0, j_s, np.nan), w, color='k', label='HFSS fit')
            top = np.nanmax(np.r_[j_c, j_m, np.where(j_s > 0, j_s, np.nan)])
            for i in bad:
                ax.text(x[i] + w, 0.5 * top, 'HFSS\nfit\nfailed', ha='center', va='center', fontsize=7, color='C3')
            ax.set_ylim(0, 1.35 * top)
            ax.set(ylabel='J / 2pi (MHz)', title='J'); ax.legend(fontsize=8, loc='upper left', ncol=3)
            ax.grid(alpha=.3, axis='y')

            ax = axs[3, j]
            for key, lab, mk, col in (('f_p_model_cin_GHz', 'f_p circuit', 'o', 'C3'), ('f_r_model_cin_GHz', 'f_r circuit', 's', 'C3'),
                                      ('f_p_hfss_GHz', 'f_p HFSS', 'o', 'k'), ('f_r_hfss_GHz', 'f_r HFSS', 's', 'k')):
                v = (np.array([mod[q][key] for q in names]) - tgt) * 1e3
                ax.plot(x, v, mk, color=col, mfc='none' if 'f_r' in key else col, ms=7, label=lab)
            ax.axhline(0, color='0.5', lw=.8)
            ax.set(ylabel='offset from target (MHz)', title='f_p, f_r offsets (HFSS 1 pass: mesh error)')
            ax.legend(fontsize=8, ncol=2); ax.grid(alpha=.3)
            for r in (1, 2, 3):
                axs[r, j].set_xticks(x, [f'{q}\n{t:.1f} GHz' for q, t in zip(names, tgt)])
            table[c['rev']] = {q: dict(kappa_hanger=k_h[i], kappa_circuit=k_m[i], kappa_hfss=k_s[i], J_closed=j_c[i],
                                       J_circuit=j_m[i], J_hfss=j_s[i], J_hfss_fit_ok=i not in bad)
                               for i, q in enumerate(names)}
        fig.suptitle('Design5 PurcellReadout: HFSS driven S21 vs Q3D-based model (no HFSS input to the model)', fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, .98))
        fig.savefig(folder / 's21_vs_q3d.png', dpi=130); plt.close(fig)
        finish_run(folder, status='completed',
                   summary=json.loads(json.dumps(dict(table=table, figure='s21_vs_q3d.png',
                                                      validation='plot of existing runs; HFSS 1 pass'), default=float)))
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
