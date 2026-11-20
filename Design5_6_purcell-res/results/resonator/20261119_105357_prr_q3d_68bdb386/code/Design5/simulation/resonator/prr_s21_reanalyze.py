"""Re-analyse a finished Design5 PurcellReadout driven S21 run (new run folder, upstream = the raw run).

    python -I -B Design5/simulation/resonator/prr_s21_reanalyze.py --upstream Design5/results/resonator/<driven run>

Reads raw/<s_file> (S_wide.csv for wide-only runs) and the Q3D-model prediction stored in the upstream run.json.
Per channel window (split halfway between neighbouring target frequencies) prr_sim.s21_channel gives:
  zeros method  f_r = max |S21| between the two dips, f_p = z1 + z2 - f_r, J^2 = f_p f_r - z1 z2;
  pole fit      S21 = (c0 + c1 x)(f - z1)(f - z2) / ((f - p1)(f - p2)) with the zeros fixed, then the two-mode
                inversion of the poles -> f_p, f_r, J, kappa_p.
Dips are found in dB; the two most prominent per window are the channel zeros, weaker extra dips are reported.
Compared with the prediction (no fitting of lengths or capacitances).
"""
from pathlib import Path
import argparse
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator')]
import numpy as np  # noqa: E402
import prr_sim as S  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

DESIGN = 'Design5'


def windows(pred, f):
    names = sorted(pred, key=lambda k: pred[k]['f_r_GHz'])
    fr = [pred[q]['f_r_GHz'] for q in names]
    edges = [f.min()] + [(a + b) / 2 for a, b in zip(fr, fr[1:])] + [f.max()]
    return {q: (edges[i], edges[i + 1]) for i, q in enumerate(names)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--upstream', type=Path, required=True)
    a = ap.parse_args()
    up = a.upstream.resolve()
    rec = json.loads((up / 'run.json').read_text(encoding='utf-8'))
    raw = json.loads((up / 'raw' / 'driven.json').read_text(encoding='utf-8'))
    s_file = up / 'raw' / raw.get('s_file', 'S_fine.csv')
    pred = rec['parameters']['prediction']
    folder = create_run(DESIGN, 'resonator', 'prr_s21_reanalysis', entrypoint=__file__,
                        inputs=[s_file, up / 'run.json', up / 'raw' / 'driven.json'], upstream=[str(up)],
                        parameters=dict(method='prr_sim.s21_channel per target window (dB dips, two most prominent)',
                                        upstream_passes=raw['passes'], upstream_mesh_elements_note='see upstream raw',
                                        s_file=s_file.name))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        f, s = S.load_s(s_file)
        out, rows = {}, {}
        fig, ax = plt.subplots(figsize=(13, 4.8))
        ax.plot(f, 20 * np.log10(abs(s)), 'k-', lw=1.5, label='HFSS driven S21')
        for q, w in windows(pred, f).items():
            r = S.s21_channel(f, s, w)
            out[q] = dict(zeros=r['zeros'], poles=r['poles'], extra_dips=r['extra_dips'], window_GHz=list(w))
            p, pr = r['poles'], pred[q]
            rows[q] = dict(target_GHz=pr['f_r_GHz'], df_p_MHz=(p['f_p_GHz'] - pr['f_p_GHz']) * 1e3,
                           df_r_MHz=(p['f_r_GHz'] - pr['f_r_GHz']) * 1e3, J_MHz=p['J_MHz'],
                           J_model_MHz=pr['J_hammer_MHz'], J_ratio=p['J_MHz'] / pr['J_hammer_MHz'],
                           kappa_p_MHz=p['kappa_p_MHz'], kappa_p_model_MHz=pr['kappa_p_MHz'],
                           kappa_ratio=p['kappa_p_MHz'] / pr['kappa_p_MHz'], J_zeros_MHz=r['zeros']['J_MHz'],
                           fit_rms=p['fit_rms'])
            ax.plot(r['f'], 20 * np.log10(abs(r['fit'])), '--', lw=1, label=f'{q} pole fit')
            ax.axvline(pr['f_r_GHz'], color='0.7', lw=.6, ls=':')
        ax.set_xlabel('Frequency (GHz)'); ax.set_ylabel('|S21| (dB)'); ax.grid(alpha=.3)
        ax.legend(fontsize=7, ncol=4)
        ax.set_title(f"Design5 prr_r1, six channels: HFSS {raw['passes']} pass, wide sweep "
                     f"{raw.get('wide_step_GHz', 0) * 1e3:g} MHz; dotted = Q3D-model target", fontsize=9)
        fig.tight_layout(); fig.savefig(folder / 's21_channels.png', dpi=150); plt.close(fig)
        summary = dict(prediction=pred, hfss=out, comparison=rows,
                       validation=f"upstream fixed budget {raw['passes']} pass (user quick check); not converged")
        for q, v in sorted(rows.items(), key=lambda kv: kv[1]['target_GHz']):
            print(f"{q} {v['target_GHz']:.1f} GHz: df_p {v['df_p_MHz']:+.1f} df_r {v['df_r_MHz']:+.1f} MHz | "
                  f"J {v['J_MHz']:.2f} (model {v['J_model_MHz']:.2f}, x{v['J_ratio']:.2f}; zeros {v['J_zeros_MHz']:.2f}) | "
                  f"kappa_p {v['kappa_p_MHz']:.1f} (model {v['kappa_p_model_MHz']:.1f}, x{v['kappa_ratio']:.2f}) "
                  f"rms {v['fit_rms']:.3f}", flush=True)
        finish_run(folder, status='completed', summary=json.loads(json.dumps(summary, default=float)))
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
