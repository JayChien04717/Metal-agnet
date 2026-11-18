"""Compare refined HFSS with frozen, previously exported Q3D circuit predictions."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import hammer_deembed as H
from hammer_fixed_s21 import model_plot_style
from run_record import create_run, finish_run


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('model_run', type=Path)
    ap.add_argument('refined_run', type=Path)
    args = ap.parse_args()
    original, refined = args.model_run.resolve(), args.refined_run.resolve()
    provenance = json.loads((original/'comparison.json').read_text())
    assert 'port_polarity' in provenance, 'Use a model run with verified HFSS port polarity'
    old = pd.read_csv(original/'model_comparison.csv')
    new = pd.read_csv(refined/'raw/S_fine.csv', index_col=0)
    f = new.index.to_numpy(float)
    assert len(f) == len(old) and np.allclose(f, old.frequency_GHz, atol=1e-10, rtol=0), 'Require identical grids; no interpolation'
    hfss = new.S21.map(complex).to_numpy()
    folder = create_run('Design4', 'resonator', 'hammer_refined_model', entrypoint=__file__,
        inputs=[original/'model_comparison.csv', original/'comparison.json', refined/'raw/S_fine.csv', refined/'raw/hfss_convergence.csv'],
        upstream=[str(original), str(refined)], parameters=dict(EM_rerun=False, fitting=False,
            method='Reuse exact circuit predictions on identical frequency grid; replace only HFSS spectrum'))
    print('RUN_FOLDER='+str(folder), flush=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    ax[0].plot(f, 20*np.log10(abs(hfss)), color='black', linestyle='-', linewidth=2, label='HFSS Driven S21 - full-wave')
    output = pd.DataFrame({'frequency_GHz': f, 'HFSS_S21': hfss})
    metrics = {}
    for col in old:
        if not col.endswith('_S21') or col == 'HFSS_S21':
            continue
        z = old[col].map(complex).to_numpy()
        output[col] = z
        label = col.removesuffix('_S21')
        metrics[label] = dict(rms_complex_S21=float(np.sqrt(np.mean(abs(z-hfss)**2))))
        ax[0].plot(f, 20*np.log10(abs(z)), **model_plot_style(label))
        ax[1].plot(f, abs(z-hfss), **model_plot_style(label))
    for a in ax:
        a.legend(); a.grid(alpha=.3)
    ax[0].set_ylabel('S21 (dB)'); ax[1].set_ylabel('|Q3D circuit - HFSS Driven|'); ax[1].set_xlabel('Frequency (GHz)')
    fig.suptitle('Design4: fixed physical line lengths\nQ3D curves: circuit predictions from extracted capacitances')
    fig.tight_layout(); fig.savefig(folder/'comparison.png', dpi=180)
    output.to_csv(folder/'model_comparison.csv', index=False)
    old_hfss = old.HFSS_S21.map(complex).to_numpy()
    summary = dict(models=metrics, HFSS_change_rms=float(np.sqrt(np.mean(abs(hfss-old_hfss)**2))),
        source_model=str(original), source_HFSS=str(refined),
        last_delta_s=float(pd.read_csv(refined/'raw/hfss_convergence.csv').iloc[-1]['Max Mag. Delta S']),
        limitations=provenance['limitations'])
    H.write(folder/'comparison.json', summary)
    finish_run(folder, status='completed', summary=summary)
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
