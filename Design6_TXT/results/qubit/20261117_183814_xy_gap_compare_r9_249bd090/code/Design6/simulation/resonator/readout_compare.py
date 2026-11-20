"""Overlay of HFSS readout-only S21 and the Q3D-parameter circuit model (re-analysis, no EM rerun).

Inputs: an HFSS readout_s21 run (raw/S_wide.csv), its pole fit (readout_poles summary.json) and a readout_model run
(run.json parameters: Q3D Cp / hammer values, geometry, lengths). The model S21 is recomputed with readout_model.port_Y
on exactly those parameters; nothing is fitted. Output: compare.png and compare.json (HFSS vs model bare fp / fr / J /
kappa_p, and the differences)."""
from argparse import ArgumentParser
from pathlib import Path
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'Design6/simulation/resonator')]
from readout_model import s21, NAMES          # noqa: E402
from readout_poles import load                # noqa: E402
from run_record import create_run, finish_run  # noqa: E402


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('hfss', type=Path); ap.add_argument('poles', type=Path); ap.add_argument('model', type=Path)
    a = ap.parse_args(); hf, po, mo = (p.resolve() for p in (a.hfss, a.poles, a.model))
    mrun = json.loads((mo / 'run.json').read_text(encoding='utf-8'))
    hrun = json.loads((hf / 'run.json').read_text(encoding='utf-8'))
    pm = mrun['parameters']
    assert abs(pm['res_extra_mm'] - hrun['parameters']['res_extra_mm']) < 1e-12, 'model and HFSS lengths differ'
    folder = create_run('Design6', 'resonator', 'readout_compare', entrypoint=__file__,
                        inputs=[hf / 'run.json', hf / 'raw/S_wide.csv', po / 'summary.json', mo / 'run.json',
                                mo / 'model.json'],
                        upstream=[str(hf), str(po), str(mo)],
                        parameters=dict(EM_rerun=False, res_extra_mm=pm['res_extra_mm'],
                                        model='readout_model.port_Y with the Q3D values of the model run, no fitting'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        lens = {q: tuple(pm['lengths_mm'][q]) for q in NAMES}
        args = (lens, pm['geometry'], pm['cp_fF'], pm['hammer_deembedded_fF'])
        f, sh = load(hf / 'raw/S_wide.csv')
        fm = np.arange(6.8, 7.2, .00002); sm = s21(fm, *args)
        hp = json.loads((po / 'summary.json').read_text(encoding='utf-8'))['channels']
        mp = json.loads((mo / 'model.json').read_text(encoding='utf-8'))['prediction_current']['bare']
        keys = ('bare_fp_GHz', 'bare_fr_GHz', 'detuning_fr_minus_fp_MHz', 'J_MHz', 'kappa_p_MHz')
        table = {q: {k: dict(hfss=hp[q]['wide_only'][k], model=mp[q][k],
                             model_minus_hfss=(mp[q][k] - hp[q]['wide_only'][k]) * (1000 if 'GHz' in k else 1),
                             diff_units='MHz') for k in keys} for q in NAMES}

        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        HF, MD = '#1f1f1f', '#2a6fdb'
        fig = plt.figure(figsize=(12, 8.2))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=.38, wspace=.22)
        for col, q in enumerate(('Q2', 'Q1')):
            ax = fig.add_subplot(gs[0, col]); h, m = table[q]['bare_fp_GHz'], table[q]['bare_fr_GHz']
            lo = min(h['hfss'], h['model']) - .045; hi = max(m['hfss'], m['model']) + .045
            ax.plot(f, 20 * np.log10(abs(sh)), color=HF, lw=2, label='HFSS driven (+0.1 mm, no claw)')
            ax.plot(fm, 20 * np.log10(abs(sm)), color=MD, lw=2, ls='--', label='Circuit model, Q3D Cp + hammer')
            for src, c, y in (('hfss', HF, -46), ('model', MD, -51)):
                for k, mk in (('bare_fp_GHz', 'v'), ('bare_fr_GHz', 'o')):
                    ax.plot(table[q][k][src], y, mk, color=c, ms=8, mec='white', mew=1.5, clip_on=False)
            ax.text(lo + .003, -40, 'bare  ▼ f_p   ● f_r', fontsize=8, color='#555')
            ax.set_xlim(lo, hi); ax.set_ylim(-54, 2); ax.set_title(f'{q} channel', fontsize=11)
            ax.set_xlabel('Frequency (GHz)'); ax.grid(alpha=.25, lw=.6)
            if col == 0:
                ax.set_ylabel('|S21| (dB)')
            for s in ('top', 'right'): ax.spines[s].set_visible(False)
        fig.legend(*fig.axes[0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(.5, .955), ncol=2,
                   fontsize=9, frameon=False)

        ax = fig.add_subplot(gs[1, 0])   # bare-frequency offsets from HFSS f_p of each channel
        for i, q in enumerate(('Q2', 'Q1')):
            ref = table[q]['bare_fp_GHz']['hfss']
            for src, c, dy in (('hfss', HF, .12), ('model', MD, -.12)):
                fp, fr = ((table[q][k][src] - ref) * 1000 for k in ('bare_fp_GHz', 'bare_fr_GHz'))
                ax.plot([fp, fr], [i + dy] * 2, color=c, lw=2)
                ax.plot(fp, i + dy, 'v', color=c, ms=9, mec='white', mew=1.5)
                ax.plot(fr, i + dy, 'o', color=c, ms=9, mec='white', mew=1.5)
                ax.annotate(f'{fr - fp:+.1f} MHz', ((fp + fr) / 2, i + dy), xytext=(0, 6 if dy > 0 else -13),
                            textcoords='offset points', ha='center', fontsize=8.5, color='#333')
        ax.axvline(0, color='#999', lw=.8)
        ax.set_yticks([0, 1], ['Q2', 'Q1']); ax.set_ylim(-.6, 1.6)
        ax.set_xlabel('Bare frequency − HFSS bare f_p (MHz)')
        ax.set_title('Bare f_p (▼) and f_r (●): HFSS (black) vs model (blue)', fontsize=10)
        ax.grid(axis='x', alpha=.25, lw=.6)
        for s in ('top', 'right'): ax.spines[s].set_visible(False)

        ax = fig.add_subplot(gs[1, 1]); ax.axis('off')
        rows = [['', 'Q1 HFSS', 'Q1 model', 'Q2 HFSS', 'Q2 model']]
        for k, lab, sc, fmt in (('bare_fp_GHz', 'bare f_p (GHz)', 1, '.4f'), ('bare_fr_GHz', 'bare f_r (GHz)', 1, '.4f'),
                                ('detuning_fr_minus_fp_MHz', 'f_r − f_p (MHz)', 1, '+.1f'),
                                ('J_MHz', 'J (MHz)', 1, '.1f'), ('kappa_p_MHz', 'κ_p (MHz)', 1, '.1f')):
            rows.append([lab] + [format(table[q][k][s] * sc, fmt) for q in NAMES for s in ('hfss', 'model')])
        t = ax.table(cellText=rows, loc='upper center', cellLoc='center', colWidths=[.28, .18, .18, .18, .18])
        t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.35)
        for (r, c), cell in t.get_celld().items():
            cell.set_edgecolor('#ddd')
            if r == 0: cell.set_text_props(weight='bold')
        cp, hm = pm['cp_fF'], pm['hammer_deembedded_fF']
        ax.text(0, .08, f"Model inputs (Q3D): Cp mutual {cp[0]:.2f} fF, ground {cp[1]:.1f} / {cp[2]:.1f} fF (not converged)\n"
                        f"hammer C_j {hm[0]:.3f} fF, excess ground {hm[1]:.2f} / {hm[2]:.2f} fF; parallel-line C_m omitted\n"
                        f"CPW Z0 {pm['geometry']['z0_ohm']:.2f} Ω, v {pm['geometry']['velocity_m_s']:.4g} m/s (conformal, no fit)",
                fontsize=8, color='#444', transform=ax.transAxes, va='bottom')
        fig.suptitle('Design6 r9 readout only (+0.1 mm L_res, no claw): HFSS vs Q3D-parameter circuit model', fontsize=12)
        fig.savefig(folder / 'compare.png', dpi=160, bbox_inches='tight')
        out = dict(hfss_run=str(hf), poles_run=str(po), model_run=str(mo), comparison=table,
                   note='HFSS values from the wide-grid 4-pole fit (one lossy mode, real J, no intrinsic loss); model values '
                        'from the same fit on the model S21. Validation status unknown.')
        (folder / 'compare.json').write_text(json.dumps(out, indent=1))
        finish_run(folder, status='completed', summary=out)
        print(json.dumps(table, indent=1))
    except Exception:
        import traceback
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()}); raise


if __name__ == '__main__':
    main()
