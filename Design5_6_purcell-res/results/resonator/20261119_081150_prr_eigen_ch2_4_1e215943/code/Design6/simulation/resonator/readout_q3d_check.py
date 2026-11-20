"""Static check of the Design6 readout circuit model against the full readout Q3D (re-analysis, no EM rerun).

The circuit model (readout_model.port_Y) reduces at f -> 0, with the shorted ends opened as in the full Q3D, to a
5-net capacitance matrix: each TL contributes C' * length to ground, each pi network its mutual and ground parts,
and (optionally) coupled-line sections c'm * overlap between filter and resonator. That matrix is compared entry by
entry with the full Q3D Maxwell matrix (readout_q3d_full.py). Differences show which lumped element or line constant
of the model is wrong, independent of HFSS.

Bookkeeping per net (lengths are the model's TL lengths, mm):
  feed      : C' (feed_in + feed_mid + feed_out + 2 stub_p) + 2 Cp ground (frame side)
  filter_Q  : C' (L_filter + spacing/2) + Cp ground (lead side) + hammer excess ground (filter side)
  res_Q     : C' (L_res + spacing/2) + hammer excess ground (res side)
  feed-filter_Q : Cp mutual;  filter_Q-res_Q : hammer mutual (+ c'm * parallel overlap)
"""
from argparse import ArgumentParser
from pathlib import Path
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code/general'))
from run_record import create_run, finish_run  # noqa: E402

ORDER = ['feed', 'filter_Q1', 'res_Q1', 'filter_Q2', 'res_Q2']
OVERLAP_MM = 1.636   # filter-resonator side-by-side run per channel at 200 um (y 0.064..1.700 mm, layout geometry)


def model_matrix(pm, c_line, cm_per_mm=0.0):
    F, sp = pm['geometry']['feeds_mm'], pm['geometry']['spacing_mm']
    cp, hm = pm['cp_fF'], pm['hammer_deembedded_fF']
    M = np.zeros((5, 5)); g = np.zeros(5); mut = {}
    g[0] = c_line * (F['feed_in'] + F['feed_mid'] + F['feed_out'] + F['stub_p_Q1'] + F['stub_p_Q2']) + 2 * cp[1]
    for k, q in ((1, 'Q1'), (3, 'Q2')):
        Lf, Lr = pm['lengths_mm'][q]
        g[k] = c_line * (Lf + sp[q] / 2) + cp[2] + hm[1]
        g[k + 1] = c_line * (Lr + sp[q] / 2) + hm[2]
        mut[(0, k)] = cp[0]; mut[(k, k + 1)] = hm[0] + cm_per_mm * OVERLAP_MM
    for (i, j), c in mut.items():
        M[i, j] = M[j, i] = -c
    for i in range(5):
        M[i, i] = g[i] - M[i].sum()
    return M


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('q3d', type=Path); ap.add_argument('model', type=Path)
    a = ap.parse_args(); q3, mo = a.q3d.resolve(), a.model.resolve()
    pm = json.loads((mo / 'run.json').read_text(encoding='utf-8'))['parameters']
    qs = json.loads((q3 / 'summary_out.json').read_text())
    assert qs['order'] == ORDER
    folder = create_run('Design6', 'resonator', 'readout_q3d_check', entrypoint=__file__,
                        inputs=[q3 / 'run.json', q3 / 'summary_out.json', mo / 'run.json'], upstream=[str(q3), str(mo)],
                        parameters=dict(EM_rerun=False, overlap_mm=OVERLAP_MM))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        Q = np.array(qs['maxwell_fF'])
        geo = pm['geometry']; c_conf = 1e12 / (geo['z0_ohm'] * geo['velocity_m_s'])
        variants = dict(model_as_run=model_matrix(pm, c_conf),
                        with_parallel_cm=model_matrix(pm, c_conf, 0.56341))
        out = dict(q3d_run=str(q3), model_run=str(mo), c_line_conformal_fF_per_mm=c_conf, order=ORDER, q3d=Q.tolist(),
                   q3d_converged=qs['solver_converged'], q3d_last_change=qs['last_pass_max_relative_change'], variants={})
        lab = lambda i, j: ORDER[i] if i == j else f'{ORDER[i]}-{ORDER[j]}'
        for name, M in variants.items():
            rows = {}
            for i in range(5):
                gq, gm = Q[i].sum(), M[i].sum()
                rows[f'{ORDER[i]} ground'] = dict(q3d=gq, model=gm, model_minus_q3d=gm - gq,
                                                  equiv_line_um=(gm - gq) / c_conf * 1000)
                for j in range(i + 1, 5):
                    if abs(Q[i, j]) > .05 or abs(M[i, j]) > 0:
                        rows[lab(i, j)] = dict(q3d=-Q[i, j], model=-M[i, j], model_minus_q3d=-M[i, j] + Q[i, j])
            out['variants'][name] = dict(model=M.tolist(), compare=rows)
        (folder / 'check.json').write_text(json.dumps(out, indent=1))
        finish_run(folder, status='completed', summary=out)
        for name, v in out['variants'].items():
            print('==', name)
            for k, r in v['compare'].items():
                extra = f"  ({r['equiv_line_um']:+.1f} um of line)" if 'equiv_line_um' in r else ''
                print(f"{k:24s} Q3D {r['q3d']:9.3f}  model {r['model']:9.3f}  diff {r['model_minus_q3d']:+8.3f} fF{extra}")
    except Exception:
        import traceback
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()}); raise


if __name__ == '__main__':
    main()
