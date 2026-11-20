"""Design6_TXT readout circuit model corrected with the full readout Q3D, then checked against HFSS S21 (no EM rerun).

Starting point: readout_model.py (uncoupled lambda/4 lines, Q3D Cp and hammer coupons). The full readout Q3D
(readout_q3d_full.py, same geometry as the HFSS reference) and its static check (readout_q3d_check.py) show:
  * filter-resonator static mutual 2.19 / 2.23 fF, vs hammer 0.66 fF (+0.92 fF from c'm * 1.636 mm parallel run);
  * a direct feed-resonator mutual of 0.235 fF (Cp frame <-> resonator tail at y = 1.70 mm), absent from the model;
  * filter ground C in the model a few fF above the resonator's, after removing a common per-length offset
    (model C' conformal vs Q3D): attributed to the Cp lead-side ground (coupon reference plane).
Every correction below is taken from Q3D only; HFSS is used only for the comparison.

Variants (cumulative unless stated):
  A  readout_model as run (asserted identical admittance for cm = 0, no extras)
  B  + coupled-line sections, c'm = 0.56341 fF/mm (Q3D pair, Design4_Single_purcell-res e1d28295), homogeneous L'/C' coupling
  C  + feed-res 0.235 fF (Cp side a <-> resonator node at the top of the parallel run) + Cp lead ground - delta_q
  D  C + the rest of the Q3D filter-res static mutual as a lumped C at the hammer (total = Q3D)
  E  C with c'm scaled so that hammer + c'm * overlap = Q3D (alternative placement of the same missing mutual)
  F  D with the Q3D pair C' (165.8 fF/mm) instead of the conformal C' and L' kept conformal (v, Z0 rescaled)
  G  F with the C_p mutual taken from the full Q3D (42.17 fF, coupon 41.25 fF); the C_p ground parts unchanged
"""
from argparse import ArgumentParser
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'Design6_TXT/simulation/resonator')]
import readout_model as M        # noqa: E402  (sets sys.path)
from readout_network import port_Y, FILTER_TAP, RES_TAP, TEE_HALF
from readout_poles import invert, fit  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402
from scipy.signal import find_peaks, peak_widths  # noqa: E402

P, NAMES = M.P, M.NAMES
assert (FILTER_TAP, RES_TAP, TEE_HALF) == (P.FILTER_TAP, P.RES_TAP, P.TEE_HALF)
PAIR_RUN = ROOT / 'Design4_Single_purcell-res/results/resonator/20261116_220124_hammer_coupled_model_e1d28295'
CM_PAIR, CSELF_PAIR = 0.56341, 165.845          # fF/mm, Q3D two-CPW pair at 200 um (mean of both ground per mm)
ORDER = ['feed', 'filter_Q1', 'res_Q1', 'filter_Q2', 'res_Q2']


def sections(geo_json, q):
    """Parallel-run lengths (mm) of channel q from the rendered geometry: straight vertices, fillets ignored."""
    from shapely import wkt
    ids = {n: i for i, n in enumerate(geo_json['options'], start=1)}
    lines = {r['component']: list(wkt.loads(r['wkt']).coords) for r in geo_json['path'] if not r['subtract']}
    ft, ro, rt = (lines[ids[f'{n}_{q}']] for n in ('filter_tail', 'res_open', 'res_tail'))
    tee_lo, tee_hi = ro[0][1], rt[0][1]
    assert np.isclose(ft[0][1], tee_lo) and np.isclose(tee_hi - tee_lo, 2 * P.TEE_HALF)
    above = rt[1][1] - tee_hi                      # filter feed side beside the resonator tail
    below = tee_lo - max(ft[1][1], ro[1][1])       # filter tail beside the resonator open side
    return dict(above_mm=above, below_mm=below, filter_free_mm=P.FILTER_TAP - P.TEE_HALF - above,
                res_open_free_mm=P.RES_TAP - P.TEE_HALF - below, overlap_mm=above + below + 2 * P.TEE_HALF)


def s21(grid, *args, **kw):
    out = []
    for x in grid:
        Yp = port_Y(x, *args, **kw); S = np.linalg.solve(np.eye(2) + 50 * Yp, np.eye(2) - 50 * Yp); out.append(S[1, 0])
    return np.array(out)


def bare_of(*args, **kw):
    """Fresh dip detection + four-pole fit on the model S21 (readout_model.poles_of procedure)."""
    G = M.GRID
    S = s21(G, *args, **kw); db = 20 * np.log10(abs(S))
    idx, _ = find_peaks(-db, prominence=.3)
    width = peak_widths(-db, idx, rel_height=.5)[0] * .00002
    assert len(idx) == 4, f'dips {G[idx]}, widths {width}'
    import itertools
    best = None      # deep dips make the dB width meaningless for hybridized pairs: try several width guesses
    for wq2, wq1 in itertools.product(((10.0, .2), (.2, 10.0), (8.0, 8.0)), repeat=2):
        g = [(G[idx][0], wq2[0]), (G[idx][1], wq2[1]), (G[idx][2], wq1[0]), (G[idx][3], wq1[1])]
        try:
            r = fit(G, S, g, (G[0], G[-1]))
        except Exception:
            continue
        if best is None or r[2] < best[2]:
            best = r
    poles, _, rms = best
    assert rms < 1e-3, f'model pole fit rms {rms}'
    out = {}
    for k, q in ((0, 'Q2'), (2, 'Q1')):
        out[q] = invert(*sorted(poles[k:k + 2], key=lambda z: -z.imag))
    return out, S


def ground_delta(Q, model, pm):
    """Filter ground excess of the model beyond the per-length offset seen on the resonator of the same channel."""
    out = {}
    for k, q in ((1, 'Q1'), (3, 'Q2')):
        Lf, Lr = pm['lengths_mm'][q]; sp = pm['geometry']['spacing_mm'][q] / 2
        dgf, dgr = model[k].sum() - Q[k].sum(), model[k + 1].sum() - Q[k + 1].sum()
        out[q] = dict(filter_excess_fF=dgf, res_excess_fF=dgr, common_per_mm_fF=dgr / (Lr + sp),
                      delta_fF=dgf - dgr * (Lf + sp) / (Lr + sp))
    return out


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('q3d', type=Path, help='readout_q3d_full run'); ap.add_argument('model', type=Path)
    ap.add_argument('poles', type=Path, help='readout_poles run of the HFSS reference')
    ap.add_argument('hfss', type=Path, help='HFSS readout_s21 run (S21 for the figure)')
    a = ap.parse_args(); q3, mo, po, hf = (p.resolve() for p in (a.q3d, a.model, a.poles, a.hfss))
    sys.path.insert(0, str(Path(__file__).parent)); import readout_q3d_check as C
    pm = json.loads((mo / 'run.json').read_text(encoding='utf-8'))['parameters']
    qs = json.loads((q3 / 'summary_out.json').read_text()); assert qs['order'] == ORDER
    geo_json = json.loads((q3 / 'raw/geometry.json').read_text())
    hp = json.loads((po / 'summary.json').read_text(encoding='utf-8'))['channels']
    folder = create_run('Design6_TXT', 'resonator', 'readout_model_q3d', entrypoint=__file__,
                        inputs=[q3 / 'run.json', q3 / 'summary_out.json', q3 / 'raw/geometry.json', mo / 'run.json',
                                po / 'summary.json', hf / 'raw/S_wide.csv', PAIR_RUN / 'summary.json'],
                        upstream=[str(q3), str(mo), str(po), str(hf), str(PAIR_RUN)],
                        parameters=dict(EM_rerun=False, cm_pair_fF_per_mm=CM_PAIR, cself_pair_fF_per_mm=CSELF_PAIR,
                                        coupling='homogeneous quasi-TEM: L matrix = (v^2 C_matrix)^-1, k_L = k_C',
                                        hfss_used_for='comparison only, no parameter is fitted to HFSS'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        geo, cp, hm = pm['geometry'], pm['cp_fF'], pm['hammer_deembedded_fF']
        lens = {q: tuple(pm['lengths_mm'][q]) for q in NAMES}
        sec = {q: sections(geo_json, q) for q in NAMES}
        c_conf = 1e12 / (geo['z0_ohm'] * geo['velocity_m_s'])
        # identity check vs readout_model
        for f in (6.9 + .001j, 7.1 + .002j):
            assert np.allclose(port_Y(f, lens, geo, cp, hm, sec, coupled=False), M.port_Y(f, lens, geo, cp, hm))
        Q = np.array(qs['maxwell_fF'])
        base = C.model_matrix(pm, c_conf)
        dlt = ground_delta(Q, base, pm)
        # stability of delta over the last Q3D passes (pass CSVs are in F)
        inv = {v: k for k, v in qs['net_object_names'].items()}; o = [inv[x] for x in ORDER]
        dpass = {}
        for k in range(10, qs['passes_saved'] + 1):
            Mk = pd.read_csv(q3 / f'raw/pass_{k}.csv', index_col=0).loc[o, o].to_numpy(float) * 1e15
            dpass[k] = {q: v['delta_fF'] for q, v in ground_delta(Mk, base, pm).items()}
        feed_res = float(np.mean([-Q[0, 2], -Q[0, 4]]))
        q3d_mut = {q: -Q[k, k + 1] for k, q in ((1, 'Q1'), (3, 'Q2'))}
        mean = lambda d: float(np.mean(list(d.values())))
        delta = mean({q: v['delta_fF'] for q, v in dlt.items()})
        ov = mean({q: s['overlap_mm'] for q, s in sec.items()})
        extra = mean(q3d_mut) - hm[0] - CM_PAIR * ov
        cm_scaled = (mean(q3d_mut) - hm[0]) / ov
        cp_c = [cp[0], cp[1], cp[2] - delta]
        cp_full = float(np.mean([-Q[0, 1], -Q[0, 3]])); cp_g = [cp_full, cp[1], cp[2] - delta]
        sc = np.sqrt(c_conf / CSELF_PAIR)
        geo_f = dict(geo, z0_ohm=geo['z0_ohm'] * sc, velocity_m_s=geo['velocity_m_s'] * sc)
        variants = dict(
            A=(dict(cp=cp, geo=geo), dict(coupled=False)),
            B=(dict(cp=cp, geo=geo), dict(cm=CM_PAIR)),
            C=(dict(cp=cp_c, geo=geo), dict(cm=CM_PAIR, feed_res=feed_res)),
            D=(dict(cp=cp_c, geo=geo), dict(cm=CM_PAIR, feed_res=feed_res, extra_cj=extra)),
            E=(dict(cp=cp_c, geo=geo), dict(cm=cm_scaled, feed_res=feed_res)),
            F=(dict(cp=cp_c, geo=geo_f), dict(cm=CM_PAIR, feed_res=feed_res, extra_cj=extra)),
            G=(dict(cp=cp_g, geo=geo_f), dict(cm=CM_PAIR, feed_res=feed_res, extra_cj=extra)))
        keys = ('bare_fp_GHz', 'bare_fr_GHz', 'detuning_fr_minus_fp_MHz', 'J_MHz', 'kappa_p_MHz')
        res, curves = {}, {}
        for name, (pa, kw) in variants.items():
            b, S = bare_of(lens, pa['geo'], pa['cp'], hm, sec, **kw)
            curves[name] = S
            res[name] = {q: {k: dict(model=b[q][k], hfss=hp[q]['wide_only'][k],
                                     diff=(b[q][k] - hp[q]['wide_only'][k]) * (1000 if k.endswith('GHz') else 1))
                             for k in keys} for q in NAMES}
            print(name, ' | '.join(f"{q}: dfp {r['bare_fp_GHz']['diff']:+6.1f} dfr {r['bare_fr_GHz']['diff']:+6.1f} "
                                   f"J {r['J_MHz']['model']:5.2f}/{r['J_MHz']['hfss']:5.2f} "
                                   f"kp {r['kappa_p_MHz']['model']:5.1f}/{r['kappa_p_MHz']['hfss']:5.1f}"
                                   for q, r in res[name].items()), flush=True)
        np.savez(folder / 'raw/model_s21.npz', f_GHz=M.GRID, **curves)
        out = dict(q3d_run=str(q3), model_run=str(mo), poles_run=str(po), hfss_run=str(hf), sections_mm=sec,
                   corrections=dict(ground_delta=dlt, ground_delta_by_pass=dpass, delta_used_fF=delta,
                                    feed_res_fF=feed_res, q3d_filter_res_mutual_fF=q3d_mut,
                                    extra_lumped_cj_fF=extra, cm_scaled_fF_per_mm=cm_scaled,
                                    line_scale_F=float(sc), cp_mutual_full_q3d_fF=cp_full),
                   q3d_converged=qs['solver_converged'], q3d_last_pass_max_relative_change=qs['last_pass_max_relative_change'],
                   variants={k: res[k] for k in variants}, units='GHz for frequencies, MHz for diff / J / kappa',
                   note='HFSS values from the wide-grid 4-pole fit of the HFSS reference; model values from the same fit '
                        'on the model S21. diff = model - HFSS (MHz). Validation status unknown.')
        (folder / 'compare.json').write_text(json.dumps(out, indent=1, default=float))
        figure(folder, hf, curves, res)
        finish_run(folder, status='completed', summary=dict(corrections=out['corrections'], variants=out['variants']))
        print(json.dumps(out['corrections'], indent=1, default=float))
    except Exception:
        import traceback
        (folder / 'error.log').write_text(traceback.format_exc())
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()}); raise


def figure(folder, hf, curves, res):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    hs = pd.read_csv(hf / 'raw/S_wide.csv', index_col=0)
    fh, sh = hs.index.to_numpy(float), hs['S21'].map(lambda x: complex(x.strip('()'))).to_numpy()
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw=dict(width_ratios=[1, 1, 1.1]))
    win = dict(Q2=(6.82, 7.00), Q1=(7.03, 7.20))
    sty = dict(A=('#9aa0a6', '--'), B=('#7a9cc6', ':'), G=('#d4731c', '-'))
    for axx, q in zip(ax[:2], ('Q2', 'Q1')):
        axx.plot(fh, 20 * np.log10(abs(sh)), color='k', lw=2, label='HFSS')
        for n, (c, ls) in sty.items():
            axx.plot(M.GRID, 20 * np.log10(abs(curves[n])), color=c, ls=ls, lw=1.4, label=f'model {n}')
        axx.set_xlim(*win[q]); axx.set_title(f'{q} channel'); axx.set_xlabel('f (GHz)'); axx.set_ylabel('|S21| (dB)')
        axx.grid(alpha=.3)
    names = list(res)
    for i, q in enumerate(NAMES):
        for j, (k, mk) in enumerate((('bare_fp_GHz', 'o'), ('bare_fr_GHz', 's'))):
            ax[2].plot(names, [res[n][q][k]['diff'] for n in names], marker=mk, ls='-' if i == 0 else ':',
                       color=('#1f6feb', '#d4731c')[j], label=f"{q} {'f_p' if j == 0 else 'f_r'}")
    ax[2].axhline(0, color='k', lw=.8); ax[2].set_ylabel('model - HFSS (MHz)'); ax[2].set_title('bare frequency error')
    ax[2].grid(alpha=.3); ax[2].legend(fontsize=8)
    tx = '\n'.join(f"{n}: J {res[n]['Q1']['J_MHz']['model']:.1f}/{res[n]['Q2']['J_MHz']['model']:.1f}" for n in names)
    ax[2].text(1.02, .5, 'J model (Q1/Q2), MHz\n' + tx + f"\nHFSS: {res['A']['Q1']['J_MHz']['hfss']:.1f}/"
               f"{res['A']['Q2']['J_MHz']['hfss']:.1f}", transform=ax[2].transAxes, fontsize=8, va='center')
    h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, loc='upper center', ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0, .93, .92)); fig.savefig(folder / 'compare.png', dpi=130); plt.close(fig)


if __name__ == '__main__':
    main()
