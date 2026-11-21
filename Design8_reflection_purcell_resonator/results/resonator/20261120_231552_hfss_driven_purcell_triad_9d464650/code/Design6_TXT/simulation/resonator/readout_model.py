"""Design6_TXT readout-only circuit model with the new Q3D Cp and hammer values; length re-tuning.

Replaces purcell_6500.tune_lengths() for Design6_TXT (that function uses the Design5_6_purcell-res C_in network and an interpolated
~2 fF hammer). Nodal TL network of the actual r9 readout: both channels on one feed line (no C_in), 50 ohm ports,
enclosing Cp as a pi network, deembedded hammer (mutual + excess grounds) at the centre plane with spacing/2 arms,
filter / resonator as uncoupled lambda/4 lines from the taps (as Design4_Single_purcell-res hammer_fixed_s21.model, 'deembedded').
No claw: the readout open end is a bare open, matching readout_s21.py.

Complex poles are found by Newton iteration on det(I + 50 Y_port(f)); each channel's pole pair is inverted to bare
fp / fr / J with the same assumptions as readout_poles.py. Outputs:
  1. prediction for the r9 lengths with --res-extra-mm (default 0.1 = CLAW_DL, i.e. the running HFSS geometry);
  2. lengths giving bare fp = fr = f_target * HFSS_CORR per channel (the tune_lengths convention).
Line parameters: quasi-static conformal CPW, eps_eff = (1 + eps_r)/2, no fitting. Parallel-line coupling is omitted.
"""
from argparse import ArgumentParser
from pathlib import Path
import json
import sys
import numpy as np
from scipy.special import ellipk

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'Design6_TXT/simulation/resonator')]
import readout_s21 as R           # noqa: E402  (sets sys.path for layout / run_record)
from readout_poles import invert, fit  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

L, P = R.L, R.L.P
CAPS_RUN = ROOT / 'Design6_TXT/results/resonator/20261116_221901_readout_caps_0e896b6c'
NAMES = ('Q1', 'Q2')


def geometry(d):
    pl = lambda n: L.d5.f3.path_length(d, n)
    th = P.TEE_HALF
    feeds = dict(feed_in=pl('feed_in') + th, feed_mid=pl('feed_mid') + 2 * th, feed_out=pl('feed_out') + th,
                 **{f'stub_p_{q}': pl(f'stub_p_{q}') + .1 for q in NAMES})
    spacing = {q: abs(d.parse_value(d.components[f'jf_{q}'].options.pos_x) -
                      d.parse_value(d.components[f'jr_{q}'].options.pos_x)) for q in NAMES}
    w, g = (d.parse_value(L.CPW_RO[k]) for k in ('trace_width', 'trace_gap'))
    eps = (1 + L.d5.SI_EPS) / 2; k = w / (w + 2 * g)
    z0 = 30 * np.pi / np.sqrt(eps) * ellipk(1 - k * k) / ellipk(k * k)
    return dict(feeds_mm=feeds, spacing_mm=spacing, z0_ohm=z0, velocity_m_s=299792458 / np.sqrt(eps))


def port_Y(f, lengths, geo, cp, hammer):
    """2x2 port admittance (S) at complex frequency f (GHz)."""
    w = 2 * np.pi * np.asarray(f, complex) * 1e9
    z0, v = geo['z0_ohm'], geo['velocity_m_s']
    # 0 in, 1 out, 2 n0_Q1, 3 n0_Q2, then per channel 4 nodes: cp_a, cp_b(=filter start), filter tap, res tap, hl, hr
    n = 4 + 6 * 2
    Y = np.zeros((n, n), complex)

    def tl(a, b, mm):
        t = w / v * mm * 1e-3; dg, off = -1j / z0 / np.tan(t), 1j / z0 / np.sin(t)
        Y[a, a] += dg; Y[b, b] += dg; Y[a, b] += off; Y[b, a] += off

    def pi(a, b, c):
        cc, ga, gb = np.asarray(c) * 1e-15
        Y[a, a] += 1j * w * (cc + ga); Y[b, b] += 1j * w * (cc + gb); Y[a, b] -= 1j * w * cc; Y[b, a] -= 1j * w * cc

    def stub(a, mm, shorted):
        t = w / v * mm * 1e-3
        Y[a, a] += -1j / z0 / np.tan(t) if shorted else 1j / z0 * np.tan(t)

    F = geo['feeds_mm']
    tl(0, 2, F['feed_in']); tl(2, 3, F['feed_mid']); tl(3, 1, F['feed_out'])
    for i, q in enumerate(NAMES):
        a, b, tf, tr, hl, hr = range(4 + 6 * i, 10 + 6 * i)
        Lf, Lr = lengths[q]
        tl(2 + i, a, F[f'stub_p_{q}']); pi(a, b, cp)
        tl(b, tf, P.FILTER_TAP); stub(tf, Lf - P.FILTER_TAP, True)
        stub(tr, P.RES_TAP, False); stub(tr, Lr - P.RES_TAP, True)
        s = geo['spacing_mm'][q] / 2
        tl(tf, hl, s); pi(hl, hr, hammer); tl(hr, tr, s)
    Yp = Y[:2, :2] - Y[:2, 2:] @ np.linalg.solve(Y[2:, 2:], Y[2:, :2])
    return Yp


def s21(f, *args):
    out = []
    for x in f:
        Yp = port_Y(x, *args); S = np.linalg.solve(np.eye(2) + 50 * Yp, np.eye(2) - 50 * Yp); out.append(S[1, 0])
    return np.array(out)


GRID = np.arange(6.6, 7.5, .00002)   # 20 kHz, same extraction as readout_poles.py on HFSS


def channel_poles(guesses, *args):
    """guesses: {'Q2': (purcell, readout), 'Q1': (...)} complex GHz, Im = half-width (GHz).
    Four-pole rational fit of the model S21 (readout_poles.fit); returns {q: (broad, narrow)}."""
    order = [('Q2', 0), ('Q2', 1), ('Q1', 0), ('Q1', 1)]
    g = [(guesses[q][i].real, max(guesses[q][i].imag * 1000, 1e-3)) for q, i in order]
    poles, _, rms = fit(GRID, s21(GRID, *args), g, (GRID[0], GRID[-1]))
    assert rms < 1e-3, f'model pole fit rms {rms}'
    out = {}
    for k, q in ((0, 'Q2'), (2, 'Q1')):
        a, b = poles[k], poles[k + 1]
        out[q] = (a, b) if a.imag >= b.imag else (b, a)
    return out


def bare(poles):
    return {q: invert(*pq) for q, pq in poles.items()}


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--res-extra-mm', type=float, default=L.CLAW_DL)
    args = ap.parse_args()
    caps = json.loads((CAPS_RUN / 'summary_out.json').read_text())
    assert caps['cp']['signals'] == ['frame_cp', 'lead_b_cp']      # frame on the n0 side, lead_b on the filter side
    cp = [caps['cp']['mutual_fF'], *caps['cp']['ground_fF']]
    hammer = [caps['hammer']['mutual_fF'], *caps['hammer_excess_ground_fF']]
    d, _, lengths = R.build(args.res_extra_mm)
    geo = geometry(d)
    cur = {q: (lengths[q]['L_filter'], lengths[q]['L_res']) for q in NAMES}
    folder = create_run('Design6_TXT', 'resonator', 'readout_model', entrypoint=__file__,
                        inputs=[CAPS_RUN / 'run.json', CAPS_RUN / 'summary_out.json', Path(L._length_file())],
                        upstream=[str(CAPS_RUN)],
                        parameters=dict(EM_rerun=False, res_extra_mm=args.res_extra_mm, cp_fF=cp, hammer_deembedded_fF=hammer,
                                        cp_q3d_converged=caps['cp']['solver_converged'], geometry=geo, lengths_mm=cur,
                                        HFSS_CORR=L.HFSS_CORR, targets_GHz={q: L.READOUT[q]['f'] for q in NAMES},
                                        omitted=['parallel-line C_m', 'claw (bare open end)', 'Cp / tee endpoint deembedding']))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        model = lambda lens: (lens, geo, cp, hammer)
        from scipy.signal import find_peaks, peak_widths

        def poles_of(lens):
            """Fresh dip detection on the model S21, then the four-pole fit (guesses never carried over)."""
            db = 20 * np.log10(abs(s21(GRID, *model(lens))))
            idx, _ = find_peaks(-db, prominence=.3)
            width = peak_widths(-db, idx, rel_height=.5)[0] * .00002
            assert len(idx) == 4, f'dips {GRID[idx]}, widths {width}'
            gu = {}      # Q2 lines are longer -> the lower pair; the wider dip of a pair is the Purcell-like mode
            for q, k in (('Q2', [0, 1]), ('Q1', [2, 3])):
                b, n = sorted(k, key=lambda i: -width[i])
                gu[q] = tuple(complex(*x) for x in sorted([(GRID[idx][b], .010), (GRID[idx][n], .0002)]))
            return channel_poles(gu, *model(lens))

        poles0 = poles_of(cur)
        pred = bare(poles0)
        targets = {q: L.READOUT[q]['f'] * L.HFSS_CORR for q in NAMES}

        def iterate(update, lens, tol=1e-5, steps=20):
            """Fixed point with f ~ 1/L per line; update(q, lens, bare) -> new (L_filter, L_res)."""
            hist = []
            for _ in range(steps):
                b = bare(poles_of(lens)); hist.append({q: dict(L=lens[q], **b[q]) for q in NAMES})
                new = {q: update(q, lens[q], b[q]) for q in NAMES}
                if max(abs(np.subtract(new[q], lens[q])).max() for q in NAMES) < tol:
                    return lens, b, hist, True
                lens = new
            return lens, b, hist, False

        tuned_lengths, tuned, hist_t, ok_t = iterate(
            lambda q, l, b: (l[0] * b['bare_fp_GHz'] / targets[q], l[1] * b['bare_fr_GHz'] / targets[q]), dict(cur))
        al_lengths, al, hist_a, ok_a = iterate(
            lambda q, l, b: (l[0], l[1] * b['bare_fr_GHz'] / b['bare_fp_GHz']), dict(cur))
        same_freq = {q: dict(L_res_mm=al_lengths[q][1], delta_L_res_mm=al_lengths[q][1] - cur[q][1], **al[q])
                     for q in NAMES}
        same_freq['converged'] = ok_a
        out = dict(
            prediction_current=dict(lengths_mm=cur, res_extra_mm=args.res_extra_mm,
                                    poles=[[q, [[z.real, z.imag * 1000] for z in pq]] for q, pq in poles0.items()],
                                    bare=pred),
            align_resonator_only=same_freq,
            retune_to_targets=dict(model_targets_GHz=targets, lengths_mm=tuned_lengths, bare=tuned,
                                   delta_mm={q: [tuned_lengths[q][0] - cur[q][0], tuned_lengths[q][1] - cur[q][1]] for q in NAMES},
                                   converged=ok_t, iterations=len(hist_t),
                                   max_residual_GHz=max(abs(tuned[q][k] - targets[q]) for q in NAMES for k in ('bare_fp_GHz', 'bare_fr_GHz')),
                                   note='Model frequency = f_target * HFSS_CORR (Design4_Single_purcell-res 6.08/5.975); the correction is not '
                                        're-derived for this model. For the fabricated chip subtract CLAW_DL from L_res again.'),
            limitations=['Cp Q3D not converged (0.14 %)', 'parallel-line coupling omitted (Design4_Single_purcell-res: ~6 MHz of J)',
                         'Design4_Single_purcell-res analogue model was 20-35 MHz low in fp / fr vs HFSS', 'no claw'])
        (folder / 'model.json').write_text(json.dumps(out, indent=1, default=float))
        finish_run(folder, status='completed', summary=out)
        print(json.dumps(out, indent=1, default=float))
    except Exception:
        import traceback
        (folder / 'error.log').write_text(traceback.format_exc())
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()}); raise


if __name__ == '__main__':
    main()
