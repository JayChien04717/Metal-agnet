"""Design4: Q3D circuit-model S21 against HFSS driven S21 for each stub variant (no EM rerun, nothing fitted).

Circuit (lossless, 16/9 um CPW preset of purcell_readout_model, e^{+jwt}): port 1 - feed line - tee n0 - port 2, with
the channel branch at the tee: tee -> C_p pin a line (drawn length), C_p pi network (frame ground at pin a, series
cp_mutual, plate ground at pin b), filter open end -> tap (1.2 mm) -> short arm; tap -> s/2 line -> head excess,
series C_j, head excess -> s/2 line -> resonator tap; resonator short arm, open arm 1.3 mm -> OpenToGround end C.
All loads from the variant's Q3D run, lengths as drawn (build_info.json of the driven run).
Both S21 traces go through the same extraction (purcell_resonator_pair_sim.s21_channel: S21 zeros + S21 = 1 point,
and a pole fit with the zeros fixed -> bare f_p, f_r, J, kappa_p). The parallel-run bound uses
purcell_readout_model.parallel_j on PurcellReadout.parallel_runs().

    python -I -B Design4/simulation/resonator/compare.py <driven run> [<driven run> ...]
"""
from pathlib import Path
import argparse
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator'), str(ROOT / 'Design4' / 'layout')]
import project_paths  # noqa: E402
project_paths.configure(__file__)
import numpy as np  # noqa: E402

import build as B  # noqa: E402
import purcell_readout_model as M  # noqa: E402
import purcell_resonator_pair_sim as S  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402


def abcd_tl(f, L_mm, c):
    th = 2 * np.pi * f * 1e9 * L_mm * 1e-3 / c['v']
    z0 = c['z0']
    return np.array([[np.cos(th), 1j * z0 * np.sin(th)], [1j * np.sin(th) / z0, np.cos(th)]])


def y_in(f, L_mm, y_load, c):
    (a, b), (cc, d) = abcd_tl(f, L_mm, c)
    return (cc + d * y_load) / (a + b * y_load)


def y_branch(f, lengths, loads, tee_to_pin_a, c):
    """Admittance (S) of the whole channel seen from the feed tee."""
    w = 2 * np.pi * f * 1e9
    g = B.GEOMETRY['module']
    s2 = lengths['stub_mm']
    Lf, Lr, tf, tr = lengths['filter_mm'], lengths['resonator_mm'], g['filter_tap_mm'], g['resonator_tap_mm']
    yc = lambda c_fF: 1j * w * c_fF * 1e-15
    y_short = lambda length: -1j / (c['z0'] * np.tan(2 * np.pi * f * 1e9 * length * 1e-3 / c['v']))
    y_tr = y_short(Lr - tr) + y_in(f, tr, yc(lengths['c_end_resonator_fF']), c)
    y_b = y_in(f, s2, y_tr, c) + yc(loads['head_excess'])
    y_a = 1 / (1 / y_b + 1 / yc(loads['cj'])) + yc(loads['head_excess'])
    y_tf = y_in(f, s2, y_a, c) + y_short(Lf - tf)
    y_plate = y_in(f, tf, y_tf, c) + yc(loads['cp_plate_ground'])
    y_frame = 1 / (1 / y_plate + 1 / yc(loads['cp_mutual'])) + yc(loads['cp_frame_ground'])
    return y_in(f, tee_to_pin_a, y_frame, c)


def s21_model(fs, lengths, loads, sections, c):
    (x1, x2), xt = sections['port_x'], sections['tee_x']
    out = np.empty(len(fs), complex)
    for i, f in enumerate(fs):
        m = abcd_tl(f, xt - x1, c) @ np.array([[1, 0], [y_branch(f, lengths, loads, sections['tee_to_cp_pin_a'], c), 1]]) \
            @ abcd_tl(f, x2 - xt, c)
        (a, b), (cc, d) = m
        out[i] = 2 / (a + b / c['z0'] + cc * c['z0'] + d)
    return out


def closed_form(lengths, loads):
    g = B.GEOMETRY['module']
    cpw = lengths['cpw']
    fp = M.tapped_lambda4_f(lengths['filter_mm'], g['filter_tap_mm'], lengths['c_end_filter_fF'], lengths['stub_mm'],
                            lengths['c_stub_fF'], cpw)
    fr = M.tapped_lambda4_f(lengths['resonator_mm'], g['resonator_tap_mm'], lengths['c_end_resonator_fF'],
                            lengths['stub_mm'], lengths['c_stub_fF'], cpw)
    j = M.j_hammer((fp + fr) / 2, loads['cj'], lengths['filter_mm'], g['filter_tap_mm'], lengths['c_end_filter_fF'],
                   lengths['resonator_mm'], g['resonator_tap_mm'], lengths['c_end_resonator_fF'], cpw)
    return dict(f_p_GHz=fp, f_r_GHz=fr, J_hammer_MHz=j, kappa_p_hanger_MHz=M.kappa_hanger(fp, loads['cp_mutual'], cpw))


def strip(r):
    return {k: v for k, v in r.items() if k not in ('f', 's', 'fit')}


def analyse(run):
    p = json.loads((run / 'run.json').read_text(encoding='utf-8'))['parameters']
    info = json.loads((run / 'build_info.json').read_text(encoding='utf-8'))
    raw = json.loads((run / 'raw' / 'driven.json').read_text(encoding='utf-8'))
    v, lengths, loads = p['variant'], p['lengths'], p['q3d_loads']
    c = M.line(lengths['cpw'])
    f, s_h = S.load_s(run / 'raw' / raw['s_file'])
    s_m = s21_model(f, lengths, loads, info['circuit_sections_mm'], c)
    sign = min((1, -1), key=lambda k: np.mean(abs(k * s_m - s_h) ** 2))      # port polarity convention
    s_m = sign * s_m
    window = (float(f.min()), float(f.max()))
    hf, mo = S.s21_channel(f, s_h, window), S.s21_channel(f, s_m, window)
    d, _ = B.build(v, lengths=lengths, wire_bonds=False)
    runs = d.components['purcell_resonator_pair'].parallel_runs()
    fmid = (mo['poles']['f_p_GHz'] + mo['poles']['f_r_GHz']) / 2
    par = M.parallel_j(fmid, runs, lengths['filter_mm'], lengths['c_end_filter_fF'], lengths['resonator_mm'],
                       lengths['c_end_resonator_fF'], lengths['cpw']) if runs else dict(J_max=0.0)
    hp, mp = hf['poles'], mo['poles']
    diff = dict(df_p_MHz=(hp['f_p_GHz'] - mp['f_p_GHz']) * 1e3, df_r_MHz=(hp['f_r_GHz'] - mp['f_r_GHz']) * 1e3,
                J_ratio_hfss_over_model=hp['J_MHz'] / mp['J_MHz'], kappa_ratio_hfss_over_model=hp['kappa_p_MHz'] / mp['kappa_p_MHz'],
                abs_s21_rms=float(np.sqrt(np.mean((abs(s_h) - abs(s_m)) ** 2))),
                complex_s21_rms=float(np.sqrt(np.mean(abs(s_h - s_m) ** 2))), model_sign=sign)
    return dict(variant=v, run=str(run.relative_to(ROOT)), line_spacing_um=p['line_spacing_um'],
                stub_metal_um=p['stub_metal_um'], q3d_loads=loads, lengths_mm=[lengths['filter_mm'], lengths['resonator_mm']],
                hfss_passes=raw['passes'], hfss_last_delta_s=raw['last_delta_s'], closed_form=closed_form(lengths, loads),
                hfss=strip(hf), model=strip(mo), difference=diff, parallel_runs=runs, parallel_J=par), (f, s_h, s_m, hf, mo)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('driven', type=Path, nargs='+')
    args = ap.parse_args()
    runs = [r.resolve() for r in args.driven]
    folder = create_run('Design4', 'resonator', 'compare', entrypoint=__file__,
                        inputs=[r / 'run.json' for r in runs], upstream=[str(r) for r in runs],
                        parameters=dict(driven_runs=[str(r) for r in runs], fitted='nothing (pole extraction only)'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        out, traces = [], []
        for r in runs:
            res, tr = analyse(r)
            out.append(res)
            traces.append(tr)
        fig, axs = plt.subplots(2, len(runs), figsize=(6.5 * len(runs), 7), squeeze=False,
                                gridspec_kw=dict(height_ratios=[3, 1.3]))
        for k, (res, (f, s_h, s_m, hf, mo)) in enumerate(zip(out, traces)):
            ax, ax2 = axs[0, k], axs[1, k]
            ax.plot(f, 20 * np.log10(abs(s_h)), 'k-', lw=1.4, label='HFSS Driven')
            ax.plot(f, 20 * np.log10(abs(s_m)), color='tab:blue', ls='--', lw=1.4, label='Q3D circuit model')
            hp, mp = res['hfss']['poles'], res['model']['poles']
            ax.set_title(f"{res['variant']}: stub {res['stub_metal_um']:g} um, spacing {res['line_spacing_um']:g} um\n"
                         f"J HFSS {hp['J_MHz']:.2f} / model {mp['J_MHz']:.2f} MHz; "
                         f"f_p {res['difference']['df_p_MHz']:+.1f}, f_r {res['difference']['df_r_MHz']:+.1f} MHz", fontsize=9)
            ax.set_ylabel('|S21| (dB)'); ax.grid(alpha=.3); ax.legend(fontsize=8)
            ax2.plot(f, abs(s_h - s_m), color='tab:red', lw=1)
            ax2.set_ylabel('|S21 HFSS - model|'); ax2.set_xlabel('Frequency (GHz)'); ax2.grid(alpha=.3)
        fig.tight_layout(); fig.savefig(folder / 'model_vs_driven.png', dpi=150); plt.close(fig)
        summary = dict(variants=out, figure='model_vs_driven.png',
                       note='black = HFSS Driven, blue dashed = circuit fed with Q3D capacitances (Q3D does not give S21 '
                            'directly); bare parameters from the same two-pole extraction for both traces')
        for res in out:
            print(res['variant'], json.dumps(dict(hfss=res['hfss']['poles'], model=res['model']['poles'],
                                                  difference=res['difference'], parallel_J=res['parallel_J']), default=float))
        finish_run(folder, status='completed', summary=summary)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
