"""Q3D-based circuit model of the whole Design5 prr_r1 feed (C_in + six channels), compared with HFSS S21.

    python -I -B Design5/simulation/resonator/prr_circuit_model.py --hfss <prr_s21_reanalysis run>

Nothing is fitted to HFSS. Inputs: the Q3D loads and channel lengths of parameters/<rev>_build_checks.json (Q3D run
20261119_005709), the C_p / C_in pi network of the same Q3D run (series cp_mutual; pin a = frame, cp_frame_ground;
pin b = plate, cp_plate_ground), lossless 16/9 um CPW (purcell_readout_model preset), and the drawn positions:
  port 1 at x = -4.22 mm (LPPin tie) ... C_in pins a / b at -4.0325 / -3.7675 ... tees at chip6_prr.ORDER x ...
  port 2 at +4.22 mm; each tee: 150 um CPW to C_p pin a (frame side), C_p, plate side = filter open end.
Each channel: filter (open end -> tap 1.2 mm -> short), hammer arm (s/2 CPW, head excess to ground, C_j, head excess,
s/2 CPW), resonator (tap 1.3 mm -> open end with C_end, tap -> short). The model S21 is analysed with the same
prr_sim.s21_channel as HFSS. Cases: 'with_cin' (as drawn) and 'no_cin' (C_in replaced by a straight CPW).
"""
from pathlib import Path
import argparse
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator'),
                str(ROOT / 'pattern design' / 'resonator')]
import numpy as np  # noqa: E402
import prr_sim as S  # noqa: E402
import purcell_readout_model as PRM  # noqa: E402
from run_record import create_run, finish_run  # noqa: E402

DESIGN = 'Design5'
ORDER = {2: -3.20, 4: -1.90, 6: -0.60, 1: 0.70, 3: 2.00, 5: 3.30}   # chip6_prr.ORDER x (mm)
PORT_X, CIN_A, CIN_B, STUB_MM = 4.22, -4.0325, -3.7675, 0.150


def tl(f, L_mm, c):
    th = 2 * np.pi * f * L_mm * 1e-3 / c['v']
    z0 = c['z0']
    return np.array([[np.cos(th), 1j * z0 * np.sin(th)], [1j * np.sin(th) / z0, np.cos(th)]])


def shunt(y):
    return np.array([[1, 0], [y, 1]], dtype=complex)


def series(z):
    return np.array([[1, z], [0, 1]], dtype=complex)


def y_through(f, L_mm, y_load, c):
    """Input admittance of a line of L_mm terminated in y_load."""
    (a, b), (cc, d) = tl(f, L_mm, c)
    return (cc + d * y_load) / (a + b * y_load)


def y_channel(f, ch, lens, L, cfg, c):
    w = 2 * np.pi * f
    m = cfg['module']
    s2 = m['line_spacing_um'] / 2000
    Lf, Lr, tf, tr = lens['L_filter'], lens['L_res'], m['filter_tap_mm'], m['resonator_tap_mm']
    y_short = lambda length: -1j / (c['z0'] * np.tan(2 * np.pi * f * length * 1e-3 / c['v']))
    y_tr = y_short(Lr - tr) + y_through(f, tr, 1j * w * L['open_end'] * 1e-15, c)
    y_b = y_through(f, s2, y_tr, c) + 1j * w * L['head_excess'] * 1e-15
    y_a = 1 / (1 / y_b + 1 / (1j * w * L['cj'] * 1e-15)) + 1j * w * L['head_excess'] * 1e-15
    y_tf = y_through(f, s2, y_a, c) + y_short(Lf - tf)
    y_open_end = y_through(f, tf, y_tf, c)                       # filter seen from its open (C_p) end
    y_plate = y_open_end + 1j * w * L['cp_plate_ground'] * 1e-15
    y_frame = 1 / (1 / y_plate + 1 / (1j * w * L['cp_mutual'] * 1e-15)) + 1j * w * L['cp_frame_ground'] * 1e-15
    return y_through(f, STUB_MM, y_frame, c)                    # at the feed tee


def s21(f, lens, L, cfg, c, with_cin=True):
    out = np.empty(len(f), complex)
    w = 2 * np.pi * f
    for i, fi in enumerate(f):
        M = np.eye(2, dtype=complex)
        x = -PORT_X
        if with_cin:
            M = M @ tl(fi, CIN_A - x, c) @ shunt(1j * w[i] * L['cp_frame_ground'] * 1e-15) \
                @ series(1 / (1j * w[i] * L['cp_mutual'] * 1e-15)) @ shunt(1j * w[i] * L['cp_plate_ground'] * 1e-15)
            x = CIN_B
        for ch, xt in sorted(ORDER.items(), key=lambda kv: kv[1]):
            M = M @ tl(fi, xt - x, c) @ shunt(y_channel(fi, ch, lens[ch], L, cfg, c))
            x = xt
        M = M @ tl(fi, PORT_X - x, c)
        (a, b), (cc, d) = M
        z0 = c['z0']
        out[i] = 2 / (a + b / z0 + cc * z0 + d)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--hfss', type=Path, required=True, help='prr_s21_reanalysis run folder')
    ap.add_argument('--rev', default='2026-11-19_prr_r1')
    a = ap.parse_args()
    hf = a.hfss.resolve()
    checks_file = ROOT / DESIGN / 'parameters' / f'{a.rev}_build_checks.json'
    checks = json.loads(checks_file.read_text(encoding='utf-8'))
    q3d_run = ROOT / checks['loads_source']
    q3d = json.loads((q3d_run / 'summary.json').read_text(encoding='utf-8'))['derived']
    cfg = json.loads((ROOT / DESIGN / 'parameters' / 'readout_prr.json').read_text(encoding='utf-8'))
    L = dict(checks['loads'], cp_frame_ground=q3d['cp_frame_ground'])
    lens = {int(k): v for k, v in checks['channels'].items()}
    c = PRM.line(cfg['cpw']['model_preset'])
    folder = create_run(DESIGN, 'resonator', 'prr_circuit_model', entrypoint=__file__,
                        inputs=[checks_file, q3d_run / 'summary.json', hf / 'summary.json'],
                        upstream=[str(q3d_run), str(hf)],
                        parameters=dict(loads_fF=L, cpw=cfg['cpw']['model_preset'], positions_mm=dict(
                            ports=[-PORT_X, PORT_X], cin_pins=[CIN_A, CIN_B], tees=ORDER, tee_to_cp_mm=STUB_MM),
                            model='lossless TL + Q3D lumped C; no HFSS input except for the comparison'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        f = np.arange(5.75, 7.25, 0.00005) * 1e9
        hs = json.loads((hf / 'summary.json').read_text(encoding='utf-8'))
        pred = hs['prediction']
        names = sorted(pred, key=lambda k: pred[k]['f_r_GHz'])
        fr = [pred[q]['f_r_GHz'] for q in names]
        edges = [5.75] + [(x + y) / 2 for x, y in zip(fr, fr[1:])] + [7.25]
        cases, curves = {}, {}
        for case, cin in (('with_cin', True), ('no_cin', False)):
            s = s21(f, lens, L, cfg, c, cin)
            np.savetxt(folder / f'model_S21_{case}.csv', np.c_[f / 1e9, s.real, s.imag], delimiter=',',
                       header='f_GHz,re_S21,im_S21', comments='')
            curves[case] = s
            cases[case] = {}
            for i, q in enumerate(names):
                r = S.s21_channel(f / 1e9, s, (edges[i], edges[i + 1]))
                cases[case][q] = dict(r['poles'], zeros=r['zeros'])
        rows = {}
        for q in names:
            h = hs['hfss'][q]['poles']
            rows[q] = dict(target_GHz=pred[q]['f_r_GHz'], kappa_hanger_MHz=pred[q]['kappa_p_MHz'],
                           J_closed_form_MHz=pred[q]['J_hammer_MHz'],
                           kappa_model_cin_MHz=cases['with_cin'][q]['kappa_p_MHz'],
                           kappa_model_nocin_MHz=cases['no_cin'][q]['kappa_p_MHz'], kappa_hfss_MHz=h['kappa_p_MHz'],
                           J_model_cin_MHz=cases['with_cin'][q]['J_MHz'], J_hfss_MHz=h['J_MHz'],
                           f_p_model_cin_GHz=cases['with_cin'][q]['f_p_GHz'], f_r_model_cin_GHz=cases['with_cin'][q]['f_r_GHz'],
                           f_p_hfss_GHz=h['f_p_GHz'], f_r_hfss_GHz=h['f_r_GHz'])
        fig, ax = plt.subplots(2, 1, figsize=(13, 8))
        fh, shf = S.load_s(next((ROOT / hs_up for hs_up in [Path(json.loads((hf / 'run.json').read_text(
            encoding='utf-8'))['upstream'][0])]), None) / 'raw' / 'S_wide.csv')
        ax[0].plot(fh, 20 * np.log10(abs(shf)), 'k-', lw=1.4, label='HFSS 1 pass')
        ax[0].plot(f / 1e9, 20 * np.log10(abs(curves['with_cin'])), 'C3-', lw=.9, label='circuit model, with C_in')
        ax[0].plot(f / 1e9, 20 * np.log10(abs(curves['no_cin'])), 'C0-', lw=.7, alpha=.7, label='circuit model, no C_in')
        ax[0].set_ylabel('|S21| (dB)'); ax[0].set_xlabel('Frequency (GHz)'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
        ax[0].set_title('Design5 prr_r1: Q3D circuit model (no HFSS input) vs HFSS S21', fontsize=9)
        xs = np.arange(len(names))
        for k, (key, lab, mk) in enumerate((('kappa_hanger_MHz', 'hanger formula (matched feed)', 's'),
                                            ('kappa_model_nocin_MHz', 'circuit, no C_in', 'v'),
                                            ('kappa_model_cin_MHz', 'circuit, with C_in', 'o'),
                                            ('kappa_hfss_MHz', 'HFSS S21 fit (1 pass)', 'D'))):
            ax[1].plot(xs + (k - 1.5) * .08, [rows[q][key] for q in names], mk, ms=7, label=lab)
        ax[1].set_xticks(xs, [f'{q}\n{rows[q]["target_GHz"]:.1f} GHz' for q in names])
        ax[1].set_ylabel('kappa_p / 2pi (MHz)'); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
        fig.tight_layout(); fig.savefig(folder / 'circuit_vs_hfss.png', dpi=140); plt.close(fig)
        for q in names:
            v = rows[q]
            print(f"{q} {v['target_GHz']:.1f}: kappa hanger {v['kappa_hanger_MHz']:.1f} | no-C_in {v['kappa_model_nocin_MHz']:.1f}"
                  f" | with-C_in {v['kappa_model_cin_MHz']:.1f} | HFSS {v['kappa_hfss_MHz']:.1f} MHz ;"
                  f" J model {v['J_model_cin_MHz']:.2f} HFSS {v['J_hfss_MHz']:.2f} ;"
                  f" f_p/f_r model {v['f_p_model_cin_GHz']:.4f}/{v['f_r_model_cin_GHz']:.4f}"
                  f" HFSS {v['f_p_hfss_GHz']:.4f}/{v['f_r_hfss_GHz']:.4f}", flush=True)
        summary = dict(comparison=rows, model=cases, validation='model only; HFSS comparison is a 1-pass run')
        finish_run(folder, status='completed', summary=json.loads(json.dumps(summary, default=float)))
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
