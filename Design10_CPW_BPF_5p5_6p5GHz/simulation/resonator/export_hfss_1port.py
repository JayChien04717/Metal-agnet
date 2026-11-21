# -*- coding: utf-8 -*-
"""Re-analyze raw HFSS 50-Ohm S-parameter CSV for Design10 (1-port 3-stage reflection BPF),
exporting normalized open-end probe transfer, unloaded 1-port S11 (.csv + .s1p), and a pure
HFSS-only 3-panel plot (no analytic curve overlay).
"""

import io
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general')]
import project_paths  # noqa: E402
from run_record import _long_path, create_run, finish_run  # noqa: E402

project_paths.configure(__file__)
DESIGN = project_paths.design_dir(__file__)
Z_PROBE_OHM = 50000.0


def parse_complex_col(series):
    if np.iscomplexobj(series):
        return series.to_numpy(complex)
    return series.map(lambda v: complex(str(v).replace('i', 'j').strip())).to_numpy(complex)


def terminate_open_probe(s_df, z_probe_ohm=Z_PROBE_OHM, z0_ohm=50.0):
    s11 = parse_complex_col(s_df['S11'])
    s21 = parse_complex_col(s_df['S21'])
    s22 = parse_complex_col(s_df['S22'])
    g2 = (z_probe_ohm - z0_ohm) / (z_probe_ohm + z0_ohm)
    denom = 1.0 - g2 * s22
    s11_open = s11 + (g2 * (s21 ** 2)) / denom
    s21_probe = (math.sqrt(1.0 - g2 ** 2) * s21) / denom
    return s11_open, s21_probe


def main():
    parent_name = sys.argv[1] if len(sys.argv) > 1 else '20261121_162252_hfss_driven_1port_3stage_9d5ee50f'
    parent = DESIGN / 'results' / 'resonator' / parent_name
    s_df = pd.read_csv(_long_path(parent / 'raw' / 'S_params_raw50.csv'), index_col=0)
    conv = pd.read_csv(_long_path(parent / 'raw' / 'hfss_convergence.csv'))

    re_folder = create_run(
        DESIGN.name,
        'resonator',
        'reanalysis_hfss_1port_s21_s11',
        entrypoint=__file__,
        inputs=[parent / 'raw' / 'S_params_raw50.csv', parent / 'raw' / 'hfss_convergence.csv'],
        upstream=[str(parent.relative_to(ROOT).as_posix())],
        parameters=dict(parent_run=parent.name, z_probe_ohm=Z_PROBE_OHM, geometry_revision='2026-11-21_r2'),
    )
    rec = json.loads(_long_path(re_folder / 'run.json').read_text(encoding='utf-8'))
    rec['geometry_revision'] = '2026-11-21_r2'
    _long_path(re_folder / 'run.json').write_text(json.dumps(rec, indent=2) + '\n', encoding='utf-8')
    print('REANALYSIS_FOLDER=' + re_folder.name, flush=True)

    freqs = s_df.index.to_numpy(float)
    s11_open, s21_probe = terminate_open_probe(s_df)
    s21_db = 20.0 * np.log10(np.maximum(np.abs(s21_probe), 1e-12))
    s21_norm = s21_db - float(np.max(s21_db))
    s11_db = 20.0 * np.log10(np.maximum(np.abs(s11_open), 1e-12))
    ph_rad = np.unwrap(np.angle(s11_open))
    tau_ns = np.abs(np.gradient(ph_rad) / (2.0 * np.pi * np.gradient(freqs) * 1e9)) * 1e9

    _long_path(re_folder / 'summaries').mkdir(parents=True, exist_ok=True)
    out_tbl = pd.DataFrame(dict(
        freq_GHz=freqs,
        S21_probe_norm_dB=np.round(s21_norm, 4),
        S21_probe_dB=np.round(s21_db, 4),
        S11_dB=np.round(s11_db, 4),
        S11_deg=np.round(np.degrees(np.angle(s11_open)), 3),
        S11_unwrapped_deg=np.round(np.degrees(ph_rad), 3),
        tau_g_ns=np.round(tau_ns, 5),
    ))
    out_tbl.to_csv(_long_path(re_folder / 'summaries' / 'hfss_1port_2_7p5GHz.csv'), index=False)

    s1p_lines = [
        f'! {DESIGN.name} HFSS Driven unloaded 1-port S11 (2026-11-21_r2)',
        f'! Parent HFSS run: {parent.name} (3 passes, 10 um mesh, 451272 tets, Delta S = 0.033366)',
        '# GHz S RI R 50',
    ]
    for f, s11v in zip(freqs, s11_open):
        s1p_lines.append(f'{f:.6f} {s11v.real:.8e} {s11v.imag:.8e}')
    _long_path(re_folder / 'summaries' / 'hfss_1port_s11.s1p').write_text('\n'.join(s1p_lines) + '\n', encoding='utf-8')

    pb = (freqs >= 5.5) & (freqs <= 6.5)
    sb = (freqs >= 2.0) & (freqs <= 4.0)
    idx3 = np.where(s21_norm >= -3.0)[0]
    f3_lo, f3_hi = float(freqs[idx3[0]]), float(freqs[idx3[-1]])
    last_ds = float(conv['Max Mag. Delta S'].dropna().iloc[-1])
    tets = int(conv['Solved Elements'].dropna().iloc[-1])

    def at_f(ft):
        return float(s21_norm[np.argmin(np.abs(freqs - ft))])

    pb_min = float(np.min(s21_norm[pb]))
    sb_max = float(np.max(s21_norm[sb]))

    metrics = dict(
        band_3dB_GHz=[round(f3_lo, 4), round(f3_hi, 4)],
        f0_center_3dB_GHz=round(0.5 * (f3_lo + f3_hi), 4),
        pb_5p5_6p5_min_norm_dB=round(pb_min, 3),
        sb_2_4_max_norm_dB=round(sb_max, 2),
        sb_2_4_meets_40dB=bool(sb_max <= -40.0),
        norm_at_2GHz_dB=round(at_f(2.0), 2),
        norm_at_3GHz_dB=round(at_f(3.0), 2),
        norm_at_4GHz_dB=round(at_f(4.0), 2),
        norm_at_5GHz_dB=round(at_f(5.0), 2),
        norm_at_6GHz_dB=round(at_f(6.0), 2),
        norm_at_7GHz_dB=round(at_f(7.0), 2),
        tau_g_5p5_6p5_median_ns=round(float(np.median(tau_ns[pb])), 3),
        tau_g_3GHz_ns=round(float(tau_ns[np.argmin(np.abs(freqs - 3.0))]), 4),
        passes=int(len(conv)),
        solved_elements=tets,
        last_delta_s=round(last_ds, 6),
        mesh_length='10um',
        strap_count=40,
        aedt_design='D10_1P_9d5ee50f_hfss',
        parent_run=parent.name,
    )

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(9.2, 9.2), sharex=True)
    ax0.axvspan(2.0, 4.0, color='#fadbd8', alpha=0.45, label='Stopband target: 2–4 GHz (≤ -40 dB)')
    ax0.axvspan(5.5, 6.5, color='#d6eaf8', alpha=0.50, label='Passband target: 5.5–6.5 GHz')
    ax0.plot(freqs, s21_norm, color='#1b4f72', lw=2.0, label='HFSS normalized |S21_probe| (50 kΩ open-end probe)')
    ax0.axhline(-40.0, color='#922b21', lw=0.9, ls=':')
    ax0.axhline(-3.0, color='#555555', lw=0.8, ls=':')
    ax0.set_ylim(-75, 2)
    ax0.set_ylabel('Normalized |S21_probe| (dB)')
    ax0.set_title(
        f'{DESIGN.name} · 1-Port 3-Stage Reflection BPF · HFSS Driven Only '
        f'(3 passes, 10 μm, ΔS={last_ds:.4f}, {tets:,} tets)\n'
        f'-3 dB: {f3_lo:.2f}–{f3_hi:.2f} GHz | 5.5–6.5 GHz min={pb_min:.2f} dB | '
        f'2–4 GHz max={sb_max:.1f} dB',
        fontsize=8.6,
        weight='bold',
    )
    ax0.grid(True, alpha=0.3)
    ax0.legend(loc='lower right', fontsize=8.3, frameon=True)

    ax1.axvspan(2.0, 4.0, color='#fadbd8', alpha=0.35)
    ax1.axvspan(5.5, 6.5, color='#d6eaf8', alpha=0.45)
    ax1.plot(freqs, s11_db, color='#c0392b', lw=1.8, label='HFSS |S11| (dB, left)')
    ax1.set_ylim(-3.0, 0.5)
    ax1.set_ylabel('HFSS |S11| (dB)', color='#c0392b')
    ax1.tick_params(axis='y', labelcolor='#c0392b')
    ax1.grid(True, alpha=0.3)
    ax1r = ax1.twinx()
    ax1r.plot(freqs, np.degrees(ph_rad), color='#7d3c98', lw=1.5, ls='--', label='HFSS ∠S11 unwrapped (deg, right)')
    ax1r.set_ylabel('∠S11 (deg)', color='#7d3c98')
    ax1r.tick_params(axis='y', labelcolor='#7d3c98')
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax1r.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc='lower left', fontsize=8.3, frameon=True)
    ax1.set_title('1-Port Reflection S11: Magnitude & Unwrapped Phase (HFSS Driven Only)', fontsize=9.2, weight='bold')

    ax2.axvspan(2.0, 4.0, color='#fadbd8', alpha=0.35)
    ax2.axvspan(5.5, 6.5, color='#d6eaf8', alpha=0.45)
    ax2.plot(freqs, tau_ns, color='#1e8449', lw=1.8, label='HFSS 1-port reflection group delay τ_g(S11)')
    ax2.set_yscale('log')
    ax2.set_xlim(freqs.min(), freqs.max())
    ax2.set_xlabel('Frequency (GHz)')
    ax2.set_ylabel('Group Delay τ_g(S11) (ns)')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend(loc='upper left', fontsize=8.3, frameon=True)

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=220, facecolor='white')
    plt.close(fig)
    png_bytes = buf.getvalue()
    notes_fig = DESIGN / 'notes' / 'fig' / 'hfss_1port_s21_s11.png'
    for p in (re_folder / 'figures' / 'hfss_1port_s21_s11.png', notes_fig):
        _long_path(p.parent).mkdir(parents=True, exist_ok=True)
        _long_path(p).write_bytes(png_bytes)

    finish_run(
        re_folder,
        status='completed',
        summary=dict(metrics=metrics, figure=str(notes_fig.relative_to(ROOT).as_posix())),
    )
    print('DONE', flush=True)


if __name__ == '__main__':
    main()
