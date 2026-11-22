# -*- coding: utf-8 -*-
"""Re-analyze raw HFSS 50-Ohm S-parameter CSV for Design10 (1-port 3-stage reflection BPF),
applying exact ideal-open termination (Gamma_2 = +1, zero resistive loading) at R3_open to
export:
  1. Unloaded 1-port reflection S11 (= VNA S21 through a circulator): magnitude, phase, tau_g
  2. Internal open-end voltage isolation |V_open / V_inc| = |2*S21 / (1 - S22)| (NOT a physical Port 2)
  3. Pure HFSS-only 3-panel plot (no analytic overlay) + .csv + .s1p.
"""

import io
import json
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


def parse_complex_col(series):
    if np.iscomplexobj(series):
        return series.to_numpy(complex)
    return series.map(lambda v: complex(str(v).replace('i', 'j').strip())).to_numpy(complex)


def terminate_ideal_open(s_df):
    """Exact ideal open-circuit termination (Gamma_2 = +1) at R3_open."""
    s11 = parse_complex_col(s_df['S11'])
    s21 = parse_complex_col(s_df['S21'])
    s22 = parse_complex_col(s_df['S22'])
    denom = 1.0 - s22
    s11_open = s11 + (s21 ** 2) / denom
    v_open_over_vinc = (2.0 * s21) / denom
    return s11_open, v_open_over_vinc


def main():
    parent_name = sys.argv[1] if len(sys.argv) > 1 else '20261121_162252_hfss_driven_1port_3stage_9d5ee50f'
    parent = DESIGN / 'results' / 'resonator' / parent_name
    s_df = pd.read_csv(_long_path(parent / 'raw' / 'S_params_raw50.csv'), index_col=0)
    conv = pd.read_csv(_long_path(parent / 'raw' / 'hfss_convergence.csv'))

    re_folder = create_run(
        DESIGN.name,
        'resonator',
        'reanalysis_hfss_1port_s11_vopen',
        entrypoint=__file__,
        inputs=[parent / 'raw' / 'S_params_raw50.csv', parent / 'raw' / 'hfss_convergence.csv'],
        upstream=[str(parent.relative_to(ROOT).as_posix())],
        parameters=dict(
            parent_run=parent.name,
            open_termination='ideal_open_gamma2_plus1',
            geometry_revision='2026-11-21_r2',
        ),
    )
    rec = json.loads(_long_path(re_folder / 'run.json').read_text(encoding='utf-8'))
    rec['geometry_revision'] = '2026-11-21_r2'
    _long_path(re_folder / 'run.json').write_text(json.dumps(rec, indent=2) + '\n', encoding='utf-8')
    print('REANALYSIS_FOLDER=' + re_folder.name, flush=True)

    freqs = s_df.index.to_numpy(float)
    s11_open, v_ratio = terminate_ideal_open(s_df)
    v_db = 20.0 * np.log10(np.maximum(np.abs(v_ratio), 1e-12))
    v_norm = v_db - float(np.max(v_db))
    s11_db = 20.0 * np.log10(np.maximum(np.abs(s11_open), 1e-12))
    ph_rad = np.unwrap(np.angle(s11_open))
    tau_ns = np.abs(np.gradient(ph_rad) / (2.0 * np.pi * np.gradient(freqs) * 1e9)) * 1e9

    _long_path(re_folder / 'summaries').mkdir(parents=True, exist_ok=True)
    out_tbl = pd.DataFrame(dict(
        freq_GHz=freqs,
        S11_dB=np.round(s11_db, 5),
        S11_deg=np.round(np.degrees(np.angle(s11_open)), 3),
        S11_unwrapped_deg=np.round(np.degrees(ph_rad), 3),
        tau_g_ns=np.round(tau_ns, 5),
        Vopen_over_Vinc_norm_dB=np.round(v_norm, 4),
        Vopen_over_Vinc_dB=np.round(v_db, 4),
    ))
    out_tbl.to_csv(_long_path(re_folder / 'summaries' / 'hfss_1port_2_7p5GHz.csv'), index=False)

    s1p_lines = [
        f'! {DESIGN.name} HFSS Driven ideal-open 1-port S11 (= Circulator + VNA S21) (2026-11-21_r2)',
        f'! Parent HFSS run: {parent.name} (3 passes, 10 um mesh, 451272 tets, Delta S = 0.033366)',
        '# GHz S RI R 50',
    ]
    for f, s11v in zip(freqs, s11_open):
        s1p_lines.append(f'{f:.6f} {s11v.real:.8e} {s11v.imag:.8e}')
    _long_path(re_folder / 'summaries' / 'hfss_1port_s11.s1p').write_text('\n'.join(s1p_lines) + '\n', encoding='utf-8')

    pb = (freqs >= 5.5) & (freqs <= 6.5)
    sb = (freqs >= 2.0) & (freqs <= 4.0)
    idx3 = np.where(v_norm >= -3.0)[0]
    f3_lo, f3_hi = float(freqs[idx3[0]]), float(freqs[idx3[-1]])
    last_ds = float(conv['Max Mag. Delta S'].dropna().iloc[-1])
    tets = int(conv['Solved Elements'].dropna().iloc[-1])

    def at_f(ft):
        return float(v_norm[np.argmin(np.abs(freqs - ft))])

    pb_min = float(np.min(v_norm[pb]))
    sb_max = float(np.max(v_norm[sb]))

    metrics = dict(
        open_termination='ideal_open_gamma2_plus1',
        s11_mag_range_dB=[round(float(np.min(s11_db)), 5), round(float(np.max(s11_db)), 5)],
        band_3dB_GHz=[round(f3_lo, 4), round(f3_hi, 4)],
        f0_center_3dB_GHz=round(0.5 * (f3_lo + f3_hi), 4),
        pb_5p5_6p5_min_vopen_dB=round(pb_min, 3),
        sb_2_4_max_vopen_dB=round(sb_max, 2),
        sb_2_4_meets_40dB=bool(sb_max <= -40.0),
        vopen_at_2GHz_dB=round(at_f(2.0), 2),
        vopen_at_3GHz_dB=round(at_f(3.0), 2),
        vopen_at_4GHz_dB=round(at_f(4.0), 2),
        vopen_at_5GHz_dB=round(at_f(5.0), 2),
        vopen_at_6GHz_dB=round(at_f(6.0), 2),
        vopen_at_7GHz_dB=round(at_f(7.0), 2),
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

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(9.2, 9.4), sharex=True)

    # Panel 1 (Top): Circulator + VNA S21 (= chip 1-port S11) magnitude & unwrapped phase
    ax0.axvspan(2.0, 4.0, color='#fadbd8', alpha=0.38, label='Stopband: 2–4 GHz')
    ax0.axvspan(5.5, 6.5, color='#d6eaf8', alpha=0.48, label='Passband: 5.5–6.5 GHz')
    ax0.plot(freqs, s11_db, color='#c0392b', lw=1.9, label='|S21_VNA| = |S11| (0 dB, lossless 1-port)')
    ax0.set_ylim(-2.0, 0.5)
    ax0.set_ylabel('|S21_VNA| = |S11| (dB)', color='#c0392b')
    ax0.tick_params(axis='y', labelcolor='#c0392b')
    ax0.grid(True, alpha=0.3)
    ax0r = ax0.twinx()
    ax0r.plot(freqs, np.degrees(ph_rad), color='#7d3c98', lw=1.7, ls='--',
              label='∠S21_VNA = ∠S11 unwrapped (deg, right)')
    ax0r.set_ylabel('∠S21_VNA = ∠S11 (deg)', color='#7d3c98')
    ax0r.tick_params(axis='y', labelcolor='#7d3c98')
    h1, l1 = ax0.get_legend_handles_labels()
    h2, l2 = ax0r.get_legend_handles_labels()
    ax0.legend(h1 + h2, l1 + l2, loc='lower left', fontsize=8.0, frameon=True)
    ax0.set_title(
        f'{DESIGN.name} · 1-Port 3-Stage Reflection BPF · HFSS Driven Only '
        f'(3 passes, 10 μm, ΔS={last_ds:.4f}, {tets:,} tets)\n'
        f'Circulator + VNA Measurement: S21_VNA = Chip S11 (Ideal Open End, No Dip Because Lossless)',
        fontsize=8.6,
        weight='bold',
    )

    # Panel 2 (Middle): Circulator + VNA reflection group delay
    ax1.axvspan(2.0, 4.0, color='#fadbd8', alpha=0.38)
    ax1.axvspan(5.5, 6.5, color='#d6eaf8', alpha=0.48)
    ax1.plot(freqs, tau_ns, color='#1e8449', lw=1.9,
             label='Reflection group delay τ_g(S11) = -(1/2π) d∠S11/df')
    ax1.set_yscale('log')
    ax1.set_ylabel('Group Delay τ_g(S11) (ns)')
    ax1.set_title(
        f'Circulator + VNA Reflection Group Delay τ_g(S11) · '
        f'5.5–6.5 GHz median={metrics["tau_g_5p5_6p5_median_ns"]:.2f} ns | @3 GHz={metrics["tau_g_3GHz_ns"]:.3f} ns',
        fontsize=8.9,
        weight='bold',
    )
    ax1.grid(True, alpha=0.3, which='both')
    ax1.legend(loc='upper left', fontsize=8.2, frameon=True)

    # Panel 3 (Bottom): Internal open-end voltage isolation |V_open / V_inc| (NOT a physical Port 2)
    ax2.axvspan(2.0, 4.0, color='#fadbd8', alpha=0.42, label='Stopband target: 2–4 GHz (≤ -40 dB)')
    ax2.axvspan(5.5, 6.5, color='#d6eaf8', alpha=0.50, label='Passband target: 5.5–6.5 GHz')
    ax2.plot(freqs, v_norm, color='#1b4f72', lw=2.0,
             label='Internal open-end voltage transfer |V_open / V_inc| (normalized, no Port 2)')
    ax2.axhline(-40.0, color='#922b21', lw=0.9, ls=':')
    ax2.axhline(-3.0, color='#555555', lw=0.8, ls=':')
    ax2.set_ylim(-75, 2)
    ax2.set_xlim(freqs.min(), freqs.max())
    ax2.set_xlabel('Frequency (GHz)')
    ax2.set_ylabel('Normalized |V_open / V_inc| (dB)')
    ax2.set_title(
        f'Internal Open-End Isolation at R3_open (Not a Physical Port 2) · '
        f'-3 dB: {f3_lo:.2f}–{f3_hi:.2f} GHz | 5.5–6.5G min={pb_min:.2f} dB | 2–4G max={sb_max:.1f} dB',
        fontsize=8.7,
        weight='bold',
    )
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='lower right', fontsize=8.1, frameon=True)

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
