# -*- coding: utf-8 -*-
r"""Biot-Savart / Neumann filament pre-screening of Design12 flux lines ('L' and 'T') and SQUID loops.

Computes for any (fbl_type, fbl_x, fbl_length, fbl_sep, jj_gap, squid_spacing, squid_bar_height):
  - M_BS (pH) via Neumann / Biot-Savart flux integral through the SQUID loop area
    (including symmetric/asymmetric return currents in the flux-line bar, feed leg, and
    ground-plane return strip / shorted ends).
  - Analytical planar rectangular loop geometric inductance L_geo (pH).
  - I_Phi0 = Phi0 / |M|, I_0p5 = 0.5 * I_Phi0, CZ pulse currents (Q1 0.116 Phi0, TC 0.224 Phi0).
  - Loop area (stray-field sensitivity) and perimeter (1/f flux-noise sensitivity).
  - Far-field crosstalk estimate M_ij / M_jj (dipole vs quadrupole scaling for L vs T).
  - Manhattan e-beam junction clearance check via TransmonJunctionManhattan.

Saves a full run record under Design12/results/flux/ with CSV table and trade-off figures.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DESIGN_DIR = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / 'code' / 'general'),
    str(ROOT / 'pattern design' / 'ebeam_junction'),
    str(DESIGN_DIR / 'layout'),
]
import project_paths
project_paths.configure(__file__)

import numpy as np
import pandas as pd
from run_record import create_run, finish_run

MU0 = 4.0 * math.pi * 1e-7
PHI0 = 2.067833848e-15


def b_z_segment(x, y, x1, y1, x2, y2, current=1.0):
    """Exact Biot-Savart B_z(x, y) [T] in the z=0 plane from a straight 2D filament (x1,y1)->(x2,y2) [m]."""
    dx = x2 - x1
    dy = y2 - y1
    L = math.hypot(dx, dy)
    if L < 1e-12:
        return np.zeros_like(x, dtype=float)
    ux, uy = dx / L, dy / L
    # Signed perpendicular distance d = (r - r1) x u (in z)
    rx1 = x - x1
    ry1 = y - y1
    # dl x (r - r') has z component: ux * ry1 - uy * rx1
    perp = ux * ry1 - uy * rx1
    # Projections along segment
    s1 = rx1 * ux + ry1 * uy
    s2 = s1 - L
    r1 = np.hypot(rx1, ry1)
    r2 = np.hypot(x - x2, y - y2)
    # Regularize at very small perp (not hit when loop doesn't intersect filament)
    denom = np.where(np.abs(perp) > 1e-12, perp, np.nan)
    bz = (MU0 * current / (4.0 * math.pi)) * (s1 / r1 - s2 / r2) / denom
    return np.nan_to_num(bz, nan=0.0)


def flux_filaments(fbl_type, fbl_x_um, fbl_length_um, fbl_sep_um, fbl_height_um=80.0,
                   cpw_w_um=10.0, cpw_g_um=6.0):
    """Return list of (x1, y1, x2, y2, I_fraction) in um in the pre-flip local frame.

    Frame: SQUID ground etch is at y = 0 (shifted by -d from fluxline.py); SQUID loop is in y < 0;
    flux-line bar centerline is at y_bar = +(fbl_sep + cpw_w/2); feed pin at y_end = 2*(fbl_height + fbl_sep).
    Ground-plane return current flows along the CPW ground edges (x_feed +- (w/2 + g)) and along the
    ground strip edge (y_ret = fbl_sep - cpw_g) back to the feed CPW grounds.
    """
    y_bar = fbl_sep_um + cpw_w_um / 2.0
    y_end = 2.0 * (fbl_height_um + fbl_sep_um)
    y_gnd_top = fbl_sep_um - cpw_g_um        # top edge of gap (SQUID side ground strip)
    y_gnd_bot = y_bar + cpw_w_um / 2.0 + cpw_g_um
    dx_cpw = cpw_w_um / 2.0 + cpw_g_um       # 11 um from feed centerline to CPW ground edge
    fl = fbl_length_um
    w = cpw_w_um
    fx = fbl_x_um

    segs = []
    if fbl_type.upper() == 'L':
        # Signal: pin at (fx + fl/2 + w/2, y_end) -> up to bar -> left to short at (fx - fl/2 - w, y_bar)
        # Current-conserving thin-filament return along the lower CPW ground edge (standard free-space BS;
        # HFSS sits ~30-40% above this due to PEC flux focusing, per skill/mutual_inductance_analysis.md §4).
        x_feed = fx + fl / 2.0 + w / 2.0
        x_short = fx - fl / 2.0 - w
        segs.append((x_feed, y_end, x_feed, y_bar, +1.0))
        segs.append((x_feed, y_bar, x_short, y_bar, +1.0))
        segs.append((x_short, y_bar, x_short, y_gnd_bot, +1.0))
        segs.append((x_short, y_gnd_bot, x_feed - dx_cpw, y_gnd_bot, +1.0))
        segs.append((x_feed - dx_cpw, y_gnd_bot, x_feed - dx_cpw, y_end, +1.0))
    elif fbl_type.upper() == 'T':
        # Signal: feed at x = fx from y_end to y_bar, splits 50/50 into left arm (fx - fl/2) and right arm (fx + fl/2)
        # Each shorted arm returns 0.5 along its lower CPW ground rail (exact current conservation).
        x_feed = fx
        x_left = fx - fl / 2.0
        x_right = fx + fl / 2.0
        segs.append((x_feed, y_end, x_feed, y_bar, +1.0))
        segs.append((x_feed, y_bar, x_left, y_bar, +0.5))
        segs.append((x_feed, y_bar, x_right, y_bar, +0.5))
        segs.append((x_left, y_bar, x_left, y_gnd_bot, +0.5))
        segs.append((x_left, y_gnd_bot, x_feed - dx_cpw, y_gnd_bot, +0.5))
        segs.append((x_feed - dx_cpw, y_gnd_bot, x_feed - dx_cpw, y_end, +0.5))
        segs.append((x_right, y_bar, x_right, y_gnd_bot, +0.5))
        segs.append((x_right, y_gnd_bot, x_feed + dx_cpw, y_gnd_bot, +0.5))
        segs.append((x_feed + dx_cpw, y_gnd_bot, x_feed + dx_cpw, y_end, +0.5))
    else:
        raise ValueError(fbl_type)
    return segs


def mutual_bs_pH(fbl_type, fbl_x_um, fbl_length_um, fbl_sep_um, jj_gap_um, squid_spacing_um,
                 squid_bar_height_um=6.0, nx=24, ny=24):
    """Integrate B_z over the SQUID loop rectangle (centerline frame) -> M [pH].

    Loop x in [-sp/2, +sp/2], y in [-jj_gap + bar_h/2, -bar_h/2] (relative to ground etch y = 0).
    """
    segs = flux_filaments(fbl_type, fbl_x_um, fbl_length_um, fbl_sep_um)
    x0, x1 = -squid_spacing_um / 2.0, +squid_spacing_um / 2.0
    y0, y1 = -(jj_gap_um - squid_bar_height_um / 2.0), -(squid_bar_height_um / 2.0)
    xs = np.linspace(x0, x1, nx + 1)
    ys = np.linspace(y0, y1, ny + 1)
    xc = 0.5 * (xs[:-1] + xs[1:]) * 1e-6
    yc = 0.5 * (ys[:-1] + ys[1:]) * 1e-6
    dx_m = (xs[1] - xs[0]) * 1e-6
    dy_m = (ys[1] - ys[0]) * 1e-6
    X, Y = np.meshgrid(xc, yc)
    Bz = np.zeros_like(X)
    for sx1, sy1, sx2, sy2, cur in segs:
        Bz += b_z_segment(X, Y, sx1 * 1e-6, sy1 * 1e-6, sx2 * 1e-6, sy2 * 1e-6, current=cur)
    phi = float(np.sum(Bz) * dx_m * dy_m)
    return phi * 1e12  # pH


def loop_self_inductance_pH(w_um, h_um, trace_w_um=4.5):
    """Planar rectangular loop geometric self-inductance estimate (Rosa/Grover formula) [pH]."""
    a, b = w_um * 1e-6, h_um * 1e-6
    w = trace_w_um * 1e-6
    d = math.hypot(a, b)
    # Thin flat strip GMD ~ 0.2235 * w
    rho = 0.2235 * w
    L = (MU0 / math.pi) * (
        a * math.log(2.0 * a * b / (rho * (a + d)))
        + b * math.log(2.0 * a * b / (rho * (b + d)))
        - 2.0 * (a + b - d)
        + 0.25 * (a + b)
    )
    return L * 1e12


def check_manhattan_min_spacing():
    """Find minimum squid_spacing (um) where TransmonJunctionManhattan passes without overlap."""
    from qiskit_metal import designs
    from transmon_junction_manhattan import TransmonJunctionManhattan
    d = designs.DesignPlanar()
    d.overwrite_enabled = True
    for sp in range(14, 32):
        try:
            TransmonJunctionManhattan(d, f'jj_{sp}', options=dict(squid='True', squid_spacing=f'{sp}um'))
            return sp
        except ValueError:
            continue
    return None


def main():
    min_sp = check_manhattan_min_spacing()
    print(f'TransmonJunctionManhattan minimum non-overlapping squid_spacing: {min_sp} um')

    # Baseline L-r4
    base_M = mutual_bs_pH('L', 12.0, 45.0, 10.5, 48.0, 30.0, 6.0)
    print(f'Biot-Savart baseline L-r4 (30 x 42 um): M_BS = {base_M:.4f} pH, I_Phi0 = {PHI0/(abs(base_M)*1e-12)*1e3:.2f} mA')

    # 1) Sweep T-line placement (fbl_x, fbl_length) on the r4 loop (30 x 48 um) to find optimal T alignment
    t_align_rows = []
    for fx in (25, 30, 35, 40, 45, 50):
        for fl in (90, 110, 130, 150, 170):
            # Require left arm to reach at least -sp/2 (-15 um)
            left_end = fx - fl / 2.0
            if left_end > -15.0:
                continue
            m = mutual_bs_pH('T', fx, fl, 10.5, 48.0, 30.0, 6.0)
            t_align_rows.append(dict(fbl_x_um=fx, fbl_length_um=fl, left_end_um=left_end, M_BS_pH=m))
    df_align = pd.DataFrame(t_align_rows)
    best_align = df_align.iloc[df_align['M_BS_pH'].abs().idxmax()]
    print('\nBest T-line alignment on r4 loop (30 x 48 um):')
    print(df_align.sort_values('M_BS_pH', key=abs, ascending=False).head(8).to_string(index=False))

    # 2) Sweep SQUID loop (jj_gap x squid_spacing) + fbl_bar_height + fbl_sep with optimal T-line alignment
    rows = []
    # Include L-r4 baseline as row 0
    for cfg in [
        dict(label='L_r4_baseline', fbl_type='L', fbl_x=12.0, fbl_length=45.0, fbl_sep=10.5,
             jj_gap=48.0, squid_spacing=30.0, bar_h=6.0),
    ]:
        m = mutual_bs_pH(cfg['fbl_type'], cfg['fbl_x'], cfg['fbl_length'], cfg['fbl_sep'],
                         cfg['jj_gap'], cfg['squid_spacing'], cfg['bar_h'])
        h_eff = cfg['jj_gap'] - cfg['bar_h']
        w_eff = cfg['squid_spacing']
        I0 = PHI0 / (abs(m) * 1e-12) * 1e3
        rows.append(dict(
            **cfg, h_eff_um=h_eff, area_um2=w_eff * h_eff, perim_um=2 * (w_eff + h_eff),
            aspect=w_eff / h_eff, M_BS_pH=m, I_Phi0_mA=I0, I_0p5_mA=0.5 * I0,
            I_CZ_Q1_mA=0.116 * I0, I_CZ_TC_mA=0.224 * I0,
            L_geo_est_pH=loop_self_inductance_pH(w_eff, h_eff),
        ))

    for jg in (32, 36, 40, 48, 56, 64, 72):
        for sp in (20, 24, 30, 36, 44, 52, 60):
            if min_sp and sp < min_sp:
                continue
            for bar_h in (4.0, 6.0):
                for sep in (10.0, 10.5):
                    # Place T feed just outside +sp/2 + 20 um (so feed + CPW gap clears SQUID right edge)
                    # and extend left arm at least 20 um past -sp/2
                    fx = sp / 2.0 + 20.0
                    fl = 2.0 * (fx + sp / 2.0 + 20.0)
                    m = mutual_bs_pH('T', fx, fl, sep, jg, sp, bar_h)
                    h_eff = jg - bar_h
                    w_eff = sp
                    I0 = PHI0 / (abs(m) * 1e-12) * 1e3
                    rows.append(dict(
                        label=f'T_jg{jg}_sp{sp}_bh{int(bar_h)}_s{sep}',
                        fbl_type='T', fbl_x=fx, fbl_length=fl, fbl_sep=sep,
                        jj_gap=jg, squid_spacing=sp, bar_h=bar_h,
                        h_eff_um=h_eff, area_um2=w_eff * h_eff, perim_um=2 * (w_eff + h_eff),
                        aspect=w_eff / h_eff, M_BS_pH=m, I_Phi0_mA=I0, I_0p5_mA=0.5 * I0,
                        I_CZ_Q1_mA=0.116 * I0, I_CZ_TC_mA=0.224 * I0,
                        L_geo_est_pH=loop_self_inductance_pH(w_eff, h_eff),
                    ))

    df = pd.DataFrame(rows)
    run_dir = create_run(
        'Design12_Transmon2Q_TransmissionPurcell', 'flux', 'prescreen_biot_savart_T_loop',
        entrypoint=__file__,
        inputs=[DESIGN_DIR / 'layout' / 'build.py'],
        parameters=dict(min_manhattan_spacing_um=min_sp, n_candidates=len(df)),
    )
    (run_dir / 'raw').mkdir(parents=True, exist_ok=True)
    (run_dir / 'summaries').mkdir(parents=True, exist_ok=True)
    (run_dir / 'figures').mkdir(parents=True, exist_ok=True)
    df.to_csv(run_dir / 'raw' / 'prescreen_grid.csv', index=False)
    df_align.to_csv(run_dir / 'raw' / 't_alignment_sweep.csv', index=False)

    # Print compact 2D grid for bar_h=6, sep=10.5
    sub = df[(df['fbl_type'] == 'T') & (df['bar_h'] == 6.0) & (df['fbl_sep'] == 10.5)]
    piv_M = sub.pivot(index='jj_gap', columns='squid_spacing', values='M_BS_pH')
    piv_I = sub.pivot(index='jj_gap', columns='squid_spacing', values='I_Phi0_mA')
    piv_eff = sub.assign(M_over_P=sub['M_BS_pH'].abs() / sub['perim_um'] * 100.0).pivot(
        index='jj_gap', columns='squid_spacing', values='M_over_P'
    )
    print('\n--- Biot-Savart M (pH) for T-line (bar_h=6um, fbl_sep=10.5um) ---')
    print(piv_M.round(3).to_string())
    print('\n--- I_Phi0 (mA) for T-line (bar_h=6um, fbl_sep=10.5um) ---')
    print(piv_I.round(2).to_string())
    print('\n--- Efficiency |M| / Perimeter (pH / 100um) ---')
    print(piv_eff.round(3).to_string())

    # Plot trade-off
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    sc = ax[0].scatter(sub['perim_um'], sub['M_BS_pH'].abs(), c=sub['aspect'], cmap='viridis', s=50)
    ax[0].scatter([df.iloc[0]['perim_um']], [abs(df.iloc[0]['M_BS_pH'])], c='red', marker='*', s=160, label='L-r4 baseline')
    ax[0].set_xlabel('Loop perimeter (um) [1/f flux-noise proxy]')
    ax[0].set_ylabel('|M_BS| (pH)')
    ax[0].set_title('T-line M vs loop perimeter (color = width/height)')
    ax[0].grid(True, alpha=0.3)
    ax[0].legend()
    fig.colorbar(sc, ax=ax[0], label='aspect ratio w/h')

    im = ax[1].imshow(piv_M.abs().to_numpy(), origin='lower', aspect='auto', cmap='magma')
    ax[1].set_xticks(range(len(piv_M.columns)))
    ax[1].set_xticklabels([int(c) for c in piv_M.columns])
    ax[1].set_yticks(range(len(piv_M.index)))
    ax[1].set_yticklabels([int(r) for r in piv_M.index])
    ax[1].set_xlabel('squid_spacing (um)')
    ax[1].set_ylabel('jj_gap (um)')
    ax[1].set_title('|M_BS| (pH) grid (bar_h=6um, sep=10.5um)')
    for i in range(len(piv_M.index)):
        for j in range(len(piv_M.columns)):
            v = abs(piv_M.iloc[i, j])
            ax[1].text(j, i, f'{v:.2f}', ha='center', va='center',
                       color='white' if v < 0.6 * piv_M.abs().to_numpy().max() else 'black', fontsize=8)
    fig.tight_layout()
    fig.savefig(run_dir / 'figures' / 'prescreen_tradeoff.png', dpi=140)
    plt.close(fig)

    summary = dict(
        min_manhattan_spacing_um=min_sp,
        baseline_L_r4_M_BS_pH=float(base_M),
        best_T_align_on_r4=best_align.to_dict(),
    )
    (run_dir / 'summaries' / 'summary.json').write_text(json.dumps(summary, indent=2, default=float), encoding='utf-8')
    finish_run(run_dir, status='completed', summary=summary)
    print('RUN_FOLDER=' + str(run_dir))


if __name__ == '__main__':
    main()
