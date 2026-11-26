# -*- coding: utf-8 -*-
r"""Aggregate L-r4 baseline and shortlisted T-line + SQUID HFSS screening runs into a table & figure."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DESIGN_DIR = Path(__file__).resolve().parents[2]
sys.path[:0] = [
    str(ROOT / 'code' / 'general'),
    str(DESIGN_DIR / 'simulation' / 'flux'),
]
import project_paths
project_paths.configure(__file__)

import numpy as np
import pandas as pd
import prescreen_biot_savart as BS
from run_record import create_run, finish_run


def main():
    flux_dir = DESIGN_DIR / 'results' / 'flux'
    tags = [
        ('L_r4_p12', 'L-r4 baseline'),
        ('T_r4loop_p8', 'T + r4 loop'),
        ('T_c36_24_p8', 'T + 24x30'),
        ('T_c36_30_p8', 'T + 30x30 (sq)'),
        ('T_c36_36_p8', 'T + 36x30'),
        ('T_c32_36_bh4_p8', 'T + 36x28 (bh4)'),
    ]
    lgeo_tags = ['lgeo_r4_30x42', 'lgeo_c36_30x30', 'lgeo_c36_36x30']
    lgeo = {}
    for lt in lgeo_tags:
        m = sorted(flux_dir.glob(f'*hfss_mutual_{lt}_*'))
        if m:
            s = json.loads((m[-1] / 'summaries' / 'summary.json').read_text(encoding='utf-8'))
            lgeo[lt] = dict(run=m[-1].name, L_geo_pH=s['L_geo_pH']['Q1'])

    rows = []
    for tag, label in tags:
        m = sorted(flux_dir.glob(f'*hfss_mutual_{tag}_*'))
        if not m:
            continue
        p = m[-1]
        s = json.loads((p / 'summaries' / 'summary.json').read_text(encoding='utf-8'))
        lg = s['loop_geom']['Q1']
        jg, sp, bh = lg['jj_gap_um'], lg['squid_spacing_um'], lg['squid_bar_height_um']
        h_eff = jg - bh
        area = sp * h_eff
        perim = 2.0 * (sp + h_eff)
        fbl = s['overrides'].get('fbl', dict(fbl_type='L', fbl_x='12um', fbl_length='45um', fbl_sep='10.5um'))
        m_bs = BS.mutual_bs_pH(
            fbl['fbl_type'], float(fbl['fbl_x'].rstrip('um')), float(fbl['fbl_length'].rstrip('um')),
            float(fbl['fbl_sep'].rstrip('um')), jg, sp, bh,
        )
        l_est = BS.loop_self_inductance_pH(sp, h_eff)
        ms = s['M_self_pH']
        i0 = s['I_phi0_mA']
        i5 = s['I_0p5_mA']
        icz = s['I_cz_mA']
        xt = s['crosstalk_pct']
        cv = s['convergence']['last_row']
        rows.append(dict(
            tag=tag, label=label, run=p.name, role=s['role'],
            fbl_type=fbl['fbl_type'], fbl_x_um=float(fbl['fbl_x'].rstrip('um')),
            fbl_length_um=float(fbl['fbl_length'].rstrip('um')),
            jj_gap_um=jg, squid_spacing_um=sp, bar_h_um=bh, h_eff_um=h_eff,
            area_um2=area, perim_um=perim, aspect=sp / h_eff,
            M_Q1_pH=abs(ms['Q1']), M_TC_pH=abs(ms['TC']), M_Q2_pH=abs(ms['Q2']),
            M_BS_pH=m_bs, HFSS_over_BS=abs(ms['Q1']) / m_bs,
            I0_Q1_mA=i0['Q1'], I0_TC_mA=i0['TC'], I0_Q2_mA=i0['Q2'],
            I0p5_Q1_mA=i5['Q1'], I0p5_TC_mA=i5['TC'],
            Icz_Q1_mA=icz['Q1'], Icz_TC_mA=icz['TC'],
            L_geo_est_pH=l_est,
            M_spread_Q1_pct=s['M_spread_pct']['Q1'],
            xt_Q1_TC=xt['Q1<-TC'], xt_Q2_TC=xt['Q2<-TC'],
            xt_TC_Q1=xt['TC<-Q1'], xt_TC_Q2=xt['TC<-Q2'],
            xt_Q1_Q2=xt['Q1<-Q2'], xt_Q2_Q1=xt['Q2<-Q1'],
            passes=s['convergence']['passes'],
            elements=int(cv['Solved Elements']),
            delta_s=float(cv['Max Mag. Delta S']),
        ))

    df = pd.DataFrame(rows)
    out_dir = create_run(
        'Design12_Transmon2Q_TransmissionPurcell', 'flux', 'T_screen_agg',
        entrypoint=__file__, parameters=dict(n_rows=len(df)),
    )
    (out_dir / 'summaries').mkdir(parents=True, exist_ok=True)
    (out_dir / 'figures').mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / 'summaries' / 'table.csv', index=False)
    (out_dir / 'summaries' / 'summary.json').write_text(
        json.dumps(dict(rows=rows, lgeo_hfss=lgeo), indent=2), encoding='utf-8'
    )

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.5))
    # Panel 1: |M| vs Perimeter (1/f flux-noise proxy)
    for _, r in df.iterrows():
        mk = '*' if r['fbl_type'] == 'L' else 'o'
        sz = 180 if r['fbl_type'] == 'L' else 90
        c = 'crimson' if r['fbl_type'] == 'L' else 'tab:blue'
        ax[0].scatter([r['perim_um']], [r['M_Q1_pH']], marker=mk, s=sz, color=c)
        ax[0].annotate(r['label'], (r['perim_um'], r['M_Q1_pH']), textcoords='offset points', xytext=(6, 4), fontsize=8)
    ax[0].axhspan(0.5, 1.0, color='green', alpha=0.10, label='Suggested M range (0.5–1.0 pH)')
    ax[0].set_xlabel('Loop perimeter (um) [1/f flux-noise proxy]')
    ax[0].set_ylabel('|M_Q1| HFSS (pH)')
    ax[0].set_title('Mutual inductance vs loop perimeter')
    ax[0].grid(True, alpha=0.3)
    ax[0].legend(fontsize=8)

    # Panel 2: Crosstalk comparison
    x = np.arange(len(df))
    w = 0.14
    for k, (col, lbl) in enumerate([
        ('xt_Q1_TC', 'Q1<-TC'), ('xt_Q2_TC', 'Q2<-TC'),
        ('xt_TC_Q1', 'TC<-Q1'), ('xt_TC_Q2', 'TC<-Q2'), ('xt_Q1_Q2', 'Q1<-Q2'),
    ]):
        ax[1].bar(x + (k - 2) * w, df[col].abs(), width=w, label=lbl)
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(df['label'], rotation=20, ha='right', fontsize=8)
    ax[1].set_ylabel('|M_ij / M_jj| (%)')
    ax[1].set_title('Flux crosstalk (3x3 local HFSS)')
    ax[1].grid(True, axis='y', alpha=0.3)
    ax[1].legend(fontsize=7, ncol=2)

    fig.tight_layout()
    fig.savefig(out_dir / 'figures' / 'tradeoff.png', dpi=140)
    plt.close(fig)

    finish_run(out_dir, status='completed', summary=dict(lgeo_hfss=lgeo, n_rows=len(df)))
    print(df[['label', 'squid_spacing_um', 'h_eff_um', 'area_um2', 'perim_um', 'M_Q1_pH', 'M_TC_pH', 'M_Q2_pH',
              'M_BS_pH', 'HFSS_over_BS', 'I0_Q1_mA', 'I0_TC_mA', 'Icz_Q1_mA', 'Icz_TC_mA', 'passes', 'elements', 'delta_s']].to_string(index=False))
    print('\nCrosstalk (%):')
    print(df[['label', 'xt_Q1_TC', 'xt_Q2_TC', 'xt_TC_Q1', 'xt_TC_Q2', 'xt_Q1_Q2', 'xt_Q2_Q1']].to_string(index=False))
    print('\nL_geo HFSS:', lgeo)
    print('OUT_DIR=' + str(out_dir))


if __name__ == '__main__':
    main()
