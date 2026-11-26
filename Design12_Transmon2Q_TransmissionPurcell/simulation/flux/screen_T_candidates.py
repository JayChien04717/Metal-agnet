# -*- coding: utf-8 -*-
r"""Sequential HFSS screening of shortlisted 'T' flux-line + SQUID loop candidates for Design12.

Runs one AEDT session at a time (8 adaptive passes at 100 MHz, 4 discrete points 10-100 MHz,
recorded as role='screening') on the 6-port 3x3 matrix fixture (Q1, TC, Q2), plus two quick
single-device L_geo measurements (rect_jj2 Lj=1e-6 nH) to measure pure geometric loop inductance.

Outputs a consolidated comparison summary and figure under Design12/results/flux/.
"""
from __future__ import annotations

import json
import subprocess
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

import prescreen_biot_savart as BS
from run_record import create_run, finish_run

PY = sys.executable
HFSS_SCRIPT = DESIGN_DIR / 'simulation' / 'flux' / 'hfss_mutual_d12.py'
BASELINE_RUN = DESIGN_DIR / 'results' / 'flux' / '20261126_220549_hfss_mutual_L_r4_p12_2a1f697d'

CANDIDATES = [
    dict(
        tag='T_r4loop_p8',
        label='T + r4 loop (30x42)',
        fbl=dict(fbl_type='T', fbl_x='40um', fbl_length='130um', fbl_sep='10.5um'),
        squid=dict(jj_gap='48um', squid_spacing='30um', squid_bar_height='6um'),
    ),
    dict(
        tag='T_c36_24_p8',
        label='T + compact 24x30',
        fbl=dict(fbl_type='T', fbl_x='32um', fbl_length='110um', fbl_sep='10.5um'),
        squid=dict(jj_gap='36um', squid_spacing='24um', squid_bar_height='6um'),
    ),
    dict(
        tag='T_c36_30_p8',
        label='T + square 30x30',
        fbl=dict(fbl_type='T', fbl_x='35um', fbl_length='120um', fbl_sep='10.5um'),
        squid=dict(jj_gap='36um', squid_spacing='30um', squid_bar_height='6um'),
    ),
    dict(
        tag='T_c36_36_p8',
        label='T + near-sq 36x30',
        fbl=dict(fbl_type='T', fbl_x='38um', fbl_length='130um', fbl_sep='10.5um'),
        squid=dict(jj_gap='36um', squid_spacing='36um', squid_bar_height='6um'),
    ),
    dict(
        tag='T_c32_36_bh4_p8',
        label='T + ultra-comp 36x28 (bh4)',
        fbl=dict(fbl_type='T', fbl_x='38um', fbl_length='130um', fbl_sep='10.5um'),
        squid=dict(jj_gap='32um', squid_spacing='36um', squid_bar_height='4um'),
    ),
]


def run_hfss(args):
    cmd = [PY, '-X', 'utf8', '-B', str(HFSS_SCRIPT)] + args
    print('>>>', ' '.join(cmd), flush=True)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stdout[-2000:], flush=True)
        print(proc.stderr[-2000:], flush=True)
        raise RuntimeError(f'HFSS failed: {args}')
    run_folder = None
    for line in proc.stdout.splitlines():
        if line.startswith('RUN_FOLDER='):
            run_folder = Path(line.split('=', 1)[1].strip())
    if not run_folder:
        raise RuntimeError('Could not parse RUN_FOLDER from output')
    summary = json.loads((run_folder / 'summaries' / 'summary.json').read_text(encoding='utf-8'))
    return run_folder, summary


def main():
    import hfss_mutual_d12 as HM
    # Pre-validate all candidates in-process (no stray dry-run folders)
    for c in CANDIDATES:
        ov = dict(fbl=c['fbl'], squid=c['squid'])
        d = HM.build_core_design(ov)
        assert all(nm in d.components for nm in ('Q1', 'TC', 'Q2'))
    print('All 5 candidates passed geometry validation!', flush=True)

    results = []
    # Include baseline L-r4
    if BASELINE_RUN.exists():
        b_sum = json.loads((BASELINE_RUN / 'summaries' / 'summary.json').read_text(encoding='utf-8'))
        results.append(dict(
            tag='L_r4_p12', label='L-r4 baseline (30x42)', run=str(BASELINE_RUN),
            fbl=dict(fbl_type='L', fbl_x='12um', fbl_length='45um', fbl_sep='10.5um'),
            squid=dict(jj_gap='48um', squid_spacing='30um', squid_bar_height='6um'),
            M_BS_pH=BS.mutual_bs_pH('L', 12.0, 45.0, 10.5, 48.0, 30.0, 6.0),
            L_geo_est_pH=BS.loop_self_inductance_pH(30.0, 42.0),
            summary=b_sum,
        ))

    for c in CANDIDATES:
        ov = dict(fbl=c['fbl'], squid=c['squid'])
        folder, s = run_hfss([
            '--tag', c['tag'], '--mode', 'matrix',
            '--min-passes', '8', '--max-passes', '8', '--min-converged', '2',
            '--max-delta-s', '0.0005', '--npts', '4', '--role', 'screening',
            '--overrides', json.dumps(ov),
        ])
        fx = float(c['fbl']['fbl_x'].rstrip('um'))
        fl = float(c['fbl']['fbl_length'].rstrip('um'))
        sep = float(c['fbl']['fbl_sep'].rstrip('um'))
        jg = float(c['squid']['jj_gap'].rstrip('um'))
        sp = float(c['squid']['squid_spacing'].rstrip('um'))
        bh = float(c['squid']['squid_bar_height'].rstrip('um'))
        m_bs = BS.mutual_bs_pH('T', fx, fl, sep, jg, sp, bh)
        l_est = BS.loop_self_inductance_pH(sp, jg - bh)
        results.append(dict(
            tag=c['tag'], label=c['label'], run=str(folder),
            fbl=c['fbl'], squid=c['squid'], M_BS_pH=m_bs, L_geo_est_pH=l_est, summary=s,
        ))
        print(f"DONE {c['tag']}: M_self={s['M_self_pH']}, xtalk={s['crosstalk_pct']}", flush=True)

    # Measure pure geometric L_loop on Q1 for r4 (30x42) and square (30x30)
    lgeo_runs = {}
    for ltag, sq in (
        ('lgeo_r4_30x42', dict(jj_gap='48um', squid_spacing='30um', squid_bar_height='6um')),
        ('lgeo_c36_30x30', dict(jj_gap='36um', squid_spacing='30um', squid_bar_height='6um')),
        ('lgeo_c36_36x30', dict(jj_gap='36um', squid_spacing='36um', squid_bar_height='6um')),
    ):
        ov = dict(fbl=dict(fbl_type='T', fbl_x='35um', fbl_length='120um', fbl_sep='10.5um'), squid=sq)
        folder, s = run_hfss([
            '--tag', ltag, '--mode', 'single', '--device', 'Q1',
            '--min-passes', '8', '--max-passes', '8', '--min-converged', '2',
            '--npts', '2', '--role', 'screening', '--measure-lgeo',
            '--overrides', json.dumps(ov),
        ])
        lgeo_runs[ltag] = dict(run=str(folder), L_geo_HFSS_pH=s['L_loop_nH']['Q1'] * 1e3)
        print(f"DONE {ltag}: L_geo_HFSS = {s['L_loop_nH']['Q1']*1e3:.1f} pH", flush=True)

    agg_dir = create_run(
        'Design12_Transmon2Q_TransmissionPurcell', 'flux', 'screen_T_candidates_summary',
        entrypoint=__file__,
        inputs=[HFSS_SCRIPT],
        parameters=dict(n_candidates=len(CANDIDATES)),
    )
    (agg_dir / 'summaries').mkdir(parents=True, exist_ok=True)
    out = dict(candidates=results, lgeo_checks=lgeo_runs)
    (agg_dir / 'summaries' / 'screening_comparison.json').write_text(
        json.dumps(out, indent=2, default=float), encoding='utf-8'
    )
    finish_run(agg_dir, status='completed', summary=dict(
        lgeo_checks=lgeo_runs,
        table=[dict(
            tag=r['tag'], label=r['label'],
            M_Q1=r['summary']['M_self_pH']['Q1'],
            M_TC=r['summary']['M_self_pH']['TC'],
            M_Q2=r['summary']['M_self_pH']['Q2'],
            M_BS=r['M_BS_pH'],
            I_Phi0_Q1=r['summary']['I_phi0_mA']['Q1'],
            I_Phi0_TC=r['summary']['I_phi0_mA']['TC'],
            max_xtalk_pct=max(abs(v) for v in r['summary']['crosstalk_pct'].values()),
        ) for r in results],
    ))
    print('SCREENING_SUMMARY_DIR=' + str(agg_dir), flush=True)


if __name__ == '__main__':
    main()
