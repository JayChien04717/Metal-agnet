# -*- coding: utf-8 -*-
r"""Re-analysis: read HFSS eigenmode adaptive-pass solutions from the AEDT results folder (text *.sd files) and
invert the Purcell filter / readout resonator mode pairs to bare (f_p, f_r, J, kappa_p). No AEDT is opened.

Usage (qmetal or scqenv; analysis only):
    python -X utf8 -B Design11_Transmon2Q_Manhattan/simulation/resonator/eigen_from_sd.py \
        --parent Design11_Transmon2Q_Manhattan/results/resonator/<hfss_eigen run> \
        --results Design11_Transmon2Q_Manhattan/simulation/aedt/D11.aedtresults/E_<hash>_hfss.results

Each DV*_SOL*_M<n>_V*.sd file holds one adaptive pass: Modes(f1, f1'', f2, f2'', ...) in Hz, complex frequency
f + i f''. kappa = 2 f'' and Q = f / (2 f''). Modes with Q < Q_MIN (port / box modes) or outside the readout band
are dropped. The remaining modes are paired in ascending frequency: the lowest pair is Q2 (6.15 GHz target) and the
next pair is Q1 (6.35 GHz). Each pair is inverted with purcell_resonator_pair_sim.invert. The mode with the larger
f'' is the filter-like mode. The source results folder is read only.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'resonator')]

import purcell_resonator_pair_sim as S
from run_record import create_run, finish_run

DESIGN = 'Design11_Transmon2Q_Manhattan'
TARGET = dict(Q1=6.35, Q2=6.15)
BAND = (5.6, 6.9)
Q_MIN = 10.0


def read_passes(results: Path) -> dict:
    out = {}
    for f in sorted(results.glob('*_SOL*_M*_V*.sd')):
        t = f.read_text(encoding='utf-8', errors='replace')
        name = re.search(r"Name='([^']+)'", t).group(1)
        npass = int(re.search(r'AdaptivePass=(\d+)', t).group(1))
        vals = [float(x) for x in re.search(r'Modes\(([^)]*)\)', t).group(1).split(',')]
        modes = [dict(f_GHz=vals[i] / 1e9, im_GHz=vals[i + 1] / 1e9) for i in range(0, len(vals), 2)]
        for m in modes:
            m['kappa_MHz'] = 2 * abs(m['im_GHz']) * 1e3
            m['Q'] = m['f_GHz'] / (2 * abs(m['im_GHz'])) if m['im_GHz'] else float('inf')
        out[name] = dict(file=f.name, adaptive_pass=npass, modes=modes)
    return out


def invert_pass(modes: list) -> dict:
    keep = sorted([m for m in modes if BAND[0] <= m['f_GHz'] <= BAND[1] and m['Q'] >= Q_MIN], key=lambda m: m['f_GHz'])
    dropped = [m for m in modes if m not in keep]
    res = dict(used_modes=keep, dropped_modes=dropped, channels={})
    if len(keep) != 4:
        res['error'] = f'expected 4 readout modes in band, found {len(keep)}'
        return res
    for ch, pair in (('Q2', keep[0:2]), ('Q1', keep[2:4])):
        z1 = complex(pair[0]['f_GHz'], abs(pair[0]['im_GHz']))
        z2 = complex(pair[1]['f_GHz'], abs(pair[1]['im_GHz']))
        inv = S.invert(z1, z2)
        inv['filter_like_mode_GHz'] = max(pair, key=lambda m: abs(m['im_GHz']))['f_GHz']
        inv['resonator_like_mode_GHz'] = min(pair, key=lambda m: abs(m['im_GHz']))['f_GHz']
        inv['df_p_MHz'] = (inv['f_p_GHz'] - TARGET[ch]) * 1e3
        inv['df_r_MHz'] = (inv['f_r_GHz'] - TARGET[ch]) * 1e3
        inv['detuning_fr_minus_fp_MHz'] = (inv['f_r_GHz'] - inv['f_p_GHz']) * 1e3
        res['channels'][ch] = inv
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--parent', required=True, help='parent hfss_eigen run folder (repo-relative)')
    ap.add_argument('--results', required=True, help='AEDT <design>.results folder holding the *.sd pass files')
    ap.add_argument('--tag', default=None)
    a = ap.parse_args()
    parent = (ROOT / a.parent).resolve()
    results = (ROOT / a.results).resolve()
    sd_files = sorted(results.glob('*_SOL*_M*_V*.sd')) + sorted(results.glob('*.EMConvData'))
    tag = a.tag or parent.name.split('_hfss_eigen_')[-1].rsplit('_', 1)[0]
    folder = create_run(DESIGN, 'resonator', f'eigen_from_sd_{tag}', entrypoint=__file__, inputs=sd_files,
                        upstream=[str(parent.relative_to(ROOT))],
                        parameters=dict(parent_run=parent.name, aedt_results=str(results.relative_to(ROOT)),
                                        band_GHz=BAND, q_min=Q_MIN))
    passes = read_passes(results)
    conv = {}
    for f in results.glob('*.EMConvData'):
        t = f.read_text(encoding='utf-8', errors='replace')
        conv = dict(is_converged=re.search(r'IsConverged=(\w+)', t).group(1),
                    max_delta_freq_pct={int(p): float(v) for p, v in
                                        re.findall(r'PassNo=(\d+).*?MaxDeltaFreq\(([\d.eE+-]+)', t, flags=re.S)},
                    solved_elements={int(p): int(v) for p, v in
                                     re.findall(r'PassNo=(\d+)\s*Intrinsics=\S*\s*SolvedElements\((\d+)', t)})
    analysis = {name: invert_pass(p['modes']) for name, p in passes.items()}
    out = dict(parent_run=parent.name, passes=passes, convergence=conv, analysis=analysis)
    (folder / 'summaries').mkdir(parents=True, exist_ok=True)
    (folder / 'summaries' / 'summary.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
    finish_run(folder, status='completed', summary=out)
    print('RUN_FOLDER=' + str(folder))
    for name in sorted(analysis):
        ch = analysis[name].get('channels', {})
        print(name, {q: {k: round(v, 5) for k, v in d.items() if isinstance(v, float)} for q, d in ch.items()},
              analysis[name].get('error', ''))


if __name__ == '__main__':
    main()
