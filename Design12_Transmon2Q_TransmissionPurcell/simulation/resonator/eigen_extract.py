# -*- coding: utf-8 -*-
r"""Read a saved HFSS eigenmode solution from simulation/aedt/D12.aedt into a new run (re-analysis of a parent run).

Why: on this machine the COM call that waits on Analyze can raise com_error -2147352567 part way through an eigenmode
solve while HFSS keeps solving and later reports normal completion (Design12 r4, 2026-11-26). The parent run is then
marked failed although the solution exists. This script opens the project, reads the convergence table and the
complex eigenfrequencies of --design, classifies filter / resonator modes (hfss_purcell_eigen.classify), and stores
them in a new run whose upstream is the parent run.

Usage (repository root, qmetal):
    python -X utf8 -B Design12_Transmon2Q_TransmissionPurcell/simulation/resonator/eigen_extract.py \
        --parent Design12_Transmon2Q_TransmissionPurcell/results/resonator/<run> --design E12_xxxxxxxx_hfss \
        --tag r4_4p_crop --f-pf 6.24
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code' / 'general'))
import project_paths
project_paths.configure(__file__)
DESIGN = project_paths.design_dir(__file__)
sys.path.insert(0, str(DESIGN / 'layout'))
sys.path.insert(0, str(Path(__file__).resolve().parent))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--parent', required=True)
    ap.add_argument('--design', required=True, help='AEDT design name inside D12.aedt')
    ap.add_argument('--tag', required=True)
    ap.add_argument('--f-pf', type=float, default=6.24, help='expected filter frequency for mode classification, GHz')
    a = ap.parse_args()
    from run_record import create_run, finish_run
    from pyEPR import ProjectInfo
    import build as B
    from hfss_purcell_eigen import classify

    parent = Path(a.parent).resolve()
    folder = create_run(DESIGN.name, 'resonator', f'eigen_extract_{a.tag}', entrypoint=__file__,
                        upstream=[str(parent.relative_to(ROOT))],
                        parameters=dict(parent_run=parent.name, aedt_design=a.design, f_pf_classify_GHz=a.f_pf,
                                        reason='parent run stopped on a COM error at Analyze; read the saved solution'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        pinfo = ProjectInfo(project_path=str(DESIGN / 'simulation' / 'aedt'), project_name='D12',
                            design_name=a.design)
        setup = pinfo.design.get_setup('Setup')
        conv = setup.get_convergence()[0]
        (folder / 'raw').mkdir(parents=True, exist_ok=True)
        conv.to_csv(folder / 'raw' / 'eigen_convergence.csv')
        f, im = setup.get_solutions().eigenmodes()
        raw = dict(frequencies_GHz=list(map(float, f)), imag_GHz=list(map(float, im)), passes=int(len(conv)),
                   last_delta_f_pct=float(conv.iloc[-1, -1]) if len(conv) else None, aedt_design=a.design)
        (folder / 'raw' / 'eigen_modes.json').write_text(json.dumps(raw, indent=2), encoding='utf-8')
        targets = {q: B.RR_CHANNELS[q]['f_GHz'] for q in ('Q1', 'Q2')}
        modes = classify(raw, targets, a.f_pf)
        summary = dict(parent_run=parent.name, hfss_eigen=raw, modes=modes)
        (folder / 'summaries').mkdir(parents=True, exist_ok=True)
        (folder / 'summaries' / 'summary.json').write_text(json.dumps(summary, indent=2, default=float),
                                                          encoding='utf-8')
        finish_run(folder, status='completed', summary=summary)
        print(conv.to_string(), flush=True)
        print(json.dumps(modes, indent=1, default=float), flush=True)
    except Exception as exc:
        finish_run(folder, status='failed', summary={'error': repr(exc)})
        raise


if __name__ == '__main__':
    main()
