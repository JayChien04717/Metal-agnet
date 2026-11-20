"""Core Q3D (Q1, TC, Q2 + readout claws + XY lines) on the current r9 layout, written to a new run folder.

Same rendering and analysis as the legacy entry point q3d_core.py (run q3d_v8, 2026-11-14, revision A):
  components Q1, Q2, TC, xy_Q1, xy_Q2, xy_end_Q1, xy_end_Q2; opens at the readout claw pins and the XY route starts;
  coupler flux line drawn as ground; SQUID electrodes flush with ground (simulation geometry); channels=False (claws keep
  their 40 um CPW stub, no Purcell channel); analysis = q3d_core.analyse (Ec, g_n, readout g, XY T1 / crosstalk).
Differences from q3d_v8: build_design defaults of the current layout (B2 flux lines, updated launch pads) and
wafer_ports=True (r9 side launches, so the XY routes follow the r9 paths).
Solver: 20 passes / 0.1 % (the Design6 readout_caps setting) instead of v8's 8 passes / 0.5 %; every pass is saved, so
pass 8 gives the like-for-like comparison with v8 and the later passes give the convergence evidence.
"""
from argparse import ArgumentParser
from pathlib import Path
import json
import sys
import time
import traceback

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(Path(__file__).parent), str(ROOT / 'code/general')]
import q3d_core as QC                          # noqa: E402  (sets the layout / note paths)
from run_record import create_run, finish_run  # noqa: E402

L = QC.L
V8 = ROOT / 'Design6/results/qubit/2026-11-14/q3d_v8'
COMPS = ['Q1', 'Q2', 'TC', 'xy_Q1', 'xy_Q2', 'xy_end_Q1', 'xy_end_Q2']
OPENS = [('Q1', 'readout'), ('Q2', 'readout'), ('xy_Q1', 'start'), ('xy_Q2', 'start')]
REVISION = '2026-11-16_wafer_r9'
SETUP = dict(min_passes=2, freq_ghz=5.0, solution_order='Medium', auto_increase_solution_order=False)


def build():
    cp = dict(flux_as_ground='True', squid_contact_overlap='0um', squid_shift_toward_stem='0um')
    return L.build_design(channels=False, wafer_ports=True, coupler=cp), cp


def brief(res):
    return dict(Ec_GHz=res['Ec_GHz'], C_sigma_fF=res['C_sigma_fF'], g_n_MHz=res['g_n_MHz'],
                C_mutual_fF=res['C_mutual_fF'],
                readout={k: {kk: v for kk, v in r.items() if isinstance(v, (int, float))} for k, r in res['readout'].items()},
                charge_lines=res['charge_lines'],
                xy_crosstalk_pct={k: v['crosstalk_pct'] for k, v in res['xy_crosstalk'].items()})


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--max-passes', type=int, default=20)
    ap.add_argument('--percent-error', type=float, default=0.1)
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()
    d, cp = build()
    names = {c: list(d.components[c].pins.keys()) for c in COMPS}
    open_pins = [
        dict(component=component, pin=pin,
             middle_mm=d.components[component].pins[pin]['middle'].tolist(),
             normal=d.components[component].pins[pin]['normal'].tolist())
        for component, pin in OPENS
    ]
    setup = dict(SETUP, max_passes=a.max_passes, percent_error=a.percent_error)
    folder = create_run('Design6', 'qubit', 'q3d_core_r9', entrypoint=__file__,
                        inputs=[ROOT / 'Design6/design.json', ROOT / 'Design6/layout/design6_layout.py',
                                V8 / 'summary.json', V8 / 'cmat_fF.csv', V8 / 'convergence.csv'],
                        upstream=[str(V8)],
                        parameters=dict(geometry='current layout defaults (B2 flux lines) + wafer_ports=True (r9), '
                                                 'channels=False', coupler_overrides=cp, components=COMPS,
                                        open_terminations=open_pins, pins=names, F_Q_GHz={k: v / 1e9 for k, v in QC.F_Q.items()},
                                        F_RO_GHz={k: v / 1e9 for k, v in QC.F_RO.items()},
                                        Q3D=setup, ground_straps='none (legacy core fixture)',
                                        resonator_model='claw plus 40 um stub; lambda/4 estimate at target frequency',
                                        comparison='q3d_v8 (rev-A, 8 passes / 0.5 %)', dry_run=a.dry_run))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        manifest = folder / 'run.json'
        record = json.loads(manifest.read_text(encoding='utf-8'))
        record['geometry_revision'] = REVISION
        manifest.write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        if a.dry_run:
            finish_run(folder, status='completed', summary=dict(dry_run=True, components=names, open_pins=open_pins))
            print(json.dumps(open_pins, indent=2))
            return
        from qiskit_metal.analyses.quantization import LOManalysis
        import save_aedt, hfss_dirs
        temp = hfss_dirs.use_d_drive()
        project = save_aedt.project_file(ROOT / 'Design6/simulation')
        save_aedt.use_project(project)
        record['environment']['solver_version'] = save_aedt._desktop().get_version()
        record['environment']['solver_temp'] = temp['after']
        record['parameters']['aedt_project'] = str(project)
        manifest.write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        lom = LOManalysis(d, 'q3d')
        lom.sim.setup.update(setup)
        t0 = time.time()
        lom.sim.run(name='D6CORE_' + folder.name[-8:], components=COMPS, open_terminations=OPENS)
        cmat = lom.sim.capacitance_matrix
        if str(lom.sim.units).lower() not in ('ff', 'femtofarad', 'femtofarads'):
            raise RuntimeError(f'unexpected units {lom.sim.units}')
        cmat.to_csv(folder / 'raw/cmat_fF.csv')
        passes = {}
        for k, m in lom.sim.capacitance_all_passes.items():
            f = (m if isinstance(m, pd.DataFrame) else pd.DataFrame(m, index=cmat.index, columns=cmat.columns))
            # Qiskit Metal's Q3D renderer returns the adaptive matrices in farads.
            f = f * 1e15
            f.to_csv(folder / f'raw/pass_{k}_fF.csv')
            try:
                passes[int(k)] = brief(QC.analyse(f, verbose=False))
            except Exception as ex:                                       # early passes may miss small entries
                passes[int(k)] = dict(error=repr(ex)[:200])
        res = QC.analyse(cmat)
        convergence_error = None
        try:
            convergence, convergence_report = lom.sim.renderer.pinfo.setup.get_convergence()
            if convergence_report:
                (folder / 'raw/convergence.txt').write_text(convergence_report, encoding='utf-8')
            if convergence is not None:
                convergence.to_csv(folder / 'raw/convergence.csv')
        except Exception as ex:
            convergence_error = repr(ex)
        ks = sorted(k for k in passes if 'Ec_GHz' in passes[k])
        conv = {}
        if len(ks) > 1:
            last, prev = passes[ks[-1]], passes[ks[-2]]
            conv = dict(Ec_rel_change={q: abs(last['Ec_GHz'][q] / prev['Ec_GHz'][q] - 1) for q in last['Ec_GHz']},
                        g_rel_change={q: abs(last['g_n_MHz'][q] / prev['g_n_MHz'][q] - 1) for q in last['g_n_MHz']})
        v8 = json.loads((folder / 'inputs/2/summary.json').read_text(encoding='utf-8'))
        cmp = {}
        for key in ('Ec_GHz', 'C_sigma_fF', 'g_n_MHz', 'C_mutual_fF'):
            cmp[key] = {q: dict(r9_final=res[key][q], r9_pass8=passes.get(8, {}).get(key, {}).get(q), v8=v8[key][q],
                                final_minus_v8=res[key][q] - v8[key][q]) for q in res[key]}
        out = dict(final=brief(res), net_mapping=res['nodes'],
                   solver_converged=lom.sim.is_converged, convergence_export_error=convergence_error,
                   passes_saved=len(passes),
                   last_pass_change=conv, by_pass=passes, vs_v8=cmp, minutes=(time.time() - t0) / 60,
                   aedt_design=lom.sim.renderer.pinfo.design.name)
        lom.sim.renderer.pinfo.project._project.Save()
        (folder / 'summary_out.json').write_text(json.dumps(out, indent=1, default=float))
        finish_run(folder, status='completed', summary=dict(final=out['final'], solver_converged=out['solver_converged'],
                                                            last_pass_change=conv, vs_v8=cmp))
        print(json.dumps(dict(vs_v8=cmp, converged=out['solver_converged'], last_pass_change=conv), indent=1, default=float))
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()}); raise


if __name__ == '__main__':
    main()
