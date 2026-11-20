"""Record mutual inductance for both sense gaps of the r9 TC fixture."""
import json
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code/general'))
from run_record import create_run, finish_run
import recheck_readout_flux as runner


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', nargs='+', choices=('tc_left', 'tc_right'), default=['tc_left', 'tc_right'])
    parser.add_argument('--fixed-passes', type=int)
    args = parser.parse_args()
    metadata = ROOT / 'Design3/design.json'
    candidate = json.loads(metadata.read_text(encoding='utf-8'))['fluxonium_charge_ground_candidate']
    options = dict(candidate['layout_options'])
    options['tmon_options'] = {'flux_bias_line_options': dict(
        fbl_type='L', fbl_sep=f'{candidate["fbl_sep_um"]}um', fbl_height='80um',
        fbl_length='45um', fbl_x='12um', cpw_width='10um', cpw_gap='6um', make_fbl=True)}
    folder = create_run('Design3', 'flux', 'r9_mutual', entrypoint=__file__,
                        inputs=[metadata, ROOT / 'Design3' / candidate['reference_gds']],
                        parameters={'revision': candidate['revision'], 'layout_options': options,
                                    'criterion_abs_M_pH': 1.0, 'fixed_passes': args.fixed_passes})
    print('RUN_FOLDER', folder, flush=True)
    runner.PROJECT = ROOT / 'Design3/simulation/aedt/Design3.aedt'
    try:
        runner.base.save_aedt.use_project(runner.PROJECT)
        if args.fixed_passes:
            runner.base.mutual.SETUP.update(max_passes=args.fixed_passes,
                                           min_passes=args.fixed_passes, min_converged=1)
        for case in args.case:
            runner.run_mutual(case, layout_options=options, output=folder / 'raw')
            runner.base.save_aedt._desktop().get_active_project()._project.Save()
        summaries = {case: json.loads((folder / 'raw' / f'm_{case}_summary.json').read_text())
                     for case in args.case}
        runner.base.save_aedt.snapshot(folder / 'raw', 'r9_mutual')
        finish_run(folder, status='completed', validation='recorded', summary=summaries)
    except Exception:
        finish_run(folder, status='failed')
        raise


if __name__ == '__main__':
    main()
