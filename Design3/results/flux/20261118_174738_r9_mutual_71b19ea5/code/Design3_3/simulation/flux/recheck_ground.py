"""Record mutual inductance for both sense gaps of the r9 TC fixture."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code/general'))
from run_record import create_run, finish_run
import recheck_readout_flux as runner


def main():
    metadata = ROOT / 'Design3_3/design.json'
    candidate = json.loads(metadata.read_text(encoding='utf-8'))['fluxonium_charge_ground_candidate']
    options = dict(candidate['layout_options'])
    options['tmon_options'] = {'flux_bias_line_options': dict(
        fbl_type='L', fbl_sep=f'{candidate["fbl_sep_um"]}um', fbl_height='80um',
        fbl_length='45um', fbl_x='12um', cpw_width='10um', cpw_gap='6um', make_fbl=True)}
    folder = create_run('Design3_3', 'flux', 'r9_mutual', entrypoint=__file__,
                        inputs=[metadata, ROOT / 'Design3_3' / candidate['reference_gds']],
                        parameters={'revision': candidate['revision'], 'layout_options': options,
                                    'criterion_abs_M_pH': 1.0})
    print('RUN_FOLDER', folder, flush=True)
    runner.PROJECT = ROOT / 'Design3_3/simulation/aedt/Design3_3.aedt'
    try:
        runner.base.save_aedt.use_project(runner.PROJECT)
        for case in ('tc_left', 'tc_right'):
            runner.run_mutual(case, layout_options=options, output=folder / 'raw')
            runner.base.save_aedt._desktop().get_active_project()._project.Save()
        summaries = {case: json.loads((folder / 'raw' / f'm_{case}_summary.json').read_text())
                     for case in ('tc_left', 'tc_right')}
        runner.base.save_aedt.snapshot(folder / 'raw', 'r9_mutual')
        finish_run(folder, status='completed', validation='recorded', summary=summaries)
    except Exception:
        finish_run(folder, status='failed')
        raise


if __name__ == '__main__':
    main()
