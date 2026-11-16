"""Compare the 190 µm straight readout hammer on Design6 with its 400 µm baseline."""

from argparse import Namespace
import json
from pathlib import Path
import sys
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / path) for path in ('code/general', 'Design6/layout',
                                           'Design5/simulation/resonator')]
import design6_layout as layout
import hfss_readout
from q3d_hammer_check import solve
from run_record import create_run, finish_run


def main():
    import hfss_dirs
    import save_aedt

    home = ROOT / 'Design6'
    filename = home / 'parameters/Design6_2026-11-14.json'
    channels = json.loads(filename.read_text(encoding='utf-8'))['channels']
    lengths = {ch['name']: dict(L_filter=ch['L_filter_mm'], L_res=ch['L_res_mm']) for ch in channels}
    config_file = home / 'design.json'
    config = json.loads(config_file.read_text(encoding='utf-8'))['wafer_candidate']
    previous_gds = (home / config['fabrication_gds']).resolve()
    baseline = {k: layout.HAMMER[k] for k in ('width_um', 'height_um', 'gap_um', 'ground_um')}
    target = dict(baseline, cpw_centres_um=190, height_um=40, gap_um=40)
    preview = layout.build_design(lengths=lengths, readout_hammer=target, wafer_ports=True)
    channel_info = preview.metadata['channels']
    limits = tuple(float(x) for x in np.mean(
        [channel['parallel_limits_from_tap_mm'] for channel in channel_info.values()], axis=0))
    # Reuse the validated 2 × 2 mm coupon domain. The selected overlap fits inside it.
    assert -.99 < limits[0] < -.15 and .15 < limits[1] < .99
    for info in channel_info.values():
        assert abs(abs(info['jr_x'] - info['jf_x']) - .190) < 1e-10
    folder = create_run('Design6', 'resonator', 'parallel_hammer_quick', entrypoint=__file__,
                        inputs=[config_file, filename, previous_gds],
                        parameters=dict(original=baseline, candidate_seed=target, q3d_limits_mm=limits,
                                        q3d_overlap_by_channel_mm={q: info['parallel_limits_from_tap_mm']
                                                                   for q, info in channel_info.items()},
                                        q3d_passes=8, q3d_pct=.3, q3d_eps_r=11.5,
                                        hfss_passes=3, hfss_mesh='25um', lengths_mm=lengths,
                                        coupling_comparison_tolerance=.05, frequency_comparison_tolerance=.01))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(home / 'simulation'))
        manifest_path = folder / 'run.json'
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        manifest['environment']['solver_version'] = save_aedt._desktop().get_version()
        manifest['geometry_revision'] = 'Design6 B2 straight 190 um hammer candidate'
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
        ident = folder.name[-8:]
        ref = solve(baseline, folder / 'raw/q3d_reference', f'D6CJ_{ident}_ref', layout.P.TEE_HALF,
                    axial_limits_mm=limits)
        goal = ref['mutual_fF']
        candidates = []
        for k in range(3):
            row = solve(target, folder / f'raw/q3d_candidate_{k}', f'D6CJ_{ident}_c{k}', layout.P.TEE_HALF,
                        axial_limits_mm=limits)
            row['raw_path'] = f'raw/q3d_candidate_{k}/caps.json'
            candidates.append(row)
            if abs(row['mutual_fF'] / goal - 1) <= .05:
                break
            height = target['height_um']
            height_next = height * goal / row['mutual_fF']
            if len(candidates) > 1:
                earlier = candidates[-2]
                slope = (row['mutual_fF'] - earlier['mutual_fF']) / (height - earlier['geometry']['height_um'])
                if slope > 0:
                    height_next = height + (goal - row['mutual_fF']) / slope
            height_next = round(float(np.clip(height_next, 15, 75)), 1)
            if height_next == height:
                break
            target = dict(target, height_um=height_next)
        chosen = min(candidates, key=lambda item: abs(item['mutual_fF'] - goal))
        print('SELECTED_HAMMER=' + json.dumps(chosen['geometry']), flush=True)
        quick = {}
        for q in ('Q1', 'Q2'):
            modes = {}
            for key, hammer in (('reference', baseline), ('parallel', chosen['geometry'])):
                output = folder / 'raw' / f'{q}_{key}'
                output.mkdir()
                args = Namespace(qubit=q, trim=0, tag=key, mesh='25um', maxp=3, quick=True,
                                 lengths=lengths, readout_hammer=hammer,
                                 design_name=f'D6RO_{ident}_{q}_{key}')
                modes[key] = hfss_readout.run(args, output)
            f0, f1 = (np.sort(modes[key]['f_GHz']) for key in ('reference', 'parallel'))
            identified = len(f0) == len(f1) == 3 and f0[1] - f0[0] > .7 and f1[1] - f1[0] > .7
            modes['readout_pair_candidate_GHz'] = dict(reference=f0[1:].tolist(), parallel=f1[1:].tolist()) if identified else None
            modes['readout_pair_candidate_shift_MHz'] = ((f1[1:] - f0[1:])*1000).tolist() if identified else None
            modes['mode_identification'] = 'lowest mode separated from upper doublet; no field participation check' if identified else 'ambiguous'
            quick[q] = modes
        summary = dict(q3d_reference=ref, q3d_selected=chosen, q3d_cases=candidates,
                       q3d_relative_delta=chosen['mutual_fF']/goal-1, eigenmode=quick,
                       q3d_scope='Local coupon at mean shared straight CPW span of the two channels',
                       eigenmode_scope='Three adaptive passes per loaded readout channel; no full-chip validation')
        save_aedt._desktop().get_active_project()._project.Save()
        finish_run(folder, status='completed', validation='rough_readout_comparison', summary=summary)
        print(json.dumps(summary, indent=2), flush=True)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
