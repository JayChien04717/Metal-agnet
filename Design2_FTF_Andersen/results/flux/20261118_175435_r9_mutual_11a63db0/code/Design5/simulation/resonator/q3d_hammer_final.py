"""Audit the delivered Design5 r5 with the full shared straight CPW span."""

from hashlib import sha256
import json
import traceback

import gdstk
import numpy as np

from q3d_hammer_check import ROOT, P, create_run, finish_run, solve


def check_mask(mask, hammer, limits_mm):
    """Verify both trace axes on the delivered mask across the modelled span."""
    from chip6_d4 import ORDER

    top, = gdstk.read_gds(str(mask)).top_level()
    metal = top.get_polygons()
    trace_um = 16
    for ch, x, up in ORDER:
        ys = (P.TAP_Y + np.linspace(limits_mm[0] + .005, limits_mm[1] - .005, 31)) * (
            1000 if up else -1000)
        for axis in (1000 * x, 1000 * x + hammer['cpw_centres_um']):
            points = [(axis + offset, y) for y in ys for offset in (-9, -7, 0, 7, 9)]
            assert gdstk.inside(points, metal) == (False, True, True, True, False) * len(ys), ch
        assert hammer['cpw_centres_um'] > trace_um


def main():
    import hfss_dirs
    import save_aedt

    home = ROOT / 'Design5'
    manifest = home / 'design.json'
    lengths = home / 'parameters/Design5_2026-11-13.json'
    config = json.loads(manifest.read_text(encoding='utf-8'))['wafer_candidate']
    mask = (home / config['delivery_fabrication_gds']).resolve()
    selected = config['hammer_geometry_um']
    canonical = home / config['fabrication_gds']
    digest = sha256(mask.read_bytes()).hexdigest()
    assert digest == config['delivery_gds_sha256'] == sha256(canonical.read_bytes()).hexdigest()
    radius_mm = float(P.d5.f3.FILLET.removesuffix('um')) / 1000
    limits = (P.WIDE['y_bottom'] + radius_mm - P.TAP_Y,
              P.TEE_HALF + P.WIDE['filter']['lead'] - radius_mm)
    check_mask(mask, selected, limits)
    prior = home / config['hammer_q3d_check']
    inputs = [manifest, lengths, mask, prior]
    folder = create_run('Design5', 'resonator', 'hammer_final_overlap_q3d', entrypoint=__file__,
                        inputs=inputs, upstream=[prior.parent.name],
                        parameters=dict(mask_sha256=digest, hammer=selected, cpw=P.CPW,
                                        axial_limits_mm=limits, straight_overlap_um=1000*(limits[1]-limits[0]),
                                        max_passes=8, percent_error=.3, freq_ghz=6.5,
                                        domain_mm=[2, 2, .5], eps_r=11.5))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        hfss_dirs.use_d_drive()
        save_aedt.use_project(save_aedt.project_file(home / 'simulation'))
        manifest_path = folder / 'run.json'
        record = json.loads(manifest_path.read_text(encoding='utf-8'))
        record['environment']['solver_version'] = save_aedt._desktop().get_version()
        record['geometry_revision'] = config['revision']
        manifest_path.write_text(json.dumps(record, indent=2), encoding='utf-8')
        original = {key: P.PAPER_HAMMER[key] for key in ('width_um', 'height_um', 'gap_um', 'ground_um')}
        base = solve(original, folder / 'raw/reference_400', 'D5FINAL_' + folder.name[-8:] + '_ref',
                     P.TEE_HALF, axial_limits_mm=limits)
        final = solve(selected, folder / 'raw/final_190', 'D5FINAL_' + folder.name[-8:] + '_r5',
                      P.TEE_HALF, axial_limits_mm=limits)
        compare = final['mutual_fF'] / base['mutual_fF'] - 1
        summary = dict(baseline=base, final=final, relative_delta=compare,
                       measured_parallel_span_um=1000*(limits[1]-limits[0]), mask_sha256=digest,
                       scope='Quasi-static open-ended coupon spanning the parallel straight CPWs. Meanders and bridges omitted.')
        save_aedt._desktop().get_active_project()._project.Save()
        finish_run(folder, status='completed', validation='final_layout_local_check', summary=summary)
        print(f'Original {base["mutual_fF"]:.5f} fF, final {final["mutual_fF"]:.5f} fF ({compare:+.2%})', flush=True)
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
