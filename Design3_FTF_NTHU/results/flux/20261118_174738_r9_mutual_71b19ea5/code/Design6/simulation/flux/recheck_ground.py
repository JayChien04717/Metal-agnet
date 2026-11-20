"""Remeasure the confirmed r11 transmon flux fixture with recorded inputs."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code/general'))
from run_record import create_run, finish_run
import hfss_flux


def main():
    metadata = ROOT / 'Design6/design.json'
    candidate = json.loads(metadata.read_text(encoding='utf-8'))['transmon_flux_ground_candidate']
    flux = {'fbl_sep': f'{candidate["fbl_sep_um"]}um'}
    overrides = {'qflux': flux, 'coupler': {'flux_bias_line_options': dict(
        fbl_type='L', fbl_sep=flux['fbl_sep'], fbl_height='80um', fbl_length='45um',
        fbl_x='12um', cpw_width='10um', cpw_gap='6um', make_fbl=True)}}
    folder = create_run('Design6', 'flux', 'r11_mutual', entrypoint=__file__,
                        inputs=[metadata, ROOT / 'Design6' / candidate['reference_gds']],
                        parameters={'revision': candidate['revision'], 'overrides': overrides,
                                    'criterion_abs_M_pH': 1.0, 'max_passes': 6, 'mesh': '4um'})
    print('RUN_FOLDER', folder, flush=True)
    args = SimpleNamespace(overrides=json.dumps(overrides), routed=False, mesh='4um',
                           maxp=6, npts=10, tag='r11_mutual')
    try:
        hfss_flux.run(args, str(folder / 'raw'))
        result = hfss_flux.analyse(str(folder / 'raw'))
        result['above_1pH'] = {key: abs(value) > 1 for key, value in result['M_self_pH'].items()}
        result['convergence_evidence'] = 'Single six-pass seeded-mesh solve; pass-count stability pending'
        finish_run(folder, status='completed', validation='recorded', summary=result)
    except Exception:
        finish_run(folder, status='failed')
        raise


if __name__ == '__main__':
    main()
