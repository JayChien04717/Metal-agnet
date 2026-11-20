"""Record both TC sense-gap mutual inductances for the confirmed r9 geometry."""
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'code/general'))
from run_record import create_run, finish_run


def main():
    metadata = ROOT / 'Design2/design.json'
    candidate = json.loads(metadata.read_text(encoding='utf-8'))['transmon_flux_ground_candidate']
    folder = create_run('Design2', 'flux', 'r9_mutual', entrypoint=__file__,
                        inputs=[metadata, ROOT / 'Design2' / candidate['reference_gds']],
                        parameters={'revision': candidate['revision'],
                                    'tc_flux_offset_um': candidate['tc_flux_offset_um'],
                                    'criterion_abs_M_pH': 1.0})
    print('RUN_FOLDER', folder, flush=True)
    env = dict(os.environ, MUTUAL_OUTPUT_DIR=str(folder / 'raw'),
               TC_FLUX_OFFSET_UM=str(candidate['tc_flux_offset_um']))
    script = Path(__file__).with_name('run_mutual.py')
    try:
        subprocess.run([sys.executable, '-B', '-X', 'utf8', str(script), 'tc_left', 'tc_right'],
                       env=env, check=True)
        import pandas as pd
        result = {}
        for case in ('tc_left', 'tc_right'):
            table = pd.read_csv(folder / 'raw' / f'D2_M_{case}.csv')
            if min(table.Z11_L_nH.iloc[0], table.Z22_L_nH.iloc[0]) <= 0:
                raise ValueError(f'{case}: low-frequency port self-inductance is nonpositive')
            result[case] = {'M_pH_fmin': float(table.M_pH.iloc[0]),
                            'above_1pH': bool(abs(table.M_pH.iloc[0]) > 1),
                            'max_reciprocity_error': float(table.reciprocity.max())}
        import save_aedt
        save_aedt.snapshot(folder / 'raw', 'r9_mutual')
        finish_run(folder, status='completed', validation='recorded', summary=result)
    except Exception:
        finish_run(folder, status='failed')
        raise


if __name__ == '__main__':
    main()
