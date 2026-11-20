"""Estimate bare transitions using the Design2 circuit-model inputs."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import json
from pathlib import Path
import sys

import numpy as np

RESULTS = _P('Design2/results/general/legacy/layout3_full')
FX_EJ_GHZ = 4.0
FX_EL_GHZ = 0.6
FX_FLUX = 0.5
TC_EJ_OVER_EC = 50


def main():
    sys.path.insert(0, _I('Design3/simulation'))
    from zz_analysis import fluxonium, transmon

    electrical = json.loads((_P('Design2/results/general/legacy/layout3_full/electrical.json')).read_text(encoding='utf-8'))
    fx_names = ('FX1', 'FX2')
    devices = {
        name: fluxonium(electrical[f'Ec_{name} [GHz]'], FX_EJ_GHZ, FX_EL_GHZ,
                        phi_ext=2 * np.pi * FX_FLUX, n_keep=2)
        for name in fx_names
    }
    ec_tc = electrical['Ec_TC [GHz]']
    devices['TC'] = transmon(ec_tc, TC_EJ_OVER_EC * ec_tc, n_keep=2)
    order = ('FX1', 'TC', 'FX2')
    n01 = {name: float(abs(devices[name][1][0, 1])) for name in order}
    result = {
        'FX_EJ_GHz': FX_EJ_GHZ,
        'FX_EL_GHz': FX_EL_GHZ,
        'FX_flux_in_Phi0': FX_FLUX,
        'TC_EJ_over_EC': TC_EJ_OVER_EC,
        'f01_GHz': {name: float(devices[name][0][1]) for name in order},
        'n01': n01,
        'transition_coupling_MHz': {
            f'{a}_{b}': abs(electrical[j]) * n01[a] * n01[b]
            for a, b, j in (
                ('FX1', 'TC', 'J_1c [MHz]'),
                ('FX2', 'TC', 'J_2c [MHz]'),
                ('FX1', 'FX2', 'J_12 [MHz]'),
            )
        },
    }
    (_P('Design2/results/general/legacy/layout3_full/transition_estimates.json')).write_text(json.dumps(result, indent=2) + '\n',
                                                       encoding='utf-8')


if __name__ == '__main__':
    main()
