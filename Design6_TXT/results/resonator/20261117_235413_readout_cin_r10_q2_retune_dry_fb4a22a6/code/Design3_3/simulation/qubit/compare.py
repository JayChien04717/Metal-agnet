"""Compare signed Q3D couplings and the full three-by-three drive matrix."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse
from pathlib import Path
import sys

import pandas as pd

HERE = _P('Design3_3')
sys.path.insert(0, _I('Design3_3'))
from validate_gate import extracted


def parameters(path):
    result, beta = extracted(path)
    matrix = pd.read_csv(path, index_col=0)
    electrodes = {
        'FX1': ['pad_bot_FX1', 'pad_top_FX1'],
        'TC': ['tmon_TC'],
        'FX2': ['pad_bot_FX2', 'pad_top_FX2'],
    }
    for i, (target, pads) in enumerate(electrodes.items()):
        for j, line in enumerate(electrodes):
            result[f'beta {line} charge -> {target}'] = beta[i, j]
            result[f'C {line} charge -> {target} [fF]'] = -matrix.loc[f'charge_line_{line}', pads].sum()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('matrix', type=Path, nargs='+')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    base = _P('Design3_3/results/general/legacy/design3e_e4_cmat_fF.csv')
    table = pd.DataFrame({'Design3e e4': parameters(base),
                          **{path.stem.removesuffix('_cmat_fF'): parameters(path) for path in args.matrix}})
    if args.output:
        table.to_csv(args.output, index_label='Parameter')
    print(table.to_string())


if __name__ == '__main__':
    main()
