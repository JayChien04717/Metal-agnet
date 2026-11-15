"""Compare signed Q3D couplings and the full three-by-three drive matrix."""
import argparse
from pathlib import Path
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
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
    base = HERE / 'results/design3e_e4_cmat_fF.csv'
    table = pd.DataFrame({'Design3e e4': parameters(base),
                          **{path.stem.removesuffix('_cmat_fF'): parameters(path) for path in args.matrix}})
    if args.output:
        table.to_csv(args.output, index_label='Parameter')
    print(table.to_string())


if __name__ == '__main__':
    main()
