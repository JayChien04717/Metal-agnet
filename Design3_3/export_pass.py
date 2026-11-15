"""Export an adaptive Q3D matrix from the saved Design3_3 project."""
import argparse
from pathlib import Path

from q3d import HERE, base


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tag')
    parser.add_argument('passes', type=int, nargs='+')
    parser.add_argument('--project', type=Path,
                        default=Path.home() / 'AppData/Local/Temp/D33Tune/D33Tune.aedt')
    args = parser.parse_args()
    base.save_aedt.use_project(args.project)
    project = base.save_aedt._desktop().get_active_project()
    setup = project.get_design(f'{args.tag}_q3d').get_setup('Setup')
    for step in args.passes:
        matrix, *_ = setup.get_matrix(pass_number=step, solution_kind='AdaptivePass')
        matrix.to_csv(HERE / 'results' / f'{args.tag}_pass{step}_cmat_fF.csv')


if __name__ == '__main__':
    main()
