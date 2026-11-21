"""Export an adaptive Q3D matrix from the saved Design3_FTF_NTHU project."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
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
        matrix.to_csv(_relocate(_P('Design3_FTF_NTHU/results/general/legacy') / f'{args.tag}_pass{step}_cmat_fF.csv'))


if __name__ == '__main__':
    main()
