"""Delete old designs (and their solver results) from a Design's shared AEDT project.

    python -X utf8 -B code/aedt/prune_aedt_designs.py <project.aedt> --keep NAME [--keep NAME ...] [--apply]

Without --apply it only lists what would be deleted (dry run). With --apply it opens the project in AEDT, deletes every
design not named by --keep (oProject.DeleteDesign), saves, and reports the .aedtresults size before / after.
Run only when no solver run is using the project (no other AEDT session holding it, no <project>.aedt.lock).
Extracted run data live in results/<category>/<run>/ and are not touched; deleted designs cannot be reopened.
The AEDT session started here is closed by code/aedt/aedt_session.py (imported via project_paths).
"""
from pathlib import Path
import argparse
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'code' / 'aedt')]
import project_paths
project_paths.configure(__file__)


def _size(folder):
    return sum(f.stat().st_size for f in Path(folder).rglob('*') if f.is_file()) if Path(folder).exists() else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('project')
    ap.add_argument('--keep', action='append', default=[], help='design name to keep (repeatable)')
    ap.add_argument('--apply', action='store_true')
    a = ap.parse_args()
    path = Path(a.project).resolve()
    results = path.with_suffix('.aedtresults')
    lock = Path(str(path) + '.lock')
    if lock.exists():
        raise SystemExit(f'{lock} exists: the project is open in AEDT (a run may be using it). Close it first.')

    import save_aedt
    save_aedt.use_project(str(path))
    desk = save_aedt._desktop()
    proj = desk._desktop.GetActiveProject()
    names = [n.split(';')[-1] for n in proj.GetTopDesignList()]
    missing = [k for k in a.keep if k not in names]
    if missing:
        raise SystemExit(f'--keep names not in the project: {missing} (designs: {names})')
    drop = [n for n in names if n not in a.keep]
    before = _size(results)
    print(f'{path.name}: {len(names)} designs, .aedtresults {before / 1e9:.2f} GB')
    print('keep  :', a.keep)
    print('delete:', drop)
    if not a.apply:
        print('dry run: nothing deleted (add --apply)')
        return
    for n in drop:
        proj.DeleteDesign(n)
    proj.Save()
    after = _size(results)
    print(f'deleted {len(drop)} designs; .aedtresults {before / 1e9:.2f} -> {after / 1e9:.2f} GB')


if __name__ == '__main__':
    main()
