"""Keep the Git (GitHub) copy of the repository small: track code, parameters and the
records of the runs each Design actually uses; everything else stays local.

Policy (see skill/parameter_provenance.md section 7):
* `.gitignore` excludes run folders (`Design*/results/*/*/`) and generated GDS
  (`Design*/gds/**/*.gds`); the layout code regenerates the GDS.
* Runs referenced by a Design's `design.json` (any key except historical / failed / removed
  lists) are force-tracked, but only their records: `run.json`, `summaries/`, `figures/`, the
  top-level `*.json` / `*.md` / `*.txt` / `*.png` files, and any file `design.json` points to.
  `code/`, `inputs/`, `raw/`, `logs/` snapshots stay local.

Run after changing `adopted_results` / candidates in a `design.json` (or after editing
`.gitignore`). Default is a dry run; `--apply` updates the Git index only (local files are
never deleted).

Usage (repository root):
    python -I -X utf8 -B code/general/sync_tracked_results.py            # dry run
    python -I -X utf8 -B code/general/sync_tracked_results.py --apply
"""

import argparse
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUN_ID = re.compile(r'\d{8}_\d{6}_[A-Za-z0-9_\-]+')
SKIP_KEYS = re.compile(r'historical|failed|removed', re.I)
RECORD_DIRS = ('summaries', 'figures')
RECORD_SUFFIXES = ('.json', '.md', '.txt', '.png')


def git(*args, stdin=None):
    out = subprocess.run(['git', '-c', 'core.quotepath=off', *args], cwd=ROOT, input=stdin,
                         capture_output=True, check=True)
    return out.stdout.decode('utf-8')


def referenced_runs(design_json):
    """Run IDs referenced by design.json, except under historical / failed / removed keys."""
    d = json.loads(design_json.read_text(encoding='utf-8'))
    text = json.dumps({k: v for k, v in d.items() if not SKIP_KEYS.search(k)})
    paths = set(re.findall(r'"((?:[^"]*/)?results/[^"]+)"', text))
    return set(RUN_ID.findall(text)), paths


def keep_files():
    """{repo-relative posix path} of run records that must stay tracked; plus missing runs."""
    keep, missing = set(), []
    for dj in sorted(ROOT.glob('Design*/design.json')):
        design = dj.parent
        ids, paths = referenced_runs(dj)
        for rid in sorted(ids):
            runs = [p for p in design.glob(f'results/*/{rid}') if p.is_dir()]
            if not runs:
                missing.append(f'{design.name}: {rid}')
            for run in runs:
                for f in run.iterdir():
                    if f.is_file() and f.suffix.lower() in RECORD_SUFFIXES:
                        keep.add(f)
                for sub in RECORD_DIRS:
                    keep.update(p for p in (run / sub).rglob('*') if p.is_file())
        for p in paths:
            f = design / p[p.index('results/'):]
            if f.is_file():
                keep.add(f)
    return {f.relative_to(ROOT).as_posix() for f in keep}, missing


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--apply', action='store_true', help='update the Git index (default: dry run)')
    args = ap.parse_args()

    keep, missing = keep_files()
    tracked = set(git('ls-files', '-z').split('\0')) - {''}
    ignored_tracked = set(git('ls-files', '-z', '-ci', '--exclude-standard').split('\0')) - {''}
    untrack = sorted(ignored_tracked - keep)
    add = sorted(keep - tracked)

    print(f'run records to keep tracked: {len(keep)} ({len(add)} not yet tracked)')
    print(f'tracked files to untrack (ignored, not run records): {len(untrack)}')
    for m in missing:
        print(f'WARNING referenced run folder not found locally: {m}')
    if not args.apply:
        print('dry run; pass --apply to update the index')
        return
    if untrack:
        git('rm', '--cached', '-q', '--pathspec-from-file=-', '--pathspec-file-nul',
            stdin='\0'.join(untrack).encode('utf-8'))
    if add:
        git('add', '-f', '--pathspec-from-file=-', '--pathspec-file-nul',
            stdin='\0'.join(add).encode('utf-8'))
    print('index updated; review with `git status` before committing')


if __name__ == '__main__':
    main()
