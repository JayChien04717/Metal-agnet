"""Minimal run records: unique folder, input/source snapshots and explicit provenance.

Use create_run() before solving, write outputs into its returned directory, then
finish_run(). Existing legacy scripts are not silently redirected by this module.
"""
from datetime import datetime
from hashlib import sha256
from importlib.metadata import distributions
from pathlib import Path
import argparse
import json
import re
import subprocess
import sys
from uuid import uuid4

ROOT = Path(__file__).resolve().parents[2]
CATEGORIES = ('general', 'qubit', 'resonator', 'flux')


def _git(*args):
    result = subprocess.run(['git', '-C', str(ROOT), *args], capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else None


def create_run(design, category, tag, *, entrypoint, inputs=(), upstream=(), parameters=None):
    """Snapshot active project Python sources conservatively, including cross-Design dependencies.

    Historical source copies under results are excluded. Solver files and GDS can
    be passed explicitly as inputs. The caller owns actual output routing.
    """
    if category not in CATEGORIES or not re.fullmatch(r'Design\d+(?:_\d+)?', design):
        raise ValueError('Expected Design<N> and general/qubit/resonator/flux')
    design_root = ROOT / design
    if not design_root.is_dir():
        raise ValueError(f'Unknown design: {design}')
    entry = Path(entrypoint).resolve()
    if not entry.is_relative_to(ROOT) or not entry.is_file():
        raise ValueError('entrypoint must be an existing project file')
    # Validate all input files before creating a run.
    inputs = [Path(p).resolve() for p in inputs]
    for source in inputs:
        if not source.is_file():
            raise FileNotFoundError(source)
    now = datetime.now().astimezone()
    safe_tag = re.sub(r'[^\w.-]+', '_', tag).strip('._') or 'run'
    run_id = now.strftime('%Y%m%d_%H%M%S') + '_' + safe_tag + '_' + uuid4().hex[:8]
    folder = design_root / 'results' / category / run_id
    folder.mkdir(parents=True, exist_ok=False)
    (folder / 'raw').mkdir()
    git_status = _git('status', '--porcelain')
    record = dict(run_id=run_id, design=design, category=category, started_at=now.isoformat(),
                  status='prepared', validation_status='unknown', geometry_revision=None,
                  entrypoint=entry.relative_to(ROOT).as_posix(), command=sys.argv,
                  working_directory=str(Path.cwd()), parameters=parameters or {},
                  upstream=list(upstream), git_commit=_git('rev-parse', 'HEAD'),
                  git_dirty=None if git_status is None else bool(git_status), inputs=[], code=[],
                  environment={'python':sys.version, 'executable':sys.executable,
                               'packages':{d.metadata['Name']:d.version for d in distributions()
                                           if d.metadata.get('Name')}, 'solver_version':None},
                  source_scope='All active project Python files; snapshots are current, not historical run evidence')
    try:
        for index, source in enumerate(inputs):
            target = folder / 'inputs' / str(index) / source.name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source.read_bytes())
            record['inputs'].append({'source':str(source), 'snapshot':target.relative_to(folder).as_posix(),
                                     'sha256':sha256(target.read_bytes()).hexdigest()})
        for base in [ROOT/'code', ROOT/'pattern design'] + sorted(ROOT.glob('Design*')):
            for source in base.rglob('*.py'):
                relative = source.relative_to(ROOT)
                if any(p in ('results', '__pycache__', '.migration-work') for p in relative.parts):
                    continue
                target = folder / 'code' / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(source.read_bytes())
                record['code'].append({'path':relative.as_posix(), 'sha256':sha256(target.read_bytes()).hexdigest()})
        # Snapshot the relocation index used by the Python sources as well.
        source = ROOT/'code/general/path_migration.json'
        target = folder/'code/code/general/path_migration.json'
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(source.read_bytes())
        record['code'].append({'path':'code/general/path_migration.json', 'sha256':sha256(target.read_bytes()).hexdigest()})
    except Exception:
        record['status'] = 'preparation_failed'
        (folder/'run.json').write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        raise
    (folder/'run.json').write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    return folder


def finish_run(folder, *, status, validation='unknown', summary=None):
    folder = Path(folder).resolve()
    if not folder.is_relative_to(ROOT) or status not in ('completed', 'failed', 'cancelled'):
        raise ValueError('Invalid run folder or terminal status')
    manifest = folder/'run.json'
    record = json.loads(manifest.read_text(encoding='utf-8'))
    if record['status'] != 'prepared':
        raise ValueError('This run is already finalized; create a new run for re-analysis')
    if summary is not None:
        with (folder/'summary.json').open('x', encoding='utf-8') as stream:
            json.dump(summary, stream, ensure_ascii=False, indent=2)
    record.update(status=status, validation_status=validation, finished_at=datetime.now().astimezone().isoformat())
    manifest.write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('design')
    parser.add_argument('category', choices=CATEGORIES)
    parser.add_argument('tag')
    parser.add_argument('--entrypoint', required=True)
    parser.add_argument('--input', action='append', default=[])
    parser.add_argument('--upstream', action='append', default=[])
    args = parser.parse_args()
    print(create_run(args.design, args.category, args.tag, entrypoint=args.entrypoint,
                     inputs=args.input, upstream=args.upstream))
