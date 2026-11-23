"""Minimal run records: unique folder, input/source snapshots and explicit provenance.

Use create_run() before solving, write outputs into its returned directory, then
finish_run(). Existing legacy scripts are not silently redirected by this module, except for the qiskit-metal
eigenmode convergence CSV, which goes to the active run (project_paths.install_convergence_redirect).
"""
from datetime import datetime
from hashlib import sha256
from importlib.metadata import distributions
from pathlib import Path
import argparse
import ast
import json
import os
import re
import subprocess
import sys
from uuid import uuid4

ROOT = Path(__file__).resolve().parents[2]
CATEGORIES = ('general', 'qubit', 'resonator', 'flux')
ACTIVE_RUN = None   # set by create_run, cleared by finish_run; project_paths routes solver side files here
_EXCLUDED_PARTS = frozenset(('results', '__pycache__', '.migration-work', '.git'))


def _long_path(p):
    """Add Windows extended-length path prefix only when approaching MAX_PATH (260 chars)."""
    s = str(Path(p).resolve())
    if os.name == 'nt' and len(s) >= 248 and not s.startswith('\\\\?\\'):
        return Path('\\\\?\\' + s)
    return Path(s)


def _git(*args):
    result = subprocess.run(['git', '-C', str(ROOT), *args], capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else None


def _is_project_source(path):
    try:
        p = Path(path).resolve()
    except OSError:
        return False
    if p.suffix != '.py' or not p.is_file() or not p.is_relative_to(ROOT):
        return False
    return not any(part in _EXCLUDED_PARTS for part in p.relative_to(ROOT).parts)


def _loaded_project_sources():
    found = []
    for mod in list(sys.modules.values()):
        mod_file = getattr(mod, '__file__', None)
        if mod_file and _is_project_source(mod_file):
            found.append(Path(mod_file).resolve())
    return found


def _default_search_dirs(entry, design_root):
    if str(ROOT / 'code' / 'general') not in sys.path:
        sys.path.insert(0, str(ROOT / 'code' / 'general'))
    import project_paths

    candidates = [
        entry.parent,
        design_root / 'layout',
        design_root / 'layout' / 'components',
        *(design_root / 'simulation' / c for c in CATEGORIES),
        *(ROOT / 'code' / c for c in CATEGORIES),
    ]
    # Place pattern design/archive last so active QComponents take precedence.
    pattern_dirs = [Path(p) for p in project_paths.pattern_design_paths()]
    candidates += [p for p in pattern_dirs if p.name != 'archive']
    candidates += [p for p in pattern_dirs if p.name == 'archive']
    candidates.append(ROOT)
    for item in sys.path:
        if not item:
            continue
        try:
            p = Path(item).resolve()
        except OSError:
            continue
        if p.is_dir() and p.is_relative_to(ROOT) and not any(
            part in _EXCLUDED_PARTS for part in p.relative_to(ROOT).parts
        ):
            candidates.append(p)
    seen, dirs = set(), []
    for d in candidates:
        try:
            rd = d.resolve()
        except OSError:
            continue
        if rd.is_dir() and rd not in seen:
            seen.add(rd)
            dirs.append(rd)
    return dirs


def _resolve_module(mod_name, search_dirs):
    parts = [p for p in mod_name.split('.') if p]
    if not parts:
        return []
    hits = []
    for base in search_dirs:
        py_file = base.joinpath(*parts).with_suffix('.py')
        init_file = base.joinpath(*parts) / '__init__.py'
        if _is_project_source(py_file):
            hits.append(py_file.resolve())
            break
        if _is_project_source(init_file):
            hits.append(init_file.resolve())
            break
    # Also include parent package __init__.py files when importing a submodule.
    if hits and len(parts) > 1:
        for i in range(1, len(parts)):
            for base in search_dirs:
                parent_init = base.joinpath(*parts[:i]) / '__init__.py'
                if _is_project_source(parent_init):
                    hits.append(parent_init.resolve())
                    break
    return hits


def _scan_imports(source_path, search_dirs):
    """Return project .py files directly imported by source_path (AST scan, including lazy imports)."""
    try:
        tree = ast.parse(_long_path(source_path).read_text(encoding='utf-8'), filename=str(source_path))
    except (OSError, SyntaxError, UnicodeDecodeError):
        return []

    # Discover legacy cross-Design import_path(...) / _I(...) calls even before execution.
    import project_paths
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and node.args and isinstance(node.args[0], ast.Constant):
            arg = node.args[0].value
            if not isinstance(arg, str):
                continue
            fn = node.func
            fn_name = fn.id if isinstance(fn, ast.Name) else (fn.attr if isinstance(fn, ast.Attribute) else None)
            if fn_name in ('import_path', '_I'):
                try:
                    resolved = Path(project_paths.import_path(arg)).resolve()
                except Exception:
                    continue
                if resolved.is_dir() and resolved.is_relative_to(ROOT) and resolved not in search_dirs:
                    if not any(part in _EXCLUDED_PARTS for part in resolved.relative_to(ROOT).parts):
                        search_dirs.append(resolved)

    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                found.extend(_resolve_module(alias.name, search_dirs))
        elif isinstance(node, ast.ImportFrom):
            if node.level and node.level > 0:
                base = source_path.parent
                for _ in range(node.level - 1):
                    base = base.parent
                rel_dirs = [base]
                if node.module:
                    found.extend(_resolve_module(node.module, rel_dirs))
                    for alias in node.names:
                        if alias.name != '*':
                            found.extend(_resolve_module(f'{node.module}.{alias.name}', rel_dirs))
                else:
                    for alias in node.names:
                        if alias.name != '*':
                            found.extend(_resolve_module(alias.name, rel_dirs))
            elif node.module:
                found.extend(_resolve_module(node.module, search_dirs))
                for alias in node.names:
                    if alias.name != '*':
                        found.extend(_resolve_module(f'{node.module}.{alias.name}', search_dirs))
    return found


def _collect_project_sources(entry, design_root, extra_code=()):
    """Collect the entrypoint plus only project Python modules it imports (AST + runtime sys.modules)."""
    search_dirs = _default_search_dirs(entry, design_root)
    queue = [entry]
    for p in extra_code:
        rp = Path(p).resolve()
        if not _is_project_source(rp):
            raise ValueError(f'extra_code must be an active project .py file: {p}')
        queue.append(rp)
    queue.extend(_loaded_project_sources())

    visited = set()
    while queue:
        src = queue.pop(0)
        if src in visited or not _is_project_source(src):
            continue
        visited.add(src)
        for dep in _scan_imports(src, search_dirs):
            if dep not in visited:
                queue.append(dep)
    return sorted(visited, key=lambda p: p.relative_to(ROOT).as_posix())


def _flat_snapshot_name(rel_path, used_names):
    """Choose a short flat filename under code/, disambiguating with parent parts only on collision."""
    parts = [re.sub(r'[^\w.-]+', '_', part) for part in rel_path.parts]
    for depth in range(1, len(parts) + 1):
        candidate = '__'.join(parts[-depth:])
        if candidate not in used_names:
            used_names.add(candidate)
            return candidate
    base = '__'.join(parts)
    suffix = sha256(rel_path.as_posix().encode('utf-8')).hexdigest()[:8]
    candidate = f'{Path(base).stem}_{suffix}{Path(base).suffix}'
    used_names.add(candidate)
    return candidate


def _snapshot_sources(folder, sources, existing_entries=()):
    """Write flat snapshots of project .py files into <folder>/code/ and return updated manifest entries."""
    code_dir = folder / 'code'
    _long_path(code_dir).mkdir(parents=True, exist_ok=True)
    entries = list(existing_entries)
    recorded_paths = {item['path'] for item in entries}
    used_names = {Path(item.get('snapshot', item['path'])).name for item in entries}

    # Pre-check basename collisions among newly added sources so colliding peers both get qualified names.
    new_rels = []
    for src in sources:
        rel = src.relative_to(ROOT)
        if rel.as_posix() not in recorded_paths:
            new_rels.append((src, rel))

    name_counts = {}
    for _, rel in new_rels:
        clean = re.sub(r'[^\w.-]+', '_', rel.name)
        name_counts[clean] = name_counts.get(clean, 0) + 1

    for src, rel in new_rels:
        clean = re.sub(r'[^\w.-]+', '_', rel.name)
        if name_counts[clean] > 1 or clean in used_names:
            # Start disambiguation from parent__name when basenames collide.
            parts = [re.sub(r'[^\w.-]+', '_', part) for part in rel.parts]
            snap_name = None
            for depth in range(2, len(parts) + 1):
                candidate = '__'.join(parts[-depth:])
                if candidate not in used_names:
                    snap_name = candidate
                    used_names.add(candidate)
                    break
            if snap_name is None:
                snap_name = _flat_snapshot_name(rel, used_names)
        else:
            snap_name = clean
            used_names.add(snap_name)

        target = code_dir / snap_name
        data = _long_path(src).read_bytes()
        _long_path(target).write_bytes(data)
        entries.append({
            'path': rel.as_posix(),
            'snapshot': f'code/{snap_name}',
            'sha256': sha256(data).hexdigest(),
        })
    entries.sort(key=lambda item: item['path'])
    return entries


def create_run(design, category, tag, *, entrypoint, inputs=(), upstream=(), parameters=None, extra_code=()):
    """Create a run folder and snapshot only the entrypoint and imported project Python modules.

    Snapshots are stored flat under <run>/code/ (original repo-relative paths and
    SHA-256 hashes are recorded in run.json) to avoid copying unrelated Designs and
    exceeding Windows MAX_PATH. Non-Python inputs (geometry JSON, GDS, etc.) should
    be passed explicitly via inputs=. The caller owns actual output routing.
    """
    design_root = ROOT / design
    if category not in CATEGORIES or not (
        re.fullmatch(r'[Dd]esign[\w-]+', design) or (design_root / 'design.json').is_file()
    ):
        raise ValueError('Expected Design<N> and general/qubit/resonator/flux')
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
                  source_scope='Entrypoint and imported project Python modules (AST + runtime sys.modules); flat snapshots under code/ with original paths in run.json')
    try:
        for index, source in enumerate(inputs):
            target = folder / 'inputs' / str(index) / source.name
            _long_path(target.parent).mkdir(parents=True, exist_ok=True)
            data = _long_path(source).read_bytes()
            _long_path(target).write_bytes(data)
            record['inputs'].append({'source':str(source), 'snapshot':target.relative_to(folder).as_posix(),
                                     'sha256':sha256(data).hexdigest()})
        sources = _collect_project_sources(entry, design_root, extra_code=extra_code)
        record['code'] = _snapshot_sources(folder, sources)
    except Exception:
        record['status'] = 'preparation_failed'
        (folder/'run.json').write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        raise
    (folder/'run.json').write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    global ACTIVE_RUN
    ACTIVE_RUN = folder
    if str(ROOT / 'code' / 'general') not in sys.path:
        sys.path.insert(0, str(ROOT / 'code' / 'general'))
    from project_paths import install_convergence_redirect
    install_convergence_redirect()
    return folder


def finish_run(folder, *, status, validation='unknown', summary=None):
    folder = Path(folder).resolve()
    if not folder.is_relative_to(ROOT) or status not in ('completed', 'failed', 'cancelled'):
        raise ValueError('Invalid run folder or terminal status')
    manifest = folder/'run.json'
    record = json.loads(manifest.read_text(encoding='utf-8'))
    if record['status'] != 'prepared':
        raise ValueError('This run is already finalized; create a new run for re-analysis')
    # Capture any project Python modules imported lazily after create_run().
    record['code'] = _snapshot_sources(folder, _loaded_project_sources(), existing_entries=record.get('code', ()))
    if summary is not None:
        with (folder/'summary.json').open('x', encoding='utf-8') as stream:
            json.dump(summary, stream, ensure_ascii=False, indent=2)
    record.update(status=status, validation_status=validation, finished_at=datetime.now().astimezone().isoformat())
    manifest.write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    global ACTIVE_RUN
    if ACTIVE_RUN is not None and Path(ACTIVE_RUN).resolve() == folder:
        ACTIVE_RUN = None


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
