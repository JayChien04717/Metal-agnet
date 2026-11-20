"""Resolve relocated project paths without changing historical solver records.

The migration map is an explicit old -> new index, not a latest-result selector.
No filesystem operations or solver connections happen here.
"""
from pathlib import Path
import json
import os
import sys
import glob as _glob

ROOT = Path(__file__).resolve().parents[2]
_MAP = json.loads((ROOT / 'code/general/path_migration.json').read_text(encoding='utf-8'))
_PREFIXES = sorted(_MAP['prefixes'].items(), key=lambda x: len(x[0]), reverse=True)


def relocate(value):
    """Map known old paths; leave unrelated paths and value types unchanged."""
    if not isinstance(value, (str, Path)):
        return value
    original = Path(value)
    absolute = original if original.is_absolute() else Path.cwd() / original
    try:
        relative = Path(os.path.normpath(absolute)).relative_to(ROOT).as_posix()
    except ValueError:
        return value
    parts = relative.split('/')
    if (len(parts) >= 3 and parts[1] == 'results'
            and parts[2] in ('general', 'qubit', 'resonator', 'flux')):
        # A mixed legacy root can contain old names of separately moved runs.
        aliases = _MAP.get('result_aliases', {})
        for old, new in sorted(aliases.items(), key=lambda x: len(x[0]), reverse=True):
            if relative == old or relative.startswith(old + '/'):
                result = ROOT / (new + relative[len(old):])
                return str(result) if isinstance(value, str) else result
        return value
    for old, new in _PREFIXES:
        if relative == old or relative.startswith(old + '/'):
            result = ROOT / (new + relative[len(old):])
            return str(result) if isinstance(value, str) else result
    return value


def path(relative):
    return relocate(ROOT / relative)


def string_path(relative):
    return str(path(relative))


def import_path(relative):
    """Resolve a former mixed script directory for an explicit sys.path entry."""
    parts = relative.split('/')
    if len(parts) == 2 and parts[1] == 'simulation':
        category = ('resonator' if parts[0] in ('Design4', 'Design5') else 'qubit')
        return str(ROOT / relative / category)
    return string_path(relative)


def glob_files(pattern, *, recursive=False):
    """Search categorized files as well as old loose-result names in the migration index."""
    hits = list(_glob.glob(str(pattern), recursive=recursive))
    absolute = os.path.normpath(os.path.abspath(pattern))
    import fnmatch
    for mapping in (_MAP['prefixes'], _MAP.get('result_aliases', {})):
        for old, new in mapping.items():
            candidate = os.path.normpath(str(ROOT / old))
            # Match direct historical file aliases; never fabricate a file or walk solver caches.
            if fnmatch.fnmatch(candidate, absolute):
                target = ROOT / new
                if target.is_file():
                    hits.append(str(target))
    return list(dict.fromkeys(hits))


def glob_paths(folder, pattern):
    return (Path(p) for p in glob_files(os.path.join(folder, pattern)))


PATTERN_GROUPS = ('fluxonium', 'transmon', 'resonator', 'marker', 'ebeam_junction')


def pattern_design_paths():
    """Shared QComponent folders: 'pattern design' (general parts) and its category subfolders.
    Modules keep flat names (``from my_tmon import ...``), so every folder goes on sys.path."""
    base = ROOT / 'pattern design'
    return [str(base)] + [str(base / group) for group in PATTERN_GROUPS]


def configure(filename):
    """Make local categories and shared utilities available for legacy bare imports."""
    file = Path(filename).resolve()
    relative = file.relative_to(ROOT)
    folders = []
    if relative.parts[0].startswith('Design'):
        design = ROOT / relative.parts[0]
        folders += [file.parent, design / 'layout']
        folders += [design / 'simulation' / k for k in ('general', 'qubit', 'resonator', 'flux')]
    folders += [ROOT / 'code' / k for k in ('general', 'qubit', 'resonator', 'flux')]
    folders += [Path(p) for p in pattern_design_paths()]
    for folder in reversed(folders):
        if folder.is_dir() and str(folder) not in sys.path:
            sys.path.insert(0, str(folder))
