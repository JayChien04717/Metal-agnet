"""Resolve relocated project paths without changing historical solver records.

The migration map is an explicit old -> new index, not a latest-result selector.
No filesystem operations or solver connections happen here (importing it only installs the AEDT session cleanup
hook from aedt_session.py, which acts when pyEPR connects to AEDT).
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
        category = ('resonator' if parts[0] in ('Design4_Single_purcell-res', 'Design5_6_purcell-res') else 'qubit')
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
    Modules keep flat names (``from my_tmon_dolan import ...``), so every folder goes on sys.path."""
    base = ROOT / 'pattern design'
    return [str(base)] + [str(base / group) for group in PATTERN_GROUPS]


def design_dir(filename):
    """Return the enclosing Design directory (containing design.json) for a script path.

    Keeps layout and simulation scripts independent of the Design folder name.
    """
    file = Path(filename).resolve()
    relative = file.relative_to(ROOT)
    for parent in (file if file.is_dir() else file.parent).parents:
        if parent == ROOT:
            break
        if (parent / 'design.json').is_file():
            return parent
    top = ROOT / relative.parts[0]
    if (top / 'design.json').is_file() or relative.parts[0].startswith('Design'):
        return top
    raise ValueError(f'{file} is not inside a Design directory under {ROOT}')


def design_name(filename):
    """Return the enclosing Design folder name for a script path."""
    return design_dir(filename).name


def all_design_dirs():
    """Return all top-level Design directories under ROOT (detected by design.json)."""
    return sorted(p for p in ROOT.iterdir() if p.is_dir() and (p / 'design.json').is_file())


CODE_GROUPS = ('general', 'aedt', 'layout', 'fab', 'ebeam', 'qubit', 'resonator', 'flux')


def code_paths():
    """All active shared code subfolders under ``code/``."""
    return [str(ROOT / 'code' / k) for k in CODE_GROUPS if (ROOT / 'code' / k).is_dir()]


def configure(filename):
    """Make local categories and shared utilities available for legacy bare imports."""
    file = Path(filename).resolve()
    relative = file.relative_to(ROOT)
    folders = []
    top = ROOT / relative.parts[0]
    if (top / 'design.json').is_file() or relative.parts[0].startswith('Design'):
        design = design_dir(file)
        folders += [file.parent, design / 'layout']
        folders += [design / 'simulation' / k for k in ('general', 'qubit', 'resonator', 'flux')]
    folders += [ROOT / 'code' / k for k in CODE_GROUPS]
    folders += [Path(p) for p in pattern_design_paths()]
    for folder in reversed(folders):
        if folder.is_dir() and str(folder) not in sys.path:
            sys.path.insert(0, str(folder))
    install_convergence_redirect(file)


# ---- HFSS eigenmode convergence CSV -------------------------------------------------------------------------------
# qiskit-metal's QHFSSRenderer.get_f_convergence exports 'hfss_eig_f_convergence.csv' into the current working
# directory (normally the repository root), overwriting the previous one. The redirect below sends it to the active
# run (run_record.create_run -> <run>/raw/) or, for legacy scripts without a run record, to
# <Design>/results/<category>/<script>_<timestamp>_hfss_eig_f_convergence.csv. Installed by configure() and by
# run_record.create_run(); the renderer module is patched when (or if) it is imported.
_ENTRY = {'file': None}
_RENDERER = 'qiskit_metal.renderers.renderer_ansys.hfss_renderer'


def convergence_target():
    """(folder, filename) for the next convergence CSV."""
    from datetime import datetime
    run = sys.modules.get('run_record') and getattr(sys.modules['run_record'], 'ACTIVE_RUN', None)
    if run:
        run = Path(run)
        return (run / 'raw' if (run / 'raw').is_dir() else run), 'hfss_eig_f_convergence.csv'
    file = _ENTRY['file'] or Path(sys.argv[0] or 'interactive').resolve()
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    name = f'{file.stem}_{stamp}_hfss_eig_f_convergence.csv'
    try:
        parts = file.relative_to(ROOT).parts
    except ValueError:
        parts = ()
    if parts and parts[0].startswith('Design'):
        cat = parts[2] if len(parts) > 3 and parts[1] == 'simulation' and parts[2] in (
            'general', 'qubit', 'resonator', 'flux') else 'general'
        return ROOT / parts[0] / 'results' / cat, name
    return file.parent, name


def _patch_renderer(module):
    cls = module.QHFSSRenderer
    if getattr(cls.get_f_convergence, '_redirected', False):
        return
    original = cls.get_f_convergence

    def get_f_convergence(self, *args, **kwargs):
        from contextlib import chdir
        folder, name = convergence_target()
        folder.mkdir(parents=True, exist_ok=True)
        with chdir(folder):
            out = original(self, *args, **kwargs)
        written = folder / 'hfss_eig_f_convergence.csv'
        if name != written.name and written.exists():
            written.replace(folder / name)
        return out
    get_f_convergence._redirected = True
    cls.get_f_convergence = get_f_convergence


def install_convergence_redirect(entry=None):
    if entry and _ENTRY['file'] is None:
        _ENTRY['file'] = Path(entry).resolve()
    if _RENDERER in sys.modules:
        _patch_renderer(sys.modules[_RENDERER])
        return
    if any(isinstance(f, _RendererHook) for f in sys.meta_path):
        return
    sys.meta_path.insert(0, _RendererHook())


class _RendererHook:
    """Meta-path finder that patches the renderer right after its module is executed."""

    def find_spec(self, name, path=None, target=None):
        if name != _RENDERER:
            return None
        import importlib.machinery
        spec = importlib.machinery.PathFinder.find_spec(name, path)
        if spec is None or spec.loader is None:
            return spec
        loader = spec.loader
        exec_module = loader.exec_module

        def exec_and_patch(module):
            exec_module(module)
            _patch_renderer(module)
        loader.exec_module = exec_and_patch
        return spec


# ---- AEDT session cleanup -----------------------------------------------------------------------------------------
# An AEDT that a run launches is saved, closed and quit when the run ends (aedt_session.py; AEDT_KEEP_OPEN=1 keeps
# it). An AEDT that was already open is left alone. Installed on import: pyEPR.ansys is patched when it is imported.
def _install_aedt_session():
    import importlib.util
    if 'aedt_session' in sys.modules:
        sys.modules['aedt_session'].install()
        return
    aedt_path = ROOT / 'code' / 'aedt' / 'aedt_session.py'
    if not aedt_path.is_file():
        aedt_path = Path(__file__).with_name('aedt_session.py')
    spec = importlib.util.spec_from_file_location('aedt_session', aedt_path)
    if spec is None:
        return
    module = importlib.util.module_from_spec(spec)
    sys.modules['aedt_session'] = module
    spec.loader.exec_module(module)
    module.install()


_install_aedt_session()
