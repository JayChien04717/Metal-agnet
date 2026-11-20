"""Where finished results are kept: one folder per run, grouped by date, never overwritten (traceability).

    results/<YYYY-MM-DD>/<HHMM>_<tag>/        HFSS / Q3D run: raw data (S.csv, cmat), meta, log, analysis
    results/<YYYY-MM-DD>/<HHMM>_<tag>/<tag>_compare_<YYYYmmdd_HHMM>.png / .json    every (re-)analysis, time-stamped

The files directly in results/ are only the current working values (chip6_lengths.json, q3d5_full_caps.json ...)
that other scripts read; their time-stamped copies go into the dated folders.
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I, glob_files as _glob_files, glob_paths as _glob_paths
_configure(__file__)
import os, glob, time
RES = _S('Design5_6_purcell-res/results/resonator/legacy')


def stamp():
    return time.strftime('%Y%m%d_%H%M')


def new_run_dir(tag, res=RES):
    """results/<YYYY-MM-DD>/<HHMM>_<tag>/ (suffix _2, _3 ... if it exists already)."""
    # Legacy working values stay readable; new runs go beside legacy, not inside it.
    if os.path.basename(res) == 'legacy':
        res = os.path.dirname(res)
    d = _relocate(os.path.join(res, time.strftime('%Y-%m-%d'), f'{time.strftime("%H%M")}_{tag}'))
    out, k = d, 2
    while os.path.exists(out):
        out = f'{d}_{k}'; k += 1
    os.makedirs(out)
    return out


def find_run_dirs(tag, res=RES):
    """All run folders of `tag`, oldest first (dated layout and the older results/<stamp>_<tag>/ layout)."""
    roots = [res]
    if os.path.basename(res) == 'legacy':
        roots.append(os.path.dirname(res))
    hits = [hit for root in roots for pattern in (os.path.join(root, '*', f'*_{tag}'),
                                                 os.path.join(root, f'*_{tag}'))
            for hit in _glob_files(_relocate(pattern))]
    hits = [h for h in hits if os.path.isdir(h)]
    return sorted(hits, key=os.path.getmtime)
