"""Where finished results are kept: one folder per run, grouped by date, never overwritten (traceability).

    results/<YYYY-MM-DD>/<HHMM>_<tag>/        HFSS / Q3D run: raw data (S.csv, cmat), meta, log, analysis
    results/<YYYY-MM-DD>/<HHMM>_<tag>/<tag>_compare_<YYYYmmdd_HHMM>.png / .json    every (re-)analysis, time-stamped

The files directly in results/ are only the current working values (chip6_lengths.json, q3d5_full_caps.json ...)
that other scripts read; their time-stamped copies go into the dated folders.
"""
import os, glob, time
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')


def stamp():
    return time.strftime('%Y%m%d_%H%M')


def new_run_dir(tag, res=RES):
    """results/<YYYY-MM-DD>/<HHMM>_<tag>/ (suffix _2, _3 ... if it exists already)."""
    d = os.path.join(res, time.strftime('%Y-%m-%d'), f'{time.strftime("%H%M")}_{tag}')
    out, k = d, 2
    while os.path.exists(out):
        out = f'{d}_{k}'; k += 1
    os.makedirs(out)
    return out


def find_run_dirs(tag, res=RES):
    """All run folders of `tag`, oldest first (dated layout and the older results/<stamp>_<tag>/ layout)."""
    hits = glob.glob(os.path.join(res, '*', f'*_{tag}')) + glob.glob(os.path.join(res, f'*_{tag}'))
    hits = [h for h in hits if os.path.isdir(h)]
    return sorted(hits, key=os.path.getmtime)
