"""Render the standard figures of every Design from its current GDS, in one style.

Each Design's ``fig/`` folder holds only:
- ``chip.png``: the whole chip, with each part window boxed and labelled by its file name. A Design with
  several chip variants has one ``chip_<variant>.png`` per variant.
- ``<part>.png``: one close-up per part (qubit, coupler, purcell_resonator_pair, ...). With variants the name
  is ``<variant>_<part>.png``.
- ``README.md``: generated. It lists each figure, its source GDS, the revision and the window.
- ``.gitignore``: generated. Only the files above are tracked; anything else that a script writes into
  ``fig/`` (check plots, design-ID previews) stays local. Superseded figures live in the ignored ``previous/``.

Figures that support a note (data plots, method figures) belong in ``notes/fig/`` next to that note, not here.

What to draw comes from ``code/general/design_figures.json``. All layer-1 datatypes are drawn as metal (grey).
Layer 50/0 junction annotations in simulation GDS files are drawn red. White is bare substrate.

Usage (repository root, qmetal):
    python -X utf8 -B code/general/plot_design_figs.py            # every Design
    python -X utf8 -B code/general/plot_design_figs.py Design6_TXT    # one Design
"""

import argparse
import hashlib
import json
from pathlib import Path

import gdstk
import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SPEC = Path(__file__).with_name('design_figures.json')
METAL, JUNCTION, EDGE, BOX = '#7d8a98', '#d7322b', '#2b3540', '#d98b00'


def load(path):
    top, = gdstk.read_gds(str(path)).top_level()
    metal = [p for p in top.get_polygons() if p.layer == 1]
    junction = top.get_polygons(layer=50, datatype=0)
    return metal, junction, np.asarray(top.bounding_box())


def clip(polygons, window):
    x0, x1, y0, y1 = window
    keep = []
    for p in polygons:
        (a, b), (c, d) = p.bounding_box()
        if a < x1 and c > x0 and b < y1 and d > y0:
            keep.append(p.points)
    return keep


def scale_bar(ax, window, unit):
    x0, x1, y0, y1 = window
    span = x1 - x0
    nice = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
    length = max(v for v in nice if v <= span / 4)
    xs, ys = x0 + .05 * span, y0 + .05 * (y1 - y0)
    ax.add_patch(Rectangle((xs, ys), length, .012 * (y1 - y0), fc='black', ec='none', zorder=5))
    label = f'{length / 1000:g} mm' if length >= 1000 else f'{length:g} µm'
    ax.text(xs + length / 2, ys + .03 * (y1 - y0), label, ha='center', va='bottom', fontsize=9, zorder=5,
            bbox=dict(fc='white', ec='none', alpha=.8, pad=1))


def draw(ax, metal, junction, window):
    ax.set_facecolor('white')
    # A hairline edge in the metal colour hides the anti-aliasing seams between touching polygons.
    ax.add_collection(PolyCollection(clip(metal, window), fc=METAL, ec=METAL, lw=.25))
    j = clip(junction, window)
    if j:
        ax.add_collection(PolyCollection(j, fc=JUNCTION, ec='none', zorder=3))
    x0, x1, y0, y1 = window
    ax.set(xlim=(x0, x1), ylim=(y0, y1), aspect='equal')
    ax.set_xticks([]); ax.set_yticks([])
    for side in ax.spines.values():
        side.set_color(EDGE); side.set_linewidth(.8)


def render_chip(path, title, metal, junction, window, parts):
    fig, ax = plt.subplots(figsize=(8, 8 * (window[3] - window[2]) / (window[1] - window[0])))
    draw(ax, metal, junction, window)
    for name, part in parts.items():
        x0, x1, y0, y1 = part['window']
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=BOX, lw=1.4, zorder=4))
        ax.text(x0, y1, f' {name}', ha='left', va='bottom', fontsize=8, color=BOX, weight='bold',
                zorder=4, bbox=dict(fc='white', ec='none', alpha=.75, pad=.5))
    scale_bar(ax, window, 'um')
    fig.suptitle(title, x=.01, ha='left', fontsize=11, weight='bold')
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor='white')
    plt.close(fig)


def fit(window, low=.75, high=2.5):
    """Widen a very narrow (or flatten a very wide) window about its centre so every close-up reads well."""
    x0, x1, y0, y1 = window
    w, h, cx, cy = x1 - x0, y1 - y0, (x0 + x1) / 2, (y0 + y1) / 2
    if w / h < low:
        w = h * low
    elif w / h > high:
        h = w / high
    return [round(cx - w / 2, 3), round(cx + w / 2, 3), round(cy - h / 2, 3), round(cy + h / 2, 3)]


def render_part(path, title, metal, junction, window):
    w, h = window[1] - window[0], window[3] - window[2]
    fig, ax = plt.subplots(figsize=(6, max(2.6, min(8.5, 6 * h / w))))
    draw(ax, metal, junction, window)
    scale_bar(ax, window, 'um')
    fig.suptitle(title, x=.01, ha='left', fontsize=11, weight='bold')
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor='white')
    plt.close(fig)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()[:12]


def render_design(name, spec):
    folder = ROOT / name / 'fig'
    folder.mkdir(exist_ok=True)
    rows, written = [], []
    for variant, chip in spec['chips'].items():
        prefix = '' if variant == 'main' else f'{variant}_'
        gds = ROOT / name / chip['gds']
        metal, junction, bbox = load(gds)
        window = chip.get('window') or [bbox[0][0], bbox[1][0], bbox[0][1], bbox[1][1]]
        label = chip.get('label', '')
        chip_png = folder / ('chip.png' if variant == 'main' else f'chip_{variant}.png')
        render_chip(chip_png, f"{name} · {chip['revision']}{' · ' + label if label else ''} · whole chip",
                    metal, junction, window, {k: dict(v, window=fit(v['window'])) for k, v in chip.get('parts', {}).items()})
        rows.append((chip_png.name, 'whole chip' + (f' ({label})' if label else ''), chip, window))
        written.append(chip_png)
        for part, info in chip.get('parts', {}).items():
            png = folder / f'{prefix}{part}.png'
            window = fit(info['window'])
            where = f' ({label})' if label else ''
            render_part(png, f"{name} · {chip['revision']} · {info['label']}{where}", metal, junction, window)
            rows.append((png.name, info['label'] + where, chip, window))
            written.append(png)
    lines = [f'# {name} figures', '',
             'Generated by [`plot_design_figs.py`](../../code/general/plot_design_figs.py) from',
             '[`design_figures.json`](../../code/general/design_figures.json). Do not edit by hand.',
             '',
             '- `chip*.png` is the whole chip; each `<part>.png` is one close-up.',
             '- Grey is metal and white is bare substrate. Red marks the junction annotations (layer 50/0) in',
             '  simulation GDS files.',
             '- The figures follow the source GDS revision below. Re-render after any geometry change.', '',
             '| Figure | Shows | Revision | Source GDS (sha256, first 12) | Window x0, x1, y0, y1 (µm) |',
             '|---|---|---|---|---|']
    for png, shows, chip, window in rows:
        gds = ROOT / name / chip['gds']
        lines.append(f"| [{png}]({png}) | {shows} | {chip['revision']} | [`{gds.name}`](../{chip['gds']}) "
                     f"`{sha(gds)}` | {', '.join(f'{v:g}' for v in window)} |")
    if spec.get('note'):
        lines += ['', spec['note']]
    (folder / 'README.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    allow = ['*', '!.gitignore', '!README.md'] + [f'!{p.name}' for p in written]
    (folder / '.gitignore').write_text('# Generated by code/general/plot_design_figs.py: only the standard set is tracked.\n'
                                       + '\n'.join(allow) + '\n', encoding='utf-8')
    return written


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('designs', nargs='*')
    args = parser.parse_args()
    spec = json.loads(SPEC.read_text(encoding='utf-8'))
    for name in args.designs or [k for k in spec if not k.startswith('_')]:
        written = render_design(name, spec[name])
        print(name, len(written), 'figures')


if __name__ == '__main__':
    main()
