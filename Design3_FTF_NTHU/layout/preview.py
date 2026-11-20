"""Check and draw the rounded brackets."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from layout import build_design
from run_q3d import verify_layout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry', type=Path, help='Q3D geometry JSON for a specific run.')
    args = parser.parse_args()
    options = {}
    output = _P('Design3_FTF_NTHU/fig/bracket_preview.png')
    if args.geometry:
        geometry = json.loads(args.geometry.read_text(encoding='utf-8'))
        options['arm_lengths'] = geometry['length'] if 'length' in geometry else (geometry['arm_length'],) * 2
        output = args.geometry.with_name(args.geometry.stem.replace('_geometry', '_layout') + '.png')
    design = build_design(for_process=True, fx_options={'array_end_metal': 'False'}, **options)
    verify_layout(design)
    table = design.qgeometry.tables['poly']
    tc = design.components['TC']
    island = table.loc[(table.component == tc.id) & (table.name == 'tmon'), 'geometry'].item()
    junctions = design.qgeometry.tables['junction']
    for fx in ('FX1', 'FX2'):
        component = design.components[fx].id
        pads = table.loc[(table.component == component) & (table.name == 'pad_bot')]
        assert abs(island.distance(pads.geometry.item()) * 1000 - 4) < .02
        gap = min(island.distance(row.geometry.buffer(row.width / 2)) * 1000
                  for row in junctions.loc[junctions.component == component].itertuples())
        assert gap > 0, f'{fx} junction overlaps the TC island'
        print(f'{fx} nearest junction sheet clearance: {gap:.2f} um')

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), gridspec_kw={'width_ratios': [1, 2.4, 1]})
    panels = [((-430, -200), (-980, -805)), ((-450, 450), (-1110, -570)), ((200, 430), (-980, -805))]
    for ax, (xlim, ylim) in zip(axes, panels):
        for device, color in [('FX1', '#287db1'), ('TC', '#d48132'), ('FX2', '#287db1')]:
            rows = table.loc[(table.component == design.components[device].id) & ~table['subtract']]
            for row in rows.itertuples():
                vertices = [(x * 1000, y * 1000) for x, y in row.geometry.exterior.coords]
                ax.add_patch(Polygon(vertices, facecolor=color, edgecolor='#253341', linewidth=.55))
        ax.set(xlim=xlim, ylim=ylim, aspect='equal', xlabel='x (µm)')
    for ax, title in zip(axes, ('FX1 detail', 'Design3_FTF_NTHU', 'FX2 detail')):
        ax.set_title(title)
    axes[0].set_ylabel('y (µm)')
    fig.tight_layout()
    fig.savefig(output, dpi=170)


if __name__ == '__main__':
    main()
