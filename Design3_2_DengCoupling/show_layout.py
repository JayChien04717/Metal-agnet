"""Plot the two coupling electrodes directly from the rounded solver geometry."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from run_q3d import build_variant, verify_layout

parser = argparse.ArgumentParser()
parser.add_argument('geometry', type=Path)
args = parser.parse_args()
config = json.loads(args.geometry.read_text(encoding='utf8'))
d = build_variant(config['gap'],config['length'],rounded=True,for_process=True,
                  tmon_options=dict(pad_head_width=f"{config['stem_width']}um",paddle_gap=f"{config['tc_etch']}um"))
verify_layout(d)
table = d.qgeometry.tables['poly']
fig, ax = plt.subplots(figsize=(10,5))
for q,color in (('FX1','#287db1'),('FX2','#287db1'),('TC','#d48132')):
    for row in table.loc[(table.component==d.components[q].id)&(~table['subtract'])].itertuples():
        part=row.geometry
        ax.add_patch(Polygon([(1000*x,1000*y) for x,y in part.exterior.coords],facecolor=color,edgecolor='black',lw=.5))
        for ring in part.interiors:
            ax.add_patch(Polygon([(1000*x,1000*y) for x,y in ring.coords],facecolor='white',edgecolor='black',lw=.5))
ax.set_aspect('equal')
ax.set(xlim=(-460,460),ylim=(-1140,-690),xlabel='x (µm)',ylabel='y (µm)',
       title='Design3_2 · FX–TC coupling electrodes')
fig.tight_layout()
fig.savefig(args.geometry.with_name(args.geometry.stem.replace('_geometry','_layout')+'.png'),dpi=160)
