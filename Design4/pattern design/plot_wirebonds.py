"""Layout with all ground bonds: qiskit-metal auto wire bonds (P.WB_OPTIONS) + hammer / tee straps.
Usage: python plot_wirebonds.py [--gap 45]  -> ../fig/purcell_wide400_gap<gap>_wirebonds.png"""
import argparse, warnings
warnings.filterwarnings('ignore')
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import purcell_6500 as P

ap = argparse.ArgumentParser(); ap.add_argument('--gap', type=float, default=45); a = ap.parse_args()
d = P.build(P.LENGTHS, paper_hammer=True, paper_hammer_gap_um=a.gap, spacing_mm=0.4)
P.set_wire_bonds(d)
fig, ax = plt.subplots(figsize=(8, 7.5))
for r in d.qgeometry.tables['path'].itertuples():
    if r.subtract:
        ax.plot(*np.asarray(r.geometry.coords).T, color='0.55', lw=2.5, solid_capstyle='butt')
for r in d.qgeometry.tables['poly'].itertuples():
    g = r.geometry; gs = getattr(g, 'geoms', [g])
    for p in gs:
        ax.fill(*p.exterior.xy, color='0.8' if r.subtract else 'goldenrod', lw=0)
span = 3 * 0.022
auto = P.wire_bond_positions(d, float(P.WB_OPTIONS['wb_threshold'][:-2]) / 1000)
for x, y, ox, oy in auto:
    ax.plot([x - ox * span / 2, x + ox * span / 2], [y - oy * span / 2, y + oy * span / 2], color='red', lw=2.5)
straps = P.ground_straps(d, 'all')
for s in straps:
    (x, y), (ox, oy) = s['pos'], s['ori']
    ax.plot([x - ox * s['span'] / 2, x + ox * s['span'] / 2], [y - oy * s['span'] / 2, y + oy * s['span'] / 2], color='blue', lw=2.5)
ax.plot([], [], 'r', lw=2.5, label=f'qiskit auto wire bonds ({len(auto)}, every {P.WB_OPTIONS["wb_threshold"]})')
ax.plot([], [], 'b', lw=2.5, label=f'hammer-stub / tee straps ({len(straps)})')
ax.set_aspect('equal'); ax.set_xlim(-1.75, 1.75); ax.set_ylim(-0.3, 3.0); ax.legend(loc='upper left', fontsize=8)
ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')
ax.set_title(f'wide layout, 400 um spacing, hammer gap {a.gap:g} um: ground bonds in HFSS / Q3D', fontsize=9)
out = f'../fig/purcell_wide400_gap{a.gap:g}_wirebonds.png'
fig.tight_layout(); fig.savefig(out, dpi=130); print(out, len(auto), len(straps))
