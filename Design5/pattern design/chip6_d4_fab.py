"""Fabrication GDS of the Design5 six-channel Purcell chip (chip6_d4.py), 10 x 10 mm.

Lab fab rules (note/gds_export_settings.md, skill/chip_drc.md C10/C13), same as Design1/3:
  - 5 x 5 um square cheese holes on a 50 um pitch in the ground (layer 1),
  - no holes within 50 um of any etched feature, no holes in the outer 0.5 mm of the chip,
  - test-junction pockets in the four corners (800 x 300 um, pads 225/295 x 129 um, gap 211 um),
  - die marks (crosses in 400 um windows, 0.35 mm in from the corners) and the die ID.
GDS written in um (chip -5000 .. 5000).

Usage: python chip6_d4_fab.py [--die-id PF600]
  -> gds/Design5_<date>.gds, fig/fab/Design5_<date>_{chip,testjj,holes}.png;
  fab settings: section `fab` of parameters/Design5_<date>.json
"""
import json, os, sys, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import chip6_d4 as C
import design5_layout as d5

f3 = d5.f3                                   # Design3 ftf_layout: TEST_JJ, marks, fab options
FAB_GDS = os.path.join(ROOT, 'gds', C.TAG + '.gds')
FIG_DIR = os.path.join(ROOT, 'fig', 'fab')


def build_fab(die_id='PF600'):
    d = C.build_chip()
    t = f3.TEST_JJ
    (wl, hl), (wr, hr) = t['pad_left'], t['pad_right']
    for name, sx, sy in (('LL', -1, -1), ('LR', 1, -1), ('UL', -1, 1), ('UR', 1, 1)):
        f3.TestJunction(d, f'testjj_{name}', options=dict(
            chip='main', pos_x=f"{sx * t['x']}mm", pos_y=f"{sy * t['y']}mm",
            pocket_width=f"{t['pocket'][0]}mm", pocket_height=f"{t['pocket'][1]}mm",
            pad_l_width=f'{wl}mm', pad_l_height=f'{hl}mm', pad_r_width=f'{wr}mm', pad_r_height=f'{hr}mm',
            pad_gap=f"{t['pad_gap']}mm"))
    f3.add_corner_marks(d, die_id=die_id, inset=f3.MARK_INSET, **f3.MARK)
    return d


def export(d, filename):
    d.qgeometry.tables['junction'] = d.qgeometry.tables['junction'].iloc[0:0]
    gds = d.renderers.gds
    f3.apply_fab_gds_options(gds)
    gds.export_to_gds(filename)
    f3._prune_gds(filename)


def plots(filename):
    import gdspy, matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    top = [c for c in gdspy.GdsLibrary(infile=filename).top_level() if c.name == 'TOP'][0]
    spec = top.get_polygons(by_spec=True)
    os.makedirs(FIG_DIR, exist_ok=True)
    for tag, win, size in (('chip', (-5100, 5100, -5100, 5100), (9, 9)),
                           ('testjj', (-4450, -3450, 3950, 4550), (8, 5)),
                           ('holes', (600, 1100, 300, 800), (8, 8))):
        fig, ax = plt.subplots(figsize=size)
        for (layer, dt), ps in sorted(spec.items()):
            if layer == 1:
                ax.add_collection(PolyCollection(ps, facecolor='#9aa3ad' if dt == 100 else '#4b5560', edgecolor='none'))
        ax.set_xlim(win[:2]); ax.set_ylim(win[2:]); ax.set_aspect('equal')
        ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)'); ax.set_title(f'{C.TAG} ({tag})')
        fig.tight_layout(); fig.savefig(os.path.join(FIG_DIR, f'{C.TAG}_{tag}.png'), dpi=100); plt.close(fig)
    bb = top.get_bounding_box()
    n = {k: len(v) for k, v in spec.items()}
    print('bbox (um):', bb.tolist(), '\npolygons per layer/datatype:', n)


if __name__ == '__main__':
    import geopandas as gpd
    if not hasattr(gpd.GeoDataFrame, 'append'):
        gpd.GeoDataFrame.append = gpd.GeoDataFrame._append
    die = sys.argv[sys.argv.index('--die-id') + 1] if '--die-id' in sys.argv else 'PF600'
    d = build_fab(die)
    export(d, FAB_GDS)
    plots(FAB_GDS)
    print('gds', FAB_GDS)
