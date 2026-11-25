"""Report metal narrower than the minimum width (default 4 um) in GDS masks.

    python -I code/general/check_metal_width.py Design6_TXT/gds/<mask>.gds [--min-width 4]
        [--layer 1 --datatype 0] [--window xmin ymin xmax ymax] [--png out.png]
        [--exempt xmin ymin xmax ymax ...]

Exempt by rule (user decision 2026-11-18): resonator parts (readout / Purcell resonators, coupled-line
tees, resonator capacitors; their coupling strips set kappa) and the etched die-ID / marker label.
With a qiskit-metal design, design_exemptions(design) returns these regions; for a GDS alone, pass
them with --exempt (um).

Units: um. Checks one metal layer (positive metal, default 1/0) of the top cell. Rule and
algorithm: pattern design/metal_width.py, skill/chip_drc.md C12a. Exit code 1 when any
violation is found. Cheesed masks are slow; use --window to inspect a region.
"""
from pathlib import Path
import argparse
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'pattern design' / g) for g in ('', 'fluxonium', 'transmon', 'resonator', 'marker', 'ebeam_junction')]

import gdstk
from shapely.geometry import Polygon, box
from shapely.ops import unary_union

from metal_width import MIN_METAL_WIDTH_UM, thin_metal


def gds_metal(path, layer=1, datatype=0, window=None):
    library = gdstk.read_gds(str(path))
    if abs(library.unit - 1e-6) > 1e-12:
        raise ValueError(f'{path}: GDS unit {library.unit} m, expected 1e-6 (um)')
    top, = library.top_level()
    polygons = top.get_polygons(layer=layer, datatype=datatype)
    if window is not None:
        polygons = gdstk.boolean(polygons, gdstk.rectangle(window[:2], window[2:]), 'and', precision=1e-3)
    return unary_union([Polygon(p.points).buffer(0) for p in polygons])


RESONATOR_CLASSES = {'RouteMeander', 'CoupledLineTee', 'HammerCap', 'EnclosingCap', 'InterdigitalCap',
                     'ResEnd', 'OpenToGround', 'ShortToGround', 'TapJunction'}


def design_exemptions(design, min_width=MIN_METAL_WIDTH_UM):
    """Exempt regions (um) of a qiskit-metal design: resonator components (RESONATOR_CLASSES or a class
    defined in pattern design/resonator) grown by min_width, and XMark labels grown by min_width."""
    from shapely.affinity import scale
    unit = design.parse_value('1um')
    regions = []
    for component in design.components.values():
        cls = type(component)
        module = getattr(sys.modules.get(cls.__module__), '__file__', '') or ''
        resonator = cls.__name__ in RESONATOR_CLASSES or Path(module).parent.name == 'resonator'
        for kind in ('poly', 'path'):
            table = design.qgeometry.tables[kind]
            rows = table[table.component == component.id]
            if not resonator:
                rows = rows[rows['name'].str.startswith('label')] if cls.__name__ == 'XMark' else rows.iloc[0:0]
            for row in rows.itertuples():
                geom = row.geometry
                if kind == 'path':      # route 'cut' rows already carry width + 2 gap
                    geom = geom.buffer(row.width / 2, cap_style=2)
                regions.append(scale(geom, 1 / unit, 1 / unit, origin=(0, 0)).buffer(min_width))
    return unary_union(regions) if regions else None


def check(path, min_width=MIN_METAL_WIDTH_UM, layer=1, datatype=0, window=None, exempt=None):
    """[(region, width_um)] below min_width; regions cut by the window edge, and regions inside
    `exempt` (shapely geometry or list of (xmin, ymin, xmax, ymax), um), are dropped."""
    metal = gds_metal(path, layer, datatype, window)
    found = thin_metal(metal, min_width)
    if window is not None:
        inner = box(*window).buffer(-min_width)
        found = [(g, w) for g, w in found if g.within(inner)]
    if exempt is not None:
        if isinstance(exempt, (list, tuple)):
            exempt = unary_union([box(*b) for b in exempt]) if exempt else None
        if exempt is not None:
            found = [(g, w) for g, w in found if not g.within(exempt)]
    return found


def plot(path, found, png, layer=1, datatype=0, pad=25):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n = len(found)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4), squeeze=False)
    for ax, (region, width) in zip(axes[0], found):
        x0, y0, x1, y1 = region.bounds
        win = (x0 - pad, y0 - pad, x1 + pad, y1 + pad)
        metal = gds_metal(path, layer, datatype, win)
        for g in getattr(metal, 'geoms', [metal]):
            ax.fill(*g.exterior.xy, color='#9aa3ad')
            for hole in g.interiors:
                ax.fill(*hole.xy, color='white')
        for g in getattr(region, 'geoms', [region]):
            ax.fill(*g.exterior.xy, color='red')
        ax.set_xlim(win[0], win[2]); ax.set_ylim(win[1], win[3]); ax.set_aspect('equal')
        ax.set_title(f'{width:.2f} um at ({region.centroid.x:.0f}, {region.centroid.y:.0f})', fontsize=9)
    fig.tight_layout()
    fig.savefig(png, dpi=120)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gds', type=Path, nargs='+')
    parser.add_argument('--min-width', type=float, default=MIN_METAL_WIDTH_UM)
    parser.add_argument('--layer', type=int, default=1)
    parser.add_argument('--datatype', type=int, default=0)
    parser.add_argument('--window', type=float, nargs=4)
    parser.add_argument('--png', type=Path, help='plot every violation (first file only)')
    parser.add_argument('--exempt', type=float, nargs=4, action='append', default=[],
                        metavar=('XMIN', 'YMIN', 'XMAX', 'YMAX'), help='exempt box in um (repeatable)')
    args = parser.parse_args()
    failed = False
    for i, path in enumerate(args.gds):
        found = check(path, args.min_width, args.layer, args.datatype, args.window, args.exempt)
        print(f'{path}: {len(found)} region(s) below {args.min_width} um')
        for region, width in sorted(found, key=lambda item: item[1]):
            x0, y0, x1, y1 = region.bounds
            print(f'  {width:6.3f} um  bounds ({x0:.1f}, {y0:.1f}) - ({x1:.1f}, {y1:.1f})')
        if found and args.png and i == 0:
            plot(path, found, args.png, args.layer, args.datatype)
        failed |= bool(found)
    sys.exit(1 if failed else 0)
