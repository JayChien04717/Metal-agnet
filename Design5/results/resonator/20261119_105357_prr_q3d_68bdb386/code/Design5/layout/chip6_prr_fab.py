"""Fabrication GDS of the Design5 PurcellReadout chip (chip6_prr.py), 10 x 10 mm wafer die.

Same flow as chip6_d4_fab.export_wafer (r7): positive metal on layer 1/0, 1 nm wafer grid, 20 um ground cheese
outside the launch keep-outs, no corner crosses, no test junctions. The chip (design) ID, PR01 by default (user
choice 2026-11-19), is etched into the ground at the upper left, glyph height xmark.LABEL_SIZE_UM, 50 um from the left
and top die edges (Design6 TXT0x convention); the wafer die number goes at the upper right. The channel lengths are
read from the revision's build checks (parameters/<rev>_build_checks.json), so the fab mask has exactly the geometry
of gds/Design5_<rev>_sim.gds. HFSS wire bonds are not drawn. Air bridges at PurcellReadout.strap_points() are not on
this mask (separate process layer).

    python -I -B Design5/layout/chip6_prr_fab.py [--rev 2026-11-19_prr_r1] [--chip-id PR01]
        -> gds/Design5_<rev>_<id>_fab.gds, fig/Design5_<rev>_<id>_fab.png, fig/..._chip_id.png
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'code' / 'general'), str(ROOT / 'Design5' / 'layout')]
import chip6_prr as C  # noqa: E402  (configures paths, SI_EPS 11.5)
from gds_tools import export_gds, plot_gds, write_wafer_gds  # noqa: E402

HOME = ROOT / 'Design5'
CHIP_ID_INSET_UM = 50.0


def chip_id_glyphs(text):
    """Etched chip-ID glyphs (shapely, um) at the upper-left die corner."""
    from shapely.affinity import translate
    from xmark import LABEL_SIZE_UM, label_geometry
    g = label_geometry(text, LABEL_SIZE_UM, 10, 400, 30)
    x0, _, _, y1 = g.bounds
    return translate(g, -5000 + CHIP_ID_INSET_UM - x0, 5000 - CHIP_ID_INSET_UM - y1)


def etch_chip_id(mask, text):
    """Subtract the chip-ID glyphs from layer 1/0 and check the etched area against the glyphs."""
    import gdstk
    from shapely import union_all
    from shapely.geometry import Polygon, box
    glyphs = chip_id_glyphs(text)
    cut = [gdstk.Polygon(list(p.exterior.coords)[:-1]) for p in getattr(glyphs, 'geoms', [glyphs])]
    holes = [gdstk.Polygon(list(h.coords)[:-1]) for p in getattr(glyphs, 'geoms', [glyphs]) for h in p.interiors]
    if holes:
        cut = gdstk.boolean(cut, holes, 'not', precision=.001)
    src = gdstk.read_gds(str(mask))
    top, = src.top_level()
    top.flatten()                      # wafer GDS keeps metal in referenced tiles
    lib = gdstk.Library(unit=src.unit, precision=src.precision)
    lib.add(top)
    window = box(*glyphs.buffer(60).bounds).intersection(box(-5000, -5000, 5000, 5000))
    x0, y0, x1, y1 = window.bounds

    def near(cell):                    # layer-1 polygons whose bounding box meets the ID window
        out = []
        for p in cell.polygons:
            if p.layer != 1:
                continue
            (a, b), (c, d) = p.bounding_box()
            if a < x1 and c > x0 and b < y1 and d > y0:
                out.append(p)
        return out

    def metal(cell):
        clip = gdstk.rectangle((x0, y0), (x1, y1))
        pieces = gdstk.boolean(near(cell), clip, 'and', precision=.001)
        return union_all([Polygon(p.points).buffer(0) for p in pieces], grid_size=.001)
    before = metal(top)
    assert window.difference(before).area < 1, 'chip-ID corner is not solid ground'
    touched = near(top)                # only these are re-cut; the rest of the chip is untouched
    new = gdstk.boolean(touched, cut, 'not', precision=.001, layer=1, datatype=0)
    top.remove(*touched)
    top.add(*new)
    lib.write_gds(str(mask))
    check, = gdstk.read_gds(str(mask)).top_level()
    etched = window.difference(metal(check))
    err = etched.symmetric_difference(glyphs).area
    assert err < 0.02 * glyphs.area, f'chip ID mismatch {err:.2f} um^2'
    print(f'chip ID {text!r} etched at upper left, inset {CHIP_ID_INSET_UM} um; glyph XOR {err:.3f} um^2', flush=True)
    return glyphs


def build_fab(rev):
    checks = json.loads((HOME / 'parameters' / f'{rev}_build_checks.json').read_text(encoding='utf-8'))
    lens = {int(ch): v for ch, v in checks['channels'].items()}
    design = C.build_chip(lens, C.config(), wire_bonds=False)
    assert design.qgeometry.tables['junction'].empty
    return design, checks


def export_wafer(output, rev, chip_id='PR01'):
    from export_wafer_masks import SIDE_PORTS, verify_output
    from wafer_cheese_grid import replace_wafer_cheese

    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    design, checks = build_fab(rev)
    with TemporaryDirectory(dir=output.parent) as temp:
        raw, mask = Path(temp) / 'raw_fab.gds', Path(temp) / output.name
        export_gds(design, str(raw), fab=True, cheese=False)
        write_wafer_gds(raw, mask)
        holes = replace_wafer_cheese(mask, launches=SIDE_PORTS)
        verify_output(mask, 'Design5', True, expect_die_id=False)   # before the ID: r7 checks need bare corners
        if chip_id:
            etch_chip_id(mask, chip_id)
        mask.replace(output)
    print(f'Fab GDS: {output}', flush=True)
    return checks, holes


def main():
    import geopandas as gpd
    if not hasattr(gpd.GeoDataFrame, 'append'):
        gpd.GeoDataFrame.append = gpd.GeoDataFrame._append
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rev', default='2026-11-19_prr_r1')
    ap.add_argument('--chip-id', default='PR01')
    a = ap.parse_args()
    stem = f'Design5_{a.rev}' + (f'_{a.chip_id.lower()}' if a.chip_id else '')
    out = HOME / 'gds' / f'{stem}_fab.gds'
    _, holes = export_wafer(out, a.rev, a.chip_id)
    print('cheese:', holes, flush=True)
    plot_gds(str(out), str(HOME / 'fig' / f'{stem}_fab.png'), window=(-5.15, 5.15, -5.15, 5.15),
             ports=__import__('export_wafer_masks').SIDE_PORTS)
    if a.chip_id:
        plot_gds(str(out), str(HOME / 'fig' / f'{stem}_fab_chip_id.png'), window=(-5.05, -4.10, 4.55, 5.05),
                 figsize=(6, 3.6))


if __name__ == '__main__':
    main()
