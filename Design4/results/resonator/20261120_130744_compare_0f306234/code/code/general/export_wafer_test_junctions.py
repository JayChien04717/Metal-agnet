"""Export wafer fabrication masks with the current side launch coordinates."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import subprocess
from tempfile import TemporaryDirectory

import gdstk
from shapely import STRtree, get_parts, union_all
from shapely.affinity import scale, translate
from shapely.geometry import Point, Polygon, box

from export_wafer_masks import SIDE_PORTS, verify_output
from gds_tools import export_gds, write_wafer_gds
from project_paths import ROOT
from wafer_cheese_grid import replace_wafer_cheese

sys.path.insert(0, str(ROOT / 'DesignTestJunction' / 'layout'))
from build import build_design, GEOMETRY, read_region
from wafer_pcbport import wafer_pcb_ports
from xmark import CORNERS, LABEL_SIZE_UM, label_geometry


DESIGNS = ('Design1', 'Design2', 'Design3', 'Design4', 'Design5', 'Design6')   # Design3 = main FTF (was Design3_3)
GRID = .001  # µm
AREA_TOL = .2  # µm², allows for vertex rounding on the 1 nm grid


def polygons(shape):
    """Convert a region in µm to positive metal, preserving holes."""
    result = []
    for part in get_parts(shape):
        if part.geom_type != 'Polygon' or not part.area:
            continue
        shell = gdstk.Polygon(part.exterior.coords[:-1], layer=1)
        holes = [gdstk.Polygon(ring.coords[:-1]) for ring in part.interiors]
        result.extend(gdstk.boolean(shell, holes, 'not', layer=1, precision=GRID) if holes else [shell])
    return result


def region_in(window, metal, index):
    fragments = [metal[i] for i in index.query(window)]
    clip = gdstk.rectangle(window.bounds[:2], window.bounds[2:])
    return read_region(gdstk.boolean(fragments, clip, 'and', precision=GRID))


def move_wafer_die_id(mask, settings, *, required=True):
    """Move the existing etched glyphs to the upper-left die margin."""
    expected = label_geometry(settings['text'], LABEL_SIZE_UM, 10, 400, 30)
    if settings['source_corner'] == 'LR':
        x0, _, x1, _ = expected.bounds
        expected = translate(expected, 4650 - x0 - x1, -4650)
    else:
        expected = translate(expected, -4650, 4650)
    old_window = expected.envelope.buffer(8, join_style=2)

    library = gdstk.read_gds(str(mask))
    original, = library.top_level()
    all_polygons = original.get_polygons()
    metal = [p for p in all_polygons if (p.layer, p.datatype) == (1, 0)]
    index = STRtree([box(*low, *high) for low, high in (p.bounding_box() for p in metal)])
    etched = old_window.difference(region_in(old_window, metal, index))
    if not required and etched.area < AREA_TOL:
        return False
    assert etched.intersection(expected).area > .85 * expected.area, (mask, 'source die ID')
    assert etched.difference(expected.buffer(4)).area < AREA_TOL, (mask, 'unrelated etch near die ID')
    etched = union_all([part for part in get_parts(etched) if part.area >= AREA_TOL], grid_size=GRID)

    margin = settings['edge_inset_um']
    xmin, _, _, ymax = etched.bounds
    label = translate(etched, -5000 + margin - xmin, 5000 - margin - ymax)
    new_window = label.envelope.buffer(5, join_style=2)
    assert new_window.difference(region_in(new_window, metal, index)).area < AREA_TOL, (mask, 'occupied die-ID location')

    result = gdstk.Library(unit=library.unit, precision=library.precision)
    top = result.new_cell(original.name)
    cuts = polygons(label)
    hits = set(index.query(new_window))
    for i, polygon in enumerate(metal):
        top.add(*(gdstk.boolean(polygon, cuts, 'not', precision=GRID, layer=1) if i in hits else [polygon]))
    top.add(*(p for p in all_polygons if (p.layer, p.datatype) != (1, 0)))
    top.add(*polygons(etched))
    result.write_gds(str(mask), max_points=199)

    top, = gdstk.read_gds(str(mask)).top_level()
    checked = [p for p in top.get_polygons() if (p.layer, p.datatype) == (1, 0)]
    index = STRtree([box(*low, *high) for low, high in (p.bounding_box() for p in checked)])
    assert old_window.difference(region_in(old_window, checked, index)).area < AREA_TOL, (mask, 'old ID remains')
    moved = new_window.difference(region_in(new_window, checked, index))
    assert moved.symmetric_difference(label).area < AREA_TOL, (mask, 'relocated die ID')
    return label      # etched die-ID glyphs (um), truthy


def relabel_wafer_die_id(parent, output, settings):
    """Replace an upper-left positive-metal mask's etched die ID in place or in a new file."""
    from os.path import commonprefix

    from wafer_region_overlay import overlay_wafer_regions

    def label(text):
        glyphs = label_geometry(text, LABEL_SIZE_UM, 10, 400, 30)
        left, _, _, top = glyphs.bounds
        margin = settings['edge_inset_um']
        return translate(glyphs, -5000 + margin - left, 5000 - margin - top)

    original = label(settings['source_text'])
    replacement = label(settings['text'])
    window = original.union(replacement).envelope.buffer(5, join_style=2)
    prefix = commonprefix((settings['source_text'], settings['text']))
    if prefix:
        # Keep the existing etched prefix, including the renderer's rounded vertices.
        _, ymin, xmax, ymax = window.bounds
        window = box(label(prefix).bounds[2] + 5, ymin, xmax, ymax)
        original = original.intersection(window)
        replacement = replacement.intersection(window)

    def etched_region(path):
        layout = gdstk.read_gds(str(path))
        top, = layout.top_level()
        metal = [p for p in top.get_polygons() if (p.layer, p.datatype) == (1, 0)]
        index = STRtree([box(*low, *high) for low, high in (p.bounding_box() for p in metal)])
        return window.difference(region_in(window, metal, index))

    etched = etched_region(parent)
    if (etched.intersection(original).area < .85 * original.area
            or etched.difference(original.buffer(4)).area > AREA_TOL):
        raise ValueError(f'Existing etch does not match {settings["source_text"]} in {parent}')
    with TemporaryDirectory(dir=Path(output).parent) as temporary:
        staged = Path(temporary) / Path(output).name
        library = gdstk.Library(unit=1e-6, precision=1e-9)
        die = library.new_cell('DIE_ID_TEMPLATE')
        die.add(*gdstk.boolean(gdstk.rectangle((-5000, -5000), (5000, 5000)), polygons(replacement),
                               'not', precision=GRID, layer=1))
        library.write_gds(str(staged))
        patch = [round(coordinate / GRID) for coordinate in window.bounds]
        overlay_wafer_regions(parent, staged, [patch])
        new = etched_region(staged)
        if new.symmetric_difference(replacement).area > AREA_TOL:
            raise ValueError(f'Etched mask does not match {settings["text"]} after GDS readback')
        staged.replace(output)
    return patch


def coupon_shapes(variant, x, y, label_um):
    label_options = {key: f'{value}um' for key, value in label_um.items()}
    design = build_design(variant, fab_contacts=True, label_options=label_options)
    table = design.qgeometry.tables['poly']
    placed = {row.name: translate(scale(row.geometry, 1000, 1000, origin=(0, 0)), x, y)
              for row in table.itertuples()}
    return placed['pocket'], placed['letter'], (placed['pad_l'], placed['pad_r'])


def add_corner_coupons(source, output, placements, label_um):
    lib = gdstk.read_gds(str(source))
    if (lib.unit, lib.precision) != (1e-6, 1e-9):
        raise ValueError('Expected parent wafer GDS with µm units and 1 nm precision')
    top, = lib.top_level()
    metal = top.get_polygons()
    assert {(p.layer, p.datatype) for p in metal} == {(1, 0)}
    index = STRtree([box(*low, *high) for low, high in (p.bounding_box() for p in metal)])

    cuts, additions, checks = [], [], []
    for corner, sx, sy in CORNERS:
        x, y = placements[corner]['centre_um']
        variant = placements[corner]['variant']
        assert x * sx > 0 and y * sy > 0, corner
        pocket, letter, pads = coupon_shapes(variant, x, y, label_um)
        etch = pocket.union(letter)
        assert box(-4500, -4500, 4500, 4500).covers(pocket), corner
        assert box(-4700, -4700, 4700, 4700).covers(letter), corner
        keepout = etch.buffer(50)
        window = etch.envelope.buffer(65, join_style=2)
        existing = region_in(window, metal, index)
        voids = get_parts(window.difference(existing))
        holes = [hole for hole in voids if hole.intersects(keepout)]
        for hole in holes:
            x0, y0, x1, y1 = hole.bounds
            assert hole.area <= 30 and max(x1 - x0, y1 - y0) <= 6, (corner, 'existing feature in coupon keep-out')
            additions.extend(polygons(hole.difference(etch)))
        cuts.extend(polygons(etch))
        for pad in pads:
            additions.extend(polygons(pad))

        # Replace the entire old mark and a 5 µm margin with one metal polygon.
        # Keeping the original cross polygon would leave its outline visible in CAD.
        cross = box(4450, 4450, 4850, 4850)
        cross = scale(cross, sx, sy, origin=(0, 0))
        marker_void = cross.difference(region_in(cross, metal, index))
        had_cross = marker_void.area > cross.area / 2
        if had_cross:
            repair = cross.buffer(5, join_style=2)
            assert not repair.intersects(keepout), corner
            margin = repair.difference(cross)
            assert margin.difference(region_in(repair, metal, index)).area < AREA_TOL, corner
            cuts.extend(polygons(repair))
            additions.extend(polygons(repair))
        checks.append((corner, variant, x, y, pocket, letter, pads, cross, had_cross, window))

    # Keep the original positive metal and cut only polygons that meet an etch.
    cut_shapes = [Polygon(cut.points) for cut in cuts]
    hits = index.query(cut_shapes)
    by_polygon = {}
    for cut_i, metal_i in zip(*hits):
        by_polygon.setdefault(metal_i, []).append(cuts[cut_i])
    result = gdstk.Library(unit=lib.unit, precision=lib.precision)
    cell = result.new_cell('TOP')
    for i, poly in enumerate(metal):
        if i in by_polygon:
            cell.add(*gdstk.boolean(poly, by_polygon[i], 'not', precision=GRID, layer=1))
        else:
            cell.add(poly)
    cell.add(*additions)
    result.write_gds(str(output), max_points=199)
    return checks


def check_corners(filename, checks, *, cheesed=False, allow=None):
    """allow: etched geometry (the relocated die ID) that may overlap a removed corner mark."""
    lib = gdstk.read_gds(str(filename))
    top, = lib.top_level()
    metal = top.get_polygons()
    index = STRtree([box(*low, *high) for low, high in (p.bounding_box() for p in metal)])
    for corner, variant, x, y, pocket, letter, pads, cross, had_cross, window in checks:
        local = region_in(window, metal, index)
        etch = pocket.union(letter)
        pads_in_mask = local.intersection(pocket.buffer(-GRID))
        assert pads_in_mask.symmetric_difference(pads[0].union(pads[1])).area < AREA_TOL, corner
        assert len(get_parts(pads_in_mask)) == 2 and not local.covers(Point(x, y)), corner
        assert local.intersection(letter).area < AREA_TOL, (corner, 'letter etch')
        if had_cross and not cheesed:
            # An etched silhouette would fail the area check. A filled silhouette
            # made of touching polygons would fail the single-polygon check.
            void = cross.difference(region_in(cross, metal, index))
            if allow is not None:
                void = void.difference(allow.buffer(GRID * 10))
            assert void.area < AREA_TOL, (corner, 'old cross')
            interiors = [Polygon(metal[i].points).buffer(0).intersection(cross)
                         for i in index.query(cross)]
            covering = [part for part in interiors if part.area > AREA_TOL]
            if allow is not None and allow.intersects(cross):
                # The die-ID cut may split the fill into several polygons; none may share an edge (seam).
                target = cross.difference(allow.buffer(GRID * 10))
                assert union_all(covering).buffer(GRID).covers(target), (corner, 'old marker seam')
                for i, a in enumerate(covering):
                    for b in covering[i + 1:]:
                        assert a.boundary.intersection(b.boundary).length < GRID, (corner, 'old marker seam')
            else:
                assert len(covering) == 1 and covering[0].covers(cross), (corner, 'old marker seam')
        keepout = etch.buffer(50).difference(etch)
        assert keepout.difference(local).area < AREA_TOL, (corner, 'holes within keep-out')
        envelope = GEOMETRY['variants'][variant].get('array_envelope_um')
        if envelope:
            array = translate(box(*envelope), x, y)
            assert not local.intersects(array) and abs(array.distance(pocket.boundary) - 55) < GRID, corner


def review_layout_options(name, review):
    if not review:
        return {'readout': 'legacy'} if name == 'Design6' else {}
    if name == 'Design2':
        return {'tc_flux_offset_um': review['tc_flux_offset_um'],
                **review.get('layout_options', {})}
    from fluxline import B2_FLUX_OPTIONS

    flux = dict(B2_FLUX_OPTIONS, make_fbl=True, fbl_sep=f'{review["fbl_sep_um"]}um')
    if name == 'Design3':
        options = {'tmon_options': {'flux_bias_line_options': flux}, **review.get('layout_options', {})}
        if 'ground_tail_um' in review:
            options['tmon_options']['squid_ground_tail'] = f'{review["ground_tail_um"]}um'
        return options
    options = {'qflux': flux, 'coupler': {'flux_bias_line_options': flux}}
    if 'ground_tail_um' in review:
        options['coupler']['squid_ground_tail'] = f'{review["ground_tail_um"]}um'
        options['qsquid'] = {'ground_tail': f'{review["ground_tail_um"]}um'}
    # Design6: candidates up to r13 were built with the TapJunction + HammerCap readout; build_design defaults to
    # readout='purcell_resonator_pair' since 2026-11-19, so a candidate must name readout='purcell_resonator_pair' in layout_options to get PurcellReadout.
    return {**options, 'readout': 'legacy', **review.get('layout_options', {})}


def design6_readout_args(candidate, layout_options):
    """Legacy hammer / lengths from wafer_candidate only for the legacy readout; purcell_resonator_pair lengths come from build_design."""
    if layout_options.get('readout', 'legacy') == 'purcell_resonator_pair':
        return dict(readout_hammer=None, lengths=layout_options.pop('lengths', None))
    return dict(readout_hammer=candidate.get('readout_hammer_um'), lengths=candidate.get('readout_lengths_mm'))


def build_launch_mask(name, staging, *, readout_hammer=None, readout_lengths=None, review=None):
    """Rebuild the base fab mask with the Design's wafer launch anchors."""
    sys.path.insert(0, str(ROOT / name / 'layout'))
    raw = staging / 'launches_raw.gds'
    rounded = staging / 'launches_rounded.gds'
    mask = staging / 'launches_wafer.gds'
    layout_options = review_layout_options(name, review)
    if name == 'Design2':
        from layout3_full import export

        source = ROOT / name / 'gds' / 'previous' / 'layout3 ftf.gds'
        export(source, raw, fab=True, standard_launches=True, port_options=wafer_pcb_ports, **layout_options)
    elif name == 'Design4':
        from purcell_6500 import LENGTHS, build

        export_gds(build(LENGTHS, wafer=True), str(raw), fab=True, cheese=False)
    else:
        options = dict(jj_metal=False, fab_contacts=True)
        if name == 'Design6':
            from design6_layout import build_design as build_chip

            options.update(wafer_ports=True, **design6_readout_args(
                dict(readout_hammer_um=readout_hammer, readout_lengths_mm=readout_lengths), layout_options))
        else:
            options.update(include_tc_landing_pads=True, include_test_junctions=False)
            if name == 'Design3':
                from layout import build_design as build_chip

                options['rounded'] = False
            else:
                from ftf_layout import build_design as build_chip
        export_gds(build_chip(**options, **layout_options), str(raw), fab=True, cheese=False)
    if name in ('Design3', 'Design6'):
        sys.path.insert(0, str(ROOT / 'Design3' / 'layout' / 'components'))
        from round_mask import round_mask

        round_mask(raw, rounded)
        raw = rounded
    write_wafer_gds(raw, mask)
    return mask


def apply_regional_overlay(home, candidate, kind, output):
    overlay = candidate.get('regional_overlay')
    if overlay:
        from wafer_region_overlay import overlay_wafer_regions

        boxes = overlay.get(f'{kind}_boxes_nm') or overlay['boxes_nm']
        overlay_wafer_regions(home / overlay[f'parent_{kind}_gds'], output, boxes.values())


def verify_ground_strips(path, review):
    """Measure the ground strip from the nearest control-line slot to its neighbouring etch."""
    library = gdstk.read_gds(str(path))
    top, = library.top_level()
    metal = top.get_polygons(layer=1, datatype=0)
    index = STRtree([box(*low, *high) for low, high in (polygon.bounding_box() for polygon in metal)])
    widths = {}
    for device, inspection in review['ground_inspection'].items():
        window = box(*inspection['window_um'])
        etches = list(get_parts(window.difference(region_in(window, metal, index))))
        nearest_slot, *other_slots = [
            next(etch for etch in etches if etch.contains(Point(*probe)))
            for probe in inspection['slot_probes_um']
        ]
        if nearest_slot in other_slots:
            raise ValueError(f'{device}: the slots meet inside the inspection window')
        slots = [nearest_slot, *other_slots]
        width = min(nearest_slot.distance(etch) for etch in etches if etch not in slots)
        if width + 1e-9 < review['minimum_ground_width_um']:
            raise ValueError(f'{device}: ground strip is {width:.3f} µm in {path.name}')
        widths[device] = round(width, 3)
    if 'tail_inspection' in review:
        from shapely.geometry import LineString

        for device, inspection in review['tail_inspection'].items():
            window = box(*inspection['window_um'])
            etch = window.difference(region_in(window, metal, index))
            xmin, ymin, xmax, ymax = window.bounds
            center = inspection['center_x_um']
            floors = {}
            for offset in (-30, 0, 30):
                x = center + offset
                section = etch.intersection(LineString([(x, ymin), (x, ymax)]))
                if section.is_empty:
                    raise ValueError(f'Missing SQUID pocket section at x = {x} µm')
                floors[offset] = section.bounds[1]
            for side in (-30, 30):
                tail = floors[0] - floors[side]
                if abs(tail - 7.0) > GRID:
                    raise ValueError(f'{device} SQUID ground tail is {tail:.3f} µm in {path.name}')
            widths[device + '_ground_tail'] = round(floors[0] - floors[30], 3)
    return widths


def export_reference(name, candidate, target, *, review=None):
    """Export the unnotched or annotated companion from the selected wafer geometry."""
    home = ROOT / name
    sys.path.insert(0, str(home / 'layout'))
    with TemporaryDirectory(dir=target.parent) as temp:
        staging = Path(temp)
        raw, checked = staging / 'model.gds', staging / target.name
        layout_options = review_layout_options(name, review)
        if name == 'Design2':
            from layout3_full import export

            export(home / 'gds/previous/layout3 ftf.gds', raw,
                   standard_launches=True, port_options=wafer_pcb_ports, **layout_options)
        elif name == 'Design3':
            from export_gds import export_simulation

            export_simulation(raw, **layout_options)
        elif name == 'Design4':
            from purcell_6500 import LENGTHS, build

            export_gds(build(LENGTHS, wafer=True), str(raw))
        elif name == 'Design5':
            from chip6_d4_fab import build_fab

            source = home / (candidate.get('channel_lengths_source') or 'parameters/Design5_2026-11-13.json')
            lengths = json.loads(source.read_text(encoding='utf-8'))['channels']
            export_gds(build_fab(lengths=lengths, hammer=candidate['hammer_geometry_um']), str(raw))
        elif name == 'Design6':
            from design6_layout import build_design, export_simulation_gds

            design = build_design(wafer_ports=True, **design6_readout_args(candidate, layout_options), **layout_options)
            export_simulation_gds(design, raw)
        else:
            from ftf_layout import build_design

            design = build_design(jj_metal=False, include_tc_landing_pads=True, include_test_junctions=False)
            export_gds(design, str(raw))
        write_wafer_gds(raw, checked)
        die_id = candidate.get('die_id') or {}
        if review and review.get('regional_overlay'):
            apply_regional_overlay(home, review, 'sim', checked)
        else:
            if 'edge_inset_um' in die_id:
                move_wafer_die_id(checked, dict(die_id, text=die_id.get('source_text', die_id['text'])), required=False)
            apply_regional_overlay(home, candidate, 'sim', checked)
        ports = SIDE_PORTS if name in ('Design4', 'Design5') else wafer_pcb_ports
        verify_output(checked, name, False, ports_override=ports, expect_die_id='edge_inset_um' not in die_id)
        if 'edge_inset_um' in die_id and die_id.get('also_relocated_in_companion'):
            verify_design_id(checked, dict(die_id, text=die_id.get('source_text', die_id['text'])))
        if review:
            print(f'{target.name} ground widths (µm): {verify_ground_strips(checked, review)}')
        checked.replace(target)


def export_ground_review(name, target, *, reference_only, review_key):
    home = ROOT / name
    metadata = json.loads((home / 'design.json').read_text(encoding='utf-8'))
    candidate, review = metadata['wafer_candidate'], metadata[review_key]
    target = target or home / review['reference_gds' if reference_only else 'fabrication_gds']
    if target.exists():
        raise FileExistsError(target)
    if reference_only:
        export_reference(name, candidate, target, review=review)
        return
    with TemporaryDirectory(dir=target.parent) as temp:
        mask = build_launch_mask(name, Path(temp), readout_hammer=candidate.get('readout_hammer_um'),
                                 readout_lengths=candidate.get('readout_lengths_mm'), review=review)
        if name == 'Design2':
            # Recompute cheese eligibility around the displaced imported CPW.
            replace_wafer_cheese(mask, launches=wafer_pcb_ports, rounded=False)
        apply_regional_overlay(home, review, 'fab', mask)
        verify_output(mask, name, True, ports_override=wafer_pcb_ports,
                      expect_die_id='edge_inset_um' not in (candidate.get('die_id') or {}))
        print(f'{target.name} ground widths (µm): {verify_ground_strips(mask, review)}')
        mask.replace(target)


def verify_design_id(path, die_id, *, png=None):
    """Assert the etched design ID (die_id['text']) at the upper-left die margin of a mask; optional crop PNG.

    Expected glyphs: xmark.label_geometry(text, LABEL_SIZE_UM) placed edge_inset_um from the left and top die
    edges (the convention of move_wafer_die_id / relabel_wafer_die_id). Fails if a glyph is missing, an extra
    etch sits within 20 um of the text, or the text is not at the expected place.
    """
    text, margin = die_id['text'], die_id['edge_inset_um']
    glyphs = label_geometry(text, LABEL_SIZE_UM, 10, 400, 30)
    left, _, _, top = glyphs.bounds
    expected = translate(glyphs, -5000 + margin - left, 5000 - margin - top)
    window = expected.envelope.buffer(20, join_style=2)
    top_cell, = gdstk.read_gds(str(path)).top_level()
    metal = [p for p in top_cell.get_polygons() if (p.layer, p.datatype) == (1, 0)]
    index = STRtree([box(*low, *high) for low, high in (p.bounding_box() for p in metal)])
    etched = window.difference(region_in(window, metal, index))
    # Rounded masks (round_mask) round every glyph corner (about 1 % of the text area); one wrong
    # character, a wrong size or a wrong position changes far more than 2 % (checked: > 7 %).
    extra = etched.difference(expected).area
    missing = expected.difference(etched).area
    limit = 0.02 * expected.area
    if extra > limit or missing > limit:
        raise ValueError(f'{Path(path).name}: design ID {text!r} mismatch: extra etch {extra:.1f} um^2, '
                         f'missing {missing:.1f} um^2 (limit {limit:.0f})')
    if png:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        x0, y0, x1, y1 = window.buffer(60).bounds
        fig, ax = plt.subplots(figsize=(6, 6 * (y1 - y0) / (x1 - x0) + .6))
        shapes = region_in(box(x0, y0, x1, y1), metal, index)
        for g in getattr(shapes, 'geoms', [shapes]):
            ax.fill(*g.exterior.xy, color='#9aa3ad', lw=0)
            for h in g.interiors:
                ax.fill(*h.xy, color='white', lw=0)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect('equal')
        ax.set_title(f'{Path(path).stem}: design ID {text} (height {LABEL_SIZE_UM:g} um, inset {margin} um)',
                     fontsize=8)
        fig.tight_layout(); fig.savefig(png, dpi=130); plt.close(fig)
    print(f'{Path(path).name}: design ID {text!r} verified at x {expected.bounds[0]:.1f}, '
          f'top {expected.bounds[3]:.1f} um', flush=True)
    return expected


def export_full_candidate(name, review_key, *, reference_only=False):
    """Rebuild a complete wafer mask (and its companion) from the current QComponents.

    Uses the Design's wafer_candidate (launches, corner coupons, die ID, cheese) with the layout
    options of design.json[review_key]; unlike export_ground_review no parent GDS regions are pasted in.
    """
    home = ROOT / name
    metadata = json.loads((home / 'design.json').read_text(encoding='utf-8'))
    candidate, review = metadata['wafer_candidate'], metadata[review_key]
    if review.get('regional_overlay'):
        raise ValueError(f'{review_key}: a full rebuild must not use a regional overlay')
    output = review['reference_gds' if reference_only else 'fabrication_gds']
    if (home / output).exists():
        raise FileExistsError(home / output)
    clean = dict(candidate, regional_overlay=None)
    override = dict(regional_overlay=None, fabrication_gds=output)
    if review.get('die_id'):          # e.g. a new design ID for a new readout revision
        override['die_id'] = dict(candidate.get('die_id') or {}, **review['die_id'])
        clean['die_id'] = override['die_id']
    if reference_only:
        export_reference(name, clean, home / output, review=review)
    else:
        export_design(name, override=override, review=review)


def export_design(name, *, override=None, review=None):
    home = ROOT / name
    candidate = json.loads((home / 'design.json').read_text(encoding='utf-8'))['wafer_candidate']
    if override:
        candidate.update(override)
    if name == 'Design5':
        from shutil import copyfile

        source = home / candidate['fabrication_gds']
        target = home / candidate['delivery_fabrication_gds']
        verify_output(source, name, True, expect_die_id=bool(candidate.get('die_id')))
        target.parent.mkdir(parents=True, exist_ok=True)
        copyfile(source, target)
        assert sha256(source.read_bytes()).digest() == sha256(target.read_bytes()).digest()
        return
    target = home / candidate['fabrication_gds']
    target.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=target.parent) as temp:
        mask = Path(temp) / target.name
        launch_mask = build_launch_mask(name, Path(temp),
                                        readout_hammer=candidate.get('readout_hammer_um'),
                                        readout_lengths=candidate.get('readout_lengths_mm'), review=review)
        placements = candidate.get('corner_test_junctions')
        if placements:
            checks = add_corner_coupons(launch_mask, mask, placements, candidate['corner_label_um'])
        else:
            launch_mask.replace(mask)
        die_id = candidate.get('die_id') or {}
        moved_id = None
        if 'edge_inset_um' in die_id:
            moved_id = move_wafer_die_id(mask, dict(die_id, text=die_id.get('source_text', die_id['text'])))
        if placements:
            check_corners(mask, checks, allow=moved_id)
        ports = SIDE_PORTS if name == 'Design4' else wafer_pcb_ports
        perforations = replace_wafer_cheese(mask, launches=ports, rounded=name in ('Design3', 'Design6'))
        apply_regional_overlay(home, candidate, 'fab', mask)
        if 'source_text' in die_id:
            relabel_wafer_die_id(mask, mask, die_id)
        if placements:
            check_corners(mask, checks, cheesed=True)
        verify_output(mask, name, True, ports_override=ports, expect_die_id='edge_inset_um' not in die_id)
        if 'edge_inset_um' in die_id:
            fig = home / 'fig' / (target.stem + '_design_id.png')
            verify_design_id(mask, die_id, png=None if fig.exists() else fig)
        if review:
            print(f'{target.name} ground widths (µm): {verify_ground_strips(mask, review)}', flush=True)
        mask.replace(target)
    stage = ' before applying the regional overlay' if candidate.get('regional_overlay') else ''
    print(f'{name}: verified {perforations["new_holes"]} perforations on the 20 µm lattice{stage}', flush=True)
    if placements:
        print(f'{name}: ' + ', '.join(f'{corner}={item["variant"]}' for corner, item in placements.items()),
              flush=True)
    return perforations


if __name__ == '__main__':
    import geopandas as gpd

    if not hasattr(gpd.GeoDataFrame, 'append'):
        gpd.GeoDataFrame.append = gpd.GeoDataFrame._append
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--design', choices=DESIGNS, nargs='+', default=DESIGNS)
    parser.add_argument('--reference-only', action='store_true', help='Export the adopted local companion GDS')
    reviews = parser.add_mutually_exclusive_group()
    reviews.add_argument('--flux-ground-review', action='store_true',
                         help='Export the transmon ground-strip candidate for one Design')
    reviews.add_argument('--charge-ground-review', action='store_true',
                         help='Export the Design3 fluxonium charge-line ground-strip candidate')
    reviews.add_argument('--candidate', metavar='KEY',
                         help='Full rebuild with the layout options of design.json[KEY] (no regional overlay)')
    reviews.add_argument('--squid-tail-review', action='store_true',
                         help='Export the seven-micrometre SQUID ground-tail candidate')
    parser.add_argument('--output', type=Path, help='Write a single-design reproduction to this GDS path')
    args = parser.parse_args()
    selected = args.design
    if args.output and len(selected) != 1:
        parser.error('--output requires exactly one design')
    if args.flux_ground_review and (len(selected) != 1 or selected[0] not in ('Design2', 'Design3', 'Design6')):
        parser.error('--flux-ground-review requires one of Design2, Design3 or Design6')
    if args.charge_ground_review and selected != ['Design3']:
        parser.error('--charge-ground-review requires --design Design3')
    if args.squid_tail_review and (len(selected) != 1 or selected[0] not in ('Design2', 'Design3', 'Design6')):
        parser.error('--squid-tail-review requires one of Design2, Design3 or Design6')
    if len(selected) == 1:
        name = selected[0]
        target = args.output.resolve() if args.output else None
        if target:
            target.parent.mkdir(parents=True, exist_ok=True)
        if args.candidate:
            export_full_candidate(name, args.candidate, reference_only=args.reference_only)
        elif args.flux_ground_review:
            export_ground_review(name, target, reference_only=args.reference_only,
                                 review_key='transmon_flux_ground_candidate')
        elif args.charge_ground_review:
            export_ground_review(name, target, reference_only=args.reference_only,
                                 review_key='fluxonium_charge_ground_candidate')
        elif args.squid_tail_review:
            export_ground_review(name, target, reference_only=args.reference_only,
                                 review_key='squid_ground_tail_candidate')
        elif args.reference_only:
            home = ROOT / name
            candidate = json.loads((home / 'design.json').read_text(encoding='utf-8'))['wafer_candidate']
            output = candidate.get('reference_gds') or candidate['simulation_gds']
            export_reference(name, candidate, target or home / output)
        elif target and name == 'Design5':
            sys.path.insert(0, str(ROOT / name / 'layout'))
            from chip6_d4_fab import export_wafer

            export_wafer(target)
        else:
            export_design(name, override={'fabrication_gds': str(target)} if target else None)
    else:
        # Several designs use modules with the same name in their own layout directories.
        flags = ['--reference-only'] if args.reference_only else []
        for name in selected:
            subprocess.run([sys.executable, '-X', 'utf8', '-B', __file__, '--design', name, *flags], check=True)
