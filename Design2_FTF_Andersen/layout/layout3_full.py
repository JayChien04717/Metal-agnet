"""Recover Design2_FTF_Andersen planar conductors directly from the layout3 FTF mask."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
from project_paths import pattern_design_paths as _pattern_paths
from pathlib import Path
import argparse
import sys
import numpy as np
import gdstk
from shapely import STRtree, get_parts, set_precision
from shapely.affinity import rotate, scale, translate
from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import unary_union
from metal_width import MIN_METAL_WIDTH_UM

PORTS = {
    "readout_in": (-2000, 4450), "readout_out": (2000, 4450),
    "flux_left": (-4450, 2000), "flux_right": (4450, 2000),
    "charge_left": (-4450, -2000), "charge_right": (4450, -2000),
    "flux_tc": (2000, -4450), "grounded_spare": (-2000, -4450),
}
PCB_ASSIGNMENTS = {
    'readout_in': 'P8', 'readout_out': 'P7',
    'flux_left': 'P1', 'flux_right': 'P6',
    'charge_left': 'P2', 'charge_right': 'P5',
    'flux_tc': 'P4', 'grounded_spare': 'P3',
}
DEVICES = {
    "fx_left_top": (-485, -880), "fx_left_bottom": (-485, -1120),
    "fx_right_top": (485, -880), "fx_right_bottom": (485, -1120),
    "transmon": (0, -1400), "res_left": (-485, -640), "res_right": (485, -640),
}
CORNER_TESTS = box(-5000,3500,-3500,4950).union(box(3500,3500,5000,4950)).union(
    box(-5000,-4950,-3500,-3500)).union(box(3500,-4950,5000,-3500))
PLANE = box(-5000,-5000,5000,5000)
SQUID_CENTER_Y = -1568.0
DIE_ID = "FTFC00"

# Trim the original route up to a straight splice, then join it from the
# standard launch. Coordinates are in µm. The last waypoint is the splice.
WAFER_ROUTES = {
    'readout_in': ([(-1100, 3270), (-2750, 3270)], (-2750, 3150, 0, 4950)),
    'readout_out': ([(1100, 3270), (2750, 3270)], (0, 3150, 2750, 4950)),
    'flux_left': ([(-3670, 1600), (-3670, 1350)], (-4950, 1350, 0, 4950)),
    'flux_right': ([(3620, 1600), (3620, 1350)], (0, 1350, 4950, 4950)),
    'charge_left': ([(-515, -1600), (-515, -1425)], (-4950, -4950, 0, -1425)),
    'charge_right': ([(515, -1600), (515, -1425)], (0, -4950, 4950, -1425)),
    'flux_tc': ([(1100, -3620), (900, -3620)], (900, -4950, 4950, -3500)),
    'grounded_spare': ([], (-4950, -4950, 0, -3650)),
}


def launch_probes(standard_launches=False, *, port_options=None):
    if not standard_launches:
        return PORTS
    from pcbport import pcb_ports

    ports = pcb_ports if port_options is None else port_options
    return {name: tuple(float(ports[port][f'pos_{axis}'][:-2]) * 1000
                        for axis in ('x', 'y'))
            for name, port in PCB_ASSIGNMENTS.items()}


def replace_launch_pads(ground, sheets, *, standard_launches=False, cheese_tree=None, port_options=None):
    """Splice the shared launch into its original signal conductor."""
    sys.path[:0] = [p for p in _pattern_paths() if p not in sys.path]
    from bondpad import add_bond_pad
    from qiskit_metal import designs

    template = designs.DesignPlanar()
    ports = launch_probes(standard_launches, port_options=port_options)
    for name, (x, y) in PORTS.items():
        vertical = abs(y) > abs(x)
        outer = y if vertical else x
        orientation = (270 if outer > 0 else 90) if vertical else (180 if outer > 0 else 0)
        nominal = 4500 if outer > 0 else -4500
        destination = ports[name]
        transverse = destination[0 if vertical else 1]
        px, py = (transverse, nominal) if vertical else (nominal, transverse)
        grounded = name.startswith('flux_') or name == 'grounded_spare'
        width = 15
        pad = add_bond_pad(template, name, dict(pos_x=f'{px}um', pos_y=f'{py}um',
                                               orientation=orientation),
                           trace_width=f'{width}um', trace_gap='10um')
        table = template.qgeometry.tables['poly']
        rows = table.loc[table.component.eq(pad.id)].set_index('name').geometry
        pad_metal = set_precision(scale(rows['launch_pad'], 1000, 1000, origin=(0, 0)), .001)
        pad_cut = set_precision(scale(rows['pocket'], 1000, 1000, origin=(0, 0)), .001)
        tie = np.round(pad.pins['tie'].middle * 1000, 3)
        neck_len = abs(tie[1 if vertical else 0]) - 3790

        def place(shape):
            return set_precision(translate(rotate(shape, orientation, origin=(0, 0)), *tie), .001)

        original = next(sheet for sheet in sheets if sheet.covers(Point(x, y)))
        sheets.remove(original)
        delta = (x - px, 0) if vertical else (0, y - py)
        shift = round(rotate(Point(*delta), -orientation, origin=(0, 0)).y, 3)
        window = place(box(-800, shift - 425, neck_len, shift + 425))
        if standard_launches:
            points, bounds = WAFER_ROUTES[name]
            if points and not vertical:
                points = [(points[0][0], py), *points[1:]]
            trim = box(*bounds)
            removed = original.intersection(trim)
            # Close the former CPW slots without filling any unrelated etch.
            window = window.union(removed.buffer(14).intersection(trim))
            retained = original.difference(trim)
            path = np.array([tie, *points])

            def trace(cpw_width):
                if not points:
                    return Polygon()
                polygon, = gdstk.FlexPath(path, cpw_width, bend_radius=75,
                                          tolerance=.005).to_polygons()
                return set_precision(Polygon(polygon.points), .001)

            neck = trace(width)
            neck_cut = trace(width + 20)
        else:
            retained = original.difference(window)
            neck = place(box(0, -width / 2, neck_len, width / 2))
            neck_cut = place(box(0, -width / 2 - 10, neck_len, width / 2 + 10))
        replacement = pad_metal.union(neck)
        ground = ground.union(window.intersection(box(-4950, -4950, 4950, 4950)))
        launch_cut = pad_cut.union(neck_cut)
        if cheese_tree is not None:
            hits = cheese_tree.query(launch_cut, predicate='dwithin', distance=50)
            ground = ground.union(unary_union(cheese_tree.geometries[hits]))
        if standard_launches:
            # Before subtracting the new etch, require a clear corridor on
            # ground. This also catches a route crossing an earlier CPW slot.
            occupied = launch_cut.difference(ground)
            if occupied.area > .05 or any(launch_cut.intersection(sheet).area > .05 for sheet in sheets):
                raise ValueError(f'{name}: launch route hits existing geometry ({occupied.area:.3f} µm²)')
        ground = ground.difference(launch_cut)
        if grounded:
            ground = ground.union(replacement)
            if not retained.is_empty:
                sheets.append(retained)
        else:
            sheets.append(retained.union(replacement))
        # Check the pad, taper and neck against all positive-metal inputs.
        neck_probe = -15 if standard_launches and name == 'grounded_spare' else 20
        for along, offset in ((-480, 210), (-170, 140), (neck_probe, width / 2 + 5)):
            for side in (-1, 1):
                side_slot = place(Point(along, side * offset))
                if ground.covers(side_slot) or any(sheet.covers(side_slot) for sheet in sheets):
                    raise ValueError(f'{name}: blocked launch slot at {along}, {side * offset} um')
        if standard_launches:
            for begin, end in zip(path[:-1], path[1:]):
                centre = (begin + end) / 2
                dx, dy = (end - begin) / np.linalg.norm(end - begin)
                for side in (-1, 1):
                    slot = Point(centre + side * (width / 2 + 5) * np.array([-dy, dx]))
                    if ground.covers(slot) or any(sheet.covers(slot) for sheet in sheets):
                        raise ValueError(f'{name}: blocked route slot at {slot}')
    return ground, sheets

def move_tc_flux_slots(metal, offset_um):
    """Move the imported TC termination south and splice back into its CPW."""
    splice, shoulder = -1760, -1670
    window = box(-85, splice, 85, -1590.5)
    cuts = window.difference(metal)
    slots = [part for part in get_parts(cuts)
             if any(part.covers(Point(0, y)) for y in (-1595, -1614))]
    if len(slots) != 2 or any(part.bounds[0] <= -85 or part.bounds[2] >= 85 for part in slots):
        raise ValueError("Expected two TC flux slots in the source mask")
    slots = unary_union(slots)
    transition = slots.intersection(box(-85, splice, 85, shoulder))
    moved = slots.difference(transition)
    shifted_slots = unary_union([
        scale(transition, 1, 1 - offset_um / (shoulder - splice), origin=(0, splice)),
        translate(moved, yoff=-offset_um),
    ])
    return metal.union(slots).difference(shifted_slots)


# Fluxonia Q1 (left) / Q2 (right): shared FluxoniumDelft (pattern design/fluxonium/fluxonium_delft.py).
# Its defaults are the lab Q1; Q2 is the mirror image with the flux signal on the other U leg.
# array_length 54 um + fl_sep 15 um keeps the hairpin end at -55 um and U tip at -70 um (NTHU fl_sep convention).
FX_COMMON = dict(array_length='54um', flux_line_options=dict(fl_sep='15um'))
FX_OPTIONS = {
    'left': dict(FX_COMMON, pos_x='-485um', pos_y='-1000um'),
    'right': dict(FX_COMMON, pos_x='485um', pos_y='-1000um', mirror='True',
                  flux_line_options=dict(fl_sep='15um', signal_leg='bottom')),
}
# Pocket region shared with the TC coupler crescent (local Q1 frame, um): kept from the source mask.
FX_TC_ZONE = box(50, -275, 450, -95)


def fx_design(fab=False, jj_metal=False):
    """Design with Q1 / Q2 as FluxoniumDelft (fab: FluxoniumDelftFab, contact notches)."""
    sys.path[:0] = [p for p in _pattern_paths() if p not in sys.path]
    from qiskit_metal import designs
    from fluxonium_delft import FluxoniumDelft
    from fluxonium_delft_fab import FluxoniumDelftFab

    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.size.update(size_x='10mm', size_y='10mm')
    cls = FluxoniumDelftFab if fab else FluxoniumDelft
    for side, options in FX_OPTIONS.items():
        cls(design, f'fx_{side}', options=dict(options, jj_metal=str(bool(jj_metal) and not fab)))
    return design


def fx_geometry(fab=False):
    """(metal, cut, keep) of both fluxonia in um: component metal, component etch (pockets and
    line gaps) and the part of the etch whose source-mask content is kept (TC crescents)."""
    table = fx_design(fab).qgeometry.tables['poly']
    um = [scale(g, 1000, 1000, origin=(0, 0)) for g in table.geometry]
    metal = unary_union([g for g, s in zip(um, table['subtract']) if not s])
    cut = unary_union([g for g, s in zip(um, table['subtract']) if s])
    zones = [translate(FX_TC_ZONE, -485, -1000), translate(scale(FX_TC_ZONE, -1, 1, origin=(0, 0)), 485, -1000)]
    keep = unary_union(zones).intersection(cut).difference(metal.buffer(5))
    return set_precision(metal, .001), set_precision(cut, .001), keep


def use_fluxonium_delft(joined, fab=False):
    """Replace the source-mask fluxonia (Manhattan JJ stubs included) by FluxoniumDelft."""
    metal, cut, keep = fx_geometry(fab)
    return joined.difference(cut).union(joined.intersection(keep)).union(metal)


def recover(source, fab=False, *, standard_launches=False, port_options=None, tc_flux_offset_um=0,
            squid_ground_tail_um=None):
    lib = gdstk.read_gds(source)
    source_ground = unary_union([Polygon(p.points).buffer(0) for p in lib["TOP_MAIN_1_CHEESE_100"].polygons])
    # Holes smaller than 30 um² are the reference cheese pattern.
    cuts = unary_union([Polygon(r) for r in source_ground.interiors if Polygon(r).area >= 30])
    ground = box(*source_ground.bounds).difference(cuts).union(CORNER_TESTS.intersection(box(*source_ground.bounds)))
    cheese_tree = None
    if fab:
        cheese_holes = [Polygon(r) for r in source_ground.interiors if Polygon(r).area < 30]
        ground = ground.difference(unary_union(cheese_holes))
        if standard_launches:
            cheese_tree = STRtree(cheese_holes)
    metal = unary_union([Polygon(p.points).buffer(0) for p in lib["TOP"].polygons if p.layer == 5])
    sheets = [p for p in get_parts(metal) if p.area > 1000 and not p.intersects(CORNER_TESTS)
              and PLANE.covers(p)]
    # The FX contact notches (fab) come from FluxoniumDelftFab in use_fluxonium_delft() below.
    tc = next(p for p in sheets if p.covers(Point(*DEVICES["transmon"])))
    # The source GDS cuts the two coupler crescents around the FX bottom pads with a straight edge at
    # x = +/-429.5 um, leaving tapering tips narrower than 4 um (skill/chip_drc.md C12a). Remove only
    # the metal of those tips that a 4 um opening cannot reach; the rest of the crescent is unchanged.
    tips = unary_union([box(sx * 400, -1230, sx * 445, -1090) for sx in (-1, 1)])
    opened = tc.buffer(-(MIN_METAL_WIDTH_UM / 2 - 1e-3), quad_segs=16).buffer(MIN_METAL_WIDTH_UM / 2 - 1e-3, quad_segs=16)
    tc_tips = tc.intersection(tips).difference(opened)
    if not 0 < tc_tips.area < 50:
        raise ValueError(f'Unexpected TC crescent tip trim: {tc_tips.area:.2f} um^2')
    sheets.remove(tc)
    tc = tc.difference(tc_tips)
    sheets.append(tc)
    cy = SQUID_CENTER_Y
    aperture = box(-10,cy-12,10,cy+12)
    island = box(-20,cy+2.5,20,cy+23).difference(aperture)
    landing = box(-20,cy-23,20,cy-2.5).difference(aperture)
    if fab:
        for x in (-15,15):
            island = island.difference(box(x-1,cy+2.5,x+1,cy+6.5))
            landing = landing.difference(box(x+1,cy-6.5,x+3,cy-2.5))
    # Remove legacy stem metal where it would bridge the rectangular aperture.
    sheets.remove(tc)
    sheets.append(tc.union(island).difference(aperture).difference(box(-20,cy-2.5,20,cy+2.5)))
    ground = ground.union(landing).difference(aperture)
    ground, sheets = replace_launch_pads(ground, sheets, standard_launches=standard_launches,
                                        cheese_tree=cheese_tree, port_options=port_options)
    if standard_launches:
        # Extend the continuous positive ground to the die edge on layer 1.
        ground = ground.union(PLANE.difference(box(-4950, -4950, 4950, 4950)))
    ports = launch_probes(standard_launches, port_options=port_options)
    joined = use_fluxonium_delft(unary_union([ground,*sheets]), fab)
    if tc_flux_offset_um:
        joined = move_tc_flux_slots(joined, tc_flux_offset_um)
    if squid_ground_tail_um is not None:
        if not 0 < squid_ground_tail_um <= 10:
            raise ValueError('SQUID ground tail must be between 0 and 10 µm')
        # The loop floor is cy - 12. Restore ground up to the requested
        # pocket edge across both sides of the ground-side landing.
        edge = cy - 12 - squid_ground_tail_um
        joined = joined.union(box(-37.5, cy - 22, 37.5, edge))
    parts = list(get_parts(joined))
    ground = max(parts,key=lambda p:p.area)
    nets = {}
    probes = dict(DEVICES, feed=ports["readout_in"],
                charge_left=ports["charge_left"],charge_right=ports["charge_right"])
    for name, xy in probes.items():
        net, = [p for p in parts if p.covers(Point(*xy))]
        if net is ground or any(net.equals(old) for old in nets.values()):
            raise ValueError(f"Shorted net {name}")
        nets[name] = net
    if not nets["feed"].covers(Point(*ports["readout_out"])):
        raise ValueError("Disconnected feedline")
    for name in ("flux_left","flux_right","flux_tc","grounded_spare"):
        if not ground.covers(Point(*ports[name])):
            raise ValueError(f"Lost grounded return at {name}")
    for p in parts:
        if p is not ground and all(not p.equals(v) for v in nets.values()):
            print("Unmapped",p.area,p.bounds,flush=True)
    if len(nets)+1 != len(parts):
        raise ValueError(f"Unmapped metal islands: {len(parts)-len(nets)-1}")
    return ground,nets

def die_id_geometry():
    """Etched ID at the Design3_FTF_NTHU lower-right label position, in µm."""
    sys.path[:0] = [p for p in _pattern_paths() if p not in sys.path]
    from xmark import LABEL_SIZE_UM, label_geometry

    label = label_geometry(DIE_ID, LABEL_SIZE_UM, 10, 400, 30)
    return translate(label, 4650 - label.bounds[0] - label.bounds[2], -4650)


def export(source, output, fab=False, *, standard_launches=False, port_options=None, tc_flux_offset_um=0,
           squid_ground_tail_um=None):
    ground,nets = recover(source, fab, standard_launches=standard_launches, port_options=port_options,
                          tc_flux_offset_um=tc_flux_offset_um, squid_ground_tail_um=squid_ground_tail_um)
    lib = gdstk.Library(unit=1e-6,precision=1e-9)
    top = lib.new_cell("TOP")
    label = die_id_geometry()
    glyphs = []
    for glyph in get_parts(label):
        outer = gdstk.Polygon(np.asarray(glyph.exterior.coords), layer=1)
        counters = [gdstk.Polygon(np.asarray(r.coords)) for r in glyph.interiors]
        glyphs.extend(gdstk.boolean(outer, counters, "not", layer=1, precision=.001))
    left, bottom, right, top_edge = label.bounds
    for poly in [ground,*nets.values()]:
        outline = gdstk.Polygon(np.asarray(poly.exterior.coords),layer=1)
        holes = [gdstk.Polygon(np.asarray(r.coords)) for r in poly.interiors]
        fragments = gdstk.boolean(outline, holes, "not", layer=1, precision=.001)
        for fragment in (part for polygon in fragments for part in polygon.fracture(max_points=199, precision=.001)):
            (x0,y0),(x1,y1) = fragment.bounding_box()
            if x1 < left or x0 > right or y1 < bottom or y0 > top_edge:
                top.add(fragment)
            else:
                top.add(*gdstk.boolean(fragment,glyphs,"not",layer=1,precision=.001))
    lib.write_gds(output,max_points=199)

def build_design(source, readout_ports=False, *, tc_flux_offset_um=0, squid_ground_tail_um=None):
    from qiskit_metal import designs
    from qiskit_metal.qlibrary.core import QComponent
    ground,nets = recover(source, tc_flux_offset_um=tc_flux_offset_um,
                          squid_ground_tail_um=squid_ground_tail_um)

    class Conductor(QComponent):
        def make(self):
            self.add_qgeometry("poly",{"metal":scale(nets[self.name],.001,.001,origin=(0,0))})
            if readout_ports and self.name == "feed":
                for side, x in (("left", -2), ("right", 2)):
                    self.add_pin(side, np.array([[x,4.849],[x,4.85]]), .4,
                                 gap=.05, input_as_norm=True)

    class DielectricCuts(QComponent):
        def make(self):
            cuts = PLANE.difference(ground)
            for i,poly in enumerate(get_parts(cuts)):
                # Boolean subtraction in AEDT expects simple outlines. Retain
                # ground islands as holes in the etched region.
                self.add_qgeometry("poly",{f"cut_{i}":scale(poly,.001,.001,origin=(0,0))},subtract=True)

    design = designs.DesignPlanar()
    design.chips.main.material = "silicon"
    design.chips.main.size.update(size_x="10mm",size_y="10mm",size_z="-525um")
    for name in nets:
        Conductor(design,name)
    DielectricCuts(design,"mask")
    return design




def build_linear_design(source, readout_ports=False, *, tc_flux_offset_um=0, squid_ground_tail_um=None):
    """Add the explicit fluxonium arrays and SQUID inductors for HFSS."""
    d = build_design(source, readout_ports, tc_flux_offset_um=tc_flux_offset_um,
                     squid_ground_tail_um=squid_ground_tail_um)
    # Add the linearized junctions to the driven device geometry so the field
    # model includes both fluxonia and the SQUID.
    # Leads, hairpin end and junction lines: FluxoniumDelft with jj_metal=True (same names as before).
    fx = fx_design(jj_metal=True)
    fx_poly, fx_junction = fx.qgeometry.tables['poly'], fx.qgeometry.tables['junction']
    for side in ('left', 'right'):
        cid = fx.components[f'fx_{side}'].id
        rows = fx_poly[fx_poly.component.eq(cid)].set_index('name').geometry
        for part, pad in (('top', 'pad_top'), ('bottom', 'pad_bot')):
            comp=d.components[f'fx_{side}_{part}']
            table=d.qgeometry.tables['poly']
            idx,=table.index[table.component.eq(comp.id)]
            table.at[idx,'geometry']=table.at[idx,'geometry'].union(rows[pad])   # stem + lead
        component=d.components[f'fx_{side}_top']
        component.add_qgeometry('poly',{'array_end':rows['array_end']})
        jrows = fx_junction[fx_junction.component.eq(cid)].set_index('name')
        for name, source_name in (('jj', 'rect_jj'), ('array_top', 'array_top'), ('array_bottom', 'array_bot')):
            row = jrows.loc[source_name]
            component.add_qgeometry('junction', {name: row.geometry}, width=row.width,
                                    hfss_inductance=row.hfss_inductance)
    for label,x in (('left',-15),('right',15)):
        y=SQUID_CENTER_Y
        d.components['transmon'].add_qgeometry('junction',{
            f'jj_{label}':LineString(np.array([(x,y+2.5),(x,y-2.5)])/1000)},
            width=.007,hfss_inductance='50nH')
    return d


def main():
    root = _P('Design2_FTF_Andersen')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=root / "gds/previous/layout3 ftf.gds")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--fab", action="store_true", help="Keep the source cheese and add SQUID recesses")
    parser.add_argument('--standard-launches', action='store_true', help='Route to the shared 8-port PCB anchors')
    args = parser.parse_args()
    if args.standard_launches and args.output is None:
        parser.error('--standard-launches requires an explicit output path')
    default_name = "ftf_layout_fab.gds" if args.fab else "ftf_layout.gds"
    output = args.output or _relocate(root / "gds" / default_name)
    if args.standard_launches:
        from wafer_pcbport import wafer_pcb_ports

        export(args.reference, output, args.fab, standard_launches=True, port_options=wafer_pcb_ports)
    else:
        export(args.reference, output, args.fab)
    if args.output is None and not args.fab:
        from gds_tools import plot_gds
        plot_gds(str(output), str(_relocate(root / 'fig/layout/ftf_layout_chip.png')))


if __name__ == "__main__":
    main()
