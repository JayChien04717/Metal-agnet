"""Recover Design2 planar conductors directly from the layout3 FTF mask."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[2]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
from pathlib import Path
import argparse
import sys
import numpy as np
import gdstk
from shapely import get_parts, set_precision
from shapely.affinity import rotate, scale, translate
from shapely.geometry import Point, Polygon, box
from shapely.ops import unary_union

PORTS = {
    "readout_in": (-2000, 4450), "readout_out": (2000, 4450),
    "flux_left": (-4450, 2000), "flux_right": (4450, 2000),
    "charge_left": (-4450, -2000), "charge_right": (4450, -2000),
    "flux_tc": (2000, -4450), "grounded_spare": (-2000, -4450),
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


def replace_launch_pads(ground, sheets):
    """Splice the shared pad into the straight neck at 3790 µm from die centre."""
    shared = _P('pattern design')
    if str(shared) not in sys.path:
        sys.path.insert(0, str(shared))
    from bondpad import add_bond_pad
    from qiskit_metal import designs

    template = designs.DesignPlanar()
    for name, (x, y) in PORTS.items():
        vertical = abs(y) > abs(x)
        outer = y if vertical else x
        transverse = x if vertical else y
        orientation = (270 if outer > 0 else 90) if vertical else (180 if outer > 0 else 0)
        nominal = 4500 if outer > 0 else -4500
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
        neck = place(box(0, -width / 2, neck_len, width / 2))
        neck_cut = place(box(0, -width / 2 - 10, neck_len, width / 2 + 10))
        window = place(box(-800, -425, neck_len, 425))
        replacement = pad_metal.union(neck)
        ground = ground.union(window.intersection(box(-4950, -4950, 4950, 4950)))
        ground = ground.difference(pad_cut.union(neck_cut))
        probe = Point(x, y)
        original = next(sheet for sheet in sheets if sheet.covers(probe))
        sheets.remove(original)
        retained = original.difference(window)
        if grounded:
            ground = ground.union(replacement)
            sheets.append(retained)
        else:
            sheets.append(retained.union(replacement))
        # Check the pad, taper and neck against all positive-metal inputs.
        for along, offset in ((-480, 210), (-170, 140), (20, width / 2 + 5)):
            for side in (-1, 1):
                side_slot = place(Point(along, side * offset))
                if ground.covers(side_slot) or any(sheet.covers(side_slot) for sheet in sheets):
                    raise ValueError(f'{name}: blocked launch slot at {along}, {side * offset} um')
    return ground, sheets

def recover(source, fab=False):
    lib = gdstk.read_gds(source)
    source_ground = unary_union([Polygon(p.points).buffer(0) for p in lib["TOP_MAIN_1_CHEESE_100"].polygons])
    # Holes smaller than 30 um² are the reference cheese pattern.
    cuts = unary_union([Polygon(r) for r in source_ground.interiors if Polygon(r).area >= 30])
    ground = box(*source_ground.bounds).difference(cuts).union(CORNER_TESTS.intersection(box(*source_ground.bounds)))
    if fab:
        ground = ground.difference(unary_union([Polygon(r) for r in source_ground.interiors if Polygon(r).area < 30]))
    metal = unary_union([Polygon(p.points).buffer(0) for p in lib["TOP"].polygons if p.layer == 5])
    sheets = [p for p in get_parts(metal) if p.area > 1000 and not p.intersects(CORNER_TESTS)
              and PLANE.covers(p)]
    if fab:
        # The two reference stems end on the same y coordinates on both FX devices.
        for x in (-485, 485):
            for part, y0, y1 in (("top", -985, -981), ("bottom", -1019, -1015)):
                stem = next(p for p in sheets if p.covers(Point(x, y1 + 5 if part == "top" else y0 - 5)))
                hook = box(x - 1, y0, x + 1, y1)
                if not stem.covers(hook):
                    raise ValueError(f"Hook outside FX stem at {x}, {part}")
                sheets.remove(stem)
                sheets.append(stem.difference(hook))
    tc = next(p for p in sheets if p.covers(Point(*DEVICES["transmon"])))
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
    ground, sheets = replace_launch_pads(ground, sheets)
    joined = unary_union([ground,*sheets])
    parts = list(get_parts(joined))
    ground = max(parts,key=lambda p:p.area)
    nets = {}
    probes = dict(DEVICES, feed=PORTS["readout_in"],
                charge_left=PORTS["charge_left"],charge_right=PORTS["charge_right"])
    for name, xy in probes.items():
        net, = [p for p in parts if p.covers(Point(*xy))]
        if net is ground or any(net.equals(old) for old in nets.values()):
            raise ValueError(f"Shorted net {name}")
        nets[name] = net
    if not nets["feed"].covers(Point(*PORTS["readout_out"])):
        raise ValueError("Disconnected feedline")
    for name in ("flux_left","flux_right","flux_tc","grounded_spare"):
        if not ground.covers(Point(*PORTS[name])):
            raise ValueError(f"Lost grounded return at {name}")
    for p in parts:
        if p is not ground and all(not p.equals(v) for v in nets.values()):
            print("Unmapped",p.area,p.bounds,flush=True)
    if len(nets)+1 != len(parts):
        raise ValueError(f"Unmapped metal islands: {len(parts)-len(nets)-1}")
    return ground,nets

def die_id_geometry():
    """Etched ID at the Design3 lower-right label position, in µm."""
    shared_patterns = _P('pattern design')
    if str(shared_patterns) not in sys.path:
        sys.path.insert(0, str(shared_patterns))
    from xmark import label_geometry

    label = label_geometry(DIE_ID, 100, 10, 400, 30)
    return translate(label, 4650 - label.bounds[0] - label.bounds[2], -4650)


def export(source, output, fab=False):
    ground,nets = recover(source,fab)
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

def build_design(source, readout_ports=False):
    from qiskit_metal import designs
    from qiskit_metal.qlibrary.core import QComponent
    ground,nets = recover(source)

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




def build_linear_design(source, readout_ports=False):
    """Add the explicit fluxonium arrays and SQUID inductors for HFSS."""
    d = build_design(source, readout_ports)
    # Add the linearized junctions to the driven device geometry so the field
    # model includes both fluxonia and the SQUID.
    from shapely.geometry import LineString,box
    from shapely.affinity import scale
    for side,sign in (('left',-1),('right',1)):
        x=sign*485
        arm=sign*539
        for part,y0,y1 in (('top',-999,-984.75),('bottom',-1015.25,-1001)):
            comp=d.components[f'fx_{side}_{part}']
            table=d.qgeometry.tables['poly']
            idx,=table.index[table.component.eq(comp.id)]
            lead=scale(box(x-1,y0,x+1,y1),.001,.001,origin=(0,0))
            table.at[idx,'geometry']=table.at[idx,'geometry'].union(lead)
        component=d.components[f'fx_{side}_top']
        component.add_qgeometry('poly',{'array_end':scale(box(arm-1,-1005,arm+1,-995),.001,.001,origin=(0,0))})
        for name,poles,inductance in [('jj',[(x,-999),(x,-1001)],'34.38nH'),
            ('array_top',[(x+sign,-996),(arm-sign,-996)],'100nH'),
            ('array_bottom',[(x+sign,-1004),(arm-sign,-1004)],'100nH')]:
            component.add_qgeometry('junction',{name:LineString(np.array(poles)/1000)},
                                    width=.002,hfss_inductance=inductance)
    for label,x in (('left',-15),('right',15)):
        y=SQUID_CENTER_Y
        d.components['transmon'].add_qgeometry('junction',{
            f'jj_{label}':LineString(np.array([(x,y+2.5),(x,y-2.5)])/1000)},
            width=.007,hfss_inductance='50nH')
    return d


def main():
    root = _P('Design2')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=_relocate(root / "gds/layout3 ftf.gds"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--fab", action="store_true", help="Keep the source cheese and add SQUID recesses")
    args = parser.parse_args()
    default_name = "ftf_layout_fab.gds" if args.fab else "ftf_layout.gds"
    output = args.output or _relocate(root / "gds" / default_name)
    export(args.reference, output, args.fab)
    if args.output is None and not args.fab:
        from gds_tools import plot_gds
        plot_gds(str(output), str(_relocate(root / 'fig/layout/ftf_layout_chip.png')))


if __name__ == "__main__":
    main()
