"""(Copied from Design3e so Design3_3 runs on its own.) Q3D tuning of the lower-pad FX--TC coupling on a copy of Design3_2."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse
import json
from math import cos, pi, sin
from qiskit_metal import Dict
from pathlib import Path
import sys

from shapely.affinity import scale
from shapely.geometry import LineString, box
from shapely.ops import unary_union

HERE = _P('Design3_3')
ROOT = _P('.')
sys.path[:0] = [_I('Design3_2/simulation'), _I('Design3_2/layout'), _I('Design3_2/layout/components')]
import run_verification as base
import ftf_layout as layout
from my_tmon import MyTmonCoupler
from my_tmon_fab import MyTmonCouplerFab

def curved_arm(p, pad_edge, pad_half_height, gap, length):
    """Sweep a constant-width arm around the inner corner of the FX pad."""
    half_width = p.pad_arm_width / 2
    corner_radius = gap + half_width
    side_x = pad_edge - corner_radius
    approach_radius = .022
    start_x = side_x - approach_radius
    tip_x = pad_edge + length - half_width
    if start_x <= p.pad_arm_length / 2 - .02 or tip_x <= pad_edge:
        raise ValueError('Curved arm needs room for its bend and rounded tip.')

    def arc(cx, cy, radius, a, b):
        return [(cx + radius * cos(a + (b-a)*i/32),
                 cy + radius * sin(a + (b-a)*i/32)) for i in range(33)]

    centerline = [(p.pad_arm_length / 2 - .018, 0)]
    centerline += arc(start_x, approach_radius, approach_radius, -pi/2, 0)
    centerline += arc(pad_edge, pad_half_height, corner_radius, pi, pi/2)
    centerline += [(tip_x, pad_half_height + corner_radius)]
    return LineString(centerline).buffer(half_width, quad_segs=24)


def build_variant(gaps, lengths, rounded=False, arm_style='e4', arm_builder=None, **options):
    if any(gap < 2 or length <= 0 for gap,length in zip(gaps,lengths)):
        raise ValueError('Use at least 2 µm gap with positive overlap on both sides.')
    options = dict(options)
    original = layout.MyTmonCoupler
    original_fab = layout.MyTmonCouplerFab
    if arm_builder is None and arm_style == 'curved':
        arm_builder = curved_arm
    def modified(parent):
        class Coupler(parent):
            def make_pocket(self):
                super().make_pocket()
                p = self.p
                fingers = []
                for side, gap, length in zip((-1, 1), gaps, lengths):
                    gap, length = gap/1000, length/1000
                    outer = layout.FX_X - layout._mm(layout.FX_OPTIONS['pad_width'])/2
                    half_height = layout._mm(layout.FX_OPTIONS['pad_height'])/2
                    if arm_builder is not None:
                        part = arm_builder(p, outer, half_height, gap, length)
                    else:
                        edge = p.pad_arm_length/2+p.paddle_width
                        neck = box(edge-.015, -.014, outer-gap, half_height+gap+.012)
                        tip = box(outer-gap-.008, half_height+gap, outer+length, half_height+gap+.012)
                        part = unary_union([neck, tip])
                    fingers.append(part if side==1 else scale(part,-1,1,origin=(0,0)))
                placed = self._place(fingers)
                table = self.design.qgeometry.tables['poly']
                for name in ('tmon','pad_etch'):
                    idx, = table.index[(table.component == self.id) & (table.name == name)]
                    original_geometry = table.at[idx,'geometry']
                    if name == 'tmon' and arm_builder is not None:
                        # Drop the original flat paddles, retaining the bar and SQUID stem.
                        x, y = p.pos_x, p.pos_y
                        half_length = p.pad_arm_length / 2
                        original_geometry = original_geometry.intersection(
                            box(x-half_length,y-1,x+half_length,y+1)).buffer(0)
                    join_style = 1 if arm_builder is not None else 2
                    addition = placed if name == 'tmon' else [g.buffer(.016,join_style=join_style) for g in placed]
                    table.at[idx,'geometry'] = unary_union([original_geometry,*addition]).simplify(1e-9)
        return Coupler
    layout.MyTmonCoupler = modified(MyTmonCoupler)
    layout.MyTmonCouplerFab = modified(MyTmonCouplerFab)
    try:
        design = (base.rounded_design if rounded else layout.build_design)(**options)
        for name in ('res_TC','tee_TC','feed_FX1_TC','feed_TC_FX2'):
            design.delete_component(name)
        table = design.qgeometry.tables['poly']
        mask = (table.component == design.components['TC'].id) & table.name.isin(('readout_line','readout_window'))
        table.drop(index=table.index[mask], inplace=True)
        design.components['TC'].pins.pop('readout_line', None)
        layout.RouteStraight(design, 'feed_FX1_FX2', Dict(chip='main', pin_inputs=Dict(
            start_pin=Dict(component='tee_FX1',pin='prime_end'),
            end_pin=Dict(component='tee_FX2',pin='prime_start'))))
        assert 'res_TC' not in design.components and 'tee_TC' not in design.components
        assert 'readout_line' not in design.components['TC'].pins
        assert 'feed_FX1_FX2' in design.components
        return design
    finally:
        layout.MyTmonCoupler = original
        layout.MyTmonCouplerFab = original_fab

def verify_layout(design):
    table = design.qgeometry.tables['poly']
    def geom(q, name):
        row, = table.loc[(table.component == design.components[q].id) & (table.name == name)].itertuples()
        return row.geometry
    island = geom('TC','tmon')
    assert island.geom_type == 'Polygon' and island.is_valid
    for q in ('FX1','FX2'):
        charge_gap = island.distance(geom(q,'charge_line'))*1000
        print(q, 'TC electrode to charge line minimum gap um', charge_gap, flush=True)
        assert charge_gap >= 6
        for pad in ('pad_top','pad_bot'):
            gap = island.distance(geom(q,pad))*1000
            print(q, pad, 'minimum gap um', gap, flush=True)
            assert gap > 0.1
    etch = unary_union(list(table.loc[table["subtract"],'geometry']))
    assert island.difference(etch).area < 1e-9, 'TC island crosses ground'
    assert island.distance(geom('TC','squid_ground')) > 0
    for row in table.loc[~table["subtract"]].itertuples():
        if row.component != design.components['TC'].id or row.name != 'tmon':
            assert island.intersection(row.geometry).area < 1e-10, row.name
    print('Planar overlap check passed',flush=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tag')
    parser.add_argument('--project', type=Path, default=Path.home()/'AppData/Local/Temp/D3E/D3E.aedt')
    parser.add_argument('--gap', type=float, nargs=2, default=(4,4))
    parser.add_argument('--arm-style', choices=('e4','curved'), default='e4')
    parser.add_argument('--length', type=float, nargs=2, default=(23.5,45))
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--local',action='store_true')
    parser.add_argument('--charge-sep', type=float, nargs=2, default=(9,8))
    parser.add_argument('--stem-width', type=float, default=15)
    parser.add_argument('--tc-etch', type=float, default=55)
    parser.add_argument('--passes',type=int,default=12)
    parser.add_argument('--min-passes',type=int,default=4)
    parser.add_argument('--error',type=float,default=.5)
    parser.add_argument('--order', choices=('Low','Medium','High'), default='Medium')
    args=parser.parse_args()
    def build(**kw):
        tmon = dict(kw.pop('tmon_options', None) or {})
        tmon.update(pad_head_width=f'{args.stem_width}um', paddle_gap=f'{args.tc_etch}um')
        for q, sep in zip(('FX1','FX2'),args.charge_sep):
            key = q.lower()+'_options'
            fx = dict(kw.pop(key,None) or {})
            fx['charge_line_options'] = Dict(layout.FX_OPTIONS['charge_line_options'])
            fx['charge_line_options']['cl_sep'] = f'{sep}um'
            kw[key] = fx
        return build_variant(args.gap,args.length,rounded=True,arm_style=args.arm_style,tmon_options=tmon,**kw)
    verify_layout(build(for_process=True,fx_options={'array_end_metal':'False'}))
    if args.check:
        return
    base.RESULTS=base.q3d.RESULTS=_P('Design3_3/results/general/legacy')
    base.RESULTS.mkdir(exist_ok=True)
    base.save_aedt.use_project(args.project)
    base.q3d.build_design=build
    class Runner(base.q3d.Q3DRunner):
        def _prepare(self, design):
            lom=super()._prepare(design)
            base.join_squid_ground(lom.sim.renderer,merge_readouts=not args.local)
            renderer = lom.sim.renderer
            def analyze(name):
                # Set the native CG order explicitly. The renderer updates only pass counts and error.
                setup = renderer.pinfo.get_setup(name)
                setup.parent._setup_module.EditSetup(name, [
                    f'NAME:{name}', 'AdaptiveFreq:=','5GHz', 'SaveFields:=',False,'Enabled:=',True,
                    ['NAME:Cap','MaxPass:=',args.passes,'MinPass:=',args.min_passes,
                     'MinConvPass:=',2,'PerError:=',args.error,'PerRefine:=',25,
                     'AutoIncreaseSolutionOrder:=',False,'SolutionOrder:=',args.order,'Solver Type:=','Iterative']])
                try:
                    setup.analyze(name)
                except Exception:
                    desktop = base.save_aedt._desktop()._desktop
                    print(desktop.GetMessages(args.project.stem,f'{args.tag}_q3d',0),flush=True)
                    raise
            renderer.analyze_setup = analyze
            run=lom.sim.run
            lom.sim.run=lambda **kw: run(name=args.tag,**kw)
            return lom
    base.q3d.OPEN_PINS = [p for p in base.q3d.OPEN_PINS if p != ('TC','readout_line')]
    if not args.local:
        base.q3d.COMPONENTS = [*base.q3d.COMPONENTS,'res_FX1','res_FX2']
        base.q3d.OPEN_PINS = [p for p in base.q3d.OPEN_PINS if p[1] != 'readout_line']
        base.q3d.OPEN_PINS += [('res_FX1','start'),('res_FX2','start')]
    runner=Runner(setup=dict(max_passes=args.passes,min_passes=args.min_passes,percent_error=args.error))
    try:
        runner.run(args.tag,with_resonators=False)
        (base.RESULTS/f'{args.tag}_geometry.json').write_text(json.dumps(vars(args),indent=2,default=str))
        setup=runner.lom.sim.renderer.pinfo.setup
        conv,_=setup.get_convergence()
        conv.to_csv(base.RESULTS/f'{args.tag}_convergence.csv')
    finally:
        base.save()

if __name__=='__main__':
    main()
