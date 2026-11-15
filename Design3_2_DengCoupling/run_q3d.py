"""Q3D tuning of the lower-pad FX--TC coupling on a copy of Design3_2."""
import argparse
from pathlib import Path
import sys

from shapely.affinity import scale
from shapely.geometry import box
from shapely.ops import unary_union

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT/'Design3_2/simulation'), str(ROOT/'Design3_2/pattern design')]
import run_verification as base
import ftf_layout as layout
from my_tmon import MyTmonCoupler
from my_tmon_fab import MyTmonCouplerFab

def build_variant(gaps, lengths, rounded=False, **options):
    original = layout.MyTmonCoupler
    original_fab = layout.MyTmonCouplerFab
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
                    edge = p.pad_arm_length/2+p.paddle_width
                    # Fingers extend the existing paddle down and wrap under the bottom pad.
                    neck = box(edge-.015, -.014, outer-gap, half_height+gap+.012)
                    tip = box(outer-gap-.008, half_height+gap, outer+length, half_height+gap+.012)
                    part = unary_union([neck, tip])
                    fingers.append(part if side==1 else scale(part,-1,1,origin=(0,0)))
                placed = self._place(fingers)
                table = self.design.qgeometry.tables['poly']
                for name in ('tmon','pad_etch'):
                    idx, = table.index[(table.component == self.id) & (table.name == name)]
                    addition = placed if name == 'tmon' else [g.buffer(.016,join_style=2) for g in placed]
                    table.at[idx,'geometry'] = unary_union([table.at[idx,'geometry'],*addition]).simplify(1e-9)
        return Coupler
    layout.MyTmonCoupler = modified(MyTmonCoupler)
    layout.MyTmonCouplerFab = modified(MyTmonCouplerFab)
    try:
        return (base.rounded_design if rounded else layout.build_design)(**options)
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
    parser.add_argument('--project', type=Path, default=HERE/'aedt/Design3_2_DengCoupling.aedt')
    parser.add_argument('--gap', type=float, nargs=2, default=(4,4))
    parser.add_argument('--length', type=float, nargs=2, default=(45,65))
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--local',action='store_true')
    parser.add_argument('--stem-width', type=float, default=15)
    parser.add_argument('--tc-etch', type=float, default=30)
    parser.add_argument('--passes',type=int,default=12)
    parser.add_argument('--min-passes',type=int,default=4)
    parser.add_argument('--error',type=float,default=.5)
    parser.add_argument('--order', choices=('Low','Medium','High'), default='Medium')
    args=parser.parse_args()
    def build(**kw):
        tmon = dict(kw.pop('tmon_options', None) or {})
        tmon.update(pad_head_width=f'{args.stem_width}um', paddle_gap=f'{args.tc_etch}um')
        return build_variant(args.gap,args.length,rounded=True,tmon_options=tmon,**kw)
    verify_layout(build(for_process=True,fx_options={'array_end_metal':'False'}))
    if args.check:
        return
    base.RESULTS=base.q3d.RESULTS=HERE/'results'
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
    runner=Runner(setup=dict(max_passes=args.passes,min_passes=args.min_passes,percent_error=args.error))
    try:
        runner.run(args.tag,with_resonators=not args.local)
        import json
        (base.RESULTS/f'{args.tag}_geometry.json').write_text(json.dumps(vars(args),indent=2,default=str))
        setup=runner.lom.sim.renderer.pinfo.setup
        conv,_=setup.get_convergence()
        conv.to_csv(base.RESULTS/f'{args.tag}_convergence.csv')
    finally:
        base.save()

if __name__=='__main__':
    main()
