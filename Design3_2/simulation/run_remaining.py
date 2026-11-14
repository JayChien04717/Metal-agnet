"""Design3_2 electrical verification using rounded Qiskit Metal polygons.

Uses the Design3 extraction conventions and the Design3_2 component geometry.
CPW paths and junction sheets retain their native renderer definitions.
"""
from pathlib import Path
import sys
import time

import pandas as pd
from shapely.affinity import scale
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union
from shapely import set_precision

DESIGN=Path(__file__).resolve().parents[1]
ROOT=DESIGN.parent
D3=ROOT/'Design3'
sys.path[:0]=[str(DESIGN/'pattern design'),str(D3/'simulation/q3d'),
              str(D3/'simulation/eigenmode'),str(D3/'simulation/driven'),
              str(ROOT/'Design5/simulation'),str(ROOT/'note')]
from ftf_layout import build_design, RESONATORS, path_length
from round_mask import rounded_ring
import q3d_ftf as q3d
import q3d_kappa as kappa
import hfss_mutual_fx as mutual
from hfss_eigen_res import SETUP, TC_LJ, EigenmodeSim
import save_aedt

RESULTS=DESIGN/'simulation/results'
PROJECT=DESIGN/'simulation/aedt/Design3_2.aedt'
q3d.RESULTS=kappa.RESULTS=mutual.RESULTS=str(RESULTS)

def rounded_design(**options):
    d=build_design(marks=False,include_test_junctions=False,**options)
    polygons=d.qgeometry.tables['poly']
    contacts={}
    for row in d.qgeometry.tables['junction'].itertuples():
        # Preserve the original landing geometry at lumped junction boundaries.
        contacts.setdefault(row.component,[]).append(row.geometry.buffer(row.width/2+.001))
    for component in d.components.values():
        for pin in component.pins.values():
            contacts.setdefault(component.id,[]).append(
                Point(pin.middle).buffer(pin.width/2+pin.gap+.005))
    keep={key:unary_union(items) for key,items in contacts.items()}
    for idx,row in polygons.iterrows():
        # Round in micrometres, the units used by the mask fillet routine.
        g=scale(row.geometry,1000,1000,origin=(0,0))
        window=keep.get(row.component)
        if window is not None:
            window=scale(window,1000,1000,origin=(0,0))
        g=g.simplify(0,preserve_topology=True)
        new=Polygon(rounded_ring(g.exterior,5,window),
                    [rounded_ring(ring,5,window) for ring in g.interiors])
        new=scale(new,.001,.001,origin=(0,0))
        new=set_precision(new, .000005).simplify(.000005, preserve_topology=True)
        if not new.is_valid or new.is_empty or new.geom_type!='Polygon':
            raise ValueError((row.component,row['name'],'invalid rounded polygon'))
        polygons.at[idx,'geometry']=new
    return d

def join_squid_ground(renderer, merge_readouts=False):
    subtract=renderer.subtract_from_ground
    def subtract_and_join():
        subtract()
        if 'squid_ground_TC' in renderer.assign_perfE:
            renderer.modeler.unite(['ground_main_plane','squid_ground_TC'])
            renderer.assign_perfE.remove('squid_ground_TC')
        if merge_readouts:
            for q in RESONATORS:
                pad, trace = f'readout_line_{q}', f'trace_res_{q}'
                if pad in renderer.assign_perfE and trace in renderer.assign_perfE:
                    renderer.modeler.unite([pad, trace])
                    renderer.assign_perfE.remove(trace)
    renderer.subtract_from_ground=subtract_and_join

def save():
    save_aedt._desktop().get_active_project()._project.Save()

class MatrixRunner(q3d.Q3DRunner):
    def _prepare(self,design):
        lom=super()._prepare(design)
        join_squid_ground(lom.sim.renderer, merge_readouts=True)
        original=lom.sim.run
        lom.sim.run=lambda **kw: original(name='D32_R5_matrix_contact',**kw)
        return lom

def run_matrix():
    q3d.build_design=rounded_design
    runner=MatrixRunner(setup=dict(max_passes=20))
    runner.run('d32_matrix',with_resonators=True)
    setup=runner.lom.sim.renderer.pinfo.setup
    conv,_=setup.get_convergence()
    conv.to_csv(RESULTS/'d32_matrix_convergence.csv')
    save()

def run_kappa():
    kappa.build_design=rounded_design
    kappa.KAPPA_SETUP.update(max_passes=24)
    class NamedLOM(kappa.LOManalysis):
        def __init__(self,*args,**kwargs):
            super().__init__(*args,**kwargs)
            original=self.sim.run
            self.sim.run=lambda **kw: original(name='D32_R5_kappa_final',**kw)
    kappa.LOManalysis=NamedLOM
    kappa.run('d32_kappa')
    project=save_aedt._desktop().get_active_project()
    conv,_=project.get_design('D32_R5_kappa_final_q3d').get_setup('Setup').get_convergence()
    conv.to_csv(RESULTS/'d32_kappa_convergence.csv')
    save()

def run_eigen(name):
    d=rounded_design(for_process=True,tmon_options={'Lj':TC_LJ})
    sim=EigenmodeSim(d,'hfss')
    join_squid_ground(sim.renderer)
    sim.setup.update(dict(SETUP, max_passes=24, pct_refinement=23))
    tee=f'tee_{name}'
    sim.run(name=f'D32_R5_eigen_{name}_contact',components=[f'res_{name}',tee,name],
            open_terminations=[(tee,'prime_start'),(tee,'prime_end')])
    conv=sim.convergence_f
    conv=conv.dropna(how='all')
    conv.to_csv(RESULTS/f'd32_eigen_{name}_convergence.csv')
    total=path_length(d,f'res_{name}')+path_length(d,tee,'second_cpw')
    print(name,'path mm',total,'last modes',list(conv.iloc[-1]),flush=True)
    save()

def run_mutual(case):
    tc=case.startswith('tc_')
    q='TC' if tc else case.upper()
    d=rounded_design(for_process=True,fx_options={'flux_as_ground':'False'},
                     tmon_options={'flux_as_ground':'False'})
    jj='rect_jj'
    if tc:
        jj='rect_jj1' if case.endswith('left') else 'rect_jj2'
        comp=d.components[q]
        table=d.qgeometry.tables['junction']
        row,=table.loc[table.component.eq(comp.id)&table.name.ne(jj)].itertuples()
        a,b=[Point(p) for p in row.geometry.coords]
        u,v=(b.x-a.x)/row.geometry.length,(b.y-a.y)/row.geometry.length
        land=.00025
        bridge=LineString([(a.x-land*u,a.y-land*v),(b.x+land*u,b.y+land*v)]).buffer(row.width/2,cap_style=2)
        comp.add_qgeometry('poly',{'mutual_loop_closure':bridge})
        d.qgeometry.tables['junction']=table.drop(index=row.Index)
    mutual.run_design('d32_mutual_fixed_'+case,d,q,'flux_bias_line' if tc else 'flux_line',jj,
                      [(q,'readout_line'),(q,'charge_line')],info={'corner_radius_target_um':5})
    save()

def main():
    cases={'matrix':run_matrix,'kappa':run_kappa}
    args=sys.argv[1:]
    if args==['check']:
        d=rounded_design(for_process=True)
        print('Valid rounded polygons',len(d.qgeometry.tables['poly']),flush=True)
        return
    RESULTS.mkdir(parents=True,exist_ok=True)
    save_aedt.use_project(PROJECT)
    failed=[]
    for case in args:
        try:
            if case in cases:
                cases[case]()
            elif case in RESONATORS:
                run_eigen(case)
            elif case in ('fx1','fx2','tc_left','tc_right'):
                run_mutual(case)
            else:
                raise ValueError(case)
        except Exception:
            import traceback
            traceback.print_exc()
            failed.append(case)
    if failed:
        raise SystemExit('Failed cases: '+', '.join(failed))

if __name__=='__main__':
    main()
