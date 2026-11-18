
# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
from pathlib import Path
import os
import sys

design_dir = _P('Design2')
root = _P('.')
out = Path(os.environ.get('MUTUAL_OUTPUT_DIR', _S('Design2/results/general/legacy/layout3_full')))
out.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, _I('Design2/layout'))

import numpy as np
import gdstk
from shapely import get_parts
from shapely.geometry import Polygon,Point,box,LineString
from shapely.ops import unary_union
from shapely.affinity import scale
from qiskit_metal import designs,Dict
from qiskit_metal.qlibrary.core import QComponent
from layout3_full import recover,PORTS,SQUID_CENTER_Y

sys.path[:0]=[_I('Design3/simulation/flux'),_I('Design5/simulation')]
import save_aedt
import hfss_mutual_fx as mutual
mutual.RESULTS=_S('Design2/results/general/legacy/layout3_full')
source=_P('Design2/gds/previous/layout3 ftf.gds')
offset_um = float(os.environ.get('TC_FLUX_OFFSET_UM', '0'))
ground,nets=recover(source, tc_flux_offset_um=offset_um)
lib=gdstk.read_gds(source)
metal=unary_union([Polygon(p.points).buffer(0) for p in lib['TOP'].polygons if p.layer==5])

def build(case):
    group='tc' if case.startswith('tc') else case
    line, = [p for p in get_parts(metal) if p.covers(Point(*PORTS['flux_'+group]))]
    if group == 'tc' and offset_um:
        from shapely.affinity import translate
        splice, shoulder = -1760, -1670
        lower = line.intersection(box(-5000, -5000, 5000, splice))
        middle = line.intersection(box(-5000, splice, 5000, shoulder))
        upper = line.intersection(box(-5000, shoulder, 5000, 5000))
        line = unary_union([lower,
            scale(middle, 1, 1-offset_um/(shoulder-splice), origin=(0,splice)),
            translate(upper, yoff=-offset_um)])
    fixed_ground=ground.difference(line)
    polys={'flux':line}
    jj=[]
    if group=='tc':
        polys['island']=nets['transmon']
        x=-15 if case.endswith('left') else 15
        y=SQUID_CENTER_Y
        jj=[('sense',[(x,y+2.5),(x,y-2.5)],7,'50nH')]
        polys['closure']=box(-x-3.5,y-2.75,-x+3.5,y+2.75).difference(fixed_ground)
        bounds=(-.5,-4.95,2.7,-.95)
        chip_center=(1.1,-2.95)
        pin=np.array([[2,-4.549],[2,-4.55]])
    else:
        sign=-1 if group=='left' else 1
        center=sign*485
        arm=sign*539
        polys.update(top=nets[f'fx_{group}_top'].union(box(center-1,-999,center+1,-984.75)),
                     bottom=nets[f'fx_{group}_bottom'].union(box(center-1,-1015.25,center+1,-1001)),
                     array_end=box(arm-1,-1005,arm+1,-995))
        jj=[('sense',[(center,-999),(center,-1001)],2,'34.38nH'),
            ('array_top',[(center+sign,-996),(arm-sign,-996)],2,'100nH'),
            ('array_bottom',[(center+sign,-1004),(arm-sign,-1004)],2,'100nH')]
        # The 2 um-wide lead continues from the main junction to the lower arm.
        bounds=(-4.95,-1.6,-.1,2.6) if sign<0 else (.1,-1.6,4.95,2.6)
        chip_center=((bounds[0]+bounds[2])/2,.5)
        pin=np.array([[sign*4.549,2],[sign*4.55,2]])
    chip=box(*(np.array(bounds)*1000))
    etched=chip.difference(fixed_ground)
    class Fixture(QComponent):
        def make(self):
            for name,p in polys.items():
                for i,s in enumerate(get_parts(p)):
                    self.add_qgeometry('poly',{f'{name}_{i}':scale(s,.001,.001,origin=(0,0))})
            for i,p in enumerate(get_parts(etched)):
                if p.area>1e-6:
                    self.add_qgeometry('poly',{f'cut_{i}':scale(p,.001,.001,origin=(0,0))},subtract=True)
            for name,points,width,inductance in jj:
                self.add_qgeometry('junction',{name:LineString(np.array(points)/1000)},
                    width=width/1000,hfss_inductance=inductance)
            self.add_pin('drive',pin,.5,gap=.1,input_as_norm=True)
    d=designs.DesignPlanar()
    d.chips.main.material='silicon'
    d.chips.main.size.update(center_x=f'{chip_center[0]}mm',center_y=f'{chip_center[1]}mm',
        size_x=f'{bounds[2]-bounds[0]}mm',size_y=f'{bounds[3]-bounds[1]}mm',size_z='-525um')
    Fixture(d,'device')
    return d

if '--check' in sys.argv:
    for case in ('tc_left','tc_right','left','right'):
        design=build(case)
        print(case,design.qgeometry.tables['poly'].shape,design.qgeometry.tables['junction'].shape,flush=True)
    sys.exit()
save_aedt.use_project(Path(os.environ.get('AEDT_PROJECT', _S('Design2/simulation/aedt/D2_Layout3_M.aedt'))))
def run(case):
    d=build(case)
    sim=mutual.ScatteringImpedanceSim(d,'hfss')
    sim.setup.update(mutual.SETUP)
    sim.setup.update(dict(pct_refinement=20))
    if case.startswith('tc'):
        # The shorted JJ and ground are one PEC connection.
        orig=sim.renderer.subtract_from_ground
        def join():
            orig()
            sim.renderer.modeler.unite(['ground_main_plane','closure_0_device','flux_0_device'])
            sim.renderer.assign_perfE.remove('closure_0_device')
            sim.renderer.assign_perfE.remove('flux_0_device')
        sim.renderer.subtract_from_ground=join
    name='D2_M_'+case
    sim.run_sim(name=name,components=['device'],open_terminations=[],
        port_list=[('device','drive',50)],jj_to_port=[('device','sense',50,False)],box_plus_buffer=False)
    f,z,_=sim.renderer.get_params(['Z11','Z21','Z12','Z22'])
    w=2*np.pi*np.asarray(f)
    a,b,c,e=map(np.asarray,z)
    zm=(b+c)/2
    import pandas as pd
    pd.DataFrame(dict(f_MHz=f/1e6,M_pH=zm.imag/w*1e12,
        Z11_L_nH=a.imag/w*1e9,Z22_L_nH=e.imag/w*1e9,
        reciprocity=np.abs(b-c)/np.maximum(np.abs(zm),1e-20))).to_csv(_relocate(out/f'{name}.csv'),index=False)
    conv,_=sim.renderer.pinfo.setup.get_convergence()
    conv.to_csv(_relocate(out/f'{name}_convergence.csv'))
    desk=save_aedt._desktop()
    desk.get_active_project()._project.Save()
    (_relocate(out/f'{name}_messages.log')).write_text('\n'.join(desk._desktop.GetMessages(
        'D2_Layout3_M',name+'_hfss',0)),encoding='utf8')

for case in sys.argv[1:]:
    try:
        run(case)
    except Exception:
        import traceback
        trace = traceback.format_exc()
        print(trace, flush=True)
        (_relocate(out/f'D2_M_{case}_failure.log')).write_text(trace, encoding='utf8')
