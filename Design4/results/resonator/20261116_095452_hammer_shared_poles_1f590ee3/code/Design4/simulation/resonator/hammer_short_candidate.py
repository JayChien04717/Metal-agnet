"""Full EM candidate: 50 um arms, 30 um head height, reduced gap; fixed main lengths."""
from pathlib import Path
from functools import partial
import argparse
import traceback
import numpy as np
import hammer_deembed as H
import hammer_search as S
import hfss_hammer_straps as E
from run_record import create_run,finish_run


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('baseline',type=Path)
    ap.add_argument('--gap',type=float,default=20.)
    ap.add_argument('--dry-run',action='store_true')
    args=ap.parse_args();baseline=args.baseline.resolve()
    assert 0<args.gap<42, 'Candidate must reduce the original pad gap'
    lead,height,diameter=50.,30.,12.
    P=H.P;P.PAPER_HAMMER.update(gap_um=args.gap,height_um=height)
    P.hammer_ground_straps=partial(P.hammer_ground_straps,diameter_um=diameter)
    spacing=2*(lead+P.PAPER_HAMMER['width_um']+args.gap/2)/1000
    E.SPACING_MM=spacing;E.GAP_UM=args.gap;E.LENGTHS=P.LENGTHS;E.NO_HAMMER=False
    folder=create_run('Design4','resonator',f'hammer_arm50_height30_gap{args.gap:g}',entrypoint=__file__,
        inputs=[baseline/'run.json'],upstream=[str(baseline)],parameters=dict(
            arm_um=lead,head_height_um=height,gap_um=args.gap,head_width_um=25.,spacing_mm=spacing,
            hammer_stub_bridge_diameter_um=diameter,main_lengths_mm=P.LENGTHS,dry_run=args.dry_run,
            target='preserve original effective Purcell-readout coupling; Q3D C alone is not acceptance',
            ground_bridge_reason='17 um stub landing interval; 12 um wire leaves 2.5 um margin per side',
            Q3D_max_passes=20,Q3D_requested_percent_error=.1,HFSS_max_passes=3,HFSS_mesh_seed_um=10))
    print('RUN_FOLDER='+str(folder),flush=True)
    try:
        d,comps=E.build(0)
        actual=[P.d5.f3.path_length(d,'filter_tail')+P.FILTER_TAP+P.TEE_HALF,
                P.d5.f3.path_length(d,'res_tail')+P.RES_TAP+P.TEE_HALF]
        assert np.allclose(actual,P.LENGTHS,atol=1e-6,rtol=0)
        straps=P.ground_straps(d,'all') # checks feet against all etches
        for j,stub in [('jf','stub_f'),('jr','stub_r')]:
            arm=d.parse_value(d.components[j].options.stub_length)+P.d5.f3.path_length(d,stub)+P.d5.CJ_LEAD/1000
            assert np.isclose(arm,lead/1000,atol=1e-8)
        H.write(folder/'raw/full_geometry.json',H.geometry_dump(d))
        H.write(folder/'geometry_audit.json',dict(actual_main_lengths_mm=actual,actual_arm_um=lead,
            spacing_mm=spacing,ground_straps=straps,auto_bond_count=len(P.wire_bond_positions(d)),
            ground_foot_etch_checks_passed=True,simulation_not_fabrication_signoff=True))
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(1,2,figsize=(12,5))
        ids={d.components[n].id for n in comps}
        for kind in ('poly','path'):
            for r in d.qgeometry.tables[kind].itertuples():
                if r.component not in ids:continue
                assert r.geometry.is_valid
                shape=r.geometry.buffer(float(r.width)/2) if kind=='path' else r.geometry
                b=shape.bounds;assert b[0]>=-2.6 and b[2]<=2.6 and b[1]>=-.8 and b[3]<=3.5
                for part in getattr(shape,'geoms',[shape]):
                    if hasattr(part,'exterior'):
                        for a in ax:a.plot(*part.exterior.xy,color='black',lw=.7)
        for s in straps:
            pos=np.array(s['pos']);v=np.array(s['ori'])*s['span']/2
            for a in ax:a.plot([pos[0]-v[0],pos[0]+v[0]],[pos[1]-v[1],pos[1]+v[1]],color='tab:blue',lw=2)
        ax[1].set_xlim(-.05,spacing+.05);ax[1].set_ylim(P.TAP_Y-.1,P.TAP_Y+.1)
        for a in ax:a.set_aspect('equal');a.set_xlabel('x (mm)');a.set_ylabel('y (mm)')
        ax[0].set_title('Fixed main-line lengths');ax[1].set_title(f'Arm50 / head30 / gap{args.gap:g} um')
        fig.tight_layout();fig.savefig(folder/'layout.png',dpi=170)
        if args.dry_run:
            finish_run(folder,status='completed',summary={'dry_run':True,'geometry_checks':'passed'});return
        import save_aedt,hfss_dirs
        hfss_dirs.use_d_drive();save_aedt.use_project(save_aedt.project_file(H.ROOT/'Design4/simulation'))
        extracted=[]
        for reference,name in [(False,'cj'),(True,'ref')]:
            local,component=H.coupon(lead,reference)
            result=H.solve_cap(local,component,folder/f'raw/{name}_50um',
                'D4SHORT_'+folder.name[-8:]+'_'+name,max_passes=20,percent_error=.1)
            extracted.append(result)
        cap,ref=extracted;half=ref['ground_fF'][0]/2
        H.plot_sweep([dict(length_um=lead,mutual_fF=cap['mutual_fF'],raw_a_fF=cap['ground_fF'][0],
            raw_b_fF=cap['ground_fF'][1],reference_half_fF=half,
            excess_a_fF=cap['ground_fF'][0]-half,excess_b_fF=cap['ground_fF'][1]-half)],folder)
        full=S.execute('hammer_fixed_s21.py',[folder,'--calibration',baseline,'--gap',args.gap,
            '--height',height,'--lead-um',lead,'--extraction-length-um',lead,
            '--hammer-bridge-diameter-um',diameter,'--straps','all'],folder/'raw/full_hfss.log')
        H.write(folder/'hfss_source.json',{'run':str(full)})
        fit=S.execute('hammer_coupling_fit.py',[full],folder/'raw/coupling_fit.log')
        finish_run(folder,status='completed',summary=dict(hfss_run=str(full),coupling_fit_run=str(fit),
            capacitance=cap,reference=ref,target_met='pending comparison to baseline effective coupling'))
    except Exception:
        error=traceback.format_exc();(folder/'error.log').write_text(error)
        finish_run(folder,status='failed',summary={'error':error});raise


if __name__=='__main__':main()
