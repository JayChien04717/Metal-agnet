"""Full EM candidate: 50 um arms, 30 um head height, reduced gap; fixed main lengths."""
from pathlib import Path
from functools import partial
import argparse
import ast
import inspect
import json
import shutil
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
    ap.add_argument('--width',type=float,default=25.)
    ap.add_argument('--straps',choices=['all','cpw'],default='all')
    ap.add_argument('--spacing-um',type=float,help='Fix main-line centres; derive each arm from gap and head width')
    ap.add_argument('--single-band',type=float,nargs=2,metavar=('LOW_GHZ','HIGH_GHZ'))
    ap.add_argument('--reference-coupon',type=Path,help='Reuse an explicitly named identical straight-CPW reference after geometry/source checks')
    ap.add_argument('--calibration',type=Path,help='Explicit refined Cin/Cp coupon run; baseline geometry reference is unchanged')
    ap.add_argument('--dry-run',action='store_true')
    ap.add_argument('--max-passes',type=int,default=3,help='Explicit HFSS validation budget, default unchanged')
    args=ap.parse_args();baseline=args.baseline.resolve()
    assert args.gap>0
    lead,height,diameter=50.,30.,12.
    if args.spacing_um is not None:
        lead=(args.spacing_um-args.gap-2*args.width)/2
    assert lead>0, 'Head and gap must fit between the fixed main lines'
    assert args.width>4, 'Head width must accommodate the existing corner radius'
    P=H.P;P.PAPER_HAMMER.update(gap_um=args.gap,height_um=height,width_um=args.width)
    original_lengths=json.loads((baseline/'run.json').read_text())['parameters']['lengths_mm']
    assert np.allclose(P.LENGTHS,original_lengths,atol=1e-12,rtol=0), 'Main-line inputs differ from the explicitly selected baseline'
    P.hammer_ground_straps=partial(P.hammer_ground_straps,diameter_um=diameter)
    spacing=2*(lead+P.PAPER_HAMMER['width_um']+args.gap/2)/1000
    E.SPACING_MM=spacing;E.GAP_UM=args.gap;E.LENGTHS=P.LENGTHS;E.NO_HAMMER=False
    inputs=[baseline/'run.json'];upstream=[str(baseline)]
    calibration=args.calibration.resolve() if args.calibration else baseline
    if calibration!=baseline:
        inputs += [calibration/'run.json',calibration/'raw/cp_coupon/caps.json']
        upstream.append(str(calibration))
    if args.reference_coupon:
        args.reference_coupon=args.reference_coupon.resolve()
        inputs += [args.reference_coupon/'geometry.json',args.reference_coupon/'caps.json']
        upstream.append(str(args.reference_coupon.parent.parent))
    folder=create_run('Design4','resonator',f'h30_g{args.gap:g}_w{args.width:g}_a{lead:g}',entrypoint=__file__,
        inputs=inputs,upstream=upstream,parameters=dict(
            arm_um=lead,head_height_um=height,gap_um=args.gap,head_width_um=args.width,spacing_mm=spacing,straps=args.straps,
            fixed_spacing_um=args.spacing_um,arm_source='(main centre spacing - gap - 2*head width)/2' if args.spacing_um is not None else '50 um input',
            single_band_GHz=args.single_band,
            reference_coupon_source=str(args.reference_coupon) if args.reference_coupon else None,
            cin_cp_calibration_source=str(calibration),
            hammer_stub_bridge_diameter_um=diameter,main_lengths_mm=P.LENGTHS,dry_run=args.dry_run,
            target='preserve original effective Purcell-readout coupling; Q3D C alone is not acceptance',
            spacing_exception='User-selected short-arm geometry overrides the older 400 um inter-line spacing preference; meander pitch unchanged',
            ground_bridge_reason='12 um local wire; actual landing and etch checks enforced for every geometry',
            Q3D_max_passes=20,Q3D_requested_percent_error=.1,HFSS_max_passes=args.max_passes,HFSS_mesh_seed_um=10))
    print('RUN_FOLDER='+str(folder),flush=True)
    try:
        d,comps=E.build(0,autobonds=args.straps!='cpw')
        actual=[P.d5.f3.path_length(d,'filter_tail')+P.FILTER_TAP+P.TEE_HALF,
                P.d5.f3.path_length(d,'res_tail')+P.RES_TAP+P.TEE_HALF]
        assert np.allclose(actual,P.LENGTHS,atol=1e-6,rtol=0)
        straps=P.ground_straps(d,args.straps) # checks feet against all etches
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
                shape=r.geometry.buffer(float(r.width)/2,cap_style=2) if kind=='path' else r.geometry
                b=shape.bounds;assert b[0]>=-2.6 and b[2]<=2.6 and b[1]>=-.8 and b[3]<=3.5
                for part in getattr(shape,'geoms',[shape]):
                    if hasattr(part,'exterior'):
                        for a in ax:a.plot(*part.exterior.xy,color='black',lw=.7)
        for s in straps:
            pos=np.array(s['pos']);v=np.array(s['ori'])*s['span']/2
            for a in ax:a.plot([pos[0]-v[0],pos[0]+v[0]],[pos[1]-v[1],pos[1]+v[1]],color='tab:blue',lw=2)
        ax[1].set_xlim(-.05,spacing+.05);ax[1].set_ylim(P.TAP_Y-.1,P.TAP_Y+.1)
        for a in ax:a.set_aspect('equal');a.set_xlabel('x (mm)');a.set_ylabel('y (mm)')
        ax[0].set_title('Fixed main-line lengths');ax[1].set_title(f'Arm{lead:g} / head30 / width{args.width:g} / gap{args.gap:g} um')
        fig.tight_layout();fig.savefig(folder/'layout.png',dpi=170)
        if args.dry_run:
            finish_run(folder,status='completed',summary={'dry_run':True,'geometry_checks':'passed'});return
        import save_aedt,hfss_dirs
        hfss_dirs.use_d_drive();save_aedt.use_project(save_aedt.project_file(H.ROOT/'Design4/simulation'))
        extracted=[]
        for reference,name in [(False,'cj'),(True,'ref')]:
            local,component=H.coupon(lead,reference)
            out=folder/f'raw/{name}_{lead:g}um'
            if reference and args.reference_coupon:
                source=args.reference_coupon
                old_geometry=json.loads((source/'geometry.json').read_text())
                current_geometry=H.geometry_dump(local)
                assert all(old_geometry[k]==current_geometry[k] for k in ('poly','path')), 'Reference geometry differs'
                snapshot=source.parent.parent/'code/Design4/simulation/resonator/hammer_deembed.py'
                old_coupon=next(n for n in ast.parse(snapshot.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='coupon')
                current_coupon=ast.parse(inspect.getsource(H.coupon)).body[0]
                assert ast.dump(old_coupon)==ast.dump(current_coupon), 'Fixed-domain coupon definition changed'
                assert P.d5.SI_MATERIAL=='silicon_eps11p5'
                shutil.copytree(source,out)
                H.write(out/'reused_source.json',dict(source=str(source),geometry_equal=True,fixed_domain_definition_equal=True,
                    material=P.d5.SI_MATERIAL,reason='Fixed main-line spacing leaves identical reference endpoints; no new Q3D reference solve'))
                result=json.loads((source/'caps.json').read_text())
            else:
                result=H.solve_cap(local,component,out,
                    'D4SHORT_'+folder.name[-8:]+'_'+name,max_passes=20,percent_error=.1)
            extracted.append(result)
        cap,ref=extracted;half=ref['ground_fF'][0]/2
        H.plot_sweep([dict(length_um=lead,mutual_fF=cap['mutual_fF'],raw_a_fF=cap['ground_fF'][0],
            raw_b_fF=cap['ground_fF'][1],reference_half_fF=half,
            excess_a_fF=cap['ground_fF'][0]-half,excess_b_fF=cap['ground_fF'][1]-half)],folder)
        hfss_args=[folder,'--calibration',calibration,'--gap',args.gap,
            '--height',height,'--width',args.width,'--lead-um',lead,'--extraction-length-um',lead,
            '--hammer-bridge-diameter-um',diameter,'--straps',args.straps,'--max-passes',args.max_passes]
        if args.single_band:
            hfss_args+=['--single-band',*args.single_band]
        full=S.execute('hammer_fixed_s21.py',hfss_args,folder/'raw/full_hfss.log')
        H.write(folder/'hfss_source.json',{'run':str(full)})
        fit=S.execute('hammer_pole_fit.py',[full],folder/'raw/coupling_fit.log')
        finish_run(folder,status='completed',summary=dict(hfss_run=str(full),coupling_fit_run=str(fit),
            capacitance=cap,reference=ref,target_met='pending comparison to baseline effective coupling'))
    except Exception:
        error=traceback.format_exc();(folder/'error.log').write_text(error)
        finish_run(folder,status='failed',summary={'error':error});raise


if __name__=='__main__':main()
