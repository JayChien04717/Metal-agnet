"""Actual 400/800 um hammer arms in full-channel HFSS; main lengths fixed."""
from pathlib import Path
import argparse
import json
import traceback
import numpy as np
import pandas as pd
import hammer_fixed_s21 as F
import hammer_search as S
import hfss_hammer_straps as E
from run_record import create_run,finish_run


def layout_check(lead):
    P=F.P
    P.PAPER_HAMMER.update(gap_um=42.,height_um=65.)
    spacing=2*(lead+P.PAPER_HAMMER['width_um']+P.PAPER_HAMMER['gap_um']/2)/1000
    E.SPACING_MM=spacing; E.GAP_UM=42.; E.LENGTHS=P.LENGTHS; E.NO_HAMMER=False
    d,comps=E.build(0)
    parse=d.parse_value
    chip=d.chips.main.size
    bounds=[parse(chip.center_x)-parse(chip.size_x)/2,parse(chip.center_y)-parse(chip.size_y)/2,
            parse(chip.center_x)+parse(chip.size_x)/2,parse(chip.center_y)+parse(chip.size_y)/2]
    checked=[]
    for kind in ('poly','path'):
        for row in d.qgeometry.tables[kind].itertuples():
            name=next(c.name for c in d.components.values() if c.id==row.component)
            if name not in comps: continue
            assert row.geometry.is_valid,(name,row.name)
            shape=row.geometry.buffer(float(row.width)/2) if kind=='path' else row.geometry
            b=shape.bounds
            assert b[0]>=bounds[0] and b[1]>=bounds[1] and b[2]<=bounds[2] and b[3]<=bounds[3],(name,b,bounds)
            checked.append(shape)
    actual_lengths=[P.d5.f3.path_length(d,'filter_tail')+P.FILTER_TAP+P.TEE_HALF,
                    P.d5.f3.path_length(d,'res_tail')+P.RES_TAP+P.TEE_HALF]
    assert np.allclose(actual_lengths,P.LENGTHS,rtol=0,atol=1e-6)
    # Tap branch + straight connector + the hammer's own 10 um lead.
    actual_arm=parse(d.components['jf'].options.stub_length)+P.d5.f3.path_length(d,'stub_f')+P.d5.CJ_LEAD/1000
    assert np.isclose(actual_arm,lead/1000,atol=1e-8)
    right_arm=parse(d.components['jr'].options.stub_length)+P.d5.f3.path_length(d,'stub_r')+P.d5.CJ_LEAD/1000
    assert np.isclose(right_arm,lead/1000,atol=1e-8)
    info=dict(physical_lead_um=lead,definition='head outer edge to main-line tap centre, per side',
        spacing_mm=spacing,tee_to_hammer_centre_mm=spacing/2,
        actual_arms_mm=[actual_arm,right_arm],actual_main_lengths_mm=actual_lengths,
        chip_bounds_mm=bounds,all_rendered_planar_geometry_within_chip=True,
        explicit_straps=P.ground_straps(d,'all'),auto_bond_count=len(P.wire_bond_positions(d)))
    return d,checked,info


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('reference_run',type=Path)
    ap.add_argument('baseline_run',type=Path)
    ap.add_argument('--dry-run',action='store_true')
    args=ap.parse_args(); reference=args.reference_run.resolve(); baseline=args.baseline_run.resolve()
    folder=create_run('Design4','resonator','hammer_physical_arms',entrypoint=__file__,
        inputs=[reference/'run.json',baseline/'run.json'],upstream=[str(reference),str(baseline)],
        parameters=dict(physical_leads_um=[400,800],dry_run=args.dry_run,main_lengths_mm=F.P.LENGTHS,
            hammer_gap_um=42,hammer_height_um=65,straps='all',solver_defaults_preserved=True))
    print('RUN_FOLDER='+str(folder),flush=True)
    records=[]
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axs=plt.subplots(1,2,figsize=(12,5))
        audits=[]
        for lead,ax in zip((400,800),axs):
            d,shapes,info=layout_check(lead);audits.append(info)
            F.write(folder/f'raw/geometry_{lead}um.json',F.geometry_dump(d))
            for shape in shapes:
                for part in getattr(shape,'geoms',[shape]):
                    if hasattr(part,'exterior'):
                        x,y=part.exterior.xy; ax.plot(x,y,lw=.45,color='black')
            ax.set_aspect('equal');ax.set_title(f'{lead} um each arm; main lengths fixed')
            ax.set_xlabel('x (mm)');ax.set_ylabel('y (mm)')
        fig.tight_layout();fig.savefig(folder/'physical_layouts.png',dpi=160)
        F.write(folder/'layout_audit.json',audits)
        if args.dry_run:
            finish_run(folder,status='completed',summary={'dry_run':True,'audits':audits});return
        for lead in (400,800):
            run=S.execute('hammer_fixed_s21.py',[reference,'--calibration',baseline,
                '--lead-um',lead,'--extraction-length-um',lead,'--gap',42,'--height',65,'--straps','all'],
                folder/f'raw/arm_{lead}um_hfss.log')
            records.append(dict(tag=f'physical_arm_{lead}um',gap_um=42.,height_um=65.,straps='all',
                physical_lead_um=lead,run=str(run)))
            F.write(folder/'progress.json',{'completed':records})
        # Arm length changes move the whole spectrum. Fine-band intersection may
        # contain neither resonance, so never rank only that intersection.
        finish_run(folder,status='completed',summary=dict(cases=records,layout_audit=audits,
            analysis_entrypoint='hammer_arms_report.py',
            scope='Two actual full-channel arm lengths, not a hammer-gap search'))
    except Exception:
        error=traceback.format_exc();(folder/'error.log').write_text(error,encoding='utf-8')
        finish_run(folder,status='failed',summary={'error':error,'completed':records});raise


if __name__=='__main__':main()
