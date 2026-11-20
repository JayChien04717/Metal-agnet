"""Extract one candidate while reusing the identical outer-endpoint CPW reference."""
from pathlib import Path
import argparse
import json
import shutil
import traceback
import hammer_deembed as H
from run_record import create_run,finish_run


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('reference_run',type=Path)
    ap.add_argument('--gap',type=float,required=True)
    ap.add_argument('--height',type=float,default=65)
    args=ap.parse_args(); reference=args.reference_run.resolve()
    original=dict(H.P.PAPER_HAMMER)
    endpoint=original['gap_um']/2+original['width_um']+800
    H.P.PAPER_HAMMER.update(gap_um=args.gap,height_um=args.height)
    lead=endpoint-args.gap/2-original['width_um']
    folder=create_run('Design4','resonator',f'hammer_gap{args.gap:g}_height{args.height:g}',
        entrypoint=__file__,inputs=[reference/'run.json',reference/'raw/ref_800um/caps.json'],
        upstream=[str(reference)],parameters=dict(hammer=H.P.PAPER_HAMMER,lead_um=lead,
            outer_endpoint_um=endpoint,max_passes=20,percent_error=.1,channel_lengths_mm=H.P.LENGTHS))
    print('RUN_FOLDER='+str(folder),flush=True)
    try:
        import save_aedt,hfss_dirs
        href,_=H.coupon(lead,True)
        current=H.geometry_dump(href)
        previous=json.loads((reference/'raw/ref_800um/geometry.json').read_text())
        # Numeric geometry is identical even when option strings differ in formatting.
        for kind in ('path','poly'):
            assert current[kind]==previous[kind], 'Reference geometry mismatch'
        out=folder/'raw'/f'ref_{lead:g}um'
        shutil.copytree(reference/'raw/ref_800um',out)
        H.write(out/'reused_source.json',dict(source=str(reference/'raw/ref_800um'),geometry_verified=True))
        hfss_dirs.use_d_drive(); save_aedt.use_project(save_aedt.project_file(H.ROOT/'Design4/simulation'))
        d,c=H.coupon(lead)
        result=H.solve_cap(d,c,folder/'raw'/f'cj_{lead:g}um','D4CAND_'+folder.name[-8:],max_passes=20,percent_error=.1)
        ref=json.loads((out/'caps.json').read_text()); half=ref['ground_fF'][0]/2
        rows=[dict(length_um=lead,mutual_fF=result['mutual_fF'],raw_a_fF=result['ground_fF'][0],
            raw_b_fF=result['ground_fF'][1],reference_half_fF=half,
            excess_a_fF=result['ground_fF'][0]-half,excess_b_fF=result['ground_fF'][1]-half)]
        H.plot_sweep(rows,folder)
        finish_run(folder,status='completed',summary=dict(rows=rows,reference_reused=True))
    except Exception:
        (folder/'error.log').write_text(traceback.format_exc(),encoding='utf-8')
        finish_run(folder,status='failed',summary={'error':traceback.format_exc()}); raise


if __name__=='__main__':main()
