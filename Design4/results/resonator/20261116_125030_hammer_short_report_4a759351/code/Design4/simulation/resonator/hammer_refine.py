"""Continue only the two longest hammer/reference pairs on their existing meshes.
Original exported results remain immutable. New run includes explicit parent data.
"""
from pathlib import Path
import argparse
import json
import shutil
import traceback
import numpy as np
import pandas as pd
import hammer_deembed as H
from run_record import create_run,finish_run


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('upstream',type=Path)
    args=ap.parse_args(); parent=args.upstream.resolve()
    folder=create_run('Design4','resonator','hammer_cpw_refine',entrypoint=__file__,
        inputs=[parent/'run.json',parent/'cpw_sweep.csv'],upstream=[str(parent)],
        parameters=dict(lengths_um=[400,800],max_passes=20,percent_error=.1,
            reason='8-pass reference CPWs did not converge; subtraction amplifies matrix errors',
            geometry_change=False,mesh_policy='continue existing setup; do not rerender'))
    print('RUN_FOLDER='+str(folder),flush=True)
    try:
        import save_aedt
        save_aedt.use_project(save_aedt.project_file(H.ROOT/'Design4/simulation'))
        project=save_aedt._desktop().get_active_project()
        rows=pd.read_csv(parent/'cpw_sweep.csv')
        for source in (parent/'raw').iterdir():
            if source.is_dir(): shutil.copytree(source,folder/'raw'/source.name)
        H.write(folder/'copied_parent_data.json',dict(source=str(parent),
            unchanged_cases=['cj_50um','ref_50um','cj_100um','ref_100um','cj_200um','ref_200um'],
            refined_cases=['cj_400um','ref_400um','cj_800um','ref_800um']))
        for length in (400,800):
            solved=[]
            for comp in ('cj','ref'):
                out=folder/'raw'/f'{comp}_{length}um'
                old=json.loads((out/'caps.json').read_text())
                design=next(d for d in project.get_designs() if d.name==old['aedt_design'])
                setup=design.get_setup('Setup')
                before=setup.get_convergence()[0]; before.to_csv(out/'convergence_before.csv')
                setup.max_pass=20
                # Setter uses the CG property server; a string preserves fractional precision.
                design._design.ChangeProperty(['NAME:AllTabs',['NAME:CG',
                    ['NAME:PropServers','AnalysisSetup:Setup'],
                    ['NAME:ChangedProps',['NAME:Percent Error','Value:=','0.1']]]])
                print('CONTINUE',design.name,flush=True)
                setup.analyze()
                after=setup.get_convergence()[0]; after.to_csv(out/'convergence_after.csv')
                cm,units,_,_=setup.get_matrix(); assert units=='fF',units
                cm.to_csv(out/'matrix.csv')
                sig=old['signals']; matrix=cm.loc[sig,sig].to_numpy(float)
                for n in range(1,len(after)+1):
                    m,u,_,_=setup.get_matrix(solution_kind='AdaptivePass',pass_number=n)
                    assert u=='fF',u
                    (m*1e-15).to_csv(out/f'pass_{n}.csv')
                result=dict(units=units,signals=sig,mutual_fF=float(-matrix[0,1]) if comp=='cj' else 0.,
                    ground_fF=matrix.sum(axis=1).tolist(),aedt_design=design.name,
                    parent_run=str(parent),max_passes=20,percent_error=.1,
                    adaptive_history_preserved=bool(after.iloc[:len(before)].equals(before)),
                    actual_passes=len(after))
                H.write(out/'caps.json',result); solved.append(result)
                print(json.dumps(result),flush=True)
            c,r=solved; baseline=r['ground_fF'][0]/2
            idx=rows.index[rows.length_um==length][0]
            rows.loc[idx,['mutual_fF','raw_a_fF','raw_b_fF','reference_half_fF','excess_a_fF','excess_b_fF']]=[
                c['mutual_fF'],*c['ground_fF'],baseline,c['ground_fF'][0]-baseline,c['ground_fF'][1]-baseline]
            H.plot_sweep(rows.to_dict('records'),folder)
        project._project.Save()
        finish_run(folder,status='completed',summary=dict(rows=rows.to_dict('records'),
            remaining='fixed-length HFSS comparison',validation='evaluate last-pass and length sensitivity'))
    except Exception:
        (folder/'error.log').write_text(traceback.format_exc(),encoding='utf-8')
        finish_run(folder,status='failed',summary={'error':traceback.format_exc()})
        raise


if __name__=='__main__':main()
