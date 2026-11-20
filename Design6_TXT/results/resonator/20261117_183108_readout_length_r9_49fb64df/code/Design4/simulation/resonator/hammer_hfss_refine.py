"""Explicit convergence validation on an existing HFSS setup; no layout rebuild."""
from pathlib import Path
import argparse
import json
import shutil
import traceback
import numpy as np
import hammer_deembed as H
from run_record import create_run,finish_run


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source',type=Path);ap.add_argument('--max-passes',type=int,default=8)
    args=ap.parse_args();source=args.source.resolve()
    original=json.loads((source/'run.json').read_text())
    hfss=json.loads((source/'raw/hfss.json').read_text())
    folder=create_run('Design4','resonator','hammer_hfss_refine',entrypoint=__file__,
        inputs=[source/'run.json',source/'raw/hfss.json',source/'raw/hfss_convergence.csv'],upstream=[str(source)],
        parameters=dict(original_geometry_parameters=original['parameters'],max_passes=args.max_passes,
            geometry_changed=False,reason='explicit validation of unconverged coupling comparison; default runner remains three passes',
            adaptive_frequency_changed=False,delta_s_target_changed=False,mesh_policy='reuse native setup; verify whether adaptive history is retained'))
    print('RUN_FOLDER='+str(folder),flush=True)
    try:
        import save_aedt
        save_aedt.use_project(save_aedt.project_file(H.ROOT/'Design4/simulation'))
        project=save_aedt._desktop().get_active_project()
        design=next(x for x in project.get_designs() if x.name==hfss['design'])
        setup=design.get_setup('Setup')
        before=setup.get_convergence()[0];before.to_csv(folder/'raw/convergence_before.csv')
        properties={k:design._design.GetPropertyValue('HfssTab','AnalysisSetup:Setup',k)
            for k in design._design.GetProperties('HfssTab','AnalysisSetup:Setup')}
        H.write(folder/'raw/setup_before.json',properties)
        setup.passes=args.max_passes
        print('REFINE',design.name,'max passes',args.max_passes,flush=True)
        setup.analyze('Setup : Fine')
        after=setup.get_convergence()[0];after.to_csv(folder/'raw/hfss_convergence.csv')
        sweep=setup.get_sweep('Fine');freq,values=sweep.get_network_data(['S11','S21'])
        H.pd.DataFrame(np.asarray(values).T,index=np.asarray(freq)/1e9,columns=['S11','S21']).to_csv(folder/'raw/S_fine.csv')
        final=setup.get_convergence()[0];final.to_csv(folder/'raw/hfss_convergence_after_fine.csv')
        preserved=after.iloc[:len(before)].equals(before)
        result=dict(design=design.name,source=str(source),geometry_changed=False,
            adaptive_history_preserved=preserved,passes=len(after),
            last_delta_s=float(after.iloc[-1]['Max Mag. Delta S']),target_delta_s=.01,
            convergence_target_met=bool(after.iloc[-1]['Max Mag. Delta S']<=.01),
            requested_max_passes=args.max_passes,fine_GHz=[float(np.min(freq))/1e9,float(np.max(freq))/1e9],
            numerical_validation='Target check is not field identification or direct-frequency validation')
        if (source/'raw/channel_geometry.json').exists():
            shutil.copy2(source/'raw/channel_geometry.json',folder/'raw/channel_geometry.json')
        project._project.Save()
        H.write(folder/'raw/hfss.json',result)
        finish_run(folder,status='completed',summary=result)
        print(json.dumps(result),flush=True)
    except Exception:
        error=traceback.format_exc();(folder/'error.log').write_text(error)
        finish_run(folder,status='failed',summary={'error':error});raise


if __name__=='__main__':main()
