"""Refine the explicitly selected unchanged Cin/Cp Q3D coupon on its native setup."""
from pathlib import Path
import argparse
import json
import shutil
import traceback
import numpy as np
import hammer_deembed as H
from run_record import create_run,finish_run


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('source',type=Path)
    args=ap.parse_args();source=args.source.resolve();coupon=source/'raw/cp_coupon'
    old=json.loads((coupon/'caps.json').read_text())
    folder=create_run('Design4','resonator','hammer_cp_refine',entrypoint=__file__,
        inputs=[source/'run.json',coupon/'caps.json',coupon/'geometry.json'],upstream=[str(source)],
        parameters=dict(geometry_changed=False,max_passes=20,percent_error=.1,
            reason='Previously reused Cin/Cp coupon reports nonconvergence',mesh_policy='Continue existing native Q3D setup'))
    print('RUN_FOLDER='+str(folder),flush=True)
    try:
        import save_aedt
        save_aedt.use_project(save_aedt.project_file(H.ROOT/'Design4/simulation'))
        project=save_aedt._desktop().get_active_project()
        design=next(d for d in project.get_designs() if d.name==old['aedt_design'])
        setup=design.get_setup('Setup');out=folder/'raw/cp_coupon';out.mkdir()
        shutil.copy2(coupon/'geometry.json',out/'geometry.json')
        before=setup.get_convergence()[0];before.to_csv(out/'convergence_before.csv')
        setup.max_pass=20
        design._design.ChangeProperty(['NAME:AllTabs',['NAME:CG',['NAME:PropServers','AnalysisSetup:Setup'],
            ['NAME:ChangedProps',['NAME:Percent Error','Value:=','0.1']]]])
        setup.analyze();after,native=setup.get_convergence()
        after.to_csv(out/'convergence_after.csv');(out/'convergence.txt').write_text(native)
        cm,units,_,_=setup.get_matrix();assert units=='fF';cm.to_csv(out/'matrix.csv')
        sig=old['signals'];matrix=cm.loc[sig,sig].to_numpy(float);passes=[]
        for n in range(1,len(after)+1):
            m,u,_,_=setup.get_matrix(solution_kind='AdaptivePass',pass_number=n);assert u=='fF'
            (m*1e-15).to_csv(out/f'pass_{n}.csv');passes.append(m.loc[sig,sig].to_numpy(float))
        delta=float(np.max(abs(passes[-1]-passes[-2])/np.maximum(abs(passes[-1]),1e-30))) if len(passes)>1 else None
        def field(label):
            lines=[line for line in native.splitlines() if label in line]
            try:
                return float(lines[0].partition(':')[2]) if len(lines)==1 else None
            except ValueError:
                return None
        target,current,minimum=field('Target'),field('Current'),field('Minimum')
        converged=(current<=target and minimum>1) if None not in (target,current,minimum) else None
        result=dict(units=units,signals=sig,mutual_fF=float(-matrix[0,1]),ground_fF=matrix.sum(axis=1).tolist(),
            aedt_design=design.name,max_passes=20,percent_error=.1,actual_passes=len(after),
            adaptive_history_preserved=bool(after.iloc[:len(before)].equals(before)),
            last_pass_max_relative_change=delta,native_convergence_last_row=after.iloc[-1].to_dict(),
            solver_converged=converged,native_target_percent=target,native_current_percent=current,
            convergence_source='convergence.txt; native Current versus Target with minimum-pass check',
            ground_bonds='none: unchanged local coupon, common electrostatic ground')
        H.write(out/'caps.json',result);project._project.Save()
        finish_run(folder,status='completed',summary=result);print(json.dumps(result),flush=True)
    except Exception:
        error=traceback.format_exc();(folder/'error.log').write_text(error)
        finish_run(folder,status='failed',summary={'error':error});raise


if __name__=='__main__':main()
