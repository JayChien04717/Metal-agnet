"""Correct recorded port conventions and audit the completed fixed-length search."""
from pathlib import Path
import argparse
import json
import inspect
import shutil
import numpy as np
import pandas as pd
import hammer_fixed_s21 as F
import hammer_search as S
from run_record import create_run, finish_run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('search_run',type=Path)
    args=parser.parse_args(); parent=args.search_run.resolve()
    records=json.loads((parent/'progress.json').read_text())['completed']
    folder=create_run('Design4','resonator','hammer_verified_conventions',entrypoint=__file__,
        inputs=[parent/'run.json',parent/'progress.json'],upstream=[str(parent),*[r['run'] for r in records]],
        parameters=dict(EM_rerun=False,fit_parameters=[],correction='port voltage polarity',
            validation_status='not converged; convention audit only'))
    print('RUN_FOLDER='+str(folder),flush=True)
    import hfss_hammer_straps as E
    from qiskit_metal.renderers.renderer_ansys.hfss_renderer import QHFSSRenderer
    from pyEPR.ansys import Rect
    (folder/'raw/port_renderer_source.txt').write_text(inspect.getsource(QHFSSRenderer.create_ports)
        +'\n'+inspect.getsource(Rect.make_lumped_port)+'\n'+inspect.getsource(Rect.make_center_line),encoding='utf-8')
    evidence=[]
    for row in records:
        run=Path(row['run'])
        F.P.PAPER_HAMMER.update(gap_um=row['gap_um'],height_um=row['height_um'])
        E.SPACING_MM=.4; E.GAP_UM=row['gap_um']; E.LENGTHS=F.P.LENGTHS; E.NO_HAMMER=False
        d,_=E.build(0,autobonds=row['straps']!='cpw')
        current=F.geometry_dump(d); old=json.loads((run/'raw/channel_geometry.json').read_text())
        assert all(current[k]==old[k] for k in ('poly','path')), 'Recorded geometry mismatch'
        polarity=F.port_polarity(d)
        assert polarity['S21_multiplier']==-1
        # Parent data predate the polarity correction. Reject already corrected data.
        comparison=json.loads((run/'comparison.json').read_text())
        assert 'port_polarity' not in comparison
        c=pd.read_csv(run/'raw/hfss_convergence.csv')
        after=pd.read_csv(run/'raw/hfss_convergence_after_fine.csv')
        raw=pd.read_csv(run/'raw/S_fine.csv',index_col=0)
        s11=raw.S11.map(complex).to_numpy(); s21=raw.S21.map(complex).to_numpy()
        evidence.append(dict(tag=row['tag'],geometry_matches=True,polarity=polarity,
            passes=len(c),elements=int(c.iloc[-1]['Solved Elements']),
            last_delta_s=float(c.iloc[-1]['Max Mag. Delta S']),target_delta_s=.01,
            delta_s_target_met=bool(c.iloc[-1]['Max Mag. Delta S']<=.01),
            fine_sweep_adaptive_history_unchanged=c.equals(after),
            maximum_two_port_power_balance_error=float(np.max(abs(abs(s11)**2+abs(s21)**2-1))),
            hammer_fF=comparison['models']['deembedded']['hammer_fF']))
    F.write(folder/'solver_and_polarity_audit.json',evidence)
    result=S.rank(records,folder,model_to_hfss_polarity=-1)
    result['best_dip_frequency_match']=min(result['ranking'],key=lambda x:x.get('modal_frequency_rmse_MHz',np.inf))
    result['all_HFSS_delta_s_targets_met']=all(r['delta_s_target_met'] for r in evidence)
    result['interpretation']='Metric-dependent exploratory candidates, not a validated best fit. No circuit parameter was fitted.'
    result['supersedes']=str(parent/'ranking.json')
    result['reason']='Parent complex-error ranking omitted the deterministic opposite port polarity; magnitude spectra and dip frequencies are unchanged.'
    F.write(folder/'ranking.json',result)
    finish_run(folder,status='completed',summary=result)


if __name__=='__main__': main()
