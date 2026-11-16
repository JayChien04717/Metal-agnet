"""Bounded physical-geometry comparison; fixed filter/resonator lengths throughout.

Compare the baseline, denser CPW straps, gap 20/height 65, and gap 42/height 100.
Rank by complex S21 error on the SAME frequency interval, without fitted lengths,
capacitances, phases or amplitudes. This identifies the best TESTED candidate,
not a global optimum or a validated fabrication release.
"""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import traceback
import numpy as np
import pandas as pd
import hammer_deembed as H
from run_record import create_run,finish_run


def execute(script, arguments, logfile):
    command=[sys.executable,str(H.ROOT/'Design4/simulation/resonator'/script),*map(str,arguments)]
    found=None
    with logfile.open('x',encoding='utf-8') as stream:
        process=subprocess.Popen(command,cwd=H.ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
            text=True,encoding='utf-8',errors='replace',env={**H.os.environ,'PYTHONIOENCODING':'utf-8'})
        for line in process.stdout:
            stream.write(line); stream.flush()
            if line.startswith('RUN_FOLDER='):
                found=Path(line.strip().split('=',1)[1]); print(line.strip(),flush=True)
            elif 'Analyzing setup' in line or line.startswith('REUSE'):
                print(line.strip(),flush=True)
        status=process.wait()
    if status or found is None:
        raise RuntimeError(f'{script} failed ({status}); see {logfile}')
    return found


def rank(records,folder):
    from scipy.signal import find_peaks
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    frames=[]
    for row in records:
        frame=pd.read_csv(Path(row['run'])/'model_comparison.csv')
        for name in ('HFSS_S21','deembedded_S21'): frame[name]=frame[name].map(complex)
        frames.append(frame)
    lo=max(float(f.frequency_GHz.min()) for f in frames)
    hi=min(float(f.frequency_GHz.max()) for f in frames)
    if hi<=lo: raise RuntimeError('No common frequency interval')
    grid=np.arange(lo,hi,.000025)
    scores=[]
    fig,axes=plt.subplots(len(records),1,figsize=(10,3*len(records)),sharex=True,squeeze=False)
    exported=pd.DataFrame({'frequency_GHz':grid})
    for row,frame,ax in zip(records,frames,axes[:,0]):
        actual=np.interp(grid,frame.frequency_GHz,frame.HFSS_S21)
        predicted=np.interp(grid,frame.frequency_GHz,frame.deembedded_S21)
        error=abs(predicted-actual)
        hfss_db=20*np.log10(np.maximum(abs(actual),1e-12))
        model_db=20*np.log10(np.maximum(abs(predicted),1e-12))
        ih,ph=find_peaks(-hfss_db,prominence=.1)
        im,pm=find_peaks(-model_db,prominence=.1)
        modal={}
        if len(ih)>=2 and len(im)>=2:
            fh=np.sort(grid[ih[np.argsort(ph['prominences'])[-2:]]])
            fm=np.sort(grid[im[np.argsort(pm['prominences'])[-2:]]])
            modal=dict(HFSS_two_dips_GHz=fh.tolist(),model_two_dips_GHz=fm.tolist(),
                modal_frequency_rmse_MHz=float(np.sqrt(np.mean((fh-fm)**2))*1000),
                splitting_error_MHz=float(abs(np.diff(fh)[0]-np.diff(fm)[0])*1000))
        scores.append(dict(**row,rms_complex_S21=float(np.sqrt(np.mean(error**2))),
            maximum_complex_error=float(error.max()),HFSS_dips_GHz=grid[ih].tolist(),model_dips_GHz=grid[im].tolist(),
            dip_prominence_threshold_dB=.1,**modal))
        exported[row['tag']+'_HFSS']=actual; exported[row['tag']+'_model']=predicted
        ax.plot(grid,hfss_db,'k',label='HFSS'); ax.plot(grid,model_db,label='Q3D circuit, fixed lengths')
        ax.set_title(f"{row['tag']}: gap {row['gap_um']:g} um, height {row['height_um']:g} um, straps {row['straps']}")
        ax.set_ylabel('S21 (dB)'); ax.grid(alpha=.25); ax.legend()
    axes[-1,0].set_xlabel('Frequency (GHz)'); fig.tight_layout(); fig.savefig(folder/'candidate_spectra.png',dpi=160)
    exported.to_csv(folder/'common_grid.csv',index=False)
    scores.sort(key=lambda r:r['rms_complex_S21'])
    result=dict(search_domain='four explicit candidates including baseline; not a global optimum',
        fixed_line_lengths_mm=H.P.LENGTHS,frequency_interval_GHz=[lo,hi],grid_step_kHz=25,
        objective='RMS complex S21 difference, no fitted phase/amplitude/length/capacitance',
        best_tested=scores[0],ranking=scores,
        limitations=['Q3D shunt subtraction still has mesh/extension sensitivity',
            'cin/cp endpoint contributions and bridge/distributed coupling are not fully represented',
            'HFSS three-pass setup is not independently converged by pass count alone',
            'No claim of global optimum, target readout performance, or fabrication release'])
    H.write(folder/'ranking.json',result)
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('baseline',type=Path)
    ap.add_argument('reference_run',type=Path)
    args=ap.parse_args(); baseline=args.baseline.resolve(); reference=args.reference_run.resolve()
    candidates=[dict(tag='dense_ground',gap_um=42.,height_um=65.,straps='cpw'),
        dict(tag='gap20_height65',gap_um=20.,height_um=65.,straps='all'),
        dict(tag='gap42_height100',gap_um=42.,height_um=100.,straps='all')]
    folder=create_run('Design4','resonator','hammer_geometry_search',entrypoint=__file__,
        inputs=[baseline/'run.json',reference/'run.json'],upstream=[str(baseline),str(reference)],
        parameters=dict(candidates=candidates,fixed_lengths_mm=H.P.LENGTHS,
            objective='complex S21 prediction error at fixed physical lengths',
            reference_reuse='identical end coordinates, verified geometry',release_selection=False))
    print('SEARCH_FOLDER='+str(folder),flush=True)
    records=[dict(tag='baseline_gap42_height65',gap_um=42.,height_um=65.,straps='all',run=str(baseline))]
    failures=[]
    try:
        for candidate in candidates:
            print('CANDIDATE',candidate,flush=True)
            try:
                if candidate['tag']=='dense_ground':
                    caprun=reference  # capacitor geometry unchanged; bridge loading is an explicit model limitation.
                else:
                    caprun=execute('hammer_candidate.py',[reference,'--gap',candidate['gap_um'],'--height',candidate['height_um']],
                        folder/'raw'/f"{candidate['tag']}_q3d.log")
                fullrun=execute('hammer_fixed_s21.py',[caprun,'--gap',candidate['gap_um'],'--height',candidate['height_um'],
                    '--straps',candidate['straps'],'--calibration',baseline],folder/'raw'/f"{candidate['tag']}_hfss.log")
                records.append(dict(**candidate,run=str(fullrun),capacitance_run=str(caprun)))
                H.write(folder/'progress.json',dict(completed=records,failures=failures))
            except Exception as exc:
                failures.append(dict(candidate=candidate,error=repr(exc)))
                H.write(folder/'progress.json',dict(completed=records,failures=failures))
        result=rank(records,folder); result['failures']=failures
        finish_run(folder,status='completed' if not failures else 'failed',summary=result)
    except Exception:
        (folder/'error.log').write_text(traceback.format_exc(),encoding='utf-8')
        finish_run(folder,status='failed',summary={'error':traceback.format_exc()}); raise


if __name__=='__main__':main()
