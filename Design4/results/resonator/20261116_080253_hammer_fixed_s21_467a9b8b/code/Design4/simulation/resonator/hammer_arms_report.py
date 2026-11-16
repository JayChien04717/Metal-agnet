"""Report actual long-arm spectra in their own bands; no overlap-only ranking."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import hammer_fixed_s21 as F
from run_record import create_run,finish_run


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('parent',type=Path)
    args=ap.parse_args();parent=args.parent.resolve()
    records=json.loads((parent/'progress.json').read_text())['completed']
    assert len(records)==2, 'Both physical-arm solves must be complete'
    folder=create_run('Design4','resonator','hammer_physical_arms_report',entrypoint=__file__,
        inputs=[parent/'run.json',parent/'progress.json',parent/'layout_audit.json'],
        upstream=[str(parent),*[r['run'] for r in records]],
        parameters=dict(EM_rerun=False,fit_parameters=[],score_windows='own fine bands, not a common-overlap ranking'))
    print('RUN_FOLDER='+str(folder),flush=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,1,figsize=(11,12))
    results=[]
    for row,ax in zip(records,axes[:2]):
        run=Path(row['run'])
        result=json.loads((run/'comparison.json').read_text())
        assert result['port_polarity']['S21_multiplier']==-1
        curve=pd.read_csv(run/'model_comparison.csv')
        freq=curve.frequency_GHz.to_numpy()
        models={name:curve[name].map(complex).to_numpy() for name in ('HFSS_S21','raw_full_arm_S21','deembedded_S21')}
        metrics={}
        for name,z in models.items():
            db=20*np.log10(np.maximum(abs(z),1e-12))
            peaks,props=find_peaks(-db,prominence=1)
            selected=peaks[np.argsort(props['prominences'])[-2:]] if len(peaks)>=2 else peaks
            dips=np.sort(freq[selected])
            ax.plot(freq,db,label=name.replace('_S21',''),lw=1.2)
            metrics[name]=dict(dips_GHz=dips.tolist(),all_dips_GHz=freq[peaks].tolist(),
                rms_complex_S21=float(np.sqrt(np.mean(abs(z-models['HFSS_S21'])**2))))
        actual=np.array(metrics['HFSS_S21']['dips_GHz'])
        for name in ('raw_full_arm_S21','deembedded_S21'):
            predicted=np.array(metrics[name]['dips_GHz'])
            if len(actual)==len(predicted)==2:
                metrics[name].update(dip_frequency_rmse_MHz=float(np.sqrt(np.mean((actual-predicted)**2))*1000),
                    splitting_error_MHz=float(abs(np.diff(actual)[0]-np.diff(predicted)[0]))*1000)
        conv=pd.read_csv(run/'raw/hfss_convergence.csv')
        after=pd.read_csv(run/'raw/hfss_convergence_after_fine.csv')
        fine=pd.read_csv(run/'raw/S_fine.csv',index_col=0)
        wide=pd.read_csv(run/'raw/S_wide.csv',index_col=0)
        combined=pd.concat([wide[(wide.index<fine.index.min())|(wide.index>fine.index.max())],fine]).sort_index()
        combined.to_csv(folder/f"raw/HFSS_{row['physical_lead_um']}um_wide_with_fine.csv")
        axes[2].plot(combined.index,20*np.log10(np.maximum(abs(combined.S21.map(complex)),1e-12)),
            label=f"{row['physical_lead_um']} um arms")
        s11=fine.S11.map(complex).to_numpy();s21=fine.S21.map(complex).to_numpy()
        audit=dict(passes=len(conv),last_delta_s=float(conv.iloc[-1]['Max Mag. Delta S']),target_delta_s=.01,
            target_met=bool(conv.iloc[-1]['Max Mag. Delta S']<=.01),
            fine_mesh_reuse_history_unchanged=conv.equals(after),
            max_power_balance_error=float(np.max(abs(abs(s11)**2+abs(s21)**2-1))))
        results.append(dict(**row,frequency_interval_GHz=[float(freq.min()),float(freq.max())],
            spacing_mm=result['spacing_mm'],main_lengths_mm=result['line_lengths_mm'],
            extraction_length_um=result['selected_length_um'],metrics=metrics,solver=audit))
        ax.set_title(f"Actual {row['physical_lead_um']} um CPW arms per side; fixed main lengths")
    axes[2].set_title('HFSS: actual 400 vs 800 um arms (wide data with fine-band replacement)')
    for ax in axes:
        ax.set_xlabel('Frequency (GHz)');ax.set_ylabel('S21 (dB)');ax.grid(alpha=.25);ax.legend()
    fig.tight_layout();fig.savefig(folder/'physical_arms_s21.png',dpi=170)
    summary=dict(cases=results,EM_rerun=False,
        interpretation='Actual geometry changes, not only an extraction-length change. Main lengths fixed.',
        limitations=['Local fine-window RMS values use different frequency windows; do not rank them directly.',
            'Sorted dip minima are not field-identified modes or a measurement of 2g.',
            'Local Q3D coupons omit tee junctions, full-channel coupling and finite bridge loading.',
            'Same adaptive mesh reused within each wide/fine pair, not between changed geometries or solvers.',
            'HFSS convergence, Q3D subtraction sensitivity and Fast-sweep interpolation limit accuracy.'])
    F.write(folder/'summary.json',summary)
    lines=['# Actual 400 / 800 um hammer-arm HFSS S21','',
        'Full Design4 hammer channels were newly simulated with each head-outer-edge to main-line tap arm set to 400 or 800 um. Main filter/resonator lengths remained 4.535956670430899 / 4.852581441726434 mm. Hammer gap/height/width remained 42 / 65 / 25 um; CPW remained 16 / 9 um.', '',
        'The straight-arm change increases main-line centre spacing to 0.892 / 1.692 mm. This also changes distributed inter-line coupling. The same grounding policy is retained (30 automatic bonds and 9 explicit straps); this is simulation geometry, not a fabrication release.', '',
        'Raw full-arm model: the entire corresponding Q3D coupon is a pi network at the two tap planes, without additional arm transmission lines. Deembedded model: subtract matched CPW reference and represent the actual tee-to-hammer-centre arms as transmission lines (0.446 / 0.846 mm each). Both use the same fixed main lengths and unfitted Cin/Cp extraction. Model S21 includes the deterministic renderer port-polarity conversion.', '',
        '| Actual arm per side (um) | HFSS dips (GHz) | Raw-coupon dips (GHz) | Deembedded-model dips (GHz) | Last Delta S |',
        '|---|---|---|---|---|']
    for case in results:
        fmt=lambda key: ', '.join(f'{v:.6f}' for v in case['metrics'][key]['dips_GHz']) or 'not resolved'
        lines.append(f"| {case['physical_lead_um']} | {fmt('HFSS_S21')} | {fmt('raw_full_arm_S21')} | {fmt('deembedded_S21')} | {case['solver']['last_delta_s']:.6f} |")
    lines += ['', 'Requested HFSS Delta S is 0.01. Pass count is not convergence. Fine sweeps reuse each case\'s adaptive setup; different geometries do not share a mesh. Q3D/HFSS cross-solver mesh reuse is not claimed. The 400/800 um Q3D extractions and unchanged Cin/Cp coupon were reused, not solved again.', '',
        'Each model comparison uses its own fine frequency band. Do not compare the two cases by an intersection-only score: the physical arm change moves the resonances to different bands. Sorted dip minima are not field-identified modes and their splitting is not automatically 2g. Q3D shunt subtraction, Fast-sweep interpolation and unmodelled tee/bridge effects remain limitations.', '',
        '[Spectra](physical_arms_s21.png) · [Exact results and raw-run sources](summary.json) · [Run provenance](run.json)', '',
        'The parent controller was started with an earlier snapshot that attempted a common-fine-window ranking after the solves. That ranking, if present, is not used here because its overlap can exclude the resonances. This report supersedes that postprocessing; the actual HFSS raw data are unchanged.']
    (folder/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    finish_run(folder,status='completed',summary=summary)


if __name__=='__main__':main()
