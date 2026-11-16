"""Compare poles of the fixed Q3D circuit and actual HFSS, optionally refining Cin/Cp inputs."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import hammer_deembed as H
import hammer_fixed_s21 as F
from hammer_pole_fit import fit_poles
from run_record import create_run,finish_run


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('model_run',type=Path);ap.add_argument('--hfss',type=Path)
    ap.add_argument('--calibration',type=Path)
    ap.add_argument('--geometry-run',type=Path,help='Explicit original geometry manifest when historical comparison lacks spacing')
    ap.add_argument('--upper-ghz',type=float,default=6.25)
    args=ap.parse_args();source=args.model_run.resolve();hfss=(args.hfss or source).resolve()
    metrics=json.loads((source/'comparison.json').read_text());assert 'port_polarity' in metrics
    cap=metrics['cin_cp_fF'];inputs=[source/'comparison.json',hfss/'raw/S_fine.csv']
    upstream=[str(source),str(hfss)]
    spacing=metrics.get('spacing_mm')
    if args.geometry_run:
        geometry_path=args.geometry_run.resolve()/'run.json'
        geometry=json.loads(geometry_path.read_text())['parameters']
        if spacing is not None:assert np.isclose(spacing,geometry['spacing_mm'])
        spacing=geometry['spacing_mm'];inputs.append(geometry_path);upstream.append(str(args.geometry_run.resolve()))
        assert np.allclose(metrics['line_lengths_mm'],geometry['lengths_mm'])
    assert spacing is not None, 'Supply --geometry-run for historical comparisons without recorded spacing'
    if args.calibration:
        cpath=args.calibration.resolve()/'raw/cp_coupon/caps.json'
        c=json.loads(cpath.read_text());cap=[c['mutual_fF'],*c['ground_fF']]
        inputs.append(cpath);upstream.append(str(args.calibration.resolve()))
    data=pd.read_csv(hfss/'raw/S_fine.csv',index_col=0)
    data=data.loc[data.index<=args.upper_ghz]
    f=data.index.to_numpy(float);values=np.column_stack([data[k].map(complex).to_numpy() for k in ('S11','S21')])
    observed,_=fit_poles(f,values,3)
    predictions={};export=pd.DataFrame({'frequency_GHz':f,'HFSS_S21':values[:,1]})
    for name,choice in metrics['models'].items():
        if name=='legacy_no_stubs':continue
        halfspan=spacing/2-choice['tee_to_cap_reference_mm']
        omit=choice['tee_to_cap_reference_mm']==0
        a,b=F.model(f,cap,choice['hammer_fF'],metrics['line_lengths_mm'],metrics['Z0_ohm'],metrics['velocity_m_s'],
            metrics['feed_lengths_mm'],spacing=spacing,cap_half_span=halfspan,omit_hammer_stubs=omit)
        b=b*metrics['port_polarity']['S21_multiplier']
        predicted,_=fit_poles(f,np.column_stack([a,b]),3)
        predictions[name]=dict(poles=predicted,rms_complex_S21=float(np.sqrt(np.mean(abs(b-values[:,1])**2))),
            J_relative_error=predicted['J_MHz']/observed['J_MHz']-1,
            bare_fp_error_MHz=(predicted['bare_fp_GHz']-observed['bare_fp_GHz'])*1000,
            bare_fr_error_MHz=(predicted['bare_fr_GHz']-observed['bare_fr_GHz'])*1000)
        export[name+'_S21']=b
    folder=create_run('Design4','resonator','hammer_model_poles',entrypoint=__file__,inputs=inputs,upstream=upstream,
        parameters=dict(EM_rerun=False,circuit_parameters_fitted=False,frequency_band_GHz=[float(f.min()),float(f.max())],
            cap_fF=cap,original_cap_fF=metrics['cin_cp_fF'],spacing_mm=spacing,
            pole_assumption='One lossy Purcell mode and negligible intrinsic readout loss; pole fits characterize each spectrum, not retune the circuit'))
    result=dict(HFSS_poles=observed,models=predictions,
        limitations=['Shared-pole inversion is conditional, not a field-based coupling extraction.',
            'EM convergence, CPW propagation approximation and local-coupon boundaries remain separate uncertainties.'])
    export.to_csv(folder/'model_comparison.csv',index=False)
    H.write(folder/'comparison.json',result);finish_run(folder,status='completed',summary=result)
    print('RUN_FOLDER='+str(folder),flush=True);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
