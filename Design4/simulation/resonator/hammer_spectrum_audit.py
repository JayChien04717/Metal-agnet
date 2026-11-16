"""Audit lossless two-port power and an explicitly archived native adaptive point."""
from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd
import hammer_deembed as H
from run_record import create_run, finish_run


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source',type=Path)
    ap.add_argument('--adaptive-point',type=Path)
    args=ap.parse_args();source=args.source.resolve()
    inputs=[source/'raw/S_fine.csv',source/'raw/hfss_convergence.csv']
    if args.adaptive_point:
        args.adaptive_point=args.adaptive_point.resolve();inputs.append(args.adaptive_point)
    folder=create_run('Design4','resonator','hammer_spectrum_audit',entrypoint=__file__,inputs=inputs,
        upstream=[str(source)],parameters=dict(EM_rerun=False,geometry_changed=False,
        power_assumption='Lossless two-port approximation; residual is a numerical diagnostic, not mesh convergence'))
    data=pd.read_csv(source/'raw/S_fine.csv',index_col=0)
    f=data.index.to_numpy(float);s11=data.S11.map(complex).to_numpy();s21=data.S21.map(complex).to_numpy()
    power=abs(s11)**2+abs(s21)**2
    result=dict(max_two_port_power_residual=float(np.max(abs(power-1))),frequency_GHz=[float(f.min()),float(f.max())])
    if args.adaptive_point:
        native=args.adaptive_point.read_text()
        frequency=float(re.search(r'Frequency=([\deE+.-]+)',native).group(1))/1e9
        pass_no=int(re.search(r'AdaptivePass=(\d+)',native).group(1))
        block=native.split("$begin 'SMatrix'",1)[1].split("$end 'SMatrix'",1)[0]
        numbers=re.findall(r'c\(([^,]+),\s*([^)]+)\)',block)
        assert len(numbers)==4
        matrix=np.array([complex(float(a),float(b)) for a,b in numbers]).reshape(2,2)
        conv=pd.read_csv(source/'raw/hfss_convergence.csv')
        assert pass_no==int(conv.iloc[-1]['Pass Number']), 'Native point must use the final adaptive mesh'
        ix=int(np.argmin(abs(f-frequency)))
        assert abs(f[ix]-frequency)<1e-8, 'Require a sampled matching frequency; do not interpolate validation data'
        result['adaptive_point']=dict(source=str(args.adaptive_point),frequency_GHz=frequency,pass_number=pass_no,
            S11_error=float(abs(s11[ix]-matrix[0,0])),S21_error=float(abs(s21[ix]-matrix[1,0])),
            native_S11=str(matrix[0,0]),native_S21=str(matrix[1,0]),
            limitation='One independently solved adaptive-frequency point only; no direct checks at both resonance minima')
    pd.DataFrame({'frequency_GHz':f,'two_port_power':power,'absolute_residual':abs(power-1)}).to_csv(folder/'power.csv',index=False)
    H.write(folder/'audit.json',result)
    finish_run(folder,status='completed',summary=result)
    print('RUN_FOLDER='+str(folder),flush=True);print(result,flush=True)


if __name__=='__main__':main()
