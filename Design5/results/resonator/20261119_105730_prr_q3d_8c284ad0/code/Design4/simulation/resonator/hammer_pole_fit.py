"""Shared two-pole fit of complex S11/S21; conditional one-lossy-mode inversion."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.optimize import least_squares
import hammer_deembed as H
from run_record import create_run,finish_run


def fit_poles(f,values,degree=2):
    peaks,props=find_peaks(-20*np.log10(abs(values[:,1])),prominence=3)
    peaks=np.sort(peaks[np.argsort(props['prominences'])[-2:]])
    assert len(peaks)==2
    centre=float(np.mean(f[peaks])); x=(f-centre)*1000
    mask=np.zeros(len(x),bool);mask[::10]=True
    for p in peaks:mask|=abs(f-f[p])<.002
    xx=x[mask];zz=values[mask]
    p0=x[peaks]
    def basis(x,q):
        poles=q[:2]+1j*np.exp(q[2:])
        return np.column_stack([*((x/100)**k for k in range(degree+1)),1/(x-poles[0]),1/(x-poles[1])])
    def residual(q):
        b=basis(xx,q);coeff=np.linalg.lstsq(b,zz,rcond=None)[0]
        r=b@coeff-zz
        return np.r_[r.real.ravel(),r.imag.ravel()]
    bounds=(np.r_[p0-15,np.log([.0001,.0001])],np.r_[p0+15,np.log([100,100])])
    best=None
    for gamma in ((.1,10),(.5,20),(10,.1),(20,.5)):
        fit=least_squares(residual,np.r_[p0,np.log(gamma)],bounds=bounds,x_scale=[1,1,1,1],max_nfev=500)
        if best is None or fit.cost<best.cost:best=fit
    q=best.x;coeff=np.linalg.lstsq(basis(xx,q),zz,rcond=None)[0]
    predicted=basis(x,q)@coeff
    a,b=q[:2]+1j*np.exp(q[2:]);s=a+b;product=a*b
    fr=product.imag/s.imag;fp=s.real-fr;j2=fp*fr-product.real
    result=dict(polynomial_degree=degree,poles_GHz=[[(centre+z.real/1000),z.imag/1000] for z in (a,b)],
        bare_fp_GHz=centre+fp/1000,bare_fr_GHz=centre+fr/1000,kappa_p_MHz=2*s.imag,
        J_MHz=float(np.sqrt(j2)) if j2>0 else None,J_squared_MHz2=float(j2),
        full_grid_rms=float(np.sqrt(np.mean(abs(predicted-values)**2))),
        S21_rms=float(np.sqrt(np.mean(abs(predicted[:,1]-values[:,1])**2))),
        fit_success=bool(best.success),sampling='every tenth background point and all within 2 MHz of dips')
    return result,predicted


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('source',type=Path)
    args=ap.parse_args();source=args.source.resolve()
    folder=create_run('Design4','resonator','hammer_shared_poles',entrypoint=__file__,
        inputs=[source/'run.json',source/'raw/S_fine.csv'],upstream=[str(source)],parameters=dict(
            geometry_changed=False,assumption='real reciprocal J, zero intrinsic readout loss, one lossy Purcell mode',
            convention='H/h frequencies in GHz; reported J in MHz',EM_rerun=False))
    print('RUN_FOLDER='+str(folder),flush=True)
    data=pd.read_csv(source/'raw/S_fine.csv',index_col=0);f=data.index.to_numpy(float)
    values=np.column_stack([data[k].map(complex).to_numpy() for k in ['S11','S21']])
    fits=[]
    for degree in (2,3):
        fit,pred=fit_poles(f,values,degree);fits.append(fit)
        pd.DataFrame({'frequency_GHz':f,'HFSS_S21':values[:,1],'fit_S21':pred[:,1]}).to_csv(folder/f'fit_degree{degree}.csv',index=False)
    summary=dict(source=str(source),fits=fits,
        J_relative_degree_sensitivity=abs(fits[1]['J_MHz']-fits[0]['J_MHz'])/fits[0]['J_MHz'],
        interpretation='Conditional effective coupling derived from fitted poles. No physical length is retuned.',
        limitation='Additional damping channels, dissipative coupling, unresolved modes and EM mesh error can invalidate the inversion.')
    H.write(folder/'pole_fit.json',summary)
    finish_run(folder,status='completed',summary=summary)
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
