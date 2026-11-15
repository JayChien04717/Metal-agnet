"""Optimize a short TC I/Q envelope including both flux ramps."""
import argparse
import csv
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize

from validate_gate import extracted, HERE
sys.path.insert(0, str(HERE/'gate_model'))
import sweep as s
import verify as ramp
from sequence import ramp_U

def run(args):
    e,beta=extracted(args.matrix)
    s.ec1,s.ect,s.ec2,s.j1,s.j2,s.j12=e.values()
    s.beta=beta
    on=ramp.TC_CZ_GHZ
    nf,nt=args.fx_levels,args.tc_levels
    energies,_,ids,ops=s.setup(on,nf,nt)
    setup=s.setup
    s.setup=lambda f:setup(f,nf,nt)
    try:
        cases=[c for c,*_ in s.candidates(on) if c['line']=='TC']
    finally:
        s.setup=setup
    cases=sorted(cases,key=lambda r:r['delta_sel_MHz'],reverse=True)
    X=np.einsum('q,qij->ij',beta[:,1],ops)
    path=ramp.ramp_path(ramp.TC_IDLE_GHZ,on,args.ramp_points,nf,nt)
    rows=[]
    for tr in args.ramp_ns:
        before=ramp_U(path,tr)[:,path[2]]
        after=ramp_U(path,tr,True)[path[2],:]
        for gate in args.gate_ns:
            T=gate-2*tr
            for case in cases[:args.candidates]:
                i,j=case['i'],case['j']
                if abs(X[j,i]) < 1e-6:
                    continue
                carrier=energies[j]-energies[i]
                order=np.rint(energies/carrier).astype(int)
                coupling=np.where(abs(order[:,None]-order)==1,X,0)
                lam,V=eigh(coupling)
                Vh=V.conj().T
                peak=2/(T*abs(X[j,i]))
                steps=max(55,int(T*args.steps_per_ns))
                t=(np.arange(steps)+.5)/steps
                env=np.sin(np.pi*t)**2
                basis=np.array([env,env*np.cos(2*np.pi*t),env*np.cos(4*np.pi*t),
                                1j*env*np.sin(2*np.pi*t),1j*env*np.sin(4*np.pi*t)])
                dt=T/steps
                def cost(z, extra=False):
                    fd=carrier+z[1]/T
                    wave=peak*(np.r_[z[0],z[2:]] @ basis)
                    d=np.exp(-1j*np.pi*(energies-order*fd)*dt)
                    state=before.copy()
                    for a in wave:
                        phase=np.exp(1j*np.angle(a)*order)
                        state=(d*phase)[:,None]*(V@(np.exp(-1j*np.pi*lam*abs(a)*dt)[:,None]*(Vh@((d*phase.conj())[:,None]*state))))
                    U=after@(np.exp(-2j*np.pi*order*fd*T)[:,None]*state)
                    phases=np.angle(np.diag(U))
                    local=np.exp(-1j*np.array([phases[0],phases[1],phases[2],phases[1]+phases[2]-phases[0]]))
                    overlap=abs(np.dot(local*np.array([1,1,1,-1]),np.diag(U)))**2
                    error=1-(np.sum(abs(U)**2)+overlap)/20
                    if extra:
                        return error,1-min(np.sum(abs(U)**2,axis=0)),max(abs(wave))*4.135667696e-3/2
                    return error
                x0=np.array([1.,0.,0.,0.,0.,0.])
                bounds=[(.65,1.45),(-.7,.7)]+[(-.65,.65)]*4
                if args.envelope is not None:
                    x0 = np.array(args.envelope)
                fit=minimize(cost,x0,method='L-BFGS-B',bounds=bounds,
                             options=dict(maxiter=85,ftol=1e-11,gtol=1e-7))
                err,leak,voltage=cost(fit.x,True)
                rows.append(dict(**case,total_ns=gate,ramp_ns=tr,mw_ns=T,error=err,worst_leak=leak,
                                 drive_GHz=carrier+fit.x[1]/T,drive_peak_mV=voltage,
                                 envelope=' '.join(f'{z:.10g}' for z in fit.x)))
        print(tr,min(rows,key=lambda r:r['error']),flush=True)
    out=args.matrix.with_name(f'{args.matrix.stem}_iq_{args.ramp_points}steps.csv')
    with out.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('matrix',type=Path)
    ap.add_argument('--ramp-points',type=int,default=240)
    ap.add_argument('--fx-levels',type=int,default=8)
    ap.add_argument('--tc-levels',type=int,default=6)
    ap.add_argument('--steps-per-ns',type=float,default=3.)
    ap.add_argument('--gate-ns',type=int,nargs='+',default=[31])
    ap.add_argument('--ramp-ns',type=int,nargs='+',default=[2,3,4,5])
    ap.add_argument('--candidates',type=int,default=4)
    ap.add_argument('--envelope',type=float,nargs=6)
    run(ap.parse_args())
