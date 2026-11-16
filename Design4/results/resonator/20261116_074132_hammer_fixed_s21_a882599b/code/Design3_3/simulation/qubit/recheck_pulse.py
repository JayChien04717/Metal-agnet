"""Re-evaluate a saved pulse in the rotating or full time-dependent frame."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse
import csv
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import solve_ivp

from validate_gate import HERE, extracted
sys.path.insert(0,_I('Design3_3/simulation/qubit'))
import sweep as s
import verify as ramp
from sequence import ramp_U

def run(args):
    electrical,beta=extracted(args.matrix)
    s.ec1,s.ect,s.ec2,s.j1,s.j2,s.j12=electrical.values()
    s.beta=beta
    nf,nt=args.fx_levels,args.tc_levels
    on=ramp.TC_CZ_GHZ
    w,_,ids,ops=s.setup(on,nf,nt)
    X=np.einsum('q,qij->ij',beta[:,1],ops)
    path=ramp.ramp_path(ramp.TC_IDLE_GHZ,on,args.ramp_points,nf,nt)
    best=sorted(csv.DictReader(args.pulses.open()),key=lambda r:(int(r['total_ns']),float(r['error'])))
    durations=sorted({int(r['total_ns']) for r in best})
    reports=[]
    for duration in durations:
        for row in [r for r in best if int(r['total_ns'])==duration][:args.count]:
            tr,T,fd=[float(row[k]) for k in ('ramp_ns','mw_ns','drive_GHz')]
            before=ramp_U(path,tr)[:,path[2]]
            after=ramp_U(path,tr,True)[path[2],:]
            order=np.rint(w/float(row['carrier_GHz'])).astype(int)
            Xr=np.where(abs(order[:,None]-order)==1,X,0)
            upper=np.where(order[:,None]-order==1,Xr,0)
            a,_,b,c,d,e=map(float,row['envelope'].split())
            amp=2/(T*float(row['beta_n']))
            def envelope(t):
                k=2*np.pi*t/T
                return amp*np.sin(k/2)**2*(a+b*np.cos(k)+c*np.cos(2*k)+1j*(d*np.sin(k)+e*np.sin(2*k)))
            def rhs(t,y):
                A=envelope(t)
                y=y.reshape(-1,4)
                if args.full:
                    H=w[:,None]*y+np.real(A*np.exp(-2j*np.pi*fd*t))*(X@y)
                else:
                    H=(w-order*fd)[:,None]*y+(A*(upper@y)+A.conjugate()*(upper.conj().T@y))/2
                return (-2j*np.pi*H).ravel()
            evolution=solve_ivp(rhs,(0,T),before.ravel(),method='DOP853',rtol=2e-9,atol=1e-11)
            last=evolution.y[:,-1].reshape(-1,4)
            if not args.full:
                last=np.exp(-2j*np.pi*fd*order*T)[:,None]*last
            U=after@last
            phases=np.angle(np.diag(U))
            local=np.exp(-1j*np.array([phases[0],phases[1],phases[2],phases[1]+phases[2]-phases[0]]))
            trace=np.dot(local*np.array([1,1,1,-1]),np.diag(U))
            error=1-(np.sum(abs(U)**2)+abs(trace)**2)/20
            reports.append(dict(**row,recheck_error=error,recheck_worst_leak=1-min(np.sum(abs(U)**2,axis=0))))
            print(duration,tr,fd,'full' if args.full else 'RWA',error,flush=True)
    out=args.pulses.with_name(args.pulses.stem+('_full.csv' if args.full else '_recheck.csv'))
    with out.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(reports[0]));writer.writeheader();writer.writerows(reports)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('matrix',type=Path)
    parser.add_argument('pulses',type=Path)
    parser.add_argument('--ramp-points',type=int,default=440)
    parser.add_argument('--fx-levels',type=int,default=9)
    parser.add_argument('--tc-levels',type=int,default=7)
    parser.add_argument('--count',type=int,default=1)
    parser.add_argument('--full',action='store_true')
    run(parser.parse_args())
