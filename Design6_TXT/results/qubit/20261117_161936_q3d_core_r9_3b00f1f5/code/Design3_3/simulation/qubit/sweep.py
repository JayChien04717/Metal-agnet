
# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import json
import sys
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linear_sum_assignment

root=_P('.')
sys.path.insert(0,_I('Design3/simulation'))
from zz_analysis import fluxonium, transmon, coupler_EJ
data=json.loads((_P('Design3_2/results/general/legacy/electrical_parameters.json')).read_text(encoding='utf-8'))
ec1,ect,ec2,j1,j2,j12=data['EcJ'].values()
beta=np.array([[data['charge'][q][l]['beta'] for l in ('FX1','TC','FX2')] for q in ('FX1','TC','FX2')])
def setup(f,nf=8,nt=6):
    e1,n1=fluxonium(ec1,4,.6,n_keep=nf)
    e2,n2=fluxonium(ec2,4,.6,n_keep=nf)
    et,nc=transmon(ect,coupler_EJ(ect,f),nt)
    i,t=np.eye(nf),np.eye(nt)
    def k(a,b,c): return np.kron(np.kron(a,b),c)
    ops=np.array([k(n1,t,i),k(i,nc,i),k(i,t,n2)])
    h=np.diag((e1[:,None,None]+et[None,:,None]+e2[None,None,:]).ravel())
    h=h+1e-3*(j1*k(n1,nc,i)+j2*k(i,nc,n2)+j12*k(n1,t,n2))
    w,v=eigh(h,subset_by_index=(0,39))
    bare=[0,1,nt*nf,nt*nf+1]
    ids=linear_sum_assignment(-abs(v[bare,:])**2)[1]
    return w-w[ids[0]],v,ids,np.array([v.conj().T@op@v for op in ops])

def metrics(f,nf=8,nt=6):
    w,v,ids,_=setup(f,nf,nt)
    energies=w[ids]
    # Effective Hamiltonian in the bare single-excitation subspace via polar projection.
    u,_,vh=np.linalg.svd(v[np.ix_([1,nt*nf],ids[1:3])])
    q=u@vh
    exchange=abs((q@np.diag(energies[1:3])@q.conj().T)[0,1])
    return dict(f_GHz=f,zz_kHz=1e6*(energies[3]-energies[2]-energies[1]+energies[0]),exchange_kHz=exchange*1e6,
                comp_overlap=min(abs(v[[0,1,nt*nf,nt*nf+1],ids])**2))

def candidates(f, floor=1e-3):
    w,v,ids,ops=setup(f)
    drive=np.einsum('ql,qij->lij',beta,ops)
    choices=[]
    for i in ids:
        for j in range(len(w)):
            fd=w[j]-w[i]
            if not 3.25<fd<4.4 or max(abs(drive[:,j,i]))<1e-4:
                continue
            for chan,X in enumerate(drive):
                if abs(X[j,i])<1e-4:
                    continue
                others=[(abs(w[b]-w[a]-fd),a,b,abs(X[b,a])) for a in ids for b in range(len(w))
                        if (a,b)!=(i,j) and abs(X[b,a])>floor*abs(X[j,i]) and w[b]>w[a]]
                det,a,b,m=min(others)
                choices.append((dict(on_GHz=f,initial=['00','01','10','11'][list(ids).index(i)],i=i,j=j,line=['FX1','TC','FX2'][chan],
                    carrier_GHz=fd,delta_sel_MHz=1000*det,nearest_strength_ratio=m/abs(X[j,i]),n_TC=abs(ops[1,j,i]),beta_n=abs(X[j,i])),w,ids,X))
    return choices
