
# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)



import numpy as np
from scipy.linalg import eigh
def ramp_U(path,T,rev=False):
    O,W,ids,on,*_=path
    if rev:
        O=[a.conj().T for a in O[::-1]]
        W=W[::-1]
    U=O[0]
    for o,w in zip(O[1:],W):
        U=o@(np.exp(-2j*np.pi*w*T/len(W))[:,None]*U)
    return U

def flat_U(w,X,carrier,T,amp):
    order=np.rint(w/carrier).astype(int)
    X=np.where(abs(order[:,None]-order)==1,X,0)
    en=w-order*carrier
    dx,Q=eigh(X)
    steps=max(80,int(T*3))
    dt=T/steps
    d=np.exp(-1j*np.pi*en*dt)
    U=np.eye(len(w),dtype=complex)
    for j in range(steps):
        a=amp*np.sin(np.pi*(j+.5)/steps)**2
        U=d[:,None]*(Q@(np.exp(-1j*np.pi*dx*a*dt)[:,None]*(Q.conj().T@(d[:,None]*U))))
    return np.exp(-2j*np.pi*order*carrier*T)[:,None]*U
