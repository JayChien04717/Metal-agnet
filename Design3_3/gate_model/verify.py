import sweep as s

from sweep import *
from functools import lru_cache


def set_deng():
    global FT
    s.ec1,s.ec2,s.ect=1.17,1.30,.21
    s.j1,s.j2,s.j12=309,-371,-68
    FT=(4.2,.6,4.,.55)

FT=(4.,.6,4.,.6)
original_fluxonium=lru_cache(None)(s.fluxonium)
def flux(e,*args,**kwargs):
    a,b=FT[:2] if abs(e-s.ec1)<1e-9 else FT[2:]
    return original_fluxonium(e,a,b,**kwargs)
s.fluxonium=flux

@lru_cache(None)
def ramp_path(off,on,points,nf=7,nt=6,states=40):
    a=coupler_EJ(s.ect,off)
    b=coupler_EJ(s.ect,on)
    jj1,jj2=(a+b)/2,(a-b)/2
    bare_n=np.arange(-30,31)
    def at(phi):
        z=jj1*np.exp(1j*np.pi*phi)+jj2*np.exp(-1j*np.pi*phi)
        E=abs(z)
        f=off if phi==0 else on if phi==.5 else np.diff(transmon(s.ect,E,2)[0])[0]
        E=coupler_EJ(s.ect,f)
        w,v,ids,ops=s.setup(f,nf,nt)
        H=np.diag(4 * s.ect * bare_n ** 2.0)-E / 2 * (np.eye(61, k=1) + np.eye(61, k=-1))
        _,C=eigh(H)
        C=np.exp(1j*np.angle(z)*bare_n[:,None])*C[:,:nt]
        return w[:states],v[:,:states],ids,C
    mid=[at(.25*(1-np.cos(np.pi*(j+.5)/points))) for j in range(points)]
    initial=at(0.)
    end=at(.5)
    def overlap(new,old):
        return new[1].conj().T@(np.einsum('kl,albs->akbs',new[3].conj().T@old[3],old[1].reshape(nf,nt,nf,-1),optimize=True).reshape(nf*nt*nf,-1))
    nodes=[initial]+mid+[end]
    overlaps=[overlap(q,p) for p,q in zip(nodes[:-1],nodes[1:])]
    return overlaps,[x[0] for x in mid],initial[2],end[2],jj1,jj2

def ramp_error(path,T):
    overlaps,w,ids,on,*_=path
    U=np.eye(len(w[0]),dtype=complex)[:,ids]
    for O,E in zip(overlaps,w):
        U=np.exp(-2j*np.pi*E*T/len(w))[:,None]*(O@U)
    U=overlaps[-1]@U
    return max(0,1-min(abs(U[on,range(4)])**2))
