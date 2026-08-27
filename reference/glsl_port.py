"""Direct port of principia-ii/src/shaders/principia/frag.glsl decodeIC + integrator."""
import numpy as np
PI=np.pi; ALPHA_MIN=0.05; MU_MAX=5.0; Q_MAX=2.0; G=1.0
sig=lambda x:1/(1+np.exp(-x))
def decode(z):
    # z: (...,10) GLSL indexing
    mu1=MU_MAX*(2*sig(z[...,8])-1); mu2=MU_MAX*(2*sig(z[...,9])-1)
    l=np.stack([np.zeros_like(mu1),mu1,mu2],-1); e=np.exp(l-l.max(-1,keepdims=True)); m=e/e.sum(-1,keepdims=True)
    M01=m[...,0]+m[...,1]
    alpha=ALPHA_MIN+(PI/2-2*ALPHA_MIN)*sig(z[...,1])
    beta=PI*sig(z[...,0])
    muR=m[...,0]*m[...,1]/M01; muL=m[...,2]*M01
    rt=np.stack([np.cos(alpha),np.zeros_like(alpha)],-1)
    lt=np.stack([np.sin(alpha)*np.cos(beta),np.sin(alpha)*np.sin(beta)],-1)
    rho=rt/np.sqrt(muR)[...,None]; lam=lt/np.sqrt(muL)[...,None]
    r01=-m[...,2:3]*lam
    r0=r01-(m[...,1:2]/M01[...,None])*rho
    r1=r01+(m[...,0:1]/M01[...,None])*rho
    r2=M01[...,None]*lam
    pR=Q_MAX*(2*sig(z[...,4:6])-1); pL=Q_MAX*(2*sig(z[...,6:8])-1)
    p0=-pR-(m[...,0:1]/M01[...,None])*pL
    p1= pR-(m[...,1:2]/M01[...,None])*pL
    p2= pL
    return m, np.stack([r0,r1,r2],-2), np.stack([p0,p1,p2],-2)

def shape_n(r,m):
    M01=m[...,0]+m[...,1]
    rho=r[...,1,:]-r[...,0,:]
    com=(m[...,0:1]*r[...,0,:]+m[...,1:2]*r[...,1,:])/M01[...,None]
    lam=r[...,2,:]-com
    muR=m[...,0]*m[...,1]/M01; muL=m[...,2]*M01
    rt=np.sqrt(muR)[...,None]*rho; lt=np.sqrt(muL)[...,None]*lam
    A=(rt*rt).sum(-1); B=(lt*lt).sum(-1); I=np.maximum(A+B,1e-30)
    p=(rt*lt).sum(-1); q=rt[...,1]*lt[...,0]-rt[...,0]*lt[...,1]
    return np.stack([(A-B)/I,2*p/I,2*q/I],-1)

def run(z, t_max=13.0, n_macro=64, n_sub=32, r_coll=1e-3, stop_escape=False):
    m,r,p=decode(z)
    r=r.copy(); p=p.copy()
    dt=t_max/(n_macro*n_sub)
    alive=np.ones(r.shape[:-2],bool); coll=np.zeros(r.shape[:-2],bool)
    tend=np.full(r.shape[:-2],t_max)
    PAIRS=[(0,1),(0,2),(1,2)]
    def acc(r):
        f=np.zeros_like(r)
        for i,j in PAIRS:
            d=r[...,j,:]-r[...,i,:]; dd=np.sqrt((d*d).sum(-1))[...,None]
            fm=np.where(dd>1e-10, G*(m[...,i]*m[...,j])[...,None]/np.maximum(dd,1e-10)**3, 0.0)
            f[...,i,:]+=fm*d; f[...,j,:]-=fm*d
        return f
    t=0.0
    for step in range(n_macro*n_sub):
        f=acc(r); p=p+f*dt*0.5
        r=r+p/m[...,None]*dt
        f=acc(r); p=p+f*dt*0.5
        t+=dt
        dmin=np.min([np.sqrt(((r[...,j,:]-r[...,i,:])**2).sum(-1)) for i,j in PAIRS],axis=0)
        hit=(dmin<r_coll)&alive
        if hit.any():
            tend=np.where(hit,t,tend); coll|=hit; alive&=~hit
        if not alive.any(): break
    return m,r,p,coll,tend
