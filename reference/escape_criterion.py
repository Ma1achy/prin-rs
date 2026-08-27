"""GLSL port + the closure-and-energy escape criterion, with stop-on-escape as a toggle."""
import numpy as np, glsl_port as GP
PAIRS=[(0,1),(0,2),(1,2)]
def _acc(r,m):
    f=np.zeros_like(r)
    for i,j in PAIRS:
        d=r[...,j,:]-r[...,i,:]; dd=np.sqrt((d*d).sum(-1))[...,None]
        fm=np.where(dd>1e-10,(m[...,i]*m[...,j])[...,None]/np.maximum(dd,1e-10)**3,0.)
        f[...,i,:]+=fm*d; f[...,j,:]-=fm*d
    return f
def _rel(r,p,m):
    """Per-body: relative energy to the other two, and separation."""
    v=p/m[...,None]; E=[];D=[]
    for b in range(3):
        o=[k for k in range(3) if k!=b]; mb=m[...,o[0]]+m[...,o[1]]
        rc=(m[...,o[0],None]*r[...,o[0],:]+m[...,o[1],None]*r[...,o[1],:])/mb[...,None]
        vc=(m[...,o[0],None]*v[...,o[0],:]+m[...,o[1],None]*v[...,o[1],:])/mb[...,None]
        dr=r[...,b,:]-rc; dv=v[...,b,:]-vc
        d=np.sqrt((dr*dr).sum(-1))
        E.append(.5*(dv*dv).sum(-1)-mb/np.maximum(d,1e-12)); D.append(d)
    return np.stack(E,-1),np.stack(D,-1)

def render(z, t_max=26.0, dt=0.002, stop_on_escape=False, tau=1e-3, win=0.4,
           r_coll=1e-3, sample_every=5):
    """
    ESCAPE  <=>  |dn| over `win` < tau   AND   E_rel > 0
    Measured 100% precision on the config chart; `receding` and `d > r_esc`
    were both redundant once these two hold (identical to the digit).

    stop_on_escape=False  -> integrate all to t_max, shape at a COMMON time
    stop_on_escape=True   -> freeze at t_end, reproduces the patchwork artefact
    """
    m,r,p=GP.decode(z); r=r.copy(); p=p.copy()
    sh=z.shape[:-1]
    nbuf=int(round(win/(dt*sample_every)))+1
    buf=[]
    alive=np.ones(sh,bool)
    esc_body=np.full(sh,-1,np.int8); t_end=np.full(sh,np.nan)
    frozen=np.zeros(sh+(3,)); coll=np.zeros(sh,bool)
    for s in range(int(t_max/dt)):
        f=_acc(r,m); p=p+f*dt*.5; r=r+p/m[...,None]*dt; f=_acc(r,m); p=p+f*dt*.5
        if s%sample_every: continue
        t=(s+1)*dt
        n=GP.shape_n(r,m); buf.append(n)
        if len(buf)>nbuf: buf.pop(0)
        ds=np.stack([np.sqrt(((r[...,j,:]-r[...,i,:])**2).sum(-1)) for i,j in PAIRS],-1)
        hit=(ds.min(-1)<r_coll)&alive
        if hit.any():
            frozen[hit]=n[hit]; t_end=np.where(hit,t,t_end); coll|=hit; alive&=~hit
        if len(buf)<nbuf: continue
        E,D=_rel(r,p,m)
        dn=np.sqrt(((buf[-1]-buf[0])**2).sum(-1))[...,None]
        fire=(dn<tau)&(E>0)                      # THE CRITERION
        any_f=fire.any(-1)&alive&(esc_body<0)
        if any_f.any():
            b=np.argmax(fire,-1)
            esc_body=np.where(any_f,b,esc_body); t_end=np.where(any_f,t,t_end)
            if stop_on_escape:
                frozen[any_f]=n[any_f]; alive&=~any_f
        if not alive.any(): break
    frozen[alive]=GP.shape_n(r,m)[alive]
    return dict(n=frozen,esc_body=esc_body,t_end=t_end,coll=coll,
                frozen_frac=float((~alive).mean()))
