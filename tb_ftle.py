"""
Benettin FTLE + diffusion accumulators on top of tb.py.

FTLE: each trajectory carries a shadow at phase-space separation d0, renormalised every
`renorm_every` steps; S accumulates log(d/d0); ftle = S / T.  Renormalisation is what
stops it saturating (see findings doc §5.1).

Diffusion: running linear regression of log(moment of inertia) against t, from O(1)
accumulators (n, St, Stt, Sy, Sty) -- a stand-in for whatever the payload spec's
`mean_y`/`C_ty` regress; the point here is the *shape* of the quantity, not its exact choice.
"""
import numpy as np
import tb


def _phase_dist(dr, dv):
    return np.sqrt(np.einsum('nki,nki->n', dr, dr) + np.einsum('nki,nki->n', dv, dv))


def inertia(r):
    """Mass-weighted moment of inertia about the centre of mass."""
    m = tb.M[None, :, None]
    com = (m * r).sum(axis=1, keepdims=True) / tb.M.sum()
    d = r - com
    return np.einsum('k,nki->n', tb.M, d * d)


def integrate_full(r0, v0, t_max, dt, eps=0.03, d0=1e-8, renorm_every=200, seed=0):
    """Main trajectory + Benettin shadow + diffusion regression. Returns dict of fields."""
    eps2 = eps * eps
    r = r0.copy(); v = v0.copy()
    n = r.shape[0]

    rng = np.random.default_rng(seed)
    pert = rng.normal(size=r.shape); pert /= np.linalg.norm(pert.reshape(n, -1), axis=1)[:, None, None]
    rs = r + d0 * pert
    vs = v.copy()

    S = np.zeros(n)                      # accumulated log stretch
    nre = np.zeros(n, dtype=np.int64)    # completed renormalisations

    # diffusion regression accumulators (O(1), no history)
    cnt = 0.0; St = 0.0; Stt = 0.0
    Sy = np.zeros(n); Sty = np.zeros(n)

    dmin = np.full(n, np.inf)
    a = tb.accel(r, eps2); as_ = tb.accel(rs, eps2)
    steps = int(round(t_max / dt))

    for s in range(steps):
        v += 0.5 * dt * a;   r += dt * v;   a = tb.accel(r, eps2);   v += 0.5 * dt * a
        vs += 0.5 * dt * as_; rs += dt * vs; as_ = tb.accel(rs, eps2); vs += 0.5 * dt * as_

        if s % renorm_every == 0 and s > 0:
            dr = rs - r; dv = vs - v
            d = _phase_dist(dr, dv)
            ok = d > 1e-300
            S[ok] += np.log(d[ok] / d0)
            nre[ok] += 1
            sc = np.where(ok, d0 / np.maximum(d, 1e-300), 1.0)[:, None, None]
            rs = r + dr * sc; vs = v + dv * sc

        if s % 25 == 0:
            t = (s + 1) * dt
            y = np.log(np.maximum(inertia(r), 1e-300))
            cnt += 1.0; St += t; Stt += t * t; Sy += y; Sty += t * y
            dmin = np.minimum(dmin, tb.pair_dists(r).min(axis=1))

    T = steps * dt
    ftle = np.where(nre > 0, S / np.maximum(T, 1e-12), np.nan)
    den = cnt * Stt - St * St
    diffusion = np.where(abs(den) > 1e-12, (cnt * Sty - St * Sy) / den, np.nan)
    return dict(r=r, v=v, ftle=ftle, diffusion=diffusion, dmin=dmin, nre=nre, S=S)
