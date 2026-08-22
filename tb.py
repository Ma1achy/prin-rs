"""
Planar three-body: vectorised KDK leapfrog over a batch of initial conditions.
Purpose: probe whether a LOCAL uncertainty-exponent estimate separates
"refining helps" from "refining is futile" regions of IC space.

Pythagorean (Burrau) problem: m=(3,4,5) at rest at the vertices of a 3-4-5 triangle.
"""
import numpy as np

G = 1.0
M = np.array([3.0, 4.0, 5.0])
R0 = np.array([[1.0, 3.0], [-2.0, -1.0], [1.0, -1.0]])
V0 = np.zeros((3, 2))

PAIRS = [(0, 1), (0, 2), (1, 2)]


def accel(r, eps2):
    """r: (N,3,2) -> a: (N,3,2). Softened Newtonian gravity."""
    a = np.zeros_like(r)
    for i, j in PAIRS:
        d = r[:, j, :] - r[:, i, :]                  # (N,2)
        d2 = np.einsum('ij,ij->i', d, d) + eps2
        inv3 = d2 ** -1.5
        f = d * inv3[:, None]
        a[:, i, :] += G * M[j] * f
        a[:, j, :] -= G * M[i] * f
    return a


def energy(r, v, eps2):
    ke = 0.5 * np.einsum('k,nki->n', M, v * v)
    pe = np.zeros(r.shape[0])
    for i, j in PAIRS:
        d = r[:, j, :] - r[:, i, :]
        dist = np.sqrt(np.einsum('ij,ij->i', d, d) + eps2)
        pe -= G * M[i] * M[j] / dist
    return ke + pe


def pair_dists(r):
    """(N,3) distances in PAIRS order."""
    out = np.empty((r.shape[0], 3))
    for k, (i, j) in enumerate(PAIRS):
        d = r[:, j, :] - r[:, i, :]
        out[:, k] = np.sqrt(np.einsum('ij,ij->i', d, d))
    return out


def integrate(r0, v0, t_max, dt, eps=1e-3, word_len=24, sample_every=25):
    """
    KDK leapfrog. Returns dict with final state, min distances, energy drift,
    and a symbolic 'word': the sequence of changes in which pair is tightest.
    """
    r = r0.copy()
    v = v0.copy()
    eps2 = eps * eps
    n = r.shape[0]
    E0 = energy(r, v, eps2)

    dmin = np.full(n, np.inf)
    word = np.zeros((n, word_len), dtype=np.int8)
    wlen = np.zeros(n, dtype=np.int32)
    prev_tight = np.argmin(pair_dists(r), axis=1).astype(np.int8)

    a = accel(r, eps2)
    steps = int(round(t_max / dt))
    for s in range(steps):
        v += 0.5 * dt * a
        r += dt * v
        a = accel(r, eps2)
        v += 0.5 * dt * a

        if s % sample_every == 0:
            pd = pair_dists(r)
            dmin = np.minimum(dmin, pd.min(axis=1))
            tight = np.argmin(pd, axis=1).astype(np.int8)
            changed = (tight != prev_tight) & (wlen < word_len)
            if changed.any():
                idx = np.nonzero(changed)[0]
                word[idx, wlen[idx]] = tight[idx] + 1     # symbols 1..3
                wlen[idx] += 1
            prev_tight = tight

    E1 = energy(r, v, eps2)
    drift = np.abs((E1 - E0) / np.abs(E0))
    return dict(r=r, v=v, dmin=dmin, drift=drift, word=word, wlen=wlen)


def classify(res):
    """
    Coarse Chazy-like outcome: identify the most-bound pair; the remaining body
    is the candidate escaper. Class = which body is escaping (0,1,2), or 3 if
    the triple still looks bound (no body with positive energy wrt the others).
    """
    r, v = res['r'], res['v']
    n = r.shape[0]
    pd = pair_dists(r)
    tight = np.argmin(pd, axis=1)                    # index into PAIRS
    third = np.array([2, 1, 0])[tight]               # body not in the tight pair

    out = np.full(n, 3, dtype=np.int8)
    for b in range(3):
        sel = third == b
        if not sel.any():
            continue
        others = [k for k in range(3) if k != b]
        m_bin = M[others].sum()
        # binary barycentre
        rc = (M[others[0]] * r[sel][:, others[0], :] +
              M[others[1]] * r[sel][:, others[1], :]) / m_bin
        vc = (M[others[0]] * v[sel][:, others[0], :] +
              M[others[1]] * v[sel][:, others[1], :]) / m_bin
        dr = r[sel][:, b, :] - rc
        dv = v[sel][:, b, :] - vc
        dist = np.sqrt(np.einsum('ij,ij->i', dr, dr))
        sp = 0.5 * np.einsum('ij,ij->i', dv, dv) - G * m_bin / np.maximum(dist, 1e-9)
        receding = np.einsum('ij,ij->i', dr, dv) > 0
        esc = (sp > 0) & receding
        cls = np.where(esc, b, 3).astype(np.int8)
        out[sel] = cls
    return out


def burrau_grid(nx, ny, cx, cy, half, body=0, ens=0, jitter_frac=0.5, seed=0):
    """
    2D slice: perturb `body`'s initial position over a box of half-width `half`
    centred on (cx,cy). If ens>0, add that many jittered copies per grid point
    (jitter ~ jitter_frac of the cell size) -> the SSAA/ensemble bundle.
    Returns r0, v0, and the grid index of each row.
    """
    xs = np.linspace(cx - half, cx + half, nx)
    ys = np.linspace(cy - half, cy + half, ny)
    X, Y = np.meshgrid(xs, ys, indexing='xy')
    base = np.stack([X.ravel(), Y.ravel()], axis=1)          # (nx*ny, 2)
    npts = base.shape[0]

    hx = (2 * half) / max(nx - 1, 1)
    rng = np.random.default_rng(seed)
    reps = 1 + ens
    pos = np.repeat(base, reps, axis=0)                      # (npts*reps, 2)
    if ens > 0:
        jit = rng.uniform(-jitter_frac * hx, jitter_frac * hx, size=pos.shape)
        mask = np.ones(pos.shape[0], dtype=bool)
        mask[::reps] = False                                  # copy 0 = nominal
        pos[mask] += jit[mask]

    n = pos.shape[0]
    r0 = np.tile(R0, (n, 1, 1))
    v0 = np.tile(V0, (n, 1, 1))
    r0[:, body, :] = pos
    gid = np.repeat(np.arange(npts), reps)
    return r0, v0, gid, (nx, ny), hx
