"""
TEST: does  ratio = between / (between + within)  separate
"refining helps" from "refining is futile"?

within  = spread inside one bundle (E+1 copies jittered within a pixel footprint)
          -> irreducible local sensitivity. Splitting cannot reduce it.
between = spread across the quad's N x N nominal samples
          -> spatial heterogeneity. Splitting CAN reduce it.

Prediction if the criterion works:
  - resolvable boundary : ratio HIGH, and FALLS as we refine (boundary gets resolved)
  - chaotic sea         : ratio ~0.5 or LOWER and FLAT (within saturates too) -> floor
  - smooth region       : both tiny -> ratio meaningless, must be gated on magnitude

Failure mode to look for (Wada): sea keeps between HIGH at every scale while
within stays low -> ratio stays high -> refine forever.
"""
import numpy as np
import tb

MT = tb.M
MTOT = MT.sum()


def shape_vec(r):
    """(N,3,2) -> (N,3) unit vectors on the shape sphere (Hopf map of mass-weighted Jacobi)."""
    m0, m1, m2 = MT
    rho = r[:, 1, :] - r[:, 0, :]
    com01 = (m0 * r[:, 0, :] + m1 * r[:, 1, :]) / (m0 + m1)
    lam = r[:, 2, :] - com01
    mu_rho = m0 * m1 / (m0 + m1)
    mu_lam = m2 * (m0 + m1) / MTOT
    rt = np.sqrt(mu_rho) * rho
    lt = np.sqrt(mu_lam) * lam
    a = np.einsum('ij,ij->i', rt, rt)
    b = np.einsum('ij,ij->i', lt, lt)
    I = a + b
    p = rt[:, 0] * lt[:, 0] + rt[:, 1] * lt[:, 1]
    q = rt[:, 1] * lt[:, 0] - rt[:, 0] * lt[:, 1]
    n = np.stack([(a - b) / I, 2 * p / I, 2 * q / I], axis=1)
    return n / np.linalg.norm(n, axis=1, keepdims=True)


def svar(n):
    """Spherical variance of unit vectors: 1 - |resultant|, in [0,1]."""
    if n.shape[0] < 2:
        return 0.0
    return float(1.0 - np.linalg.norm(n.mean(axis=0)))


def disagree(c):
    """Fraction of samples not in the modal class, in [0,1]."""
    if len(c) < 2:
        return 0.0
    _, cnt = np.unique(c, return_counts=True)
    return float(1.0 - cnt.max() / len(c))


def probe(cx, cy, half, N=8, ens=3, t_max=25.0, dt=1e-3, eps=1e-2, body=0, seed=0):
    """One quad: return within/between/ratio for both the continuous and outcome measures."""
    r0, v0, gid, _, hx = tb.burrau_grid(N, N, cx, cy, half, body=body,
                                        ens=ens, jitter_frac=0.5, seed=seed)
    res = tb.integrate(r0, v0, t_max=t_max, dt=dt, eps=eps)
    n = shape_vec(res['r'])
    cls = tb.classify(res)
    reps = 1 + ens

    # within: per bundle, then averaged over the quad's samples
    w_shape, w_out = [], []
    for g in range(N * N):
        sel = np.nonzero(gid == g)[0]
        w_shape.append(svar(n[sel]))
        w_out.append(disagree(cls[sel]))
    within = float(np.mean(w_shape))
    within_out = float(np.mean(w_out))

    # between: across the nominal sample of each bundle (copy 0)
    nom = np.arange(0, len(gid), reps)
    between = svar(n[nom])
    between_out = disagree(cls[nom])

    ratio = between / (between + within) if (between + within) > 1e-12 else float('nan')
    ratio_out = (between_out / (between_out + within_out)
                 if (between_out + within_out) > 1e-12 else float('nan'))
    return dict(within=within, between=between, ratio=ratio,
                within_out=within_out, between_out=between_out, ratio_out=ratio_out,
                spacing=hx, drift=float(np.median(res['drift'])))


def refine_series(name, cx, cy, half0, levels=4, **kw):
    """Follow one quad down `levels` of refinement, halving the half-width each time."""
    print(f"\n=== {name}  centre=({cx},{cy}) ===")
    print(f"{'lvl':>3} {'half':>8} {'spacing':>9} | {'within':>8} {'between':>8} {'ratio':>7} "
          f"| {'w_out':>6} {'b_out':>6} {'r_out':>6} | {'drift':>8}")
    rows = []
    h = half0
    for L in range(levels):
        r = probe(cx, cy, h, **kw)
        rows.append((L, h, r))
        print(f"{L:>3} {h:>8.4f} {r['spacing']:>9.2e} | {r['within']:>8.4f} {r['between']:>8.4f} "
              f"{r['ratio']:>7.3f} | {r['within_out']:>6.3f} {r['between_out']:>6.3f} "
              f"{r['ratio_out']:>6.3f} | {r['drift']:>8.1e}")
        h /= 2.0
    return rows
