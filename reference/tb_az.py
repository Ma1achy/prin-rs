"""
Aarseth-Zare GLOBAL regularisation (planar three-body).

LC regularises ONE pair, so a close approach of either other pair still enters through an
unregularised term -- measured: near-field drift 5.0e-04 unregularised -> 1.7e+01 with LC on the
wrong pair, and switching still fails when two pairs are close at once (findings doc §7.7).

AZ regularises TWO pairs simultaneously. Pick the reference body `a` to be the one NOT in the
longest side; then both regularised pairs share `a`:

    R1 = r_b - r_a        R2 = r_c - r_a        R3 = R2 - R1 = r_c - r_b   (unregularised)

Because (b,c) is the LONGEST side, |R3| >= max(|R1|,|R2|), so R3 -> 0 only when all three
separations -> 0, i.e. genuine TRIPLE collision -- which is provably non-regularisable anyway.
Everything else is covered.

KINETIC ENERGY.  R1,R2 are relative-to-a-common-body coordinates, not Jacobi, so the kinetic
energy carries a cross term:

    T = P1^2/(2 mu1) + P1.P2/m_a + P2^2/(2 mu2),    1/mu1 = 1/m_a + 1/m_b,  1/mu2 = 1/m_a + 1/m_c

(verified by inverting the velocity-space mass matrix, whose determinant is m_a m_b m_c / M).

LEVI-CIVITA on both:  R1 = L(u1)u1,  R2 = L(u2)u2,  A = |u1|^2 = |R1|,  B = |u2|^2 = |R2|
TIME:  dt = |R1||R2| dtau = A B dtau      -- vanishes when EITHER pair closes

REGULARISED HAMILTONIAN  Gamma = A B (H - E):

    Gamma = B|p1|^2/(8 mu1) + A|p2|^2/(8 mu2) + (L1 p1).(L2 p2)/(4 m_a)
            - m_a m_b B - m_a m_c A - A B m_b m_c/|R3| - E A B

Both  -m_a m_b/|R1|  and  -m_a m_c/|R2|  have become CONSTANT-order terms (linear in B and A).
Nothing is singular unless |R3| -> 0.
"""
import numpy as np
import tb
from tb_lc import _Lmat_apply, _LmatT_apply, rho_of_u, u_of_rho

TRIPLES = [(0, 1, 2), (1, 0, 2), (2, 0, 1)]   # (reference a, b, c) for each choice of a
PAIRS = [(0, 1), (0, 2), (1, 2)]


class AZSystem:
    """Reference body a; regularised pairs (a,b) and (a,c); unregularised side (b,c)."""

    def __init__(self, a, b, c, M=None):
        M = tb.M if M is None else M
        self.a, self.b, self.c = a, b, c
        self.ma, self.mb, self.mc = M[a], M[b], M[c]
        self.M = M
        self.Mtot = M.sum()
        self.mu1 = self.ma*self.mb/(self.ma + self.mb)
        self.mu2 = self.ma*self.mc/(self.ma + self.mc)
        # velocity-space mass matrix for (Rdot1, Rdot2)
        Mt = self.Mtot
        self.K11 = self.mb*(self.ma + self.mc)/Mt
        self.K22 = self.mc*(self.ma + self.mb)/Mt
        self.K12 = -self.mb*self.mc/Mt

    # ---------- conversions ----------
    def to_reg(self, r, v):
        a, b, c = self.a, self.b, self.c
        R1 = r[:, b, :] - r[:, a, :]
        R2 = r[:, c, :] - r[:, a, :]
        V1 = v[:, b, :] - v[:, a, :]
        V2 = v[:, c, :] - v[:, a, :]
        P1 = self.K11*V1 + self.K12*V2
        P2 = self.K12*V1 + self.K22*V2
        u1 = u_of_rho(R1); u2 = u_of_rho(R2)
        p1 = 2*_LmatT_apply(u1, P1)
        p2 = 2*_LmatT_apply(u2, P2)
        E = self.energy_phys(R1, R2, P1, P2)
        s = np.concatenate([u1, p1, u2, p2, np.zeros((len(u1), 1))], axis=1)
        return s, E

    def phys_from_state(self, s):
        u1, p1, u2, p2 = s[:, 0:2], s[:, 2:4], s[:, 4:6], s[:, 6:8]
        A = np.maximum(np.einsum('sk,sk->s', u1, u1), 1e-300)
        B = np.maximum(np.einsum('sk,sk->s', u2, u2), 1e-300)
        R1 = rho_of_u(u1); R2 = rho_of_u(u2)
        P1 = _Lmat_apply(u1, p1)/(2*A[:, None])
        P2 = _Lmat_apply(u2, p2)/(2*B[:, None])
        return R1, R2, P1, P2

    def energy_phys(self, R1, R2, P1, P2):
        r1 = np.maximum(np.sqrt(np.einsum('sk,sk->s', R1, R1)), 1e-300)
        r2 = np.maximum(np.sqrt(np.einsum('sk,sk->s', R2, R2)), 1e-300)
        R3 = R2 - R1
        r3 = np.maximum(np.sqrt(np.einsum('sk,sk->s', R3, R3)), 1e-300)
        T = (np.einsum('sk,sk->s', P1, P1)/(2*self.mu1)
             + np.einsum('sk,sk->s', P1, P2)/self.ma
             + np.einsum('sk,sk->s', P2, P2)/(2*self.mu2))
        U = self.ma*self.mb/r1 + self.ma*self.mc/r2 + self.mb*self.mc/r3
        return T - U

    def to_cartesian(self, s):
        """Back to (S,3,2) in the COM frame."""
        R1, R2, P1, P2 = self.phys_from_state(s)
        det = self.K11*self.K22 - self.K12**2
        V1 = (self.K22*P1 - self.K12*P2)/det
        V2 = (self.K11*P2 - self.K12*P1)/det
        n = len(R1)
        r = np.zeros((n, 3, 2)); v = np.zeros((n, 3, 2))
        r[:, self.a, :] = 0.0;  r[:, self.b, :] = R1;  r[:, self.c, :] = R2
        v[:, self.a, :] = 0.0;  v[:, self.b, :] = V1;  v[:, self.c, :] = V2
        Mv = self.M[None, :, None]
        r -= (Mv*r).sum(axis=1, keepdims=True)/self.Mtot
        v -= (Mv*v).sum(axis=1, keepdims=True)/self.Mtot
        return r, v

    # ---------- regularised equations of motion ----------
    def deriv(self, s, E):
        u1, p1, u2, p2 = s[:, 0:2], s[:, 2:4], s[:, 4:6], s[:, 6:8]
        A = np.einsum('sk,sk->s', u1, u1)
        B = np.einsum('sk,sk->s', u2, u2)
        R1 = rho_of_u(u1); R2 = rho_of_u(u2)
        R3 = R2 - R1
        r3 = np.maximum(np.sqrt(np.einsum('sk,sk->s', R3, R3)), 1e-300)
        L1p1 = _Lmat_apply(u1, p1)
        L2p2 = _Lmat_apply(u2, p2)
        n1 = np.einsum('sk,sk->s', p1, p1)
        n2 = np.einsum('sk,sk->s', p2, p2)
        mbc = self.mb*self.mc

        du1 = B[:, None]*p1/(4*self.mu1) + _LmatT_apply(u1, L2p2)/(4*self.ma)
        du2 = A[:, None]*p2/(4*self.mu2) + _LmatT_apply(u2, L1p1)/(4*self.ma)

        # dGamma/du1
        g1 = (u1*(n2/(4*self.mu2))[:, None]
              + _LmatT_apply(p1, L2p2)/(4*self.ma)
              - 2*self.ma*self.mc*u1
              - mbc*(2*u1*(B/r3)[:, None] + 2*(A*B/r3**3)[:, None]*_LmatT_apply(u1, R3))
              - 2*(E*B)[:, None]*u1)
        # dGamma/du2
        g2 = (u2*(n1/(4*self.mu1))[:, None]
              + _LmatT_apply(p2, L1p1)/(4*self.ma)
              - 2*self.ma*self.mb*u2
              - mbc*(2*u2*(A/r3)[:, None] - 2*(A*B/r3**3)[:, None]*_LmatT_apply(u2, R3))
              - 2*(E*A)[:, None]*u2)

        dt = A*B
        return np.concatenate([du1, -g1, du2, -g2, dt[:, None]], axis=1)


def choose_reference(r):
    """Reference body = the one NOT in the longest side, so the unregularised side is longest."""
    pd = tb.pair_dists(r)                 # order PAIRS = [(0,1),(0,2),(1,2)]
    longest = np.argmax(pd, axis=1)
    ref = np.array([2, 1, 0])[longest]    # body not in that pair
    return ref


# How `dtau` is sized within a sync interval. Ported from the corrected Rust
# (`src/integrate/az/driver.rs::DtauMode`), which is the considered implementation.
#
#   'fixed'              -- once at interval entry, never updated. What this file did until now,
#                           and therefore what EVERY committed NumPy number in the corpus used.
#                           dt = A*B*dtau, so the physical step is eta*dt_left only while A*B
#                           stays near its entry value; a trajectory at a close encounter AT a
#                           sync boundary has a tiny A0*B0, so dtau is enormous, and as the
#                           bodies separate through the interval dt grows with A*B. Giant
#                           physical steps immediately after an encounter.
#   'per-step-remaining' -- eta*(dt_left - s.t)/(A*B) per step. Zeno by arithmetic:
#                           dt ~ eta*rem gives rem_{n+1} = rem_n (1-eta), so the interval is
#                           approached geometrically and never completed. Kept because it looks
#                           right, not because it is a candidate.
#   'per-step-interval'  -- min(eta*dt_left/(A*B), dtau_entry). dt ~ eta*dt_left throughout, so
#                           the step count stays at ~1/eta. The cap is ONE-SIDED in the right
#                           direction: when A*B grows the recomputed value falls and the blow-up
#                           goes; when A*B falls at a close approach the cap holds dtau at
#                           nominal so dt shrinks with the separation -- which is what
#                           regularisation buys and what the comment below always wanted.
DTAU_MODES = ('fixed', 'per-step-remaining', 'per-step-interval')
DTAU_MODE_DEFAULT = 'per-step-interval'


def _dtau_for_step(mode, eta, dt_left, dtau_entry, s):
    """(A*B) recomputed from the CURRENT state. Mirrors `dtau_for_step` in driver.rs."""
    A = np.maximum(np.einsum('sk,sk->s', s[:, 0:2], s[:, 0:2]), 1e-300)
    B = np.maximum(np.einsum('sk,sk->s', s[:, 4:6], s[:, 4:6]), 1e-300)
    if mode == 'fixed':
        return dtau_entry
    if mode == 'per-step-remaining':
        return eta*np.maximum(dt_left - s[:, 8], 0.0)/(A*B)
    if mode == 'per-step-interval':
        return np.minimum(eta*dt_left/(A*B), dtau_entry)
    raise ValueError(f"dtau_mode: expected one of {DTAU_MODES}, got {mode!r}")


def integrate_az(r0, v0, t_max, n_sync=32, eta=0.02, max_steps=8000, M=None,
                 dtau_mode=DTAU_MODE_DEFAULT):
    Mv = tb.M if M is None else M
    n = r0.shape[0]
    r = r0.copy(); v = v0.copy()
    t = np.zeros(n)
    dmin = np.full(n, np.inf)
    switches = np.zeros(n, dtype=np.int64)
    prev_ref = np.full(n, -1)
    E0 = tb.energy(r, v, 0.0)

    for kk in range(n_sync):
        t_target = (kk + 1)*t_max/n_sync
        ref = choose_reference(r)
        switches += (prev_ref >= 0) & (ref != prev_ref)
        prev_ref = ref.copy()

        for a, b, c in TRIPLES:
            sel = np.nonzero((ref == a) & (t < t_target - 1e-15))[0]
            if len(sel) == 0:
                continue
            sysm = AZSystem(a, b, c, M=Mv)
            s, E = sysm.to_reg(r[sel], v[sel])
            dt_left = t_target - t[sel]

            # The interval's ENTRY sizing. Under 'fixed' it is the whole step control; under
            # 'per-step-interval' it is the one-sided cap that keeps the physical step shrinking
            # at a close approach instead of being pinned at eta*dt_left.
            A0 = np.maximum(np.einsum('sk,sk->s', s[:, 0:2], s[:, 0:2]), 1e-300)
            B0 = np.maximum(np.einsum('sk,sk->s', s[:, 4:6], s[:, 4:6]), 1e-300)
            dtau = eta*dt_left/(A0*B0)

            for _ in range(max_steps):
                # NaN >= x is False, so a non-finite trajectory never satisfies `done` and the
                # loop would burn the entire max_steps budget (measured: 354 s vs ~3 s nominal).
                # Treat non-finite as done: it is excluded by the drift gate downstream anyway.
                bad = ~np.isfinite(s).all(axis=1)
                done = (s[:, 8] >= dt_left) | bad
                if done.all():
                    break
                h = (_dtau_for_step(dtau_mode, eta, dt_left, dtau, s)*(~done)).astype(float)[:, None]
                k1 = sysm.deriv(s, E)
                k2 = sysm.deriv(s + 0.5*h*k1, E)
                k3 = sysm.deriv(s + 0.5*h*k2, E)
                k4 = sysm.deriv(s + h*k3, E)
                s = s + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)
                R1, R2, _, _ = sysm.phys_from_state(s)
                d1 = np.sqrt(np.einsum('sk,sk->s', R1, R1))
                d2 = np.sqrt(np.einsum('sk,sk->s', R2, R2))
                dmin[sel] = np.minimum(dmin[sel], np.minimum(d1, d2))

            rc, vc = sysm.to_cartesian(s)
            r[sel] = rc; v[sel] = vc
            t[sel] = t[sel] + np.minimum(s[:, 8], dt_left)

    drift = np.abs((tb.energy(r, v, 0.0) - E0)/np.maximum(np.abs(E0), 1e-30))
    return dict(r=r, v=v, drift=drift, dmin=dmin, t=t, switches=switches)
