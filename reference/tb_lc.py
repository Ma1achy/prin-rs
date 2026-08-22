"""
Levi-Civita regularisation of the planar three-body problem.

The close pair (i,j) is regularised; the third body k is carried as a perturber.

COORDINATES
    rho = r_j - r_i                     (close-pair separation)
    lam = r_k - COM(i,j)                (third body relative to the pair barycentre)
    mu_rho = m_i m_j / (m_i + m_j)      mu_lam = m_k (m_i + m_j) / M

LEVI-CIVITA (planar, so the KS map is just complex squaring)
    rho = L(u) u,   L(u) = [[u1, -u2], [u2, u1]],   i.e. rho = (u1^2 - u2^2, 2 u1 u2)
    |rho| = |u|^2                       d rho = 2 L(u) du      p_u = 2 L(u)^T p_rho

TIME TRANSFORMATION
    dt = |rho| dtau = |u|^2 dtau

REGULARISED HAMILTONIAN   Gamma = |rho| (H - E)
    Gamma = |p_u|^2/(8 mu_rho) + |u|^2 [ |p_lam|^2/(2 mu_lam) + V_pert - E ] - G m_i m_j

The singular term -G m_i m_j / |rho| has become a CONSTANT. Gamma and its derivatives are
finite at u = 0, so the pair may pass through exact collision without the integrator failing.

Integrated with RK4 in fictitious time tau. Gamma is not separable (the |u|^2 |p_lam|^2 term
couples u and p_lam), so leapfrog does not apply; RK4 is used and energy error is monitored.
"""
import numpy as np
import tb


def _Lmat_apply(u, w):
    """L(u) @ w  for batches. u,w: (S,2) -> (S,2)"""
    return np.stack([u[:, 0]*w[:, 0] - u[:, 1]*w[:, 1],
                     u[:, 1]*w[:, 0] + u[:, 0]*w[:, 1]], axis=1)


def _LmatT_apply(u, w):
    """L(u)^T @ w"""
    return np.stack([u[:, 0]*w[:, 0] + u[:, 1]*w[:, 1],
                     -u[:, 1]*w[:, 0] + u[:, 0]*w[:, 1]], axis=1)


def rho_of_u(u):
    return np.stack([u[:, 0]**2 - u[:, 1]**2, 2*u[:, 0]*u[:, 1]], axis=1)


def u_of_rho(rho):
    """Inverse LC map (branch choice is irrelevant: rho(u) = rho(-u))."""
    r = np.sqrt(np.einsum('sk,sk->s', rho, rho))
    u1 = np.sqrt(np.maximum(0.5*(r + rho[:, 0]), 0.0))
    u2 = np.where(u1 > 1e-300, rho[:, 1]/(2*np.maximum(u1, 1e-300)),
                  np.sqrt(np.maximum(0.5*(r - rho[:, 0]), 0.0)))
    return np.stack([u1, u2], axis=1)


class LCSystem:
    """Regularises pair (i,j); k is the perturber."""

    def __init__(self, i, j, k, M=None):
        M = tb.M if M is None else M
        self.i, self.j, self.k = i, j, k
        self.mi, self.mj, self.mk = M[i], M[j], M[k]
        self.mij = self.mi + self.mj
        self.mu_rho = self.mi*self.mj/self.mij
        self.mu_lam = max(self.mk*self.mij/(self.mij + self.mk), 1e-300)  # guard m_k -> 0
        self.bi = self.mi/self.mij
        self.bj = self.mj/self.mij
        self.Gmimj = self.mi*self.mj

    # ---- perturbation from the third body ----
    def _sep(self, rho, lam):
        d1 = lam + self.bj*rho          # r_k - r_i
        d2 = lam - self.bi*rho          # r_k - r_j
        return d1, d2

    def V_pert(self, rho, lam):
        d1, d2 = self._sep(rho, lam)
        n1 = np.sqrt(np.einsum('sk,sk->s', d1, d1))
        n2 = np.sqrt(np.einsum('sk,sk->s', d2, d2))
        return -(self.mk*self.mi/np.maximum(n1, 1e-300)
                 + self.mk*self.mj/np.maximum(n2, 1e-300))

    def dV(self, rho, lam):
        """(dV/drho, dV/dlam)"""
        d1, d2 = self._sep(rho, lam)
        n1 = np.maximum(np.sqrt(np.einsum('sk,sk->s', d1, d1)), 1e-300)[:, None]
        n2 = np.maximum(np.sqrt(np.einsum('sk,sk->s', d2, d2)), 1e-300)[:, None]
        g1 = self.mk*self.mi*d1/n1**3
        g2 = self.mk*self.mj*d2/n2**3
        return self.bj*g1 - self.bi*g2, g1 + g2

    # ---- energy in physical variables ----
    def energy(self, rho, p_rho, lam, p_lam):
        return (np.einsum('sk,sk->s', p_rho, p_rho)/(2*self.mu_rho)
                + np.einsum('sk,sk->s', p_lam, p_lam)/(2*self.mu_lam)
                - self.Gmimj/np.maximum(np.sqrt(np.einsum('sk,sk->s', rho, rho)), 1e-300)
                + self.V_pert(rho, lam))

    # ---- regularised equations of motion in tau ----
    def deriv(self, s, E):
        u, pu, lam, plam = s[:, 0:2], s[:, 2:4], s[:, 4:6], s[:, 6:8]
        r2 = np.einsum('sk,sk->s', u, u)                 # |u|^2 = |rho|
        rho = rho_of_u(u)
        dVrho, dVlam = self.dV(rho, lam)
        Tlam = np.einsum('sk,sk->s', plam, plam)/(2*self.mu_lam)
        brack = Tlam + self.V_pert(rho, lam) - E         # the [.] in Gamma

        du = pu/(4*self.mu_rho)
        dVu = 2*_LmatT_apply(u, dVrho)                   # dV/du = 2 L^T dV/drho
        dpu = -2*u*brack[:, None] - r2[:, None]*dVu
        dlam = r2[:, None]*plam/self.mu_lam
        dplam = -r2[:, None]*dVlam
        dt = r2
        return np.concatenate([du, dpu, dlam, dplam, dt[:, None]], axis=1)

    # ---- conversions ----
    def to_reg(self, r, v):
        i, j, k = self.i, self.j, self.k
        rho = r[:, j, :] - r[:, i, :]
        vrho = v[:, j, :] - v[:, i, :]
        com = (self.mi*r[:, i, :] + self.mj*r[:, j, :])/self.mij
        vcom = (self.mi*v[:, i, :] + self.mj*v[:, j, :])/self.mij
        lam = r[:, k, :] - com
        vlam = v[:, k, :] - vcom
        p_rho = self.mu_rho*vrho
        p_lam = self.mu_lam*vlam
        u = u_of_rho(rho)
        pu = 2*_LmatT_apply(u, p_rho)
        E = self.energy(rho, p_rho, lam, p_lam)
        s = np.concatenate([u, pu, lam, p_lam, np.zeros((len(u), 1))], axis=1)
        return s, E

    def to_phys(self, s):
        u, pu, lam, plam = s[:, 0:2], s[:, 2:4], s[:, 4:6], s[:, 6:8]
        r2 = np.maximum(np.einsum('sk,sk->s', u, u), 1e-300)
        rho = rho_of_u(u)
        p_rho = _Lmat_apply(u, pu)/(2*r2[:, None])
        return rho, p_rho, lam, plam

    def to_cartesian(self, s):
        """Back to (S,3,2) positions/velocities with the pair barycentre at the origin."""
        rho, p_rho, lam, plam = self.to_phys(s)
        n = len(rho)
        r = np.zeros((n, 3, 2)); v = np.zeros((n, 3, 2))
        vrho = p_rho/self.mu_rho; vlam = plam/self.mu_lam
        r[:, self.i, :] = -self.bj*rho
        r[:, self.j, :] = self.bi*rho
        r[:, self.k, :] = lam
        v[:, self.i, :] = -self.bj*vrho
        v[:, self.j, :] = self.bi*vrho
        v[:, self.k, :] = vlam
        return r, v


def integrate_lc(r0, v0, t_max, pair=(1, 2), dtau=None, n_steps=200000, M=None, eta=0.01):
    """RK4 in fictitious time until every trajectory passes t_max."""
    i, j = pair
    k = [q for q in range(3) if q not in pair][0]
    sysm = LCSystem(i, j, k, M=M)
    s, E = sysm.to_reg(r0, v0)
    n = len(s)

    if dtau is None:
        # Uniform dtau is the point of LC: the regularised pair is a harmonic oscillator, so
        # constant steps in tau resolve the collision. Set it from the INITIAL separation:
        #   dt = |rho| dtau  and  dt ~ eta sqrt(|rho|^3/Gm)  =>  dtau ~ eta sqrt(|rho|/Gm)
        # which scales as sqrt(alpha) under r -> alpha r, i.e. scale-covariant.
        rho0 = rho_of_u(s[:, 0:2])
        r0n = np.sqrt(np.einsum('sk,sk->s', rho0, rho0))
        dtau = float(eta*np.median(np.sqrt(np.maximum(r0n, 1e-30)/sysm.mij)))

    dmin = np.full(n, np.inf)
    for _ in range(n_steps):
        if (s[:, 8] >= t_max).all():
            break
        act = (s[:, 8] < t_max)[:, None].astype(float)
        h = dtau*act
        k1 = sysm.deriv(s, E)
        k2 = sysm.deriv(s + 0.5*h*k1, E)
        k3 = sysm.deriv(s + 0.5*h*k2, E)
        k4 = sysm.deriv(s + h*k3, E)
        s = s + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        rho = rho_of_u(s[:, 0:2])
        dmin = np.minimum(dmin, np.sqrt(np.einsum('sk,sk->s', rho, rho)))

    rho, p_rho, lam, plam = sysm.to_phys(s)
    E1 = sysm.energy(rho, p_rho, lam, plam)
    drift = np.abs((E1 - E)/np.maximum(np.abs(E), 1e-30))
    r, v = sysm.to_cartesian(s)
    return dict(r=r, v=v, drift=drift, dmin=dmin, t=s[:, 8], sys=sysm, dtau=dtau)
