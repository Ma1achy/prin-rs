//! The shared decoder `D` and canonicaliser `C`: **every chart ends here.**
//!
//! ```text
//! (u,v) --Phi--> chart space --D--> (m, r, p) --C--> canonical IC --> integrator
//! ```
//!
//! Adding a chart means adding a `Phi`, nothing else. The integrator never learns which chart
//! produced its input.
//!
//! **Transcribed from the chart reference §0, not re-derived.** The algebra here is the kind
//! that fails silently — the crossed mass factors in [`momenta`] are the named hazard — and the
//! project's standing rule is to port rather than derive.
//!
//! # A correction to the reference, measured
//!
//! The reference says of the crossed mass factors: *"Transcribe it, then verify by asserting
//! `sum p_i = 0` to machine precision — the test that catches a swap."* **It does not catch it.**
//! Both forms sum to zero identically:
//!
//! ```text
//! crossed:    sum p = p_lam * (1 - (m0+m1)/M01) = 0
//! uncrossed:  sum p = p_lam * (1 - (m1+m0)/M01) = 0
//! ```
//!
//! Measured at Burrau's masses: `|sum p|` is `7.9e-17` crossed and `5.6e-17` uncrossed. The
//! assertion is decoration — a test that cannot fire.
//!
//! Two things do catch it, both hard:
//!
//! - **the Jacobi round-trip** `p_rho == (m0*p1 - m1*p0)/M01`: `1.1e-16` crossed against
//!   `6.8e-2` uncrossed;
//! - **the kinetic-energy identity** `K == |p_rho|^2/(2 mu_rho) + |p_lam|^2/(2 mu_lam)`:
//!   `4.4e-16` crossed against `2.6e-1` uncrossed.
//!
//! Both require `m0 != m1` — at equal masses the two forms coincide and neither can fire. Burrau
//! has `m0 = 3, m1 = 4`, so they do. `tests/charts.rs` runs them with the swap as a negative
//! control, and states the equal-mass exclusion.
//!
//! `sum p_i = 0` is still asserted, because it catches a different family of errors (a dropped
//! term, a sign flip on `p_lam`). It is just not the one the reference says it is.
//!
//! # The GLSL reference is the pin, and it carries ten slots for eight coordinates
//!
//! `Ma1achy/principia-ii`, `src/shaders/principia/frag.glsl:19-59`, is the validated
//! implementation and settles the constants the LaTeX chart reference could only guess at. Its
//! `decodeIC` takes `z0..z9` but **never reads `z2` or `z3`** — dead slots from before the chart
//! was known to be 8D. They are dropped here, and the two angle coordinates are stored in the
//! *spec's* order rather than the GLSL's:
//!
//! ```text
//!   GLSL      this module        meaning
//!   z0     -> index 1  z_beta    beta  = PI * sigmoid(z)
//!   z1     -> index 0  z_alpha   alpha = ALPHA_MIN + (PI/2 - 2*ALPHA_MIN)*sigmoid(z)
//!   z2, z3 -> --                 DROPPED, never read by decodeIC
//!   z4, z5 -> index 2, 3         p_rho    = Q_MAX*(2*sigmoid(z) - 1), per component
//!   z6, z7 -> index 4, 5         p_lambda = same
//!   z8, z9 -> index 6, 7         mu1, mu2 = MU_MAX*(2*sigmoid(z) - 1)
//! ```
//!
//! **The alpha/beta order flips.** The GLSL puts beta at index 0; the spec names the chart
//! `(z_alpha, z_beta)` and that is the order used here. A consequence is that the GLSL's `shape`
//! preset — `q1 = e0, q2 = e1`, which in *its* indexing is `beta x alpha` — becomes `alpha x
//! beta` here, so the rendered image is **transposed relative to the GLSL**. That is the port
//! being faithful to the spec, not a bug.

use crate::physics::{shape, Cart, Ic};
use crate::Vec2;

/// Logit saturation for the mass coordinates. Pinned from `frag.glsl:21`; the LaTeX chart
/// reference's guessed `4.0` was wrong.
pub const MU_MAX: f64 = 5.0;
/// Saturation for the free Jacobi momentum coordinates. Pinned from `frag.glsl:22`.
pub const Q_MAX: f64 = 2.0;
/// Buffer keeping `‖rho‖` away from zero. Note the orientation: `‖rho~‖ = cos(alpha)`, so
/// **small alpha is a LARGE inner-pair separation** and `alpha -> pi/2` is a tight inner pair
/// with a distant third body. That is the reference's own "easy to get backwards" note, and
/// `tests/charts.rs` asserts the direction rather than the formula.
pub const ALPHA_MIN: f64 = 0.05;
/// Mirror deadband in the canonicaliser.
pub const DELTA_LAM: f64 = 1e-12;
/// `M01` below this is degenerate.
pub const M01_EPS: f64 = 1e-12;
/// Minimum mass-weighted norm for a momentum seed direction to be usable.
pub const SEED_EPS: f64 = 1e-10;

/// Why a pixel could not be decoded. **No pixel is ever rejected** — the label travels with it,
/// the same way a non-finite copy is a measurement outcome rather than missing data.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Degenerate {
    /// `m0 + m1` underflowed.
    M01Tiny,
    /// The configuration has zero moment of inertia.
    InertiaZero,
    /// A target kinetic energy below the rigid-rotation minimum `Lz^2/(2I)`.
    BelowKMin,
    /// All four momentum seed directions degenerated.
    SeedsExhausted,
    /// Energy normalisation asked for `E* < U`.
    EnergyInfeasible,
}

impl Degenerate {
    pub fn name(self) -> &'static str {
        match self {
            Degenerate::M01Tiny => "M01_TINY",
            Degenerate::InertiaZero => "INERTIA_ZERO",
            Degenerate::BelowKMin => "BELOW_K_MIN",
            Degenerate::SeedsExhausted => "SEEDS_EXHAUSTED",
            Degenerate::EnergyInfeasible => "ENERGY_INFEASIBLE",
        }
    }
}

/// A decode result: always an initial condition, sometimes with a label attached.
#[derive(Clone, Copy, Debug)]
pub struct Decoded {
    pub ic: Ic<f64>,
    pub flag: Option<Degenerate>,
}

fn sigmoid(z: f64) -> f64 {
    1.0 / (1.0 + (-z).exp())
}

/// The 8D latent chart's coordinate. **No coordinate is spent on gauge** — rotation and scale
/// are fixed by the canonical-frame decode, which is why this is 8D and not 10D.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Latent {
    pub z_alpha: f64,
    pub z_beta: f64,
    pub z_q: [f64; 4],
    pub z_mu: [f64; 2],
}

impl Latent {
    /// Index into the 8 coordinates, in the reference's order:
    /// `(z_alpha, z_beta | z_q0..z_q3 | z_mu1, z_mu2)`.
    ///
    /// **Out of range panics rather than aliasing.** These previously used an irrefutable `_`
    /// arm, so index 8 or 9 silently read `z_mu[1]` — and the GLSL's ten slots collapsing to
    /// eight is exactly the situation in which an out-of-range index gets written by hand.
    pub fn get(&self, i: usize) -> f64 {
        match i {
            0 => self.z_alpha,
            1 => self.z_beta,
            2..=5 => self.z_q[i - 2],
            6 => self.z_mu[0],
            7 => self.z_mu[1],
            _ => panic!("latent coordinate {i} out of range (8 coordinates, 0..=7)"),
        }
    }
    pub fn set(&mut self, i: usize, v: f64) {
        match i {
            0 => self.z_alpha = v,
            1 => self.z_beta = v,
            2..=5 => self.z_q[i - 2] = v,
            6 => self.z_mu[0] = v,
            7 => self.z_mu[1] = v,
            _ => panic!("latent coordinate {i} out of range (8 coordinates, 0..=7)"),
        }
    }
    pub const AXIS_NAMES: [&'static str; 8] =
        ["z_alpha", "z_beta", "z_q0", "z_q1", "z_q2", "z_q3", "z_mu1", "z_mu2"];
}

// ---------------------------------------------------------------------------------------------
// 0.1 Masses
// ---------------------------------------------------------------------------------------------

/// `mu_k = MU_MAX*(2*sigmoid(z_k) - 1)`, then `softmax(0, mu1, mu2)`. Normalised so `M = 1`.
///
/// **This is `MU_MAX*tanh(z/2)` — HALF the gain of the LaTeX reference's `mu_max*tanh(z)`.**
/// `frag.glsl:35-36` is the pin. The half-gain form is shared with the momentum coordinates,
/// which is why both read `2*sigmoid(z) - 1` and neither reads `tanh`.
pub fn masses(z_mu: [f64; 2]) -> ([f64; 3], Option<Degenerate>) {
    let mu1 = MU_MAX * (2.0 * sigmoid(z_mu[0]) - 1.0);
    let mu2 = MU_MAX * (2.0 * sigmoid(z_mu[1]) - 1.0);
    // Softmax with the maximum subtracted, the same conditioning the vMF weights use.
    let mx = 0f64.max(mu1).max(mu2);
    let e = [(-mx).exp(), (mu1 - mx).exp(), (mu2 - mx).exp()];
    let z: f64 = e.iter().sum();
    let m = [e[0] / z, e[1] / z, e[2] / z];
    let flag = if m[0] + m[1] < M01_EPS { Some(Degenerate::M01Tiny) } else { None };
    (m, flag)
}

// ---------------------------------------------------------------------------------------------
// 0.2 Configuration — hyperspherical mass-weighted Jacobi, canonical frame
// ---------------------------------------------------------------------------------------------

/// `(mu_rho, mu_lam)` for the `(0,1)` inner pair with body 2 outer. Total mass need not be 1.
pub fn reduced(m: &[f64; 3]) -> (f64, f64) {
    let m01 = m[0] + m[1];
    let mtot = m01 + m[2];
    (m[0] * m[1] / m01, m[2] * m01 / mtot)
}

/// `(alpha, beta)` from the two configuration coordinates.
///
/// `alpha in [ALPHA_MIN, pi/2 - ALPHA_MIN]`, `beta in [0, pi]` — the half-range on `beta` is
/// what fixes the mirror gauge, so the canonicaliser is a no-op away from the seam.
pub fn angles(z_alpha: f64, z_beta: f64) -> (f64, f64) {
    let hi = std::f64::consts::FRAC_PI_2 - 2.0 * ALPHA_MIN;
    (ALPHA_MIN + hi * sigmoid(z_alpha), std::f64::consts::PI * sigmoid(z_beta))
}

/// Positions from `(alpha, beta)` at unit hyperradius, COM at the origin by construction.
///
/// `‖rho~‖ = cos(alpha)` — **small alpha is a wide inner pair**, not a tight one.
pub fn config(alpha: f64, beta: f64, m: &[f64; 3]) -> [Vec2<f64>; 3] {
    let (mu_rho, mu_lam) = reduced(m);
    let rho_t = Vec2::new(alpha.cos(), 0.0);
    let lam_t = Vec2::new(alpha.sin() * beta.cos(), alpha.sin() * beta.sin());

    let rho = rho_t / mu_rho.sqrt();
    let lam = lam_t / mu_lam.sqrt();

    let m01 = m[0] + m[1];
    let mtot = m01 + m[2];
    let r01 = lam * (-m[2] / mtot);
    [r01 - rho * (m[1] / m01), r01 + rho * (m[0] / m01), r01 + lam]
}

// ---------------------------------------------------------------------------------------------
// 0.3 Momentum — free Jacobi momenta
// ---------------------------------------------------------------------------------------------

/// Particle momenta from the four free Jacobi momentum coordinates.
///
/// **The transcription hazard.** The `m0` and `m1` factors are *crossed* relative to the
/// position reconstruction: positions take `-m1/M01` on `r0`, momenta take `-m0/M01` on `p0`.
/// See the module note — `sum p = 0` does **not** catch a swap here; the Jacobi round-trip and
/// the kinetic-energy identity do.
pub fn momenta(z_q: [f64; 4], m: &[f64; 3]) -> [Vec2<f64>; 3] {
    let q: Vec<f64> = z_q.iter().map(|&z| Q_MAX * (2.0 * sigmoid(z) - 1.0)).collect();
    let p_rho = Vec2::new(q[0], q[1]);
    let p_lam = Vec2::new(q[2], q[3]);
    from_jacobi_momenta(p_rho, p_lam, m)
}

/// `(p_rho, p_lam)` to particle momenta. Split out so the crossed factors have exactly one home.
pub fn from_jacobi_momenta(p_rho: Vec2<f64>, p_lam: Vec2<f64>, m: &[f64; 3]) -> [Vec2<f64>; 3] {
    let m01 = m[0] + m[1];
    [
        -p_rho - p_lam * (m[0] / m01),
        p_rho - p_lam * (m[1] / m01),
        p_lam,
    ]
}

/// Particle momenta back to `(p_rho, p_lam)`. The inverse of [`from_jacobi_momenta`].
///
/// This is the test the reference asked `sum p = 0` to do. Round-tripping through it separates
/// the crossed form from the swapped one by `6.8e-2` at Burrau's masses, where `sum p` separates
/// them by nothing at all.
pub fn to_jacobi_momenta(p: &[Vec2<f64>; 3], m: &[f64; 3]) -> (Vec2<f64>, Vec2<f64>) {
    let m01 = m[0] + m[1];
    ((p[1] * m[0] - p[0] * m[1]) / m01, p[2])
}

// ---------------------------------------------------------------------------------------------
// 0.4 / 0.5 Canonicalisation and scale gauge
// ---------------------------------------------------------------------------------------------

/// Rotate so `rho` lies along `+x`, then mirror if `lam_y < -DELTA_LAM`.
///
/// A no-op away from the seam under the canonical-frame decode, which is why it is implemented
/// anyway: charts that bypass that decode — the Burrau family, the mass simplex — do need it,
/// and having it in one place is what keeps them comparable to the latent chart.
pub fn canonicalise(r: &mut [Vec2<f64>; 3], p: &mut [Vec2<f64>; 3], m: &[f64; 3]) {
    let m01 = m[0] + m[1];
    let rho = r[1] - r[0];
    let phi = rho.y.atan2(rho.x);
    let (c, s) = ((-phi).cos(), (-phi).sin());
    let rot = |v: Vec2<f64>| Vec2::new(c * v.x - s * v.y, s * v.x + c * v.y);
    for k in 0..3 {
        r[k] = rot(r[k]);
        p[k] = rot(p[k]);
    }
    let com01 = (r[0] * m[0] + r[1] * m[1]) / m01;
    let lam = r[2] - com01;
    if lam.y < -DELTA_LAM {
        for k in 0..3 {
            r[k] = Vec2::new(r[k].x, -r[k].y);
            p[k] = Vec2::new(p[k].x, -p[k].y);
        }
    }
}

/// `l = sqrt(I)`, then `r /= l` and `p *= sqrt(l)`.
///
/// **Note the asymmetric powers** — that is what makes the transformation canonical, and getting
/// them equal would leave the Hamiltonian's scaling wrong while every configuration still looked
/// right. `I` here is the mass-weighted Jacobi inertia, matching [`shape::inertia`].
pub fn scale_gauge(
    r: &mut [Vec2<f64>; 3],
    p: &mut [Vec2<f64>; 3],
    m: &[f64; 3],
) -> Option<Degenerate> {
    let i = shape::inertia(r, m);
    if !(i > 0.0) || !i.is_finite() {
        return Some(Degenerate::InertiaZero);
    }
    let l = i.sqrt();
    for k in 0..3 {
        r[k] = r[k] / l;
        p[k] = p[k] * l.sqrt();
    }
    None
}

/// `eta_E = sqrt((E* - U)/K0)`, applied to the momenta. Feasible only when `E* >= U`.
///
/// **Forbidden on `(Lz,E)` and `(Lz,K)`**, where energy is a chart coordinate or is enforced by
/// the momentum construction — applying it there collapses the energy axis. Enforced in code
/// rather than prose: see `Chart::forbids_energy_normalisation` and the validation in
/// `grid::validate_chart`.
pub fn normalise_energy(
    r: &[Vec2<f64>; 3],
    p: &mut [Vec2<f64>; 3],
    m: &[f64; 3],
    e_star: f64,
) -> Option<Degenerate> {
    let v: [Vec2<f64>; 3] = [p[0] / m[0], p[1] / m[1], p[2] / m[2]];
    let k0 = crate::physics::energy::kinetic(&v, m);
    let u = crate::physics::energy::potential(r, m, 0.0);
    if e_star < u || !(k0 > 0.0) {
        return Some(Degenerate::EnergyInfeasible);
    }
    let eta = ((e_star - u) / k0).sqrt();
    for x in p.iter_mut() {
        *x = *x * eta;
    }
    None
}

// ---------------------------------------------------------------------------------------------
// The whole of D, for the latent chart
// ---------------------------------------------------------------------------------------------

/// `D` applied to an 8D latent coordinate.
pub fn decode(z: &Latent) -> Decoded {
    let (m, mut flag) = masses(z.z_mu);
    let (alpha, beta) = angles(z.z_alpha, z.z_beta);
    let mut r = config(alpha, beta, &m);
    let mut p = momenta(z.z_q, &m);
    canonicalise(&mut r, &mut p, &m);
    if let Some(f) = scale_gauge(&mut r, &mut p, &m) {
        flag = flag.or(Some(f));
    }
    Decoded { ic: to_ic(&r, &p, &m), flag }
}

/// Momenta to velocities, and the pair into an [`Ic`]. The integrator's interface is velocities.
pub fn to_ic(r: &[Vec2<f64>; 3], p: &[Vec2<f64>; 3], m: &[f64; 3]) -> Ic<f64> {
    Ic { m: *m, s: Cart { r: *r, v: [p[0] / m[0], p[1] / m[1], p[2] / m[2]] } }
}

// ---------------------------------------------------------------------------------------------
// 2.2 The deterministic momentum construction, for the invariant charts
// ---------------------------------------------------------------------------------------------

fn ang_mom(r: &[Vec2<f64>; 3], w: &[Vec2<f64>; 3], m: &[f64; 3]) -> f64 {
    (0..3).map(|k| m[k] * (r[k].x * w[k].y - r[k].y * w[k].x)).sum()
}

fn mass_norm_sq(w: &[Vec2<f64>; 3], m: &[f64; 3]) -> f64 {
    (0..3).map(|k| m[k] * w[k].norm_sq()).sum()
}

/// Momenta realising a target `Lz` and kinetic energy `K*`, deterministically.
///
/// Three steps, per the reference §2.2: the minimal-energy rigid rotation that realises `Lz`;
/// a direction field that adds energy without changing `Lz`, projected free of COM drift and of
/// angular momentum; then a mix to hit `K*` exactly.
///
/// **The seed family is tried in order and the best-conditioned qualifying one is chosen**, not
/// the first: `(rho, 0)`, `(0, lam)`, `(J rho, 0)`, `(0, J lam)`. Returning
/// [`Degenerate::SeedsExhausted`] only when all four fail is what keeps this a decode rather than
/// a rejection.
///
/// Three constraints, all hit exactly: `sum p = 0`, `Lz(p) = lz`, `K(p) = k_star`. All three are
/// asserted over a deterministic spread including the parabola boundary where `K* -> K_min`,
/// which is the conditioning-sensitive corner.
///
/// **`lz` is about the centre of mass**, and the configuration is centred internally. See the
/// note in the body: the reference assumes a COM-centred input and returns a drifting system
/// without one.
pub fn momenta_for(
    lz: f64,
    k_star: f64,
    r_in: &[Vec2<f64>; 3],
    m: &[f64; 3],
) -> Result<[Vec2<f64>; 3], Degenerate> {
    // **Everything below is in the centre-of-mass frame, and `lz` is about the centre of mass.**
    // The reference's §2.2 assumes the COM is at the origin, which every decoded configuration
    // satisfies -- but the rigid-rotation step is `v = omega J r`, whose total momentum is
    // `omega J (M R_com)`. Handed an off-centre configuration it silently returns momenta whose
    // sum is not zero, and the caller gets a drifting system that conserves everything else.
    // Caught by `a_degenerate_primary_seed_falls_back_rather_than_giving_up`, whose test
    // configuration is not COM-centred.
    let mtot = m[0] + m[1] + m[2];
    let c = (r_in[0] * m[0] + r_in[1] * m[1] + r_in[2] * m[2]) / mtot;
    let r = &[r_in[0] - c, r_in[1] - c, r_in[2] - c];

    let i = (0..3).map(|k| m[k] * r[k].norm_sq()).sum::<f64>();
    if !(i > 0.0) || !i.is_finite() {
        return Err(Degenerate::InertiaZero);
    }
    let j = |v: Vec2<f64>| Vec2::new(-v.y, v.x);

    // (i) minimal-energy rigid rotation realising Lz
    let omega = lz / i;
    let v_l: [Vec2<f64>; 3] = [j(r[0]) * omega, j(r[1]) * omega, j(r[2]) * omega];
    let k_min = lz * lz / (2.0 * i);
    if k_star < k_min {
        return Err(Degenerate::BelowKMin);
    }

    // (ii) a direction field that adds energy without changing Lz
    let m01 = m[0] + m[1];
    let com01 = (r[0] * m[0] + r[1] * m[1]) / m01;
    let rho = r[1] - r[0];
    let lam = r[2] - com01;
    let seeds: [(Vec2<f64>, Vec2<f64>); 4] =
        [(rho, Vec2::zero()), (Vec2::zero(), lam), (j(rho), Vec2::zero()), (Vec2::zero(), j(lam))];

    let mut best: Option<([Vec2<f64>; 3], f64)> = None;
    for (rd, ld) in seeds {
        let mut w: [Vec2<f64>; 3] = [
            ld * (-m[2] / mtot) - rd * (m[1] / m01),
            ld * (-m[2] / mtot) + rd * (m[0] / m01),
            ld * (m01 / mtot),
        ];
        // project out COM drift
        let c = (w[0] * m[0] + w[1] * m[1] + w[2] * m[2]) / mtot;
        for x in w.iter_mut() {
            *x = *x - c;
        }
        // project out angular momentum
        let beta_l = ang_mom(r, &w, m) / i;
        for k in 0..3 {
            w[k] = w[k] - j(r[k]) * beta_l;
        }
        let n2 = mass_norm_sq(&w, m);
        if n2 > SEED_EPS && n2.is_finite() {
            let n = n2.sqrt();
            let unit = [w[0] / n, w[1] / n, w[2] / n];
            // Largest ‖w‖ wins, for conditioning — not the first that qualifies.
            if best.as_ref().map_or(true, |(_, b)| n > *b) {
                best = Some((unit, n));
            }
        }
    }
    let (wu, _) = best.ok_or(Degenerate::SeedsExhausted)?;

    // (iii) mix to the target kinetic energy
    let a = (2.0 * (k_star - k_min)).max(0.0).sqrt();
    let mut p = [Vec2::zero(); 3];
    for k in 0..3 {
        p[k] = (v_l[k] + wu[k] * a) * m[k];
    }
    Ok(p)
}

// ---------------------------------------------------------------------------------------------
// 4 The Burrau family
// ---------------------------------------------------------------------------------------------

/// Euclid's parametrisation continued to real `nu = n/m in (0,1)`: the bifurcation strip.
///
/// Returns `(masses, positions)` in the reference's convention — right angle at the origin,
/// normalised by the hypotenuse, each mass equal to its opposite side.
///
/// **This is not the repo's Burrau convention**, and the difference is a body relabelling plus
/// the scale gauge, not a different system. The repo uses `MASSES = [3,4,5]` at
/// `[(1,3), (-2,-1), (1,-1)]`; the reference gives `(c,b,a)/(a+b+c) = (5,4,3)/12` at
/// `[(0,0), (a/c,0), (0,b/c)]`. Both are "mass equals opposite side". `tests/charts.rs` pins the
/// correspondence against `tests/burrau_constants.rs`'s gates (`M = 12`, `R = 2.2361`,
/// `E = -12.8167`) rather than assuming it, and every dump states which convention it is in.
pub fn burrau_family(nu: f64) -> ([f64; 3], [Vec2<f64>; 3]) {
    let (a, b, c) = (1.0 - nu * nu, 2.0 * nu, 1.0 + nu * nu);
    let s = a + b + c;
    let m = [c / s, b / s, a / s];
    let r = [Vec2::new(0.0, 0.0), Vec2::new(a / c, 0.0), Vec2::new(0.0, b / c)];
    (m, r)
}

/// The `nu` that reproduces a primitive triple from Euclid's `(m, n)`. `nu = n/m`.
pub fn nu_of(m: u32, n: u32) -> f64 {
    n as f64 / m as f64
}
