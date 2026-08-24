//! Benettin FTLE and the diffusion regression — a port of `reference/tb_ftle.py`.
//!
//! **Ported, not re-derived.** The renormalisation bookkeeping is exactly the class of algebra
//! that fails silently: a shadow that is rescaled at the wrong moment, or normalised against the
//! wrong separation, still produces a smooth plausible field.
//!
//! # What this is for
//!
//! The production colour scheme is bivariate: hue from the shape sphere, **lightness from a
//! scalar**. That makes it a criterion question rather than a presentation one — the criterion
//! asks "would splitting change what we display?", so what is displayed decides what the
//! criterion should measure. `spread_shape` maps to hue, so that half is aligned; if lightness
//! carries diffusion or FTLE, the criterion is currently blind to changes in it.
//!
//! # Renormalisation is what stops it saturating
//!
//! Without it the shadow separates until it fills the accessible space and `log(d/d0)/T` decays
//! toward zero — reporting `lambda ~ 0` for the *most* chaotic regions, which is the inversion
//! this project has now met three times. `n_renorm` is returned so a caller can assert it is
//! nonzero; an FTLE built from zero renormalisations is the saturated case in new clothing, and
//! `ftle` is **NaN** in that case rather than a number.
//!
//! # Two honest limitations, stated because they bound what can be concluded
//!
//! 1. **This sits on the plain leapfrog, not on Aarseth-Zare.** `tb_ftle.py` is built on
//!    `tb.py`, so that is the pair with a reference and the pair the cross-check compares. The
//!    unregularised integrator fails on close encounters by construction, so an FTLE from this
//!    path is trustworthy only where the trajectory is not near a collision — which is exactly
//!    where it is least interesting. Carrying a shadow through AZ is a separate step with no
//!    reference behind it, and it is validated against this one where **both** resolve.
//! 2. **The perturbation direction is a parameter, not an RNG.** `tb_ftle.py` draws it from
//!    `numpy.random.default_rng(seed).normal`, whose Ziggurat is not ported. Reproducing the
//!    stream is not required: the direction only seeds the shadow, and a comparison needs both
//!    sides to use the *same* one. So the direction is passed in, and the cross-check hands an
//!    analytic vector to both sides. Nothing about the RNG enters the validated path.

use crate::physics::{energy, newton, Cart};
use crate::{Real, Vec2};

/// Knobs, with the reference's defaults.
#[derive(Clone, Copy, Debug)]
pub struct FtleOpts {
    /// Plummer softening `eps`, squared. The reference's default `eps = 0.03`.
    pub eps2: f64,
    /// Shadow separation held by renormalisation.
    pub d0: f64,
    /// Steps between renormalisations. The reference's 200.
    pub renorm_every: usize,
    /// Steps between diffusion-regression samples. The reference's 25.
    pub sample_every: usize,
}

impl Default for FtleOpts {
    fn default() -> Self {
        Self { eps2: 0.03 * 0.03, d0: 1e-8, renorm_every: 200, sample_every: 25 }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FtleOut<T> {
    pub state: Cart<T>,
    /// `S / T`. **NaN when no renormalisation completed** — never 0, which would be
    /// indistinguishable from a perfectly regular trajectory.
    pub ftle: T,
    /// OLS slope of `log(inertia)` against `t`. NaN on a degenerate design matrix.
    pub diffusion: T,
    pub d_min: T,
    pub n_renorm: u64,
    /// Accumulated `log(d/d0)`, returned so `ftle` can be recomputed at another horizon.
    pub s_accum: T,
    pub steps: usize,
    pub finite: bool,
}

/// Phase-space separation of the shadow: position and velocity together.
fn phase_dist<T: Real>(dr: &[Vec2<T>; 3], dv: &[Vec2<T>; 3]) -> T {
    let mut acc = T::zero();
    for k in 0..3 {
        acc += dr[k].norm_sq() + dv[k].norm_sq();
    }
    acc.sqrt()
}

/// A deterministic unit direction in the 6-dimensional position space.
///
/// **Not numpy's stream**, and not claimed to be. Used where the direction only needs to be
/// reproducible, never where two implementations are being compared.
pub fn unit_perturbation<T: Real>(seed: u64) -> [Vec2<T>; 3] {
    let mut z = crate::rng::SplitMix64::new(seed);
    let mut raw = [0.0f64; 6];
    for x in raw.iter_mut() {
        // Box-Muller from two uniforms; the direction is all that survives normalisation.
        let u1 = ((z.next_u64() >> 11) as f64 / (1u64 << 53) as f64).max(1e-300);
        let u2 = (z.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
        *x = (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
    }
    normalise(raw)
}

/// Normalise a flat 6-vector into the `[Vec2; 3]` layout, matching the reference's
/// `pert /= norm(pert.reshape(n, -1))` — one norm over **all six components**, not per body.
pub fn normalise<T: Real>(raw: [f64; 6]) -> [Vec2<T>; 3] {
    let nrm = raw.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-300);
    [
        Vec2::new(T::lit(raw[0] / nrm), T::lit(raw[1] / nrm)),
        Vec2::new(T::lit(raw[2] / nrm), T::lit(raw[3] / nrm)),
        Vec2::new(T::lit(raw[4] / nrm), T::lit(raw[5] / nrm)),
    ]
}

/// Main trajectory, Benettin shadow and diffusion regression, at **fixed** `dt`.
///
/// Fixed step, not the adaptive one in `integrate::leapfrog`: `tb_ftle.integrate_full` takes
/// `steps = round(t_max/dt)` and the renormalisation cadence is counted in *steps*, so an
/// adaptive step would make `renorm_every` mean a different interval on every trajectory and
/// the comparison would not be against the same estimator.
pub fn integrate_full<T: Real>(
    s0: Cart<T>,
    m: &[T; 3],
    t_max: T,
    dt: T,
    opts: &FtleOpts,
    pert: &[Vec2<T>; 3],
) -> FtleOut<T> {
    let eps2 = T::lit(opts.eps2);
    let d0 = T::lit(opts.d0);
    let half = T::lit(0.5);

    let mut r = s0.r;
    let mut v = s0.v;
    // The shadow is displaced in POSITION only and shares the velocity, as the reference has it
    // (`rs = r + d0*pert; vs = v.copy()`), even though the separation is measured in full phase
    // space. Transcribed rather than tidied: changing where the shadow starts changes the
    // transient, and the transient is inside the short horizon this project runs at.
    let mut rs = [
        r[0] + pert[0] * d0,
        r[1] + pert[1] * d0,
        r[2] + pert[2] * d0,
    ];
    let mut vs = v;

    let mut s_accum = T::zero();
    let mut n_renorm = 0u64;

    // O(1) regression accumulators — no history retained, which is what makes this affordable
    // per footprint.
    let (mut cnt, mut st, mut stt) = (T::zero(), T::zero(), T::zero());
    let (mut sy, mut sty) = (T::zero(), T::zero());

    let mut d_min = T::infinity();
    let mut a = newton::accel(&r, m, eps2);
    let mut a_s = newton::accel(&rs, m, eps2);

    let steps = (t_max / dt).to_f64().unwrap().round().max(0.0) as usize;
    let mut finite = true;

    for s in 0..steps {
        for k in 0..3 {
            v[k] += a[k] * (half * dt);
            r[k] += v[k] * dt;
        }
        a = newton::accel(&r, m, eps2);
        for k in 0..3 {
            v[k] += a[k] * (half * dt);
        }

        for k in 0..3 {
            vs[k] += a_s[k] * (half * dt);
            rs[k] += vs[k] * dt;
        }
        a_s = newton::accel(&rs, m, eps2);
        for k in 0..3 {
            vs[k] += a_s[k] * (half * dt);
        }

        // Test explicitly: NaN >= x is false, so a diverged trajectory satisfies no guard and
        // burns its whole budget. Measured elsewhere in this project at 354 s against 3 s.
        if !r[0].x.is_finite() || !rs[0].x.is_finite() {
            finite = false;
            break;
        }

        if s % opts.renorm_every == 0 && s > 0 {
            let dr = [rs[0] - r[0], rs[1] - r[1], rs[2] - r[2]];
            let dv = [vs[0] - v[0], vs[1] - v[1], vs[2] - v[2]];
            let d = phase_dist(&dr, &dv);
            if d > T::TINY {
                s_accum += (d / d0).ln();
                n_renorm += 1;
                let sc = d0 / d.max(T::TINY);
                for k in 0..3 {
                    rs[k] = r[k] + dr[k] * sc;
                    vs[k] = v[k] + dv[k] * sc;
                }
            }
        }

        if s % opts.sample_every == 0 {
            let t = T::lit((s + 1) as f64) * dt;
            let y = energy::inertia(&r, m).max(T::TINY).ln();
            cnt += T::one();
            st += t;
            stt += t * t;
            sy += y;
            sty += t * y;
            let d = newton::pair_dists(&r);
            let step_min = d[0].min(d[1]).min(d[2]);
            if step_min < d_min {
                d_min = step_min;
            }
        }
    }

    let total_t = T::lit(steps as f64) * dt;
    let ftle = if n_renorm > 0 {
        s_accum / total_t.max(T::lit(1e-12))
    } else {
        T::nan()
    };
    let den = cnt * stt - st * st;
    let diffusion = if den.abs() > T::lit(1e-12) {
        (cnt * sty - st * sy) / den
    } else {
        T::nan()
    };

    FtleOut {
        state: Cart { r, v },
        ftle,
        diffusion,
        d_min,
        n_renorm,
        s_accum,
        steps,
        finite,
    }
}
