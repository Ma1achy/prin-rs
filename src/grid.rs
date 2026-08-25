//! The initial-condition slice: a 2D box of chart coordinates, decoded to full states.
//!
//! Row order matches `tb.burrau_grid`: `meshgrid(indexing='xy')` then C-order `ravel()`,
//! i.e. **`index = jy*nx + jx`, x fastest**. The cross-check compares row by row, so this
//! ordering is load-bearing and is asserted in a test against a hardcoded case.
//!
//! Initial conditions are always built in `f64` and cast down. Generating them separately
//! per precision would make an IC difference indistinguishable from a genuine f32
//! arithmetic effect — the decomposition that made the earlier f32 investigation
//! interpretable at all.
//!
//! **The chart is a parameter.** Every experiment before the vertical slice varied one
//! body's position in the plane, which is an *affine* decode: `J_D` is constant, so the
//! linearised path `x = x0 + J_D.delta` is exact rather than approximate and "where does
//! the linearisation start to matter" answers "never" at every depth. [`Chart::Shape`]
//! exists so that question has something to measure. [`Chart::BodyPlane`] is the old
//! behaviour, preserved **bitwise** and asserted so in `tests/vertical_slice.rs`.

use crate::physics::{burrau, decoder, shape, Cart, Ic};
use crate::{Real, Vec2};

/// How a chart coordinate `(u, v)` becomes a three-body state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Chart {
    /// One body's position over a 2D box; the other two fixed at their Burrau values.
    /// **Affine** — this is every result before the vertical slice.
    BodyPlane,
    /// An oblique 2-plane in the 6D position space: `origin + u*U + v*V`, where `U` and `V`
    /// may move more than one body. Still affine, but no longer axis-aligned, which is what
    /// §3.5 asks about. [`Chart::BodyPlane`] is the special case
    /// `U = e(body, x)`, `V = e(body, y)`.
    Plane {
        origin: Cart<f64>,
        u: [Vec2<f64>; 3],
        v: [Vec2<f64>; 3],
    },
    /// The shape sphere: `(u, v)` is a tangent step from `n0`, mapped by the exponential
    /// map and inverted through the Hopf map at fixed scale and fibre phase.
    ///
    /// **Nonlinear.** The exp map's `cos|t|`, `sin|t|/|t|` is where the curvature lives, and
    /// it is the only chart here on which the linearised decoder is an approximation at all.
    Shape {
        n0: [f64; 3],
        e1: [f64; 3],
        e2: [f64; 3],
        inertia: f64,
        phase: f64,
        /// The masses the Hopf inverse is taken with. Was hard-wired to `burrau::MASSES`
        /// inside the decode; it is a parameter now because the mass-varying charts made a
        /// global mass wrong in general, and because the Hopf inverse is mass-dependent.
        m: [f64; 3],
    },

    /// **The 8D latent chart** — the reference coordinate system, of which the axis-aligned
    /// planes and the oblique ones are both instances.
    ///
    /// `z(u, v) = z0 + u*q1 + v*q2`. The reference writes `(2u-1)*s_u*q1` over `[0,1]^2`;
    /// here the slice already supplies a signed box in chart coordinates, so `(u, v)` are the
    /// offsets directly and the box's half-width plays the part of `s_u`. Same family, one
    /// fewer place for a factor of two to hide.
    ///
    /// **Report the basis pair, or the slice is not reproducible** — [`Chart::params`] writes
    /// it into every dump header.
    Latent {
        z0: decoder::Latent,
        q1: [f64; 8],
        q2: [f64; 8],
    },

    /// **The Burrau family `(nu, K)`** — the bifurcation strip.
    ///
    /// `u` sweeps Euclid's `nu = n/m` continued to the reals, so the primitive Pythagorean
    /// triples sit at a countable set of points on a continuous curve; `v` sweeps kinetic
    /// energy through the warp `K = k_max * v^gamma_k`, at `Lz = 0`. At `v = 0` this is the
    /// classical rest start.
    ///
    /// Answers directly whether the right-angle constraint is dynamically special or merely
    /// convenient: the Burrau configurations are a **1D curve inside a 2D map**.
    BurrauFamily {
        nu_lo: f64,
        nu_hi: f64,
        k_max: f64,
        gamma_k: f64,
    },

    /// **The invariant-momentum chart `(Lz, K)`** — geometry and mass fixed, both axes
    /// conserved quantities, with the feasibility warp that makes every pixel valid.
    ///
    /// `K(t) = k_max * t^gamma_k`, `L_max(t) = sqrt(2 I K(t))`, `Lz(s,t) = (2s-1) L_max(t)`.
    /// The warp maps the unit square onto the interior of the feasible parabola, which is why
    /// it exists rather than clamping.
    ///
    /// **`(Lz, E)` is not a second chart.** The reference lists the two separately, but its own
    /// construction parameterises both by `K(t) = K_max t^gamma`: for `(Lz,E)` it then sets
    /// `E = U + K(t)`, which is a *relabelling of the axis*, not a different sweep. The two
    /// produce bitwise identical initial conditions at equal `gamma_k`, and `tests/charts.rs`
    /// asserts it. `report_e` carries the label so a dump says which axis was intended;
    /// the only thing that makes the two genuinely differ is `gamma_k`.
    ///
    /// `forbids_energy_normalisation` is **true** here: energy is a chart coordinate, and
    /// normalising it would collapse the axis. Enforced by [`Chart::validate`], not by prose.
    Invariant {
        base: decoder::Latent,
        k_max: f64,
        gamma_k: f64,
        report_e: bool,
    },

    /// **The mass simplex** at fixed geometry and fixed momentum coordinates.
    ///
    /// Barycentric coordinates onto `[0,1]^2` with a shear, then blended `margin` of the way
    /// toward the centroid so no mass reaches zero — a zero mass is not a three-body system and
    /// would be a degenerate label on every edge pixel rather than a measurement.
    MassSimplex {
        z_alpha: f64,
        z_beta: f64,
        z_q: [f64; 4],
        margin: f64,
    },
}

/// What region of the plane a chart's coordinates are defined on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Domain {
    /// Any `(u, v)`.
    Free,
    /// `[0,1]^2`. The feasibility warps and the simplex map are only onto inside it, so a slice
    /// that strays outside is not sampling what it claims to.
    Unit,
}

impl Default for Chart {
    fn default() -> Self {
        Chart::BodyPlane
    }
}

impl Chart {
    pub fn name(&self) -> &'static str {
        match self {
            Chart::BodyPlane => "body_plane",
            Chart::Plane { .. } => "plane",
            Chart::Shape { .. } => "shape",
            Chart::Latent { .. } => "latent",
            Chart::BurrauFamily { .. } => "burrau",
            Chart::Invariant { report_e, .. } => {
                if *report_e {
                    "invariant_lz_e"
                } else {
                    "invariant_lz_k"
                }
            }
            Chart::MassSimplex { .. } => "mass_simplex",
        }
    }

    /// Every parameter of this chart, for the dump header.
    ///
    /// `Chart::name()` alone is not enough: a `Plane`'s basis, a `Shape`'s `(n0, I, phase)` and
    /// a `Latent`'s `(z0, q1, q2)` are all free, so two dumps with the same name can be
    /// different configurations. The reference's rule is *report the pair used, or the slice is
    /// not reproducible*, and an oblique latent chart makes that binding.
    pub fn params(&self) -> String {
        match self {
            Chart::BodyPlane => "-".into(),
            Chart::Plane { origin, u, v } => format!(
                "origin_r={:?} u={:?} v={:?}",
                origin.r.map(|p| (p.x, p.y)),
                u.map(|p| (p.x, p.y)),
                v.map(|p| (p.x, p.y))
            ),
            Chart::Shape { n0, e1, e2, inertia, phase, m } => {
                format!("n0={n0:?} e1={e1:?} e2={e2:?} I={inertia:?} phase={phase:?} m={m:?}")
            }
            Chart::Latent { z0, q1, q2 } => format!("z0={z0:?} q1={q1:?} q2={q2:?}"),
            Chart::BurrauFamily { nu_lo, nu_hi, k_max, gamma_k } => {
                format!("nu=[{nu_lo:?},{nu_hi:?}] k_max={k_max:?} gamma_k={gamma_k:?}")
            }
            Chart::Invariant { base, k_max, gamma_k, report_e } => {
                format!("base={base:?} k_max={k_max:?} gamma_k={gamma_k:?} report_e={report_e}")
            }
            Chart::MassSimplex { z_alpha, z_beta, z_q, margin } => {
                format!("z_alpha={z_alpha:?} z_beta={z_beta:?} z_q={z_q:?} margin={margin:?}")
            }
        }
    }

    /// The natural box half-width for this chart family.
    ///
    /// **Not one number.** A `BodyPlane` coordinate is a body position in Burrau units; a
    /// `Latent` coordinate is a sigmoid pre-image. `0.05` of one is nothing like `0.05` of the
    /// other, and a single shared default therefore *silently means two different things*. That
    /// is exactly how the four GLSL presets shipped at a 3x crop: the reference UI reads
    /// `Slice +/- 3.0e+0` and the port rendered `half = 1.0`, which spans 46% of the azimuth
    /// against 90%.
    ///
    /// `0.45` for the `Domain::Unit` charts is the value already in use at centre `(0.5, 0.5)`;
    /// it is recorded here rather than changed.
    pub fn default_half(&self) -> f64 {
        match self {
            Chart::BodyPlane | Chart::Plane { .. } | Chart::Shape { .. } => 0.05,
            Chart::Latent { .. } => 3.0,
            Chart::BurrauFamily { .. } | Chart::Invariant { .. } | Chart::MassSimplex { .. } => {
                0.45
            }
        }
    }

    /// The coordinate region this chart is defined on.
    pub fn domain(&self) -> Domain {
        match self {
            Chart::BurrauFamily { .. } | Chart::Invariant { .. } | Chart::MassSimplex { .. } => {
                Domain::Unit
            }
            _ => Domain::Free,
        }
    }

    /// Whether applying `D`'s energy normalisation to this chart would destroy its own axis.
    ///
    /// True for the invariant-momentum chart, where energy is a coordinate. **Enforced**, not
    /// documented: [`Chart::validate`] refuses the combination and a test asserts the refusal.
    pub fn forbids_energy_normalisation(&self) -> bool {
        matches!(self, Chart::Invariant { .. })
    }

    /// Refuse a configuration that cannot mean what it says.
    ///
    /// Two things it catches: a non-zero `E*` on a chart whose axis *is* energy, and a slice box
    /// that leaves a `Unit`-domain chart's square, where the feasibility warp's guarantee does
    /// not hold and the chart is sampling something it does not describe.
    pub fn validate(&self, e_star: f64, cx: f64, cy: f64, half: f64) -> Result<(), String> {
        if self.forbids_energy_normalisation() && e_star != 0.0 {
            return Err(format!(
                "chart `{}` carries energy as a coordinate; energy normalisation to E* = {e_star} \
                 would collapse that axis",
                self.name()
            ));
        }
        if self.domain() == Domain::Unit {
            // Copies span the whole cell edge to edge at `jitter_frac = 0.5`, so the box has to
            // sit inside the square with room for them; the jitter itself reflects at the
            // boundary rather than clamping (see `jitter::reflect_into_unit`).
            for (name, c) in [("u", cx), ("v", cy)] {
                if c - half < 0.0 || c + half > 1.0 {
                    return Err(format!(
                        "chart `{}` has domain [0,1]^2 but the slice spans {name} in \
                         [{}, {}]",
                        self.name(),
                        c - half,
                        c + half
                    ));
                }
            }
        }
        Ok(())
    }

    /// Is the decode affine? If so the linearised path is exact and its curvature term is
    /// **structurally zero** — a fact to report as structural, never as a measurement.
    pub fn is_affine(&self) -> bool {
        matches!(self, Chart::BodyPlane | Chart::Plane { .. })
    }

    /// The `Plane` that reproduces `BodyPlane` for `body`, used to assert the two agree.
    pub fn plane_for_body(body: usize) -> Chart {
        let mut u = [Vec2::zero(); 3];
        let mut v = [Vec2::zero(); 3];
        u[body] = Vec2::new(1.0, 0.0);
        v[body] = Vec2::new(0.0, 1.0);
        // The origin's varying body sits at zero; `u`, `v` carry the whole position.
        let mut origin = burrau::state::<f64>();
        origin.r[body] = Vec2::zero();
        Chart::Plane { origin, u, v }
    }

    /// The shape chart through Burrau's own configuration: same point on the sphere, same
    /// scale, deterministic tangent frame. A slice of it is a slice of the shape sphere
    /// around the reference triangle.
    pub fn shape_at_burrau(phase: f64) -> Chart {
        let m = burrau::MASSES;
        let r: [Vec2<f64>; 3] = [
            Vec2::new(burrau::R0[0][0], burrau::R0[0][1]),
            Vec2::new(burrau::R0[1][0], burrau::R0[1][1]),
            Vec2::new(burrau::R0[2][0], burrau::R0[2][1]),
        ];
        let n0 = shape::shape_vec(&r, &m);
        let (e1, e2) = shape::tangent_frame(n0);
        Chart::Shape { n0, e1, e2, inertia: shape::inertia(&r, &m), phase, m }
    }

    /// An axis-aligned latent plane: the two named coordinates of the reference's table.
    ///
    /// `shape` is `(z_alpha, z_beta)`, `inner momentum` is `(z_q0, z_q1)`, `outer momentum` is
    /// `(z_q2, z_q3)`, `mass` is `(z_mu1, z_mu2)`, `mixed` is `(z_alpha, z_q2)`. There are
    /// `C(8,2) = 28` such planes and every one is an instance of [`Chart::Latent`].
    pub fn latent_axes(z0: decoder::Latent, i: usize, j: usize) -> Chart {
        assert!(i < 8 && j < 8 && i != j, "latent axes must be two distinct coordinates");
        let mut q1 = [0.0; 8];
        let mut q2 = [0.0; 8];
        q1[i] = 1.0;
        q2[j] = 1.0;
        Chart::Latent { z0, q1, q2 }
    }

    /// An oblique latent plane from two seed vectors, orthonormalised by Gram-Schmidt.
    ///
    /// Deterministic in the seeds so the plane is reproducible, and [`Chart::params`] writes the
    /// resulting pair into the header — a plane whose basis is not recorded is not a slice
    /// anyone can repeat.
    pub fn latent_oblique(z0: decoder::Latent, a: [f64; 8], b: [f64; 8]) -> Chart {
        let dot = |x: &[f64; 8], y: &[f64; 8]| (0..8).map(|k| x[k] * y[k]).sum::<f64>();
        // A zero seed, or a `b` parallel to `a`, divides by zero and hands back a NaN basis. A
        // NaN basis decodes every pixel identically, `ensemble_spread` reads exactly zero, and
        // the criterion reports the quad perfectly resolved — a collapsed decode wearing the
        // face of a tidy answer. Refuse it here instead.
        const SEED_EPS: f64 = 1e-12;
        let mut q1 = a;
        let n1 = dot(&q1, &q1).sqrt();
        assert!(n1 > SEED_EPS, "latent_oblique: first seed has norm {n1:e}, cannot be normalised");
        for x in q1.iter_mut() {
            *x /= n1;
        }
        let mut q2 = b;
        let p = dot(&q2, &q1);
        for k in 0..8 {
            q2[k] -= p * q1[k];
        }
        let n2 = dot(&q2, &q2).sqrt();
        assert!(
            n2 > SEED_EPS,
            "latent_oblique: second seed is parallel to the first (residual norm {n2:e}); \
             the two seeds span a line, not a plane"
        );
        for x in q2.iter_mut() {
            *x /= n2;
        }
        Chart::Latent { z0, q1, q2 }
    }

    // ---- The GLSL reference's four default slices ---------------------------------------
    //
    // `Ma1achy/principia-ii`, `src/state.ts:71-76`. Constructed here rather than at each call
    // site because the basis was wrong in one of them and the literal appeared three times: the
    // gallery and two tests. A correction has to land once.
    //
    // All four sit at `z0 = 0`, which decodes to the equilateral Lagrange configuration -- a
    // named physical state at the centre of every one of these images. Their natural extent is
    // [`Chart::default_half`], `3.0`, from the reference UI's `Slice +/- 3.0e+0`.

    /// `shape`: the two configuration coordinates.
    pub fn preset_shape() -> Chart {
        Chart::latent_axes(decoder::Latent::default(), 0, 1)
    }

    /// `prho`: the inner momentum pair. A constant-**configuration** slice -- positions in this
    /// decode do not depend on the momentum coordinates at all, so every pixel is the same
    /// triangle released with a different initial velocity, and `spread_shape` at `t = 0` is
    /// identically zero across the whole slice. Any structure in it is purely momentum-driven,
    /// which makes it the control that separates configuration effects from momentum effects.
    pub fn preset_prho() -> Chart {
        Chart::latent_axes(decoder::Latent::default(), 2, 3)
    }

    /// `plambda`: the outer momentum pair. Constant-configuration, exactly as [`Self::preset_prho`].
    pub fn preset_plambda() -> Chart {
        Chart::latent_axes(decoder::Latent::default(), 4, 5)
    }

    /// `shape_pl`: **the only preset with a cross-coupling**, and the only one that can be got
    /// wrong in this particular way.
    ///
    /// Constructed directly and **not** through [`Chart::latent_oblique`]: the reference's basis
    /// is un-normalised (each direction has norm `sqrt 2`) and Gram-Schmidt would quietly render
    /// a different slice while looking like a tidy-up. `tests/charts.rs` pins the norm.
    ///
    /// **The pairing is by GLSL SLOT.** The reference is `q1 = e0 + e6`, `q2 = e1 + e7`, and in
    /// *its* indexing (`z0 = beta`, `z1 = alpha`, `z6/z7 = pLambda.x/y`) that pairs beta with
    /// `pLambda.x` and alpha with `pLambda.y`. This module renumbers alpha and beta into the
    /// spec's order, and **must carry their momentum partners with them** -- so the *pair
    /// assignment* transposes and each pair stays intact:
    ///
    /// ```text
    ///   q1 (horizontal) = e_alpha + e_pLambda_y = e0 + e5      // the GLSL's q2
    ///   q2 (vertical)   = e_beta  + e_pLambda_x = e1 + e4      // the GLSL's q1
    /// ```
    ///
    /// Pairing alpha with `pLambda.x` instead is a **genuinely different 2-plane** through the 8D
    /// space, not a reorientation of the same one, and transposing `q1`/`q2` does not recover it
    /// -- that gives `e_beta + e_pLy`, `e_alpha + e_pLx`, still crossed. It renders as *twisted*
    /// rather than tilted, because the coupling sets how momentum co-varies with configuration
    /// across the slice and the two pairings give different shears.
    pub fn preset_shape_pl() -> Chart {
        Chart::Latent {
            z0: decoder::Latent::default(),
            q1: [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            q2: [0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Slice {
    pub nx: usize,
    pub ny: usize,
    pub cx: f64,
    pub cy: f64,
    pub half: f64,
    /// Which body's position varies under [`Chart::BodyPlane`]: 0, 1 or 2. Retained on every
    /// chart because it names the slice family in dumps and headers.
    pub body: usize,
    pub chart: Chart,
}

impl Slice {
    /// The pre-vertical-slice constructor. Every existing result is a `BodyPlane` slice and
    /// this is the only way one is built, so nothing can acquire a different chart silently.
    pub fn body_plane(nx: usize, ny: usize, cx: f64, cy: f64, half: f64, body: usize) -> Self {
        Slice { nx, ny, cx, cy, half, body, chart: Chart::BodyPlane }
    }

    pub fn with_chart(mut self, chart: Chart) -> Self {
        self.chart = chart;
        self
    }

    pub fn npix(&self) -> usize {
        self.nx * self.ny
    }

    /// Cell widths, **per axis**.
    ///
    /// The reference computes only `hx = 2*half/max(nx-1,1)` and uses it for *both* axes.
    /// Every experiment run against it used a square grid, so it never bit; on any
    /// non-square grid the y-jitter is silently mis-scaled. Fixed here.
    pub fn cell_widths(&self) -> (f64, f64) {
        let hx = 2.0 * self.half / (self.nx.max(2) - 1) as f64;
        let hy = 2.0 * self.half / (self.ny.max(2) - 1) as f64;
        (hx, hy)
    }

    /// One axis of the grid, matching `numpy.linspace` including its exact endpoint.
    fn axis(c: f64, half: f64, n: usize, i: usize) -> f64 {
        let (a, b) = (c - half, c + half);
        if n <= 1 {
            return a;
        }
        if i == n - 1 {
            return b; // numpy sets the last sample exactly, rather than accumulating
        }
        a + (i as f64) * ((b - a) / (n - 1) as f64)
    }

    /// Chart position of pixel `idx`, with `idx = jy*nx + jx`.
    pub fn decode_pos(&self, idx: usize) -> (f64, f64) {
        let jx = idx % self.nx;
        let jy = idx / self.nx;
        (
            Self::axis(self.cx, self.half, self.nx, jx),
            Self::axis(self.cy, self.half, self.ny, jy),
        )
    }

    /// **The decode**: chart coordinate `(u, v)` to a full state, in `f64`.
    ///
    /// Public because the deep-zoom ladder decodes directly, without a grid.
    pub fn decode_state(&self, u: f64, v: f64) -> Ic<f64> {
        decode_state(&self.chart, self.body, u, v)
    }

    /// Nominal (un-jittered) initial condition for pixel `idx`.
    ///
    /// This is copy 0 of the reference's ensemble — `mask[::reps] = False` leaves it
    /// un-jittered and therefore completely seed-independent. That is what makes a
    /// nominal-only cross-check possible with no RNG on either side.
    pub fn nominal<T: Real>(&self, idx: usize) -> Cart<T> {
        let (x, y) = self.decode_pos(idx);
        self.decode_state(x, y).s.cast::<T>()
    }

    /// Nominal initial condition **with its masses**.
    ///
    /// [`Self::nominal`] drops them, which is correct at the many call sites that predate the
    /// chart families and are all Burrau. Anything that decodes a chart which can vary mass must
    /// use this instead, or it integrates the right configuration with the wrong bodies.
    pub fn nominal_ic<T: Real>(&self, idx: usize) -> Ic<T> {
        let (x, y) = self.decode_pos(idx);
        self.decode_state(x, y).cast::<T>()
    }
}

/// The decode, free of any grid. `body` is read only by [`Chart::BodyPlane`].
pub fn decode_state(chart: &Chart, body: usize, u: f64, v: f64) -> Ic<f64> {
    match *chart {
        // Bitwise what this function did before the chart existed. Asserted in a test.
        Chart::BodyPlane => {
            let mut s = burrau::state::<f64>();
            s.r[body] = Vec2::new(u, v);
            Ic { m: burrau::MASSES, s }
        }
        Chart::Plane { origin, u: uu, v: vv } => {
            let mut s = origin;
            for k in 0..3 {
                s.r[k] = s.r[k] + uu[k] * u + vv[k] * v;
            }
            Ic { m: burrau::MASSES, s }
        }
        Chart::Shape { n0, e1, e2, inertia, phase, m } => {
            let t = [
                u * e1[0] + v * e2[0],
                u * e1[1] + v * e2[1],
                u * e1[2] + v * e2[2],
            ];
            let n = shape::exp_map(n0, t);
            let r = shape::from_shape(n, inertia, phase, &m);
            // Released from rest, like every configuration in this project.
            Ic { m, s: Cart { r, v: [Vec2::zero(); 3] } }
        }

        // ---- the reference's five families, all through the shared decoder ----
        Chart::Latent { z0, q1, q2 } => {
            let mut z = z0;
            for k in 0..8 {
                z.set(k, z0.get(k) + u * q1[k] + v * q2[k]);
            }
            decoder::decode(&z).ic
        }

        Chart::BurrauFamily { nu_lo, nu_hi, k_max, gamma_k } => {
            let nu = (nu_lo + (nu_hi - nu_lo) * u.clamp(0.0, 1.0)).clamp(1e-9, 1.0 - 1e-9);
            let (m, mut r) = decoder::burrau_family(nu);
            // Lz = 0: the family is defined by its geometry, and the second axis is energy
            // alone. At v = 0 the momenta vanish and this is the classical rest start.
            let k = k_max * v.clamp(0.0, 1.0).powf(gamma_k);
            let mut p = decoder::momenta_for(0.0, k, &r, &m).unwrap_or([Vec2::zero(); 3]);
            decoder::canonicalise(&mut r, &mut p, &m);
            let _ = decoder::scale_gauge(&mut r, &mut p, &m);
            decoder::to_ic(&r, &p, &m)
        }

        Chart::Invariant { base, k_max, gamma_k, .. } => {
            let (m, _) = decoder::masses(base.z_mu);
            let (alpha, beta) = decoder::angles(base.z_alpha, base.z_beta);
            let r = decoder::config(alpha, beta, &m);
            // The feasibility warp. `K >= 0` and `|Lz| <= sqrt(2 I K)` bound the feasible set
            // to the interior of a parabola with its apex at the rest start; this maps the unit
            // square onto that interior, so **no pixel is infeasible by construction** -- which
            // is why the warp exists rather than a clamp.
            let i: f64 = (0..3).map(|k| m[k] * r[k].norm_sq()).sum();
            let t = v.clamp(0.0, 1.0);
            let k = k_max * t.powf(gamma_k);
            let l_max = (2.0 * i * k).max(0.0).sqrt();
            let lz = (2.0 * u.clamp(0.0, 1.0) - 1.0) * l_max;
            let mut p = decoder::momenta_for(lz, k, &r, &m).unwrap_or([Vec2::zero(); 3]);
            let mut r = r;
            decoder::canonicalise(&mut r, &mut p, &m);
            decoder::to_ic(&r, &p, &m)
        }

        Chart::MassSimplex { z_alpha, z_beta, z_q, margin } => {
            // Square to triangle with a shear, then blended toward the centroid so no mass
            // reaches zero. A zero mass is not a three-body system, and without the margin every
            // edge pixel would carry a degenerate label rather than a measurement.
            let (uu, vv) = (u.clamp(0.0, 1.0), v.clamp(0.0, 1.0));
            let raw = [uu * (1.0 - vv), (1.0 - uu) * (1.0 - vv), vv];
            let g = margin.clamp(0.0, 1.0 / 3.0);
            let m = [
                (1.0 - 3.0 * g) * raw[0] + g,
                (1.0 - 3.0 * g) * raw[1] + g,
                (1.0 - 3.0 * g) * raw[2] + g,
            ];
            let (alpha, beta) = decoder::angles(z_alpha, z_beta);
            let mut r = decoder::config(alpha, beta, &m);
            let mut p = decoder::momenta(z_q, &m);
            decoder::canonicalise(&mut r, &mut p, &m);
            let _ = decoder::scale_gauge(&mut r, &mut p, &m);
            decoder::to_ic(&r, &p, &m)
        }
    }
}

/// BRIEF §2.2's named regions, all Burrau: `(name, cx, cy, body)`.
pub const REGIONS: [(&str, f64, f64, usize); 8] = [
    ("near-field", 1.0, 3.0, 0),
    ("mid-field", 1.0, 6.0, 0),
    ("far", 1.0, 13.0, 0),
    ("body2 core", 1.0, -1.0, 2),
    ("body2 mid", 1.0, -5.0, 2),
    ("body1 slice", -2.0, -1.0, 1),
    ("body1 far", -2.0, -7.0, 1),
    // Pathological: drives all three bodies together. Not regularisable; expected to hit
    // the triple-collision outcome. That is correct behaviour, not a bug (BRIEF §2.6).
    ("deep interior", 0.0, 0.0, 0),
];

pub fn region(name: &str, nx: usize, ny: usize, half: f64) -> Option<Slice> {
    REGIONS
        .iter()
        .find(|r| r.0 == name)
        .map(|&(_, cx, cy, body)| Slice::body_plane(nx, ny, cx, cy, half, body))
}
