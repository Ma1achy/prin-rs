//! Tests for the shared decoder and the chart families.
//!
//! The reference names `sum p_i = 0` as the test that catches the crossed-mass swap in §0.3.
//! **It does not** — both forms sum to zero identically. That is the first thing tested here,
//! with the swap built explicitly as a negative control, because a test that cannot fire is
//! indistinguishable from a test that passes.

use prin_rs::physics::decoder::{self, Degenerate, Latent};
use prin_rs::physics::{burrau, energy, shape};
use prin_rs::Vec2;

const M: [f64; 3] = [3.0 / 12.0, 4.0 / 12.0, 5.0 / 12.0];

/// The swapped form: `m1` on `p0` and `m0` on `p1`, the mistake §0.3 warns about.
fn swapped(p_rho: Vec2<f64>, p_lam: Vec2<f64>, m: &[f64; 3]) -> [Vec2<f64>; 3] {
    let m01 = m[0] + m[1];
    [-p_rho - p_lam * (m[1] / m01), p_rho - p_lam * (m[0] / m01), p_lam]
}

// ---------------------------------------------------------------------------------------------
// The crossed-mass hazard
// ---------------------------------------------------------------------------------------------

#[test]
fn sum_p_is_zero_for_both_the_correct_form_and_the_swap() {
    // The reference calls this "the test that catches a swap". It is not: both forms give
    // `sum p = p_lam * (1 - (m0+m1)/M01) = 0`. Asserted here so the claim is on record as
    // refuted rather than quietly not relied on.
    let (pr, pl) = (Vec2::new(0.31, -0.72), Vec2::new(-0.44, 0.19));
    for p in [decoder::from_jacobi_momenta(pr, pl, &M), swapped(pr, pl, &M)] {
        let t = p[0] + p[1] + p[2];
        assert!(t.norm() < 1e-15, "sum p = {} for a form this test cannot distinguish", t.norm());
    }
}

#[test]
fn the_jacobi_round_trip_catches_the_crossed_mass_swap() {
    let (pr, pl) = (Vec2::new(0.31, -0.72), Vec2::new(-0.44, 0.19));

    let good = decoder::from_jacobi_momenta(pr, pl, &M);
    let (br, bl) = decoder::to_jacobi_momenta(&good, &M);
    assert!((br - pr).norm() < 1e-15, "the correct form does not round-trip: {}", (br - pr).norm());
    assert!((bl - pl).norm() < 1e-15);

    // The negative control. If this ever passes, the round-trip is not sensitive to the swap
    // and the arm above is decoration.
    let bad = swapped(pr, pl, &M);
    let (br2, _) = decoder::to_jacobi_momenta(&bad, &M);
    assert!(
        (br2 - pr).norm() > 1e-3,
        "the swap round-trips too ({}), so this test cannot fire",
        (br2 - pr).norm()
    );
}

#[test]
fn the_kinetic_energy_identity_catches_the_crossed_mass_swap() {
    // K == |p_rho|^2/(2 mu_rho) + |p_lam|^2/(2 mu_lam) is the physical statement that these are
    // the Jacobi momenta. Measured separation at Burrau's masses: 4.4e-16 against 2.6e-1.
    let (pr, pl) = (Vec2::new(0.31, -0.72), Vec2::new(-0.44, 0.19));
    let (mu_rho, mu_lam) = decoder::reduced(&M);
    let k_jacobi = pr.norm_sq() / (2.0 * mu_rho) + pl.norm_sq() / (2.0 * mu_lam);

    let k_of = |p: [Vec2<f64>; 3]| {
        let v = [p[0] / M[0], p[1] / M[1], p[2] / M[2]];
        energy::kinetic(&v, &M)
    };
    let good = (k_of(decoder::from_jacobi_momenta(pr, pl, &M)) - k_jacobi).abs();
    let bad = (k_of(swapped(pr, pl, &M)) - k_jacobi).abs();
    assert!(good < 1e-14, "the correct form misses the Jacobi kinetic energy by {good:e}");
    assert!(bad > 1e-3, "the swap also satisfies it ({bad:e}), so this test cannot fire");
}

#[test]
fn neither_catch_fires_at_equal_masses_and_that_is_stated() {
    // The exclusion both tests above carry. At m0 == m1 the crossed and swapped forms are the
    // same expression, so a mass-simplex chart passing through the equal-mass line has a point
    // at which neither test can distinguish them. Burrau has m0 = 3, m1 = 4, so they do fire
    // there -- but a future test written on an equal-mass configuration would be empty.
    let eq = [0.25, 0.25, 0.5];
    let (pr, pl) = (Vec2::new(0.31, -0.72), Vec2::new(-0.44, 0.19));
    let a = decoder::from_jacobi_momenta(pr, pl, &eq);
    let b = swapped(pr, pl, &eq);
    for k in 0..3 {
        assert!((a[k] - b[k]).norm() < 1e-18, "the two forms differ at equal masses");
    }
}

// ---------------------------------------------------------------------------------------------
// The decoder's invariants
// ---------------------------------------------------------------------------------------------

fn latents() -> Vec<Latent> {
    // A deterministic spread over the 8 coordinates -- no RNG, so a failure is reproducible.
    let mut out = Vec::new();
    for i in 0..40 {
        let t = |k: usize| ((i * 7 + k * 13) % 23) as f64 / 11.0 - 1.0;
        out.push(Latent {
            z_alpha: 3.0 * t(0),
            z_beta: 3.0 * t(1),
            z_q: [2.0 * t(2), 2.0 * t(3), 2.0 * t(4), 2.0 * t(5)],
            z_mu: [1.5 * t(6), 1.5 * t(7)],
        });
    }
    out
}

#[test]
fn the_decoder_conserves_total_momentum_and_puts_the_com_at_the_origin() {
    for z in latents() {
        let d = decoder::decode(&z);
        assert!(d.flag.is_none(), "unexpected degenerate: {:?}", d.flag);
        let p: Vec<Vec2<f64>> =
            (0..3).map(|k| d.ic.s.v[k] * d.ic.m[k]).collect();
        let tot = p[0] + p[1] + p[2];
        assert!(tot.norm() < 1e-14, "sum p = {:e}", tot.norm());
        let c = energy::com(&d.ic.s.r, &d.ic.m);
        assert!(c.norm() < 1e-14, "com at {:e}", c.norm());
    }
}

#[test]
fn the_masses_are_a_normalised_simplex_point() {
    for z in latents() {
        let (m, flag) = decoder::masses(z.z_mu);
        assert!(flag.is_none());
        let s: f64 = m.iter().sum();
        assert!((s - 1.0).abs() < 1e-15, "masses sum to {s}");
        assert!(m.iter().all(|&x| x > 0.0), "a mass is non-positive: {m:?}");
    }
}

#[test]
fn the_canonical_frame_decode_gives_unit_inertia() {
    // I = 1 is what the scale gauge enforces, and what makes a latent slice comparable across
    // its own extent rather than a picture of the overall size varying.
    for z in latents() {
        let d = decoder::decode(&z);
        let i = shape::inertia(&d.ic.s.r, &d.ic.m);
        assert!((i - 1.0).abs() < 1e-12, "I = {i}, want 1");
    }
}

#[test]
fn small_alpha_is_a_wide_inner_pair_not_a_tight_one() {
    // The reference's "easy to get backwards" note, asserted as a DIRECTION rather than by
    // re-checking the formula against itself. Getting it inverted would put every hierarchical
    // configuration where an anti-hierarchical one belongs, and every picture would still look
    // like physics.
    let m = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
    let sep = |alpha: f64| {
        let r = decoder::config(alpha, 1.0, &m);
        (r[1] - r[0]).norm()
    };
    let small = sep(0.1);
    let large = sep(std::f64::consts::FRAC_PI_2 - 0.1);
    assert!(
        small > large,
        "small alpha gave inner separation {small} against {large} at large alpha -- inverted"
    );
}

#[test]
fn the_canonicaliser_is_idempotent_and_a_no_op_on_an_already_canonical_state() {
    for z in latents() {
        let d = decoder::decode(&z);
        let mut r = d.ic.s.r;
        let mut p: [Vec2<f64>; 3] =
            [d.ic.s.v[0] * d.ic.m[0], d.ic.s.v[1] * d.ic.m[1], d.ic.s.v[2] * d.ic.m[2]];
        let (r0, p0) = (r, p);
        decoder::canonicalise(&mut r, &mut p, &d.ic.m);
        for k in 0..3 {
            assert!(
                (r[k] - r0[k]).norm() < 1e-13 && (p[k] - p0[k]).norm() < 1e-13,
                "canonicalise moved an already-canonical state"
            );
        }
    }
}

#[test]
fn the_scale_gauge_uses_the_asymmetric_powers() {
    // Positions divide by l, momenta multiply by sqrt(l). Equal powers would leave every
    // configuration looking right and the Hamiltonian's scaling wrong.
    let m = [0.25, 0.35, 0.40];
    let mut r = [Vec2::new(2.0, 0.0), Vec2::new(-1.0, 1.0), Vec2::new(0.5, -2.0)];
    let mut p = [Vec2::new(0.1, 0.2), Vec2::new(-0.3, 0.05), Vec2::new(0.2, -0.25)];
    let (r0, p0) = (r, p);
    let i = shape::inertia(&r, &m);
    let l = i.sqrt();
    assert!(decoder::scale_gauge(&mut r, &mut p, &m).is_none());
    for k in 0..3 {
        assert!((r[k] - r0[k] / l).norm() < 1e-14, "positions did not divide by l");
        assert!((p[k] - p0[k] * l.sqrt()).norm() < 1e-14, "momenta did not multiply by sqrt(l)");
    }
    // And the control: l is not 1, or the test is vacuous.
    assert!((l - 1.0).abs() > 0.1, "l = {l}: the gauge is a no-op here and this proves nothing");
}

// ---------------------------------------------------------------------------------------------
// 2.2 The deterministic momentum construction
// ---------------------------------------------------------------------------------------------

#[test]
fn the_momentum_construction_hits_lz_and_k_exactly() {
    // Three independent constraints; the construction should hit all three to machine
    // precision, including near the parabola boundary where K* -> K_min and the mix
    // coefficient `a` goes to zero.
    let m = [3.0 / 12.0, 4.0 / 12.0, 5.0 / 12.0];
    let r = decoder::config(0.7, 1.1, &m);
    let i: f64 = (0..3).map(|k| m[k] * r[k].norm_sq()).sum();

    for &lz in &[-0.7f64, -0.05, 0.0, 0.31, 1.2] {
        let k_min = lz * lz / (2.0 * i);
        for &over in &[0.0f64, 1e-12, 1e-6, 0.05, 2.0] {
            let k_star = k_min + over;
            let p = decoder::momenta_for(lz, k_star, &r, &m).expect("should be feasible");
            let tot = p[0] + p[1] + p[2];
            assert!(tot.norm() < 1e-13, "sum p = {:e} at lz={lz} over={over}", tot.norm());

            let got_lz: f64 = (0..3).map(|k| r[k].x * p[k].y - r[k].y * p[k].x).sum();
            assert!(
                (got_lz - lz).abs() < 1e-12,
                "Lz = {got_lz} want {lz} (over = {over})"
            );

            let v = [p[0] / m[0], p[1] / m[1], p[2] / m[2]];
            let got_k = energy::kinetic(&v, &m);
            assert!(
                (got_k - k_star).abs() < 1e-12 * k_star.max(1.0),
                "K = {got_k} want {k_star} (lz = {lz})"
            );
        }
    }
}

#[test]
fn a_target_below_k_min_is_refused_rather_than_approximated() {
    let m = [3.0 / 12.0, 4.0 / 12.0, 5.0 / 12.0];
    let r = decoder::config(0.7, 1.1, &m);
    let i: f64 = (0..3).map(|k| m[k] * r[k].norm_sq()).sum();
    let lz = 1.0;
    let k_min = lz * lz / (2.0 * i);
    assert_eq!(decoder::momenta_for(lz, k_min * 0.5, &r, &m), Err(Degenerate::BelowKMin));
    // The control: just above K_min is accepted.
    assert!(decoder::momenta_for(lz, k_min * 1.000001, &r, &m).is_ok());
}

#[test]
fn a_degenerate_primary_seed_falls_back_rather_than_giving_up() {
    // Force rho = 0 -- bodies 0 and 1 coincident -- so the primary seed `(rho, 0)` has zero
    // norm. A fallback must be chosen, not `SeedsExhausted`.
    let m = [3.0 / 12.0, 4.0 / 12.0, 5.0 / 12.0];
    let r = [Vec2::new(0.2, 0.3), Vec2::new(0.2, 0.3), Vec2::new(-0.4, -0.6)];
    let i: f64 = (0..3).map(|k| m[k] * r[k].norm_sq()).sum();
    assert!(i > 0.0, "the test configuration must have nonzero inertia");
    let p = decoder::momenta_for(0.2, 0.2 * 0.2 / (2.0 * i) + 0.5, &r, &m)
        .expect("a fallback seed should have been chosen");
    let tot = p[0] + p[1] + p[2];
    // This configuration is deliberately NOT centred on its COM, which is how it caught that
    // `momenta_for` returned a drifting system for an off-centre input: the rigid-rotation step
    // is `v = omega J r`, whose total momentum is `omega J (M R_com)`.
    assert!(tot.norm() < 1e-13, "sum p = {:e} on an off-centre configuration", tot.norm());
}

// ---------------------------------------------------------------------------------------------
// 4 The Burrau family
// ---------------------------------------------------------------------------------------------

#[test]
fn burrau_at_nu_one_half_is_the_three_four_five_triangle() {
    let (m, r) = decoder::burrau_family(decoder::nu_of(2, 1));
    // (m,n) = (2,1) -> a=3, b=4, c=5, masses (c,b,a)/12 = (5,4,3)/12.
    for (got, want) in m.iter().zip([5.0 / 12.0, 4.0 / 12.0, 3.0 / 12.0].iter()) {
        assert!((got - want).abs() < 1e-15, "masses {m:?}");
    }
    // Sides: |r1-r0| = a/c = 3/5, |r2-r0| = b/c = 4/5, |r1-r2| = 1.
    let d = prin_rs::physics::newton::pair_dists(&r);
    assert!((d[0] - 0.6).abs() < 1e-15, "|r1-r0| = {}", d[0]);
    assert!((d[1] - 0.8).abs() < 1e-15, "|r2-r0| = {}", d[1]);
    assert!((d[2] - 1.0).abs() < 1e-15, "|r2-r1| = {}", d[2]);
}

#[test]
fn the_burrau_family_reproduces_the_repos_configuration_up_to_the_gauge() {
    // The convention clash, pinned rather than assumed. The repo uses MASSES = [3,4,5] at
    // [(1,3), (-2,-1), (1,-1)]; the reference gives (5,4,3)/12 at [(0,0), (3/5,0), (0,4/5)].
    // Both are "mass equals opposite side", so they are the same system up to the scale gauge
    // and a body relabelling -- which is exactly the kind of thing that is obvious until it is
    // wrong. `tests/burrau_constants.rs` pins M = 12, R = 2.2361, E = -12.8167; the invariant
    // that survives both the gauge and the relabelling is the SHAPE.
    let (mf, rf) = decoder::burrau_family(decoder::nu_of(2, 1));

    let repo_r = {
        let s = burrau::state::<f64>();
        s.r
    };
    let repo_m = burrau::MASSES;

    // The repo's body k corresponds to the reference's body `perm[k]`. Recover it from the
    // side lengths rather than asserting it: mass equals the opposite side in both, so the
    // ordering of the masses fixes the ordering of the bodies.
    let mut perm = [0usize; 3];
    for k in 0..3 {
        let target = repo_m[k] / 12.0;
        perm[k] = (0..3)
            .min_by(|&a, &b| {
                (mf[a] - target).abs().partial_cmp(&(mf[b] - target).abs()).unwrap()
            })
            .unwrap();
    }
    assert_eq!(
        {
            let mut p = perm;
            p.sort();
            p
        },
        [0, 1, 2],
        "the mass correspondence is not a permutation: {perm:?}"
    );

    // Side-length ratios are invariant to the scale gauge, so they are what compares.
    let ratio = |r: &[Vec2<f64>; 3], p: [usize; 3]| {
        let d = [
            (r[p[1]] - r[p[0]]).norm(),
            (r[p[2]] - r[p[0]]).norm(),
            (r[p[2]] - r[p[1]]).norm(),
        ];
        [d[0] / d[2], d[1] / d[2]]
    };
    let a = ratio(&repo_r, [0, 1, 2]);
    let b = ratio(&rf, perm);
    for k in 0..2 {
        assert!(
            (a[k] - b[k]).abs() < 1e-12,
            "side ratio {k}: repo {} vs reference {} under permutation {perm:?}",
            a[k],
            b[k]
        );
    }
}

#[test]
fn the_burrau_family_is_smooth_in_nu_and_the_triples_sit_on_it() {
    // The bifurcation strip's premise: the primitive triples are a countable set of points on a
    // continuous curve, so `nu` sweeps between them.
    for (mm, n, sides) in [(2u32, 1u32, [3.0, 4.0, 5.0]), (3, 2, [5.0, 12.0, 13.0]), (4, 3, [7.0, 24.0, 25.0])] {
        let (_, r) = decoder::burrau_family(decoder::nu_of(mm, n));
        let d = prin_rs::physics::newton::pair_dists(&r);
        let want = [sides[0] / sides[2], sides[1] / sides[2], 1.0];
        for k in 0..3 {
            assert!(
                (d[k] - want[k]).abs() < 1e-14,
                "({mm},{n}) side {k}: {} want {}",
                d[k],
                want[k]
            );
        }
    }
    // And smoothness: a small step in nu moves the configuration a small amount.
    let step = |nu: f64| {
        let (_, a) = decoder::burrau_family(nu);
        let (_, b) = decoder::burrau_family(nu + 1e-6);
        (0..3).map(|k| (a[k] - b[k]).norm()).fold(0.0, f64::max)
    };
    for nu in [0.1, 0.3, 0.5, 0.7, 0.9] {
        assert!(step(nu) < 1e-5, "nu = {nu} is not smooth: step {}", step(nu));
    }
}

// ---------------------------------------------------------------------------------------------
// The chart families
// ---------------------------------------------------------------------------------------------

use prin_rs::ensemble::jitter::{self, Scheme};
use prin_rs::grid::{self, Chart, Domain, Slice};

fn base_latent() -> Latent {
    Latent { z_alpha: 0.3, z_beta: -0.4, z_q: [0.2, -0.1, 0.35, 0.05], z_mu: [0.15, -0.25] }
}

#[test]
fn body_plane_still_decodes_to_burrau_masses_and_every_old_chart_does_too() {
    // The guard that made the `Ic` change safe. If a pre-existing chart ever returned anything
    // else, every result in the repo would be about a different three-body system while looking
    // identical.
    for chart in [
        Chart::BodyPlane,
        Chart::plane_for_body(0),
        Chart::plane_for_body(2),
        Chart::shape_at_burrau(0.0),
        Chart::shape_at_burrau(1.3),
    ] {
        for (u, v) in [(0.0, 0.0), (1.0, 3.0), (-0.4, 0.7)] {
            let ic = grid::decode_state(&chart, 0, u, v);
            assert_eq!(ic.m, burrau::MASSES, "chart {} moved the masses", chart.name());
        }
    }
}

#[test]
fn an_axis_aligned_latent_slice_equals_a_direct_sweep_of_those_two_coordinates() {
    // `Chart::latent_axes(z0, i, j)` must be exactly "vary coordinates i and j". If it were not,
    // the named planes in the reference's table would be describing something else.
    let z0 = base_latent();
    for (i, j) in [(0usize, 1usize), (2, 3), (6, 7), (0, 4)] {
        let chart = Chart::latent_axes(z0, i, j);
        for (u, v) in [(0.0, 0.0), (0.4, -0.9), (-1.2, 0.3)] {
            let got = grid::decode_state(&chart, 0, u, v);
            let mut z = z0;
            z.set(i, z0.get(i) + u);
            z.set(j, z0.get(j) + v);
            let want = decoder::decode(&z).ic;
            for k in 0..3 {
                assert!(
                    (got.s.r[k] - want.s.r[k]).norm() < 1e-15,
                    "latent_axes({i},{j}) at ({u},{v}) body {k}"
                );
                assert!((got.m[k] - want.m[k]).abs() < 1e-15);
            }
        }
    }
}

#[test]
fn an_oblique_latent_basis_is_orthonormal() {
    let a = [0.3, -1.2, 0.5, 0.9, -0.4, 0.1, 0.7, -0.6];
    let b = [1.1, 0.2, -0.8, 0.3, 0.6, -0.9, 0.15, 0.4];
    let Chart::Latent { q1, q2, .. } = Chart::latent_oblique(base_latent(), a, b) else {
        panic!("not a latent chart")
    };
    let dot = |x: &[f64; 8], y: &[f64; 8]| (0..8).map(|k| x[k] * y[k]).sum::<f64>();
    assert!((dot(&q1, &q1) - 1.0).abs() < 1e-14);
    assert!((dot(&q2, &q2) - 1.0).abs() < 1e-14);
    assert!(dot(&q1, &q2).abs() < 1e-14, "basis is not orthogonal: {}", dot(&q1, &q2));
    // The control: the seeds were NOT already orthonormal, so Gram-Schmidt did work.
    assert!((dot(&a, &b)).abs() > 0.1, "the seeds were already orthogonal; this proves nothing");
}

#[test]
fn the_mass_simplex_varies_the_masses_and_keeps_them_positive() {
    let chart = Chart::MassSimplex {
        z_alpha: 0.2,
        z_beta: -0.3,
        z_q: [0.1, 0.2, -0.1, 0.05],
        margin: 0.02,
    };
    let mut seen: Vec<[f64; 3]> = Vec::new();
    for i in 0..=10 {
        for j in 0..=10 {
            let (u, v) = (i as f64 / 10.0, j as f64 / 10.0);
            let ic = grid::decode_state(&chart, 0, u, v);
            let s: f64 = ic.m.iter().sum();
            assert!((s - 1.0).abs() < 1e-14, "masses sum to {s} at ({u},{v})");
            assert!(
                ic.m.iter().all(|&x| x >= 0.02 - 1e-15),
                "a mass fell below the margin at ({u},{v}): {:?}",
                ic.m
            );
            seen.push(ic.m);
        }
    }
    // The point of the chart: the masses must actually move. A simplex chart that returned one
    // mass everywhere would be an expensive way to draw a constant.
    let m0: Vec<f64> = seen.iter().map(|m| m[0]).collect();
    let span = m0.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - m0.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(span > 0.5, "m0 spans only {span} across the simplex");
}

#[test]
fn the_burrau_chart_puts_the_primitive_triples_on_a_continuous_curve() {
    let chart = Chart::BurrauFamily { nu_lo: 0.05, nu_hi: 0.95, k_max: 1.0, gamma_k: 1.5 };
    // nu = 1/2 is the (3,4,5) triangle; find the u that reaches it.
    let u = (0.5 - 0.05) / (0.95 - 0.05);
    let ic = grid::decode_state(&chart, 0, u, 0.0);
    // The scale gauge has been applied, so compare side RATIOS.
    let d = prin_rs::physics::newton::pair_dists(&ic.s.r);
    let mut sorted = d;
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (a, b, c) = (sorted[0], sorted[1], sorted[2]);
    assert!((a / c - 0.6).abs() < 1e-12, "shortest/longest = {}", a / c);
    assert!((b / c - 0.8).abs() < 1e-12, "middle/longest = {}", b / c);

    // v = 0 is the rest start; v > 0 adds kinetic energy and nothing else.
    let rest = grid::decode_state(&chart, 0, u, 0.0);
    assert!(
        rest.s.v.iter().all(|x| x.norm() < 1e-14),
        "v = 0 is not a rest start: {:?}",
        rest.s.v
    );
    let hot = grid::decode_state(&chart, 0, u, 0.8);
    let k = energy::kinetic(&hot.s.v, &hot.m);
    assert!(k > 0.0, "v = 0.8 added no kinetic energy");
}

#[test]
fn the_invariant_chart_is_feasible_at_every_pixel_of_the_unit_square() {
    // The warp's whole purpose: no `(u,v)` in [0,1]^2 produces K* < K_min. If any did, the
    // chart would need a clamp, and a clamped region is a flat region that looks like physics.
    let chart = Chart::Invariant {
        base: base_latent(),
        k_max: 2.0,
        gamma_k: 1.5,
        report_e: false,
    };
    let mut lz_span = (f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..=20 {
        for j in 0..=20 {
            let (u, v) = (i as f64 / 20.0, j as f64 / 20.0);
            let ic = grid::decode_state(&chart, 0, u, v);
            let p: Vec<Vec2<f64>> = (0..3).map(|k| ic.s.v[k] * ic.m[k]).collect();
            let tot = p[0] + p[1] + p[2];
            assert!(tot.norm() < 1e-12, "sum p = {:e} at ({u},{v})", tot.norm());
            let lz: f64 = (0..3).map(|k| ic.s.r[k].x * p[k].y - ic.s.r[k].y * p[k].x).sum();
            assert!(lz.is_finite(), "Lz is not finite at ({u},{v})");
            lz_span = (lz_span.0.min(lz), lz_span.1.max(lz));
        }
    }
    // And the axis is used: a chart whose Lz never moved would be feasible trivially.
    assert!(
        lz_span.1 - lz_span.0 > 0.5,
        "Lz spans only {} across the chart",
        lz_span.1 - lz_span.0
    );
}

#[test]
fn lz_e_and_lz_k_are_the_same_map_which_the_reference_lists_as_two_charts() {
    // The reference's §2 gives (Lz,E) and (Lz,K) as separate charts, #4 in its implementation
    // order, "most machinery, most verification". But its own construction parameterises BOTH
    // by K(t) = K_max t^gamma; for (Lz,E) it then reports E = U + K(t), which is a relabelling
    // of the axis rather than a different sweep. Asserted bitwise, so the claim is measured.
    let a = Chart::Invariant { base: base_latent(), k_max: 2.0, gamma_k: 1.5, report_e: true };
    let b = Chart::Invariant { base: base_latent(), k_max: 2.0, gamma_k: 1.5, report_e: false };
    for i in 0..=7 {
        for j in 0..=7 {
            let (u, v) = (i as f64 / 7.0, j as f64 / 7.0);
            let (x, y) = (grid::decode_state(&a, 0, u, v), grid::decode_state(&b, 0, u, v));
            for k in 0..3 {
                assert_eq!(x.s.r[k].x.to_bits(), y.s.r[k].x.to_bits());
                assert_eq!(x.s.v[k].y.to_bits(), y.s.v[k].y.to_bits());
            }
        }
    }
    // What DOES make them differ is gamma_k, and only that.
    let c = Chart::Invariant { base: base_latent(), k_max: 2.0, gamma_k: 1.0, report_e: true };
    let (x, y) = (grid::decode_state(&a, 0, 0.6, 0.4), grid::decode_state(&c, 0, 0.6, 0.4));
    let d: f64 = (0..3).map(|k| (x.s.v[k] - y.s.v[k]).norm()).sum();
    assert!(d > 1e-6, "gamma_k did not change the chart either; then it has one degree of freedom fewer than stated");

    // The names still differ, so a dump records which axis was intended.
    assert_eq!(a.name(), "invariant_lz_e");
    assert_eq!(b.name(), "invariant_lz_k");
}

#[test]
fn energy_normalisation_is_refused_on_a_chart_whose_axis_is_energy() {
    let inv = Chart::Invariant { base: base_latent(), k_max: 1.0, gamma_k: 1.5, report_e: true };
    assert!(inv.forbids_energy_normalisation());
    assert!(inv.validate(0.5, 0.5, 0.5, 0.4).is_err(), "E* != 0 was accepted");
    // The controls: E* = 0 is fine, and a chart that does not carry energy accepts any E*.
    assert!(inv.validate(0.0, 0.5, 0.5, 0.4).is_ok());
    assert!(!Chart::BodyPlane.forbids_energy_normalisation());
    assert!(Chart::BodyPlane.validate(-12.8, 1.0, 3.0, 0.05).is_ok());
}

#[test]
fn a_unit_domain_chart_refuses_a_slice_box_that_leaves_the_square() {
    let chart = Chart::BurrauFamily { nu_lo: 0.05, nu_hi: 0.95, k_max: 1.0, gamma_k: 1.5 };
    assert_eq!(chart.domain(), Domain::Unit);
    assert!(chart.validate(0.0, 0.5, 0.5, 0.4).is_ok());
    assert!(chart.validate(0.0, 0.5, 0.5, 0.6).is_err(), "a box leaving [0,1]^2 was accepted");
    assert!(chart.validate(0.0, 0.05, 0.5, 0.1).is_err());
    // A free-domain chart has no such constraint.
    assert_eq!(Chart::BodyPlane.domain(), Domain::Free);
    assert!(Chart::BodyPlane.validate(0.0, 1.0, 3.0, 5.0).is_ok());
}

#[test]
fn jitter_reflects_into_the_unit_square_and_does_not_collapse_the_copies() {
    // Clamping would put several copies on the boundary, and a collapsed decode gives
    // `ensemble_spread` exactly zero -- which reads as "perfectly resolved" and stops the
    // descent with a small tidy tree built from nothing.
    let chart = Chart::BurrauFamily { nu_lo: 0.05, nu_hi: 0.95, k_max: 1.0, gamma_k: 1.5 };
    // An EDGE cell: half a cell width from 0, so at jitter_frac = 0.5 the copies straddle it.
    let s = Slice::body_plane(4, 4, 0.0625, 0.5, 0.0625, 0).with_chart(chart);
    let c = jitter::copies_with::<f64>(&s, 0, 15, 0.5, 0, Scheme::Halton);
    assert_eq!(c.len(), 16);

    let ics: Vec<prin_rs::physics::Cart<f64>> = c.iter().map(|x| x.s).collect();
    let distinct = prin_rs::decode::distinct(&ics);
    assert_eq!(
        distinct,
        16,
        "the copies collapsed to {distinct} distinct ICs at an edge cell -- that reads as \
         perfectly resolved and is not"
    );
}

#[test]
fn every_chart_records_enough_to_be_reproduced() {
    // `Chart::name()` alone is not enough: a Plane's basis and a Latent's (z0, q1, q2) are free,
    // so two dumps with one name can be different configurations.
    let charts = [
        Chart::plane_for_body(1),
        Chart::shape_at_burrau(0.4),
        Chart::latent_axes(base_latent(), 0, 1),
        Chart::BurrauFamily { nu_lo: 0.05, nu_hi: 0.95, k_max: 1.0, gamma_k: 1.5 },
        Chart::Invariant { base: base_latent(), k_max: 2.0, gamma_k: 1.5, report_e: false },
        Chart::MassSimplex { z_alpha: 0.2, z_beta: -0.3, z_q: [0.1; 4], margin: 0.02 },
    ];
    for c in charts {
        let p = c.params();
        assert!(p.len() > 8, "chart {} records nothing reproducible: {p:?}", c.name());
    }
    // BodyPlane genuinely has no parameters, and says so rather than inventing some.
    assert_eq!(Chart::BodyPlane.params(), "-");
}

#[test]
fn only_the_affine_charts_claim_to_be_affine() {
    // `is_affine` gates the linearised decoder's exactness. Claiming a nonlinear chart is affine
    // would report a curvature term as structurally zero when it is merely unmeasured -- the
    // "a difference can be small because both sides are dead" failure at chart level.
    assert!(Chart::BodyPlane.is_affine());
    assert!(Chart::plane_for_body(0).is_affine());
    assert!(!Chart::shape_at_burrau(0.0).is_affine());
    assert!(!Chart::latent_axes(base_latent(), 0, 1).is_affine());
    assert!(!Chart::BurrauFamily { nu_lo: 0.1, nu_hi: 0.9, k_max: 1.0, gamma_k: 1.5 }.is_affine());
    assert!(!Chart::MassSimplex { z_alpha: 0.0, z_beta: 0.0, z_q: [0.0; 4], margin: 0.02 }
        .is_affine());
}

// ---------------------------------------------------------------------------------------------
// The GLSL reference's pinned decode, and its four default presets
//
// `Ma1achy/principia-ii`, `src/shaders/principia/frag.glsl:19-59` and `src/state.ts:71-76`.
// ---------------------------------------------------------------------------------------------

/// **The landmark.** `z = 0` decodes to the equilateral Lagrange configuration. That is a named
/// physical configuration, which makes it a stronger check on the reconstruction algebra than any
/// invariant — it can be verified by eye in the rendered image.
///
/// **What it cannot see, stated so it is not mistaken for a fuller check than it is:** at `z = 0`
/// the momentum coordinates and the mass logits are all zero, so `Q_MAX`, `MU_MAX` and the choice
/// between `tanh(z)` and `2*sigmoid(z)-1` every one of them drops out of the arithmetic. This test
/// passes unchanged under all three of the constants this port corrected.
/// `the_pinned_saturation_constants_are_the_glsls_not_the_latex_reference` covers those.
#[test]
fn the_origin_of_the_latent_chart_is_the_equilateral_lagrange_configuration() {
    let d = decoder::decode(&Latent::default());
    assert_eq!(d.flag, None, "the origin should decode cleanly");
    let ic = d.ic;

    for k in 0..3 {
        assert!((ic.m[k] - 1.0 / 3.0).abs() < 1e-15, "mass {k} = {}", ic.m[k]);
    }

    let want = [Vec2::new(-0.8660254037844386, -0.5), Vec2::new(0.8660254037844386, -0.5),
                Vec2::new(0.0, 1.0)];
    let mut worst = 0f64;
    for k in 0..3 {
        worst = worst.max((ic.s.r[k] - want[k]).norm());
    }
    assert!(worst < 1e-14, "positions differ from Lagrange by {worst:e}");

    // Equilateral: all three separations equal. This is the part a person can check in the image.
    let d01 = (ic.s.r[0] - ic.s.r[1]).norm();
    let d02 = (ic.s.r[0] - ic.s.r[2]).norm();
    let d12 = (ic.s.r[1] - ic.s.r[2]).norm();
    println!("Lagrange separations: {d01:.15} {d02:.15} {d12:.15}  (sqrt 3 = {:.15})", 3f64.sqrt());
    assert!((d01 - d02).abs() < 1e-14 && (d01 - d12).abs() < 1e-14, "not equilateral");
    assert!((d01 - 3f64.sqrt()).abs() < 1e-14);

    // Released from rest at the origin: every momentum coordinate saturates to zero.
    for k in 0..3 {
        assert!(ic.s.v[k].norm() < 1e-15, "body {k} is not at rest: {}", ic.s.v[k].norm());
    }

    // `I = 1` and `COM = 0` are ALGEBRAIC IDENTITIES of the canonical-frame decode
    // (`I = cos^2 a + sin^2 a`; `m0 r0 + m1 r1 = -M01 m2 lam` cancels `m2 r2`), so they hold under
    // any mass factors at all and cannot fail from a physics error. Kept as wiring guards, and
    // labelled as such rather than quoted as evidence.
    let i = shape::inertia(&ic.s.r, &ic.m);
    assert!((i - 1.0).abs() < 1e-15, "I = {i}");
    let com = ic.s.r[0] * ic.m[0] + ic.s.r[1] * ic.m[1] + ic.s.r[2] * ic.m[2];
    assert!(com.norm() < 1e-15, "COM at {com:?}");
}

/// The three constants the landmark is blind to, pinned against the GLSL so reverting any of them
/// fails here. Values recomputed from `frag.glsl:21-22, 35-36, 53-54`, not copied from the tree.
#[test]
fn the_pinned_saturation_constants_are_the_glsls_not_the_latex_reference() {
    assert_eq!(decoder::MU_MAX, 5.0, "frag.glsl:21");
    assert_eq!(decoder::Q_MAX, 2.0, "frag.glsl:22");

    let sig = |z: f64| 1.0 / (1.0 + (-z).exp());

    // Masses: `MU_MAX*(2*sigmoid(z)-1)`, which is `MU_MAX*tanh(z/2)` -- HALF the LaTeX
    // reference's `mu_max*tanh(z)`. At z = (1.0, -0.5) the two forms are far apart, so the
    // negative control below is not a rounding argument.
    let z_mu = [1.0f64, -0.5];
    let (m, flag) = decoder::masses(z_mu);
    assert_eq!(flag, None);
    let softmax = |l: [f64; 3]| {
        let e = l.map(f64::exp);
        let s: f64 = e.iter().sum();
        [e[0] / s, e[1] / s, e[2] / s]
    };
    let want = {
        let mu: Vec<f64> = z_mu.iter().map(|&z| 5.0 * (2.0 * sig(z) - 1.0)).collect();
        softmax([0.0, mu[0], mu[1]])
    };
    for k in 0..3 {
        assert!((m[k] - want[k]).abs() < 1e-15, "mass {k}: {} vs {}", m[k], want[k]);
    }
    // The negative control: the LaTeX form would give visibly different masses here.
    let latex = {
        let mu: Vec<f64> = z_mu.iter().map(|&z| 5.0 * z.tanh()).collect();
        softmax([0.0, mu[0], mu[1]])
    };
    let gap = (0..3).map(|k| (m[k] - latex[k]).abs()).fold(0.0, f64::max);
    println!("half-gain vs tanh masses at z_mu = {z_mu:?}: worst |dm| = {gap:.4}");
    assert!(gap > 0.05, "the two saturation forms agree here; this test cannot fire");

    // Momenta: `Q_MAX*(2*sigmoid(z)-1)` per component, and `p2 == p_lambda` exactly.
    let z_q = [0.7f64, -1.1, 0.3, 0.9];
    let p = decoder::momenta(z_q, &M);
    let pl = Vec2::new(2.0 * (2.0 * sig(z_q[2]) - 1.0), 2.0 * (2.0 * sig(z_q[3]) - 1.0));
    assert!((p[2] - pl).norm() < 1e-15, "p2 should be p_lambda: {:?} vs {pl:?}", p[2]);
    // The control: the assertion above is sensitive to Q_MAX only if the old value gives a
    // measurably different answer here. It does -- the gain is linear, so Q_MAX = 1 halves it.
    let old = pl * 0.5;
    println!("Q_MAX 2 vs 1 at z_q = {z_q:?}: |dp_lambda| = {:.4}", (pl - old).norm());
    assert!((pl - old).norm() > 0.1, "Q_MAX = 1 would give the same p_lambda; this cannot fire");
}

/// **The gauge is inert on the latent chart, and that is a claim worth checking.**
///
/// `decode` applies `canonicalise` and `scale_gauge`; the GLSL's `decodeIC` applies neither. They
/// should be no-ops here: `rho~ = (cos a, 0)` already sits on `+x` so the rotation angle is zero,
/// and `lam~_y = sin a sin b >= 0` for `b in [0, pi]` so the mirror never fires. But `I = 1` only
/// in exact algebra — `scale_gauge` divides by `sqrt(1 +- eps)` — so this is agreement to ~1e-15
/// with the residual printed, **not** a bitwise claim.
#[test]
fn the_canonicaliser_and_scale_gauge_are_no_ops_on_the_latent_chart() {
    let mut worst_r = 0f64;
    let mut worst_p = 0f64;
    for z in [
        Latent::default(),
        base_latent(),
        Latent { z_alpha: -2.1, z_beta: 1.7, z_q: [0.9, -1.4, 0.2, 1.1], z_mu: [-0.8, 1.3] },
        Latent { z_alpha: 3.0, z_beta: -3.0, z_q: [-1.5, 1.5, -1.5, 1.5], z_mu: [2.0, -2.0] },
    ] {
        // The raw GLSL path: decode with no gauge applied at all.
        let (m, _) = decoder::masses(z.z_mu);
        let (a, b) = decoder::angles(z.z_alpha, z.z_beta);
        let raw_r = decoder::config(a, b, &m);
        let raw_p = decoder::momenta(z.z_q, &m);

        let got = decoder::decode(&z).ic;
        for k in 0..3 {
            worst_r = worst_r.max((got.s.r[k] - raw_r[k]).norm());
            worst_p = worst_p.max((got.s.v[k] * got.m[k] - raw_p[k]).norm());
        }
    }
    println!("gauge residual over 4 latent points: positions {worst_r:.3e}, momenta {worst_p:.3e}");
    assert!(worst_r < 1e-14, "canonicaliser/scale_gauge moved positions by {worst_r:e}");
    assert!(worst_p < 1e-14, "canonicaliser/scale_gauge moved momenta by {worst_p:e}");
}

/// The reference's `shape_pl` basis is **un-normalised** — each direction has norm `sqrt(2)`.
/// Routing it through `latent_oblique` would orthonormalise it and render a different slice while
/// looking like a tidy-up. Pinned so that fails here instead.
#[test]
fn the_shape_pl_preset_basis_is_not_orthonormalised() {
    let Chart::Latent { z0, q1, q2 } = Chart::preset_shape_pl() else { panic!() };
    let dot = |x: &[f64; 8], y: &[f64; 8]| (0..8).map(|k| x[k] * y[k]).sum::<f64>();
    assert!((dot(&q1, &q1) - 2.0).abs() < 1e-15, "|q1|^2 should be 2, not 1");
    assert!((dot(&q2, &q2) - 2.0).abs() < 1e-15, "|q2|^2 should be 2, not 1");
    assert!(dot(&q1, &q2).abs() < 1e-15, "the two directions are already orthogonal");
    assert_eq!(z0, Latent::default(), "the presets sit at the origin");

    // And what `latent_oblique` would have done instead: scaled both by 1/sqrt(2), which halves
    // the extent of the slice at a fixed camera half-width.
    let Chart::Latent { q1: o1, .. } = Chart::latent_oblique(Latent::default(), q1, q2) else {
        panic!()
    };
    assert!((dot(&o1, &o1) - 1.0).abs() < 1e-14, "latent_oblique should normalise");
}

/// **The alpha-varying direction carries `pLambda.y`, not `pLambda.x`.**
///
/// The GLSL preset is `q1 = e0 + e6`, `q2 = e1 + e7`, and it pairs **by slot**: in its own
/// indexing that is `beta` with `pLambda.x` and `alpha` with `pLambda.y`. This module renumbers
/// alpha and beta into the spec's order and must carry their momentum partners with them. The
/// port did not, and paired alpha with `pLambda.x`.
///
/// `shape_pl` is the **only** preset with a cross-coupling — the other three are pure-config or
/// pure-momentum — so it is the only one that can fail this way, and it rendered as *twisted*
/// rather than tilted.
///
/// **The index assertions alone would not be a test.** They would pass on a "fix" that transposes
/// `q1` and `q2`, and the whole finding is that transposing does not work: that gives
/// `e_beta + e_pLy`, `e_alpha + e_pLx`, still crossed. It is a genuinely different 2-plane
/// through the 8D space, not a reorientation of the same one. The second arm below is what makes
/// this able to fire — it decodes the crossed form *and its transposition* and asserts both
/// disagree with the correct plane far outside rounding.
#[test]
fn the_shape_pl_preset_pairs_alpha_with_p_lambda_y() {
    let Chart::Latent { q1, q2, .. } = Chart::preset_shape_pl() else { panic!() };
    assert_eq!(q1[0], 1.0, "q1 must vary alpha");
    assert_eq!(q1[5], 1.0, "q1 must carry p_lambda.y");
    assert_eq!(q1[4], 0.0, "q1 must NOT carry p_lambda.x");
    assert_eq!(q2[1], 1.0, "q2 must vary beta");
    assert_eq!(q2[4], 1.0, "q2 must carry p_lambda.x");
    assert_eq!(q2[5], 0.0, "q2 must NOT carry p_lambda.y");

    // The crossed form that shipped, and its transposition -- the "fix" that is not one.
    let crossed = Chart::Latent {
        z0: Latent::default(),
        q1: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        q2: [0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    };
    let Chart::Latent { q1: c1, q2: c2, .. } = crossed else { panic!() };
    let transposed = Chart::Latent { z0: Latent::default(), q1: c2, q2: c1 };

    // Off-axis on purpose: on either axis one of the two coupled coordinates is zero and the
    // three planes coincide, so a test sampled there could not fire.
    let (u, v) = (0.7f64, -1.3f64);
    let good = grid::decode_state(&Chart::preset_shape_pl(), 0, u, v);
    for (label, other) in [("crossed", crossed), ("crossed, transposed", transposed)] {
        let got = grid::decode_state(&other, 0, u, v);
        let d = prin_rs::decode::max_abs_diff(&good.s, &got.s);
        println!("shape_pl vs {label:>19} at ({u},{v}): max |dIC| = {d:.4e}");
        assert!(
            d > 1e-3,
            "{label} decodes to the same slice as the correct pairing (max |dIC| = {d:e}); \
             this test cannot fire"
        );
    }

    // And the plane really is a plane in both cases -- so the separation above is a different
    // 2-plane rather than a degenerate basis.
    let dot = |x: &[f64; 8], y: &[f64; 8]| (0..8).map(|k| x[k] * y[k]).sum::<f64>();
    assert!(dot(&q1, &q2).abs() < 1e-15 && dot(&c1, &c2).abs() < 1e-15);
}

/// The four presets in one place: each is the reference's basis under the spec's index order, and
/// each decodes cleanly across its **whole** box.
///
/// The extent is [`Chart::default_half`] rather than a literal. It was `1.0` here and `1.0` in the
/// gallery, two copies of a number that was wrong in both: the reference UI reads
/// `Slice +/- 3.0e+0`, so every committed preset image was a 3x crop on the middle of the picture.
/// Driving both from the chart means the test cannot pass on a box the render does not use.
#[test]
fn the_four_default_presets_decode_across_their_whole_box() {
    let cases: [(&str, Chart); 4] = [
        ("shape", Chart::preset_shape()),
        ("prho", Chart::preset_prho()),
        ("plambda", Chart::preset_plambda()),
        ("shape_pl", Chart::preset_shape_pl()),
    ];
    let half = Chart::preset_shape().default_half();
    assert_eq!(half, 3.0, "the preset window is the reference UI's `Slice +/- 3.0e+0`");
    for (name, chart) in cases {
        let mut distinct = std::collections::HashSet::new();
        let mut positions = std::collections::HashSet::new();
        for iu in 0..9 {
            for iv in 0..9 {
                let (u, v) =
                    (-half + 0.25 * half * iu as f64, -half + 0.25 * half * iv as f64);
                let ic = grid::decode_state(&chart, 0, u, v);
                assert!(ic.is_finite(), "{name} at ({u},{v}) decoded non-finite");
                let s: f64 = ic.m.iter().sum();
                assert!((s - 1.0).abs() < 1e-14, "{name}: masses sum to {s}");
                // **The key must be the whole IC, not the configuration.** Positions in this
                // decode do not depend on the momentum coordinates at all, so `prho` and
                // `plambda` are constant-CONFIGURATION slices: every pixel is the same triangle
                // released with a different initial velocity. Keying on body 0's position reads
                // 1 distinct of 81 there and looks like a collapsed decode when nothing has
                // collapsed. Two different faults give the same count, and the fix is to measure
                // the quantity the chart actually moves.
                distinct.insert(format!(
                    "{:.12e},{:.12e},{:.12e},{:.12e}",
                    ic.s.r[0].x, ic.s.r[0].y, ic.s.v[0].x, ic.s.v[0].y
                ));
                positions.insert(format!("{:.12e},{:.12e}", ic.s.r[0].x, ic.s.r[0].y));
            }
        }
        // The guard against a collapsed decode: identical footprints give `ensemble_spread`
        // exactly zero, which reads as perfectly resolved and stops the descent on nothing.
        println!(
            "{name:>10}: {:>2} distinct ICs of 81, over {:>2} distinct configurations",
            distinct.len(),
            positions.len()
        );
        assert!(distinct.len() > 40, "{name} decode is collapsing: {} distinct", distinct.len());
    }
}

/// The GLSL's ten slots collapse to eight here, which is exactly when an out-of-range index gets
/// written by hand. It must panic, not alias onto `z_mu[1]`.
#[test]
#[should_panic(expected = "out of range")]
fn a_latent_index_past_seven_panics_rather_than_aliasing() {
    let _ = base_latent().get(9);
}

#[test]
#[should_panic(expected = "out of range")]
fn setting_a_latent_index_past_seven_panics_rather_than_aliasing() {
    let mut z = base_latent();
    z.set(8, 1.0);
}

/// `latent_oblique` divided by zero for a degenerate seed pair and handed back a NaN basis, which
/// decodes every pixel identically. A collapsed decode makes the criterion maximally confident, so
/// this is refused rather than propagated.
#[test]
#[should_panic(expected = "parallel to the first")]
fn latent_oblique_refuses_two_parallel_seeds() {
    let a = [0.3, -1.2, 0.5, 0.9, -0.4, 0.1, 0.7, -0.6];
    let b = a.map(|x| -2.5 * x);
    let _ = Chart::latent_oblique(base_latent(), a, b);
}

/// The two saved Config-chart slices reproduce the coordinate ranges their configs were quoted
/// with.
///
/// **The arm with teeth is the transposition.** GLSL `dimH = 0` is beta and `dimV = 1` is alpha,
/// and this port renumbers those two into the spec's order — so the horizontal basis vector is
/// spec index 1 and the vertical is spec index 0. Swapping them gives ranges that are *plausible*
/// and wrong, which is exactly how `shape_pl` and the 3x preset crop both got through. The
/// quoted ranges separate the two: `z_0` (beta) spans 0.0728 and `z_1` (alpha) spans the same
/// width about a centre two units away, so a transposition moves both by ~4.
#[test]
fn config_slices_reproduce_their_quoted_windows() {
    // (chart, cx, cy, half), the GLSL's z0 (beta) range, its z1 (alpha) range.
    let cases: [((prin_rs::grid::Chart, f64, f64, f64), [f64; 2], [f64; 2]); 2] = [
        (prin_rs::grid::Chart::config_basin(), [-1.4705, -1.3977], [2.5994, 2.6721]),
        (prin_rs::grid::Chart::config_stability(), [-0.7478, 0.5275], [-0.5355, 0.7398]),
    ];
    for (k, ((chart, cx, cy, half), want_beta, want_alpha)) in cases.iter().enumerate() {
        let prin_rs::grid::Chart::Latent { z0, q1, q2 } = chart else {
            panic!("config slice {k} is not a latent chart");
        };
        // beta is spec index 1, carried by the HORIZONTAL basis; alpha is index 0, vertical.
        let beta = [z0.get(1) + (cx - half) * q1[1], z0.get(1) + (cx + half) * q1[1]];
        let alpha = [z0.get(0) + (cy - half) * q2[0], z0.get(0) + (cy + half) * q2[0]];
        for i in 0..2 {
            assert!(
                (beta[i] - want_beta[i]).abs() < 5e-4,
                "slice {k} beta[{i}]: {} against the config's {}",
                beta[i], want_beta[i]
            );
            assert!(
                (alpha[i] - want_alpha[i]).abs() < 5e-4,
                "slice {k} alpha[{i}]: {} against the config's {}",
                alpha[i], want_alpha[i]
            );
        }
        // The negative control: with the basis transposed, both ranges are wrong by ~4 on the
        // first slice. Without this the test passes on a swapped basis whenever the two axes
        // happen to sit at similar values.
        let beta_swapped = z0.get(1) + cy * q2[0];
        assert!(
            (beta_swapped - 0.5 * (want_beta[0] + want_beta[1])).abs() > 1e-3,
            "slice {k}: the transposed basis is indistinguishable here, so this test cannot fire"
        );
    }
}
