//! `src/chazy.rs`: one constructed state per motion class, each with the arm that says the
//! fixture actually exercises the branch it is named for.
//!
//! Masses are equal thirds and `G = M = 1`. A two-body ratio `q = K/|U|` sets the outer motion:
//! `q > 1 + eps` hyperbolic, within `eps` of 1 parabolic, below bound.

use prin_rs::chazy::{self, ChazyCfg, Gate, Latch, Motion};
use prin_rs::physics::{energy, shape, Cart};
use prin_rs::Vec2;

const M: [f64; 3] = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];

fn v(x: f64, y: f64) -> Vec2<f64> {
    Vec2::new(x, y)
}

/// Shift to the centre-of-mass frame, so `E` is the internal energy.
fn com_frame(mut s: Cart<f64>) -> Cart<f64> {
    let c = (s.r[0] * M[0] + s.r[1] * M[1] + s.r[2] * M[2]) / (M[0] + M[1] + M[2]);
    let w = (s.v[0] * M[0] + s.v[1] * M[1] + s.v[2] * M[2]) / (M[0] + M[1] + M[2]);
    for k in 0..3 {
        s.r[k] = s.r[k] - c;
        s.v[k] = s.v[k] - w;
    }
    s
}

/// Bodies 0 and 1 on a Kepler orbit of semimajor axis `a` and eccentricity `e`, placed at
/// pericentre; body 2 at distance `d` on the +x axis, receding radially with outer ratio `q_out`.
fn binary_plus_escaper(a: f64, e: f64, d: f64, q_out: f64) -> Cart<f64> {
    let mb = M[0] + M[1];
    let rp = a * (1.0 - e);
    let vp = (mb * (1.0 + e) / rp).sqrt(); // vis-viva at pericentre
    let mu_out = M[2] * mb / (M[2] + mb);
    let v_out = (2.0 * q_out * (M[2] * mb / d) / mu_out).sqrt();
    let f0 = M[1] / mb;
    let f1 = M[0] / mb;
    com_frame(Cart::new(
        [v(0.0, -f0 * rp), v(0.0, f1 * rp), v(d, 0.0)],
        [v(f0 * vp, 0.0), v(-f1 * vp, 0.0), v(v_out, 0.0)],
    ))
}

struct Probe {
    s: Cart<f64>,
    e: f64,
    i: f64,
    n: [f64; 3],
}

fn probe(s: Cart<f64>) -> Probe {
    Probe { e: energy::energy(&s.r, &s.v, &M, 0.0), i: energy::inertia(&s.r, &M), n: shape::shape_vec(&s.r, &M), s }
}

fn classify(p: &Probe, latch: &Latch, cfg: &ChazyCfg) -> Motion {
    chazy::classify(&p.s, &M, p.e, 1.0, p.i, &p.n, latch, cfg)
}

fn plain(p: &Probe) -> Motion {
    classify(p, &Latch::new(p.i), &ChazyCfg::default())
}

#[test]
fn a_bound_pair_and_a_hyperbolic_escaper_is_he() {
    let p = probe(binary_plus_escaper(1.0, 0.0, 50.0, 2.0));
    assert!(p.e < 0.0, "fixture must be a bound system overall, E = {}", p.e);
    assert_eq!(plain(&p), Motion::HE);
}

#[test]
fn the_parabolic_band_and_the_bound_side_of_it() {
    let pe = probe(binary_plus_escaper(1.0, 0.0, 50.0, 1.0));
    assert_eq!(plain(&pe), Motion::PE);
    let ex = probe(binary_plus_escaper(1.0, 0.0, 50.0, 0.5));
    assert_eq!(plain(&ex), Motion::Excursion);
    // The band is scored with its hyperbolic neighbour, and only the parabolic labels move.
    assert_eq!(Motion::PE.coarse(), Motion::HE);
    assert_eq!(Motion::HP.coarse(), Motion::H);
    assert_eq!(Motion::P.coarse(), Motion::H);
    assert_eq!(Motion::Excursion.coarse(), Motion::Excursion);
}

#[test]
fn an_eccentric_pair_at_pericentre_is_still_bound() {
    // K/U = (1+e)/2 = 0.975 at pericentre: inside the parabolic band, so a banded pair test would
    // call the pair marginal and the system a two-stage H. The pair's binding is a sign.
    let p = probe(binary_plus_escaper(1.0, 0.95, 50.0, 2.0));
    let (k, u, _, _) = chazy::pair_kup(&p.s, &M, 0, 1);
    assert!((k / u - 1.0).abs() < ChazyCfg::default().eps, "control: the pair must sit inside the band, K/U = {}", k / u);
    assert!(k < u);
    assert_eq!(plain(&p), Motion::HE);
}

#[test]
fn an_escaper_leaving_an_unbound_pair_is_a_two_stage_h() {
    // Pair (0, 1) separating with K/U = 2, escaper far and hyperbolic.
    let mb = M[0] + M[1];
    let mu = M[0] * M[1] / mb;
    let vrel = (2.0 * 2.0 * (M[0] * M[1] / 1.0) / mu).sqrt();
    let mu_out = M[2] * mb;
    let v_out = (2.0 * 2.0 * (M[2] * mb / 50.0) / mu_out).sqrt();
    let s = com_frame(Cart::new(
        [v(0.0, -0.5), v(0.0, 0.5), v(50.0, 0.0)],
        [v(0.0, -0.5 * vrel), v(0.0, 0.5 * vrel), v(v_out, 0.0)],
    ));
    assert_eq!(plain(&probe(s)), Motion::H);
}

#[test]
fn three_bodies_flying_apart_is_h() {
    let dir = [v(1.0, 0.0), v(-0.5, 0.866_025_403_784_438_6), v(-0.5, -0.866_025_403_784_438_6)];
    let s = com_frame(Cart::new(dir.map(|d| d * 50.0), dir.map(|d| d * 2.0)));
    let p = probe(s);
    assert!(p.e > 0.0);
    assert_eq!(plain(&p), Motion::H);
    assert_eq!(chazy::fewbody_resolved(&p.s, &M, p.e, 1e-5), Motion::H);
}

#[test]
fn an_interacting_positive_energy_system_is_transient_never_bounded() {
    let dir = [v(1.0, 0.0), v(-0.5, 0.866_025_403_784_438_6), v(-0.5, -0.866_025_403_784_438_6)];
    let s = com_frame(Cart::new(dir.map(|d| d * 0.577), dir.map(|d| d * 2.0)));
    let p = probe(s);
    assert!(p.e > 0.0, "control: E must be positive, got {}", p.e);
    assert_eq!(plain(&p), Motion::Transient);
    // The same configuration at rest is bound, and is NOT transient: the certificate is E's sign.
    let rest = probe(Cart::new(s.r, [v(0.0, 0.0); 3]));
    assert!(rest.e < 0.0);
    assert_eq!(plain(&rest), Motion::BSoFar);
}

#[test]
fn the_os_candidate_needs_the_latch_and_nothing_else() {
    let dir = [v(1.0, 0.0), v(-0.5, 0.866_025_403_784_438_6), v(-0.5, -0.866_025_403_784_438_6)];
    let p = probe(Cart::new(dir.map(|d| d * 0.577), [v(0.0, 0.0); 3]));
    let cfg = ChazyCfg::default();
    assert_eq!(classify(&p, &Latch::new(p.i), &cfg), Motion::BSoFar);
    let far_out_once = Latch { i_max: p.i * 2.0 * cfg.os_ratio, n_prev: None };
    assert_eq!(classify(&p, &far_out_once, &cfg), Motion::OsCandidate);
}

#[test]
fn a_mardling_aarseth_stable_hierarchy_is_b_stable_and_only_with_the_test_on() {
    // Circular inner pair a = 0.1, third body on a circular prograde outer orbit at 10.
    let mb = M[0] + M[1];
    let a_in = 0.1;
    let vin = (mb / a_in).sqrt();
    let d: f64 = 10.0;
    let vout = (1.0 / d).sqrt();
    let f0 = M[1] / mb;
    let f1 = M[0] / mb;
    let s = com_frame(Cart::new(
        [v(0.0, -f0 * a_in), v(0.0, f1 * a_in), v(d, 0.0)],
        [v(f0 * vin, 0.0), v(-f1 * vin, 0.0), v(0.0, vout)],
    ));
    let p = probe(s);
    assert_eq!(chazy::ma_stable(&p.s, &M, 0, 1, 2), Some(true));
    assert_eq!(plain(&p), Motion::BStable);
    let off = ChazyCfg { mardling_aarseth: false, ..ChazyCfg::default() };
    assert_eq!(classify(&p, &Latch::new(p.i), &off), Motion::BSoFar);
    assert_eq!(chazy::fewbody_resolved(&p.s, &M, p.e, 1e-5), Motion::BStable);
}

#[test]
fn the_reference_resolves_he_only_once_tidally_isolated() {
    // Tidal ratio 2 (m3/Mb) (a(1+e)/d)^3 = (1/d)^3: 8e-6 at d = 50, 1.25e-4 at d = 20.
    let far = probe(binary_plus_escaper(1.0, 0.0, 50.0, 2.0));
    let near = probe(binary_plus_escaper(1.0, 0.0, 20.0, 2.0));
    assert_eq!(chazy::fewbody_resolved(&far.s, &M, far.e, 1e-5), Motion::HE);
    assert_eq!(chazy::fewbody_resolved(&near.s, &M, near.e, 1e-5), Motion::Unresolved);
    // And the tidal gate in the classifier follows it.
    let tidal = ChazyCfg { gate: Gate::Tidal, ..ChazyCfg::default() };
    assert_eq!(classify(&far, &Latch::new(far.i), &tidal), Motion::HE);
    assert_ne!(classify(&near, &Latch::new(near.i), &tidal), Motion::HE);
    // The hierarchy gate commits on both: it is the faster, less conservative test.
    assert_eq!(plain(&near), Motion::HE);
}

#[test]
fn the_settling_gate_withholds_until_the_shape_has_stopped() {
    let p = probe(binary_plus_escaper(1.0, 0.0, 50.0, 2.0));
    let cfg = ChazyCfg { tau: 1e-3, ..ChazyCfg::default() };
    let moved = Latch { i_max: p.i, n_prev: Some([p.n[0] + 0.1, p.n[1], p.n[2]]) };
    let still = Latch { i_max: p.i, n_prev: Some(p.n) };
    assert_ne!(classify(&p, &moved, &cfg), Motion::HE);
    assert_eq!(classify(&p, &still, &cfg), Motion::HE);
}
