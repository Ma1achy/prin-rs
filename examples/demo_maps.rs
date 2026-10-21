//! Scratch, for the trace demo: outcome maps whose decode is the demo's own (its fromJacobi, equal masses),
//! so the background and a hovered trace are the same system. Heggie, t = 20, r_coll = 1e-3, closure escape.
//! Run: cargo run --release --example demo_maps -- sphere|momenta out.bin
use prin_rs::integrate::heggie::{integrate_hg, HgOpts};
use prin_rs::outcome::{self, EscapeRule, CLOSURE_TAU};
use prin_rs::physics::decoder;
use prin_rs::vec2::Vec2;
use rayon::prelude::*;
use std::f64::consts::PI;

const QMAX: f64 = 2.0;
fn sig(z: f64) -> f64 { 1.0 / (1.0 + (-z).exp()) }

// the demo's fromJacobi, transcribed: mass-weighted Jacobi (rt, lt), momenta (prho, plam), equal masses
fn from_jacobi(rt: [f64; 2], lt: [f64; 2], prho: [f64; 2], plam: [f64; 2]) -> ([Vec2<f64>; 3], [Vec2<f64>; 3], [f64; 3]) {
    let m: [f64; 3] = [1.0 / 3.0; 3];
    let (m0, m1, m2) = (m[0], m[1], m[2]);
    let m01 = m0 + m1;
    let (mur, mul) = (m0 * m1 / m01, m2 * m01);
    let rho = [rt[0] / mur.sqrt(), rt[1] / mur.sqrt()];
    let lam = [lt[0] / mul.sqrt(), lt[1] / mul.sqrt()];
    let r01 = [-m2 * lam[0], -m2 * lam[1]];
    let r = [
        Vec2::new(r01[0] - m1 / m01 * rho[0], r01[1] - m1 / m01 * rho[1]),
        Vec2::new(r01[0] + m0 / m01 * rho[0], r01[1] + m0 / m01 * rho[1]),
        Vec2::new(m01 * lam[0], m01 * lam[1]),
    ];
    let p = [
        Vec2::new(-prho[0] - m0 / m01 * plam[0], -prho[1] - m0 / m01 * plam[1]),
        Vec2::new(prho[0] - m1 / m01 * plam[0], prho[1] - m1 / m01 * plam[1]),
        Vec2::new(plam[0], plam[1]),
    ];
    (r, p, m)
}
fn decode_shape(th: f64, ph: f64) -> ([Vec2<f64>; 3], [Vec2<f64>; 3], [f64; 3]) {
    let n = [ph.sin() * th.cos(), ph.sin() * th.sin(), ph.cos()];
    let (a, b, pp, q) = ((1.0 + n[0]) / 2.0, (1.0 - n[0]) / 2.0, n[1] / 2.0, n[2] / 2.0);
    let psi = if pp == 0.0 && q == 0.0 { 0.0 } else { q.atan2(pp) };
    from_jacobi([a.sqrt(), 0.0], [b.sqrt() * psi.cos(), b.sqrt() * psi.sin()], [0.0; 2], [0.0; 2])
}
fn decode_shape_mom(th: f64, ph: f64, prho: [f64; 2], plam: [f64; 2]) -> ([Vec2<f64>; 3], [Vec2<f64>; 3], [f64; 3]) {
    let n = [ph.sin() * th.cos(), ph.sin() * th.sin(), ph.cos()];
    let (a, b, pp, q) = ((1.0 + n[0]) / 2.0, (1.0 - n[0]) / 2.0, n[1] / 2.0, n[2] / 2.0);
    let psi = if pp == 0.0 && q == 0.0 { 0.0 } else { q.atan2(pp) };
    from_jacobi([a.sqrt(), 0.0], [b.sqrt() * psi.cos(), b.sqrt() * psi.sin()], prho, plam)
}
// the shape sphere's folds: any theta, any phi (periods and pole crossings)
fn fold(th: f64, ph: f64) -> (f64, f64) {
    let tau = 2.0 * PI;
    let mut f = ph.rem_euclid(tau);
    let mut t = th;
    if f > PI { f = tau - f; t += PI; }
    (t.rem_euclid(tau), f)
}
// homologous expansion at rate h: v_i = h r_i, so p_i = m_i h r_i (centre of mass and Lz stay 0)
fn dilate(r: [Vec2<f64>; 3], m: [f64; 3], h: f64) -> ([Vec2<f64>; 3], [Vec2<f64>; 3], [f64; 3]) {
    (r, [r[0] * (m[0] * h), r[1] * (m[1] * h), r[2] * (m[2] * h)], m)
}
// a right triangle of the Burrau family at nu, equal masses, right angle at body 1 (0-based 0), centred, I = 1
fn right_triangle(nu: f64) -> [Vec2<f64>; 3] {
    let (a, b) = (1.0 - nu * nu, 2.0 * nu);
    let m = 1.0 / 3.0;
    let mut r = [Vec2::new(0.0, 0.0), Vec2::new(b, 0.0), Vec2::new(0.0, a)];
    let c = (r[0] + r[1] + r[2]) / 3.0;
    for k in 0..3 { r[k] = r[k] - c; }
    let i: f64 = (0..3).map(|k| m * (r[k].x * r[k].x + r[k].y * r[k].y)).sum();
    for k in 0..3 { r[k] = r[k] / i.sqrt(); }
    r
}
fn decode_momenta(x: f64, y: f64) -> ([Vec2<f64>; 3], [Vec2<f64>; 3], [f64; 3]) {
    let amin = 0.05;
    let a = amin + (PI / 2.0 - 2.0 * amin) * sig(0.0);
    let b = PI * sig(0.0);
    from_jacobi([a.cos(), 0.0], [a.sin() * b.cos(), a.sin() * b.sin()],
        [0.0, QMAX * (2.0 * sig(x) - 1.0)], [0.0, QMAX * (2.0 * sig(y) - 1.0)])
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mode = a.get(1).cloned().unwrap_or_else(|| "sphere".into());
    let out = a.get(2).cloned().unwrap_or_else(|| format!("demo_{mode}.bin"));
    let (nx, ny) = match mode.as_str() { "sphere" => (768, 384), "momenta" => (400, 400), "thetah" => (512, 256), _ => (384, 192) };
    let tilt = 30f64.to_radians();
    let kappa = 0.6 / PI;   // the target coordinate per unit of s along the tilted axis: 0.3 at s = pi
    let opts = HgOpts::<f64> { r_coll_frac: 1e-3, stop_on_event: true, stop_on_escape: true,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU), closure_k: 1, ..Default::default() };
    let t_max = 20.0;
    let px: Vec<u8> = (0..nx * ny).into_par_iter().map(|idx| {
        let (i, j) = ((idx % nx) as f64 + 0.5, (idx / nx) as f64 + 0.5);
        let (fx, fy) = (i / nx as f64, j / ny as f64);
        let (r, p, m) = match mode.as_str() {
            "sphere" => decode_shape(2.0 * PI * fx, PI * fy),          // row 0 = phi 0, the pole
            "momenta" => decode_momenta(-0.48 + 0.96 * fx, 0.48 - 0.96 * fy),   // row 0 = top = +0.48
            // the sphere tilted 30 deg toward a momentum: s in [-1.6pi, 1.6pi], phi in [pi/2 - 0.8pi, pi/2 + 0.8pi] (row 0 = top)
            "tiltz" | "tiltp" => {
                let s = -1.6 * PI + 3.2 * PI * fx;
                let ph = PI / 2.0 - 0.8 * PI + 1.6 * PI * fy;
                let (t, f) = fold(PI + s * tilt.cos(), ph);
                let d = kappa * s * tilt.sin();
                let q3 = if mode == "tiltz" { QMAX * (2.0 * sig(d) - 1.0) } else { d };
                decode_shape_mom(t, f, [0.0; 2], [0.0, q3])
            }
            // theta (phi fixed at 1.1) x dilation rate h in [+0.96 (top), -0.96]
            "thetah" => {
                let (r, _, m) = decode_shape(2.0 * PI * fx, 1.1);
                dilate(r, m, 0.96 - 1.92 * fy)
            }
            // right triangles nu in [1/32, 31/32] x dilation rate h in [+0.96 (top), -0.96]
            _ => dilate(right_triangle(1.0 / 32.0 + (30.0 / 32.0) * fx), [1.0 / 3.0; 3], 0.96 - 1.92 * fy),
        };
        let ic = decoder::to_ic(&r, &p, &m);
        let o = integrate_hg(ic.s, &ic.m, t_max, 50, 1e-2, 20_000_000, &opts);
        outcome::classify(&o.events, &o.state, &ic.m, o.finite, o.budget_exhausted).pack()
    }).collect();
    let mut hist = std::collections::BTreeMap::new();
    for p in &px { *hist.entry(*p).or_insert(0usize) += 1; }
    println!("{mode} {nx}x{ny} {hist:?}");
    std::fs::write(&out, &px).unwrap();
}
