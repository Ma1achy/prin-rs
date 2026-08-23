//! The **PCG** offset path: per-pixel seeding, reproducible in isolation, and pinned.
//!
//! This is no longer the default — the spec's fixed Halton (2,3) prefix is, and it is
//! reproducible in isolation for a stronger reason: it does not depend on the pixel at all.
//! The PCG path is kept so every result measured before the switch reproduces, so it is still
//! tested, and every test here names `Scheme::Pcg` explicitly rather than relying on a default
//! that has moved.

use prin_rs::ensemble::jitter::{self, Scheme};
use prin_rs::grid;

fn pcg(s: &grid::Slice, idx: usize, seed: u64) -> Vec<prin_rs::physics::Cart<f64>> {
    jitter::copies_with::<f64>(s, idx, 7, 0.5, seed, Scheme::Pcg)
}

fn slice() -> grid::Slice {
    grid::region("near-field", 5, 5, 0.05).unwrap()
}

/// Load-bearing: it is what makes a nominal-only cross-check possible, and it is the copy
/// whose reference-body choices the shared policy propagates.
#[test]
fn copy_zero_is_always_the_un_jittered_nominal() {
    let s = slice();
    for idx in 0..s.npix() {
        let cs = pcg(&s, idx, 12345);
        let nom = s.nominal::<f64>(idx);
        for k in 0..3 {
            assert_eq!(cs[0].r[k], nom.r[k], "pixel {idx} copy 0 was jittered");
            assert_eq!(cs[0].v[k], nom.v[k]);
        }
        // and the others are not
        assert!(cs[1..].iter().all(|c| c.r[s.body] != nom.r[s.body]), "pixel {idx}: a copy was not jittered");
    }
}

/// BRIEF §7: "never from a global RNG, so any pixel is reproducible in isolation."
///
/// Both schemes satisfy this; they differ in what else they give. PCG makes each pixel's offsets
/// *independent* of every other's — which is exactly the property the fixed prefix deliberately
/// gives up, in exchange for a parent and its children sharing one perturbation pattern.
#[test]
fn pixels_are_independent_and_order_free() {
    let s = slice();
    let forward: Vec<_> = (0..s.npix()).map(|i| pcg(&s, i, 99)).collect();
    // Same pixels, evaluated in reverse order. A global stream would give different answers.
    let backward: Vec<_> = (0..s.npix()).rev().map(|i| pcg(&s, i, 99)).collect();
    for (i, f) in forward.iter().enumerate() {
        let b = &backward[s.npix() - 1 - i];
        for (cf, cb) in f.iter().zip(b.iter()) {
            assert_eq!(cf.r[s.body], cb.r[s.body], "pixel {i} depends on evaluation order");
        }
    }

    // Distinct pixels must not share a stream — under PCG.
    let a = pcg(&s, 0, 99);
    let c = pcg(&s, 1, 99);
    assert_ne!(a[1].r[s.body] - a[0].r[s.body], c[1].r[s.body] - c[0].r[s.body]);

    // And under Halton they must, which is the whole point of the change. Asserted here so the
    // contrast is a test rather than a comment.
    let ha = jitter::copies_with::<f64>(&s, 0, 7, 0.5, 99, Scheme::Halton);
    let hc = jitter::copies_with::<f64>(&s, 1, 7, 0.5, 99, Scheme::Halton);
    assert_eq!(
        ha[1].r[s.body] - ha[0].r[s.body],
        hc[1].r[s.body] - hc[0].r[s.body],
        "the fixed prefix must give distinct pixels the same offset"
    );
}

/// Jitter scales with the cell, **per axis**. The reference computes only `hx` and uses it
/// for both — latent on square grids, wrong on any other.
#[test]
fn jitter_is_bounded_by_the_per_axis_cell_width() {
    let s = grid::Slice::body_plane(5, 3, 1.0, 3.0, 0.05, 0);
    let (hx, hy) = s.cell_widths();
    assert!((hx - 0.025).abs() < 1e-15 && (hy - 0.05).abs() < 1e-15, "hx={hx} hy={hy}");
    let frac = 0.5;
    for idx in 0..s.npix() {
        // Both schemes must respect the per-axis bound, so both are checked.
        for scheme in [Scheme::Halton, Scheme::Pcg] {
        let cs = jitter::copies_with::<f64>(&s, idx, 7, frac, 7, scheme);
        let nom = s.nominal::<f64>(idx);
        for c in cs.iter().skip(1) {
            let d = c.r[s.body] - nom.r[s.body];
            assert!(d.x.abs() <= frac * hx, "x jitter {} exceeded {}", d.x, frac * hx);
            assert!(d.y.abs() <= frac * hy, "y jitter {} exceeded {}", d.y, frac * hy);
        }
        }
    }
}

/// Pinned, so a refactor or a dependency change cannot move the initial conditions
/// underneath a measurement without the test noticing.
#[test]
fn golden_jitter_values() {
    let s = slice();
    let cs = pcg(&s, 12, 0);
    let nom = s.nominal::<f64>(12);
    let got: Vec<(f64, f64)> = cs
        .iter()
        .skip(1)
        .take(3)
        .map(|c| {
            let d = c.r[s.body] - nom.r[s.body];
            (d.x, d.y)
        })
        .collect();
    for (i, (x, y)) in got.iter().enumerate() {
        println!("copy {}: dx = {x:.17e}, dy = {y:.17e}", i + 1);
    }
    // Measured, then pinned. Bounded by 0.5 * hx = 0.0125 for this slice.
    let want = [
        (-6.04696524908354682e-3, 9.09101887490315619e-3),
        (6.33724372630850574e-3, 1.17315358118474933e-2),
        (3.73278530955589716e-4, 2.70285819262516824e-3),
    ];
    for (i, ((gx, gy), (wx, wy))) in got.iter().zip(want.iter()).enumerate() {
        assert!((gx - wx).abs() < 1e-18 && (gy - wy).abs() < 1e-18,
                "copy {} moved: got ({gx:.17e}, {gy:.17e})", i + 1);
    }
}
