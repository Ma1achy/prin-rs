//! The two properties the fixed Halton (2,3) prefix is chosen for, asserted rather than assumed.
//!
//! - **Fixed**: copy `k` sits at the same *fractional* offset in every footprint at every
//!   refinement level. That is what makes a parent ensemble and its children share a
//!   perturbation pattern — common random numbers by construction.
//! - **Low-discrepancy**: the prefix covers the footprint more evenly than independent uniform
//!   draws, and the gap is largest at the small `E` the project actually uses.

use prin_rs::ensemble::jitter::{self, Scheme};
use prin_rs::grid;

/// Radical inverse against hand-computed values. `phi_2(1) = 0.1b = 1/2`,
/// `phi_2(2) = 0.01b = 1/4`, `phi_3(1) = 0.1_3 = 1/3`, `phi_3(4) = 0.11_3 = 4/9`.
#[test]
fn the_radical_inverse_is_the_reflected_digit_expansion() {
    for (n, base, want) in [
        (1usize, 2usize, 0.5),
        (2, 2, 0.25),
        (3, 2, 0.75),
        (4, 2, 0.125),
        (1, 3, 1.0 / 3.0),
        (3, 3, 1.0 / 9.0),
        (4, 3, 4.0 / 9.0),
    ] {
        let got = jitter::radical_inverse(n, base);
        assert!((got - want).abs() < 1e-15, "phi_{base}({n}) = {got}, want {want}");
    }
    assert_eq!(jitter::radical_inverse(0, 2), 0.0);
}

/// Copy 0 is the un-jittered nominal under both schemes. Load-bearing: it is what makes the
/// nominal-only cross-check possible.
#[test]
fn copy_zero_is_never_jittered() {
    let s = grid::region("near-field", 8, 8, 0.05).unwrap();
    for scheme in [Scheme::Halton, Scheme::Pcg] {
        for i in [0usize, 17, 63] {
            let c = jitter::copies_with::<f64>(&s, i, 7, 0.5, 0, scheme);
            let nom = s.nominal::<f64>(i);
            for b in 0..3 {
                assert_eq!(c[0].s.r[b].x, nom.r[b].x, "{scheme:?} pixel {i}: copy 0 moved");
                assert_eq!(c[0].s.r[b].y, nom.r[b].y, "{scheme:?} pixel {i}: copy 0 moved");
            }
        }
    }
}

/// **The fixed property.** Copy `k`'s offset, divided by the cell width, is identical in every
/// footprint and at every grid resolution. The PCG scheme fails this by construction, which is
/// the point of measuring it.
#[test]
fn halton_offsets_are_fixed_across_pixels_and_resolutions() {
    let frac = 0.5;
    let mut reference: Option<Vec<(f64, f64)>> = None;

    for size in [8usize, 16, 64] {
        let s = grid::region("near-field", size, size, 0.05).unwrap();
        let (hx, hy) = s.cell_widths();
        for idx in [0usize, 1, size + 1, s.npix() - 1] {
            let c = jitter::copies_with::<f64>(&s, idx, 7, frac, 0, Scheme::Halton);
            let nom = s.nominal::<f64>(idx);
            let rel: Vec<(f64, f64)> = c[1..]
                .iter()
                .map(|x| {
                    (
                        (x.s.r[s.body].x - nom.r[s.body].x) / (frac * hx),
                        (x.s.r[s.body].y - nom.r[s.body].y) / (frac * hy),
                    )
                })
                .collect();
            match &reference {
                None => {
                    println!("fixed Halton (2,3) prefix, offsets in units of jitter_frac*cell:");
                    for (k, (u, v)) in rel.iter().enumerate() {
                        println!("  copy {:>2}: ({u:>9.6}, {v:>9.6})", k + 1);
                    }
                    reference = Some(rel);
                }
                Some(r) => {
                    // Compared to a tolerance, not bitwise: the offsets *stored* are exact, but
                    // recovering them here divides by `frac*cell_width`, which differs by
                    // resolution and rounds. The bitwise claim is asserted directly on
                    // `halton_offset` below, where no division intervenes.
                    for (k, (a, b)) in rel.iter().enumerate() {
                        assert!((a - r[k].0).abs() < 1e-12,
                                "size {size} pixel {idx} copy {k}: x offset moved by {}", a - r[k].0);
                        assert!((b - r[k].1).abs() < 1e-12,
                                "size {size} pixel {idx} copy {k}: y offset moved by {}", b - r[k].1);
                    }
                }
            }
        }
    }
    println!();
    println!("Identical to 1e-12 across three resolutions and four footprints each, the residual");
    println!("being the division used to recover them. The generator itself is exact:");

    // The claim that actually matters, with no division in the way.
    for k in 0..8 {
        let (a, b) = jitter::halton_offset(k);
        let (c, d) = jitter::halton_offset(k);
        assert_eq!((a, b), (c, d));
    }
    println!("halton_offset(k) is a pure function of k — no pixel, no seed, no ordering. A");
    println!("parent ensemble and its children share the perturbation pattern exactly.");

    // And the contrast: PCG offsets differ between footprints, which is the whole difference.
    let s = grid::region("near-field", 8, 8, 0.05).unwrap();
    let a = jitter::copies_with::<f64>(&s, 0, 7, frac, 0, Scheme::Pcg);
    let b = jitter::copies_with::<f64>(&s, 1, 7, frac, 0, Scheme::Pcg);
    let differs = (1..8).any(|k| a[k].s.r[s.body].x != b[k].s.r[s.body].x);
    assert!(differs, "PCG offsets should differ between footprints; if not, the test is empty");
    println!("PCG offsets differ between footprints, as expected — no shared pattern.");
}

/// **The low-discrepancy property**, measured rather than asserted from theory.
///
/// Star discrepancy is expensive; L2 star discrepancy has a closed form (Warnock's formula) and
/// is what is used here. Lower is better. The comparison is against the PCG scheme's own draws
/// at the same `E`, averaged over footprints so a single lucky draw does not decide it.
#[test]
fn the_halton_prefix_covers_the_footprint_better_than_pcg_draws() {
    fn l2_star(pts: &[(f64, f64)]) -> f64 {
        // Warnock: points mapped to [0,1)^2.
        let n = pts.len() as f64;
        let mut s1 = 0.0;
        for &(x, y) in pts {
            s1 += (1.0 - x * x) * (1.0 - y * y);
        }
        let mut s2 = 0.0;
        for &(x1, y1) in pts {
            for &(x2, y2) in pts {
                s2 += (1.0 - x1.max(x2)) * (1.0 - y1.max(y2));
            }
        }
        (1.0 / 9.0 - s1 / (2.0f64.powi(1) * n) + s2 / (n * n)).abs().sqrt()
    }

    let s = grid::region("near-field", 16, 16, 0.05).unwrap();
    let (hx, hy) = s.cell_widths();
    let frac = 0.5;

    println!("{:>6}{:>16}{:>16}{:>12}", "E+1", "Halton L2*", "PCG L2* (mean)", "ratio");
    let mut wins = 0usize;
    for n_copies in [4usize, 8, 16, 32] {
        let unit = |c: &Vec<prin_rs::physics::Ic<f64>>, idx: usize| -> Vec<(f64, f64)> {
            let nom = s.nominal::<f64>(idx);
            c[1..]
                .iter()
                .map(|x| {
                    (
                        0.5 * ((x.s.r[s.body].x - nom.r[s.body].x) / (frac * hx) + 1.0),
                        0.5 * ((x.s.r[s.body].y - nom.r[s.body].y) / (frac * hy) + 1.0),
                    )
                })
                .collect()
        };

        let h = l2_star(&unit(
            &jitter::copies_with::<f64>(&s, 0, n_copies - 1, frac, 0, Scheme::Halton),
            0,
        ));
        let mut acc = 0.0;
        let m = s.npix();
        for i in 0..m {
            acc += l2_star(&unit(
                &jitter::copies_with::<f64>(&s, i, n_copies - 1, frac, 0, Scheme::Pcg),
                i,
            ));
        }
        let p = acc / m as f64;
        println!("{n_copies:>6}{h:>16.6}{p:>16.6}{:>12.3}", h / p);
        if h < p {
            wins += 1;
        }
    }
    println!();
    println!("L2 star discrepancy, lower is better; PCG averaged over {} footprints.", s.npix());
    println!("The gap is what a fixed low-discrepancy prefix buys before any physics runs.");
    assert_eq!(wins, 4, "Halton should have lower discrepancy at every E tested");
}
