//! **Does `spatial::straightness` measure straightness?** The controls, before any field uses it.
//!
//! The metric exists to separate a *decision boundary* from a *fractal boundary* in the drift
//! field. Both raise the density of sharp neighbour steps; only the first is straight. So a
//! number from it is worth nothing until it is shown to score a line low and a wandering curve of
//! the **same extent and the same pixel count** high — controlling for size and count is what
//! makes it a test of shape rather than of length.
//!
//! Four properties, and the third is the one a plausible-looking wrong implementation fails:
//!
//! 1. a straight line scores ~0 and a disc scores ~1;
//! 2. a wandering curve scores far from a line **at matched extent and count**;
//! 3. it is **rotation-invariant** — a diagonal line must score like a horizontal one, or the
//!    metric is reading axis alignment, which on a chart plane whose artefacts radiate at
//!    arbitrary angles would be a confound rather than a measurement;
//! 4. fewer than three points is `NaN`, because two points are collinear by construction and a
//!    number there describes the count and not the shape.
use prin_rs::spatial::{components, straightness};

/// A line, a disc, and the two ends of the scale.
#[test]
fn a_line_scores_zero_and_a_disc_scores_one() {
    let line: Vec<(usize, usize)> = (0..100).map(|x| (x, 50)).collect();
    let s_line = straightness(&line);

    let mut disc = Vec::new();
    for y in 0..60usize {
        for x in 0..60usize {
            let (dx, dy) = (x as f64 - 29.5, y as f64 - 29.5);
            if dx * dx + dy * dy <= 29.0 * 29.0 {
                disc.push((x, y));
            }
        }
    }
    let s_disc = straightness(&disc);
    println!("line {s_line:.6}   disc ({} px) {s_disc:.6}", disc.len());

    assert!(s_line < 1e-12, "a perfect line must be 0, got {s_line}");
    assert!(s_disc > 0.9, "a disc must be near isotropic, got {s_disc}");
}

/// **The arm with teeth: matched extent, matched count, different shape.**
///
/// Both structures span `x in 0..160` and hold exactly 160 pixels. Only the perpendicular
/// wander differs. A metric that scored them alike would be measuring length.
#[test]
fn a_wandering_curve_is_separated_from_a_line_at_matched_extent_and_count() {
    let n = 160usize;
    let line: Vec<(usize, usize)> = (0..n).map(|x| (x, 80)).collect();
    // Deterministic wander, no RNG: a fixed integer hash so the fixture is reproducible.
    let wander: Vec<(usize, usize)> = (0..n)
        .map(|x| {
            let h = (x as u64).wrapping_mul(6_364_136_223_846_793_005).rotate_left(17);
            let dy = (h % 61) as i64 - 30;
            (x, (80 + dy) as usize)
        })
        .collect();

    assert_eq!(line.len(), wander.len(), "the two fixtures must have equal counts");
    let ext = |v: &[(usize, usize)]| {
        let (lo, hi) = (v.iter().map(|p| p.0).min().unwrap(), v.iter().map(|p| p.0).max().unwrap());
        hi - lo
    };
    assert_eq!(ext(&line), ext(&wander), "the two fixtures must have equal x-extent");

    let (a, b) = (straightness(&line), straightness(&wander));
    println!("matched extent {} and count {}: line {a:.6}, wander {b:.6}", ext(&line), line.len());
    assert!(a < 1e-12, "line: {a}");
    assert!(b > 0.15, "the wander is not separated from the line: {b}");
    assert!(b / a.max(1e-12) > 1e11, "line and wander score alike, so this measures length");
}

/// **Rotation invariance.** A diagonal must score like a horizontal, or the metric reads axis
/// alignment — and the artefact it is aimed at radiates at arbitrary angles.
#[test]
fn straightness_does_not_depend_on_orientation() {
    let horiz: Vec<(usize, usize)> = (0..120).map(|x| (x, 60)).collect();
    let vert: Vec<(usize, usize)> = (0..120).map(|y| (60, y)).collect();
    let diag: Vec<(usize, usize)> = (0..120).map(|i| (i, i)).collect();
    let anti: Vec<(usize, usize)> = (0..120).map(|i| (i, 119 - i)).collect();
    for (name, v) in [("horiz", &horiz), ("vert", &vert), ("diag", &diag), ("anti", &anti)] {
        let s = straightness(v);
        println!("{name:>6}: {s:.3e}");
        assert!(s < 1e-12, "{name} is not read as straight: {s}");
    }
}

/// Fewer than three points is `NaN`, never `0.0` — a number there would describe the count.
#[test]
fn two_points_are_undetermined_and_not_perfectly_straight() {
    assert!(straightness(&[]).is_nan());
    assert!(straightness(&[(1, 1)]).is_nan());
    assert!(straightness(&[(1, 1), (5, 9)]).is_nan(), "two points are collinear by construction");
    assert!(straightness(&[(0, 0), (1, 0), (2, 0)]) < 1e-12, "three collinear points are straight");
}

/// `components` must split what `layout` counts, and 4-connectivity must hold: a checkerboard is
/// scatter, not one structure.
#[test]
fn components_are_four_connected_and_ordered_largest_first() {
    let n = 8usize;
    let mut mask = vec![false; n * n];
    for x in 0..n {
        mask[3 * n + x] = true; // a full row
    }
    mask[0] = true; // an isolated pixel, diagonally adjacent to nothing on the row
    let cs = components(&mask, n);
    println!("{} components, sizes {:?}", cs.len(), cs.iter().map(|c| c.len()).collect::<Vec<_>>());
    assert_eq!(cs.len(), 2, "the isolated pixel must not join the row");
    assert_eq!(cs[0].len(), n, "largest first");
    assert_eq!(cs[1].len(), 1);
    assert!(straightness(&cs[0]) < 1e-12, "the row component is a line");

    let mut checker = vec![false; n * n];
    for y in 0..n {
        for x in 0..n {
            checker[y * n + x] = (x + y) % 2 == 0;
        }
    }
    let cc = components(&checker, n);
    assert_eq!(cc.len(), n * n / 2, "a checkerboard is {} scattered cells under 4-connectivity", n * n / 2);
}

// -------------------------------------------------------------------------------------------
// Boundary straightness: the wedge measurement, generalised off the hand-drawn mask.
// -------------------------------------------------------------------------------------------

/// **A straight-edged region reads straight; a jagged-edged region of the same area does not.**
///
/// This is the arm that matters, because a *global* fit to any closed boundary is isotropic --
/// a square outline and a blob outline both score ~1 -- so without local fitting the metric would
/// read the same for a wedge and for chaos, and would look like a null rather than a broken test.
///
/// The two fixtures are matched on **area** and differ only in edge roughness.
#[test]
fn a_straight_edge_and_a_jagged_edge_of_the_same_area_are_separated() {
    use prin_rs::spatial::boundary_straightness;
    let n = 128usize;

    // A half-plane: y < 64. Perfectly straight edge.
    let mut flat = vec![false; n * n];
    for y in 0..64 {
        for x in 0..n {
            flat[y * n + x] = true;
        }
    }
    // Same area, edge displaced by a deterministic wander of +/- up to 12 px.
    let mut jag = vec![false; n * n];
    for x in 0..n {
        let h = (x as u64).wrapping_mul(6_364_136_223_846_793_005).rotate_left(23);
        let dy = (h % 25) as i64 - 12;
        let top = (64 + dy).clamp(1, (n - 1) as i64) as usize;
        for y in 0..top {
            jag[y * n + x] = true;
        }
    }
    let (a_flat, a_jag) = (flat.iter().filter(|&&b| b).count(), jag.iter().filter(|&&b| b).count());
    println!("areas: flat {a_flat}, jagged {a_jag}");

    let (sf, sj) = (
        boundary_straightness(&flat, n, 8, 20),
        boundary_straightness(&jag, n, 8, 20),
    );
    println!("boundary straightness r=8: flat {sf:.5}, jagged {sj:.5}");
    assert!(sf < 0.05, "a half-plane edge must read straight, got {sf}");
    assert!(sj > 4.0 * sf.max(1e-6), "jagged is not separated from flat: {sj} vs {sf}");
    assert!(sj > 0.15, "the jagged edge reads too straight to discriminate: {sj}");
}

/// **The stated limit, asserted so it cannot be rediscovered as a surprise.** A large circle is
/// locally straight, and that is what "local" means -- not a defect.
#[test]
fn a_large_circle_reads_locally_straight_and_that_is_the_documented_limit() {
    use prin_rs::spatial::boundary_straightness;
    let n = 128usize;
    let mut disc = vec![false; n * n];
    for y in 0..n {
        for x in 0..n {
            let (dx, dy) = (x as f64 - 63.5, y as f64 - 63.5);
            if dx * dx + dy * dy <= 55.0 * 55.0 {
                disc[y * n + x] = true;
            }
        }
    }
    let s = boundary_straightness(&disc, n, 6, 20);
    println!("disc of radius 55, local straightness at r=6: {s:.5}");
    assert!(s < 0.15, "a curve of radius 55 must be locally straight at r=6, got {s}");
}
