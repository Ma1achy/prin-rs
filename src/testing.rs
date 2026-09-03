//! **Synthetic footprint fields, for tests that need a known answer.**
//!
//! Every scheduler test before this ran a real integration, which is a fine way to check that
//! the code runs and a poor way to check what it decides: on a real field nobody knows what the
//! right tree is. These fields are analytic, so the tree a policy ought to build is known in
//! advance and a test can assert it — a smooth field must be `Keep` everywhere above the
//! bootstrap, a step must be refined along the step and nowhere else, white noise must refine
//! uniformly (or stop as stationary), and a filament through a sea must separate the two.
//!
//! A field is a `Fn(&Slice, usize) -> PixelOut` and is handed to `scheduler::descend_with` in
//! place of the integrator. It fills the fields the tolerance policy and the legacy policy read:
//! the nominal `shape_vec` and `event_class`, the within-footprint spreads, and the copy count.
//! **The within spread is derived from the cell width**, because the copies span the whole cell
//! (`jitter_frac = 0.5`, offsets in `[-1, 1)²`): a smooth field's spread scales with the cell,
//! a step's spread is 1 exactly on the cells that straddle it, a sea's is 0.5 everywhere.
//!
//! Nothing here is a footprint. The states are `Bounded` at the horizon, the drift is zero, the
//! copy records are empty. It is the projection the decision reads, and no more.

use crate::ensemble::pixel::PixelOut;
use crate::grid::Slice;
use crate::outcome::State;

/// A footprint with the given nominal shape, event class and packed outcome, everything else
/// at its default. The colour-relevant projection, for tests that build reductions by hand.
pub fn fp(shape: [f64; 3], event_class: u8, outcome: u8) -> PixelOut {
    PixelOut {
        shape_vec: shape,
        event_class,
        outcome,
        state: outcome >> 2,
        detail: outcome & 0b11,
        ..Default::default()
    }
}

fn bounded(shape: [f64; 3], class: u8, spread_shape: f64, spread_event: f64, t_max: f64) -> PixelOut {
    PixelOut {
        shape_vec: shape,
        event_class: class,
        spread_shape,
        spread_event,
        ensemble_spread: spread_shape.max(spread_event),
        n_nonfinite: 0,
        state: State::Bounded as u8,
        detail: 0,
        outcome: (State::Bounded as u8) << 2,
        t_end: t_max,
        // A nominal cost, so the live descent's catch-up accounting has something to scale.
        total_substeps: 1000,
        total_force_evals: 4000,
        ..Default::default()
    }
}

/// A deterministic hash of a position, in `[0, 1)`.
fn hash01(x: f64, y: f64, seed: u64) -> f64 {
    let mut h = seed ^ 0x9e3779b97f4a7c15;
    for v in [x.to_bits(), y.to_bits()] {
        h ^= v;
        h = h.wrapping_mul(0xbf58476d1ce4e5b9);
        h ^= h >> 31;
    }
    (h >> 11) as f64 / (1u64 << 53) as f64
}

/// **Smooth**: the shape rotates slowly with `x` at rate `g` radians per unit; one class.
///
/// The copies span a cell of width `hx`, so the shape sweeps an angle `g * hx` across them and
/// the mean chord from the centroid, halved, is close to `g * hx / 4`. Halving the cell halves
/// it: `alpha = 1`, and at any tolerance above the cell's own spread the quad is resolved.
pub fn smooth(g: f64, t_max: f64) -> impl Fn(&Slice, usize) -> PixelOut + Sync {
    move |sl: &Slice, k: usize| {
        let (x, _y) = sl.decode_pos(k);
        let (hx, _hy) = sl.cell_widths();
        let a = g * x;
        bounded([a.cos(), a.sin(), 0.0], 0, 0.25 * g * hx, 0.0, t_max)
    }
}

/// **Step**: class 0 with shape `+x` for `x < x0`, class 1 with shape `-x` for `x >= x0`. A cell
/// whose copies straddle `x0` reads a spread of exactly 1 on both arms; every other cell reads
/// exactly 0. The tree a tolerance policy ought to build is the column of cells along `x0`,
/// refined to the floor, and nothing else — `O(2^L)` leaves against `4^L`.
pub fn step(x0: f64, t_max: f64) -> impl Fn(&Slice, usize) -> PixelOut + Sync {
    move |sl: &Slice, k: usize| {
        let (x, _y) = sl.decode_pos(k);
        let (hx, _hy) = sl.cell_widths();
        let straddles = (x - 0.5 * hx) < x0 && x0 <= (x + 0.5 * hx);
        let (class, shape) = if x < x0 { (0u8, [1.0, 0.0, 0.0]) } else { (1u8, [-1.0, 0.0, 0.0]) };
        let s = if straddles { 1.0 } else { 0.0 };
        bounded(shape, class, s, s, t_max)
    }
}

/// **Sea**: white noise at every scale. The shape is a hash of the position, the class one of
/// three by the same hash, and every cell's copies disagree: spread 0.5 on the shape arm and 1
/// on the event arm. No depth resolves it; a tolerance policy without a stationarity stop
/// refines it uniformly to the floor, and with one it stops.
pub fn sea(seed: u64, t_max: f64) -> impl Fn(&Slice, usize) -> PixelOut + Sync {
    move |sl: &Slice, k: usize| {
        let (x, y) = sl.decode_pos(k);
        let u = hash01(x, y, seed);
        let v = hash01(y, x, seed.wrapping_add(1));
        let th = 2.0 * std::f64::consts::PI * u;
        let z = 2.0 * v - 1.0;
        let r = (1.0 - z * z).max(0.0).sqrt();
        let class = (hash01(x + 1.0, y - 1.0, seed) * 3.0) as u8;
        bounded([r * th.cos(), r * th.sin(), z], class, 0.5, 1.0, t_max)
    }
}

/// **A filament through a sea**: [`sea`] for `x < x0`, a flat basin (one shape, one class,
/// spread 0) for `x >= x0`, and the cells straddling `x0` read spread 1. The boundary is the
/// filament; the sea is what a stationarity stop must recognise and the filament what it must
/// not.
///
/// The basin's class is a **terminal** class (collided), which the sea never takes. With the
/// basin on one of the sea's three pair classes, a quad holding the filament in its edge column
/// had a mixture that matched its parent's within noise and read as stationary -- which says
/// something true about the stop: a single hot column of eight footprints in sixty-four is
/// invisible to the coherence arm at `N = 8`, and the mixture-against-parent arm is what catches
/// an edge filament, only when the classes differ.
pub fn filament_in_sea(x0: f64, seed: u64, t_max: f64) -> impl Fn(&Slice, usize) -> PixelOut + Sync {
    let s = sea(seed, t_max);
    let basin_class = crate::ensemble::stats::TERMINAL_TAG + ((State::Collision as u8) << 2);
    move |sl: &Slice, k: usize| {
        let (x, _y) = sl.decode_pos(k);
        let (hx, _hy) = sl.cell_widths();
        let straddles = (x - 0.5 * hx) < x0 && x0 <= (x + 0.5 * hx);
        if straddles {
            bounded([0.0, 0.0, 1.0], basin_class, 1.0, 1.0, t_max)
        } else if x < x0 {
            s(sl, k)
        } else {
            bounded([0.0, 0.0, 1.0], basin_class, 0.0, 0.0, t_max)
        }
    }
}

/// **A pulse**: the [`step`] at `x0`, plus an unresolved band around it whose half-width rises
/// to `w_max` at mid-march and falls back to nothing by the horizon, `w(t) = w_max sin(pi t /
/// t_max)`. The live series carries the band at each of `n_b` boundaries; the terminal fields
/// carry the last boundary, where only the straddling column is unresolved.
///
/// The control for the live descent: a static tree at the horizon sees only the step, a static
/// tree mid-march sees the whole band, and the two are not nested — the treadmill's rise and
/// collapse, in a field where it is known. The live tree, which can only add, must contain both.
pub fn pulse(x0: f64, w_max: f64, n_b: usize, t_max: f64) -> impl Fn(&Slice, usize) -> PixelOut + Sync {
    move |sl: &Slice, k: usize| {
        let (x, _y) = sl.decode_pos(k);
        let (hx, _hy) = sl.cell_widths();
        let straddles = (x - 0.5 * hx) < x0 && x0 <= (x + 0.5 * hx);
        let (class, shape) = if x < x0 { (0u8, [1.0, 0.0, 0.0]) } else { (1u8, [-1.0, 0.0, 0.0]) };
        let mut p = bounded(shape, class, 0.0, 0.0, t_max);
        for j in 0..n_b {
            let t = t_max * (j + 1) as f64 / n_b as f64;
            let w = w_max * (std::f64::consts::PI * t / t_max).sin().max(0.0);
            let hot = straddles || (x - x0).abs() <= w + 0.5 * hx;
            let s = if hot { 1.0 } else { 0.0 };
            p.live_t.push(t);
            p.live_spread_shape.push(s);
            p.live_spread_event.push(s);
            p.live_shape.push(shape);
            p.live_class.push(class);
        }
        // Terminal = the last boundary.
        let s = *p.live_spread_shape.last().unwrap();
        p.spread_shape = s;
        p.spread_event = s;
        p.ensemble_spread = s;
        p
    }
}
