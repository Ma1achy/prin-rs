//! Is the concentric banding a SYNC-CADENCE artefact, and does it reach the criterion?
//!
//! # The suspected artefact
//!
//! Committed charts (`preset_plambda_uniform.png`, `preset_prho_uniform.png`,
//! `preset_shape_pl_h1_uniform.png`) carry fine dotted arcs across their smooth regions, with
//! spacing that widens outward -- the signature of contours of a smooth function at uniform
//! value intervals. If `t_end` is quantised, everything derived from it renders those steps.
//!
//! # The mechanism, read out of the driver rather than guessed
//!
//! **Collision is sampled inside the RK4 loop** (`driver.rs`, `tc = t + s.t`) and carries
//! step resolution. **Escape is sampled only at sync boundaries**, where the state is already
//! Cartesian and every trajectory shares a playhead -- the reference's cadence, transcribed. So
//! `t_end` is quantised to `n_sync` values **exactly where escape is the terminating event**,
//! and is continuous where collision is.
//!
//! That makes a prediction sharp enough to be wrong: Burrau's `near-field` at `t = 13` has a
//! silent escape arm and terminates by collision, so its `t_end` should be continuous; the
//! latent charts run `escape_fraction` 0.9894-1.0000, so theirs should take about `n_sync`
//! values. **The banding should appear on the second set and not the first.**
//!
//! # The test that was proposed cannot fire, and this one can
//!
//! *"Recount `frac_hot_between`'s distinct values: 45 -> thousands means quantisation"* -- but
//! `frac_hot_between` is `frac_above_tau_between`, a fraction of the quad's `N^2` footprints. At
//! `N = 8` it can take **at most 65 distinct values by construction**, whatever the cadence
//! does. Its 45 is 45 of an arithmetic ceiling of 65, and the corpus's own `31 / 65 / 64`
//! is that ceiling showing. A test whose outcome is capped below its own decision threshold
//! reports the same answer under both hypotheses.
//!
//! So the count that decides it is **`t_end` itself**: unbounded, and the quantity actually
//! being quantised. It is reported here alongside how many `t_end` values land *exactly* on a
//! sync boundary, which is the mechanism rather than a proxy for it. `frac_hot_between` is
//! recounted too, with its ceiling printed beside it so the number is read correctly.
//!
//! # Writes
//!
//! stdout only.

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::metric::{self, Key};
use prin_rs::quad::{Agg, Criterion};

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

const SHIPPING: metric::Colouring =
    metric::Colouring::Bivariate(prin_rs::output::colour::Scalar::ShapeSpread);

/// Distinct values and the modal share, over the bit patterns so `NaN` counts as itself.
fn distinct(v: &[f64]) -> (usize, f64) {
    let mut b: Vec<u64> = v.iter().map(|x| x.to_bits()).collect();
    b.sort_unstable();
    let (mut d, mut best, mut i) = (0usize, 0usize, 0usize);
    while i < b.len() {
        let mut j = i;
        while j < b.len() && b[j] == b[i] {
            j += 1;
        }
        d += 1;
        best = best.max(j - i);
        i = j;
    }
    (d, best as f64 / b.len().max(1) as f64 * 100.0)
}

struct Target {
    name: &'static str,
    chart: Chart,
    cx: f64,
    cy: f64,
    half: f64,
    body: usize,
}

fn main() {
    let levels: u32 = arg(1, 4);
    let n: usize = arg(2, 8);
    let t_max: f64 = arg(3, 13.0);
    // Stride of the in-loop escape test. 1 is every RK4 step.
    let fine: usize = arg(4, 1);
    // The persistence guard. `1` requires an in-loop escape to still hold at the next boundary
    // before it is accepted; `0` accepts on first sight. Both arms are runnable because the
    // difference between them IS the measurement -- unguarded, `deep interior`'s escape fraction
    // reads 0.5494 against a guarded 0.1564, and 895 of 895 of the difference re-bind.
    let confirm: bool = arg::<usize>(5, 1) != 0;
    let res = (1usize << levels) * n;

    let base = EnsembleCfg::default();
    let n_sync = ((base.n_sync as f64) * t_max / base.t_max).round().max(2.0) as usize;
    let dt_sync = t_max / n_sync as f64;

    let mut targets: Vec<Target> = Vec::new();
    for &(region, cx, cy, body) in grid::REGIONS.iter() {
        if ["near-field", "deep interior"].contains(&region) {
            targets.push(Target {
                name: region,
                chart: Chart::BodyPlane,
                cx,
                cy,
                half: 0.05,
                body,
            });
        }
    }
    // The two the banding was read off. `preset_shape_pl_h1` also carries the crisp polygonal
    // edges, which is a separate question and gets a separate render.
    for (nm, ch, cx, cy, half) in grid::gallery_cases() {
        if ["preset_plambda", "preset_shape_pl_h1"].contains(&nm) {
            targets.push(Target { name: nm, chart: ch, cx, cy, half, body: 0 });
        }
    }

    println!(
        "level {levels}, N={n}, res {res}^2, t={t_max}, n_sync={n_sync} \
         (boundaries {dt_sync:.6} apart)\n\
         escape_every: 0 = the reference's boundary-only cadence; {fine} = every {fine} RK4 step(s)\n\
         escape_confirm: {confirm} -- an in-loop escape must still hold at the next boundary\n\
         \n\
         COLLISION is sampled inside the RK4 loop and is already continuous. ESCAPE is sampled\n\
         at sync boundaries. So t_end is quantised exactly where ESCAPE terminates, and the\n\
         prediction is that Burrau (collision-terminated) is clean and the latent charts are not.\n"
    );

    for tg in &targets {
        println!("=== {} ===", tg.name);
        let mut prev: Option<(usize, usize)> = None;
        for ev in [0usize, 32, 4, fine] {
            let ens = EnsembleCfg {
                refine_flagged: false,
                t_max,
                n_sync,
                keep_boundary_shapes: true,
                escape_every: ev,
                escape_confirm: confirm,
                ..Default::default()
            };
            let t0 = std::time::Instant::now();
            let (mut caches, px_of) = metric::build_multi_with_footprints(
                tg.name, tg.cx, tg.cy, tg.half, tg.body, tg.chart, levels, n, res, 1e-4, &ens,
                &[SHIPPING],
            );
            let cache = caches.pop().unwrap();
            let secs = t0.elapsed().as_secs_f64();

            let keys: Vec<Key> = {
                let mut v: Vec<Key> = cache.quads.keys().copied().collect();
                v.sort_unstable();
                v
            };
            let all: Vec<&PixelOut> = keys.iter().flat_map(|k| px_of[k].iter()).collect();
            let nf = all.len() as f64;

            let t_end: Vec<f64> = all.iter().map(|p| p.t_end).collect();
            let (dt, mt) = distinct(&t_end);

            // **The mechanism, not a proxy for it.** How many `t_end` land exactly on one of the
            // `n_sync` boundary times. A quantised `t_end` is not merely coarse -- its values
            // ARE the boundaries, and that is checkable.
            let on_boundary = t_end
                .iter()
                .filter(|x| x.is_finite())
                .filter(|x| {
                    let k = (*x / dt_sync).round();
                    (*x - k * dt_sync).abs() <= 1e-9 * t_max
                })
                .count();

            // Only escape-terminated footprints CAN be quantised, so the fraction bounds the
            // effect before any count is read.
            let esc = all.iter().filter(|p| p.state == 0).count();
            let coll = all.iter().filter(|p| p.state == 2).count();
            let run = all.iter().filter(|p| p.state == 3).count();

            let quad = |c: Criterion| -> (usize, f64) {
                let v: Vec<f64> =
                    keys.iter().map(|k| cache.get(*k).red.signal(c, Agg::Median)).collect();
                distinct(&v)
            };
            let (dfhb, mfhb) = quad(Criterion::FracHotBetween);
            let (dfhw, mfhw) = quad(Criterion::FracHotWithin);
            let (dtg, mtg) = quad(Criterion::TerminationGradient);
            let (dw, mw) = quad(Criterion::Within);
            let sev: Vec<f64> = all.iter().map(|p| p.spread_event).collect();
            let ssh: Vec<f64> = all.iter().map(|p| p.spread_shape).collect();
            let (dsev, msev) = distinct(&sev);
            let (dssh, mssh) = distinct(&ssh);

            println!(
                "  escape_every={ev:<3} built {secs:>7.1}s   footprints {}   \
                 escape {:.4}  collision {:.4}  running {:.4}",
                all.len(),
                esc as f64 / nf,
                coll as f64 / nf,
                run as f64 / nf
            );
            println!(
                "      t_end          distinct {dt:>7}  modal {mt:>6.2}%   \
                 ON A SYNC BOUNDARY: {on_boundary} of {} ({:.2}%)",
                all.len(),
                on_boundary as f64 / nf * 100.0
            );
            println!(
                "      spread_event   distinct {dsev:>7}  modal {msev:>6.2}%      \
                 spread_shape distinct {dssh:>7}  modal {mssh:>6.2}%"
            );
            println!(
                "      per-quad: frac_hot_between {dfhb:>4} (ceiling {}) modal {mfhb:>5.1}%   \
                 frac_hot_within {dfhw:>4} modal {mfhw:>5.1}%",
                n * n + 1
            );
            println!(
                "                term_grad {dtg:>8} modal {mtg:>5.1}%   \
                 within {dw:>8} modal {mw:>5.1}%"
            );
            // **Is a finer escape test a FIX or a BUG?** `escape_candidate` asks whether a body
            // is unbound and receding. The reference tests it only at sync boundaries, where the
            // state is Cartesian and every trajectory shares a playhead. Tested inside the RK4
            // loop it is also being asked DURING a close encounter, where a pair's instantaneous
            // two-body energy can transiently read positive -- so a finer cadence can fire
            // spuriously rather than merely fire earlier.
            //
            // The discriminator is `d_min`: a genuine escape does not coincide with an ultra-close
            // approach. If the newly-escaping footprints carry small `d_min_true`, they fired
            // mid-encounter and the finer test is wrong, not finer.
            {
                let qs = |v: &mut Vec<f64>| -> (f64, f64, f64) {
                    if v.is_empty() {
                        return (f64::NAN, f64::NAN, f64::NAN);
                    }
                    (
                        prin_rs::stats::quantile(v, 0.1),
                        prin_rs::stats::quantile(v, 0.5),
                        prin_rs::stats::quantile(v, 0.9),
                    )
                };
                let mut de: Vec<f64> = all
                    .iter()
                    .filter(|p| p.state == 0)
                    .map(|p| p.d_min_true)
                    .filter(|x| x.is_finite())
                    .collect();
                let mut dc: Vec<f64> = all
                    .iter()
                    .filter(|p| p.state == 2)
                    .map(|p| p.d_min_true)
                    .filter(|x| x.is_finite())
                    .collect();
                let (e1, e5, e9) = qs(&mut de);
                let (c1, c5, c9) = qs(&mut dc);
                println!(
                    "                d_min_true | escaped  p10 {e1:>9.3e} p50 {e5:>9.3e} \
                     p90 {e9:>9.3e}  (n {})",
                    de.len()
                );
                println!(
                    "                           | collided p10 {c1:>9.3e} p50 {c5:>9.3e} \
                     p90 {c9:>9.3e}  (n {})",
                    dc.len()
                );
            }
            if let Some((p_dt, p_ob)) = prev {
                println!(
                    "      DELTA vs the reference cadence: t_end distinct {p_dt} -> {dt} \
                     ({:.1}x), on-boundary {p_ob} -> {on_boundary}",
                    dt as f64 / p_dt.max(1) as f64
                );
            }
            prev = Some((dt, on_boundary));
        }
        println!();
    }

    println!(
        "HOW TO READ THIS\n\n\
         `t_end distinct` is the decisive number, and `ON A SYNC BOUNDARY` is the mechanism.\n\
         A quantised t_end does not merely take few values -- its values ARE the boundary\n\
         times, so a high on-boundary fraction at escape_every=0 that collapses at the fine\n\
         cadence is the artefact, demonstrated rather than inferred.\n\n\
         `frac_hot_between` is printed WITH ITS CEILING because it is a fraction over N^2\n\
         footprints and cannot exceed N^2+1 distinct values however fine the cadence gets.\n\
         Reading a saturation there as evidence about the physics is reading the resolution\n\
         of a 64-sample fraction.\n\n\
         `escape` vs `collision` bounds the whole effect before any count is read: only an\n\
         escape-terminated footprint can carry a boundary time, so a region terminating by\n\
         collision cannot be contaminated by this at all -- whatever its images look like.\n\n\
         `d_min_true` by terminal state decides whether the finer cadence is a FIX or a BUG.\n\
         `escape_candidate` asks whether a body is unbound and receding, and the reference asks\n\
         it only at boundaries, where the state is Cartesian and shared. Asked every RK4 step it\n\
         is also asked DURING a close encounter, where a pair's instantaneous two-body energy can\n\
         read positive transiently. If the escaped population's d_min collapses toward the\n\
         collided population's as the stride tightens, the test is firing mid-encounter and the\n\
         finer cadence is wrong rather than finer."
    );
}
