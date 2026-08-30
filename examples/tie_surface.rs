//! **Do the wedge boundaries lie on the chart-switching surfaces `d_i = d_j`?**
//!
//! The clean diagnostic to run before changing any code. `choose_reference` is
//! `THIRD[argmax_k d_k]`, so it puts a Voronoi-like partition of state space into the numerical
//! method: the representation changes abruptly wherever two separations tie, and nothing in the
//! *dynamics* changes there at all. If the observed cell boundaries coincide with those loci to
//! numerical resolution, the wedges are chart-selection artefacts rather than dynamical structure.
//!
//! # Two arms, and only the second can settle it
//!
//! **`t = 0`** — the tie loci in the initial conditions, free of any integration. But the
//! reference is re-chosen at all `n_sync` boundaries, so the `t = 0` pinwheel is one surface out
//! of 125 and can only ever explain a slice of the structure. Drawn because it is nearly free and
//! because a partial match is still informative; **not** because it is the test.
//!
//! **Along the trajectory** — for each neighbouring pair whose reference itinerary differs, find
//! the FIRST boundary at which it differs and read `ref_tie` there. If the cell boundary is a tie
//! surface, that ratio sits at ~1. Compared against the ratio's distribution over *all*
//! `(pixel, boundary)` pairs, which is the null: without it, "the ratios near a switch are close
//! to 1" would be unreadable, since a chaotic slice may sit near ties everywhere.
//!
//! `ref_tie` is second-**longest** over **longest**. Not `tie_ratio`, which is the two *tightest*
//! and decides which binary is which — reading that for this question would measure the wrong
//! pair entirely.
//!
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::physics::newton;

const WINDOW: f64 = 0.4;

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

fn main() {
    let res: usize = arg(1, 256);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let dir = format!("{root}/step_control/tie_surface");
    let _ = std::fs::create_dir_all(&dir);

    let (chart, cx, cy, half) = Chart::config_stability();
    let cfg = EnsembleCfg {
        refine_flagged: false,
        t_max: 50.0,
        n_sync: (50.0f64 / WINDOW).round() as usize,
        r_coll_frac: 0.005,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU),
        closure_k: 1,
        stop_on_escape: false,
        keep_ref_path: true,
        ..Default::default()
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    println!("config_stability {res}^2\nconfig: {}\n", cfg.provenance());

    // --- arm 1: the t = 0 loci, no integration -------------------------------------------
    let t0tie: Vec<f64> = (0..sl.npix())
        .into_par_iter()
        .map(|k| {
            let (x, y) = sl.decode_pos(k);
            let st = grid::decode_state(&chart, 0, x, y);
            let mut d = newton::pair_dists(&st.s.r);
            d.sort_by(|a, b| a.partial_cmp(b).unwrap());
            d[1] / d[2].max(1e-300)
        })
        .collect();
    let mut t0s = t0tie.clone();
    println!("== ARM 1: the t = 0 tie loci (no integration) ==");
    println!(
        "  `ref_tie` at t = 0: p50 {:.4} p90 {:.4} max {:.6}\n  \
         fraction within 1% of a tie: {:.4}\n",
        q(&mut t0s.clone(), 0.5),
        q(&mut t0s.clone(), 0.9),
        q(&mut t0s, 1.0),
        t0tie.iter().filter(|&&x| x > 0.99).count() as f64 / t0tie.len() as f64
    );

    // --- arm 2: along the trajectory ------------------------------------------------------
    let t = std::time::Instant::now();
    let px: Vec<PixelOut> = (0..sl.npix())
        .into_par_iter()
        .map(|k| pixel::evaluate::<f64>(&sl, k, &cfg))
        .collect();
    println!("integrated pass {:.1}s\n", t.elapsed().as_secs_f64());

    // For each neighbour pair whose itinerary differs, the ratio at the FIRST differing boundary.
    let mut at_switch: Vec<f64> = Vec::new();
    let mut first_k: Vec<f64> = Vec::new();
    for i in 0..px.len() {
        let (x, y) = (i % res, i / res);
        for j in [
            if x + 1 < res { Some(i + 1) } else { None },
            if y + 1 < res { Some(i + res) } else { None },
        ]
        .into_iter()
        .flatten()
        {
            let (a, b) = (&px[i].ref_path, &px[j].ref_path);
            if let Some(k) = (0..a.len().min(b.len())).find(|&k| a[k] != b[k]) {
                if k < px[i].ref_tie_path.len() {
                    // Both sides of the boundary; the surface lies between them, so either
                    // pixel's ratio is an estimate of how close that boundary came to a tie.
                    at_switch.push(px[i].ref_tie_path[k].max(px[j].ref_tie_path[k]));
                    first_k.push(k as f64);
                }
            }
        }
    }
    // The null: every (pixel, boundary) ratio, which is what "near a tie" has to beat.
    let all: Vec<f64> = px.iter().flat_map(|p| p.ref_tie_path.iter().copied()).collect();

    println!("== ARM 2: at the FIRST boundary where two neighbours' itineraries differ ==");
    println!(
        "  {:>28} {:>10} {:>10} {:>10} {:>12}",
        "population", "p10", "p50", "p90", "frac > 0.99"
    );
    for (name, v) in [("at a cell boundary", &at_switch), ("all (pixel, boundary)", &all)] {
        let mut s = v.clone();
        println!(
            "  {name:>28} {:>10.4} {:>10.4} {:>10.4} {:>12.4}",
            q(&mut s.clone(), 0.1),
            q(&mut s.clone(), 0.5),
            q(&mut s, 0.9),
            v.iter().filter(|&&x| x > 0.99).count() as f64 / v.len().max(1) as f64
        );
    }
    let mut fk = first_k.clone();
    println!(
        "\n  n = {} neighbour pairs differ; first divergence at boundary p10 {:.0} p50 {:.0} \
         p90 {:.0} of {}",
        at_switch.len(),
        q(&mut fk.clone(), 0.1),
        q(&mut fk.clone(), 0.5),
        q(&mut fk, 0.9),
        cfg.n_sync
    );

    println!(
        "\nHOW TO READ THIS\n\n\
         **If `at a cell boundary` concentrates near 1 while the null does not, the wedges are\n\
         chart-selection artefacts** -- the representation changing where nothing in the dynamics\n\
         does -- and Heggie's global regularisation is the principled cure, because it removes the\n\
         argmax rather than smoothing it.\n\n\
         If the two populations agree, the itinerary differences are NOT tie crossings: they are\n\
         neighbouring trajectories that have genuinely diverged by that boundary, and the cells\n\
         are a chaos structure the reference choice merely reports. **The first-divergence\n\
         boundary index decides which reading is available**: early divergence across the frame\n\
         means ties; late means chaos had time to separate them first."
    );
}
