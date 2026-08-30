//! **Is the drift field partly a picture of `t_end`?**
//!
//! `energy_drift_max` is accumulated to whenever the run *stopped*. Under `stop_on_event` a
//! trajectory that collides at `t = 5` has had a tenth as long to accumulate error as one that
//! reaches `t = 50`. So a region of common termination time is a region of common accumulated
//! drift **for a reason that has nothing to do with the integrator's local behaviour**, and the
//! bright wedges could be an artefact of the diagnostic rather than of the physics.
//!
//! That is the last untested rival, and it is the only one with any support: `collided` and
//! `terminated early` were the two candidates with lift above 1.3 in `wedge_id`, where chart
//! switching, close approach and cost all came back **depleted**.
//!
//! # Three arms, and the second is the real test
//!
//! ```text
//!   1  drift / t_end     a rate rather than a total. Cheap, but only a first-order fix:
//!                        error does not accumulate linearly in t.
//!   2  stop_on_event = false   EVERY trajectory integrates to t_max, so drift is accumulated
//!                        over the SAME interval for every pixel. This is the controlled version.
//!   3  t_end itself      rendered, so "the wedges are the termination map" can be judged by eye
//!                        as well as by a lift.
//! ```
//!
//! **The test:** does the hot set survive arm 2? `P(hot_common | hot_stopped) / P(hot_common)`.
//! A lift near 1 means the wedges are entirely a `t_end` artefact and vanish once the exposure is
//! equalised. A lift far above 1 means they are a property of the trajectories and the diagnostic
//! was not inventing them.
//!
//! # What could make arm 2 misleading, stated first
//!
//! Running past a collision is legitimate — AZ regularises binary collisions, that is the point
//! of it — but it is **not the same experiment**: the post-collision trajectory is a different
//! physical history, and its drift includes intervals the production run never integrates. So
//! arm 2 answers "is the wedge structure an artefact of unequal exposure", not "what would
//! production look like". `nonfin` and `budget` are printed per arm because a trajectory forced
//! through a close approach it would otherwise have stopped at is where an integrator fails.
//!
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};

const WINDOW: f64 = 0.4;
const DLO: f64 = 1e-12;
const DHI: f64 = 1e2;

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

fn ramp(x: f64) -> [u8; 3] {
    const S: [[f64; 3]; 5] = [
        [0.0, 0.0, 0.015], [0.34, 0.06, 0.43], [0.72, 0.21, 0.33],
        [0.98, 0.55, 0.04], [0.99, 1.0, 0.64],
    ];
    let t = x.clamp(0.0, 1.0) * 4.0;
    let i = (t.floor() as usize).min(3);
    let f = t - i as f64;
    let mut o = [0u8; 3];
    for k in 0..3 {
        o[k] = (255.0 * (S[i][k] * (1.0 - f) + S[i + 1][k] * f)).clamp(0.0, 255.0) as u8;
    }
    o
}

fn main() {
    let res: usize = arg(1, 256);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let dir = format!("{root}/step_control/drift_time");
    let _ = std::fs::create_dir_all(&dir);

    let (chart, cx, cy, half) = Chart::config_stability();
    let base = EnsembleCfg {
        refine_flagged: false,
        t_max: 50.0,
        n_sync: (50.0f64 / WINDOW).round() as usize,
        r_coll_frac: 0.005,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU),
        closure_k: 1,
        stop_on_escape: false,
        ..Default::default()
    };
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    println!("config_stability {res}^2\nconfig: {}\n", base.provenance());

    let run = |cfg: &EnsembleCfg| -> (Vec<PixelOut>, f64) {
        let t = std::time::Instant::now();
        let v: Vec<PixelOut> = (0..sl.npix())
            .into_par_iter()
            .map(|k| pixel::evaluate::<f64>(&sl, k, cfg))
            .collect();
        (v, t.elapsed().as_secs_f64())
    };

    let (a, ta) = run(&base);
    let common = EnsembleCfg { stop_on_event: false, ..base };
    let (b, tb) = run(&common);
    let n = a.len();
    println!(
        "  stopped-at-event pass {ta:.1}s;  common-time pass {tb:.1}s\n\
         {:>26} {:>10} {:>10} {:>12} {:>12}\n\
         {:>26} {:>10} {:>10} {:>12.4} {:>12.4}\n\
         {:>26} {:>10} {:>10} {:>12.4} {:>12.4}\n",
        "arm", "nonfin", "budget", "t_end p50", "drift p50",
        "stop_on_event (production)",
        a.iter().filter(|p| p.n_nonfinite > 0).count(),
        a.iter().filter(|p| p.state == 6).count(),
        q(&mut a.iter().map(|p| p.t_end).collect(), 0.5),
        q(&mut a.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect(), 0.5),
        "common time (to t_max)",
        b.iter().filter(|p| p.n_nonfinite > 0).count(),
        b.iter().filter(|p| p.state == 6).count(),
        q(&mut b.iter().map(|p| p.t_end).collect(), 0.5),
        q(&mut b.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect(), 0.5),
    );

    let hot = |v: &[PixelOut], f: &dyn Fn(&PixelOut) -> f64| -> Vec<bool> {
        let mut d: Vec<f64> = v.iter().map(|p| f(p)).filter(|x| x.is_finite()).collect();
        let cut = q(&mut d, 0.75);
        v.iter().map(|p| f(p) > cut).collect()
    };
    let hot_a = hot(&a, &|p| p.energy_drift_max);
    let hot_b = hot(&b, &|p| p.energy_drift_max);
    let hot_rate = hot(&a, &|p| p.energy_drift_max / p.t_end.max(1e-12));

    println!("== DOES THE WEDGE SET SURVIVE WHEN EXPOSURE IS EQUALISED? ==");
    println!(
        "  Base rates are 0.25 by construction (a p75 cut), so a lift of 1 is chance and 4 is\n\
         perfect agreement. **Near 1 means the wedges were a `t_end` artefact.**\n"
    );
    println!("  {:>34} {:>12} {:>8}", "set", "P(. | hot_A)", "lift");
    for (name, m) in [
        ("hot at common time (arm 2)", &hot_b),
        ("hot by RATE drift/t_end (arm 1)", &hot_rate),
    ] {
        let n_a = hot_a.iter().filter(|x| **x).count().max(1) as f64;
        let base = m.iter().filter(|x| **x).count() as f64 / n as f64;
        let p = (0..n).filter(|&i| hot_a[i] && m[i]).count() as f64 / n_a;
        println!("  {name:>34} {p:>12.4} {:>8.3}", p / base.max(f64::MIN_POSITIVE));
    }

    // And the reverse framing: how much of the drift field is explained by t_end alone?
    let mut te: Vec<f64> = a.iter().map(|p| p.t_end).collect();
    let (t25, t50, t75) = (q(&mut te.clone(), 0.25), q(&mut te.clone(), 0.5), q(&mut te, 0.75));
    println!("\n  drift p50 by `t_end` quartile (production arm):");
    for (nm, lo, hi) in [
        ("Q1 (earliest)", f64::NEG_INFINITY, t25),
        ("Q2", t25, t50),
        ("Q3", t50, t75),
        ("Q4 (latest)", t75, f64::INFINITY),
    ] {
        let mut d: Vec<f64> = a
            .iter()
            .filter(|p| p.t_end >= lo && p.t_end < hi)
            .map(|p| p.energy_drift_max)
            .filter(|x| x.is_finite())
            .collect();
        println!("  {nm:>18} n = {:>7}   drift p50 {:.3e}", d.len(), q(&mut d, 0.5));
    }

    let save = |nm: &str, buf: &[u8], note: &str| {
        let p = format!("{dir}/{nm}.png");
        let _ = prin_rs::output::adaptive::save_rect(&p, res, res, buf);
        let _ = prin_rs::output::provenance_sidecar(
            &p, &base, &format!("res={res}x{res}\ncase=config_stability\npanel={nm}\n{note}\n"),
        );
    };
    let dbuf = |v: &[PixelOut], f: &dyn Fn(&PixelOut) -> f64| -> Vec<u8> {
        v.iter()
            .flat_map(|p| {
                let x = f(p);
                if x.is_finite() && x > 0.0 {
                    ramp((x.log10() - DLO.log10()) / (DHI.log10() - DLO.log10()))
                } else if x == 0.0 {
                    ramp(0.0)
                } else {
                    [255, 0, 255]
                }
            })
            .collect()
    };
    save("drift_stopped", &dbuf(&a, &|p| p.energy_drift_max), "production, drift to t_end");
    save("drift_common_time", &dbuf(&b, &|p| p.energy_drift_max),
         "stop_on_event=false: every pixel integrated to t_max");
    save("drift_rate", &dbuf(&a, &|p| p.energy_drift_max / p.t_end.max(1e-12)),
         "drift / t_end -- a rate, first-order fix only");
    // `t_end` on its own ramp: it spans 0..t_max linearly, not decades.
    save(
        "t_end",
        &a.iter().flat_map(|p| ramp(p.t_end / base.t_max)).collect::<Vec<u8>>(),
        "t_end on a LINEAR ramp 0..t_max -- it does not span decades and a log ramp would lie",
    );
    println!("\nWrote {dir}/ -- drift_stopped, drift_common_time, drift_rate, t_end.");
}
