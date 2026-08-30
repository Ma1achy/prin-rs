//! **Where does this integrator advance with a step it knows it cannot resolve, and does that
//! coincide with the artefact?**
//!
//! The hypothesis under test: the pale, cut-out regions on `config_stability` are a *saturation*
//! boundary — inside it the substep control has hit a ceiling and the march advances anyway with
//! an under-resolved step, outside it the physics is fine. Cut out, resumes. Chaos does not cut
//! and resume; a threshold does.
//!
//! # The three stop reasons, and only one of them is terminal
//!
//! ```text
//!   budget_exhausted   max_steps per interval    TERMINAL -- the march breaks, the run ends
//!   ab_floored         A or B clamped to TINY    ADVANCE ANYWAY, fabricated denominator
//!   n_cap_hits         PerStepInterval's cap     ADVANCE ANYWAY, refused the step it wanted
//! ```
//!
//! `ab_floored` and `ab_min` were written by `AzOut` on every march since they were added and
//! read by **nothing** — they stopped one layer below `PixelOut`, so no render, dump, criterion
//! or test could see the floor fire. A sticky bit nothing reads is indistinguishable from one
//! that never fires. They are plumbed now, along with `dt_max` (the largest *physical* step any
//! copy took — `ab_min` says the denominator was small, `dt_max` says how far the step went
//! because of it) and `n_cap_hits`.
//!
//! # The controls, because a mask that covers everything and one that covers nothing look
//! equally convincing
//!
//! - **Per-box fractions** over the sixteen marked regions, with the six the report calls BROKEN
//!   and the ten it calls SOUND labelled in the table. **`B4` is the negative control**: its pale
//!   wedge *survives refinement* and is real structure. A mask firing there too is reading the
//!   chart, not the fault.
//! - **The lift**, which is the number the image cannot give: `P(mask | error_ratio > 10)`
//!   against `P(mask | error_ratio <= 10)`. `error_ratio` is the project's own objective flag for
//!   *this pixel is not data*, so it stands in for "the artefact" without my hand-drawing it.
//!   A lift near 1 means the mask and the fault are unrelated however alike the pictures look.
//! - **The frame base rate is printed first.** A mask at 90% has a lift of ~1 by arithmetic.
//!
//! # Settings
//!
//! `closure_render`'s exactly — `t_max = 50`, `r_coll = 0.005`, `n_sync = round(t_max/0.4)`,
//! `EscapeRule::Closure`, `closure_k = 1`, `stop_on_escape = false` — and **`refine_flagged`
//! OFF**, deliberately: that is the state the committed panel is in, and repairing it first would
//! switch off the diagnostic being measured.
//!
//! # Writes
//!
//! `<root>/saturation/`, both arguments, defaulting to `results/saturation`. Never into a
//! committed render directory.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour};

const WINDOW: f64 = 0.4;

/// `(name, u, v, half)` as fractions of the panel, origin top-left, digitised from the mark-up,
/// and `verdict` from `BUG_REPORT.md` §4. The verdict is carried so the table cannot be read
/// without it.
const BOXES: [(&str, f64, f64, f64, &str); 16] = [
    ("B1", 0.890, 0.245, 0.0517, "BROKEN"), ("B2", 0.862, 0.332, 0.0571, "BROKEN"),
    ("B3", 0.379, 0.446, 0.0326, "sound"),  ("B4", 0.476, 0.497, 0.0294, "SOUND-CTRL"),
    ("B5", 0.590, 0.437, 0.0294, "marginal"), ("B6", 0.608, 0.528, 0.0337, "sound"),
    ("B7", 0.419, 0.664, 0.0381, "sound"),  ("B8", 0.383, 0.748, 0.0403, "BROKEN"),
    ("B9", 0.807, 0.838, 0.0566, "BROKEN"), ("B10", 0.942, 0.789, 0.0522, "BROKEN"),
    ("P1", 0.539, 0.199, 0.0354, "sound"),  ("P2", 0.447, 0.426, 0.0354, "sound"),
    ("P3", 0.533, 0.457, 0.0305, "sound"),  ("P4", 0.335, 0.682, 0.0408, "sound"),
    ("P5", 0.428, 0.718, 0.0381, "BROKEN"), ("P6", 0.510, 0.742, 0.0397, "marginal"),
];

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
    let res: usize = arg(1, 512);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let dir = format!("{root}/saturation");
    let _ = std::fs::create_dir_all(&dir);

    let (chart, cx, cy, half) = Chart::config_stability();
    let (t_max, r_coll) = (50.0, 0.005);
    let n_sync = (t_max / WINDOW).round().max(4.0) as usize;

    // `refine_flagged: false` -- the committed panel's state, deliberately. See the module docs.
    let ens = EnsembleCfg {
        refine_flagged: false,
        t_max,
        n_sync,
        r_coll_frac: r_coll,
        escape_rule: EscapeRule::Closure(CLOSURE_TAU),
        closure_k: 1,
        stop_on_escape: false,
        ..Default::default()
    };

    println!(
        "config_stability, {res}^2, closure_render's settings, refine_flagged=OFF.\n\
         t_max={t_max} n_sync={n_sync} r_coll={r_coll} eta={} max_steps={} E+1={}\n\
         config: {}\n",
        ens.eta, ens.max_steps, ens.n_extra + 1, ens.provenance()
    );

    let t0 = std::time::Instant::now();
    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let px: Vec<PixelOut> = (0..sl.npix())
        .into_par_iter()
        .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl, k, &ens))
        .collect();
    println!("{:.1}s\n", t0.elapsed().as_secs_f64());

    // --- the masks -----------------------------------------------------------------------
    let floored: Vec<bool> = px.iter().map(|p| p.ab_floored).collect();
    let capped: Vec<bool> = px.iter().map(|p| p.n_cap_hits > 0).collect();
    let budget: Vec<bool> = px.iter().map(|p| p.budget_exhausted).collect();
    let union: Vec<bool> = (0..px.len()).map(|i| floored[i] || capped[i] || budget[i]).collect();
    // The stand-in for "the artefact", objective and already the project's own flag.
    let bad: Vec<bool> = px.iter().map(|p| p.error_ratio > ens.refine_threshold).collect();

    let frac = |m: &[bool]| m.iter().filter(|x| **x).count() as f64 / m.len() as f64;
    println!("== FRAME BASE RATES ==");
    println!("  a mask at 90% has a lift of ~1 by arithmetic. Read these before the lifts.\n");
    println!("  {:>16} {:>10}", "mask", "fraction");
    for (n, m) in [
        ("ab_floored", &floored), ("n_cap_hits>0", &capped),
        ("budget_exhausted", &budget), ("union", &union),
        ("error_ratio>10", &bad),
    ] {
        println!("  {n:>16} {:>10.6}", frac(m));
    }

    println!("\n== THE LIFT -- does the mask pick out the fault? ==");
    println!(
        "  P(mask | err>10) against P(mask | err<=10). **A lift near 1 means the mask and the\n\
         fault are unrelated however alike the pictures look.** `n` is the conditioning count.\n"
    );
    let n_bad = bad.iter().filter(|x| **x).count();
    let n_ok = bad.len() - n_bad;
    println!("  {:>16} {:>12} {:>12} {:>10}", "mask", "P(m|bad)", "P(m|ok)", "lift");
    for (nm, m) in [
        ("ab_floored", &floored), ("n_cap_hits>0", &capped),
        ("budget_exhausted", &budget), ("union", &union),
    ] {
        let pb = (0..px.len()).filter(|&i| bad[i] && m[i]).count() as f64 / n_bad.max(1) as f64;
        let po = (0..px.len()).filter(|&i| !bad[i] && m[i]).count() as f64 / n_ok.max(1) as f64;
        println!("  {nm:>16} {pb:>12.6} {po:>12.6} {:>10.3}", pb / po.max(f64::MIN_POSITIVE));
    }
    println!("  (n_bad = {n_bad}, n_ok = {n_ok})");

    // --- dt_max and ab_min, the continuous reads ------------------------------------------
    println!("\n== THE CONTINUOUS READS ==");
    println!(
        "  `dt_max` is the largest physical step any copy took; nominal is eta*dt_sync = {:.3e}.\n\
         A step-control cliff shows here and in no other recorded quantity.\n",
        ens.eta * t_max / n_sync as f64
    );
    for (nm, sel) in [("err>10", true), ("err<=10", false)] {
        let mut dt: Vec<f64> = (0..px.len()).filter(|&i| bad[i] == sel)
            .map(|i| px[i].dt_max).filter(|x| x.is_finite()).collect();
        let mut ab: Vec<f64> = (0..px.len()).filter(|&i| bad[i] == sel)
            .map(|i| px[i].ab_min).filter(|x| x.is_finite()).collect();
        let mut ch: Vec<f64> = (0..px.len()).filter(|&i| bad[i] == sel)
            .map(|i| px[i].n_cap_hits as f64).collect();
        println!(
            "  {nm:>8}  dt_max p50 {:.3e} p99 {:.3e} max {:.3e}   ab_min p50 {:.3e} p01 {:.3e}   \
             cap_hits p50 {:.0} p99 {:.0}",
            q(&mut dt.clone(), 0.5), q(&mut dt.clone(), 0.99), q(&mut dt, 1.0),
            q(&mut ab.clone(), 0.5), q(&mut ab, 0.01),
            q(&mut ch.clone(), 0.5), q(&mut ch, 0.99),
        );
    }

    // --- per box -------------------------------------------------------------------------
    println!("\n== PER MARKED BOX ==");
    println!(
        "  Verdicts are `BUG_REPORT.md` §4's, carried so the table cannot be read without them.\n\
         **B4 is the negative control** -- its pale wedge survives refinement and is real\n\
         structure. If the mask fires there too it is reading the chart, not the fault.\n"
    );
    println!(
        "  {:>4} {:>11} {:>10} {:>10} {:>10} {:>10} {:>11} {:>11}",
        "box", "verdict", "floored", "capped", "budget", "err>10", "dt_max p50", "ab_min p50"
    );
    for (nm, u, v, h, verdict) in BOXES.iter().chain(std::iter::once(&(
        "FRAME", 0.5, 0.5, 0.5, "baseline",
    ))) {
        let lo = |c: f64| (((c - h) * res as f64).floor().max(0.0)) as usize;
        let hi = |c: f64| (((c + h) * res as f64).ceil().min(res as f64)) as usize;
        // PNG row 0 is the MINIMUM v -- `Slice::axis` runs low-to-high with index and
        // `save_rect` writes rows in buffer order, so a fraction from the top maps straight
        // through with no flip. (A harness that flipped it rendered half a frame away.)
        let idx: Vec<usize> = (lo(*v)..hi(*v))
            .flat_map(|r| (lo(*u)..hi(*u)).map(move |c| r * res + c))
            .collect();
        if idx.is_empty() {
            continue;
        }
        let f = |m: &[bool]| idx.iter().filter(|&&i| m[i]).count() as f64 / idx.len() as f64;
        let mut dt: Vec<f64> =
            idx.iter().map(|&i| px[i].dt_max).filter(|x| x.is_finite()).collect();
        let mut ab: Vec<f64> =
            idx.iter().map(|&i| px[i].ab_min).filter(|x| x.is_finite()).collect();
        println!(
            "  {nm:>4} {verdict:>11} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>11.3e} {:>11.3e}",
            f(&floored), f(&capped), f(&budget), f(&bad),
            q(&mut dt, 0.5), q(&mut ab, 0.5)
        );
    }

    // --- the images ----------------------------------------------------------------------
    // Every panel carries its settings. See `output::provenance_sidecar`.
    let extra = format!(
        "res={res}x{res}\ncase=config_stability\nt_max={t_max}\nn_sync={n_sync}\n\
         r_coll={r_coll}\nwindow=({cx},{cy},{half})\n"
    );
    let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);
    let m_here = grid::decode_state(&chart, 0, cx, cy).m;
    let sites = colour::landmarks(&m_here);
    let base: Vec<u8> = px
        .iter()
        .flat_map(|p| colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi))
        .collect();
    let _ = adaptive::save_rect(&format!("{dir}/config_stability_uniform.png"), res, res, &base);
    let _ = prin_rs::output::provenance_sidecar(
        &format!("{dir}/config_stability_uniform.png"), &ens, &extra,
    );

    // Each mask twice: standalone binary, and overlaid on the panel it is meant to explain.
    // The overlay is what answers "do these coincide"; the binary is what says what the mask
    // actually is, without the panel's own structure suggesting a shape to the eye.
    for (nm, m) in [
        ("floored", &floored), ("capped", &capped), ("budget", &budget),
        ("union", &union), ("errhot", &bad),
    ] {
        let mut bin = Vec::with_capacity(res * res * 3);
        let mut ov = base.clone();
        for (i, &on) in m.iter().enumerate() {
            let c: [u8; 3] = if on { [255, 255, 255] } else { [0, 0, 0] };
            bin.extend_from_slice(&c);
            if on {
                ov[3 * i] = 0;
                ov[3 * i + 1] = 255;
                ov[3 * i + 2] = 255;
            }
        }
        let _ = adaptive::save_rect(&format!("{dir}/mask_{nm}.png"), res, res, &bin);
        let _ = adaptive::save_rect(&format!("{dir}/overlay_{nm}.png"), res, res, &ov);
        for f in [format!("{dir}/mask_{nm}.png"), format!("{dir}/overlay_{nm}.png")] {
            let _ = prin_rs::output::provenance_sidecar(&f, &ens, &extra);
        }
    }

    // `dt_max` as a field, log-ramped on a FIXED window: auto-ranging it would stretch whatever
    // range happened to be present to full scale, which on a question about whether a cliff
    // exists manufactures the answer.
    const DLO: f64 = 1e-6;
    const DHI: f64 = 1e0;
    let mut dbuf = Vec::with_capacity(res * res * 3);
    for p in &px {
        let c = if !p.dt_max.is_finite() || p.dt_max <= 0.0 {
            [255, 0, 255]
        } else {
            let x = (p.dt_max.log10() - DLO.log10()) / (DHI.log10() - DLO.log10());
            ramp(x)
        };
        dbuf.extend_from_slice(&c);
    }
    let _ = adaptive::save_rect(&format!("{dir}/field_dt_max.png"), res, res, &dbuf);
    let _ = prin_rs::output::provenance_sidecar(&format!("{dir}/field_dt_max.png"), &ens, &extra);

    println!(
        "\nWrote {dir}/: config_stability_uniform.png, mask_*.png, overlay_*.png (cyan),\n\
         field_dt_max.png (inferno, FIXED window {DLO:e}..{DHI:e}, magenta = non-finite).\n\n\
         HOW TO READ IT. If the saturated mask coincides with the cut-out regions the mechanism\n\
         is confirmed in one image. If the lift is ~1 and the overlay looks alike anyway, the\n\
         resemblance is the base rate and the mask is not the cause -- say so and keep the\n\
         plumbing, because an advance-anyway site with no telemetry is a gap either way."
    );
}
