//! **The four step-control candidates, as pictures.**
//!
//! `results/step_control/README.md` has the numbers; this is what they look like. One directory
//! per region, one panel set per mode, and a `README.md` written alongside so the folder explains
//! itself without the table.
//!
//! # The one thing that would invalidate this
//!
//! **Every colour window is computed once, from the BASELINE, and reused for every mode.** An
//! auto-ranged ramp per panel would stretch each mode's own range to full scale, which on a
//! before/after comparison manufactures or hides the very thing being shown — this project has
//! been caught by that twice already. The `spread` window comes from the baseline's own p1–p99;
//! the drift ramp is a fixed constant shared by the whole gallery; the `error_ratio` mask is a
//! fixed threshold of 10, the project's own flag for *this pixel is not data*.
//!
//! # Four panels, and why the diagnostic one is the point
//!
//! ```text
//!   _uniform    the shipping bivariate colouring -- what actually gets displayed
//!   _outcome    terminal class -- categorical, so a colouring artefact cannot draw it
//!   _drift      energy drift on an inferno ramp -- THE DIAGNOSTIC FIELD
//!   _errhot     error_ratio > 10, binary -- the artefact itself
//! ```
//!
//! *When a numerical defect is suspected, render the diagnostic field, not the science field*:
//! the science panels show a defect only after it has propagated into a spread or a label, while
//! `_drift` and `_errhot` show it at source. The white wedges in `_errhot` **are** the pale
//! regions that started this investigation.
//!
//! # What is deliberately NOT rendered
//!
//! A per-pixel difference panel against the baseline. Over `t = 50` any change of step size gives
//! a different trajectory through a chaotic region — the measured label-flip rate against the
//! `eta/256` ground truth is **1.0000 for a correct mode and a broken one alike**. Such a panel
//! would be a picture of chaos wearing the caption of a comparison.
//!
//! Args: `res root`. Writes `<root>/step_control/gallery/<region>/`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart, Slice};
use prin_rs::integrate::az::StepLimit;
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour, png};

const WINDOW: f64 = 0.4;
/// The drift ramp, **fixed** and shared by every panel in the gallery. See the module docs.
const DLO: f64 = 1e-12;
const DHI: f64 = 1e2;

/// `(directory name, label, mode, f)`. The parameters are the ones the measurement chose.
const MODES: [(&str, &str, StepLimit, f64); 5] = [
    ("0_baseline", "None (baseline)", StepLimit::None, 0.0),
    ("1_predictive", "B  Predictive f=0.02  <- SHIPPED", StepLimit::Predictive, 0.02),
    ("2_reject", "A  Reject f=0.02", StepLimit::Reject, 0.02),
    ("3_abgrowth", "C  AbGrowth f=2  (bitwise inert)", StepLimit::AbGrowth, 2.0),
    ("4_global", "D  Global f=0.25  (the dumb control)", StepLimit::Global, 0.25),
];

fn arg<T: std::str::FromStr>(i: usize, d: T) -> T {
    std::env::args().nth(i).and_then(|s| s.parse().ok()).unwrap_or(d)
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

fn targets(n: usize) -> Vec<(&'static str, Slice, EnsembleCfg, [f64; 3])> {
    let mut out = Vec::new();
    let (chart, cx, cy, half) = Chart::config_stability();
    out.push((
        "config_stability",
        grid::Slice::body_plane(n, n, cx, cy, half, 0).with_chart(chart),
        EnsembleCfg {
            refine_flagged: false,
            t_max: 50.0,
            n_sync: (50.0f64 / WINDOW).round() as usize,
            r_coll_frac: 0.005,
            escape_rule: EscapeRule::Closure(CLOSURE_TAU),
            closure_k: 1,
            stop_on_escape: false,
            ..Default::default()
        },
        [cx, cy, half],
    ));
    if let Some(c) = grid::gallery_cases().into_iter().find(|c| c.0 == "preset_shape") {
        out.push((
            "preset_shape",
            grid::Slice::body_plane(n, n, c.2, c.3, c.4, 0).with_chart(c.1),
            EnsembleCfg { refine_flagged: false, ..Default::default() },
            [c.2, c.3, c.4],
        ));
    }
    out.push((
        "deep_interior",
        grid::region("deep interior", n, n, 0.05).unwrap().with_chart(Chart::BodyPlane),
        EnsembleCfg { refine_flagged: false, ..Default::default() },
        [0.0, 0.0, 0.05],
    ));
    out
}

fn main() {
    let res: usize = arg(1, 384);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let base_dir = format!("{root}/step_control/gallery");

    println!(
        "STEP-CONTROL GALLERY, {res}^2. Every colour window is taken from the BASELINE and\n\
         reused for every mode -- an auto-ranged ramp per panel would manufacture the\n\
         difference it is meant to show.\n"
    );

    let mut index: Vec<String> = Vec::new();

    for (name, sl, base, win) in targets(res) {
        let dir = format!("{base_dir}/{name}");
        let _ = std::fs::create_dir_all(&dir);
        println!("\n================ {name} ================");
        println!("config: {}\n", base.provenance());
        println!(
            "  {:>34} {:>9} {:>11} {:>9} {:>11} {:>10} {:>9}",
            "mode", "secs", "err p99", "err>10", "steps p50", "overshoot", "nonfin"
        );

        // The shared windows, from the baseline. Computed before the loop on purpose.
        let mut shared: Option<(f64, f64, colour::SiteSet)> = None;
        index.push(format!("\n## {name}\n"));

        for (slug, label, mode, f) in MODES {
            let cfg = EnsembleCfg { step_limit: mode, step_limit_f: f, ..base };
            let t0 = std::time::Instant::now();
            let px: Vec<PixelOut> = (0..sl.npix())
                .into_par_iter()
                .map(|k| pixel::evaluate::<f64>(&sl, k, &cfg))
                .collect();
            let secs = t0.elapsed().as_secs_f64();

            let m_here = grid::decode_state(&sl.chart, 0, win[0], win[1]).m;
            let sites = colour::landmarks(&m_here);
            if shared.is_none() {
                let (lo, hi) = colour::range(&px, Scalar::ShapeSpread);
                shared = Some((lo, hi, sites.clone()));
            }
            let (lo, hi, ref sites) = *shared.as_ref().unwrap();

            let mut uni = Vec::with_capacity(px.len() * 3);
            let mut out = Vec::with_capacity(px.len() * 3);
            let mut dft = Vec::with_capacity(px.len() * 3);
            let mut hot = Vec::with_capacity(px.len() * 3);
            for p in &px {
                uni.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, sites, lo, hi));
                out.extend_from_slice(&png::outcome_rgb(p));
                dft.extend_from_slice(&if p.energy_drift_max.is_finite()
                    && p.energy_drift_max > 0.0
                {
                    let x = (p.energy_drift_max.log10() - DLO.log10())
                        / (DHI.log10() - DLO.log10());
                    ramp(x)
                } else if p.energy_drift_max == 0.0 {
                    ramp(0.0)
                } else {
                    [255, 0, 255]
                });
                // Binary and a FIXED threshold, so the four panels are directly comparable.
                hot.extend_from_slice(&if p.error_ratio > 10.0 {
                    [255, 255, 255]
                } else {
                    [12, 12, 16]
                });
            }
            let extra = format!(
                "res={res}x{res}\ncase={name}\nmode={label}\nstep_limit={mode:?}\n\
                 step_limit_f={f}\nwindow=({},{},{})\n\
                 spread ramp=({lo:.6e},{hi:.6e})  <- taken from the BASELINE and shared\n\
                 drift ramp=({DLO:e},{DHI:e})  <- fixed constant\n\
                 errhot threshold=10  <- fixed\n",
                win[0], win[1], win[2]
            );
            for (suffix, buf) in
                [("uniform", &uni), ("outcome", &out), ("drift", &dft), ("errhot", &hot)]
            {
                let path = format!("{dir}/{slug}_{suffix}.png");
                let _ = adaptive::save_rect(&path, res, res, buf);
                let _ = prin_rs::output::provenance_sidecar(&path, &cfg, &extra);
            }

            let n = px.len() as f64;
            let hot_frac = px.iter().filter(|p| p.error_ratio > 10.0).count() as f64 / n;
            let mut st: Vec<f64> = px.iter().map(|p| p.total_substeps as f64).collect();
            st.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let mut er: Vec<f64> =
                px.iter().map(|p| p.error_ratio).filter(|x| x.is_finite()).collect();
            er.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let ov: u64 = px.iter().map(|p| p.n_overshoot).sum();
            let nf = px.iter().filter(|p| p.n_nonfinite > 0).count();
            println!(
                "  {label:>34} {secs:>9.1} {:>11.3e} {hot_frac:>9.4} {:>11.3e} {ov:>10} {nf:>9}",
                er[((er.len() - 1) as f64 * 0.99).round() as usize],
                st[st.len() / 2],
            );
            index.push(format!(
                "| `{slug}` | {label} | {:.3e} | {hot_frac:.4} | {:.3e} | {ov} |",
                er[((er.len() - 1) as f64 * 0.99).round() as usize],
                st[st.len() / 2]
            ));
        }
    }

    // The folder explains itself.
    let mut md = String::from(
        "# Step-control gallery\n\n\
         Generated by `cargo run --release --example step_limit_gallery <res> <root>`.\n\
         Numbers and method: `../README.md`.\n\n\
         One directory per region; five modes per region; four panels per mode.\n\n\
         ```text\n  _uniform    the shipping bivariate colouring -- what gets displayed\n  \
         _outcome    terminal class, categorical\n  _drift      energy drift, inferno, \
         THE DIAGNOSTIC FIELD\n  _errhot     error_ratio > 10, binary -- the artefact itself\n\
         ```\n\n\
         **Every colour window is taken from the baseline and reused for every mode.** The spread\n\
         ramp is the baseline's own p1-p99, the drift ramp is a fixed constant, and the `errhot`\n\
         threshold is a fixed 10. An auto-ranged ramp per panel would stretch each mode's own\n\
         range to full scale and manufacture the difference it is meant to show. Each panel has a\n\
         `.cfg.txt` sidecar naming its full config.\n\n\
         **`_errhot` is the panel to look at first.** Its white wedges are the pale regions that\n\
         started this investigation, and `1_predictive` is black.\n\n\
         **No difference panel is rendered.** Over `t = 50` any change of step size gives a\n\
         different trajectory through a chaotic region -- the label-flip rate against the\n\
         `eta/256` ground truth is 1.0000 for a correct mode and a broken one alike, so such a\n\
         panel would be a picture of chaos with the caption of a comparison.\n\n\
         | dir | mode | err p99 | err>10 | steps p50 | overshoot |\n|---|---|---|---|---|---|\n",
    );
    for line in &index {
        if line.starts_with("\n## ") {
            md.push_str(&format!(
                "\n{}\n| dir | mode | err p99 | err>10 | steps p50 | overshoot |\n|---|---|---|---|---|---|\n",
                line.trim()
            ));
        } else {
            md.push_str(line);
            md.push('\n');
        }
    }
    let _ = std::fs::write(format!("{base_dir}/README.md"), md);
    println!("\nWrote {base_dir}/ -- one directory per region, plus README.md and per-panel sidecars.");
}
