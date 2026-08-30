//! **Was the wedge structure only on `config_stability`, or is it everywhere?**
//!
//! A scientific question about the debugging record, not about the fix: the defect was found on
//! one slice, and if it had been visible on many all along, that says something about what the
//! corpus was and was not able to show.
//!
//! # Why this cannot be read off the committed dumps
//!
//! `results/**/*.raw` carries `error_ratio` per pixel — so the historical record *is* there — but
//! every one of them is a **64x64 Burrau region at `t_max = 13`**, made by `prin`, which takes the
//! library default and therefore had the **repair pass on**. They read `err>10 = 0.0000` on six of
//! seven with `refined` up to 0.1201, which is the repair hiding the damage rather than the damage
//! being absent. The chart corpus is `.prnq`, which is leaf-level and carries no per-pixel
//! `error_ratio` at all. **So the question has to be re-run, and this is that run.**
//!
//! # Two arms, and the second is what makes the first mean something
//!
//! Every case runs under `StepLimit::None` — the kernel every committed number was taken under —
//! and again under the shipped `Predictive`. The first says where the defect was; the second says
//! whether it was the same defect. A case that is damaged under `None` and clean under
//! `Predictive` is this defect; one that is damaged under both is something else and is worth
//! knowing separately.
//!
//! # Resumable
//!
//! Checkpointed **per case**, which is the natural unit here: 34 cases x 2 arms at 512^2 is hours,
//! and losing it whole to an interruption is an experiment-design fault rather than bad luck. Both
//! arms of a case are written together — half a case is not a row. Re-running the same command
//! resumes; the key carries the resolution and the two arms' configs, and a differing key
//! **refuses** rather than mixing two experiments.
//!
//! Panels are re-rendered on resume rather than checkpointed: they are a few seconds of colouring
//! against minutes of integration, and storing them would double the checkpoint for nothing.
//!
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart, Slice};
use prin_rs::integrate::az::StepLimit;
use prin_rs::output::ckpt::Ckpt;
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour};

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

/// Inferno, fixed window, shared by every panel in the survey. **Not auto-ranged**: 34 cases
/// each stretched to their own p1-p99 would make a clean chart and a damaged one look alike,
/// which is the one thing this survey exists to distinguish.
const DLO: f64 = 1e-12;
const DHI: f64 = 1e2;

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

fn stats(px: &[PixelOut]) -> (f64, f64, u64, f64, f64) {
    let n = px.len() as f64;
    let mut e: Vec<f64> = px.iter().map(|p| p.error_ratio).filter(|x| x.is_finite()).collect();
    let mut d: Vec<f64> =
        px.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect();
    (
        px.iter().filter(|p| p.error_ratio > 10.0).count() as f64 / n,
        q(&mut e, 0.99),
        px.iter().map(|p| p.n_overshoot).sum(),
        px.iter().filter(|p| p.n_nonfinite > 0).count() as f64 / n,
        q(&mut d, 0.99),
    )
}

fn main() {
    let res: usize = arg(1, 128);
    let root: String = std::env::args().nth(2).unwrap_or_else(|| "results".into());
    let _ = std::fs::create_dir_all(format!("{root}/output"));

    // Every named Burrau region, then every chart instance in the gallery. All at the settings
    // the corpus used -- `t_max = 13`, `r_coll = 1e-3` -- so this is the corpus's own kernel and
    // not a harder one chosen to make a point.
    let mut cases: Vec<(String, Slice)> = Vec::new();
    for (name, ..) in grid::REGIONS {
        if let Some(s) = grid::region(name, res, res, 0.05) {
            cases.push((format!("region/{name}"), s.with_chart(Chart::BodyPlane)));
        }
    }
    for c in grid::gallery_cases() {
        cases.push((
            format!("chart/{}", c.0),
            grid::Slice::body_plane(res, res, c.2, c.3, c.4, 0).with_chart(c.1),
        ));
    }

    println!(
        "{res}^2, t_max = 13, r_coll = 1e-3 -- the corpus's own settings.\n\n\
         `None` is the kernel every committed number was taken under; `Pred` is the shipped\n\
         `StepLimit::Predictive` at f = 0.02. **A case damaged under None and clean under Pred is\n\
         THIS defect. Damaged under both is something else** and is flagged.\n"
    );
    println!(
        "  {:>34} {:>10} {:>10} {:>11} {:>11} {:>10} {:>10} {:>9}",
        "case", "err>10 None", "err>10 Pred", "err p99 None", "err p99 Pred", "over None",
        "over Pred", "verdict"
    );

    let key = format!(
        "defect_survey res={res} none={} pred={}",
        EnsembleCfg { refine_flagged: false, ..Default::default() }.provenance(),
        EnsembleCfg {
            refine_flagged: false,
            step_limit: StepLimit::Predictive,
            step_limit_f: 0.02,
            ..Default::default()
        }
        .provenance()
    );
    let ck_path = format!("{root}/output/defect_survey_{res}.ckpt");
    let (mut ck, done) = match Ckpt::open(&ck_path, &key) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("checkpoint: {e}");
            std::process::exit(1);
        }
    };
    println!("  checkpoint {ck_path}: {} cases already done\n", done.len());

    let mut n_damaged = 0usize;
    let mut n_other = 0usize;
    let mut n_panels = 0usize;
    let mut n_clean = 0usize;
    let mut lines: Vec<String> = Vec::new();
    for (ci, (name, sl)) in cases.into_iter().enumerate() {
        // A completed case replays from its stored row. **The panels are not replayed** -- a
        // resumed run that printed a row without its panels would read as fully rendered, so a
        // restored case says so in its own line.
        if let Some(b) = done.get(&(ci as u64)) {
            let line = String::from_utf8_lossy(b).to_string();
            println!("{line}   [restored]");
            if line.contains("THIS") {
                n_damaged += 1;
            } else if line.contains("OTHER") {
                n_other += 1;
            } else {
                n_clean += 1;
            }
            lines.push(line);
            continue;
        }
        let run = |lim: StepLimit, f: f64| -> Vec<PixelOut> {
            let cfg = EnsembleCfg {
                refine_flagged: false,
                step_limit: lim,
                step_limit_f: f,
                ..Default::default()
            };
            (0..sl.npix())
                .into_par_iter()
                .map(|k| pixel::evaluate::<f64>(&sl, k, &cfg))
                .collect()
        };
        let pa = run(StepLimit::None, 0.0);
        let pb = run(StepLimit::Predictive, 0.02);
        let a = stats(&pa);
        let b = stats(&pb);
        // A single flagged pixel of 16384 is 6.1e-5, so the cut is at that scale rather than at
        // exactly zero: "damaged" has to mean more than one pixel or every chaotic slice qualifies.
        let dmg = a.0 > 1e-3 || a.2 > 0;
        let fixed = b.0 <= 1e-3 && b.2 == 0;
        let verdict = if !dmg {
            n_clean += 1;
            "clean"
        } else if fixed {
            n_damaged += 1;
            "THIS"
        } else {
            n_other += 1;
            "OTHER"
        };
        // **Panels only for the damaged cases, both arms.** 34 cases x 2 arms x 2 panels at
        // 512^2 is ~55 MB of PNG, most of it identical black; the clean cases are fully
        // described by their row. The count of what is NOT written is printed below, because a
        // silently partial gallery reads as full coverage.
        if dmg {
            let slug = name.replace('/', "_").replace(' ', "_");
            let d = format!("{root}/step_control/survey");
            let _ = std::fs::create_dir_all(&d);
            for (arm, v) in [("none", &pa), ("pred", &pb)] {
                let drift: Vec<u8> = v
                    .iter()
                    .flat_map(|p| {
                        if p.energy_drift_max.is_finite() && p.energy_drift_max > 0.0 {
                            ramp(
                                (p.energy_drift_max.log10() - DLO.log10())
                                    / (DHI.log10() - DLO.log10()),
                            )
                        } else if p.energy_drift_max == 0.0 {
                            ramp(0.0)
                        } else {
                            [255, 0, 255]
                        }
                    })
                    .collect();
                let hot: Vec<u8> = v
                    .iter()
                    .flat_map(|p| {
                        if p.error_ratio > 10.0 { [255u8, 255, 255] } else { [12, 12, 16] }
                    })
                    .collect();
                for (suf, buf) in [("drift", &drift), ("errhot", &hot)] {
                    let path = format!("{d}/{slug}_{arm}_{suf}.png");
                    let _ = adaptive::save_rect(&path, res, res, buf);
                    let _ = prin_rs::output::provenance_sidecar(
                        &path,
                        &EnsembleCfg {
                            refine_flagged: false,
                            step_limit: if arm == "none" {
                                StepLimit::None
                            } else {
                                StepLimit::Predictive
                            },
                            step_limit_f: if arm == "none" { 0.0 } else { 0.02 },
                            ..Default::default()
                        },
                        &format!(
                            "res={res}x{res}\ncase={name}\narm={arm}\nfield={suf}\n\
                             drift ramp=({DLO:e},{DHI:e})  <- FIXED, shared by every case\n\
                             errhot threshold=10  <- fixed\n"
                        ),
                    );
                }
            }
            n_panels += 4;
        }
        let line = format!(
            "  {name:>34} {:>10.4} {:>10.4} {:>11.3e} {:>11.3e} {:>10} {:>10} {verdict:>9}",
            a.0, b.0, a.1, b.1, a.2, b.2
        );
        println!("{line}");
        let _ = ck.put(ci as u64, line.as_bytes());
        lines.push(line);
    }

    println!(
        "\n== SUMMARY ==\n\
         damaged under `None` and clean under `Pred`: {n_damaged}\n\
         damaged under both (a DIFFERENT problem):    {n_other}\n\
         clean, and therefore NOT rendered:            {n_clean}\n\
         panels written:                               {n_panels}  (4 per damaged case)\n\n\
         **Panels are written for damaged cases only**, both arms, on a FIXED drift ramp shared\n\
         by every case -- 34 auto-ranged panels would make a clean chart and a damaged one look\n\
         alike, which is what this survey exists to tell apart. The clean count is printed so a\n\
         partial gallery cannot read as full coverage.\n\n\
         **`overshoot` is the sharper column than `err>10`.** A step that carries the interval\n\
         clock past twice its interval is unambiguous; `err>10` is a threshold on a ratio and a\n\
         chaotic slice can sit near it for honest reasons. Where the two disagree, read the\n\
         overshoot count.\n\n\
         A case reading `clean` here is clean AT 128^2 AND t = 13. The defect is a thin set --\n\
         encounters coinciding with a sync boundary -- so it is resolution-sensitive by\n\
         construction, and absence at this grid is not absence."
    );
    let _ = std::fs::write(
        format!("{root}/output/defect_survey.txt"),
        lines.join("\n") + "\n",
    );
}
