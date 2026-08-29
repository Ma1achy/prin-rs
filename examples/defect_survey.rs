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
//! Args: `res root`.

use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart, Slice};
use prin_rs::integrate::az::StepLimit;

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

    let mut n_damaged = 0usize;
    let mut n_other = 0usize;
    let mut lines: Vec<String> = Vec::new();
    for (name, sl) in cases {
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
        let a = stats(&run(StepLimit::None, 0.0));
        let b = stats(&run(StepLimit::Predictive, 0.02));
        // A single flagged pixel of 16384 is 6.1e-5, so the cut is at that scale rather than at
        // exactly zero: "damaged" has to mean more than one pixel or every chaotic slice qualifies.
        let dmg = a.0 > 1e-3 || a.2 > 0;
        let fixed = b.0 <= 1e-3 && b.2 == 0;
        let verdict = if !dmg {
            "clean"
        } else if fixed {
            n_damaged += 1;
            "THIS"
        } else {
            n_other += 1;
            "OTHER"
        };
        let line = format!(
            "  {name:>34} {:>10.4} {:>10.4} {:>11.3e} {:>11.3e} {:>10} {:>10} {verdict:>9}",
            a.0, b.0, a.1, b.1, a.2, b.2
        );
        println!("{line}");
        lines.push(line);
    }

    println!(
        "\n== SUMMARY ==\n\
         damaged under `None` and clean under `Pred`: {n_damaged}\n\
         damaged under both (a DIFFERENT problem):    {n_other}\n\n\
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
