//! **Does Heggie remove the wedges?** The claim the port was FOR, on the object that was circled.
//!
//! # What a wedge is, from the record and not from a summary
//!
//! `INVESTIGATION.md` §0 and `results/circled/README.md`: on
//! `config_stability_stop0_uniform.png` at 1024^2, **pale, low-chroma patches with straight,
//! sharp edges**, interrupting otherwise continuous ribbon structure, with magenta speckle at
//! their cores. The pale is `spread_shape` saturating because the ensemble copies diverged to
//! garbage; the magenta is non-finite pixels.
//!
//! **The magenta is already gone.** The `dtau` fix took non-finite `30109 -> 178` and the
//! predictive step limit took `err>10` `0.1110 -> 0.0000`, both before Heggie existed. What is
//! left to ask about is the **pale saturated region with the straight edge**.
//!
//! # And it is NOT the reference-body argmax
//!
//! `INVESTIGATION.md` hypothesis 8 -- *"they are the reference-body argmax"* -- is **REFUTED four
//! ways**: switches *depleted* at 0.648, transverse to the wedges, and drift moves `7.5e-5`
//! decades under hysteresis. What survived is the **re-registration count**:
//!
//! ```text
//!   LC branch unconditioned          2.5e-6 decades
//!   hysteresis, switches 17.8 -> 6.7 7.5e-5 decades
//!   re-registration x2, fixed step   4.4e-1 decades   <- 6000x / 175000x
//! ```
//!
//! > It is not *which* chart is chosen. It is *how often* the state is passed through one.
//!
//! Heggie has no reference body and therefore **no re-registration at all**, which is why §5
//! calls it "the *specific* cure rather than a general one".
//!
//! # The mask is the committed one, not a new one
//!
//! `examples/circled_ics.rs:152` selects pale as OKLab `l > 0.86 && chroma < 0.045` on the
//! rendered panel, and its `dense` rule at `:288` keeps pixels with **>=25% pale in a 9x9** --
//! because the raw pale mask "includes the fine striation everywhere, which is not what was
//! circled". Both are reused verbatim. Inventing a mask here would make every number
//! incomparable with the investigation that defined the object.
//!
//! **The colour window is taken from the AZ arm and SHARED**, as the gallery does, or "pale"
//! would mean a different thing per arm and the census would score the auto-range.
//!
//! # The control that decides whether any of this means anything
//!
//! `az_prefix` runs AZ with the step-control fixes **off** -- `FixedPerInterval`, no landing
//! clamp, no step limit -- which is the kernel the wedges were first seen under. **If it does not
//! show a large dense component with a straight edge, the metric cannot see the artefact and
//! every other row is void.** A null against a metric that never had a subject is the failure
//! this project catalogues most.
//!
//! ```text
//! cargo run --release --example wedge_census -- [res] [out_dir]
//! ```
use rayon::prelude::*;

use prin_rs::ensemble::pixel::{self, EnsembleCfg, PixelOut};
use prin_rs::ensemble::provenance::Override;
use prin_rs::grid::{self, Chart};
use prin_rs::integrate::az::{DtauMode, StepLimit};
use prin_rs::integrate::Integrator;
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour, oklab};
use prin_rs::spatial;

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

/// `circled_ics.rs:152`, verbatim: pale is light and low-chroma on the rendered panel.
fn pale_mask(rgb: &[u8], n: usize) -> Vec<bool> {
    (0..n * n)
        .map(|i| {
            let lab = oklab::srgb_to_oklab([rgb[3 * i], rgb[3 * i + 1], rgb[3 * i + 2]]);
            lab[0] > 0.86 && lab[1].hypot(lab[2]) < 0.045
        })
        .collect()
}

/// **The same object taken from the FIELD instead of the render.**
///
/// The committed mask is defined on the rendered panel, which is right for comparability with
/// what was circled and carries two things the field does not: the pale test is a *compound* of
/// lightness (from `spread_shape`) and chroma (from the shape vector's own magnitude), and the
/// sRGB round trip quantises to 8 bits before the threshold is applied.
///
/// So `spread_shape` is masked directly, at the **shared upper end of the AZ colour window** --
/// which is exactly what "the ramp is saturated" means, and is shared across arms for the same
/// reason the window is. Neither mask is the primary: **they are each other's control.** If they
/// agree the finding is about the field; if they disagree the render is adding or hiding
/// something and that is worth knowing before any row is quoted. It is the discipline
/// `circled_ics.rs` used on its own hand-drawn and property-selected masks.
fn saturated_mask(px: &[PixelOut], cut: f64) -> Vec<bool> {
    px.iter().map(|p| p.spread_shape.is_finite() && p.spread_shape >= cut).collect()
}

/// The record's own number: *"the pale patches are `spread_shape` **saturating at 0.39**"*.
///
/// **Absolute, not a quantile.** `colour::range` is p1..p99, so a mask at its top selects 1% of
/// the frame *by construction* whatever the field does -- a quantile rule makes the count a
/// non-signal, which this project already has on record for `n_hot`. An absolute cut can read
/// zero, and reading zero is then a fact about the field rather than about the rule.
const SATURATED: f64 = 0.39;

/// `circled_ics.rs:288`, verbatim: **>=25% pale in a 9x9**. The raw pale mask carries the fine
/// striation everywhere, which is not what was circled; the circled features are the dense,
/// contiguous patches of it.
fn dense_mask(pale: &[bool], n: usize, half: usize) -> Vec<bool> {
    let w = 2 * half + 1;
    let area = w * w;
    let mut d = vec![false; n * n];
    for y in half..n.saturating_sub(half) {
        for x in half..n.saturating_sub(half) {
            if !pale[y * n + x] {
                continue;
            }
            let mut c = 0usize;
            for dy in 0..w {
                for dx in 0..w {
                    if pale[(y + dy - half) * n + (x + dx - half)] {
                        c += 1;
                    }
                }
            }
            d[y * n + x] = c * 4 > area;
        }
    }
    d
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(512);
    let dir: String = std::env::args().nth(2).unwrap_or_else(|| "results/wedge".into());
    let _ = std::fs::create_dir_all(&dir);

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0],
        z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0;
    q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (cx, cy, half) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);

    let sl = grid::Slice::body_plane(res, res, cx, cy, half, 0).with_chart(chart);
    let m_here = grid::decode_state(&chart, 0, cx, cy).m;
    let sites = colour::landmarks(&m_here);

    // `radius` scales with the raster: the hand fit ran over 621 rows of a 1024^2 panel.
    let radius = (res / 32).max(4);
    // **The dense window must cover the same AREA of the slice at any raster.** `circled_ics.rs`
    // used 9x9 at 1024^2, and the dense-pale wedges are 1935 px of 1048576 -- 0.18% of the frame.
    // Held at 9x9 on a 96^2 raster that window spans 11x more of the slice and the rule selects
    // nothing, which is exactly what the first run of this census did: `az_prefix`, the subject
    // control, read zero dense and the whole table was void.
    let dense_half = ((4 * res) / 1024).max(1);

    // `integrators` (default) compares occupants; `ablation` turns the three step-control fixes
    // off ONE AT A TIME. The first census turned all three off together, so it established *the
    // step-control work* removed the wedges and could not say which of the three did -- stated
    // as a limitation at the time rather than discovered later. This is that question.
    let mode: String = std::env::args().nth(3).unwrap_or_else(|| "integrators".into());
    let off_dtau = Override::DtauMode(DtauMode::FixedPerInterval);
    let off_clamp = Override::ClampFinalStep(false);
    let off_limit = Override::StepLimit(StepLimit::None);
    let az = Override::Integrator(Integrator::Az);
    let (arms, order): (Vec<(&str, Vec<Override>)>, Vec<&str>) = if mode == "ablation" {
        (
            vec![
                // `all` runs first so it sets the shared colour window, as `az` does above.
                ("all", vec![az.clone()]),
                ("none", vec![az.clone(), off_dtau.clone(), off_clamp.clone(), off_limit.clone()]),
                ("dtau_only", vec![az.clone(), off_clamp.clone(), off_limit.clone()]),
                ("clamp_only", vec![az.clone(), off_dtau.clone(), off_limit.clone()]),
                ("limit_only", vec![az.clone(), off_dtau.clone(), off_clamp.clone()]),
            ],
            vec!["all", "none", "dtau_only", "clamp_only", "limit_only"],
        )
    } else {
        (
            vec![
                ("az", vec![az.clone()]),
                ("az_prefix", vec![az.clone(), off_dtau, off_clamp, off_limit]),
                ("heggie", vec![Override::Integrator(Integrator::Heggie)]),
                ("logh_rk4", vec![Override::Integrator(Integrator::LogHRk4)]),
                ("plain_rk4", vec![Override::Integrator(Integrator::PlainRk4)]),
            ],
            vec!["az", "az_prefix", "heggie", "logh_rk4", "plain_rk4"],
        )
    };

    println!("WEDGE CENSUS on config_stability at {res}^2, t_max = 50, science pass.");
    println!("Mask is `circled_ics.rs`'s own: OKLab l > 0.86 && chroma < 0.045, then >=25% pale");
    println!("in a {}x{} window -- 9x9 at 1024^2, scaled to hold the same AREA of the slice.", 2 * dense_half + 1, 2 * dense_half + 1);
    println!("Colour window taken from the `az` arm and SHARED. The field mask is an ABSOLUTE");
    println!("`spread_shape >= {SATURATED}`, the value the record names, never a quantile.");
    println!("Boundary straightness is LOCAL at radius {radius}: 0.0 is a straight edge, ~0.7 a");
    println!("jagged one (tests/straightness.rs measures 0.00000 and 0.71743 at matched area).");
    println!();
    println!("{:>28}{:>27} | {:>27} |", "", "-- RENDERED PANEL --", "-- spread_shape FIELD --");
    println!(
        "{:>11} {:>7} {:>8} {:>8} {:>8} {:>8} | {:>8} {:>8} {:>8} | {:>6} {:>9}",
        "arm", "nonfin", "pale", "dense", "lgst", "bnd str", "dense", "lgst", "bnd str",
        "jacc", "shp p99"
    );

    let mut window: Option<(f64, f64)> = None;
    let mut rows: Vec<(String, String)> = Vec::new();

    for want in order.iter().copied() {
        let (label, ov) = arms.iter().find(|(l, _)| *l == want).unwrap();
        let mut o = ov.clone();
        o.push(Override::TMax(50.0));
        o.push(Override::NSync(125));
        o.push(Override::RefineFlagged(false));
        o.push(Override::MaxSteps(2_000_000));
        let cfg = EnsembleCfg::production().with_overrides(&o);

        let px: Vec<PixelOut> =
            (0..sl.npix()).into_par_iter().map(|k| pixel::evaluate::<f64>(&sl, k, &cfg)).collect();
        let (lo, hi) = *window.get_or_insert_with(|| colour::range(&px, Scalar::ShapeSpread));

        let mut rgb = Vec::with_capacity(px.len() * 3);
        for p in &px {
            rgb.extend_from_slice(&colour::rgb(p, Scalar::ShapeSpread, &sites, lo, hi));
        }
        let _ = adaptive::save_rect(&format!("{dir}/{label}_uniform.png"), res, res, &rgb);

        // Two masks of the same object, each other's control.
        let measure = |m: &[bool]| {
            let dense = dense_mask(m, res, dense_half);
            let comps = spatial::components(&dense, res);
            let lgst = comps.first().map(|c| c.len()).unwrap_or(0);
            let mut lmask = vec![false; res * res];
            if let Some(c) = comps.first() {
                for &(x, y) in c {
                    lmask[y * res + x] = true;
                }
            }
            let b = spatial::boundary_straightness(&lmask, res, radius, 40);
            (dense.iter().filter(|&&x| x).count(), lgst, b)
        };
        let pale = pale_mask(&rgb, res);
        let sat = saturated_mask(&px, SATURATED);
        let (dense_n, lgst, bstr) = measure(&pale);
        let (fdense_n, flgst, fbstr) = measure(&sat);
        // Agreement between the two masks, as Jaccard over the raw (undensified) sets.
        let (mut inter, mut uni) = (0usize, 0usize);
        for i in 0..res * res {
            if pale[i] || sat[i] {
                uni += 1;
                if pale[i] && sat[i] {
                    inter += 1;
                }
            }
        }
        let jac = if uni == 0 { f64::NAN } else { inter as f64 / uni as f64 };
        let nonfin = px.iter().filter(|p| !p.spread_shape.is_finite()).count();
        let mut ss: Vec<f64> = px.iter().map(|p| p.spread_shape).filter(|x| x.is_finite()).collect();
        ss.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let p99 = if ss.is_empty() { f64::NAN } else { ss[(ss.len() - 1) * 99 / 100] };

        let f = |v: usize| v as f64 / (res * res) as f64;
        eprintln!("  [{label}] dense {:.4} lgst {} bnd {:.4} | field dense {:.4} lgst {} bnd {:.4}",
                  f(dense_n), lgst, bstr, f(fdense_n), flgst, fbstr);
        rows.push((
            label.to_string(),
            format!(
                "{:>11} {:>7} {:>8.4} {:>8.4} {:>8} {:>8.4} | {:>8.4} {:>8} {:>8.4} | {:>6.3} {:>9.4}",
                label,
                nonfin,
                f(pale.iter().filter(|&&b| b).count()),
                f(dense_n),
                lgst,
                bstr,
                f(fdense_n),
                flgst,
                fbstr,
                jac,
                p99
            ),
        ));
    }
    // Reporting order puts the all-off arm first: it is the subject control and nothing below it
    // means anything until it separates.
    let report: Vec<&str> = if mode == "ablation" {
        vec!["none", "dtau_only", "clamp_only", "limit_only", "all"]
    } else {
        vec!["az_prefix", "az", "heggie", "logh_rk4", "plain_rk4"]
    };
    for want in report {
        if let Some((_, line)) = rows.iter().find(|(l, _)| l == want) {
            println!("{line}");
        }
    }

    println!();
    println!("HOW TO READ THIS");
    println!();
    println!("**`az_prefix` FIRST.** It is AZ with the step-control fixes off -- the kernel the");
    println!("wedges were first seen under. If it does not carry a large dense component with a");
    println!("STRAIGHT boundary, this metric cannot see the artefact and every row below it is");
    println!("void. A null from a metric that never had a subject is this project's most");
    println!("frequently catalogued failure.");
    println!();
    println!("**Then `az` against `heggie`.** `az` is the post-fix kernel, so what it still");
    println!("carries is the RESIDUAL the step-control work did not remove -- which is what the");
    println!("re-registration mechanism predicts Heggie should not have.");
    println!();
    println!("**The two mask columns are each other's control.** The rendered one is the mask");
    println!("`circled_ics.rs` defined and is comparable with what was circled; the field one is");
    println!("`spread_shape` above the shared AZ window's top, which is what saturation means,");
    println!("with no 8-bit round trip and no chroma term. `jacc` is their agreement. If they");
    println!("disagree, the render is adding or hiding something and no row should be quoted.");
    println!();
    println!("**`plain_rk4` is the no-reference-body control.** It has no chart to re-register");
    println!("into, so it cannot have the artefact -- but its field is wrecked, which is what");
    println!("says a low `bnd str` means straight and not merely empty. Read `dense` beside it.");
}
