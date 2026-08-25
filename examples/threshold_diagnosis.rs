//! **The threshold, measured against the distribution it is supposed to cut.** No integration —
//! this reads the committed `.prnq` dumps, which carry every quad's spread and decision already.
//!
//! # What it says
//!
//! `tau_display = 1e-4` in every committed run. Against the pooled leaf-spread distribution that
//! is the **0.4th percentile**, so the split predicate is true for essentially every quad. The
//! criterion is an always-split rule wearing a threshold's clothing, and every spatial field
//! built on the hot mask is saturated for the same reason.
//!
//! # The part that was missed, and it is the argument for rank
//!
//! A fixed threshold fails on **both** sides, and the corpus contains both:
//!
//! - `tau` **below** the bulk -> everything splits -> uniform at **max** depth. Most charts.
//! - `tau` **above** the bulk -> everything keeps -> uniform at **depth 2**, 16 leaves against a
//!   complete 4096. `far`, `deep interior`, and every deep zoom step: **16 of the 18 trees the
//!   veto does not bind** are stopped this way, with leaf-spread medians running 9.45e-5 down to
//!   4.26e-8 against `tau = 1e-4`, and their leaf decisions are `keep` almost to the last one.
//!
//! **Selectivity requires the threshold to cut through the bulk**, and the chart's own dynamic
//! range decides whether any fixed value can. A fixed threshold therefore has a narrow window of
//! usefulness that varies per chart. **A ranking cannot land above or below a distribution; it
//! always cuts through it.** That is a stronger argument for rank than the treadmill one.
//!
//! # `preset_shape` is a THIRD mode, and it is not the one it looks like
//!
//! It was tempting to file `preset_shape` under the upper-side failure — 16 leaves, depth
//! variance 0, the widest dynamic range among the charts. **The decision column says otherwise.**
//! Its leaf-spread median is `2.86e-1`, three and a half thousand times *above* `tau`: it clears
//! the spread gate everywhere and is stopped by **`alpha`** — 8 `floor` (below `alpha_lo`) and 8
//! `keep` (between the thresholds). Not one leaf failed the spread test.
//!
//! So `preset_shape` is the only tree in the corpus where the **alpha gate** is what is being
//! exercised, which makes it the cleanest instance of the standing result that `alpha_hi` does
//! more work than the criterion. Quoting it as a `tau` failure would have been a mechanism read
//! off a shape.
//!
//! # And the stop-reason column is the headline
//!
//! `Decision::MaxRelDepth` — a **camera veto** — stops 95%+ of leaves on 23 of 26 charts. So the
//! observed uniformity is not the criterion saying "split" and being obeyed; it is the criterion
//! never saying **stop**, with something else terminating the descent. `preset_shape` is the one
//! chart where the veto does not bind at all, which makes it the only chart in the set where the
//! criterion's own decisions are visible.

use std::collections::BTreeMap;

fn read(path: &str) -> Option<(String, Vec<String>, Vec<Vec<f64>>)> {
    let d = std::fs::read(path).ok()?;
    if d.len() < 12 || &d[..4] != b"PRNQ" {
        return None;
    }
    let hl = u32::from_le_bytes(d[8..12].try_into().ok()?) as usize;
    let hdr = String::from_utf8_lossy(&d[12..12 + hl]).into_owned();
    let mut off = 12 + hl;
    let n = u64::from_le_bytes(d[off..off + 8].try_into().ok()?) as usize;
    off += 8;
    let nf = u32::from_le_bytes(d[off..off + 4].try_into().ok()?) as usize;
    off += 4;
    let fields: Vec<String> = hdr
        .lines()
        .find(|l| l.starts_with("fields="))?
        .trim_start_matches("fields=")
        .split(',')
        .map(str::to_string)
        .collect();
    let rows = (0..n)
        .map(|i| {
            let base = off + i * nf * 8;
            (0..nf)
                .map(|k| f64::from_le_bytes(d[base + k * 8..base + k * 8 + 8].try_into().unwrap()))
                .collect()
        })
        .collect();
    Some((hdr, fields, rows))
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let k = ((v.len() - 1) as f64 * p).round() as usize;
    v[k]
}

/// Spearman by rank correlation, NaN on a constant input rather than 0 or 1.
fn spearman(x: &[f64], y: &[f64]) -> f64 {
    let rank = |v: &[f64]| {
        let mut idx: Vec<usize> = (0..v.len()).collect();
        idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap());
        let mut r = vec![0.0f64; v.len()];
        let mut i = 0;
        while i < idx.len() {
            let mut j = i;
            while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
                j += 1;
            }
            let avg = (i + j) as f64 / 2.0 + 1.0;
            for &k in &idx[i..=j] {
                r[k] = avg;
            }
            i = j + 1;
        }
        r
    };
    let (a, b) = (rank(x), rank(y));
    let n = a.len() as f64;
    if n < 3.0 {
        return f64::NAN;
    }
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for k in 0..a.len() {
        let (da, db) = (a[k] - ma, b[k] - mb);
        sab += da * db;
        saa += da * da;
        sbb += db * db;
    }
    if saa == 0.0 || sbb == 0.0 {
        f64::NAN
    } else {
        sab / (saa * sbb).sqrt()
    }
}

struct Row {
    name: String,
    leaves: usize,
    pct_max: f64,
    depth_var: f64,
    dyn_range: f64,
    spread_med: f64,
    veto_pct: f64,
    dec: [usize; 11],
}

fn dec_name(c: usize) -> &'static str {
    match c {
        0 => "pending", 1 => "split", 2 => "floor", 3 => "keep", 4 => "prec_floor",
        5 => "max_level", 6 => "budget", 7 => "screen", 8 => "max_rel_depth", 9 => "collapsed",
        _ => "?",
    }
}

fn main() {
    let dirs = ["results/charts", "results/vertical", "results/criterion"];
    let mut pooled: Vec<f64> = Vec::new();
    let mut hot_sat = (0usize, 0usize);
    let mut comp_one = (0usize, 0usize);
    let mut rows: Vec<Row> = Vec::new();
    let mut n_dumps = 0usize;
    let mut taus: BTreeMap<String, usize> = BTreeMap::new();

    for dir in dirs {
        let Ok(rd) = std::fs::read_dir(dir) else { continue };
        let mut paths: Vec<String> = rd
            .filter_map(|e| e.ok())
            .map(|e| e.path().to_string_lossy().into_owned())
            .filter(|p| p.ends_with(".prnq"))
            .collect();
        paths.sort();
        for p in paths {
            let Some((hdr, fields, recs)) = read(&p) else { continue };
            n_dumps += 1;
            let col = |k: &str| fields.iter().position(|f| f == k);
            let (Some(ci_leaf), Some(ci_lv), Some(ci_sp)) =
                (col("is_leaf"), col("level"), col("spread_median"))
            else {
                continue;
            };
            for tok in hdr.replace('\n', " ").split_whitespace() {
                if let Some(v) = tok.strip_prefix("tau_display=") {
                    *taus.entry(v.to_string()).or_default() += 1;
                }
            }
            let leaves: Vec<&Vec<f64>> = recs.iter().filter(|r| r[ci_leaf] > 0.5).collect();
            if leaves.is_empty() {
                continue;
            }
            let lv: Vec<f64> = leaves.iter().map(|r| r[ci_lv]).collect();
            let mx = lv.iter().cloned().fold(f64::MIN, f64::max);
            let mean = lv.iter().sum::<f64>() / lv.len() as f64;
            let var = lv.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / lv.len() as f64;
            let mut sp: Vec<f64> =
                leaves.iter().map(|r| r[ci_sp]).filter(|x| x.is_finite() && *x > 0.0).collect();
            pooled.extend(sp.iter().cloned());

            if let (Some(ci_h), Some(ci_c), Some(ci_nf)) =
                (col("n_hot_within"), col("n_components_within"), col("n_footprints"))
            {
                for r in &leaves {
                    hot_sat.1 += 1;
                    comp_one.1 += 1;
                    if r[ci_h] == r[ci_nf] {
                        hot_sat.0 += 1;
                    }
                    if r[ci_c] <= 1.0 {
                        comp_one.0 += 1;
                    }
                }
            }

            let veto = col("decision")
                .map(|c| {
                    leaves.iter().filter(|r| r[c] == 7.0 || r[c] == 8.0).count() as f64
                        / leaves.len() as f64
                })
                .unwrap_or(f64::NAN);

            let mut dec = [0usize; 11];
            if let Some(c) = col("decision") {
                for r in &leaves {
                    let k = r[c] as usize;
                    if k < 11 {
                        dec[k] += 1;
                    }
                }
            }
            let stem = p.rsplit('/').next().unwrap().trim_end_matches(".prnq").to_string();
            let dr = if sp.len() > 10 { q(&mut sp, 0.99) / q(&mut sp, 0.01) } else { f64::NAN };
            rows.push(Row {
                name: format!("{}/{stem}", dir.rsplit('/').next().unwrap()),
                leaves: leaves.len(),
                pct_max: 100.0 * lv.iter().filter(|&&x| x == mx).count() as f64 / lv.len() as f64,
                depth_var: var,
                dyn_range: dr,
                spread_med: q(&mut sp, 0.5),
                veto_pct: 100.0 * veto,
                dec,
            });
        }
    }

    println!("{n_dumps} committed .prnq dumps, {} leaves pooled", pooled.len());
    println!("tau_display values in the corpus: {taus:?}");
    println!();

    // ---- 1. where the threshold sits in its own distribution -------------------------
    let mut p = pooled.clone();
    println!("pooled leaf spread: p1 {:.3e}  median {:.3e}  p99 {:.3e}",
             q(&mut p, 0.01), q(&mut p, 0.5), q(&mut p, 0.99));
    println!();
    println!("  {:>8}  {:>10}  {:>12}", "tau", "% exceeding", "percentile");
    for t in [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1] {
        let above = pooled.iter().filter(|&&x| x > t).count() as f64 / pooled.len() as f64;
        println!("  {:>8.0e}  {:>10.1}  {:>12.1}", t, 100.0 * above, 100.0 * (1.0 - above));
    }
    println!("  ^ the shipped value is the first row. A predicate true for 99.6% of quads is not");
    println!("    selecting; it is an always-split rule. The interesting range starts at 1e-3.");
    println!();

    // ---- 2. the saturation this causes ----------------------------------------------
    // **The denominator is printed with the statistic, not once at the top.** Before this PR the
    // corpus was mixed-version -- `vertical/` was still PRNQ v1 and carries none of the hot-mask
    // columns -- so this ran over a strict subset of the leaves the pooled spread above used, and
    // two numbers side by side had different denominators without saying so.
    println!("hot mask, ABSOLUTE rule, over {} of {} pooled leaves:", hot_sat.1, pooled.len());
    println!("  n_hot == n_footprints : {:.1}%", 100.0 * hot_sat.0 as f64 / hot_sat.1 as f64);
    println!("  n_components <= 1     : {:.1}%", 100.0 * comp_one.0 as f64 / comp_one.1 as f64);
    println!("  ^ one blob covering the whole quad, everywhere. The spatial fields cannot");
    println!("    discriminate as thresholded, so measuring them before this is fixed returns");
    println!("    a null that was guaranteed.");
    println!();

    // ---- 3. the two-sided failure ----------------------------------------------------
    rows.sort_by(|a, b| b.dyn_range.partial_cmp(&a.dyn_range).unwrap_or(std::cmp::Ordering::Equal));
    println!("{:<34} {:>7} {:>7} {:>8} {:>10} {:>10} {:>7}",
             "dump", "leaves", "%max", "depthvar", "p99/p1", "spread med", "veto%");
    for r in &rows {
        println!("{:<34} {:>7} {:>7.1} {:>8.3} {:>10.1} {:>10.2e} {:>7.1}",
                 r.name, r.leaves, r.pct_max, r.depth_var, r.dyn_range, r.spread_med, r.veto_pct);
    }
    println!();

    let dr: Vec<f64> = rows.iter().filter(|r| r.dyn_range.is_finite()).map(|r| r.dyn_range).collect();
    let dv: Vec<f64> = rows.iter().filter(|r| r.dyn_range.is_finite()).map(|r| r.depth_var).collect();
    let pm: Vec<f64> = rows.iter().filter(|r| r.dyn_range.is_finite()).map(|r| r.pct_max).collect();
    println!("spearman(dynamic range, depth variance) = {:+.3}  over {} dumps", spearman(&dr, &dv), dr.len());
    println!("spearman(dynamic range, % at max depth) = {:+.3}", spearman(&dr, &pm));
    println!();
    println!("THE TWO-SIDED FAILURE. tau = 1e-4 sits below the bulk on almost every chart, so");
    println!("everything splits and the tree is uniform AT MAX DEPTH. Where the bulk sits BELOW");
    println!("tau instead -- far, deep interior, every deep zoom step -- everything keeps and the");
    println!("tree is uniform AT DEPTH 2, 16 leaves against a complete 4096. Same failure, both");
    println!("sides: a fixed level that does not cut through the distribution. A ranking cannot");
    println!("land above or below one. That is the argument for rank.");
    println!();
    println!("preset_shape is a THIRD mode and not the one it looks like: its spread median is");
    println!("3400x ABOVE tau, so it clears the spread gate on every leaf and is stopped by ALPHA");
    println!("(8 floor + 8 keep). It is the only tree here exercising the alpha gate. See the");
    println!("decision breakdown below before attributing any of these to tau.");
    println!();
    let bound: Vec<&Row> = rows.iter().filter(|r| r.veto_pct >= 95.0).collect();
    println!("AND THE STOP REASON. {} of {} dumps have >=95% of leaves stopped by a CAMERA VETO",
             bound.len(), rows.len());
    println!("(ScreenFloor or MaxRelDepth). On those the criterion decides almost nothing, so the");
    println!("uniformity is the criterion never saying STOP while something else terminates the");
    println!("descent -- not the criterion saying SPLIT and being obeyed. Never quote a leaf count");
    println!("without its stop-reason breakdown.");
    let free: Vec<&Row> = rows.iter().filter(|r| r.veto_pct < 5.0).collect();
    println!();
    println!("Dumps where the veto binds on <5% of leaves -- the only trees that are entirely");
    println!("their own decisions. WHICH GATE stopped each, from the dump's own decision column:");
    println!("  {:<34} {:>7} {:>10}  {}", "dump", "leaves", "spread med", "leaf decisions");
    for r in &free {
        let br: Vec<String> = (0..11)
            .filter(|&k| r.dec[k] > 0)
            .map(|k| format!("{}={}", dec_name(k), r.dec[k]))
            .collect();
        println!("  {:<34} {:>7} {:>10.2e}  {}", r.name, r.leaves, r.spread_med, br.join(" "));
    }
    println!();
    println!("  tau_display = 1e-4. A `keep` leaf failed the SPREAD gate (spread <= tau); a `floor`");
    println!("  leaf passed it and failed the ALPHA gate. Read the spread median against tau to see");
    println!("  which side of the distribution the threshold landed on for that dump.");
}
