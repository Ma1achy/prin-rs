//! **The gallery table, re-derived from the committed `.prnq` dumps** — no integration.
//!
//! It exists because `chart_gallery`'s own printed table carried a `bound` column that read
//! `crit` for every row whose budget was not exhausted, and that hid the largest fact about the
//! set: on most charts the tree is stopped by `Decision::MaxRelDepth`, a **camera veto**, and the
//! leaf count is a fact about the cap rather than about the criterion. A row stopped by a veto
//! has not exercised the criterion at all, and its `alpha` describes quads the cap forced rather
//! than quads the criterion chose. That is the screen-floor lesson at a second stop condition.
//!
//! Rather than pay a three-hour regeneration for a label, the dumps are the authority: they carry
//! `decision` per quad and always did. This reads them and prints what the run actually did.
//!
//! It also diagnoses the **mechanism test** — leaf depth against `terminated_fraction`. A Spearman
//! over that pair has three distinct ways to be uninformative and they are not distinguishable
//! from the number alone:
//!
//!   1. **x constant** — every leaf at one depth, so there is no depth axis. `spearman` is NaN.
//!   2. **y constant** — `terminated_fraction` takes one distinct value. `spearman` is NaN.
//!   3. **y saturated** — one value holds most of the leaves, so the correlation is read off the
//!      thin remainder. `spearman` is a finite number and means very little.
//!
//! Count the distinct values before reading the curve. That is the project's own standing rule
//! and it applies here exactly.

use std::collections::HashMap;

fn dec_name(c: u8) -> &'static str {
    match c {
        0 => "pending", 1 => "split", 2 => "floor", 3 => "keep", 4 => "prec_floor",
        5 => "max_level", 6 => "budget", 7 => "screen", 8 => "max_rel_depth", 9 => "collapsed",
        _ => "?",
    }
}

fn read(path: &str) -> Option<(Vec<String>, Vec<Vec<f64>>)> {
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
    let mut rows = Vec::with_capacity(n);
    for i in 0..n {
        let base = off + i * nf * 8;
        rows.push(
            (0..nf)
                .map(|k| f64::from_le_bytes(d[base + k * 8..base + k * 8 + 8].try_into().unwrap()))
                .collect(),
        );
    }
    Some((fields, rows))
}

const ORDER: &[&str] = &[
    "body_plane", "plane_00deg", "shape_sphere",
    "latent_shape", "latent_inner_p", "latent_outer_p", "latent_mass", "latent_mixed",
    "latent_oblique_a", "latent_oblique_b", "burrau_nu_k", "invariant_lz_k", "mass_simplex",
    "preset_shape", "preset_prho", "preset_plambda", "preset_shape_pl",
    "preset_shape_h1", "preset_prho_h1", "preset_plambda_h1", "preset_shape_pl_h1",
    "latent_shape_h3", "latent_inner_p_h3", "latent_outer_p_h3", "latent_mass_h3",
    "latent_mixed_h3",
];

fn main() {
    let dir = std::env::args().nth(1).unwrap_or_else(|| "results/charts".into());

    println!(
        "What actually stopped each descent, read from the committed .prnq dumps.\n\
         `max_rel_depth` and `screen` are CAMERA VETOES -- a quad stopped by one has not\n\
         exercised the criterion, and the row's leaf count and alpha are facts about the cap.\n\
         `floor` (spread below tau) and `keep` (alpha says refinement does not pay) are the\n\
         criterion's own two answers.\n"
    );
    println!(
        "{:>20} {:>7} {:>7} {:>8} {:>8} {:>8} {:>8}   stopped by",
        "case", "leaves", "depths", "veto%", "floor%", "keep%", "other%"
    );

    let mut rows_out: Vec<(String, usize, f64, f64)> = Vec::new();
    for case in ORDER {
        let path = format!("{dir}/{case}.prnq");
        let Some((f, rows)) = read(&path) else {
            println!("{case:>20}   MISSING ({path})");
            continue;
        };
        let idx = |name: &str| f.iter().position(|x| x == name).expect("field");
        let (i_leaf, i_lv, i_dec) = (idx("is_leaf"), idx("level"), idx("decision"));
        let (i_tf, i_ef) = (idx("terminated_fraction"), idx("escape_fraction"));

        let leaves: Vec<&Vec<f64>> = rows.iter().filter(|r| r[i_leaf] == 1.0).collect();
        let n = leaves.len().max(1);
        let mut by_dec: HashMap<u8, usize> = HashMap::new();
        for r in &leaves {
            *by_dec.entry(r[i_dec] as u8).or_default() += 1;
        }
        let pc = |c: u8| 100.0 * *by_dec.get(&c).unwrap_or(&0) as f64 / n as f64;
        let (veto, fl, kp) = (pc(8) + pc(7), pc(2), pc(3));
        // Clamp: the three arms are exhaustive, so any residual is float noise and
        // printing it as `-0.0%` reads as a bug rather than as rounding.
        let other = (100.0 - veto - fl - kp).max(0.0);
        let mut top: Vec<(u8, usize)> = by_dec.into_iter().collect();
        top.sort_by_key(|&(_, k)| std::cmp::Reverse(k));
        let told: Vec<String> = top
            .iter()
            .take(3)
            .map(|&(c, k)| format!("{}={}", dec_name(c), k))
            .collect();
        println!(
            "{case:>20} {:>7} {:>7} {veto:>7.1}% {fl:>7.1}% {kp:>7.1}% {other:>7.1}%   {}",
            leaves.len(),
            {
                let mut d: Vec<u32> = leaves.iter().map(|r| r[i_lv] as u32).collect();
                d.sort_unstable();
                d.dedup();
                d.len()
            },
            told.join(", ")
        );

        // The mechanism-test diagnosis.
        let depths = {
            let mut d: Vec<u32> = leaves.iter().map(|r| r[i_lv] as u32).collect();
            d.sort_unstable();
            d.dedup();
            d.len()
        };
        let tf: Vec<f64> = leaves.iter().map(|r| r[i_tf]).collect();
        let mut keys: Vec<u64> = tf.iter().map(|x| x.to_bits()).collect();
        keys.sort_unstable();
        let distinct = {
            let mut k = keys.clone();
            k.dedup();
            k.len()
        };
        let modal = {
            let mut best = 0usize;
            let (mut i, mut run) = (0usize, 0usize);
            while i < keys.len() {
                run = if i > 0 && keys[i] == keys[i - 1] { run + 1 } else { 1 };
                best = best.max(run);
                i += 1;
            }
            100.0 * best as f64 / n as f64
        };
        let ef: f64 = leaves.iter().map(|r| r[i_ef]).sum::<f64>() / n as f64;
        rows_out.push((case.to_string(), distinct, modal, ef));
        let _ = depths;
    }

    println!(
        "\n\nThe mechanism test -- leaf depth against `terminated_fraction` -- and whether it can\n\
         be read at all. Three ways for it to be uninformative, and the Spearman alone cannot\n\
         tell them apart:\n\
         \x20 x CONSTANT   every leaf at one depth; there is no depth axis\n\
         \x20 y CONSTANT   terminated_fraction takes one value\n\
         \x20 y SATURATED  one value holds >90% of leaves; the correlation reads the remainder\n"
    );
    println!(
        "{:>20} {:>9} {:>9} {:>9}   verdict",
        "case", "tf values", "modal%", "mean esc"
    );
    for (case, distinct, modal, ef) in &rows_out {
        let path = format!("{dir}/{case}.prnq");
        let Some((f, rows)) = read(&path) else { continue };
        let idx = |name: &str| f.iter().position(|x| x == name).unwrap();
        let (i_leaf, i_lv) = (idx("is_leaf"), idx("level"));
        let leaves: Vec<&Vec<f64>> = rows.iter().filter(|r| r[i_leaf] == 1.0).collect();
        let mut d: Vec<u32> = leaves.iter().map(|r| r[i_lv] as u32).collect();
        d.sort_unstable();
        d.dedup();
        let mut why: Vec<&str> = Vec::new();
        if d.len() < 2 {
            why.push("x CONSTANT");
        }
        if *distinct < 2 {
            why.push("y CONSTANT");
        }
        if *modal > 90.0 && *distinct >= 2 {
            why.push("y SATURATED");
        }
        println!(
            "{case:>20} {distinct:>9} {modal:>8.1}% {ef:>9.4}   {}",
            if why.is_empty() { "READABLE".to_string() } else { why.join(" + ") }
        );
    }
}
