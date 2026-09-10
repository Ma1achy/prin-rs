//! The parity run: CPU-f64, CPU-f32 and the GPU, over the ceil-boundary sweep.
//!
//! Every arm is reported, including the negative controls. A zero in a control column is a FAILED
//! MEASUREMENT, not a clean result -- the harness would be insensitive and the clean arms' agreement
//! would mean nothing.

use lowering_spike::*;

/// CPU-f64: the same decisions computed in double precision, then compared as integers. This is the
/// Tier-B question in its sharpest form -- the words must match the f32 arms bit for bit.
fn eval_f64(d2: f32, thr_f64: &[f64; N_MAX as usize], r_coll: f64) -> Discrete {
    let d2 = d2 as f64;
    let mut n: u32 = 1;
    for i in 0..N_MAX as usize {
        if d2 < thr_f64[i] { n = i as u32 + 2; }
    }
    let bf = if n > N_MAX { N_MAX } else { n };

    let c = (R_SUB / d2.sqrt()).powf(GAMMA).ceil().clamp(1.0, N_MAX as f64);
    let bt = c as u32;

    let packed = pack(bf % 6, bf & 3);
    let packed_ctl = pack(bf % 6, bt & 3);
    Discrete {
        bucket_frozen: bf,
        bucket_trans: bt,
        coll_sq: (d2 < r_coll * r_coll) as u32,
        coll_sqrt: (d2.sqrt() < r_coll) as u32,
        packed,
        roundtrip: from_bits(packed),
        pairsum_flag: pairsum_flag(bf),
        pairsum_brk: pairsum_brk(bf),
        packed_ctl,
        roundtrip_ctl: from_bits(packed_ctl),
    }
}

const FIELDS: [&str; 10] = [
    "bucket_frozen", "bucket_trans", "coll_sq", "coll_sqrt",
    "packed", "roundtrip", "pairsum_flag", "pairsum_brk", "packed_ctl", "roundtrip_ctl",
];
const RULE: [&str; 10] = [
    "clean", "CONTROL", "clean", "CONTROL", "clean", "clean", "clean", "CONTROL",
    "propagate", "propagate",
];

fn field(d: &Discrete, i: usize) -> u32 {
    match i {
        0 => d.bucket_frozen, 1 => d.bucket_trans, 2 => d.coll_sq, 3 => d.coll_sqrt,
        4 => d.packed, 5 => d.roundtrip, 6 => d.pairsum_flag, 7 => d.pairsum_brk,
        8 => d.packed_ctl, _ => d.roundtrip_ctl,
    }
}

fn main() {
    let thr = thresholds();
    // The f64 arm's thresholds are the SAME BITS widened, not recomputed in f64 -- recomputing them
    // would make the table itself a second source of disagreement and confound the measurement.
    let mut thr_f64 = [0f64; N_MAX as usize];
    for i in 0..N_MAX as usize { thr_f64[i] = thr[i] as f64; }

    let d2s = boundary_sweep();
    let r_coll = R_COLL;
    println!("sweep: {} states, +/-4 ulp around each of the {} ceil boundaries", d2s.len(), N_MAX);

    let cpu32: Vec<Discrete> = d2s.iter().map(|&d| eval_cpu(d, &thr, r_coll)).collect();
    let cpu64: Vec<Discrete> = d2s.iter().map(|&d| eval_f64(d, &thr_f64, r_coll as f64)).collect();

    let gpu = gpu::run(&d2s, &thr);
    match &gpu {
        Some(g) => println!("gpu adapter: backend={} name={:?} driver={:?}", g.backend, g.adapter, g.driver),
        None => println!("gpu adapter: NONE -- no arm to compare, this run is inconclusive"),
    }

    println!("\n{:<15} {:<8} {:>12} {:>12}", "field", "rule", "f32-vs-f64", "gpu-vs-f32");
    let mut clean_forks = 0usize;
    let mut control_forks = 0usize;
    let mut per_control = Vec::new();
    for i in 0..10 {
        let a = (0..d2s.len()).filter(|&k| field(&cpu32[k], i) != field(&cpu64[k], i)).count();
        let b = match &gpu {
            Some(g) => (0..d2s.len()).filter(|&k| field(&cpu32[k], i) != field(&g.out[k], i)).count() as i64,
            None => -1,
        };
        let bs = if b < 0 { "n/a".into() } else { format!("{}", b) };
        println!("{:<15} {:<8} {:>12} {:>12}", FIELDS[i], RULE[i], a, bs);
        let tot = a + b.max(0) as usize;
        match RULE[i] {
            "clean" => clean_forks += tot,
            "CONTROL" => { control_forks += tot; per_control.push((FIELDS[i], tot)); }
            _ => {}
        }
    }

    // The loop-shape control is silent only if the shapes both compute the RIGHT answer. Check it
    // against the closed form: agreement on a wrong value would be the miscompile, everywhere.
    let mut wrong_flag = 0usize;
    let mut wrong_brk = 0usize;
    if let Some(g) = &gpu {
        for k in 0..d2s.len() {
            let want = pairsum_expected(cpu32[k].bucket_frozen);
            if g.out[k].pairsum_flag != want { wrong_flag += 1 }
            if g.out[k].pairsum_brk != want { wrong_brk += 1 }
        }
        println!(
            "\nloop shapes against the closed form, on the GPU: flag-tested {} wrong, break-form {} wrong of {}",
            wrong_flag, wrong_brk, d2s.len()
        );
    }

    // How near does the sweep bring sqrt(d2) to r_coll? WGSL specifies sqrt to 1 ULP rather than
    // correctly-rounded, so it CAN fork -- but only on a state within a ULP of the threshold. If
    // the sweep never gets there, a zero says nothing about sqrt.
    let mut within = 0usize;
    let mut nearest = f32::INFINITY;
    for &d2 in &d2s {
        let r = d2.sqrt();
        let ulps = ((r.to_bits() as i64) - (r_coll.to_bits() as i64)).abs();
        if ulps <= 1 { within += 1 }
        nearest = nearest.min(ulps as f32);
    }
    println!(
        "coll_sqrt sensitivity: {} of {} states put sqrt(d2) within 1 ULP of r_coll (nearest {} ULP)",
        within, d2s.len(), nearest as i64
    );

    println!("\ncontrols, individually:");
    for (f, t) in &per_control {
        println!("  {:<15} {:>5} forks {}", f, t, if *t == 0 { "<- SILENT: see the notes above before reading it as clean" } else { "" });
    }
    println!("\nclean arms:   {} forks", clean_forks);
    println!("control arms: {} forks", control_forks);
    println!(
        "\nverdict: {}",
        if control_forks == 0 {
            "INCONCLUSIVE -- the controls did not fork, so the harness is not sensitive enough \
             and the clean arms' agreement is not evidence."
        } else if clean_forks == 0 {
            "PASS -- the discrete surface survives lowering bitwise, and the controls forked, \
             so the harness could have seen it if it had not."
        } else {
            "FAIL -- a comparison-only arm forked."
        }
    );
}
