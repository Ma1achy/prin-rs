//! Rust side of the cross-check. Emits the same TSV as `tools/xcheck/dump_ref.py`.
//!
//!     cargo run --release --bin xcheck -- --case algebra --out xcheck_out/rs_algebra.tsv
//!
//! Case parameters are duplicated here rather than read from the Python side on purpose:
//! `compare.py` asserts the emitted headers match, so a drift between the two definitions
//! is caught by the comparison instead of silently agreeing on the wrong thing.

use std::fs::File;
use std::io::{BufWriter, Write};

use prin_rs::physics::{burrau, energy, newton, shape};
use prin_rs::rng::SplitMix64;
use prin_rs::Vec2;

const ALGEBRA_N: usize = 256;
const ALGEBRA_SEED: u64 = 20260822;
const ALGEBRA_COLUMNS: &str = "idx,energy,d01,d02,d12,inertia,hyperradius,n0,n1,n2";

/// Mirrors `cases.random_configs`. Positions in [-4,4), velocities in [-1,1); deliberately
/// unconstrained, so the algebra is exercised over a wide spread of geometries including
/// near-degenerate ones.
fn random_configs(n: usize, seed: u64) -> Vec<([Vec2<f64>; 3], [Vec2<f64>; 3])> {
    let mut rng = SplitMix64::new(seed);
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let mut r = [Vec2::zero(); 3];
        for k in 0..3 {
            r[k] = Vec2::new(rng.range(-4.0, 4.0), rng.range(-4.0, 4.0));
        }
        let mut v = [Vec2::zero(); 3];
        for k in 0..3 {
            v[k] = Vec2::new(rng.range(-1.0, 1.0), rng.range(-1.0, 1.0));
        }
        out.push((r, v));
    }
    out
}

fn dump_algebra(path: &str) -> std::io::Result<()> {
    let m = burrau::masses::<f64>();
    let cfgs = random_configs(ALGEBRA_N, ALGEBRA_SEED);

    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "# case=algebra kind=algebra n={ALGEBRA_N} seed={ALGEBRA_SEED}"
    )?;
    writeln!(w, "# masses=3,4,5 G=1 eps2=0")?;
    writeln!(w, "# columns={ALGEBRA_COLUMNS}")?;

    for (i, (r, v)) in cfgs.iter().enumerate() {
        let e = energy::energy(r, v, &m, 0.0);
        let d = newton::pair_dists(r);
        let ine = energy::inertia(r, &m);
        let hyp = energy::hyperradius(r, &m);
        let n = shape::shape_vec(r, &m);
        let row = [e, d[0], d[1], d[2], ine, hyp, n[0], n[1], n[2]];
        write!(w, "{i}")?;
        for x in row {
            write!(w, "\t{x:.17e}")?;
        }
        writeln!(w)?;
    }
    w.flush()?;
    eprintln!("wrote {path}: {ALGEBRA_N} rows");
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |flag: &str| -> Option<String> {
        args.iter().position(|a| a == flag).and_then(|i| args.get(i + 1)).cloned()
    };
    let case = get("--case").unwrap_or_else(|| {
        eprintln!("usage: xcheck --case <name> --out <path>");
        std::process::exit(2);
    });
    let out = get("--out").unwrap_or_else(|| {
        eprintln!("usage: xcheck --case <name> --out <path>");
        std::process::exit(2);
    });
    if let Some(dir) = std::path::Path::new(&out).parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let res = match case.as_str() {
        "algebra" => dump_algebra(&out),
        other => {
            eprintln!("unknown case: {other}");
            std::process::exit(2);
        }
    };
    if let Err(e) = res {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}
