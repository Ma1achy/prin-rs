//! Rust side of the cross-check. Emits the same TSV as `tools/xcheck/dump_ref.py`.
//!
//!     cargo run --release --bin xcheck -- --case algebra --out xcheck_out/rs_algebra.tsv
//!
//! Case parameters are duplicated here rather than read from the Python side on purpose:
//! `compare.py` asserts the emitted headers match, so a drift between the two definitions
//! is caught by the comparison instead of silently agreeing on the wrong thing.

use std::fs::File;
use std::io::{BufWriter, Write};

use prin_rs::grid::Slice;
use prin_rs::integrate::az::{self, DtauMode};
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

/// The step-sizing mode, from `PRIN_DTAU_MODE`, defaulting to the shipped one.
///
/// **Both sides carry the same `dtau` defect and both were fixed, so the comparison that means
/// anything is old-against-old AND new-against-new** -- either alone would pass while the two
/// transcriptions diverged. It is an environment variable rather than a flag so
/// `tests/xcheck.rs` needs no argument threading: the test spawns both processes and they
/// inherit it. It is also emitted into the header, and `compare.py` asserts headers match, so a
/// mode mismatch between the two sides is a hard failure rather than a quiet disagreement.
fn dtau_mode_from_env() -> (DtauMode, &'static str) {
    match std::env::var("PRIN_DTAU_MODE").unwrap_or_default().as_str() {
        "fixed" => (DtauMode::FixedPerInterval, "fixed"),
        "per-step-remaining" => (DtauMode::PerStepRemaining, "per-step-remaining"),
        "" | "per-step-interval" => (DtauMode::PerStepInterval, "per-step-interval"),
        other => panic!("PRIN_DTAU_MODE: expected fixed|per-step-remaining|per-step-interval, got {other}"),
    }
}

/// Fixed physical length of one sync sub-interval; `n_sync` is derived from it so every
/// horizon in the sweep runs at the same per-interval resolution. Mirrors `cases.py`.
const SYNC_INTERVAL: f64 = 13.0 / 32.0;

fn n_sync_for(t_max: f64) -> usize {
    ((t_max / SYNC_INTERVAL).round() as usize).max(1)
}

/// Nominal copies only — `ens=0`, no RNG on either side. See `dump_ref.py:dump_az`.
fn dump_az(name: &str, t_max: f64, path: &str, lc_stable: bool) -> std::io::Result<()> {
    let m = burrau::masses::<f64>();
    let n_sync = n_sync_for(t_max);
    let eta = 0.01f64;
    let max_steps = 30_000usize;
    let (dtau_mode, mode_name) = dtau_mode_from_env();
    let s = Slice::body_plane(3, 3, 1.0, 3.0, 0.05, 0);

    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "# case={name} kind=az nx={} ny={} cx={:?} cy={:?} half={:?} body={}",
        s.nx, s.ny, s.cx, s.cy, s.half, s.body
    )?;
    writeln!(
        w,
        "# t_max={t_max:?} n_sync={n_sync} eta={eta:?} max_steps={max_steps} masses=3,4,5 G=1 ens=0 dtau_mode={mode_name}"
    )?;
    let refs: Vec<String> = (0..n_sync).map(|k| format!("ref{k:02}")).collect();
    writeln!(
        w,
        "# columns=idx,r0x,r0y,r1x,r1y,r2x,r2y,v0x,v0y,v1x,v1y,v2x,v2y,t,dmin_ref,drift,switches,{}",
        refs.join(",")
    )?;

    for i in 0..s.npix() {
        let o = az::integrate_az_opts(
            s.nominal::<f64>(i), &m, t_max, n_sync, eta, max_steps,
            &az::reference_opts(None, lc_stable, dtau_mode),
        );
        write!(w, "{i}")?;
        for k in 0..3 {
            write!(w, "\t{:.17e}\t{:.17e}", o.state.r[k].x, o.state.r[k].y)?;
        }
        for k in 0..3 {
            write!(w, "\t{:.17e}\t{:.17e}", o.state.v[k].x, o.state.v[k].y)?;
        }
        write!(w, "\t{:.17e}\t{:.17e}\t{:.17e}\t{:.17e}", o.t, o.d_min_ref, o.drift, o.switches as f64)?;
        for k in 0..n_sync {
            write!(w, "\t{:.17e}", o.refs[k] as f64)?;
        }
        writeln!(w)?;
    }
    w.flush()?;
    eprintln!("wrote {path}: {} rows, n_sync={n_sync}", s.npix());
    Ok(())
}

/// Benettin FTLE + the diffusion regression, against `tb_ftle.integrate_full`.
///
/// The perturbation direction is pinned analytically on both sides rather than drawn: numpy's
/// Ziggurat is not ported, reproducing it is not required, and the direction only seeds the
/// shadow. Nothing about an RNG enters the validated path.
///
/// `sample_every` is 25 on both sides because the reference **hardcodes** `s % 25 == 0` and
/// takes no parameter for it. It appears in the header so the two sides are asserted to agree
/// on it rather than coincidentally sharing a literal.
fn dump_ftle(path: &str) -> std::io::Result<()> {
    use prin_rs::physics::ftle;

    const NX: usize = 4;
    const NY: usize = 4;
    const CX: f64 = 1.0;
    const CY: f64 = 3.0;
    const HALF: f64 = 0.05;
    const BODY: usize = 0;
    const T_MAX: f64 = 2.0;
    const DT: f64 = 1e-4;
    const EPS: f64 = 0.03;
    const D0: f64 = 1e-8;
    const RENORM: usize = 200;
    const SAMPLE: usize = 25;
    const PERT: [f64; 6] = [1.0, -2.0, 3.0, -4.0, 5.0, -6.0];

    let m = burrau::masses::<f64>();
    let slice = Slice::body_plane(NX, NY, CX, CY, HALF, BODY);
    let pert = ftle::normalise::<f64>(PERT);
    let opts = ftle::FtleOpts {
        eps2: EPS * EPS,
        d0: D0,
        renorm_every: RENORM,
        sample_every: SAMPLE,
    };

    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "# case=ftle kind=ftle nx={NX} ny={NY} cx={CX:?} cy={CY:?} half={HALF:?} body={BODY}"
    )?;
    writeln!(
        w,
        // `d0` is written as the literal Python renders: `repr(1e-8)` is "1e-08" and Rust's
        // `{:?}` is "1e-8". `compare.py` asserts the headers match, so a formatting difference
        // fails the case exactly as a parameter difference would -- which is the intent, and is
        // why this is pinned rather than the assertion relaxed.
        "# t_max={T_MAX:?} dt={DT:?} eps={EPS:?} d0=1e-08 renorm_every={RENORM} \
sample_every={SAMPLE} masses=3,4,5 G=1 pert={}",
        PERT.iter().map(|x| format!("{x:?}")).collect::<Vec<_>>().join(",")
    )?;
    writeln!(w, "# columns=idx,ftle,diffusion,dmin,nre,S,rx0,ry0")?;

    for i in 0..slice.npix() {
        let s0 = slice.nominal::<f64>(i);
        let o = ftle::integrate_full(s0, &m, T_MAX, DT, &opts, &pert);
        writeln!(
            w,
            "{i}\t{:?}\t{:?}\t{:?}\t{}\t{:?}\t{:?}\t{:?}",
            o.ftle, o.diffusion, o.d_min, o.n_renorm, o.s_accum, o.state.r[0].x, o.state.r[0].y
        )?;
    }
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
    // Both sides now default to the conditioned inverse LC branch: `reference/tb_lc.py`
    // sets `USE_STABLE_LC = True`. Pass --lc-unstable on both sides to compare the
    // historical path instead.
    let lc_stable = !args.iter().any(|a| a == "--lc-unstable");
    if let Some(dir) = std::path::Path::new(&out).parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let res = match case.as_str() {
        "algebra" => dump_algebra(&out),
        "ftle" => dump_ftle(&out),
        "az_t0p5" => dump_az("az_t0p5", 0.5, &out, lc_stable),
        "az_t1" => dump_az("az_t1", 1.0, &out, lc_stable),
        "az_t2" => dump_az("az_t2", 2.0, &out, lc_stable),
        "az_t4" => dump_az("az_t4", 4.0, &out, lc_stable),
        "az_t8" => dump_az("az_t8", 8.0, &out, lc_stable),
        "az_t13" => dump_az("az_t13", 13.0, &out, lc_stable),
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
