//! Gate (c): cross-check against the Python reference at f64.
//!
//! `#[ignore]`d because it needs python3 + numpy on PATH. Run with:
//!     cargo test --release --test xcheck -- --ignored --nocapture

use std::process::Command;

fn cargo() -> String {
    std::env::var("CARGO").unwrap_or_else(|_| {
        format!("{}/.cargo/bin/cargo", std::env::var("HOME").unwrap())
    })
}

fn run_case(name: &str, rel_tol: &str, abs_tol: &str) -> bool {
    let refp = format!("xcheck_out/ref_{name}.tsv");
    let rsp = format!("xcheck_out/rs_{name}.tsv");

    let ok = Command::new("python3")
        .args(["tools/xcheck/dump_ref.py", "--case", name, "--out", &refp])
        .status()
        .expect("python3 not available")
        .success();
    assert!(ok, "reference dump failed for {name}");

    let ok = Command::new(cargo())
        .args(["run", "--release", "--quiet", "--bin", "xcheck", "--", "--case", name, "--out", &rsp])
        .status()
        .expect("cargo not available")
        .success();
    assert!(ok, "rust dump failed for {name}");

    let out = Command::new("python3")
        .args(["tools/xcheck/compare.py", &refp, &rsp, "--tol", rel_tol, "--atol", abs_tol])
        .output()
        .expect("compare failed");
    print!("{}", String::from_utf8_lossy(&out.stdout));
    out.status.success()
}

/// The algebra case: no integration, no chaos, nowhere to hide.
#[test]
#[ignore]
fn algebra_matches_the_reference() {
    assert!(run_case("algebra", "1e-15", "1e-15"), "algebra cross-check failed");
}

/// Short horizons, where chaotic amplification is negligible. A transcription error cannot
/// hide here — and in fact `t = 0.5` and `t = 1.0` come out **bitwise identical**.
#[test]
#[ignore]
fn short_horizons_match_to_1e13() {
    for c in ["az_t0p5", "az_t1", "az_t2"] {
        println!("\n=== {c} ===");
        assert!(run_case(c, "1e-13", "1e-13"), "{c} cross-check failed");
    }
}

/// The §9 horizon.
///
/// **This now passes the brief's ~1e-10 comfortably**, at 2.7e-13, with both sides running
/// the conditioned inverse LC branch. In PR #2 it sat at 1.9e-10 and I proposed amending §9
/// on the grounds that two f64 implementations cannot do better through 13 time units of
/// chaos. That was premature: the shortfall was not chaotic amplification of ulp noise, it
/// was branch-cut error injected at every one of the 32 registrations. Conditioning both
/// sides removed three orders.
///
/// The divergence-vs-horizon table is still the better evidence than any single threshold —
/// growth from an O(1e-16) intercept is what says the port is right — but §9's number stands
/// as written and needs no amendment.
#[test]
#[ignore]
fn long_horizon_t13() {
    println!("\n=== az_t13 ===");
    let ok = run_case("az_t13", "1e-7", "1e-11");
    println!("\nJudged on absolute agreement; relative figures are inflated where a");
    println!("coordinate passes near zero. See tools/xcheck/horizon.py for the growth curve.");
    assert!(ok, "az_t13 cross-check failed");
}

/// Benettin FTLE and the diffusion regression, against `tb_ftle.integrate_full`.
///
/// This is the pair with a reference: `tb_ftle.py` sits on `tb.py`'s fixed-step leapfrog, not
/// on Aarseth-Zare, and there is **no AZ+shadow reference** anywhere. So this validates the
/// estimator, and carrying a shadow through AZ is a separate step validated against this one
/// where both resolve.
///
/// The perturbation direction is pinned analytically on both sides. numpy's Ziggurat is not
/// ported and reproducing it is not required — the direction only seeds the shadow, and a
/// comparison needs both sides to use the same one. The Python side substitutes it into the
/// reference without reimplementing any of the reference's algebra.
#[test]
#[ignore]
fn ftle_matches_the_reference() {
    assert!(run_case("ftle", "1e-10", "1e-12"), "ftle cross-check failed");
}
