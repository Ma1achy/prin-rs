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
/// **The tolerance here is not the evidence — the growth curve is.** See
/// `tools/xcheck/horizon.py`: `t = 0.5` and `t = 1.0` are bitwise identical, `t = 2` agrees
/// to 8.9e-16, and divergence then grows exponentially. By `t = 13` the two trajectories
/// have separated by ~2e-10 in state through ordinary chaotic amplification of ulp-level
/// differences, which is a property of the system, not of the port.
///
/// Every *derived* quantity inherits that separation. `drift` is the clearest case: it
/// differs by 2.7e-11 absolute, which is 1.4e-2 in relative terms only because `drift`
/// itself is ~1.7e-8. Gating it tightly would be gating the divergence.
///
/// So `atol` is set to 1e-9, roughly 5x the observed state separation, and the honest
/// statement of what passed is printed alongside.
#[test]
#[ignore]
fn long_horizon_t13() {
    println!("\n=== az_t13 ===");
    let ok = run_case("az_t13", "1e-7", "1e-9");
    println!("\nJudged on absolute agreement. The relative figures are inflated where a");
    println!("coordinate passes near zero; max |dr| is the meaningful number, and the");
    println!("divergence-vs-horizon table is the actual evidence of correctness.");
    assert!(ok, "az_t13 cross-check failed");
}
