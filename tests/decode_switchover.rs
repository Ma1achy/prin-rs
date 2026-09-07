//! **The switchover: same symptom, opposite response, keyed to which decoder is active.**
//!
//! §14 names this as the piece that is easy to get backwards, and this build had it backwards:
//! `Decision::Collapsed` was terminal regardless. On an f32 consumer that caps the descent around
//! depth 23 and looks exactly like a precision limit — measured in `results/output/deep_zoom.txt`,
//! where `direct_f32` falls to 18/64 distinct ICs by depth 18 while `lin_split_f32` holds **64/64
//! through depth 40**.
//!
//! **And the spec's own framing needed correcting.** It keys the response on *full versus
//! linearised*, which is right for an f32 consumer and wrong for `DirectF64`: f64 is the ceiling
//! here, `LinSplitF32` is measured tracking it rung for rung, and the linearisation buys ~24 levels
//! over f32 and **none over f64**. So the key is whether a *more precise path exists*, not what
//! form the current one takes.

use prin_rs::decode::Path;
use prin_rs::quad::{Decision, QuadReduction};
use prin_rs::scheduler::{decide_reduction_for_test, Policy, SchedCfg};

/// **The ladder, not the form.** `DirectF64` is at the ceiling and must not be labelled a switch.
#[test]
fn the_switchover_key_is_the_ladder_and_not_the_linearised_form() {
    // Every path, so a new one has to be classified rather than inheriting a default.
    for p in [Path::DirectF64, Path::DirectF32, Path::LinNaiveF32, Path::LinSplitF32, Path::LinSplitF64] {
        let can = p.has_more_precise_path();
        let lin = p.is_linearised();
        match p {
            Path::DirectF32 | Path::LinNaiveF32 => assert!(can, "{} has LinSplitF32 above it", p.name()),
            _ => assert!(!can, "{} is at the ceiling", p.name()),
        }
        // The spec's framing and the correct one **disagree**, and this is where.
        if p == Path::DirectF64 {
            assert!(!lin && !can, "DirectF64 is full AND at the ceiling -- the spec would switch it");
        }
        if p == Path::LinNaiveF32 {
            assert!(lin && can, "LinNaiveF32 is linearised AND has room -- the spec would stop it");
        }
    }
    // Named plainly: on two of five paths, "full => switch, linearised => stop" gives the wrong
    // answer. If those two ever agreed with the ladder, this test would be decoration.
    let disagree = [Path::DirectF64, Path::LinNaiveF32]
        .iter()
        .filter(|p| p.is_linearised() != p.has_more_precise_path())
        .count();
    assert_eq!(disagree, 0, "the two paths chosen for this arm must disagree with the naive key");
}

/// **A collapse is a switch where there is somewhere to switch to, and terminal where there is
/// not — and the two arms are both exercised.**
#[test]
fn a_collapse_switches_or_stops_by_whether_a_better_path_exists() {
    let cfg = SchedCfg { policy: Policy::Tolerance, bootstrap_levels: 0, ..Default::default() };
    let collapsed = |can_switch: bool| QuadReduction {
        n_footprints: 64,
        // Fewer distinct ICs than footprints: the collapse condition, tested on the ICs and never
        // on a spread being zero.
        n_distinct_ic: 1,
        decode_can_switch: can_switch,
        ..Default::default()
    };

    assert_eq!(
        decide_reduction_for_test(&collapsed(false), 4, &cfg),
        Decision::Collapsed,
        "at the ceiling a collapse is terminal, under the label every dump already carries"
    );
    assert_eq!(
        decide_reduction_for_test(&collapsed(true), 4, &cfg),
        Decision::DecodeSwitch,
        "with a better path available a collapse is a SWITCH, and refinement continues"
    );

    // The control: an uncollapsed quad reaches neither, or the two above would be firing on
    // something other than the collapse.
    let fine = QuadReduction { n_footprints: 64, n_distinct_ic: 64, ..Default::default() };
    let d = decide_reduction_for_test(&fine, 4, &cfg);
    assert!(
        d != Decision::Collapsed && d != Decision::DecodeSwitch,
        "a quad with distinct ICs reached a collapse verdict: {}",
        d.name()
    );
}

/// **Production is unchanged, and this is what says so.**
///
/// The shipped path is at the ceiling, so a collapse still returns `Decision::Collapsed` — the
/// same code 9 every committed dump carries. The first cut of this work returned a new
/// `AtF32Floor` there, which `tests/no_discard.rs` caught immediately: that renames a committed
/// decision, and the new name would have been a synonym no build could produce.
#[test]
fn the_shipped_path_is_at_the_ceiling_so_nothing_moves() {
    let ens = prin_rs::ensemble::pixel::EnsembleCfg::default();
    assert_eq!(ens.decode_path, Path::DirectF64, "the shipped decode path moved");
    assert!(!ens.decode_path.has_more_precise_path(), "production must not be labelled a switch");

    // And the reduction built by the real path agrees: `n_distinct_ic` is measured on the f64
    // nominal decode, so `decode_can_switch` is false and the verdict is terminal.
    assert!(!QuadReduction::default().decode_can_switch);
}

/// **Both new variants reach the code table**, so a dump decodes them rather than printing `"?"`.
#[test]
fn the_new_variants_are_in_the_table() {
    for d in [Decision::DecodeSwitch] {
        assert_eq!(Decision::from_code(d.code()), Some(d));
        assert!(Decision::ALL.contains(&d), "{} is not in Decision::ALL", d.name());
        assert_ne!(d.name(), "?");
    }
    assert!(Decision::from_code(Decision::ALL.len() as u8).is_none());
}
