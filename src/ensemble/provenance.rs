//! One source of truth for the ensemble config, and a record of every departure from it.
//!
//! # Why this exists
//!
//! `prin` and every harness used to construct [`EnsembleCfg`] **independently**: production took
//! `::default()`, harnesses wrote their own struct literals. There was no single source of truth,
//! so "the production config" existed nowhere — it was whatever `default()` happened to return,
//! and 111 literal sites could disagree with it silently.
//!
//! They did. `refine_flagged: false` was introduced at `c03fc85` **correctly**, in experiment and
//! precision harnesses where measuring the repaired kernel would hide the thing being measured,
//! and the same commit wrote down the invariant: *the `render-*.txt` runs have it ON*. Over six
//! days the line was copied into render harnesses one file at a time, no commit message arguing
//! for it, until every committed panel carried it. Second instance of the pattern after
//! `k_frac = 1.0` shipping as the scheduler default.
//!
//! **The failure was never that someone chose `false`. It is that nothing recorded the choice.**
//!
//! # The shape of the fix
//!
//! - [`EnsembleCfg::production`] is the one literal. `Default` delegates to it.
//! - [`Override`] is a **named** value per field, so a harness declares what it changes rather
//!   than writing an anonymous struct field: `production().with_overrides(&[Override::RefineFlagged(false)])`.
//! - [`EnsembleCfg::overrides_vs_production`] **derives** the list by diffing, so a config gets a
//!   correct declaration however it was built — including the 111 existing literals, and
//!   including any future one whose author forgets. A hand-maintained list can itself go stale,
//!   which is the same class of failure this module exists to stop; the derived diff cannot.
//! - [`EnsembleCfg::provenance`] renders it for a header. **A convention cannot fail silently if
//!   the value is in the log.**
//!
//! # The two compile-time guards
//!
//! Both the diff and [`Override::apply`] **destructure or match exhaustively, with no `..` and no
//! `_` arm**. Adding a field to `EnsembleCfg` breaks the build here until it is handled, so the
//! record cannot silently fall behind the struct it describes. That is the joint that carries the
//! whole idea: a mechanism which reports "no overrides" because it does not know about a field is
//! exactly the failure being fixed, one level up.

use super::jitter::Scheme;
use super::pixel::EnsembleCfg;
use crate::decode::Path;
use crate::integrate::az::{driver::DtauMode, driver::StepLimit, reference_body::RefPolicy};
use crate::outcome::EscapeRule;
use crate::physics::ftle::FtleOpts;

/// A named change to one field of [`EnsembleCfg::production`].
///
/// Total over the struct: one variant per field, so `with_overrides` has no gaps and
/// [`Override::apply`]'s match is exhaustive by construction.
#[derive(Clone, Copy, Debug)]
pub enum Override {
    NExtra(usize),
    JitterFrac(f64),
    JitterScheme(Scheme),
    Seed(u64),
    TMax(f64),
    NSync(usize),
    EscapeRule(EscapeRule<f64>),
    ClosureK(usize),
    StopOnEscape(bool),
    EscapeEvery(usize),
    EscapeConfirm(bool),
    DtauMode(DtauMode),
    ClampFinalStep(bool),
    StepLimit(StepLimit),
    StepLimitF(f64),
    Eta(f64),
    MaxSteps(usize),
    RefPolicy(RefPolicy),
    LcStable(bool),
    RCollFrac(f64),
    StopOnEvent(bool),
    RefineFlagged(bool),
    RefineThreshold(f64),
    RefineEtaFactor(f64),
    RefineMaxPasses(u8),
    DecodePath(Path),
    KeepCopyOutcomes(bool),
    KeepCopyShapes(bool),
    KeepBoundaryShapes(bool),
    KeepDriftHist(bool),
    Ftle(Option<FtleOpts>),
    FtleDt(f64),
}

impl Override {
    /// Apply to a config in place. **Exhaustive by design** — see the module docs.
    pub fn apply(self, c: &mut EnsembleCfg) {
        match self {
            Override::NExtra(v) => c.n_extra = v,
            Override::JitterFrac(v) => c.jitter_frac = v,
            Override::JitterScheme(v) => c.jitter_scheme = v,
            Override::Seed(v) => c.seed = v,
            Override::TMax(v) => c.t_max = v,
            Override::NSync(v) => c.n_sync = v,
            Override::EscapeRule(v) => c.escape_rule = v,
            Override::ClosureK(v) => c.closure_k = v,
            Override::StopOnEscape(v) => c.stop_on_escape = v,
            Override::EscapeEvery(v) => c.escape_every = v,
            Override::EscapeConfirm(v) => c.escape_confirm = v,
            Override::DtauMode(v) => c.dtau_mode = v,
            Override::ClampFinalStep(v) => c.clamp_final_step = v,
            Override::StepLimit(v) => c.step_limit = v,
            Override::StepLimitF(v) => c.step_limit_f = v,
            Override::Eta(v) => c.eta = v,
            Override::MaxSteps(v) => c.max_steps = v,
            Override::RefPolicy(v) => c.ref_policy = v,
            Override::LcStable(v) => c.lc_stable = v,
            Override::RCollFrac(v) => c.r_coll_frac = v,
            Override::StopOnEvent(v) => c.stop_on_event = v,
            Override::RefineFlagged(v) => c.refine_flagged = v,
            Override::RefineThreshold(v) => c.refine_threshold = v,
            Override::RefineEtaFactor(v) => c.refine_eta_factor = v,
            Override::RefineMaxPasses(v) => c.refine_max_passes = v,
            Override::DecodePath(v) => c.decode_path = v,
            Override::KeepCopyOutcomes(v) => c.keep_copy_outcomes = v,
            Override::KeepCopyShapes(v) => c.keep_copy_shapes = v,
            Override::KeepBoundaryShapes(v) => c.keep_boundary_shapes = v,
            Override::KeepDriftHist(v) => c.keep_drift_hist = v,
            Override::Ftle(v) => c.ftle = v,
            Override::FtleDt(v) => c.ftle_dt = v,
        }
    }
}

impl EnsembleCfg {
    /// Apply named overrides to a config, returning it. The idiom every harness should use.
    pub fn with_overrides(mut self, ov: &[Override]) -> Self {
        for o in ov {
            o.apply(&mut self);
        }
        self
    }

    /// Every field of `self` that differs from [`EnsembleCfg::production`], as
    /// `(field, this value, production's value)`.
    ///
    /// **Derived, not declared.** Comparison is on each field's `Debug` rendering, which is
    /// uniform across `f64`, the enums and `Option<FtleOpts>` alike, and which treats two `NaN`s
    /// as equal — the right answer for a settings diff, where `!=` is not.
    ///
    /// The destructuring is exhaustive with no `..`, so adding a field to `EnsembleCfg` fails to
    /// compile here until it is listed. See the module docs.
    pub fn overrides_vs_production(&self) -> Vec<(&'static str, String, String)> {
        let p = Self::production();
        // Exhaustive on purpose. Do not add `..`.
        let EnsembleCfg {
            n_extra, jitter_frac, jitter_scheme, seed, t_max, n_sync, escape_rule, closure_k,
            stop_on_escape, escape_every, escape_confirm, dtau_mode, clamp_final_step,
            step_limit, step_limit_f, eta,
            max_steps, ref_policy, lc_stable, r_coll_frac, stop_on_event, refine_flagged,
            refine_threshold, refine_eta_factor, refine_max_passes, decode_path,
            keep_copy_outcomes, keep_copy_shapes, keep_boundary_shapes, keep_drift_hist, ftle,
            ftle_dt,
        } = self;

        let mut out = Vec::new();
        macro_rules! cmp {
            ($name:literal, $mine:expr, $theirs:expr) => {
                let (a, b) = (format!("{:?}", $mine), format!("{:?}", $theirs));
                if a != b {
                    out.push(($name, a, b));
                }
            };
        }
        cmp!("n_extra", n_extra, p.n_extra);
        cmp!("jitter_frac", jitter_frac, p.jitter_frac);
        cmp!("jitter_scheme", jitter_scheme, p.jitter_scheme);
        cmp!("seed", seed, p.seed);
        cmp!("t_max", t_max, p.t_max);
        cmp!("n_sync", n_sync, p.n_sync);
        cmp!("escape_rule", escape_rule, p.escape_rule);
        cmp!("closure_k", closure_k, p.closure_k);
        cmp!("stop_on_escape", stop_on_escape, p.stop_on_escape);
        cmp!("escape_every", escape_every, p.escape_every);
        cmp!("escape_confirm", escape_confirm, p.escape_confirm);
        cmp!("dtau_mode", dtau_mode, p.dtau_mode);
        cmp!("clamp_final_step", clamp_final_step, p.clamp_final_step);
        cmp!("step_limit", step_limit, p.step_limit);
        cmp!("step_limit_f", step_limit_f, p.step_limit_f);
        cmp!("eta", eta, p.eta);
        cmp!("max_steps", max_steps, p.max_steps);
        cmp!("ref_policy", ref_policy, p.ref_policy);
        cmp!("lc_stable", lc_stable, p.lc_stable);
        cmp!("r_coll_frac", r_coll_frac, p.r_coll_frac);
        cmp!("stop_on_event", stop_on_event, p.stop_on_event);
        cmp!("refine_flagged", refine_flagged, p.refine_flagged);
        cmp!("refine_threshold", refine_threshold, p.refine_threshold);
        cmp!("refine_eta_factor", refine_eta_factor, p.refine_eta_factor);
        cmp!("refine_max_passes", refine_max_passes, p.refine_max_passes);
        cmp!("decode_path", decode_path, p.decode_path);
        cmp!("keep_copy_outcomes", keep_copy_outcomes, p.keep_copy_outcomes);
        cmp!("keep_copy_shapes", keep_copy_shapes, p.keep_copy_shapes);
        cmp!("keep_boundary_shapes", keep_boundary_shapes, p.keep_boundary_shapes);
        cmp!("keep_drift_hist", keep_drift_hist, p.keep_drift_hist);
        cmp!("ftle", ftle, p.ftle);
        cmp!("ftle_dt", ftle_dt, p.ftle_dt);
        out
    }

    /// One line naming every departure from production, for an output header.
    ///
    /// Reads `production` when there are none — an explicit statement rather than an empty
    /// string, because a blank field and an absent field look the same in a log and the whole
    /// point is that the choice is recorded either way.
    pub fn provenance(&self) -> String {
        let ov = self.overrides_vs_production();
        if ov.is_empty() {
            return "production".into();
        }
        let body: Vec<String> =
            ov.iter().map(|(k, a, b)| format!("{k}={a} (production {b})")).collect();
        format!("production + {} override(s): {}", ov.len(), body.join(", "))
    }

    /// Whether this config departs from production at all.
    pub fn is_production(&self) -> bool {
        self.overrides_vs_production().is_empty()
    }
}
