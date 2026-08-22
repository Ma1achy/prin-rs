# Build TODO

Tracks the gated build order. Each step must pass before the next begins.
Plan: `BRIEF.md` (spec), `CLAUDE.md` (agreement), `NOTES.md` (findings).

Status: `[ ]` not started · `[~]` in progress · `[x]` done · `[!]` blocked

---

## Step 1 — reproduce the reference

- [x] Run `reference/README.md` smoke test
  - Got `median |dE/E| = 3.892633125701676e-09`, expected ~3.9e-09
- [x] **Finding:** this number is *not* a Rust target — it is a median over 36 rows
      including jittered ones, and the seeding conflict makes it unreachable

## Step 2 — physics core, f64 · *no PR, lands with Step 4*

- [x] Crate skeleton: `Cargo.toml`, `src/lib.rs`
- [x] `src/real.rs` — `Real` supertrait + per-precision floors (`TINY`, `DRIFT_FLOOR`,
      `SYNC_EPS`, `DIST_FLOOR`). f64 gets the reference's literals
- [x] `src/vec2.rs` — `Vec2<T>` with operator impls
- [x] `src/physics/newton.rs` — `accel` (takes **eps squared**), `pair_dists`
      (`PAIRS = [(0,1),(0,2),(1,2)]` ordering is load-bearing)
- [x] `src/physics/energy.rs` — `energy`, `inertia`, hyperradius `R` at t=0
- [x] `src/physics/shape.rs` — `shape_vec` (Hopf map)
- [x] `src/physics/burrau.rs` — masses, `R0`, `V0`
- [x] **Harness now, not later:** `tools/xcheck/{cases,dump_ref,compare}.py`, `algebra` case
      — validated on functions with no chaos, so at Step 4 the harness is not a suspect
- [x] `tests/burrau_constants.rs` — **gate:** `M=12`, `R=2.2361`, `E=-12.8167`, crossing `0.9652`
- [x] `tests/shape_vec.rs` — unit norm; invariant under translation, rotation, scale
- [x] **gate:** `algebra` cross-check ≤1e-15

## Step 3 — adaptive-dt integrator, unregularised · *no PR*

- [x] `src/integrate/leapfrog.rs` — `dt = eta·min_pairs(r^{3/2}/√(G(mi+mj)))`, KDK,
      explicit `is_finite` in the loop condition
- [x] `tests/leapfrog_energy.rs` — near-field t=13
- [x] **gate:** drift **falls with `eta`** on tame pixels. Failure on close encounters is
      expected — record which pixels and how many, don't tune around it
- [x] If drift is insensitive to `eta`, the equations are wrong. **Do not begin Step 4**

## Step 4 — Aarseth–Zare · **FIRST PR**

- [x] `src/integrate/az/state.rs` — `[u1,p1,u2,p2,t]`, `to/from_array9`
- [x] `src/integrate/az/system.rs` — `to_reg`, `phys_from_state`, `energy_phys`, `to_cartesian`
- [x] `src/integrate/az/hamiltonian.rs` — **`gamma()` and `deriv()` adjacent**
      (Γ exists only in the reference's docstring; write it)
- [x] `src/integrate/az/rk4.rs` — fixed `dtau`, `h = dtau*(!done)` **transcribed as-is**
- [x] `src/integrate/az/driver.rs` — sync loop, `dmin`, `switches`
- [x] `src/integrate/az/reference_body.rs` — `choose_reference`, `RefPolicy`, `ref_disagree`
- [x] `src/compat.rs` — `ReferenceQuirks` named booleans

Transcription hazards (all silent):
- [x] `g1` carries `ma*mc`, `g2` carries `ma*mb` — cross mass pairs
- [x] `R3` sub-term **sign flips** between `g1` and `g2`
- [x] `_LmatT_apply(p1, L2p2)` — `p1` in the **matrix** slot, deliberate
- [x] `energy_phys` keeps the `P1·P2/ma` cross term
- [x] `dtau` fixed per sub-interval, **never** shrunk at close approach
- [x] `deriv` leaves `A`,`B` unfloored; `phys_from_state` floors them. Asymmetric — keep
- [x] Overshoot: clip the time accumulator only, leave state overshot

Validation chain — **order matters, each link anchors the next**:
- [x] 1. Burrau constants anchor Cartesian `energy()` (Step 2)
- [x] 2. `energy_phys(to_reg(r,v)) == energy(r,v)`; `to_cartesian(to_reg(·))` round-trips
- [x] 3. **`gamma(s,E) == A·B·(energy_phys(s) − E)`** — a sign error in Γ cannot survive this
- [x] 4. FD: central-difference Γ vs analytic `deriv`, all 8 components
  - [x] assert error falls as **h²** (two step sizes)
  - [x] `#[should_panic]` sign-flipped variant, to prove the test has teeth
- [x] 5. Runtime monitor: `Γ ≡ 0` along the trajectory. Track `max|Γ|`, dump it

Gates (paste raw output into the PR):
- [x] (a) the chain above
- [x] (b) two-body radial collision: `d_min < 1e-10`, `|dE/E| < 1e-12`
- [x] (c) cross-check vs Python at f64, **nominal copies only**
  - [x] pin `n_sync=32, eta=0.01, max_steps=30000` explicitly
  - [x] row order `index = jy·nx + jx`, x fastest — unit-test against a hardcoded 3×2
  - [x] dump from `tb_az.integrate_az` directly, **not** `tb_all_az`
  - [x] **log the chosen reference body per sync as a compared column** — a near-tie broken
        differently by `argmax` fails the cross-check looking exactly like a transcription error
  - [x] hard-assert ≤1e-13 at `t≤2` — **got bitwise identity at t≤1, 8.9e-16 at t=2**
- [!] ≤1e-10 at `t=13` — **got 1.9e-10 absolute.** Exceeds the target, for reasons
      unrelated to correctness; see the horizon table. Proposed as a BRIEF amendment
  - [x] **divergence-vs-horizon table** `t ∈ {0.5,1,2,4,8,13}` — the curve shape is stronger
        evidence than a single pass/fail
- [!] Record the healthy-f64 `error_ratio` distribution — **deferred to 5a**: error_ratio is
      a cross-copy statistic and there is no ensemble until 5a. The drift distribution that
      feeds it is recorded (near-field t=13 nominal: max |dE/E| = 1.7e-08)
- [x] Benchmark `deriv`, extrapolate to 1024²×8 at t=13

## Step 5a — ensemble, fields, outputs · **SECOND PR**

- [x] `src/grid.rs` — slice decode, per-axis cell widths
- [x] `src/ensemble/jitter.rs` — per-pixel seeding `(i,j,seed)` via SplitMix64 → `Pcg64Mcg`
  - [x] **jitter uses per-axis cell width** — the reference uses `hx` for *both* axes
        (not x-only as first described; latent on square grids, wrong on any other)
  - [x] copy 0 always un-jittered — load-bearing, assert it
  - [x] `tests/seeding_golden.rs` pins first N values so a dep bump can't move ICs
- [x] `src/ensemble/stats.rs` — median, MAD, spreads
- [x] `error_ratio`: MAD inside, **max** aggregation over pixels
- [!] **MAD defeats the field's purpose.** Damaged-vs-healthy separation is 1.06 with MAD
      and 59.51 with max-deviation. Both computed and dumped; contradicts a CLAUDE.md
      non-negotiable, so the choice is the user's
- [x] `spread_shape`: dump **both** BRIEF §4's definition and `refine_test.svar` — different
      statistics; the brief's is spec, `svar` is the one with a reference
- [x] `d_min_ref`, `d_min_true`, `d_min_gap` as three fields
- [x] per-copy `finite` flag; never update `d_min` from non-finite state; never discard a copy
- [x] `src/output/raw.rs`, `src/output/png.rs`, `src/render.rs`, `src/config.rs`, `src/bin/prin.rs`

Validation (no oracle — invariants only):
- [x] `error_ratio` == 1.0 exactly at t=0
- [x] shift/scale invariance of the ratio
- [x] exactly-integrable control (two-body + distant third) sits at 1
- [x] **step-size convergence: `error_ratio − 1` falls as `eta → 0`** — the check that makes a
      field with no oracle trustworthy
- [x] **gate:** gauge invariance, `α ∈ {0.25,1,4}`, `t` by `α^{3/2}`, ~10 decimals
- [x] **gate:** `max(error_ratio)` under the Step-4-derived bound; report max/median/p99/argmax

Confounds to dump, not hide:
- [x] `sigma_E_0` and `sigma_E_t` as separate fields — `σ_E(0) ∝ cell width`, so the ratio
      inflates with resolution for a trivial reason. Threatens §8 exp 3, and exp 1 too
      (2×2 aggregation compares across different effective cell widths)
- [x] **confirmed with a number:** cell width falls 9x, `sigma_E(0)` falls 8.6x
      (proportional), `ensemble_spread` falls only 2.1x — not proportional, so free of it
- [x] `d_min_gap` **identically zero in all five regions** — NOTES §2.1 settled
- [x] `ref_disagree` vs `error_ratio` correlation, unshared run — settles NOTES §1
- [!] Predicted `error_ratio ≈ 1 + 1e-3`; **measured max 2.449, p99 1.433, median 1.000000**.
      The prediction was wrong because it assumed no damaged pixels; 23 of 1024 have
      `drift_max > 1e-3`. See the MAD finding below

## Step 5b — termination and outcome encoding · **THIRD PR**

- [x] `src/outcome.rs` — `classify_legacy()` (reference-exact, xcheck only) and `classify()`
      (BRIEF §2.4 3-bit state + 2-bit detail)
- [x] escape arm is a **port** (`tb_all_az` sync-boundary test); collision/triple arms are
      **invented** — marked as such in the source, in one block
- [x] ≥2-pair rule encoded as `detail = 3`, both arms
- [x] `r_coll` canonical, fixed at t=0, from each trajectory's own `R`, never co-moving
- [x] collision detection samples **inside** the AZ inner loop (sync boundaries are 0.4 apart)
- [x] property test: two pairs below `r_coll` is **never** labelled an ordinary binary collision
- [x] property test: the `(R1,R2,R3)` → pair-index mapping, for every reference body
- [x] scale invariance of outcomes and `t_end` under `α` — 25/25 labels held, `t_end` bitwise
      for `α ∈ {0.25,4}` and `9.9e-16` relative for `{3.7, 1/3}`
- [x] the reference-matching entry points keep `stop_on_event = false`, so the cross-check is
      unchanged — verified, `az_t13` still passes
- [!] `deep interior` terminates in bounded wall-clock — but as a **binary** collision, not a
      triple. BRIEF §2.6 is wrong about this pixel; see NOTES §2.4
- [!] the **escape arm never fires at `t = 13`** (0 of 1024 near-field pixels; 109 at `t=20`),
      so at the project horizon the outcome fields are driven entirely by the invented arm.
      NOTES §2.5
- [x] `r_coll` sensitivity sweep `{1e-4, 1e-3, 1e-2}·R` on 64×64
- [!] the default is on a **slope, not a plateau**: collision fraction 0.0000 → 0.0242 → 1.0000
      across the three decades. Reported, not tuned
- [x] outcome fractions vs `lc_stable`: **0 label flips of 4096** at every `r_coll` — the
      branch cut does not reach the outcome encoding, unlike `spread_shape`

## Spec corrections carried into Step 6's PR

- [x] `spread_event` := disagreement over the **event class** (currently-tightest pair at each
      sync boundary, joined with the terminal class), not the terminal `(state, detail)`.
      Measured: 110 vs 0 nonzero pixels at `t_max = 8`; 165 vs 22 at `t = 13`, strictly nested
- [!] the "~4 time units earlier" framing does **not** reproduce as a lead time — zero on the
      22 pixels both flag. The gain is coverage (7.5x) and horizon-independence. NOTES §2.8
- [!] the playhead value can **un-fire** (the tightest pair fluctuates), so it is non-monotone
      in the horizon. `spread_event_max` dumped alongside; the spec one stays the default
- [x] `t_spread_event` dumped, **NaN** when the copies never disagree
- [x] BRIEF §2.6 rewritten: `deep interior` is a binary collision (landed in PR #4)
- [x] `d_min_true` primary, `r_coll` a recorded parameter, in BRIEF §4 and CLAUDE.md

## Step 6 — f32 · **FIFTH PR**

- [x] `impl Real for f32` — floors asserted, not assumed: `tests/f32_precision.rs` checks that
      `1e-300` casts to exactly 0 and that `1e-15` is a no-op against `ulp(13) = 9.537e-7`
- [!] `TINY*TINY` underflows at f32 and `A*B` is a product of two floored quantities, so a
      doubly-degenerate state gives `dtau = inf`. Caught by the explicit `is_finite` test, not
      by the floor — the guard doing the work is not the one it looks like
- [x] `Precision` wired through the single dispatch in `render.rs`
- [x] ICs generated once in f64 and cast down (`ensemble/jitter.rs`) — never per precision
- [x] acceptance thresholds parameterised by precision. Gate (b): f64 `1.2881e-11` / `2.9473e-14`,
      f32 `1.2218e-11` / `2.8553e-6`. The `d_min` half meets the f64 bound as written; the
      `1e-12` energy bound is five orders below f32 eps and is not asserted there
- [x] all four {f32,f64} × shared-reference {on,off} reported — NOTES §4
- [x] **Q1** the conditioned branch fixes f32 `spread_shape`: reference branch inflates the
      median 32x (6.1582e-2 against the f64 truth 1.9095e-3); conditioned tracks to 1.2%
- [!] this run produced **no** NaN pixels at f32 on either branch, where PR #3 saw them — the
      NaN observation is configuration-dependent and must not be quoted as general
- [x] **Q2** shared reference does not help: `spread_shape` median moves 1%, and at f32 the
      worst pixel gets **56% worse**. Default stays unshared. Demoted, not eliminated
- [!] the aggregate hides it — sharing changes 268/1024 pixels' `spread_shape`, worst **1.86x**,
      while the median moves 1%. Exactly the NOTES §1 failure mode
- [x] **Q3** the branch cut **does** reach the outcome encoding at f32: **152 of 1024** label
      flips at the default `r_coll`, against 0 at f64. The zero at `r_coll = 1e-2` is
      saturation, not reassurance
- [x] bug found by the shared-reference path: the nominal can now terminate early, so its
      `refs` is shorter than `n_sync` and the policy indexed off the end. Fixed, with a
      regression test
- [!] f32 tail: 2 pixels of 1024 have `|dE/E| > 1`. Not data. Flagged, never discarded

---

## Spec amendments accepted, not yet written into BRIEF.md

- [x] all of the above are written into `BRIEF.md`, `CLAUDE.md` and `NOTES.md` as they land

## Step 7 — BRIEF §8, the reason for the build · **SIXTH PR**

- [x] **A.** refinement criterion by 2x2 aggregation, no quadtree. `alpha` distribution,
      region discrimination, `alpha_sigma_E(0)` as the control (truth exactly 1.0)
- [!] the control does **not** return 1.0 per quad: interdecile width **0.48** at `E+1 = 8`,
      falling as `1/sqrt(E)`. The criterion resolves **regions**, not individual quads
- [!] **sample-size bias**: a parent pools `4(E+1)` copies against a child's `E+1`, so the
      exponent reads +7.6% high before any physics. Subsampling to matched counts gives
      1.0191 at `E+1 = 8` and 0.9999 at 16. Use an rms deviation here, not `max_dev`
- [x] the criterion discriminates: `spread_shape` alpha median 0.17 (near-field), 0.34
      (body2 core), 1.18 (mid-field), 1.17 (far) — chaotic regions do not pay for refinement
- [x] `sigma_E_0` and `sigma_E_t` reported separately, per PR #3
- [x] **B.** statistical convergence by subsampling a fixed 128x128 grid
- [!] at `n = 64` every fraction has interdecile scatter **2.0-4.4x the quantity itself**.
      That is the 1.2x-to-18.8x mechanism, quantified
- [!] p99 of anything heavy-tailed is **not estimable** at any `n` tested (3.3-6.8 for
      `spread_shape p99`, flat)
- [!] `drift max` reads 0.000 scatter at `n <= 256` and 1.000 at `n = 4096` — apparent
      stability is the statistic never drawing the tail
- [!] 128x128 found 7 pixels of 16384 with `|dE/E| > 1` that 64x64 could not. `error_ratio`
      flags 7 of 7. Drift falls 13 orders for a 3.3x change in `eta`: resolution, not a wrong
      equation. **`eta = 1e-2` is not sufficient at this resolution** — written into BRIEF §2.5
- [x] no scheduler machinery built, per the instruction

## Consolidation · **SEVENTH PR**

- [x] `RESULTS.md` — a findings document readable cold, not a PR summary
- [x] `results/` committed: 256x256 images per region, 64x64 raw dumps, captured stdout for
      every example
- [x] three new standing rules in `NOTES.md` §3a and `CLAUDE.md`
- [x] BRIEF §2.5's remedy implemented: flag on `error_ratio`, re-integrate at `eta/4`, up to
      3 passes, recording coarse and refined values and the `eta` used
- [!] one pass is **not** always enough — `deep interior` needs all three (`1.10e12` -> `1.99e1`
      after one, `1.146e-1` after three at 128x128). My PR #6 claim that one step suffices was
      near-field only
- [!] the pass budget is calibrated on **f64**: f32 needs ~6 at 128x128 and leaves 1578 of 65536
      still flagged at 256x256 with the default 3. Read `n_still_flagged`
- [!] `error_ratio` detects **spread, not drift**: after refinement `deep interior` has 0 pixels
      flagged and a worst drift of `1.15e-1`. An ensemble that drifts *together* is invisible to
      it. Threshold on `energy_drift_max` if absolute conservation is what matters
- [x] the spread image is log scaled over the grid's own p1-p99; a linear ramp hid the structure
      entirely
- [x] `RESULTS.md` §5 lists what this invalidates in the prior numpy work
- [x] comparison tests and experiment examples pin `refine_flagged: false`, so they measure
      arithmetic rather than which pixels got a second pass
