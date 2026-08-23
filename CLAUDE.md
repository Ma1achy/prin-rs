# Working agreement — prin-rs

`BRIEF.md` is the authoritative spec. Read it in full before writing code or restructuring
anything. This file is the working agreement: what holds, what must not be broken, and how work
gets reviewed. It exists so the agreement survives context resets.

---

## SCOPE

Uniform grid, one pass, PNG + raw dump, exit.

**No** quadtree. **No** scheduler. **No** GUI. **No** interaction. **No** streaming.

These omissions are deliberate, not deferred. **If the program grows a scheduler, that is a bug,
not progress.** The refinement criterion of BRIEF.md §8 is testable by aggregating 2×2 blocks of a
fine uniform grid — a scheduler buys nothing it does not already have.

---

## PORT, DON'T DERIVE

`reference/tb_az.py` is a validated Aarseth–Zare implementation. **Port it. Transcribe the
algebra; do not re-derive it.**

The algebra is error-prone and **fails silently**. Two sign errors in it were invisible until
someone finite-differenced the Hamiltonian and compared against the analytic derivatives. There is
no symptom to notice — wrong AZ algebra produces trajectories that look like physics.

Finite-difference the Hamiltonian against your analytic derivatives as a test. Not once during
development — as a committed test.

---

## THE DIAGNOSTIC SIGNATURE

**Energy drift that does NOT fall when you shrink the step size means a WRONG EQUATION, never a
step-size problem.**

Internalise this. It has caught three separate bugs in this project. Measured example: a co-moving
softening length gave `|dE/E| = 3.06e-02`, *identical* at `dt=1e-4` and `dt=2e-5`. Insensitivity to
step size is the tell.

When drift is bad, the first question is never "smaller step?" — it is "which equation is wrong?"

**The same signature, one level up: a diagnostic that gets WORSE as resolution improves is
measuring the wrong thing.** This has now caught two separate bugs.

- A co-moving softening length gave drift identical at `dt=1e-4` and `dt=2e-5`.
- The `Gamma` residual, first normalised by `A*B`, read `3.9e-1` on a trajectory whose actual
  drift was `1.1e-13`, and got *worse* as `eta` fell — because `A*B -> 0` at a collision, so
  it blew up exactly where the integrator works hardest.

Whenever a number moves the wrong way under refinement, suspect the definition before the
implementation.

**And do not read a scaling law off a single trajectory.** `d_min` near a collision is set by
where a sample happens to land relative to the crossing. Its *scale* falls as `eta^2`, but
the realisation scatter across phases is as wide as the shift between decades, so one
trajectory can show a ratio of 1.0 across a 10x change in `eta` and mean nothing by it.
Measure scalings over an ensemble.

---

## BUILD ORDER

1. **CPU first. f64 first.** Rayon over pixels.
2. **Match the Python reference to ~1e-10 on a small grid.** Nothing else is verifiable until this
   holds. Do not proceed past it.
3. **Only then** go generic over `f32`/`f64` from one source. Both precisions exercised from that
   point on — the f32 question (BRIEF.md §8 experiment 2) is a reason this exists.

GPU comes later. Correctness now.

---

## NON-NEGOTIABLES

**`error_ratio` uses the maximum deviation internally, aggregates by max.**
Internally: `max|x − median|`, with non-finite treated as an infinite deviation. Not a standard
deviation, which returns NaN on precisely the pathological pixel the statistic exists to flag —
and **not MAD**, which was the earlier answer to that and overshot. Robustness is the wrong
property here: with 8 copies, one wild value sits above the median of eight deviations and is
arithmetically invisible. Measured damaged/healthy separation **1.06 with MAD, 59.51 with max
deviation**; a pixel whose worst copy drifted 120× the total energy read 1.1369 under MAD, inside
the healthy p99 of 1.0756. Keep `error_ratio_mad` dumped; never gate on it.

Aggregate over pixels by **max**, not median — max tracks damage at Spearman +0.956 against +0.599
for median. **That is a separate decision from the one above**, at a different level: it compares
`error_ratio_max` to `error_ratio_median` across footprints. A per-pixel correlation between the
two within-footprint estimators is a different measurement and is not evidence about it. Both
decisions land on the word "max" for unrelated reasons; do not collapse them.

Treat the result as a **boolean flag**; its magnitude is unstable.

**`spread_event` is over the EVENT CLASS, never the terminal outcome.**
The event class is which pair is *currently the tightest binary*, evaluated at every sync
boundary and joined with the terminal `(state, detail)` for copies that have terminated. The
terminal outcome was explicitly rejected as the contributor and reinstating it is a regression,
not a simplification: it is terminal-grain and inverts under lockstep — early in the march
nothing has terminated, every copy agrees, and the field reports maximum confidence at exactly
the playhead where least is known. Measured, near-field 32×32, nonzero pixels of 1024: at
`t_max = 8`, **110 against 0**; at `t_max = 13`, 165 against 22 and strictly nested. The gain is
coverage and horizon-independence, not lead time — on pixels both flag, the lead is zero.
The `(state, detail)` encoding stays as the **outcome**, for classification and rendering.

**`d_min` is primary; `r_coll` is a recorded parameter, not a physical constant.**
The collision label is *derived* from `d_min`. Measured: the collision fraction runs
0.0000 → 0.0242 → 1.0000 across `r_coll/R ∈ {1e-4, 1e-3, 1e-2}` while the grid's `d_min/R`
spans less than one decade. No threshold in that range is a physical event boundary. `r_coll`
appears in every output header. The default `1e-3` separates tail from bulk on this slice and
claims nothing more.

**Ensemble offsets are a FIXED Halton (2,3) prefix indexed by copy index.**
Not per-pixel, not pseudo-random. Fixed, so copy `k` sits at the same offset in every footprint
at every refinement level and a parent shares its children's perturbation pattern — common
random numbers by construction. Measured: the control's per-quad noise floor falls from 0.4796
to **0.0010** at `E+1 = 8` and the parent/child correlation from 0.175 to **0.9998**.
`Scheme::Pcg` reproduces the reference's stream and every result measured before the switch;
never make it the default again.

**But it does not reduce the scatter in `alpha_shape`, and more copies will not either.**
Sampling noise is only ~7% of that scatter (`var` falls 5.725e-1 -> 5.331e-1, against
`var(alpha_E) = 3.75e-2`). The other 93% is chaotic divergence. The per-quad scatter is **0.63**
under either scheme against a region separation of ~1.0 — so the criterion resolves regions, not
quads, for a reason compute cannot buy off.

**Never compute a refinement exponent by pooling a 2x2 block. Render at two resolutions.**
"A fine grid contains every coarser scale by aggregation" is true of the positions and false of
the ensemble: with fixed offsets a pooled block is four *exact repeats* of one pattern at four
cell centres, while a true parent carries offsets scaled to its own, wider cell. Measured on a
control whose true value is exactly 1.0, the pooled exponent is **+38.6% at `E+1 = 8`** (falling
as `1/E`) against a true two-resolution error **flat at +2.3%**; pooling also understates the
per-quad scatter by ~2x. Not a correction factor — a different measurement.

**In tame regions the criterion resolves individual quads; in chaotic ones nothing does.**
True two-resolution `alpha_shape`: `mid-field` and `far` sit at 1.023 with an interdecile of
**0.0004-0.001**, `near-field` and `body2 core` at 0.04-0.18 with **1.1-1.3**. The region
separation is 0.986. "Not resolvable per quad" *is* the answer for a chaotic quad, and the
scatter is the measurement rather than an error bar around one.

**Read the interdecile, never the variance, for `alpha_shape`.**
Excess kurtosis is **110**: the variance lives in the tails, the interdecile describes the bulk
(interdecile/sd = 0.866 against a normal 2.563). A scheduler decides per typical quad. The
Halton switch cut sampling variance 267,000x on the control and moved the interdecile not at all;
quoting the 6.9% variance reduction as the improvement would be quoting the tail.

**Never discard an ensemble copy.**
Every pixel carries exactly `E+1` copies, always. A badly-integrated trajectory is a *measurement
outcome* — "this could not be determined" — not missing data. Discarding biases the sample toward
tame trajectories, which on a chaos instrument is backwards.

**`r_coll` and `epsilon` are canonical and fixed at `t=0`.**
Expressed as a fraction of the initial hyperradius `R`, evaluated once at `t=0`, **never
co-moving**. A co-moving length makes the Hamiltonian time-dependent and destroys energy
conservation. An absolute length breaks the scale invariance the project deliberately quotients
out — measured, the same physical system gave answers differing by 1.66× purely from an arbitrary
overall size.

**Triple collision/ejection uses the ≥2-pair rule, encoded as `detail = 3`.**
Two or more pairs below `r_coll`, not all three. By the triangle inequality `|AB| < r_coll` and
`|AC| < r_coll` forces `|BC| < 2 r_coll`, so "exactly two pairs" is reachable and is already a
near-triple. Requiring all three would silently misclassify it as an ordinary binary collision.
`detail = 3` means "all three" for both arms — triple collision and triple ejection.

**`deep interior` is a binary collision, not a triple.**
BRIEF §2.6's warning predates AZ. Measured in both implementations: pair (0,2) closes to
`2.28e-5`, pairs (0,1) and (1,2) never register even at `r_coll = R`, and it reaches `t = 13`
in about a second with `|dE/E| ~ 1.4e-7`. The 190 s failure was the *unregularised* integrator.
Do not restore a "this region is intractable" assumption — that has already produced one false
finding in this project.

**The escape arm contributes nothing at `t = 13`.**
Burrau's escape happens later than the horizon: zero of 1024 near-field pixels fire at
`t_max = 13`, 109 at `t_max = 20`. So `spread_event` and the outcome image at the project
horizon are driven entirely by the collision arm — the one arm *with* a reference is the one
that is silent. Escape is also sampled at sync boundaries and latches on first firing, both
transcribed from the reference.

**f32 is usable, with two named caveats — and one of them is the branch cut.**
Median drift 9.3e-6 against f64's 2.8e-9, which is what `eps ~ 1.19e-7` over ~5000 RK4 steps
predicts; outcome labels agree with f64 on 1022 of 1024 pixels. Caveat one: 2 pixels of 1024
have `|dE/E| > 1` and are not data — flagged by `error_ratio`, never discarded. Caveat two: the
**unconditioned** LC branch at f32 inflates `spread_shape` 32x and flips **152 of 1024 outcome
labels**, where at f64 it flipped none. Never run f32 on the reference branch.

**The floors are the one place f32 and f64 are not the same algorithm.**
`1e-300` casts to exactly zero (f32 uses `1e-37`); `ulp(13)` at f32 is `9.537e-7`, so a `1e-15`
sync slack is no slack at all (f32 uses `1e-6`). Also: `TINY*TINY` underflows at f32 and `A*B`
is a product of two floored quantities, so a doubly-degenerate state gives `dtau = inf` — caught
by the explicit `is_finite` test, **not** by the floor. Know which guard is doing the work.

**`eta = 1e-2` is not sufficient above ~64x64, and the failure is a cliff.**
At 128x128 near-field, 7 pixels of 16384 have `|dE/E| > 1` (worst `1.49e4`), all finite, all at
`d_min ~ 2e-3` against `r_coll = 2.21e-3`. Drift falls thirteen orders for a 3.3x change in
`eta`, so it is resolution, not a wrong equation. `error_ratio` flags 7 of 7. Re-integrate
flagged pixels at finer `eta` rather than lowering `eta` globally — no scheduler needed.

**The refinement exponent is a region-level statistic, not a per-quad one.**
`alpha = log2(spread_parent/spread_child)` on a control whose true value is exactly 1.0 has an
interdecile width of **0.48 at `E+1 = 8`**, falling as `1/sqrt(E)`. The measured region
separation is about 1.0. So it resolves regions and not individual quads. Also: a parent pools
`4(E+1)` copies against a child's `E+1`, and a spread estimator's expectation depends on sample
size — **match the counts** or the exponent is biased (+7.6% at `E+1 = 8`) before any physics
enters.

**A statistic can report maximum confidence precisely when it is least informed.**
The inversion, not mere noise. `drift max` scatter reads 0.000 at `n <= 256` because small
samples never draw the one bad pixel of 16384 — stable at the wrong answer. Terminal-outcome
purity reads pure under lockstep because nothing has terminated yet. Ask what the statistic
would say about a system nothing is known about; if the answer is "confident", it is wrong.

**`n_sync` fixed while `t_max` varies compares different discretisations.**
`dtau = eta*dt_left/(A0*B0)`, so changing `t_max` at fixed `n_sync` changes the step size and
the rows are not one trajectory at different playheads. Scale `n_sync` with `t_max`, or run
once and evaluate at each boundary.

**Never compute a refinement exponent by pooling — and in the scheduler, never aggregate silently.**
A quad holds `N x N` footprint spreads and needs one number. Measured: **half the shared decisions
flip** between mean, median and p90 (54.1% and 49.1% against median in near-field), and the trees
overlap by 3-13%. median under-refines thin structure (blind to a filament crossing a quad); mean
over-refines and blows the budget; p90 refines deepest, narrowest, and floors 55% of leaves. Three
schedulers wearing one name. State the aggregation wherever a tree is quoted.

**Coarse `N` OVER-refines. The cheaper quad is a false economy.**
The concern on record was that a low `N` makes a quad call itself *coherent* by undersampling its
area. Measured, the opposite: leaf count falls monotonically with `N` (near-field 106, 31, 19, 16
at `N = 4, 7, 8, 16`), because a noisy spread estimate biases toward *refine* — the conservative
failure direction. `N = 4` spends 4x the quads of `N = 16` on the same region.

**In the scheduler, `alpha_hi` does more work than the criterion, and `tau` is often inert.**
`tau = 1e-8` and `1e-6` give identical trees in near-field (the spread never falls that low), while
`alpha_hi` from 0.20 to 0.50 collapses the tree **80x**. The alpha median is +0.389, so the
threshold sits inside the distribution. Sweep both before quoting any tree; never pick either to
make a picture look right.

**Never conclude "no effect" from an aggregate without the per-pixel distribution.**
An aggregate can only say the distribution did not move; it cannot say the pixels did not.
Measured twice in one PR: LC-branch `spread_shape` rows identical to five digits while **all
1024** pixels moved, worst 6.7%; shared references moving the median 1% while 268 of 1024
pixels moved, worst **1.86x**. Both would have been written up as "inert".

**A test that cannot fail is indistinguishable from a test that passes.**
Ask what would have to be true for the test to fire. Three catches: a label-flip count of zero
at `r_coll = 1e-2`, where every pixel collides anyway so the label is *saturated*; a
scale-invariance test at `t_max = 6` where nothing terminated, so `t_end` was the horizon and
the invariance was the rescaling's own arithmetic; and an FD test on `Gamma` that a sign error
shared by `Gamma` and `deriv` would have passed. If the answer is "nothing in this
configuration", the test is decoration.

**Test `is_finite` explicitly.**
`NaN >= x` is `false`, so a diverged trajectory never satisfies a loop exit and burns its entire
step budget. Measured: 354 s against 3 s nominal.

**Do not shrink the fictitious step at close approach.**
AZ's time transformation `dt = |R1||R2| dtau` *already* shrinks the physical step — that is what
regularisation buys. Shrinking `dtau` too drove `dt → 1e-13` and produced a false "this region is
intractable".

**Do not build an `L_z` version of `error_ratio`.**
Released from rest, `L_z = 0` for every copy, so `sigma_Lz(0) = 0` and the ratio is `0/0`.
Structurally undefined for this whole configuration family.

---

## SMOKE TEST

`reference/README.md` carries one. Expected **median `|dE/E| ≈ 3.9e-09`** (measured here:
`3.892633125701676e-09`). **Reproduce it before porting**, so you know the reference runs on your
machine and you are comparing against a live number rather than a quoted one.

---

## WORKING WITH REVIEW

**Open a PR per meaningful chunk of work.**

In every PR description, paste the **actual acceptance test output** (BRIEF.md §5) — not a
summary. The numbers are the review.

| test | expected |
|---|---|
| two-body radial collision | `d_min < 1e-10` with `\|dE/E\| < 1e-12` |
| gauge invariance: ICs by `alpha ∈ {0.25, 1, 4}`, `t` by `alpha^{3/2}` | `shape_vec` spread identical to ~10 decimals |
| `error_ratio` at `t=13`, near-field | `1.0000` |
| Burrau constants | `M=12`, `R=2.2361`, `E=-12.8167` |
| cross-check vs Python reference at f64 | agreement ~`1e-10` on a small grid |

Review is against **the physics and the reference implementation**: whether the AZ algebra
transcribes correctly, whether the outcome encoding is right, whether an invariant has been quietly
broken. It is **not** for compilation or Rust idiom — the compiler is better at that.

**Report negative and messy results.** A clean summary that hides scatter is worse than useless
here. Several conclusions in this project have been overturned by whoever was closest to the code,
and that has been the most valuable part of the process. **If a framing in the brief is wrong, say
so.**

---

## OPEN DESIGN QUESTION — shared reference body

BRIEF.md §3. AZ picks a reference body per trajectory and it can change mid-run (in the reference,
at each of the `n_sync` sync points — `reference/tb_az.py:165`). Should all `E+1` copies of a pixel
share the nominal copy's reference?

**Implement it as a flag. Measure both.** Notes toward settling it are in `NOTES.md`.

The flag governs **cross-copy sharing only** — not freezing the reference across time. Freezing
across time would break AZ outright.
