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

**The screen floor is the everyday stop, and without it the criterion is measured wrong.**
`tile_size = quad_width * zoom / N`; at `N = 8` on a 512x512 viewport samples stop being
displayable at **level 6**. PR #11 descended to **12**. Under the veto, near-field goes 4617 quads
to **549** and **61.2% of leaves are stopped by the view, not the criterion**. It is a **veto,
never a trigger** — `Camera::veto` returns `Option<Decision>` and cannot return `Split` — and it is
**view-relative and never cached on a quad**: a floored quad must refine again when zoomed into.
`deep interior` is byte-identical with it on and off, so the veto neither causes nor fixes that
region's bad tree.

**Under the screen floor, `tau` is the dominant knob and `alpha_hi` is demoted — the opposite of
PR #11.** `alpha` is a *rate* statistic and needs levels to express itself; with `bootstrap = 2`
and a floor at 6 there are **four** discretionary levels against twelve. `tau` is a *level*
statistic and keeps all its room. Measured: `alpha_hi` 0.20 -> 0.50 collapses near-field **21.7x**
(was 80x) and `far` **not at all**; `tau` moves `far` **64x** and near-field **16x**, where it was
called inert. Sweep both, and say which regime the tree was measured in.

**A difference can be small because both sides are right or because both are dead.**
Before reading any agreement number, assert each side still resolves what it is supposed to. Three
catches in one planning session: a curvature term on an affine chart (identically zero at every
depth), a linearised f32 sum whose samples all collapse to `x0` and agree perfectly with a direct
path that collapsed too, and an `E` null that a veto-capped tree would have produced whatever `E`
did. `decode::distinct` is the guard: count distinct ICs first, read divergence second.

**A collapsed decode makes the criterion maximally confident.** Identical footprints give
`ensemble_spread` exactly zero, which reads as "perfectly resolved" and stops the descent with a
small tidy tree built from nothing. Treat a collapsed quad as **undetermined**, the same way a
non-finite copy is a measurement outcome rather than missing data.

**The deep-zoom floor is a property of where you zoom.** PR #11's level 45.87 is conditional on the
chart coordinate being O(1), and the condition was never stated. The same box at the chart origin
has **no cell-width floor at all** in the tested range, on either precision, because there is no
O(1) neighbour for the increment to be absorbed into. Quote the coordinate magnitude with any floor
depth.

**The linearised decoder buys ~24 levels over f32 and none over f64.** `L-split` (x0 in f64 on the
CPU, `delta` and `J_D.delta` in f32, summed in f64) tracks `direct_f64` rung for rung: both hold
64/64 samples to depth 44 and reach 1 by 50. The literal formula `x0 + J_D.delta` **all in f32**
collapses on exactly the same curve as plain f32. The ICs must be formed as absolute O(1) numbers
before integration — three-body separations are O(1) — so the linearisation escapes the
chart-coordinate floor and not the IC-magnitude one.

**`N` and `E` fail in opposite directions and are not interchangeable.** `N` controls how well a
quad knows its **area**; undersampling inflates the between-footprint variation that drives
`alpha`, so coarse `N` **over**-refines. `E` controls how well a footprint knows its **value**;
undersampling deflates the within-footprint spread compared against `tau`, so low `E`
**under**-refines. Measured without the veto, near-field leaf count runs 742 -> 2713 -> 3463 at
`E+1 = 2, 4, 8`. Never trade one against the other as if they were the same knob.

**The jitter is in CHART coordinates, not one body's Cartesian position.** The original form added
the offset to `c.r[slice.body]`, which is right only because `Chart::BodyPlane` writes `(u, v)`
straight into that slot. On an oblique or shape chart it perturbs a body instead of taking a
sub-cell sample of the chart. Bitwise identical for `BodyPlane`; `tests/seeding_golden.rs` holds
that.

**A ranking is invariant to a monotone rescaling of the signal; a threshold is not.**
The between-footprint arm runs **1.17x** the within arm in `near-field` and **9.56x** in `far`,
so swapping criterion at a fixed `tau` changes the effective threshold by up to 8x, region by
region. Compare criteria as **orderings**, and a rescaling costs nothing; compare them against
`tau` and the measurement scores the rescaling instead of the signal.

**Two different faults give the same flat error curve: a BAD ordering and NO ordering.**
Count the signal's distinct values before reading any curve. `within/median` is flat to `B=767`
in near-field with **5418 distinct values of 5461** (modal 0.3%) — a fine-grained ordering that
is actively bad, beaten by random at every budget past 383. `frac_hot_within` and `layout` are
flat with **58 and 78** distinct (modal 40.8%) — no ordering at all, and their curve is the
tie-break's scan order. Different faults, different fixes, and `error(B)` alone cannot tell them
apart.

**But signal resolution is not what makes a ranking good.** `frac_hot_between` is the **best**
criterion in `deep interior` on **65** distinct values, beating a 4994-valued one. And
`term_grad` is **NaN on 97.1%** of near-field yet reaches the oracle's zero by `B = 383`: the
2.9% it scores are the right quads. A high `nan%` is a property to read, not a defect to hide.

**Draw the tree, not only the image — and never over a uniform base.**
The adaptive render says *what is displayed*; the wireframe says *where the tree cut*. A coarse
texel tells you a leaf is coarse; only the wire tells you whether the structure around it was
subdivided *around* it or straight *through* it. PR #11 drew boundaries over a **uniform** base,
conflating the two, and `deep interior`'s bad tree survived a whole build unnoticed. Measured:
the wire pair at `B = 682` shows `within/median` shredding near-field's **top-left corner** into
a fine mesh while the collision region sits in **two level-1 leaves** — so it is not noisy and
not failing to order, it is systematically refining the wrong corner. No table showed that.

**Greedy is a strong reference, never a ceiling — and this has now happened.**
Greedy on immediate `Δerror` is optimal only when gains are independent and immediately
available. On a tree they are neither — a quad whose own split gains little may unlock children
with large gains two levels down, and greedy declines it. **A criterion beating `greedy_oracle`
indicates lookahead value, not a bug**, and there must be no test asserting it dominates.
Measured at `t = 20`, near-field: `greedy_oracle` **plateaus at 0.00048** from `B = 383` through
`B = 3071` while `first_div` reaches **0.00000 at `B = 1535`**. A dominance test would have
fired on correct behaviour.

**What is displayed decides what the criterion should measure, and it is measurably blind to
half of it.** The production colouring is bivariate: hue from the shape sphere (aligned with
`spread_shape` by construction), lightness from a scalar. Measured on near-field at `B = 341`,
the gap between the best criterion and random runs **total** under `outcome`, **6.6x** under
lightness=spread, **1.5x** under diffusion and **1.2x** under FTLE. Under a lightness field the
criterion does not read, it is barely better than spending the budget at random — and the best
criterion changes identity too. Choose a criterion under the colouring that will ship.

**`error = 0` against a finite reference means "matches this sampling", not "correct".**
The reference is the fully-refined tree at one sample per pixel. At the screen floor sub-pixel
structure is sampled arbitrarily — which side of a filament a pixel lands on is an accident of
where its sample fell. The exactly-locatable zero is a virtue for *comparing* criteria; it is
not a statement about image quality, and every table quoting the curve says so.

**The within/between distinction is one of SCALE and AGGREGATION, not of kind.**
Matched for extent and sample count, the two arms agree to **1.01** in every region: they are
the same estimator. The brief's premise — "the ICs there are identical up to perturbation" —
does not describe this implementation: `jitter_frac` is 0.5 and `halton_offset` returns
`[-1,1)^2` scaled by cell width, so the copies span the **whole cell, edge to edge**, and the
Halton control's true `alpha` is exactly **1.0**, which an irreducible within-point statistic
cannot be. But `rho` on quads containing a transition is only **0.58-0.64**, so at their actual
settings they rank quads differently and the practical conclusion survives its mechanism.

**A control with no randomness in it cannot measure sampling noise.**
`sigma_E(0)` looks like the perfect control for `alpha_sibling_spread` — true `alpha` exactly
1.0, true sibling range exactly 0, no integration. It reads **0.003 and does not move with `N`
or `E+1` at all**, which is the tell: under the fixed Halton prefix the offsets and footprints
are both fixed, so the whole quantity is deterministic and the residual is geometry. Keep it as
a floor, label it, and vary a `Pcg` seed for a real draw. Measured that way, sampling noise p90
is **0.21-0.36** against `sib_tau = 0.5`, and the sibling median is 0.45 (near-field) and
0.79-1.05 (`deep interior`) — the threshold sits inside the noise-broadened bulk in both.

**`t_end` termination is not escape, and conflating them contradicts a standing result while
appearing to agree with it.** `t_end` is set by whichever terminating event came first.
`deep interior` reads **terminated = 0.99 with the escape arm silent** — those are collisions.
Carry `terminated_fraction` and `escape_fraction` separately.

**The camera has no position term, so panning changes no scheduling decision.**
`Camera::veto` reads `tile_size_px`, which depends on the quad's width and the camera's
`half_world` and `viewport` — and **not** on `cx`/`cy`. A pan study that reported "the tree
persists perfectly" would be reporting an identity. `Camera::covers` exists to *measure* what
would be evictable and is deliberately **not** consulted by `veto`: adding a position term to the
floor would make a quad's decision depend on where the camera points, which is what
"never cached as a quad fact" exists to keep out.

**The temporal accumulators' event arm already existed.** `spread_event_max` is a running max
over boundaries, `t_spread_event` a first-divergence time that is NaN rather than `t_max`, and
`spread_event_latched` the persistence-guarded latch. Only the **continuous** arm was missing,
and it is not a null at `t = 13`: `running_max` reaches 0.00158 in `deep interior` where
`within/median` sits at 0.01509.

**A map can be continuous and still lose an axis — and "seam" was the wrong diagnosis.**
The shipped hue map `chroma*(cos h, sin h)` with `h = atan2(n2,n1)` and
`chroma = C_MAX*hypot(n1,n2)` is identically `C_MAX*(n1,n2)` — agreement `4.2e-17` over a sphere
sweep — so it is linear and has no branch cut at all. Its fault is that it is exactly **2-to-1**:
it discards `n0`, so a tight binary with a distant third body and a wide pair with a close third
render bitwise the same. Before writing a continuity test, check whether the map is even
injective; a discontinuity is the failure that is easy to name and not the one that was there.

**An auto-ranged ramp cannot tell "no signal" from "signal", and a ratio threshold is not
enough.** The lightness window is each region's own p1–p99, so a field with no dynamic range has
its **noise** stretched to full scale. `far` reads `error(root) = 0.60` under the shipping
colouring against `0.00000` under `outcome`, which looks like a rescue and is not: its window is
`(1.3e-9, 1.1e-8)`. A ratio test missed it — span `x8`, above any sensible ratio bound. The
second arm compares the window against the region's **own median energy drift**: a field whose
whole range sits within two orders of the integrator's arithmetic is not physics.

**The p99 of a composite field can be set by a different estimator than the one being ramped.**
`ensemble_spread = max(spread_shape, spread_event)`. The event arm has **5 distinct values**
(modal 98.2%) and dominates only **1.7%** of near-field footprints — all in the top tail. So it
sets the p99 and nothing else: the window ran to `2.857e-1` (exactly `2/7`) where the continuous
arm's own p99 is `2.244e-2`, **12.7x narrower**. A linear ramp over a window an order of
magnitude too wide, set by a staircase describing 1.7% of the region. Colour on `spread_shape`,
and print `quantisation` and `event_arm_fraction` before any image.

**`ensemble_spread` carries a scale term, so a multi-resolution render is partly a picture of the
tree.** It is a spread over copies jittered within the **cell**. Measured per level: cell width
halves each level (a proportional field would show `2.000`) and the median spread ratio runs
**1.19–1.62, falling with depth to 1.048** — sub-linear and saturating, the chaotic-divergence
signature. About 12% of the lightness range across five levels. `t_end` is the scale-free control
and is flat at 0.998–1.009. Stated rather than corrected: normalising by cell width would change
what the field means.

**`sum p_i = 0` does NOT catch the crossed-mass swap in the decoder, whatever the reference
says.** Both forms give `p_lam*(1 - (m0+m1)/M01) = 0`; measured `7.9e-17` crossed and `5.6e-17`
uncrossed. What catches it: the **Jacobi round-trip** `p_rho == (m0 p1 - m1 p0)/M01` (`1.1e-16`
against `6.8e-2`) and the **kinetic-energy identity** (`4.4e-16` against `2.6e-1`). Both are
empty at `m0 == m1`, where the two forms are the same expression — Burrau has `m0 = 3, m1 = 4`,
but the mass simplex passes through that line.

**`(Lz,E)` and `(Lz,K)` are one chart with two labels.** The reference lists them separately as
its most-machinery item, but its own warp parameterises both by `K(t) = K_max t^gamma` and then
reports `E = U + K(t)` — a relabelling of the axis, not a different sweep. Bitwise identical over
the unit square. Only `gamma_k` makes them differ.

**A construction that assumes a COM-centred input returns a drifting system without one.**
`momenta_for`'s rigid-rotation step is `v = omega J r`, whose total momentum is
`omega J (M R_com)`. Every decoded configuration is COM-centred, so it would never have fired in
production — until a chart handed it something else. Centre internally; state which frame `Lz` is
about.

**The criterion question is settled by the AGGREGATION, not by the within/between arm.** Measured
at `B = 383` at 1024² under the shipping colouring: `frac_hot_between/median` is the **only**
criterion that beats the random band, in **both** measurable regions — near-field `0.12486`
against a band of `0.13059–0.14726`, `deep interior` `0.08075` against `0.08334–0.10058`. And
`between/median` is beaten by random too (`0.15243`, `0.10205`), so swapping the arm is not the
fix. Counting how many footprints are hot beats taking their median.

**Greedy is beaten by an arbitrary scan order where the field is featureless.** In `far` at
`B = 12287`, `greedy_oracle` reads `0.26391` while every other criterion reads `0.10885` — 2.4x
worse. Immediate `delta-error` is noise there and greedy chases it. Every non-greedy row is
identical to five digits, including ones with a single distinct value, which is what says the
ranking is irrelevant rather than good.

**Softness in an image is a raster size, not a rendering fault — measure the dimensions first.**
Wireframe lines are written with integer pixel `set` calls and adaptive texels are
nearest-neighbour, so neither can be soft in the file. If an image looks blurry, a viewer is
upscaling a small raster. Corollary, and it has now cost a round trip twice: **never let a
validation run write into `results/`.** A `criterion_metric -- 3 8` pass overwrote committed 512²
artefacts with 128x64 ones, and a small raster reads as a rendering fault rather than a stale
file.

**An argument hardcoded past is worse than an argument missing.** `pan_sequence` took a
`viewport` argument that set the camera while `frame_res` stayed hardcoded at 384, and
`between_vs_within` took one while its render was a literal `512` — so asking for a larger raster
produced the same small one and looked like a rendering limit. Both are fixed; the remaining
`Camera::framing(..., 512)` sites use the camera for **scheduling** only and write no images,
which is why they stay.

**`prin --size` drives the image AND the per-pixel dump, and they want opposite sizes.** The
images want 1024²; the raw dumps are documented as 64×64 and at 1024² they are **320 MB per
region**. Run it twice — once large for the images, once small for the dump — and never let the
large run's `.raw` land in `results/`. One knob for two artefacts with opposite requirements is
the shape of the problem; the workaround is in `results/README.md`.

**The GLSL is the pin for the latent chart; the LaTeX reference guessed and was wrong.**
`Ma1achy/principia-ii`, `src/shaders/principia/frag.glsl:19-59`. `MU_MAX = 5.0` (the chart
reference said 4.0), `Q_MAX = 2.0` (said 1.0), and the mass saturation is
`MU_MAX*(2*sigmoid(z) - 1)` — algebraically `MU_MAX*tanh(z/2)`, **half** the gain of the spec's
`mu_max*tanh(z)`. The reconstruction algebra was already right, crossing included. Its `decodeIC`
carries `z0..z9` but never reads `z2`/`z3`, and puts beta at index 0 where the spec names the chart
`(z_alpha, z_beta)` — so the preset images are **transposed relative to the GLSL**, which is
faithfulness to the spec and not a bug. `decoder.rs`'s module header holds the index table.

**The `z = 0` Lagrange landmark is blind to every constant this project has got wrong there.**
It is a named physical configuration — masses `(1/3,1/3,1/3)`, separations all `sqrt 3`,
`I = 1` — checkable by eye in the render, and a strong gate on the reconstruction algebra. But at
`z = 0` the momentum coordinates and mass logits are all zero, so `MU_MAX`, `Q_MAX` and the choice
between `tanh(z)` and `2*sigmoid(z)-1` every one of them drops out. `I = 1` and `COM = 0` are
worse: they are algebraic identities (`I = cos^2 a + sin^2 a`; `m0r0 + m1r1 = -M01 m2 lam` cancels
`m2 r2`) that hold under **any** mass factors. Keep all three as wiring guards, label them as such,
and put the constants under a test that fires when they move.

**A chart's tameness is set by WHICH COORDINATES it varies, not by where it is centred.**
The standing result — chart families sit at `alpha` 0.99-1.01 and do not exercise the criterion —
held for twelve of thirteen instances and was read as a property of the base point. It is not.
`preset_shape` and `preset_prho` share `z0 = 0` **exactly** and differ **5.7x** in leaf count (577
against 3268), with `alpha` interdecile 3.07 against 0.12 and 1.237% undetermined pixels against
0.007%. Only the first sweeps `(alpha, beta)`, the configuration coordinates, and so passes through
collision-adjacent shapes; the momentum slices hold one configuration and vary its initial velocity.
A mixed basis (`preset_shape_pl`) lands with the momentum slices, so the configuration sweep does
not dominate one. Moving a base point to find chaos is tuning; changing which coordinates the plane
spans is not.

**A default that spans two coordinate systems silently means two different things.**
`half = 0.05` is a body position in Burrau units on `BodyPlane` and a sigmoid pre-image on
`Latent`. One shared default shipped every GLSL preset at `half = 1.0` against the reference UI's
`Slice +/- 3.0e+0` — a **3x crop**, `alpha in [0.446,1.125]` against `[0.120,1.451]`, 46% of the
azimuth against 90%. The structure was right and the window was wrong, which reads as
*similar but not the same* and sends you looking in the physics. `Chart::default_half()` is the
one table now. And the crop control is the `_h1` twin: same chart, same basis, one number changed.

**Comparing across colour modes is how a rendering choice gets mistaken for a physics bug.**
The GLSL image first compared against was a **continuous** field; the reference's own WebGPU panel
reads `Colour mode: Event class, Palette: viridis` — **discrete** — and that one looks far closer.
A continuous field and a categorical map cannot look alike even when both are correct. Render
prin-rs under a mode matched to the reference (`Colouring::EventClass`, `<case>_event.png`), and
state the mode of every committed image. Its alphabet is **fixed at 27 slots, not data-derived**,
so the same class is the same colour in two slices — which makes adjacent ordinals close in colour
by construction, so the legend and the per-class histogram are the instrument, not the image.

**A basis that pairs coordinates pairs them BY SLOT, and renumbering must carry the partner.**
The GLSL's `shape_pl` is `q1 = e0 + e6`, `q2 = e1 + e7` — in its indexing, `beta` with
`pLambda.x` and `alpha` with `pLambda.y`. This port renumbers alpha and beta into the spec's order
and the first cut did **not** carry their momentum partners across. **No transposition of
`q1`/`q2` repairs it**: that gives `e_beta + e_pLy`, `e_alpha + e_pLx`, still crossed. It is a
genuinely different 2-plane through the 8D space, not a reorientation, which is why it rendered as
*twisted* rather than tilted — the coupling sets how momentum co-varies with configuration across
the slice. Measured `max |dIC|` against the correct plane: **5.4483** crossed, **1.2042**
transposed. An index assertion alone would have passed on the transposition, so the test carries
both as negative controls. `shape_pl` is the **only** preset with a cross-coupling, which is why it
was the only one that looked wrong in this particular way — that consistency is itself evidence.

**A tree can be set by a CAMERA VETO while the table calls it criterion-bound.**
`chart_gallery`'s `bound` column read `crit` unless the *budget* was exhausted and never asked what
actually stopped the descent. Measured from the `.prnq` dumps, which carried `decision` all along:
**`Decision::MaxRelDepth` stops 95%+ of leaves on 23 of 26 charts**, and 100% on three, where every
leaf sits at one depth — complete capped trees whose leaf counts are facts about the cap. This is
the mechanism behind the standing "the chart families do not exercise the criterion": the criterion
decides **under 1% of leaves** there, and the `alpha` near 1.0 describes quads a cap forced. Same
lesson as the screen floor, at a second stop condition, unnoticed for the same reason — nothing
printed which one fired. Never quote a leaf count without the stop-reason breakdown.

**`preset_shape` at the corrected window is where the criterion fails outright — 16 leaves.**
0% veto, 8 `Floor` and 8 `Keep`, depth 2, against a complete 4096 — the only case in the set whose
tree is entirely its own decisions. And not on a featureless field: ramp span **18966.6**, four
decades, the widest in the set, `alpha` interdecile **13.47** against 0.04-0.99 for the tame rows.
The wide window admits the smooth surroundings, their spread falls below `tau`, `Agg::Median` reads
the quad resolved, and the fractal core sits unrefined inside a level-2 quad. That is
"median under-refines thin structure" at full strength, and `half = 1.0` hid it behind a
plausible-looking 577-leaf tree.

**"The escape arm contributes nothing at `t = 13`" is about Burrau's near-field and does not
generalise.** On the latent charts `escape_fraction` is **0.9894-1.0000** — essentially everything
escapes. `preset_shape` is the counter-example inside that family: escape **0.0547** with
`terminated_fraction` median 0.984, so its terminations are collisions (event histogram dominated
by `collision d0=361886`). Carrying `terminated_fraction` and `escape_fraction` separately is what
makes the difference legible.

**A finding read off a wireframe is a finding about an appearance. Test the cause — then check the
test can fire.** "Refinement goes to smooth regions" was read off a wireframe at the wrong window.
The mechanism (terminated regions are absorbing, so copies agree, `spread_event` collapses and the
quad reads resolved) predicts leaf depth **anti-correlated** with `terminated_fraction`. Measured on
26 charts, it is **readable on 2**: three distinct failure modes that a Spearman cannot tell apart —
**x constant** (all leaves at one depth), **y constant** (`terminated_fraction` takes one value, ten
charts), **y saturated** (modal share >90%, twelve charts). The two readable charts **disagree**:
`shape_sphere` -0.2245 with medians 0.766 -> 0.000, `preset_shape_h1` +0.3756 with medians
0.812 -> 1.000. Neither established nor refuted. And **read the per-depth medians, not the pooled
correlation** — `shape_sphere` puts 908 of 970 leaves at the two deepest levels, so one number over
that design understates a four-step fall.

**`preset_prho` and `preset_plambda` are a control, not just a curiosity.**
Every pixel of both is the *same triangle* at a different initial velocity, so `spread_shape` at
`t = 0` is identically zero across the whole slice. Any structure in them is purely
momentum-driven — which makes the pair what separates configuration effects from momentum effects.

**A distinctness check keyed on the wrong quantity reads as a collapsed decode.**
Positions in the latent decode do not depend on the momentum coordinates at all, so the `prho` and
`plambda` presets are constant-**configuration** slices: every pixel is the same triangle released
with a different initial velocity. Keying `decode::distinct` on a body position reads **1 of 81**
there. Nothing has collapsed — the chart simply does not move that quantity. Two different faults
give the same count, and the guard has to measure what the chart actually varies.

**A fixed threshold fails on BOTH sides, and that is the argument for rank.** `tau_display = 1e-4`
sits at the **0.4th percentile** of the `charts/` leaf-spread distribution (75,359 leaves, median
`6.61e-4`) and the **4.3rd** of the whole corpus (92,880 leaves) — **state which scope**, because
the first write-up quoted the first under a heading naming the second. So 95.7-99.6% of quads
clear it and the tree is uniform at **max depth**. The other side is
in the same corpus: where the bulk sits *below* `tau` everything keeps and the tree is uniform at
**depth 2** — 16 leaves against a complete 4096. **16 of the 18 trees the camera veto does not bind
are stopped that way** (`far`, `deep interior`, every deep zoom step; medians `9.45e-5` down to
`4.26e-8`). Selectivity requires the threshold to *cut through the bulk*, and dynamic range decides
whether any fixed value can — `spearman(p99/p1, depth variance) = +0.727`. **A ranking cannot land
above or below a distribution.**

**The camera veto is doing the stopping, not the criterion.** `ScreenFloor` or `MaxRelDepth` stops
≥95% of leaves on **21 of 69 dumps**. So "13 of 17 charts at 99%+ max depth" is not the criterion
saying *split* and being obeyed — it is the criterion never saying *stop*. **The observed
uniformity is what a permissive criterion looks like when a veto terminates the descent.** Never
quote a leaf count without its stop-reason breakdown.

**`preset_shape` is ALPHA-bound, not `tau`-bound, and the shape of its tree does not tell you
which.** 16 leaves, depth variance 0, widest dynamic range among the charts — it reads exactly like
the upper-side `tau` failure. Its leaf-spread median is `2.86e-1`, **3400x above `tau`**: it clears
the spread gate on every leaf and is stopped by `alpha` (8 `floor` + 8 `keep`, zero spread-gate
failures). It is the only tree in the corpus exercising the alpha gate. Read the decision column;
a mechanism read off a tree shape is a guess.

**A quantile hot rule makes `n_hot` a non-signal, so BOTH masks are kept.** Under any quantile rule
the count above the cut is set by the rule, not the field (31 of 64 at `N = 8, q = 0.5`), so
`frac_hot` carries essentially nothing. `frac_hot_between/median` is the **best criterion measured
on this project**, so replacing the absolute mask would have deleted the best-performing signal and
read as an improvement. `spatial::HotRule` selects which mask the *shape* criteria read;
`frac_above_tau_*` is untouched. The exception: on a **tied** field the tie structure sets the
count — a two-valued field reads the same at `q = 0.5, 0.75, 0.9`, which is the case when the
five-valued event arm dominates.

**And the obvious desaturation test cannot fail.** *"Assert `n_hot < N^2` for a stated majority"*
passes trivially and unconditionally under a quantile rule. The form with teeth runs both masks
over **one descent** and asserts they disagree — measured 100% saturated / 100% single-blob under
the absolute rule against 0% / 59.7% under `q = 0.5`, with 40.3% resolving more than one component.
The absolute arm is the control.

**The mask saturation is REGIONAL, and in `far` the mask is EMPTY rather than full.** `far`'s median
is `4.26e-8` against `tau = 1e-4`, so nothing clears the cut: `n_hot == 0`, `perimeter_ratio` is
`NaN`, and every mask-derived criterion takes **one** distinct value over all 16 leaves. `deep
interior` already resolves a median of 2 components unmasked. It is near-field and the latent
charts where it is full. A full mask and an empty mask are the same threshold landing on either
side of the distribution.

**Desaturating the mask COARSENS the ordering.** Near-field's median component count runs 1 -> 5
from absolute to `q[0.50]`, while `Criterion::LayoutRel`'s distinct-value count falls
**78 -> 26 -> 17 -> 9** across `abs, q[0.50], q[0.75], q[0.90]` against `Criterion::Layout`'s
steady 58. With `n_hot` pinned, `largest/n_hot` has only as many values as there are component
sizes. Not disqualifying — signal resolution is not what makes a ranking good — but a criterion
whose ordering coarsens as its input improves is worth watching, and `error(B)` decides.

**The criterion is rarely what decides anything, and there are FOUR stop reasons.** Measured across
the corpus and the re-run sweep: the **camera veto** stops >=95% of leaves on 21 of 69 dumps; the
**budget** stops up to 91% at low `tau` with no veto (near-field 869 of 1498, `deep interior` 1357
of 1498); the **spread gate** stops everything at high `tau` (16 leaves, all `keep`); and the
**alpha gate** stops `preset_shape`, the only tree in the corpus that exercises it. "The criterion
refines too much" is the weaker statement. Print the stop-reason breakdown with every leaf count.

**`alpha_hi` collapses near-field 79x without the veto and `tau` is inert above it.** 1498 -> 19
leaves between `alpha_hi` 0.20 and 0.50; at `alpha_hi >= 0.5`, `tau` changes nothing at any rung.
`tau` is live in **one row of thirty-six**. Under the veto the ratio falls to 21.68x, which is the
screen floor demoting `alpha_hi` and promoting `tau`, as recorded. And within the live row `1e-8`,
`1e-6`, `1e-4` and `3e-4` give a **bitwise identical** tree -- the whole live range of `tau` is
`3e-4 .. 3e-3`, one decade, the 11.6th to 96.9th percentile.

**WHICH RUNG OF A SWEEP IS DEGENERATE IS A FACT ABOUT THE REGION.** `1e-8` and `1e-6` were dropped
from the `tau` ladder as "measuring the always-split regime twice" -- true in near-field, where
`1e-8`, `1e-6` and `1e-4` are bitwise identical, and **false in `far`**, whose median is `4.26e-8`
so `1e-8` is the only rung below its bulk. Without it `far` read 16 leaves in all 32 cells and the
sweep said "`tau` is inert here", silently contradicting the standing 64x result. The regional
spread medians span **six orders** -- `4.26e-8`, `9.45e-5`, `9.75e-4` -- which `sched_sweep`'s own
module header has said since it was written. Restoring the rung recovered the 64x exactly.

**A span quoted between two NAMED rungs is an argument hardcoded past.** `sweep_screen` printed the
`tau` span as `at(1e-8)/at(1e-6)`, which read `x0.00` once `1e-8` was removed and would have been
reporting whichever pair happened to straddle the bulk in one region even when present. Take
max/min over the whole ladder at the `alpha_hi` where the knob is live. Same defect as
`pan_sequence`'s hardcoded viewport, at a different site.

**THE TREADMILL DOES NOT HAPPEN; THE OPPOSITE DOES.** The standing argument for rank was that
"spread grows with `t` everywhere", so any fixed threshold must eventually fire on every quad.
Measured on a fixed tree across `t in {4..20}`: near-field's median leaf spread **peaks at
`t = 10` and then falls 81x**, ending below `tau = 1e-4`; `deep interior` falls **31x
monotonically** and is under `tau` from `t = 13`. The mechanism is already on record one level
down -- **terminal states are absorbing**, so as termination saturates the copies share an outcome
and the disagreement collapses. At large `t` a fixed `tau` fires **nowhere**, the spread gate keeps
everything and the tree **shrinks**: near-field 256 -> 40 leaves, `keep` 0 -> 20. This
*strengthens* the case for rank -- a rise-then-collapse has no correct fixed value at either end
and no monotone schedule that tracks it, where a monotone rise at least would.

**`order_queue` never read `cfg.criterion`.** It sorted on `red.spread(agg)`, so every
`--order spread` run in the corpus ordered by the within arm whatever its header said. It
reproduces the corpus exactly only because every committed run has `criterion=within` and
`signal(Within, agg)` *is* `spread(agg)`. Fixed, asserted, and recorded rather than quietly
corrected -- a prior `order` result compared the budget-truncation point under one signal while
naming another.

**A deferred quad is `Keep`, not `BudgetExhausted`.** Under `k_frac < 1` a quad that falls down
the ranking was outranked, not refused for want of budget. Conflating them would hide the ranking
inside the stop-reason column that exists to expose it.

**A control that is budget-bound is not a control.** `Mode::Uniform` proves the depth-variance test
discriminates only if the **veto** stops it -- then it reads exactly 0.0000 at every `t`. At a
viewport where the budget runs out first it reads 0.18 and proves nothing. Print
`budget_exhausted` per row. Related: the first cut of the mode test ran at `t = 2`, where
near-field is tame enough that the criterion floors nothing, so **both** arms hit the veto at 256
leaves and variance 0 -- and it read as "balanced degenerated".

**Churn must be read over SHARED quads, and the count printed.** A quad present at one playhead
and not the other has not changed its decision; counting it folds the tree's size change into a
statistic about its stability. Near-field at `t = 16` shares 14 quads, so its churn of 0.4286 is
6 of 14 -- thin, and labelled thin.

**The structure term needs THREE factors and a test found the third.** Connectedness x thinness
scored a **single isolated hot cell at 1.0**: it is trivially the largest component, and
`perimeter_ratio == 4` saturates thinness. Maximum structure, for one cell. **Extent**
(`largest_component / N`) is the graded form of `looks_like_boundary`'s `largest >= N/2`. Measured:
fully hot 0.0000, filament 1.0000, checkerboard 0.0039, isolated cell 0.1250 — each factor catches
a case the other two score at maximum. `NaN` on an empty mask, never 0, because `far`'s mask is
empty on every leaf.

**STRUCTURE NEITHER REPLACES NOR MULTIPLIES — and it wrecks the best criterion.** §2.2's
recommendation on record was multiply. Measured by `error(B)` at levels 6 under the shipping
colouring, on near-field, `deep interior` and `preset_shape`: multiply **never helps**, is a wash
to five digits on the `within` arm (already the worst criterion tested), and takes
`frac_hot_between` on `preset_shape` from **0.07038 to 0.13133** — from greedy's neighbourhood to
worse than the random band's upper edge. `replace` is worse still and is **not a second data
point**: `signal_with(_, _, Replace)` discards both arguments, so `replace x within`,
`replace x between` and `structure_only` are one expression printed three times. As
`structure_only` it is the **worst row in the `preset_shape` table**, beaten by random *high*.

**`frac_hot_between/median` with structure OFF is the answer, on 31 distinct values.** It beats the
random band at nearly every budget in all three targets, and on `preset_shape` — the only tree the
criterion actually controls — it reaches **0.07038 against greedy's 0.06881**. It does this on 31 /
65 / 64 distinct values with modal shares of **83.1% / 33.9% / 40.4%**, while `within/median`'s 5418
distinct values are beaten by random at every budget. *Signal resolution is not what makes a ranking
good*, at full strength.

**And it reads the ABSOLUTE mask** — the one the "make the threshold relative" instruction would
have replaced. The relative mask desaturated the spatial fields exactly as intended, and every
criterion built on it still loses to a saturated 31-valued count. **Desaturating was necessary and
is not sufficient.** `grad_rms`, threshold-free with every quad distinct, sits mid-pack.

**A test whose subject never executes is decoration, and `t = 2` is where that keeps happening.**
Three times in this build a scheduler test was written at a short horizon and read as a failure of
the thing under test: near-field at `t = 2` is tame enough that every leaf reads `Keep`, so the
ranking never runs, no tree ever becomes unbalanced, and a pan cannot change anything. Pin
scheduler tests at `t = 13`. Same family: under the veto near-field reaches a **complete tree at
one depth**, where 2:1 holds trivially — use `deep interior` for the balance test.

**`band_of` conflated `NaN` with `+inf` under one `is_finite` guard.** Undetermined and
maximally-important into the same bucket. Nothing in the current signal produces `+inf`, which is
exactly why it would have sat there unnoticed. `NaN` to the bottom, `+inf` to the top.

**The coarse-ancestor fill was a missing FILTER, not a missing feature.** `adaptive::render` drew
only leaves, so an uncomputed leaf left raw background — a hole, which reads as "nothing here"
rather than "not yet resolved". Drawing every node with samples, coarsest first, IS §4.5's option
1, and it is **bitwise identical wherever the tree is complete** because leaves tile the root. Keep
`LeafTexel` leaves-only: including the fill would double its rows and halve the apparent texel size
at every level.

**The zoom-out assertion cannot be made in the form the brief states.** *"Newly-computed quads
after a zoom-out is ~0"* presupposes a tree persisting across frames, and the scope discipline is
*no eviction, no caching, no async, no promotion*. What is measurable is the arithmetic underneath:
a zoomed-out descent computes 537 quads against a zoomed-in 597, and **zero** of its boxes are
absent from the zoomed-in run — so a persistent tree would compute none of them. Say that, rather
than the claim the build cannot support.

**A DOCUMENTED REPRODUCTION COMMAND CAN BE WRONG, AND ONLY RUNNING IT FINDS OUT.** RESULTS §13's
lines for `pan_sequence` and `slice_gallery` named `9 2000 512` and `4000 ... 512` where the
committed dumps were made at `9 20000 1024` and `40000 ... 1024`. Nineteen dumps failed to
reproduce: nine at a tenth of the budget, and ten with an **identical tree and a different
`decision` column** -- 252 leaves moved from `MaxRelDepth` to `ScreenFloor` purely by the viewport.
Same leaf count, different stop reason. **Verify a regeneration over the WHOLE corpus, never a
sample** -- eleven dumps had already reported "reproduces bitwise" -- and diff `decision`
specifically, because it is the column that moves when a parameter is wrong and the tree is not.

**The corpus was mixed-version and a corpus-wide statistic silently ran on a subset.** `vertical/`
was PRNQ **v1** (24 columns, no hot-mask block) while `charts/` and `criterion/` were v2, so any
statistic over the mask columns covered the v2 dumps only -- and two numbers printed side by side
had different denominators without saying so. Measured both ways: `tau = 1e-4` is at the **0.4th
percentile of `charts/`** (75,359 leaves) and the **4.3rd of the whole corpus** (92,880); the mask
saturates on **98.8%** of chart leaves and **87.1%** corpus-wide. The conclusion survives either
scope, which is exactly why the mislabelling was survivable and had to be caught by arithmetic
rather than by a wrong answer. **State the scope with every figure.**

**THE MECHANISM SHIPPED DISABLED TWICE, AND THE SECOND TIME IT WAS THE DEFAULT.** `k_frac = 1.0`
takes the top 100% of the ranked frontier: `Mode::Balanced` computes the priority, sorts the queue,
and refines all of it. The ranking runs and changes nothing. It was `SchedCfg::default()` through
PR #21, so all 69 dumps in `results/charts`, `results/criterion` and `results/vertical` are the
uniform-mode control with new columns attached, and so was every render made from them. The default
is now `K_FRAC_RANKED = 0.25` and `K_FRAC_UNRANKED = 1.0` is named and kept. **A configuration that
silently reproduces the old behaviour needs a guard, not a convention** —
`scheduler::assert_not_uniform_in_disguise` refuses a `results/` path from the degenerate cell, and
its test asserts all four cells because a guard that always fires passes just as easily as one that
never does.

**The uniform arm was being truncated by the knob it exists to isolate.** `decide` short-circuits
`Mode::Uniform` to `Split` — the criterion is *off*, not permissive — but the rank truncation in
`descend` ran afterwards and demoted the outranked quads to `Keep` anyway. Measured: near-field at
`t = 4` gave **40 leaves and depth variance 0.6900 under both arms, to four digits**, budget never
exhausted. Two arms agreeing to the digit is the same tell as three unrelated charts agreeing.
Exempted, and with it the §5 acceptance test discriminates for the first time: uniform reads
**256 leaves, variance exactly 0.0000, all 256 stopped by the veto**, balanced 0.23–0.69 with churn
0.08–0.53.

**`rho(depth, spread)` is confounded TWICE and both obvious repairs leave one.** Against the leaf's
own spread it reads `-0.817` at `k = 1` and `+0.821` at `k = 0.25` — but refining a quad *reduces*
the spread of the pieces it becomes. Substituting the **parent's** spread removes that arm and
leaves the scale term: `ensemble_spread` is over copies jittered within the cell, the cell halves
each level, and the measured inter-level ratio is 1.19–1.62 — so it comes out **negative at every
`k`, including the ranked ones where the ranking demonstrably works**. The form with neither
confound is **blocked by level**: among the quads at level `L`, did the ones that got split have the
higher spread? Measured near-field at `tau = 1e-4`: **-0.295, -0.028, +0.265, +0.137, +0.108** at
`k = 1, 0.5, 0.25, 0.1, 0.05`. The sign flip survives, at a third of the naive magnitude, and it is
`NaN` on every degenerate row rather than 0 — a level with one outcome has no correlation, and
saying so is the point.

**`k_frac` has an over-sparse end, and 0.25 is where two independent statistics peak.** Near-field
depth variance runs **1.015, 1.900, 2.053, 2.046, 1.334** across `k = 1, 0.5, 0.25, 0.1, 0.05`, and
the tree loses a whole level at `0.05` (4 distinct levels against 5). The level-blocked `rho` peaks
at the same rung. `tau` over a whole decade at fixed `k` moves the variance **1.900 -> 1.866**:
`tau` is a threshold and cannot land inside a distribution whose median moves six orders between
regions, where a rank always cuts through it.

**`k_frac` is a BUDGET-QUALITY TRADE, not an improvement at fixed cost — and the depth-variance
table alone does not say so.** Scored by `Cache::error_of` at each tree's own leaf count against a
five-seed random band, near-field: **no rung is sparse** — every ranked tree is at or below the
band, and the gap to `greedy_oracle` is small (0.07539 against 0.07117 at `B = 40`). But the tree
error rises monotonically as `k` falls, **0.05841 -> 0.07667**, because the tree displays less. The
selectivity is in the *shape* — depth variance, criterion-bound stops, the level-blocked `rho`
turning positive — bought at a real cost in displayed error.

**A leaf set scored against a finite reference must TILE the root, and the mapping between the two
trees is where that breaks.** `Cache::error_of` over a set with a hole is an average over the hole.
Two guards earn their place in `equal_budget`: the descent is capped at the reference's own depth
(without it, 170 of 223 leaves fell outside and the run scored nothing rather than dropping them),
and the recovered `(level, ix, iy)` is checked back against the quad centre. The second caught a
half-cell error — a centre sits at `(2i+1)h`, so dividing by the cell width `2h` gives `i + 0.5` and
`.round()` lands on `i + 1`, mapping every quad to its right/upper neighbour. **Without the check it
would have scored a perfectly coherent leaf set belonging to a shifted tree.**

**A self-describing filename has to carry every setting that is swept.** `criterion_sweep`'s stem
held tau/k/struct/crit and not `alpha_hi`, so stage 3's six alpha rows all landed on the one stem
stage 1 had written and the last writer won. Caught by re-running stage 1 over the committed corpus
and diffing: the `k0.25` dumps came back `alpha_hi=0.2 quads=29` where the committed ones read
`alpha_hi=-1 quads=53`. Same family as the wrong reproduction command, at the filename.

**A pixel count of a debug colour is a fact about the texel size.** `results/glsl/shape.png` carries
1046 magenta pixels, which reads as scattered failure. The adaptive render is nearest-neighbour, so
one footprint of a level-2 quad paints ~16x16: measured, the magenta is **3 axis-aligned blocks** in
`shape` and **1** in `plambda` — four footprints. The census over `preset_shape` at N=64 and N=128
gives **0.244% and 0.201%** undetermined, **zero `SimFailed` and zero `DecodeFailed`** at both, all
of it non-finite copies and non-finite `shape_vec` (a triple collision). Stable across resolution is
the tell that it is the chart, not the grid. It is the instrument reporting; a decode failure would
have been the fault.

**`greedy_oracle` WAS NEVER A BOUND, AND THE NAME COST TWO PRs.** Renamed
`greedy_lookahead_1`. On `far` at `B = 1535` it read **0.54760** against a random band of
**0.48550-0.52047** and every criterion at **0.36557** -- the worst strategy in the table, under a
name asserting it could not be. The real ceiling is `Cache::dp_optimal`, the exact minimum over all
tree-shaped leaf sets, and **it is cheap**: 5461 splits at `levels = 7` in **0.01-0.03 s** against
296-1487 s to build the cache. The naive `O(quads x B^2)` reading is wrong twice -- the 4-way merge
is three 2-way convolutions, and each node's cap is bounded by its own subtree, so only the top two
levels see the full budget. State the three roles wherever a curve is quoted: **floor** = random
band, **reference** = `greedy_lookahead_1`, **ceiling** = `dp_optimal`.

**THE ACCOUNTING IDENTITY IS A TAUTOLOGY, AND IT WAS PROPOSED AS THE FIX.**
`error_of(leaves) == err_sum(root) - sum(gains)` telescopes directly from `error_of` being a sum of
`err_sum` and `gain` being parent-minus-children. It holds for **any** ranking, any sequence, and
any values `err_sum` holds -- random numbers included. So does "assert greedy picked the argmax",
which re-runs `replay_with_leaves`'s own argmax over a pure, static `gain`. Both report PASS and
read as clearing the metric. *A test that cannot fail is indistinguishable from a test that passes*
applies to the **metric**, not only to the physics. The form with teeth is the exact optimum: *no
ranking may beat it*. Measured worst margin `+0.0e0`, `-1.4e-16`, `-1.4e-17` across three regions
-- **the replay is sound and no `error(B)` number is suspect.**

**ON A SMOOTH FIELD, RANKING BY SPREAD *IS* BREADTH-FIRST — so `far` degenerating is CORRECT.**
A quad of width `w` over a gradient `g` has variation `~ g*w`, so spread tracks cell size and
argmax-on-spread picks the shallowest quad. A constant signal has no argmax and falls through to a
tie-break that is lexicographic on `(level, ix, iy)` -- level first, also shallowest. Both are
breadth-first, which on a smooth field is within **2e-5 of exactly optimal**. That is why `far`'s
thirteen non-greedy rows agree to five digits: `within/median` on **21845** distinct values and
`frac_hot_between` on **1** produce the **identical leaf set** `5:982 6:168`. One allocation, two
routes -- not thirteen criteria agreeing. **A criterion can only differ from breadth-first where
the field is not smooth, and smoothness is where there is nothing to find.** `far` is the control
that shows what a featureless field looks like.

**A LEVEL BARRIER IN `err_sum` DEFEATS GREEDY, AND THE SAME STATISTIC EXPLAINS ITS ABSENCE.**
`far`'s `err_sum` per pixel is **flat through level 3** (1.00000, 1.00067, 1.00016, 0.99884) and its
gains there are noise -- `-3.0e-7`, `6.3e-7`, `-6.4e-9` with **13 of 16 negative at level 2** --
against `8.6e2` at level 3. Greedy will not pay a negative-gain split to unlock a gain five orders
larger beneath it, so it opens one subtree and descends to the bottom inside it: **14 leaves left at
level 2** while 1024 reach level 7. `near-field` and `deep interior` fall gradually at every level,
have no barrier, and greedy is **near-optimal** there (gap 0.0004 and 0.0301 against `far`'s
**0.7693**). Read `err_sum` by level before diagnosing any allocation failure.

**A SPLIT CAN MAKE THE IMAGE WORSE; TAKE THE PREFIX-MIN.** A parent's `N x N` sample grid and its
children's are **different approximation families**, not a nested refinement of one. Measured: the
root's own gain on `far` is `-3.022e-7` -- one 8x8 grid over the whole box beats four 8x8 grids over
quarters. Negative gains number 14 / 14 / **102** in `far` / `near-field` / `deep interior`. So the
ceiling at budget `B` is `min over s <= (B-1)/4` of `f_root(s)`, never `f_root(S)`, and where the
prefix-min binds is a direct measurement of negative gain.

**READ THE LEAF HISTOGRAM, NOT ONLY `error(B)`.** The error cell says a criterion lost; the
level histogram says how. `within/median` in `near-field` at `B = 1535` leaves **2 leaves at level
1** and drives **984 to level 7**, against an optimum with nothing below level 3 and only 16 at
level 7 -- its allocation is *inverted*, not merely worse. No error cell showed that, and the
headroom it implies is real: the best criterion leaves **5.1% at `B = 1535` and 16.3% at
`B = 12287`** of achievable improvement in `near-field`. Against greedy that figure read
*negative*, which is why the question was unanswerable by inspection. **This is the strongest
finding in that run and the strongest argument for changing the default -- stronger than any error
digit.** Those percentages are `frac_hot_between` against `dp`; with `Rank::Uniform` in the table the
best *row* forfeits 3.4% and 12.4% at the same budgets, and it is uniform.

**"BEATS RANDOM" IS THE WRONG BAR; THE BASELINE IS BREADTH-FIRST.** `Rank::Uniform` was missing from
every table in the corpus, and against it most of the standing comparisons change sign.
`frac_hot_between/median` -- the best criterion measured on this project -- **never beats uniform in
`near-field` at any budget** (`0.11148` against `0.10984` at `B = 1535`), while sitting well clear of
the random band. In `far` uniform **is** the exact optimum. Only `deep interior` at `B >= 6143` shows
a criterion decisively ahead (`0.04035` against `0.04813`). Random is a floor no strategy should be
below; it is not the thing a criterion has to beat.

**THE HEADROOM RISES WITH STRUCTURE, AND A CRITERION ONLY EARNS ITS KEEP WHERE STRUCTURE IS
LOCALISED.** Share of achievable improvement the best row forfeits at `B = 1535`: `far` **0.0002**
(nothing varies), `near-field` **0.0336** (structure localised), `deep interior` **0.0999**
(structure everywhere). Where nothing varies there is nothing to rank; where everything varies the
budget must go everywhere and breadth-first is right again. So `far` degenerating is **correct
behaviour** -- it is the control that shows what a featureless field looks like, not a region where
the criteria mysteriously tie.

**QUOTE THE GAP AS A CURVE; ONE BUDGET SUPPORTS A DIFFERENT STORY IN EVERY REGION.** `near-field`
rises `0.0336 -> 0.1106 -> 0.1241` over `B = 1535, 6143, 12287` -- the criterion is adequate early
and progressively worse late, which fits the inverted histogram: the failure is **late-stage
allocation**. `deep interior` peaks at `B = 3071` and falls; `far` collapses to zero by `B = 1535`.
Same statistic, three shapes.

**THE CRITERION-TO-UNIFORM GAP RISES WITH DEPTH -- AND THAT IS THE ARGUMENT AGAINST A DEPTH
PARAMETER, NOT FOR ONE.** `captured = (uniform - row)/(uniform - dp)` by level in `deep interior`:
**-9.38, -1.97, -0.45, -0.04, +0.53** at levels 2-6. Shallow splits the criterion chooses are
actively worse than raster order. But the crossing sits at level 6 there, is **never reached** in
`near-field`, and is undefined in `far` -- so any fixed depth is a tunable constant that must be
retuned per region, the same defect that killed a fixed `tau`. And the ranking **already
self-adapts**: on a flat signal it ties, ties fall to a level-first tie-break, and that *is*
breadth-first. Protect that property rather than overriding it.

**AMPLITUDE CANNOT TELL A SMALL REAL SIGNAL FROM NOISE; COHERENCE CAN.** The AUTO-RANGED OVER NOISE
guard's two arms both read amplitude, and its absolute arm compares against the region's own median
energy drift -- so in a tame region the floor falls with the field it is meant to bound, a ratio in
disguise. `far` cleared it at `ramp.1 = 1.064e-8` against a floor of `4.478e-9`. The discriminator
is **spatial coherence**: lag-1 neighbour correlation of the ramped scalar reads **0.9434 -> 1.0000**
across levels on `far`, and its p1/p99 halve **exactly** (ratio 2.000) per level -- `spread ~ g*w`
measured, a real gradient of tiny magnitude. `far` is a smooth field, not an amplified noise floor,
and the level-2 `err_sum` barrier is a property of the field. A third guard arm now reads coherence.

**AN OUTPUT ROOT IS AN ARGUMENT, NOT A CONSTANT.** `criterion_metric` wrote to a hardcoded
`results/`, so the reduced-`levels` validation pass that exists to fire the `dp_optimal` assert would
overwrite the committed 512^2 artefacts with a small raster -- the failure already on record, at a
second site. The root is now argument five, defaulting to `results`. **And an assert that has never
run is a gap that survives indefinitely**: the bound assertion was executed at scale, worst margins
`-1.4e-17`, `-5.6e-17`, `+0.0e0` across three regions.

**A DEGENERATE DENOMINATOR IS THE FINDING, NOT A RATIO.** Where uniform *is* the optimum -- what a
smooth field looks like -- `uniform - dp` is a rounding epsilon of either sign, and `captured` printed
**-0.6786** on `far` from a denominator of `-1.55e-15`. It reads as "the criterion is 68% worse" and
means "the two are identical to machine precision". Print the identity.

**`t_end` IS QUANTISED WHERE ESCAPE TERMINATES, AND IT RENDERS AS CONCENTRIC CONTOUR BANDS.**
Collision is sampled inside the RK4 loop (`tc = t + s.t`) and carries step resolution; **escape is
sampled only at sync boundaries**, the reference's cadence. So `t_end` takes `n_sync` values across
a whole chart wherever escape is the terminating event, and every derived field draws those steps.
Measured at `escape_every` 0 -> 1: `preset_plambda` **16 -> 2623** distinct with **99.52% -> 0.26%**
landing exactly on a boundary; `preset_shape_pl_h1` **41 -> 2316**. **`near-field` is bitwise
unaffected at every stride** -- its escape arm is silent at `t = 13` and it terminates by collision,
so there is nothing to quantise. The mechanism predicts which images band, which is what makes it a
mechanism. `AzOpts::escape_every` defaults to `0`; turning it on is a **spec change**, because the
cross-check and the horizon table were both measured coarse.

**A COUNT CAPPED BELOW ITS OWN DECISION THRESHOLD CANNOT DECIDE ANYTHING.** The proposed test for
that artefact was *"recount `frac_hot_between`: 45 -> thousands means quantisation, 45 -> 45 means
Wada"*. But `frac_hot_between` is `frac_above_tau_between`, a fraction over `N^2` footprints, so at
`N = 8` it has **at most 65 distinct values by construction** -- the corpus's own `31 / 65 / 64` is
that ceiling showing. It would have reported "Wada" under either hypothesis. **Check a statistic's
arithmetic range against the threshold before using it as a test.** The count that can fire is
`t_end` itself, which is unbounded and is the quantity actually being quantised. Measured anyway,
`frac_hot_between` moves the *wrong way* -- more saturated under the finer cadence, 31 -> 12
distinct with modal 41.2% -> **83.5%** -- so the saturation is not quantisation.

**"ON A BOUNDARY" CONFLATES QUANTISED WITH FINISHED.** `near-field` reads **97.85%** of `t_end`
landing on a sync boundary while being completely clean, because 97.8% of its footprints are
*Bounded* and sit at `t_end = t_max = 13`, itself a multiple of the sync interval. The horizon is a
boundary time. **Read the delta across cadences, never the level.**

**THE CADENCE'S REAL EFFECT IS AN OUTCOME RE-LABELLING, NOT A RESOLUTION GAIN.** `deep interior`'s
`t_end` was already continuous (2983 -> 3303, 1.1x) but its terminal class moves **escape
0.0945 -> 0.5494, collision 0.8965 -> 0.4482** -- half the region. Under `stop_on_event` a genuine
early escape is not *noticed* until the next boundary; the run keeps integrating, dips below
`r_coll`, and a collision wrongly wins precedence. A different defect from the one the banding
pointed at, found by the same measurement.

**AND THE `d_min` DISCRIMINATOR SAYS IT IS A FIX, NOT SPURIOUS FIRING.** Testing escape inside the
RK4 loop also asks it *during* a close encounter, where a pair's instantaneous two-body energy can
read positive transiently -- which would re-label the footprints with the **smallest** separations.
The re-labelled ones have the **largest**: escaped `d_min` p50 **1.063e-3 -> 4.419e-3** and p10
1.891e-4 -> 1.164e-3, while collided p10 rises 2.734e-4 -> 7.016e-4. Tightest approaches stay
collisions; marginal ones just under `r_coll` become escapes. It converges between strides 4 and 1.
**Not tested:** whether a fired escape persists to the following boundary.

**THE BANDING IS A COLOURING ARTEFACT; THE CRISP EDGES ARE NOT.** Under outcome-class colouring
(`preset_shape_pl_h1_uniform_outcome.png`, already committed) the arcs vanish entirely while the
straight edges survive and sharpen -- including a genuinely **circular** boundary around the
central fan. Outcome-class boundaries are real regime structure. A circle plus radiating wedges is
**polar structure in the chart plane**: a radius threshold and an angle threshold, with saturation
the candidate. Stated and stopped there -- two image diagnoses on this project were settled by one
targeted measurement each, and speculating past the render is how the earlier ones went wrong.

**`level` NEARLY SOLVES THE DP-LABEL TASK, SO EVERY CORRELATION AGAINST IT MUST BE BLOCKED.** The
optimum splits shallow quads and keeps deep ones, so `level` scores **|rho| = 0.993** against the
label pooled and every signal tracking cell width inherits it. Blocked within a level,
`spread_median` **flips sign**: pooled +0.225 -> +0.488 across the budget ladder, blocked
**-0.173, -0.098, -0.134, -0.257**. Same confound and same repair as `rho(depth, spread)`. A
logistic fit shows it too -- held-out AUC **0.88 -> 0.37** when `level` and `cell_width` are
removed, so a fit carrying them is a depth model wearing 55 names.

**THE DP's LABELS ARE MONOTONE IN BUDGET, so a fixed ordering is not structurally excluded.**
Measured on `near-field` at `B = 169, 681, 2729, 10921`: **`1 -> 0` flips are exactly zero** at
every rung, with 91 / 241 / 1049 flips all `0 -> 1`. The optimum never un-splits, so its split sets
are strictly nested -- and nested sets are precisely what a single fixed priority order produces as
prefixes. The per-level `captured` result is **not** a symptom of budget-dependent labels. The
shared population is the whole smaller tree, not a thin slice.

**A QUAD OUTSIDE THE OPTIMAL TREE WAS NEVER DECIDED.** `Dp::labels` maps internal -> split,
leaf -> keep, and **absent -> absent**. Labelling an unreached quad `keep` invents tens of thousands
of labels: at `B = 2729` the tree holds 2729 nodes of 21845. Report the population with any
statistic taken over it.

**THE ESCAPE CONDITION IS NOT ABSORBING, AND A FINER TEST LATCHES TRANSIENTS.** `escape_candidate`
is relative energy `> 0` and receding, which during a close encounter is transiently true. Of the
**895** `deep interior` trajectories that escape under an in-loop test and not at the reference
cadence, **0.000 are still unbound one boundary later** -- and 0.000 at +2, +3, +4 and +8. All 895
are transients, and latching them took the escape fraction from 0.0947 to **0.5494**. Where escape
genuinely terminates (`preset_plambda`) the finer stride adds **zero** new escapes and only sharpens
the time. `AzOpts::escape_confirm` holds an in-loop detection provisional until the next boundary
and commits the **first crossing** time; guarded, `deep interior` reads 0.1564 and converges at
stride 4, and `preset_shape_pl_h1`'s labels are **stride-invariant** while its `t_end` resolution
improves 56x. **That makes the stride a COST knob, not a correctness one.**

**A `d_min` READ OVER THE WHOLE RUN CANNOT TEST WHETHER THE RUN STOPPED TOO EARLY.** The first
attempt at fix-or-bug compared `d_min_true` by terminal state and found the re-labelled footprints
carried *larger* separations (p50 `1.063e-3 -> 4.419e-3`), read as "not mid-encounter firing". But a
run stopped early by a spurious escape never reaches its close approach, so its `d_min` is larger
**because** it terminated early. The statistic was confounded by the effect it was measuring, and
its answer was the opposite of the truth. **Test the mechanism directly; a summary statistic over a
truncated run inherits the truncation.**

**AND THE FIRST PERSISTENCE SWEEP WAS CONFOUNDED THE SAME WAY.** It re-ran to `t_e + w` with
`n_sync` rescaled per window, so every window was a different discretisation, and produced
`0.162, 0.219, 0.011, 0.083, 0.335` -- read as "the condition flickers". Recording candidacy at
every boundary of **one** run at **one** step size (`AzOut::escape_flags`) gives a flat zero.
*`n_sync` fixed while `t_max` varies compares different discretisations* is the standing form; this
was its inverse, inside a diagnostic written to catch exactly this.

**"COLLISION IS SAMPLED CONTINUOUSLY SO IT FIRES FIRST" CONFUSES DETECTED WITH OCCURRED.**
`classify` ranked collision above escape unconditionally and **discarded both times**. Collision is
sampled inside the RK4 loop and escape only where the state is Cartesian, so an escape at `t = 5.0`
noticed at `t = 5.28` lost to a collision at `t = 5.1`. Deciding by `min(t)`, with `t_end` set the
same way so state and time cannot disagree, removes the dependence on *when each arm is sampled*.
Measured alone it moves **one footprint of 5440** in production, because `stop_on_event` breaks on
the first *detected* event so only one is ever recorded. On the reference path it is large:
`preset_plambda` **990 of 996** trajectories that fired both arms escaped first and were labelled
collision -- **42.97% of the whole slice** -- median lead 6.45 sync intervals. **The ordering
guarantee comes from sampling both arms at the same cadence; `min(t)` is what stops the state and
the time disagreeing.**

**A GUARD NEEDS THE ARM THAT SAYS IT DID NOT CUT TOO MUCH.** `escape_confirm`'s test asserts it
**reduces** the escape count in `deep interior` *and* leaves `preset_plambda`'s **unchanged**. A
guard that rejects everything passes the first arm exactly as well as a correct one.

**THE CLOSURE CRITERION PASSES THE CHECK THAT KILLED THE LAST ONE, AND ITS HEADLINE GAP DOES NOT
REPRODUCE.** `|dn| over a window < tau AND E_rel > 0`, transcribed from
`reference/escape_criterion.py`. Of the trajectories it fires on, **1.0000 are still unbound at +1,
+2, +4, +8 boundaries and at `3 t_max`** in `deep interior`, `preset_plambda` and
`config_stability` -- against **0 of 895** for `spec > 0 && receding`. But the reference's 383x
separation between escapers and bound trajectories is not reachable here: the best cell is **6.8x**
(`deep interior`, `2 t_max`, `k = 2`) and **`near-field` shows none at all** (0.9-1.1 at both
horizons, 0 fires of 576). Three separable causes, none decided: maturity (`|dn/dt| ~ 1/t^3`, the
reference quotes `t = 25-30` and this ships at 13, and every region that separates separates better
at `2 t_max`), population (a geometric ground truth counts a **triple dispersal** as an escape and
nothing converges to a pole in one), and `near-field` simply not having escaped yet.

**A WIDER CLOSURE WINDOW MAKES THE GAP WORSE.** `k` from 1 to 4 runs `deep interior` 6.1 -> 4.0 and
`config_stability` 0.6 -> 0.4; `k = 1` is best or joint-best in every region that separates at all.
The window is a **TIME**, so `n_sync` must scale with `t_max` -- at `t_max = 50, n_sync = 32` it is
1.5625 against the reference's 0.4, a different criterion wearing the same name.

**THE CLOSURE WINDOW CANNOT RESOLVE INNER-BINARY PHASE, ANYWHERE.** `t_close = 2 pi sqrt(d_min^3/M)`
runs **17x to 274x below the window** in every region, and a two-end chord cannot tell a full
revolution from stationarity -- the reference buffers `nbuf` samples and reads only `buf[-1]` and
`buf[0]`. So the closure arm is structurally blind to a tight bound pair and rejecting one rests
**entirely on the energy arm**. Transcribed, not a defect, and the reason neither arm is redundant.
`tests/outcome_encoding.rs::a_full_revolution_aliases_to_zero_closure` holds it as a property.

**PERSISTENCE IS THE ENERGY ARM, NOT FULL CANDIDACY.** Closure is a difference of neighbouring
samples, so it jitters above `tau` on a perfectly settled escape: `preset_plambda` reads **0.4777**
still-candidate at the last boundary against **1.0000** still-unbound. Reading the 0-of-895 test off
candidacy would have scored ordinary jitter as a re-binding and reported a correct criterion as
broken. `AzOut::unbound_flags` exists for that question; `escape_flags` is rule-aware and answers a
different one.

**AND PRECISION MUST BE READ AGAINST IT.** `deep interior` reads precision 0.3214 while **every** one
of its fires is still unbound at `3 t_max`. The geometric ground truth demands 3x separation growth
by then, so a slow genuine escape is not certified -- a precision shortfall with persistence at
1.0000 is the **ground truth missing them**, not the criterion inventing them.

**THE `t_end` REPLAY REFINEMENT IS MEASURABLY DECORATION -- `at entry` is 1.0000 everywhere.** The
energy arm always already holds when closure settles, so there is never a crossing inside the firing
interval to find. `t_end` under this criterion is irreducibly quantised to the boundary cadence,
because closure is *defined* from the boundary series and has no finer resolution. Kept with the
counter, because the counter is what says it is inert -- and it returns the **boundary** time, not
the entry time, since the entry time would claim an escape at a playhead the criterion had not yet
concluded one at.

**A COLLIDED RUN FREEZES AND ITS CLOSURE READS EXACTLY ZERO.** Under `stop_on_event` the shape stops
changing, so a collided trajectory is the most "settled" thing in the sample -- and it lands in the
**bound** population and destroys the gap it is supposed to measure. Measure the signal with nothing
terminal and split collisions out. Exactly-zero closure survives even then: `deep interior` 10 of 63
bound, `config_stability` 4 of 248, **0 of every escaper population**. Count it; do not let a
percentile absorb it.

**THE ESCAPING-BODY LABEL IS THE LOWEST FIRING INDEX, NOT THE ESCAPING BODY.** The reference's
`b = np.argmax(fire, -1)` returns index order, and on a dispersing system all three bodies read
unbound -- checked against the reference itself, `E = [9.474, 12.178, 31.825]`, `argmax = 0`. So
`detail`, which classification and rendering read, says 0 where the physics says 2. Transcribed, not
corrected. **14 of the 40 golden rows discriminate it**, and a tightest-pair-first ordering fails on
them.

**THE TOGGLE'S MEDIAN IS EXACTLY ZERO AND 392,466 OF 1,048,576 PIXELS MOVED.** `stop_on_escape` on
`preset_plambda` at 1024^2: median `0.000e0`, **37.43% of pixels move, worst 5.993e-1** -- a third of
the shape sphere's diameter. *Never conclude "no effect" from an aggregate without the per-pixel
distribution*, now three times, and this time on the run written to **confirm** the prediction. And
**the small grid understated it by eight times**: at `n = 24` the worst was 7.4e-2. Read the max and
the moved count at the resolution that ships.

**CLOSURE DOES NOT CERTIFY THAT THE DISPLAYED QUANTITY HAS SETTLED, AND THE RENDER SAYS SO.** Under
`stop_on_escape = false` `preset_plambda`'s ribbons run continuously to every frame edge; under
`true`, same criterion and same physics, **two large smooth arcs sweep up from the lower corners and
cut through the ribbon structure**. The domes regrow. The criterion fires at a median `t = 11.8` of
13 with persistence 1.0000, and the shape still moves by up to 0.6 in the remaining 1.2 time units.
`|dn/dt| ~ 1/t^3` bounds the **rate**; a small rate integrated over 1.2 time units is not a small
displacement. **The criterion is right about what escaped and silent about whether the shape has
stopped** -- those are different questions and one window cannot answer both. `stop_on_escape` stays
**off**, and the `_outcome.png` pair is **bitwise identical**, so the toggle changes only *when*
`shape_vec` is read and the comparison is clean by construction.

**`escape_confirm` AND `escape_every` ARE VACUOUS UNDER `Closure`.** Closure is only defined where
the state is Cartesian, so there is nothing to test between boundaries, and the window already is
the persistence guard. Stated, not silently dropped -- and the replay refinement is scoped to
`Closure` alone, because under `Reference`/`Distance` both arms are continuous and refining on
energy would return a time at which the full condition may not hold.


**`dtau` SIZED ONCE PER SYNC INTERVAL BLOWS UP AFTER A BOUNDARY-COINCIDENT ENCOUNTER.**
`dt = A*B*dtau`, so the physical step is `eta*dt_left` only while `A*B` stays near its entry value.
A trajectory at a close encounter **at a boundary** has a tiny `A0*B0`, so `dtau` is enormous and
`dt` grows by orders as the bodies separate. Not "close encounters" — encounters *coinciding with a
boundary*, a thin set, which is why the damage clusters spatially and a `d_min` correlation is flat
(`1.04e-04` in both populations). `DtauMode::PerStepInterval` recomputes `A*B` per step with
`dt_left` **held fixed** and caps at the entry value; the cap is one-sided in the right direction —
when `A*B` grows the value falls and the blow-up goes, when `A*B` falls at a close approach the cap
holds `dtau` at nominal so `dt` shrinks with the separation, which is what regularisation buys.
Cost is **~10% more steps** and the budget count does **not** rise (`deep interior` 3 -> 1). It
removes the over-correction without reintroducing what was over-corrected for. `tb_az.py` carries
the same three modes; **fix both sides or the cross-check disagrees for the right reason.**

**PUTTING THE REMAINING TIME IN THE NUMERATOR IS ZENO BY ARITHMETIC.** `dtau = eta*(dt_left-s.t)/(A*B)`
gives `dt ~ eta*rem`, so `rem_{n+1} = rem_n (1-eta)` and the interval is approached geometrically and
**never completed** — measured `t/t_max` of 0.0833, 0.0303, 0.0303, **0.0080**, with the whole budget
burned on 2304 of 2304. **And its drift is the best in the table by ten orders** (`1.3e-14`) because
it went nowhere. Kept as `DtauMode::PerStepRemaining`, a named axis and never a candidate. Print
`t/t_max` before any drift column: *a difference can be small because both sides are right or
because one side is dead*, and this arrived inside the diagnostic written to catch it.

**THE 109 NEAR-FIELD ESCAPES AT `t = 20` WERE THE `dtau` BUG.** The standing finding "zero of 1024
fire at `t_max = 13`, 109 at `t_max = 20`" is corrected: under per-step stepping **zero** fire at
`t = 20`. The discriminator is on the trajectories — the 109 that fire have median energy drift
**1.147**, which is 115% of the total energy, against **6.2e-5** for the 915 that stay silent, and
under the fix those same pixels sit at 1.6e-3 and do not escape. A giant post-encounter step throws
a body outward, it reads as unbound and receding at the next boundary, and the arm latches. Genuine
escape is **not** suppressed: at `t = 40` the modes give 280 and 308 of 1024. Burrau's escape is
simply later than 20.

**`TINY*TINY` UNDERFLOWS AT f64 TOO, NOT ONLY f32.** `1e-300 * 1e-300 = 1e-600` is zero at f64, so
the doubly-degenerate hole the standing note describes as an f32 property is open at **both**
precisions. Measured: `ab_floored` fires on one trajectory of 2304 with raw `A*B` reaching exactly
`0.000e0`. Under `PerStepInterval` the `.min(dtau_entry)` cap absorbs the resulting `inf`; under
`FixedPerInterval` there is no cap and nothing to absorb it. Know which guard is doing the work.

**THE CLUSTERING RATIO ROSE WHERE THE PREDICTION SAID IT WOULD FALL, AND THE FIX STILL HELD.**
Observed neighbour fraction over chance `1-(1-base)^4` went 1.164 -> 1.650, 1.000 -> 1.164 and
12.133 -> 16.154 under the fix. The counts fell with it — `n_hot` by 15-47%, non-finite 11 -> 2,
42 -> 15, 7 -> 5 — so what is removed is the *diffuse* high-drift population and what survives is
the genuinely clustered core, a higher ratio on a smaller set. **A ratio and its base rate move
together; read both or the ratio reads backwards.** And in `deep interior` it went *down*
(1.000 -> 0.982) for a fourth reason: at **92% hot** the mask is saturated, chance is ~1 and the
ratio cannot say anything — the standing regional mask-saturation result, at a third statistic.
`nf w/ hot nbr` is worse still: **1.0000 in every cell of the table, before and after**, so the
6-of-6 the mechanism was first quoted from is saturated and could not have come out otherwise. Same defect one level down: ranking the worsened trajectories by `after/before`
returned pixels whose baseline was accidentally excellent (`9.2e-10 -> 4.5e-4`, ratio `4.9e5`, still
the region's best) while the region's `drift max` *improved* 36x. Rank by the absolute rise.

**A MEDIAN CONDITIONED ON EACH ARM'S OWN SELECTION IS A DIFFERENT STATISTIC IN EACH ROW.**
`med(drift | drift > 1e-6)` read the `dtau` fix **backwards** — 3.53e-4 -> 4.44e-4 — because the
fix removes pixels from the selection. Over the pixels hot under **either** arm, the same set in
every row, it falls **236x, 2400x, 1275x and 123x**.

**THE `dtau` FIX REMOVES 169 OF EVERY 170 MAGENTA PIXELS AND TRADES NOTHING FOR THEM.**
`config_stability` at 1024^2: non-finite **30109 -> 178** (2.87% of the frame to 0.017%) with
`budget_exhausted` **0 on both sides**, so the failure was not swapped for a differently-coloured
one. 348,314 of 1,048,576 outcome labels flip, concentrated where the speckle halos were; the
bounded regions become solid and the ribbons continuous. What remains -- red regions stippled with
blue and green -- is genuine fractal mixing: no magenta in it, not ringing a core, stable across
resolution. **`shape d` between the two modes is median 0.111 with a max of 2.000, the full
diameter of the shape sphere, over 909,184 pixels. That is not evidence of anything** -- the two
modes integrate different trajectories through a chaotic region and must diverge. The images
answer whether the structure got cleaner, not whether the pixels agree.

**WHEN A NUMERICAL DEFECT IS SUSPECTED, RENDER THE DIAGNOSTIC FIELD, NOT THE SCIENCE FIELD.**
`Scalar::Drift` / `colour::drift_rgb`: `energy_drift_max` on an inferno ramp, magenta for the same
veto set `colour::rgb` applies, auto-ranged p2-p98. It is already in the payload, so it is a
colouring and not a computation. The science fields show a defect only *after* it has propagated
into a spread or a label; the drift map shows it at source, as coherent arcs with the non-finite
pixels sitting inside them. Write it for **both** arms of a before/after — the "before" map is the
artefact worth keeping, because it is what the signature looks like.

**A FIXTURE THAT HAS TO BE MEASURED ONCE HAS TO BE MEASURED EVERY TIME THE PHYSICS MOVES.** The 2:1
balance test's region has now swapped **twice** — `deep interior`, then `near-field` when the escape
distance gate flattened `deep interior`, then back to `deep interior` under the `dtau` fix (near-field
is now gap 1 at **all twenty-four** swept cells of `alpha_hi x tau x n`). Each time it was the
*control* arm — "the unbalanced tree must actually violate 2:1" — that caught it, not the property
under test. Same session: `escape_matches_the_legacy_classifier` went vacuous because its arm stopped
firing at `t = 20`. **The assertion that the test is exercised is the part that keeps working.**

**THE `dtau` FIX SHIPPED WITHOUT ITS PARTNER, AND ON ITS OWN IT MADE THE IMAGES WORSE.** The march
exits a sync interval by **overshooting** it and only the *clock* was corrected -- the Cartesian
state written back was the overshot one. A first-order error at every boundary inside an RK4 march.
Under `FixedPerInterval` `dtau` is constant across the interval, so the overshoot is a fixed slice
of time, neighbouring trajectories overshoot alike and the error is large but spatially **smooth**.
Under `PerStepInterval` the last step's size is a function of local `A*B`, so the overshoot becomes
a function of local state and neighbouring pixels overshoot by different amounts -- a
spatially-varying error. `AzOpts::clamp_final_step`
(default **on**) lands the final step on the boundary; it is applied **after** `dtau_for_step`
returns, so it composes with every mode rather than being a fourth one, and it reuses the **same
floored `A*B`** the mode used. Measured at **1024^2**, the shipping resolution: under the clamp,
switching `dtau_mode` moves the field **2.5x** less (chord p50 `1.113e-1` -> `4.367e-2`).
**Neither change ships alone.**

**AND THE COARSE GRID OVERSTATED THAT BY 26x -- the same defect as understating it, in the other
direction.** On 48x48 the same ratio reads **66x** in `config_stability` and **316x** in
`near-field`, because 2304 samples over the window are dominated by the tame majority while a
million land in the chaotic population, where any step-control change diverges regardless. *Read
it at the resolution that ships* now has both signs: §23's coarse grid understated a max eightfold;
this one overstated a median twenty-six-fold. Per-trajectory statistics (drift, steps) do not have
this problem; **no chord ratio may be quoted from a coarse grid.**

**AND `moved` ORDERS THE PAIRS BACKWARDS.** At 1024^2 `B->D` moves the **most** pixels (0.9343) and
displaces them the **least** (1.242e-2), while `A->B` moves the fewest (0.8671) and displaces them
the most (1.113e-1); the whole table sits in 0.867-0.934. It counts pixels differing in the last
bit, which on a chaotic field is a fact about the field. §23's 87% is that statistic -- correct, and
answering a different question than it was read as answering. `chord p50` is the discriminator, and
`chord max` is **2.000, antipodal, in every pair**.

**AND THE APPEARANCE THAT PROMPTED IT WAS NOT CAUSED BY IT.** The nested-arc banding in
`config_stability` is present in **all four arms**, including `fixed`+overshoot, which predates both
changes; under outcome-class colouring arm D's arcs **vanish** while the region boundaries sharpen.
So it is §21's standing result at a new site -- *the banding is a colouring artefact; the crisp
edges are not* -- and the spatially-varying-overshoot mechanism, though real and measured, is not
what draws them. What the two changes **do** remove is the magenta: 30109 -> 2071 (clamp alone) ->
**178** (both), with `simfail` 0 throughout. *A finding read off a render is a finding about an
appearance*: taking it seriously found a genuine first-order defect, and the defect was not the
cause. Record both halves.

**READ THE ORDER, NOT THE ERROR -- AND THE FIGURE-EIGHT IS THE INSTRUMENT.** Chenciner-Montgomery
is exactly periodic, so `|state(T) - state(0)|` is a pure error with no reference trajectory and no
chaos, in under a second. Convergence order across `eta in [0.02, 0.001]` at `n_sync = 32`:
`fixed+overshoot` **1.13**, `perstep+overshoot` **1.06**, `fixed+clamp` **3.06**,
`perstep+clamp` **2.08**; the error at `eta = 1e-3` falls **827,000x**. The per-rung two-point
estimates are noise -- `fixed+clamp` runs 2.34, 2.58, 1.36, 6.49 -- so quote the endpoint slope.
**`perstep+clamp` lands at 2, not 3**, because the clamp sizes the last step from the
*instantaneous* `A*B`, a first-order predictor of the time increment, so the landing residual is
`O(h^2)`; stated rather than smoothed.

**A DIAGNOSTIC FIELD IS SPECIFIC TO A CLASS OF DEFECT, AND ENERGY DRIFT IS BLIND TO THIS ONE.**
The clamp buys 24,000x on the figure-eight while moving `near-field`'s median drift **37x the WRONG
way** (1.5e-9 -> 5.6e-8) and the NumPy smoke median 3.197e-9 -> 4.047e-9. The overshoot displaces
the state in *time* and the AZ energy is nearly stationary along the flow. *Render the diagnostic
field, not the science field* is right and incomplete: ask what the diagnostic would say about the
defect you are hunting **before** reading it as clean. What the drift table does confirm is the
prediction's last clause -- D beats B, `drift p99` 6.077e5 -> 2.881e4 in `deep interior`.

**AN ABSOLUTE TOLERANCE IN A SCALE-INVARIANT CODE IS A BUG, AND THE BITWISE GAUGE TEST IS WHAT
CATCHES IT.** The clamp's landing test first used `T::SYNC_EPS` absolute; all times rescale by
`alpha^{3/2}`, so at `alpha = 0.25` the same slack is 8x wider in relative terms and a rescaled twin
lands one step earlier. `shape_spread_is_invariant_under_the_scale_symmetry` asserts **bitwise**
equality and fired at `4.24e-15`. `Real::LAND_EPS_REL` is relative to `dt_left`. Related: the
relative tolerance also gives the Zeno mode a floor -- `PerStepRemaining` now *completes*, after
`ln(1/eps)/eta ~ 3200` steps per interval against a nominal ~100, and the test that asserted it
stalls read as a failure until it was pinned to `clamp_final_step: false`.

**`config_stability` IS RESOLVED AT THE SHIPPING SETTINGS; `deep interior` IS NOT AND NEVER WILL
BE.** Arm D at `eta`, `eta/2`, `eta/4`: `config_stability`'s median inter-rung displacement falls
**4.98x per halving** from `2.1e-5` on a diameter-2 sphere, consistent with order 2.08, and
`near-field` 4.01x. Horizon 50 was a live worry and is not the answer. `deep interior` falls 8.6x
but from `4.7e-2` with `chord max = 1.999` -- antipodal -- at every rung: chaotic divergence over
`t = 13`, which no step size buys off. A difference can be large because the physics is.

**A FIXTURE THAT HAS TO BE MEASURED ONCE HAS TO BE MEASURED EVERY TIME THE PHYSICS MOVES -- now
four times.** The clamp flattened `near-field`'s scheduler tree (184 leaves, no `Split`, levels 2-4
`Keep` and level 5 `ScreenFloor`), so `a_pan_is_an_identity_until_camera_bias_is_switched_on` stopped
firing: a pan of `0.04` against `half_world = 0.05` no longer moved any decision. Measured across
`{0.01 .. 0.10}`, the tree first differs at **0.08**, where the pan exceeds `half_world`. Re-pinned.
**The control arm caught it, not the property under test** -- as at every previous move.

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

**THE BLEACHING IS A LOSS OF TEXTURE, AND EVERY OBVIOUS STATISTIC READS IT BACKWARDS.** "The pale
regions grew" measures as median OKLab lightness **falling** (0.796 -> 0.843 -> 0.614 across the
walk) and the strictly white population falling monotonically 0.0179 -> 0.0150 — the panel gets
*darker and more saturated*. What the eye reads is those regions going **flat**. The statistic
that matches is local contrast, the 5x5 s.d. of chroma, and it collapses **3.1x at `5cc8dec` and
nowhere else** (0.01542, 0.01726, 0.01735, **0.00557**, 0.00461). **And what left was the
incoherent component**: lag-1 chroma coherence *rises* 0.62 -> 0.70 in the same step, and on the
55.5% of the frame at least 8 px from any pre-fix magenta the fall is the same 2.9x. A texture
that vanishes while what remains becomes more coherent was noise, not structure — *amplitude
cannot tell a small real signal from noise; coherence can*, now at a third site.

**`nonfin` IS THE WRONG DENOMINATOR FOR "TEN TIMES MORE CHANGED THAN WAS BROKEN".** Median energy
drift before the `dtau` fix is **2.2905 — 229% of the total energy** at the median pixel. The
frame was NaN on 2.87% and quantitatively meaningless on most of the rest, so against a **682x**
fall in median drift, 33% of labels changing is a small number. The flips are enriched 2.22x on
the pre-fix magenta pixels and 1.36x at 16 px, but that set is 2.87% of the frame — **~94% of the
flips are outside it**, which is what a correct fix looks like here.

**THE TIMELINE HAS ONE STATE FROM 08-25 TO 08-27, AND THE BIGGEST MOVE IS THE ESCAPE COMMIT.**
`f4084de` through `483b630` are **bitwise identical** over 1048576 pixels — 34 hours, the whole
08-26 scheduler run, zero pixels. `077b092 -> e53223d` then flips **39.95%** of labels against the
`dtau` fix's 33.22%, and the two escape commits take escape 0.676 -> 0.262 and bounded
0.054 -> 0.369. The screenshots being reasoned from sit at `4b26466`, **after** both, so the state
treated as original already contains the largest label move in the range. Four null controls pass
bitwise, including `4b26466 == 71de13f` (the named "true before" touches no `src/`) and
`FixedPerInterval` reproducing `71de13f` exactly — so the flag *is* a faithful reconstruction of
the pre-fix commit and those comparisons were not post-fix against post-fix.

**`hot` FALLS WITH `eta`; THE LABELS DO NOT CONVERGE.** `hot` 0.8750 -> 0.4871 and median drift
9.027e-3 -> 7.281e-7 over an 8x refinement — **12400x**, order ~4.5, so truncation and not a wrong
equation, and the 86% is the threshold sitting far below this slice's bulk. But escape fraction
moves monotonically `0.2016 -> 0.1171` with no sign of settling and `chord max` is **2.000,
antipodal, at every rung**. At horizon 50 the *shape* field converges (`chord p50` falls ~6x per
halving) and the *classification* does not.

**IT IS NOT HOW OFTEN THE REFERENCE SWITCHES, IT IS HOW OFTEN NEIGHBOURS SWITCH DIFFERENTLY.**
`config_stability` switches **21x** as often as `preset_prho` (5.717 against 0.270 over 32
boundaries), but count alone does not order the slices — `preset_plambda` switches *more* than
`preset_shape` and carries **70x fewer** sharp gradients, because its switches all fall at one
boundary and neighbours switch together. The ordering quantity is the fraction of neighbour pairs
with differing switch history: **58.4% / 13.6% / 5.1% / 1.6%** against gradient densities 0.0537 /
0.0212 / 0.0003 / 0.0006. The alignment test reads **3.43x** (against 1.69x from the `t = 0`
proxy) with a **shifted control at 0.65x** — none of it is smoothness. And the paired increment
inside one trajectory has the switch boundary above the hold boundary on **82.7-99.0%** of pixels
in all four slices. The confound is that a switch coincides with fast deformation; the argument
against is `preset_prho`, quiescent, where switches still cost four orders more than holds.
**The reference is already chosen once per SYNC BOUNDARY, never per step**, so that remedy is the
shipped behaviour.

**THE COLLISION-CADENCE CANDIDATE POINTS THE OTHER WAY.** The median physical RK4 step on this
slice is **9.6e-3** against the reference's `dtMacro = 0.002`, so **the reference samples collision
about five times more often in physical time than prin-rs does**. `total_substeps` is summed over
all `E+1` copies, and dividing by the raw sum understates the step **eightfold** — which is the
whole difference between "prin-rs tests far more often" and the truth.

**A BISECT HARNESS MUST BE REGENERATED PER COMMIT, AND AN ABSENT FIELD IS THE SIGNAL.** Defaults
are exactly what changed across a range like this, so a harness checked out with the code measures
the defaults. `results/timeline/harness/run.sh` emits the `EnsembleCfg` literal field by field from
one template and **logs the fields that do not exist yet** (22, 23, 24, 26, 27, 28, 29 across the
walk). Both colour windows are fixed constants shared by the strip: an auto-ranged ramp per panel
would stretch each commit's own p1-p99 to full scale, which on a question about bleaching
manufactures or hides the thing being measured.

**THE MASS GATE FIRED, AND WHAT IT CAUGHT WAS A DIFFERENT PHYSICAL SYSTEM WEARING THE SAME
WINDOW.** Extending the bisect to 08-24 stops for two distinct reasons and the gate separates
them: `Chart::Latent` and `decoder::Latent` **do not exist** at `961a313` or `be478e1`, so the
two leading candidates cannot have altered a chart that was not there — and `be478e1` is the
commit that *creates* the decoded-mass path, so there is no earlier state in which this slice
exists to be broken. At `30d713f` and `45e7dcb` the chart exists but decodes to
`(0.31628, 0.48444, 0.19928)` against the expected `(0.32735, 0.42763, 0.24502)`; rendering them
would have produced a plausible panel of the wrong problem. `030de1a` (08-25 13:39) is the
earliest renderable state and is **bitwise identical** to everything through `483b630`.

**AND THE MASS PATH ITSELF IS CLEAN — 8388608 SYSTEMS, AND ONLY ONE ROW OF THE TABLE COULD HAVE
FAILED.** `examples/mass_audit.rs` on every pixel and all `E+1` copies, built by
`jitter::copies_with_path` with `evaluate_at`'s own arguments: `max|sum m_i v_i| = 2.9e-17` and
`max|sum m_i r_i| = 2.1e-16` on `config_stability`. The two are asserted **separately**, because
zero momentum does not imply zero first moment and a construction that assumes a COM-centred
input returns a drifting system without one. The presets are in the table and carry nothing: on
an equal-mass slice a mass error is invisible by construction, so an equal-mass control could
not have produced this result. `m spread` across a footprint is exactly zero, which is what makes
`copies[0].m` exact on a configuration chart rather than an approximation.

**THE HISTORY IS LINEAR OVER THE INTEGRATOR, AND TWO "DIVERGENT BRANCHES" ARE ONE COMMIT.**
`git log --all --graph` over `driver.rs` and `outcome.rs` is a single straight line; no merge in
`f4084de..f7d2a31` touches `driver.rs`, `outcome.rs` or `pixel.rs`. `dtau-step-control` and
`overshoot-clamp` are **the same commit** `f7d2a31`. `closure-criterion` and
`escape-distance-gate` are ancestors of it and are already in `origin/main` (`66639b2`, PR #25).
The one non-ancestor tip, `criterion-sweep`, has a `src/` tree **identical** to `84f9cbd`. Local
`main` at `0c070e4` is a stale ref, not a fork — **check `origin/main` before reading a branch
list as a divergence.** Triple collision is likewise not a missing prototype: `State::is_triple`,
`triple_ejection` and the `>=2-pair` rule are in `src/outcome.rs` at `0114be4` and earlier.

**`R = 1` ON THE LATENT CHARTS IS AN ALGEBRAIC IDENTITY, NOT A `z0 = 0` COINCIDENCE — AND THE
UNITS QUESTION IS LIVE IN THE OPPOSITE PLACE.** `decoder::config` writes `rho~ = (cos a, 0)` and
`lam~ = sin a (cos b, sin b)` in **mass-weighted** Jacobi coordinates, so
`I = mu_rho|rho|^2 + mu_lam|lam|^2 = cos^2 a + sin^2 a` with the mass factors cancelling, and
`sum m = 1`. Measured over 8388608 systems, `config_stability` reads `max|R-1| = 4.441e-16` —
**exactly the same residual as `preset_shape`**, which has `z0 = 0`; the residual tracks whether
the chart sweeps alpha and beta (trig round-off), not whether `z0` is zero, and the two
constant-configuration presets read exactly 1.0. So `EscapeRule::Distance(12)` computes an
absolute gate of exactly 12, which is what the app's `rEsc = 12` means. Where the literals do
**not** transfer is Burrau: `escape_gate.txt` §0 has `near-field` at `R = 2.236` and
`deep interior` at `1.369`, so `r_esc = 5` means 11.18 and 6.85 absolute.

**A COMMIT MESSAGE'S COVERAGE GAP IS NOT ALWAYS THE HARNESS'S.** `e53223d`'s write-up quotes only
`deep interior`, `near-field` and `preset_plambda`, but `examples/escape_gate.rs` already carries
`config_stability` at its own settings and `results/output/escape_gate.txt` reports it: at
`r_esc = 5` the gate takes persistence at +1/+2/+4/+8 from `0.784/0.769/0.753/0.734` to
`0.968/0.958/0.944/0.923`. The gate works on this slice too. Its `r_esc` sensitivity is strong and
monotone here (escape `0.6944 -> 0.4358` across the ladder) where `near-field` is flat at 0.0000
at every rung — **read the committed output before re-running the sweep.** And the 24x24 sweep
validates against the 1024^2 panel to **0.3%** at the app's own `r_esc = 12` (escape 0.4913
against 0.4899), which is what a per-trajectory statistic is expected to do.

**`escape_confirm` IS STRUCTURALLY INERT IN THE RENDER PATH.** Read from the diff, not the
message: it holds a detection provisional only when the detection came from an **in-loop** test,
and the render path runs at `escape_every = 0`, so there are none. The strip agrees — `077b092`
moves **347 labels of 1048576 and not one shape vector** (`moved` 0, `chord p50` and `chord max`
both exactly 0).

**AN ARTEFACT CAN FAIL TO REPRODUCE AT ITS OWN COMMIT WITHOUT THE BUILD BEING
NON-DETERMINISTIC — CHECK THE FILE'S mtime AGAINST THE REFLOG BEFORE SUSPECTING THE CODE.**
`results/closure/config_stability_stop0_uniform.png` re-rendered at `220d928`, the commit that
adds it, differs on **84.12% of pixels**; re-rendered at `4b26466` it is **bitwise identical**,
both panels, every printed number to the digit. The reflog and the file mtime say why: HEAD sat
at `4b26466` from 13:52 to 21:19 on 08-27, the PNG was written at **16:38**, and `220d928`
committed it at 21:36 — **after `5cc8dec`, the `dtau` fix**. A pre-fix render committed into a
post-fix tree. The magenta fraction is 0.0029 committed against 0.0001 at `220d928`, a factor of
29, so anything using this file as a *before* was comparing pre-`dtau` against post-`dtau` and
measuring §23 over again. **Commit renders in the same commit as the code that made them, or
name the commit in the filename.**

**AND THE HARNESS DIFF IS THE FIRST THING TO CLEAR, NOT THE LAST.** `220d928` also changed
`closure_render.rs` by 13 lines in the same commit — `sub` becomes argument 4 and `tau` moves to
5. That looks like exactly the kind of argument renumbering that silently changes a run, and it
is **not** the cause: the committed invocation (`closure_render 1024 results`, from the log's own
header) passes neither argument, so both take their defaults either way. Checking it first cost
one `git show` and removed a whole candidate before any render was spent on it.

**THE IMAGE THE BLEACHING THREAD IS ABOUT IS 256x256, AND IT REPRODUCES BYTE FOR BYTE AT
`4b26466`.** `results/closure/esgate_fixed/config_stability_stop0_uniform.png` — mtime 08-27
15:20, inside the window where HEAD was `4b26466` (13:52 to 21:19) — is reproduced **bitwise, both
panels**, by `closure_render 256 <root> config_stability` at that commit. **There is no
uncommitted working-tree difference to hunt for**, and the build is deterministic across a fresh
checkout and rebuild. It was committed at `220d928` (21:36), *after* `5cc8dec` landed at 21:19, so
it is a pre-`dtau`-fix render sitting in a post-fix tree; the same command at `220d928` gives
escape 0.0433 against the committed 0.0588. Its magenta fraction is 0.0030 against the 1024^2
render's 0.0029 — same physics, quarter the linear raster — so *softness in an image is a raster
size* now has a fourth site, and a magenta cluster of a few footprints reads as a **blob** at
4x4 screen pixels per footprint.

**THE CIRCLED WEDGES ARE HIERARCHICAL ICs WITH THE HEAVIEST AND LIGHTEST BODIES FAR APART — AND
THEY ARE *NOT* NEAR AN ARGMAX DEGENERACY.** Measured on the ICs alone, which need no integration:
the magenta and dense-pale populations both order `d(0,1) < d(0,2) < d(1,2)` with **`d(1,2)`'s
10th percentile at 1.84 against the frame's 1.42** — the small-separation tail is gone — and
`tightest == (1,2)` occurs on **0.07% of the magenta against 28.9% of the frame, a 413x
depletion**. With `m = (0.32735, 0.42763, 0.24502)` that is the heaviest and lightest body as the
wide pair. Alongside: aspect 1.57 against 1.43, larger `alpha` (tighter inner pair), and **|Lz|
about 60% of the frame median**.

**The near-degeneracy hypothesis is refuted by the sign, not merely unsupported.** If these were
the pixels where AZ's `argmax` is a coin flip the tie statistics would be enriched; they are
**depleted ~2.5x** (`d[2nd]/d[longest] > 0.95` on 0.0825 against 0.2132). And `ic_class.png` — the
reference body and tightest pair over the whole frame — is a **smooth six-sector pinwheel meeting
at one point**, which is where the straight edges in the chart plane come from and which does
**not** draw the wedges. *A finding read off a wireframe is a finding about an appearance*, at the
IC layer: the straight edges suggested a discrete IC boundary, and rendering the boundary showed
it is somewhere else.

**And the class constrains without drawing.** `P(magenta | reference=0 and tightest=(0,1))` is
**0.0070** against a base rate of 0.0029 — a 2.4x lift on a class holding 22% of the frame. The
initial conditions bound where the artefact can appear; the fine structure inside them is
dynamical. **Read a hand-drawn mask's area before its statistics**: the digitised circles cover
25.4% of the frame, so that row is mostly background and only the property-selected populations
carry the signal.

**WITHDRAWN, AND THE REPLACEMENT IS BELOW — the zoom it rested on was half a frame from where it
was said to be.** ~~**THE PALE WEDGES HAVE NO INTERIOR, AND DO NOT ACQUIRE ONE UNDER 6x
MAGNIFICATION.**~~ Euclidean
distance transform on `config_stability_stop0_uniform.png`: the pale class's median inscribed
radius is **1.00 px** and **0.0%** survives a 5x5 opening, against the red band's **32.25 px** and
**93.6%**. Re-rendered over a 1/8 window at 6x finer cell width the pale class reads **1.00 px**
again and its *maximum* radius **falls** 2.24 -> 1.41, while the control shrinks the way a solid
region does under magnification. **A fixed-size geometric artefact would have grown 6x in
pixels.** The zoom render shows why: the "solid wedge with a sharp edge" is **hundreds of parallel
laminae** with smooth curved boundaries and no polygonal edge. The appearance is sub-pixel
aliasing — a mixture whose *density* changes over a few pixels, not a boundary — which is the
raster lesson one level down, at the sampling of the field rather than the size of the image.

**AND THE MASK THAT FAILED FAILED SILENTLY.** The first threshold (`C < 0.045`) selected dust:
11640 components, largest 211 px, run lengths of 2. A box-counting dimension from that reads near
1 whatever the truth is, because a scatter of isolated pixels is not a boundary — it would have
"confirmed" a geometric edge. **Component sizes and run lengths caught it; the dimension number
did not.** Check that a mask selects regions before measuring their boundaries.

**A ZOOM HARNESS FLIPPED THE VERTICAL AXIS AND RENDERED HALF A FRAME AWAY FROM WHERE IT SAID.**
`Slice::axis` runs low-to-high with index and `save_rect` writes rows in buffer order, so PNG row
0 is the **minimum** `v`: fractional `v` from the top of a saved panel maps straight to the axis
with **no flip**. `wedge_zoom.rs` wrote `1 - fv` and landed at 0.75 instead of 0.25. **Caught by
the pixels, not by re-reading the code** — cropping the panel at both candidates and comparing
against the render gives RMS 39.6 against 85.5. Draw the box on the source image and check it
before reading anything off a sub-window render.

**A DISTANCE TRANSFORM ON A PERFORATED REGION READS THE SAME AS ON DUST.** The corrected zoom puts
the wedge core at one connected component of **178486 px, 30% of the frame, 99.3% pale in a 192^2
interior window** — and its median inscribed radius is still **1.00 px** with 0% surviving a 5x5
opening, because isolated non-pale pixels riddle it. The statistic cannot tell a sponge from a
cloud; **component size can, and it was in the same output**. So *the white class has no interior*
is **withdrawn for the wedge cores** and stands only for the striated halo around them.

**THE WEDGE EDGE IS STRAIGHT AND STAYS STRAIGHT AT 6x — TEN TIMES STRAIGHTER THAN THE RIBBON
BESIDE IT.** Total-least-squares fit to the boundary: pale edge **RMS 4.28 px over 621 rows**
(max deviation 17.7) against the red band's **42.05 px over 186 rows** (max 178.4) in the same
image. The striated band dissolved into hundreds of filaments under the same magnification; the
wedge core did not. **A straight edge in the chart plane is not by itself a bug** — `ic_class.png`
shows the AZ reference-body partition is a six-sector pinwheel with straight edges — but the
observation survives its test, which the lamination reading did not.

**AND THE FOURIER TEST ON THE INTERIOR WAS VACUOUS AS RUN.** The window is 99.3% pale, so the
transform measures the sparse perforation rather than the structure. Reported rather than dropped;
it needs a window straddling the edge.

**THE RENDER HARNESSES DISABLE THE REPAIR PASS, AND `results/README.md` ASSERTS THE OPPOSITE.**
`EnsembleCfg::default()` has `refine_flagged: true` — BRIEF §2.5's remedy for the `eta` cliff,
re-integrating pixels whose `error_ratio` exceeds 10 at `eta/4`. **62 files under `examples/` set
it to `false`**, including every render harness. Measured on `config_stability`, one field changed:
`error_ratio` p99 **1.039e10 -> 35.6**, drift p99 **8.87e6 -> 1.64e-2**, drift max **1.97e12 ->
6.74e-2**, **non-finite 30109 -> 0**, escape fraction **0.0403 -> 0.0067**, with **11.1% of the
slice re-integrated**. The pale patches are `spread_shape` saturating at 0.39 because the copies
diverged to garbage. `src/bin/prin.rs` is **unaffected** — it takes the default and prints the
before/after drift.

**THE OVERRIDE WAS CORRECT WHERE IT WAS BORN AND SPREAD BY COPY.** `c03fc85` introduced it in the
same commit as the pass itself, with a rationale that still stands: experiments must characterise
the *unrepaired* kernel, and f32/f64 comparisons must not flag different pixel sets. That commit
also wrote **"the `render-*.txt` runs have it ON"** — the invariant. Over six days the line was
copied from experiment harnesses into render harnesses (`closure_render` at `71de13f`, 08-27
13:02), no commit message arguing for it, and `results/README.md:190` still asserts renders have it
on. **A convention that outlived its justification, never re-examined at the boundary it was meant
to stop at.** Same shape as `k_frac = 1.0` shipping as the default: *a configuration that silently
reproduces the old behaviour needs a guard, not a convention* — and none of these harnesses even
print the flag.

**THE MEDIAN IS BLIND TO IT; `error_ratio`'s TAIL IS NOT.** Slice-wide drift p50 moves
`4.251e-7 -> 2.560e-7`, essentially nothing, while p99 moves eight orders. A render checked on a
median would pass. That is why it survived six days and why the statistic that exists to flag
undetermined pixels is the one to read.

**AND NOT EVERY MARKED REGION IS THE BUG.** Of sixteen sampled: six read `error_ratio` p50
`9.8e5`-`4.9e7` and are the fault; ten read `1.001`-`2.0` and their pale structure is **real**.
The broken set runs at **4.7x-11.3x the nominal step rate** — working harder, not stepping coarser
— with `d_min` `3.1e-3`-`4.4e-3` against `r_coll = 5e-3`. `steps/copy` alone conflates a big step
with a short run; **the step RATE is the honest measure** and it reverses the reading.

**A STICKY BIT THAT NOTHING READS IS INDISTINGUISHABLE FROM ONE THAT NEVER FIRES.**
`AzOut::ab_floored` and `ab_min` were computed on every march since they were added and read by
**nothing** — `pixel.rs` never touched either, so they stopped one layer below `PixelOut` and no
render, dump, criterion or test could see the `T::TINY` floor fire. The floor is a genuine
*advance-anyway* site: unlike `budget_exhausted`, which is terminal, it divides by a fabricated
denominator and the run continues. Now plumbed with `dt_max` and `n_cap_hits`. The test that
earns its place asserts the floor **fires** on a constructed degenerate state and that the step is
finite — `(1e-200)^2` underflows at f64, so the doubly-degenerate hole is open at both precisions
and the floor is load-bearing rather than decorative.

**THE SATURATION HYPOTHESIS IS REFUTED IN ALL THREE FORMS THIS PORT HAS, AND ONE OF THE THREE
COULD NOT HAVE ANSWERED.** `config_stability` at 1024²-equivalent settings, 512²: `ab_floored`
**0.000000**, `budget_exhausted` **0.000000**, `n_cap_hits > 0` on **every pixel of 262144**. The
third is saturated, so its lift against `error_ratio > 10` is **exactly 1.000 by arithmetic** —
which is why the frame base rate is printed above the lift table. `capped` firing everywhere is
routine and not a fault: it fires whenever `A*B` falls below its interval-entry value, i.e. every
time bodies approach mid-interval. `principia_integrator_contract.md`'s `substep_bucket`/`N_sub`/
`N_max`/descriptor bit 5 **do not exist in this port** — that is the GLSL app's contract, and the
question has to be asked in this codebase's own terms.

**THE CLIFF IS A SLOPE, AND `error_ratio p99 = 35.6` IS THE PASS COUNT.** Over four decades of
`eta` the flagged population converges completely: median `error_ratio` **2.13e5 -> 1.000**, drift
**8.6e1 -> 3.9e-14**, p90 **3.39e9 -> 1.000**, with the control flat at 1.000 throughout. **0 of
128 fail to clear.** The `+0.00` slope at the last rung is not a floor — 1.000 is `error_ratio`'s
*converged* value, since `sigma_E(t) -> sigma_E(0)` under exact dynamics, and a statistic
normalised to 1 can never reach an arithmetic floor beneath it. **82.0% clear by rung 3**, the
shipped `refine_max_passes`, so ~18% survive — that tail *is* the p99. The pass converges and is
stopped one rung early; "the repair does not repair" is a different claim and the run was built to
be able to make either.

**A SINGLE RK4 STEP ADVANCED THE PHYSICAL CLOCK BY 2.209e128 AGAINST A SYNC INTERVAL OF 0.4, AND
THE MARCH RECORDED A CLEAN LANDING.** `dt_max` — the largest physical step as an actual `s.t`
difference across one step — reads p99 **1.263e43** and max **2.209e128** on `error_ratio > 10`
pixels against a nominal `4.0e-3`, and p99 **1.874e-2** on the rest. The acceptance path is code,
not inference: `1e128` is finite so the `is_finite` guard passes; `s.t >= dt_left - land_tol` is
satisfied by 128 orders so `landed = true`; and `t += dt_left` then corrects the **clock** to the
boundary while keeping the **state** reached at `s.t = 1e128`. **The clamp corrects the clock and
cannot un-take the step**, and `t` is clamped on both branches — so the overshoot was invisible in
every recorded quantity until `dt_max` existed. An *unbounded step with no acceptance test*, not a
cap, and curable by `eta`. The remedy shape is a step-acceptance test — reject and retry when the
taken increment exceeds its own remaining interval — which is local and needs no re-integration,
so unlike `refine_flagged` it has a live-playhead analogue.

**A SETTING COPIED FROM A CONTEXT WHERE IT WAS CORRECT INTO ONE WHERE IT IS NOT, INVISIBLE BECAUSE
NOTHING PRINTS IT.** Second instance after `k_frac = 1.0`. `refine_flagged: false` was introduced
at `c03fc85` **correctly**, in experiment and precision harnesses, in the same commit that wrote
the invariant *"the `render-*.txt` runs have it ON"*; over six days it was copied into render
harnesses one file at a time with no commit message arguing for it, reaching `closure_render` at
`71de13f`. **The failure was never the choice; it is that nothing recorded the choice.** The
architectural cause is that there was no single source of truth at all — `prin` and 111 literal
sites each constructed `EnsembleCfg` independently, so "the production config" existed nowhere.
`EnsembleCfg::production` is now the one literal, `Override` is a named value per field, and
`overrides_vs_production` **derives** the declaration by diffing so a config declares itself
however it was built — a hand-maintained list can go stale, which is the same failure one level
up. Both the diff and `Override::apply` are exhaustive with no `..` and no `_` arm, so adding a
field breaks the build until it is handled. `output::provenance_sidecar` puts it beside every
panel: the `.raw` and `.prnq` dumps have carried a full settings header since they were written —
**the PNGs were the blind spot, and so was harness stdout.**

**`refine_flagged` IS A BATCH-ONLY WORKAROUND AND HAS NO LIVE-PATH ANALOGUE.** It re-integrates a
whole trajectory from `t = 0` at finer `eta` after the fact. Under a live playhead there is
nothing to re-integrate *from* — a pixel bad at `t = 30` cannot be repaired at `t = 30.1` without
redoing thirty time units. So the reading of the propagation bug inverts: the render harnesses
were showing the **unmasked kernel** and `prin.rs`, on the default, has been hiding the same
defect behind the repair. Closing the propagation gap **closes a discrepancy, not the defect** —
and doing it first switches off the diagnostic that revealed it. The renders are still not valid
*science images*, because `spread_shape` saturates on diverged copies and the terminal class is
reclassified; `_drift.png` is where the unmasked kernel is legible and `_uniform.png` is not.
Both halves are true and the report states both.

**B WINS: A PREDICTIVE, BRANCH-FREE STEP LIMIT FIXES THE DEFECT FOR +1.9% OF THE STEPS.**
`dtau <= f*d_min/(|v_rel|_max*A*B)`, one divide from values `phys_from_state` already returns, no
trial step and no retry. On `config_stability` at 192²: `error_ratio` p99 **7.108e9 -> 1.109**,
the fraction above the flag threshold **0.1110 -> 0.0000**, and the overshoot count **634 -> 0**,
on `steps p50` 1.033e5 -> 1.053e5. `preset_shape` +0.5%, `deep interior` +30%.
`StepLimit::Predictive` at `f = 0.02` is now production; `None` stays named because every
committed number was taken under it, and `reference_opts` pins it so xcheck is 4/4.
**Read `steps`, not `secs`** — under load 85-100 the winning row timed *faster* than the baseline
while doing 1.7% more work, and `total_substeps` is the machine-independent cost.

**AND THE DUMB CONTROL DOES NOT FIX IT AT FOUR TIMES THE COST.** `Global f=0.25` — a uniform `eta`
cut, the stand-in for "widen the substep table" — leaves **153 overshoots** and `err>10 = 0.0767`
on `config_stability` at `steps p50` 4.134e5, four times the baseline. A uniform refinement buys
accuracy everywhere and still cannot bound a step whose size is set by *local* geometry. That is
the argument for a per-step limit over a global one, and it comes from the control rather than
from the candidate.

**A IS NOT VIABLE ON A GPU, AND `warps hit` IS THE NUMBER THAT SAYS SO — NOT THE DIVERGENCE
FACTOR.** Reject-and-retry reaches the same quality as B and costs +77% of the steps on CPU. The
divergence factor `mean(max per warp)/mean(per lane)` rises 1.577 -> 2.557 linear and
1.432 -> 1.975 tiled — but **the absolute level is the FIELD's, not the mode's**: step counts vary
lane to lane with no retries at all, which is why the control row reads 1.577 and quoting it alone
would condemn a mode for the field's own structure. The killer is that at the parameter A needs,
**every warp contains a retrying lane (1.0000, both dispatch shapes)** and the worst lane retried
**5.2 million** times. A also **plateaus above 1.0**: on `preset_shape`, `f` 0.1 -> 0.02 moves
`err p99` 2.087e5 -> 2.337e5, *up*, at 4x the cost, because 39 of 96 sampled pixels exhaust the
retry budget. A halving ladder bounded at 8 cannot reach where one well-chosen step goes directly.

**C WAS ALREADY SHIPPED UNDER ANOTHER NAME — `DtauMode::PerStepInterval` IS AN `A*B` GROWTH CLAMP
AT `C = 1`.** `StepLimit::AbGrowth` is **bitwise identical to the baseline** on every region at
every parameter. The brief's `dt = min(A*B, ab_entry*C)*dtau` assumes `dtau` fixed across the
interval, which is `FixedPerInterval`; under the shipped mode `dtau = eta*dt_left/(A*B)` is
recomputed per step, so `dt ~ eta*dt_left` however much `A*B` grows. Held as a two-armed test —
inert under `PerStepInterval`, active under `FixedPerInterval` — so if the older mechanism changes,
the thing that silently starts mattering announces itself. **And with B in force the cap is
redundant**: `deep interior` is identical to four digits *including the step count*,
`config_stability` moves `err p99` 1.555 -> 1.686 on 6% fewer steps. Reported, not removed —
that is a second corpus-invalidating change. `clamp_final_step` is **not** a removal candidate: it
is a correctness property, taking the convergence order 1.06 -> 2.08.

**THE CHART FAMILIES ARE TAME IN `alpha` AND NOT IN THE STEP CONTROL.** `preset_shape` was chosen
as the *clean* control and its baseline carries `err>10 = 0.0824` with **584 overshoots** at 192².
The standing "chart families sit at `alpha` 0.99-1.01 and do not exercise the criterion" is a
statement about the refinement criterion and says nothing about integration. A control has to be
clean *in the quantity being measured*, and this one was not.

**A GROUND-TRUTH COMPARISON CAN BE SATURATED AND SAY NOTHING.** Chord and label flips against the
`eta/256` reference read **1.0000 flips for a correct mode and a broken one alike** — over
`t = 50` any change of step size gives a different trajectory through a chaotic region. Reported as
saturated rather than quoted. What discriminates is `error_ratio`, because it is normalised to
exactly 1.0 under exact dynamics and therefore has an absolute meaning that a difference does not.

**THREE TESTS FAILED WHEN THE LIMIT SHIPPED, AND EVERY ONE FAILED CORRECTLY.**
`refined_pixels_are_repaired` fell over on its **own** guard — *nothing was flagged, so this test
has no subject*; `mad_based_error_ratio_cannot_separate_damaged_pixels` had no damaged population
left; and `error_ratio_minus_one_falls_with_step_size` went non-monotone because the residual is
now **round-off, not truncation** (`3.1e-9 -> 1.6e-9` across a decade). All three are pinned to
`StepLimit::None`, the kernel they are about. **That the fix deletes the subject of three
characterisation tests is the strongest corroboration in the suite**, and a pin needs the arm that
justifies it or it reads as a tolerance loosened to keep green —
`the_shipped_limit_is_already_at_the_floor_at_the_coarsest_step` asserts an eightfold step cut
moves the residual by *less* than an order.

**`eta/256` AND A FOURTH REFINEMENT PASS ARE CHARACTERISATION, NOT REMEDIES.** Both repair the
damage and neither survives a live playhead: `refine_flagged` re-integrates from `t = 0`, and a
global `eta/256` pays 256x everywhere for a local failure. Their value is diagnostic and it is
large — *`eta/256` brings every flagged pixel to `error_ratio` 1.000* is what proves this is
ordinary under-resolution rather than a wrong equation, a saturating cap or a threshold, which is
why it is the **ground truth** in the comparison and not a candidate in it.

**`Gamma*` IS DEGREE SIX JOINTLY AND AT MOST QUADRATIC IN EACH SINGLE COMPONENT, SO THE `h^2` FD
TEST HAS NO SUBJECT.** `R_i = Q_ix^2 + Q_iy^2` enters every term of Heggie's `Gamma*` at most
linearly and `W_i = L(Q_i)P_i` is linear in each argument, so the third derivative with respect to
any one coordinate is **identically zero** — and that is exactly what central differencing
truncates at. Measured: median FD error **3.3e-16 at `h = 1.0`**, a hundred-per-cent perturbation
of every component, **rising** as `h` falls at the `1/h` roundoff law (29x over 128x). A floor at
every reachable step, not a slope. AZ's version of the same test works only because its `Gamma`
carries `A B m_b m_c/|R3|`, which has third derivatives in abundance. **Asserting "no truncation"
alone would be asserting the harness produces small numbers**, which a broken harness does equally
well — so the replacement carries a control arm that adds a deliberately cubic term and requires
the same harness over the same states to recover `0.25`, which it does **at every rung**. The
transferable form: when a test ported from a sibling cannot be written, ask what property of the
new object removed its subject before weakening the threshold.

**THE `1/m` RISK HEGGIE WARNS ABOUT IS NOT LIVE ON THIS CORPUS, AND THE BOUND OVERSTATES IT BY
ORDERS.** `Gamma*` carries `1/m_1, 1/m_2, 1/m_3` and Heggie §4 says Eq. (21) is inapplicable if a
mass vanishes — a risk the Heggie port introduces that AZ never had. `Chart::Latent` saturates its
logits at `MU_MAX = 5`, which admits masses of order `e^-10`. Measured over 2,097,152 systems at
512^2, every pixel and all `E+1` copies: `config_stability` reads **`min m = 0.245`, `max 1/m =
4.08`, mass ratio 1.745**, presets exactly `1/3` and `3`. The saturation bound is a statement about
the chart family; the windows actually rendered sit nowhere near it. `examples/heggie_mass_floor.rs`
is read-only and integrates nothing.

**HEGGIE'S Eq. (18) IS VERIFIED AS THE JACOBIAN OF Eq. (17) IN FULL 3D, NOT ON THE PLANAR SLICE.**
The planar reduction `A_i = 2 L(Q_i)^T` is the one step of the transcription the paper does not
state, and checking it only where `Q_3 = Q_4 = 0` would leave the two out-of-plane columns untested
— which is what makes the reduction a claim rather than a definition. Finite-differenced in four
dimensions it agrees to **4.1e-10**, and restricted it equals `2 L(Q)^T` **exactly** while the
transposed block differs by **8.0**. Same shape as the standing `shape_pl` result: an index
assertion alone passes on a transposition, so the negative control is the test.
