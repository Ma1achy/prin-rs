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
