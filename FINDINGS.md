# What this kernel measured

A coherent account of the work, written to be read cold by someone who has not seen the
repository. It is not a changelog and not a spec: it is what was measured, what it overturned, and
what the machinery finally settled into.

`BRIEF.md` is the authoritative spec. `CLAUDE.md` is the working agreement and carries the full
findings record in the order it was discovered, which is long and lumpy by design. This file is the
same material arranged as an argument. `RESULTS.md` and `NOTES.md` hold the per-experiment tables.

Every number below was measured in this repository, and every table is reproducible from a
committed harness under `examples/`.

---

## 1. The object, and what a pixel is

The planar three-body problem, released from rest, with the scale and orientation symmetries
quotiented out. A **slice** is a 2-plane through the space of initial conditions; a **pixel** is one
point on that plane; **one pixel is one full simulation** integrated to a playhead time `t`.

**Units.** `G = 1` and the scale symmetry is quotiented out, so `t` is in units of the system's own
crossing time `τ = √(R³/GM)` and has no duration in seconds until a mass and a length are pinned.
On the latent charts `Σm = 1` and `R = 1` by algebraic identity, so `t` *is* the crossing time
directly; Burrau in the repo's units (`M = 12`, `R = 2.2361`) has `τ = 0.965`. For scale: a
Sun-mass system an AU across puts `t = 52` at 8.3 years, a wide 3 M☉ triple at 100 AU puts it at
4780. **The kernel does not choose, and that is the symmetry being factored out** — an absolute
length would break it, which is why `r_coll` and `epsilon` are fractions of the initial hyperradius
fixed at `t = 0`.

Each pixel carries an **ensemble** of `E+1` copies, jittered inside its own cell. The copies are not
noise: their *disagreement* is the measurement. Two things are read off them.

- **`spread_shape`** — how far apart the copies' final configurations are on the shape sphere.
- **`spread_event`** — whether they agree about which pair is the tightest binary, joined with the
  terminal outcome for copies that have terminated.

The ensemble offsets are a **fixed Halton (2,3) prefix indexed by copy index** — not per-pixel, not
pseudo-random. Copy `k` sits at the same offset in every footprint at every refinement level, so a
parent shares its children's perturbation pattern and the two are compared under common random
numbers by construction. Measured: the control's per-quad noise floor falls from 0.4796 to
**0.0010** at `E+1 = 8`, and the parent/child correlation from 0.175 to **0.9998**.

That buys less than it looks like it should. Sampling noise is only ~7% of the scatter in the
refinement exponent (`var` 5.725e-1 → 5.331e-1); the other 93% is chaotic divergence. **More copies
cannot buy that off**, and the per-quad scatter (0.63) against a region separation of ~1.0 is why
the criterion resolves *regions* and not individual quads.

---

## 2. Regularisation: three implementations, one physics

Close approaches make the Newtonian `1/r²` force unbounded. Every integrator here handles that by a
coordinate transformation plus a time transformation, and the differences between them turned out
to matter more than expected — and for reasons other than the obvious one.

### 2.1 Aarseth–Zare

Ported from `reference/tb_az.py`, transcribed rather than re-derived. Two pairs are regularised
through a Kustaanheimo–Stiefel-like map about a **reference body**, chosen from the geometry, with
the time transformation `dt = A·B·dτ` where `A`, `B` are the two regularised separations.

The reference body is re-chosen **at every sync boundary**. That re-registration is the whole story
of this integrator's behaviour, in both directions:

- It makes AZ **label-insensitive**. Over all five non-identity label permutations of Burrau to
  `t = 6`, free AZ reads `3.227e-15`; with the reference frozen, `3.411e-6` — **a factor of 1.06e9**
  on a run that switches reference only 4 times in 64 boundaries.
- It makes AZ **spatially discontinuous**. Neighbouring pixels that switch reference at different
  boundaries diverge; the ordering quantity is not how often a slice switches but the fraction of
  neighbour pairs with *differing* switch history — **58.4% / 13.6% / 5.1% / 1.6%** across four
  slices, against gradient densities 0.0537 / 0.0212 / 0.0003 / 0.0006.

`xcheck` runs AZ against the NumPy reference and reads `0.000e+00` over all 32 reference points. It
is the only integrator here with an independent implementation.

### 2.2 Heggie

Heggie's global regularisation removes all three pairs at once — no reference body, so **the state
is never re-registered mid-run**. `Γ*` is degree six jointly in the coordinates, and the time
transformation is `dt = R₁R₂R₃ dτ`.

Verified by five mutation-armed anchors, not by inspection:

- Eq. (18) checked as the **Jacobian of Eq. (17) in full 3D**, not on the planar slice — the planar
  reduction `A_i = 2 L(Q_i)ᵀ` is the one step the paper does not state. Agrees to `4.1e-10`;
  restricted it equals `2 L(Q)ᵀ` exactly while the transposed block differs by **8.0**.
- The two-body radial collision run **three times, once per pair**: `d_min = 5.422e-27` and
  **39455 steps, identical for all three**. That pair-independence *is* the globality, as a number.
  AZ's version of the same test must first assert the colliding pair is one of the two regularised
  ones.

**Cost does not follow Heggie's own quoted 1.6×.** That figure is per *step* — `Γ*` has more terms.
Per `deriv` call, measured: `AZ 20.29 ns`, `HG Eq.(20)/(21) 17.82 ns` (**0.88×**), `HG
Eqs.(22)–(24) 27.40 ns` (**1.36×**), the planar reduction being why. And **up to 28% of AZ's
per-call cost is ulp fidelity, not algebra** — its `deriv` calls `r3.powf(3.0)` deliberately, to
route through the same libm call NumPy does (`6.42 ns` against `0.43 ns` for `r3*r3*r3`). Step
*counts* are equal: with the step limit off, `AZ 1.228e5, HG 1.224e5`.

### 2.3 The comparison, and the trap inside it

The first Heggie-versus-AZ tables were wrong in two ways, both instructive.

**`HgOut::drift` was the energy of the regularised state and `AzOut::drift` the energy of the
returned Cartesian state.** Different quantities, differing by **280×** on the collision fixture.
Fixed: `drift` is Cartesian on both, `drift_reg` carried beside it — because the *gap* is the
finding. Across a 32× refinement, `drift_reg` is flat at 4.4e-15 while `drift` **rises** 1.245e-12
→ 5.4e-10 with `d_min = 5.422e-27` at every rung. The integration is at its round-off floor and the
**readout** degrades: after an exact collision the Cartesian energy is `kinetic − potential`
cancelling from enormous terms.

**`stop_on_event` contaminated the diagnostic and produced five wrong conclusions in a row.** A
trajectory stopped by an event is parked *at* a close approach, where that same cancellation is
worst. Termination off → on moves Heggie's `far` drift `2.925e-12 → 7.872e-8` (**27000×**) while
its `drift_reg` sits still. What that one flag corrupted: "Heggie is 560× worse on `deep interior`"
(untruncated: **14× better**), "worse on three of five Burrau regions" (it wins **five of six**),
a spurious `near-field` tail, and `config_stability`'s magnitude understated **40×**.

**The science field and the diagnostic field want opposite settings and one run cannot serve
both.** `_uniform`/`_outcome` need termination on, because the outcome class *is* the physics;
`_drift`/`_gain` need it off.

### 2.4 The clean result

256², diagnostic pass, `refine_flagged` off, same predictive limit and step budget both sides, 32
cases — `config_stability`, five Burrau regions, and the full chart gallery:

```
              case   AZ drift p50   HG drift p50   AZ err>10  HG err>10   gain p50  frac better
  config_stability       4.209e-8      9.038e-10         423          0      +1.67        0.908
        near-field      6.446e-11      1.451e-11           0          0      +0.71        0.666
         mid-field      6.000e-10      1.203e-11           0          0      +1.70        0.992
               far      2.824e-13      2.189e-12           0          0      -0.89        0.000
        body2_core      7.316e-11      1.459e-11           0          0      +0.71        0.649
     deep_interior       1.345e-9      9.298e-10         438          0      +0.15        0.571
```

**Heggie wins 31 of 32 cases**, and `err>10` — the project's own flag for *this pixel is not data* —
runs **3915 against 74** across the gallery. The 74 are not Heggie's: every one falls in a case
where AZ reports the same count, and those are chart properties (non-finite `shape_vec` from a
genuine triple collision), not integration failures.

**The whole-frame median hides it completely.** Conditioned on AZ's own drift decile:

```
    d0 [2.9e-11,2.9e-9]   gain p50 -1.21   frac better 0.232
    d4 [2.5e-8, 4.8e-8]   gain p50 -0.06   frac better 0.487
    d9 [9.5e-6, 9.3e-4]   gain p50 +3.28   frac better 1.000
```

Two real, opposite, spatially separated effects: **1900× on the damaged decile with zero
exceptions**, and ~16× the other way where drift is already 1e-10. A median over both reads
4.240e-8 against 4.824e-8 — nothing. The statistic that answers it is the **gain map**, an image:
`log10(drift_AZ/drift_HG)` on a fixed diverging ramp, 59.9% of the frame better, 42.4% by over a
decade, in coherent arms.

`far` is the one clean AZ win — all 65536 pixels, a flat 0.7–0.9 decades. The mechanism story first
written for it ("smooth and wide, nothing to regularise") is **refuted**: `latent_shape` is tamer
still and Heggie takes it on 100% of pixels. The distinguishing feature left is **scale** — `far`
spans body positions to 13 units where the latent charts sit at `R = 1` by algebraic identity, and
`Γ*` is degree six in the coordinates where AZ's `Γ` is linear in `A` and `B`. A conditioning story,
labelled as a guess.

### 2.5 logH, and why losing does not settle anything

logH (a chartless logarithmic-Hamiltonian time transformation) was built as the falsification test:
if Heggie's advantage is *no re-registration*, a method with no reference body should inherit it.

logH loses decisively — **29× to 1494×** (RK4) and **16,000× to 215,000×** (KDK) against Heggie at
matched force evaluations. **And that does not establish the claim**, because logH differs from
Heggie in *two* ways: no re-registration *and* no coordinate transformation. The second is fatal
here. logH reaches `d_min < 1e-10` and **never** `|dE/E| < 1e-12`, its drift *rising* with
penetration depth where Heggie's `drift_reg` is flat at `4.4e-15`. A KS map removes the `1/r` from
the Hamiltonian; a time transformation only slows the clock. Every region tested is
collision-dominated, so logH was graded almost entirely on the one thing it is known not to do.

**A chartless method is not a coordinate regularisation with the coordinates left out.**

### 2.6 The default, and its actual justification

**Heggie is the default.** Not because of the scoreboard — a default justified by a win count does
not survive the next investigation — but because of the cadence measurement. Unconfounded at 256²
with step control held constant:

```
    AZ  n=250 CONTROLLED   lift 2.061   chord 5.162e-1
    HG  n=250 CONTROLLED   lift 3.824   chord 4.832e-2      <- 10.7x smaller
    HG  n=250 confounded   lift 2.327   chord 4.185e-1      <- the harness sees change
    HG  refresh_h n=125    lift 2.909   chord 5.570e-1      <- and sees BOUNDARY change
```

Removing re-registration removed the sync-cadence sensitivity. Three guards make the null mean
something: `steps p50` flat to 1.3% (so `eta` held the step size), a confounded arm at 10× (so the
harness *can* see a large chord in Heggie), and `refresh_h` at 13× (so it sees **boundary**
sensitivity specifically). **Without that last arm the null could have been the instrument.**

AZ is kept. Default does not mean only.

---

## 3. Step control: three fixes doing three different jobs

This is where the visible artefacts actually lived.

### 3.1 The `dτ` sizing

`dt = A·B·dτ`, so sizing `dτ` once per sync interval keeps the physical step at `eta·dt_left` only
while `A·B` stays near its entry value. A trajectory at a close encounter **at a boundary** has a
tiny `A₀·B₀`, so `dτ` is enormous and `dt` grows by orders as the bodies separate. Not "close
encounters" — encounters *coinciding with a boundary*, a thin set, which is why the damage clusters
spatially and a `d_min` correlation is flat.

`DtauMode::PerStepInterval` recomputes `A·B` per step with `dt_left` **held fixed**, capped at the
entry value. The cap is one-sided in the right direction: when `A·B` grows the value falls and the
blow-up goes; when `A·B` falls at a close approach the cap holds `dτ` at nominal so `dt` shrinks
with the separation — which is what regularisation buys. Cost: ~10% more steps.

A named non-candidate is kept beside it. Putting the remaining time in the numerator,
`dτ = eta·(dt_left − s.t)/(A·B)`, gives `rem_{n+1} = rem_n(1−eta)`: **Zeno by arithmetic**, the
interval approached geometrically and never completed. Measured `t/t_max` of 0.0833, 0.0303,
0.0303, **0.0080** — **and its drift is the best in the table by ten orders** (`1.3e-14`) because it
went nowhere. *Print `t/t_max` before any drift column.*

### 3.2 The landing clamp

The march exited a sync interval by **overshooting** it, and only the *clock* was corrected — the
Cartesian state written back was the overshot one. A first-order error at every boundary inside an
RK4 march. `clamp_final_step` lands the final step on the boundary, applied *after* `dtau_for_step`
returns so it composes with every mode.

The figure-eight (Chenciner–Montgomery) is the instrument: exactly periodic, so
`|state(T) − state(0)|` is a pure error with no reference trajectory and no chaos.

```
    fixed + overshoot   order 1.13        AZ                 Heggie
    perstep + overshoot order 1.06        1.06 unclamped     1.03 unclamped
    fixed + clamp       order 3.06        2.08 clamped       2.40 clamped
    perstep + clamp     order 2.08
```

Error at `eta = 1e-3` falls **827,000×**. `perstep+clamp` lands at 2 and not 3 because the clamp
sizes the last step from the *instantaneous* `A·B`, a first-order predictor, so the landing residual
is `O(h²)` — stated rather than smoothed. **The unclamped arm is the control that makes the clamped
number mean anything.**

**Neither change ships alone.** Under `FixedPerInterval` the overshoot is a fixed slice of time and
neighbouring trajectories overshoot alike — a large but spatially *smooth* error. Under
`PerStepInterval` the last step's size is a function of local `A·B`, so the overshoot becomes a
function of local state and neighbouring pixels overshoot by different amounts. Measured at 1024²:
under the clamp, switching `dtau_mode` moves the field **2.5×** less.

### 3.3 The predictive step limit — the one that removed the wedges

`config_stability` carried pale wedges with straight edges severing the ribbons, and magenta
speckle inside them. A single RK4 step was advancing the physical clock by **2.209e128** against a
sync interval of 0.4 — finite, so the `is_finite` guard passed; and `s.t >= dt_left − tol` was
satisfied by 128 orders, so the march recorded a **clean landing**. **The clamp corrects the clock
and cannot un-take the step.** An unbounded step with no acceptance test, invisible in every
recorded quantity until `dt_max` existed to record it.

The cure is a **pericentre-resolution condition**:

```
    dτ  ≤  f · d_min / ( |v_rel|_max · A · B ) ,        f = 0.02
```

One divide from values `phys_from_state` already returns. No trial step, no retry, branch-free.
On `config_stability` at 192²: `error_ratio` p99 **7.108e9 → 1.109**, fraction above the flag
threshold **0.1110 → 0.0000**, overshoot count **634 → 0**, for **+1.9% of the steps**.

Rauch & Holman (1999) show the Wisdom–Holman mapping goes artificially chaotic through overlap of
step-size resonances unless the step resolves periapse, and quote `T_p/20`; this ships at 1/50.

Three alternatives were measured and all lose:

- **Reject-and-retry** reaches the same quality at +77% of the steps on CPU and is **not viable on
  a GPU** — at the parameter it needs, *every warp contains a retrying lane* (1.0000, both dispatch
  shapes) and the worst lane retried **5.2 million** times. It also plateaus above 1.0: a halving
  ladder bounded at 8 cannot reach where one well-chosen step goes directly.
- **A uniform `eta` cut** (`Global f=0.25`) leaves **153 overshoots** at **four times** the cost. A
  uniform refinement buys accuracy everywhere and still cannot bound a step whose size is set by
  *local* geometry. That is the argument for a per-step limit, and it comes from the control.
- **An `A·B` growth clamp** is **bitwise identical to the baseline** — `DtauMode::PerStepInterval`
  already *is* one at `C = 1`.

### 3.4 Which fix did what

Ablated one at a time on `config_stability` at 512², `t_max = 50`:

```
        arm    nonfin  x none |    dense  x none | field-dense
       none      1153   1.000 |   0.0026   1.000 |     0.2058
  dtau_only        32   0.028 |   0.0026   1.000 |     0.1918
 clamp_only      1531   1.328 |   0.0024   0.923 |     0.1886
 limit_only         0   0.000 |   0.0001   0.038 |     0.1537
        all         0   0.000 |   0.0001   0.038 |     0.1524
```

**`limit_only` reproduces `all` on every column.** The `dtau` fix removes the **magenta** and leaves
the wedges untouched; the clamp alone makes non-finite *worse*. Two artefacts, never one defect —
and the clamp is not removable regardless, because it is what takes the marched order 1.06 → 2.08.

Post-fix, this slice is **integrator-insensitive at the median**: Heggie, AZ and even plain RK4 with
no regularisation at all give `chord p50` of `1.2e-6` and `4.1e-6`, invisible at 8 bits — while
*all 9216 pixels move* and `chord max` reaches 1.907 of a possible 2.0.

### 3.5 The live-playhead criterion

The target design marches trajectories forward and renders as it goes. At playhead `t` the only
things known are the initial conditions and the trajectory **up to `t`**. For any proposed remedy:

> Could this decision be made by a trajectory that has only reached `t`, without re-running it from
> `t = 0` and without knowing anything after `t`?

**Live-compatible:** the predictive step limit, `PerStepInterval`, the landing clamp, in-loop event
detection with a *bounded* confirmation lag, any per-trajectory constant drawn at `t = 0` from the
pixel index, choosing a different integrator.

**Not live-compatible, whatever the number says:** `refine_flagged` (re-integrates from `t = 0`
after seeing `error_ratio` — a pixel bad at `t = 30` cannot be repaired at `t = 30.1` without
redoing thirty time units); any rule keyed on `t_end`; any two-pass render; a global `eta` chosen
after seeing the result.

**Live-compatible and still not a fix** — the category that is easy to miss: a per-pixel dither of
the step phase would break a spatial beat by decohering it, and passes the criterion cleanly. But it
converts a *coherent* artefact into *incoherent noise of the same amplitude*. It removes the
evidence rather than the error. **A remedy that only changes the spatial correlation of an error is
cosmetic.**

---

## 4. Termination: what counts as an outcome

`d_min` is primary; `r_coll` is a recorded parameter, **not a physical constant**. The collision
fraction runs 0.0000 → 0.0242 → 1.0000 across `r_coll/R ∈ {1e-4, 1e-3, 1e-2}` while the grid's
`d_min/R` spans less than one decade. No threshold in that range is a physical event boundary.

**Escape went through three rules.** The reference's `spec > 0 && receding` is **not absorbing**:
of the 895 `deep interior` trajectories that escape under an in-loop test and not at the reference
cadence, **0.000 are still unbound one boundary later** — all transients. The shipped criterion is
`|Δn| over a window < τ AND E_rel > 0`, transcribed from `reference/escape_criterion.py`. Of the
trajectories it fires on, **1.0000 are still unbound at +1, +2, +4, +8 boundaries and at `3·t_max`**
in every region tested.

Its headline gap does not reproduce here: the reference's 383× separation between escapers and
bound trajectories is at best **6.8×** in this build, and `near-field` shows **none at all**. Three
separable causes, none decided — maturity (`|dn/dt| ~ 1/t³`; the reference quotes `t = 25–30` and
this ships at 13), population (a geometric ground truth counts a **triple dispersal** as an escape
and nothing converges to a pole in one), and `near-field` simply not having escaped yet.

Two structural facts about it are worth carrying:

- **The closure window cannot resolve inner-binary phase, anywhere.**
  `t_close = 2π√(d_min³/M)` runs **17× to 274× below the window** in every region, and a two-end
  chord cannot tell a full revolution from stationarity. So the closure arm is structurally blind to
  a tight bound pair and rejecting one rests **entirely on the energy arm**. Neither arm is
  redundant.
- **Persistence is the energy arm, not full candidacy.** Closure is a difference of neighbouring
  samples, so it jitters above `τ` on a perfectly settled escape — 0.4777 still-candidate against
  **1.0000** still-unbound. Reading the persistence test off candidacy would have scored ordinary
  jitter as a re-binding.

**Precedence is by `min(t)`, not by which arm is sampled more often.** `classify` ranked collision
above escape unconditionally and discarded both times; collision is sampled inside the RK4 loop and
escape only where the state is Cartesian, so an escape at `t = 5.0` noticed at `t = 5.28` lost to a
collision at `t = 5.1`. On the reference path this is large: `preset_plambda` had **990 of 996**
trajectories that fired both arms escape first and be labelled collision — **42.97% of the slice**.

**The escaping-body label is the lowest firing index, not the escaping body.** The reference's
`argmax` over a boolean returns index order, and on a dispersing system all three read unbound.
Transcribed, not corrected; 14 of the 40 golden rows discriminate it.

**Sampling cadence is a cost knob, not a correctness one** — once the confirmation guard exists.
Escape sampled only at sync boundaries quantises `t_end` to `n_sync` values wherever escape
terminates, and that draws visible concentric bands. But the labels are stride-invariant under the
guard while the `t_end` resolution improves 56×.

---

## 5. The refinement mechanism

This is the part the whole build was for, and the part that was wrong longest.

### 5.1 What it is for

A quadtree over the slice. Each quad holds an `N×N` grid of footprints, each footprint an ensemble
of `E+1` trajectories. The scheduler decides, per quad: **split, keep, or stop** — and under a
budget, *in what order*. The goal is an image at a given cost that is as close as possible to the
fully-resolved one.

### 5.2 The criterion was inverted, and the metric could not see it

The original policy split a quad where `alpha = log2(spread_parent/spread_child) ≥ alpha_hi` — that
is, where halving the cell had *already* halved the spread. **That is the smooth, converging
regions.** It floored or kept exactly the quads whose spread does not fall: discontinuities and
fractal mixing at every scale coarser than their filaments.

It survived because of the metric. Every `error(B)` curve scored OKLab distance under the shipping
colouring, whose lightness is auto-ranged to each region's own p1–p99 — so a smooth region's `1e-8`
residual counted as error at every depth, and breadth-first came out near-optimal **by construction
of the metric**.

Scored on the payload instead — the nominal `shape_vec` and event class against a fixed
`eps = 0.01` in `spread_shape`'s own units — the same footprints give a completely different
answer: `far` resolved at the root; `near-field` zero unresolved at **137 quads** against uniform's
5449; `deep interior` **329** against 5377.

### 5.3 What it landed on: `Policy::Tolerance`

**Split iff any footprint in the quad is unresolved. No exponent in the split test.**

```
    unresolved(f)  ⟺  spread_shape(f) > eps                     (the payload has not settled)
                   ∨  the footprint's copies disagree on event class
                   ∨  the footprint is undetermined

    decide(q):
        if every footprint of q is resolved            → Keep         (nothing to buy)
        if q is undetermined                           → Undetermined (finer eta, not finer cells)
        if the decode collapsed                        → Collapsed    (refining makes it worse)
        if camera veto fires                           → ScreenFloor  (a veto, never a trigger)
        if alpha_area(q) < alpha_lo and spread is flat → Floor        (no gain from splitting)
        else                                           → Split
```

It reaches zero unresolved on both `near-field` and `deep interior` within **1.12× of the exact
optimum**, where the alpha policy stops at the bootstrap.

### 5.4 The area floor — the no-gain test the alpha policy was trying to be

A tolerance policy alone descends forever on a *sea* — a region unresolvable at every level. The
floor is an exponent, but on the quantity the policy actually optimises:

```
    alpha_area = log2( unresolved_area(coarse) / unresolved_area(children) )
```

with edge footprints weighted by the share of their cell inside the box, so the cells tile. This
has a clean reading: **a line reads 1, a sea reads 0, a boundary of box dimension `d` reads
`2 − d`** — so `alpha_lo` is a *dimension threshold*.

Two details are load-bearing.

- **It is judged on the children once computed, never predicted**, and floors only when the
  *spread* exponent is flat too. A smooth field under a tolerance below its cell spread gains no
  area at any level — it is unresolved everywhere until the level at which it resolves everywhere at
  once — and its spread halving per level is the gain that saves it.
- **A per-parent exponent of a thin structure is off by a factor of two at a quad boundary, in both
  directions.** A footprint cell on a boundary is shared by both quads at half weight while the
  finer grid below locates it on one side at full weight. So the exponent is judged over **two
  levels**, from the grandparent's quadrant — whose cells split cleanly at the midlines — and
  declines where more than half a sibling set's structure sits on the parent's outer edges.
  Measured on a synthetic step: min 1.00, max 1.05 over 15 splits, against 0 at exactly one level
  under the one-level form.

**The default is `alpha_lo = 0.005`, and it is a "no gain at a noise margin" test rather than a
dimension threshold.** At `alpha_lo = 0.2` — a genuine dimension cut at `d = 1.8` — the floor costs
**11% of `config_stability`'s resolvable pixels** and puts its tree *above* uniform at its own
error. At 0.005 no chart is above uniform and the savings mostly survive.

**And no threshold separates a sea from a fat fractal.** At `alpha_lo = 0.001` — floors only where
children resolve under 0.14% of their structured area, i.e. saturated to the sampling — 206 boxes
still floor on `config_stability` and 2.6% of resolvable pixels sit inside them. A sponge that thins
only below the coarse sampling scale is exactly saturated seen from above. **No statistic of the
levels computed can see the levels not computed**: the floor is a bet on the depth it has declined
to buy, and the ladder's flatness from 0.001 to 0.1 is that fact as a number.

### 5.5 Telling noise from structure

An unresolved footprint is **structure** iff at least two of its eight neighbours share its class
and a nominal shape within `STRUCTURE_AGREE = 0.1` (chord/2, ~11°). A sea footprint's neighbours are
independent draws on the sphere, so two agreeing by chance is a few in ten thousand; a filament's
neighbours along it agree exactly.

Two earlier forms failed, in opposite directions. A class-conditional coherence statistic read
**0.27 against a bar of 0.3** for a one-column filament between sea and resolved basin, and
**0.3–0.5** for a sea class confined to the unresolved third of a mixed quad — missing real
structure and counting a sea as structure. One-of-four agreement let a few percent of sea through,
diluting the edge share to exactly one half.

**A quadrant-mixture stationarity stop was built, measured, and defaults to off.** It is calibrated
on white noise, and real seas are not white: on the sea chart it fires on 34 of 2064 floor leaves
and where it fires it is wrong (3585 quads reaching 0.02% unresolved with it off, against 3397 quads
leaving 2.7% with it on). It is worse on every chart where it fires. **A stop calibrated on white
noise is a stop for white noise.**

### 5.6 Stop reasons — never quote a leaf count without them

There are **four** things that stop a descent, and for most of this project's history the tables
attributed all of them to the criterion:

| stop | what it means |
|---|---|
| `ScreenFloor` | the quad's texels are below display resolution — **a veto, never a trigger**, view-relative and never cached on a quad |
| `MaxLevel` / `MaxRelDepth` | a depth cap |
| `BudgetExhausted` | the frame or run budget ran out |
| the criterion | `Keep`, `Floor`, `Undetermined`, `Collapsed` |

Measured across the corpus: the camera veto stops **≥95% of leaves on 21 of 69 dumps**; the budget
stops up to 91% at low `τ`; the spread gate stops *everything* at high `τ` (16 leaves, all `keep`).
**"The criterion refines too much" was the weaker statement** — most of the time the criterion was
not deciding at all.

A resolved or stationary quad is now decided **ahead of** the caps: reporting a resolved quad at the
cap as `MaxLevel` attributes the stop to the cap when the criterion had already decided.

### 5.7 Ranking under a budget

`k_frac` takes the top fraction of the ranked frontier each round. **It shipped disabled twice**,
the second time as the default: `k_frac = 1.0` computes the priority, sorts the queue, and refines
all of it — the ranking runs and changes nothing. All 69 dumps in the scheduler corpus are the
uniform-mode control with new columns attached. The default is now `0.25`, and
`scheduler::assert_not_uniform_in_disguise` refuses a `results/` path from the degenerate cell —
**a configuration that silently reproduces the old behaviour needs a guard, not a convention.**

What the ranking is worth, scored by `Cache::error_of` against a five-seed random band:

- **"Beats random" is the wrong bar; the baseline is breadth-first.** `Rank::Uniform` was missing
  from every table in the corpus, and against it most standing comparisons change sign. The best
  criterion measured here never beats uniform in `near-field` at any budget, while sitting well
  clear of the random band. In `far`, uniform **is** the exact optimum.
- **On a smooth field, ranking by spread *is* breadth-first** — a quad of width `w` over a gradient
  `g` has variation `~g·w`, so spread tracks cell size and argmax picks the shallowest quad; a
  constant signal ties and falls to a level-first tie-break, which is also breadth-first. So `far`
  degenerating is **correct behaviour**, not a failure. Thirteen criteria there produce the
  *identical* leaf set — one allocation by two routes, not thirteen criteria agreeing.
- **The headroom rises with structure.** Share of achievable improvement the best row forfeits at
  `B = 1535`: `far` **0.0002** (nothing varies), `near-field` **0.0336** (structure localised),
  `deep interior` **0.0999** (structure everywhere). Where nothing varies there is nothing to rank;
  where everything varies the budget must go everywhere.
- **`greedy_lookahead_1` was never a bound**, and the old name (`greedy_oracle`) cost two PRs. On
  `far` at `B = 1535` it reads **0.54760** against a random band of 0.48550–0.52047 — the worst
  strategy in the table, under a name asserting it could not be. The real ceiling is
  `Cache::dp_optimal`, the exact minimum over all tree-shaped leaf sets, and it is **cheap**: 5461
  splits in 0.01–0.03 s. Quote three roles: **floor** = random band, **reference** = greedy,
  **ceiling** = dp.
- **A split can make the image worse.** A parent's `N×N` sample grid and its children's are
  different approximation families, not a nested refinement of one — the root's own gain on `far` is
  `−3.022e-7`. So the ceiling at budget `B` is a **prefix-min**, never `f_root(S)`.

### 5.8 The knobs, and which of them are real

- **`N` and `E` fail in opposite directions.** `N` controls how well a quad knows its *area*;
  undersampling inflates between-footprint variation, so coarse `N` **over**-refines (leaf count
  falls monotonically with `N`: 106, 31, 19, 16 at `N = 4, 7, 8, 16`). `E` controls how well a
  footprint knows its *value*; undersampling deflates the within-footprint spread against `τ`, so
  low `E` **under**-refines. Never trade one against the other.
- **A fixed threshold fails on both sides.** `τ = 1e-4` sits at the 0.4th percentile of the chart
  corpus and the 4.3rd of the whole corpus, so 95.7–99.6% of quads clear it and the tree is uniform
  at max depth; where the bulk sits *below* `τ` everything keeps and the tree is uniform at depth 2.
  **A ranking cannot land above or below a distribution.** The regional spread medians span **six
  orders** (`4.26e-8`, `9.45e-5`, `9.75e-4`).
- **The treadmill does not happen; the opposite does.** The standing argument for rank was that
  spread grows with `t` everywhere. Measured on a fixed tree across `t ∈ {4..20}`: near-field's
  median leaf spread **peaks at `t = 10` then falls 81×**, ending below `τ`; `deep interior` falls
  31× monotonically. Terminal states are absorbing, so as termination saturates the copies share an
  outcome and the disagreement collapses. At large `t` a fixed `τ` fires **nowhere** and the tree
  *shrinks*. This **strengthens** the case for rank — a rise-then-collapse has no correct fixed
  value at either end.
- **2:1 balance is a rendering requirement paid for in physics.** Every forced split integrates
  `N²(E+1) = 512` trajectories, 1.03–1.76× the tree. But leaves tile the root exactly at any depth
  difference (zero gaps, zero overlaps measured), so an unbalanced tree produces a resolution
  *step*, not a hole. It becomes load-bearing only under interpolation across leaves or when quads
  are drawn as GPU geometry — where the remedy is render-side stitching, which costs no integration.
- **`balance_forced/split` is not a quantity.** The same six trees report 0.113 or 0.007 for it
  depending on the throttle — a factor of fifteen — because at `k_frac = 1.0` the criterion reaches
  those quads first and takes the credit. Spearman between the two arms is +0.771: the statistic
  disagrees with itself. Quote `quad ×` (1.03–1.76) instead — but note it is throttle-*dependent*
  under a binding frame quota, which is the regime that matters.

### 5.9 The live march

Children requested at a boundary can first be decided at the next one, because they have to catch up
to the playhead. Measured: **84–91% of the work is catch-up**. `near-field` splits nothing until
`t = 9.75`, reaches 37 quads at the horizon and the same 153-quad tree as the static descent
seventeen rounds later.

**Merging is the split rule read backwards at a later boundary.** A parent whose four children are
leaves that did not split this boundary merges back when it has become resolved or its split shows
no gain; the children become `Merged`, and `QuadTree::resident` is what a live design holds against
`quads_computed`. The no-gain memory must **expire**, and expiring it by *deleting* it is the
opposite of expiring it — `decide` reads `map_or(true, …)`, so no memory means *never merged for no
gain* and the floor stands on its own merits, flooring **more** (421 quads against 645). A separate
`no_gain_expired` flag is the correct form (741 against 645). A time-to-live was built and **loses
at every rung**, and `ttl = 1` and `ttl = 4` are bitwise identical because a lapsed memory splits,
the split still shows no gain, and the merge pass records a fresh memory.

**A capped leaf was terminal in the live frontier**, so a resolved parent could never merge it — 24
parents reading zero unresolved footprints while holding 96 `MaxLevel` children. The caps are
re-tested every boundary now.

**A live view inherited a verdict on the whole march.** `PixelOut::n_nonfinite` counts copies the
driver flagged over the march to `t_max`; projecting it into a frame at `t = 0.8` painted the future
into the past — 19 footprints magenta at every one of 16 boundaries with **zero** non-finite at any
boundary in the live series. `live_nonfinite` records the first boundary at which each copy became
unusable, monotone by construction.

---

## 6. What the images are, and what they are not

- **Draw the tree, not only the image — and never over a uniform base.** The adaptive render says
  *what is displayed*; the wireframe says *where the tree cut*. A coarse texel tells you a leaf is
  coarse; only the wire tells you whether the structure around it was subdivided *around* it or
  straight *through* it.
- **The shipped hue map is exactly 2-to-1.** `chroma·(cos h, sin h)` with `h = atan2(n₂,n₁)` is
  identically `C_MAX·(n₁,n₂)` — linear, no branch cut at all — and it discards `n₀`, so a tight
  binary with a distant third body and a wide pair with a close third render bitwise the same.
  *Before writing a continuity test, check whether the map is even injective.*
- **An auto-ranged ramp cannot tell "no signal" from "signal".** The lightness window is each
  region's own p1–p99, so a field with no dynamic range has its **noise** stretched to full scale.
  `far` reads `error(root) = 0.60` under the shipping colouring against `0.00000` under outcome —
  and its window is `(1.3e-9, 1.1e-8)`. A ratio test misses it (span ×8). The guard needs three
  arms: span, a comparison against the region's own median energy drift, and **spatial coherence**,
  because amplitude cannot tell a small real signal from noise and coherence can.
- **`ensemble_spread` carries a scale term**, so a multi-resolution render is partly a picture of
  the tree: the copies are jittered within the *cell*, and the median spread ratio per level runs
  1.19–1.62 against the 2.000 a proportional field would show. About 12% of the lightness range
  across five levels. Stated rather than corrected — normalising by cell width would change what the
  field means.
- **Softness in an image is a raster size, not a rendering fault.** Wireframe lines are integer
  pixel writes and adaptive texels are nearest-neighbour; neither can be soft in the file. This has
  cost round trips five times.

---

## 7. The measurement discipline

The findings above are mostly *corrections*, and they were found by a small number of habits that
generalise past this project.

**A test that cannot fail is indistinguishable from a test that passes.** Ask what would have to be
true for it to fire. Catches: a label-flip count of zero at `r_coll = 1e-2` where every pixel
collides anyway (saturated); a scale-invariance test at `t_max = 6` where nothing terminated, so the
invariance was the rescaling's own arithmetic; an FD test on `Γ` that a sign error shared by `Γ` and
`deriv` would have passed; and an accounting identity that telescopes for *any* ranking, random
numbers included — which was proposed as the fix.

**Never conclude "no effect" from an aggregate without the per-pixel distribution.** Four sites,
one of them inside the run written to *confirm* a prediction: rows identical to five digits while
**all 1024** pixels moved; a median of exactly 0.000 while **37.43% of a million pixels moved, worst
5.993e-1**; a whole-frame drift median of 1.14× hiding a 1900× fix on the damaged decile.

**A statistic can report maximum confidence precisely when it is least informed.** `drift max`
scatter reads 0.000 at `n ≤ 256` because small samples never draw the one bad pixel of 16384 —
stable at the wrong answer. A collapsed decode gives `ensemble_spread` exactly zero, which reads as
*perfectly resolved*. A starved footprint reads `error_ratio` exactly **1.0000**, its converged
value, because every copy stopped at the same early point and so agrees perfectly. **Ask what the
statistic would say about a system nothing is known about; if the answer is "confident", it is
wrong.**

**Insensitivity to step size means a wrong equation, never a step-size problem** — and one level up,
**a diagnostic that gets worse as resolution improves is measuring the wrong thing**. Both have
caught multiple bugs here.

**A difference can be small because both sides are right or because both are dead.** Assert each
side still resolves what it is supposed to before reading any agreement number: a curvature term on
an affine chart is identically zero at every depth; a linearised f32 sum whose samples all collapse
to `x0` agrees perfectly with a direct path that collapsed too.

**The control arm is what keeps working.** Every fixture that had to be measured once had to be
measured again when the physics moved — five times here, and *every* time it was the control arm
that noticed, not the property under test. A guard needs the arm saying it did not cut too much; a
null needs the arm saying the instrument could have seen the effect.

**The fix for a class of defect is never the instance; it is the column.** `refine_flagged` shipped
disabled in every render harness for six days because nothing *printed* it. `k_frac = 1.0` shipped
as the default. A stop-reason column did not exist, so leaf counts were quoted as criterion
decisions. The repairs are structural: `EnsembleCfg::production` as the one literal,
`overrides_vs_production` **deriving** the declaration by diffing so a config declares itself however
it was built, an absolute kernel stamp beside the diff, `Decision::ALL` as the one decision table,
and a provenance sidecar beside every panel.

---

## 8. What is open

- **`config_slice` carries a corner window convention** (`2·pan − 1 + zoom`) while
  `latent_ui_slice` centres. Whether that is a second UI convention or the same defect at an older
  site wants its own measurement.
- **The 25–26 August scheduler corpus predates every integrator fix.** The arithmetic survives — the
  accounting identity, `level` at `|ρ| = 0.993` against the DP label, nested split sets,
  tile-the-root, breadth-first as the baseline, the four stop reasons. **Every field-conditional
  number does not**: all `error(B)` curves, the criterion verdict, every sweep, the mask saturation
  fractions. Re-measuring it is a re-derivation, not a regeneration — the criterion's input field
  was re-*ordered*, not rescaled (`ρ` 0.59–0.84 across regions), and `Decision::Floor` collapses
  17 → 0.
- **`f64::max` ignores `NaN`, so `ensemble_spread = max(spread_shape, spread_event)` swallows an
  undetermined shape spread** — 11 footprints on `deep interior` invisible in the `nonfin` column.
  Not repaired, because propagating the `NaN` moves every tree and every render.
- **`error_ratio` never consults the driver's usability flag**, so a starved footprint reads its
  converged value. Recorded, not repaired: moving it moves every `error_ratio` in the corpus.
- **`far` is the only AZ win and the mechanism is a guess** — a conditioning story about `Γ*` being
  degree six in coordinates that run to 13 units. Unmeasured, and labelled as such.
