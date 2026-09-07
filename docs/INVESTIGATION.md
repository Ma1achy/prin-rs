# The wedges: a debugging record

A four-day investigation into pale, sharp-edged patches on one slice of a three-body renderer.
Written up because the **process** turned out to be worth more than the answer: eleven hypotheses
were tested, nine were refuted, three of the refutations came from the reviewer rather than the
code, and at least four measurements were confounded in ways that would have produced confident
wrong answers.

The chronology is kept in the order it happened, wrong turns included. A write-up that presents
only the surviving line is a worse document than one that shows where the effort went, because
the failure modes recur and the conclusions do not.

---

## 0. The complaint

A rendered slice of initial-condition space (`config_stability`, latent chart, `t_max = 50`,
1024²) carried pale patches with **straight, sharp edges** interrupting otherwise continuous
ribbon structure, with magenta speckle at their cores. The reported intuition was that this looked
like a bug rather than fractal structure: *chaos does not cut and resume; a threshold does.*

That intuition drove the whole investigation and was **broadly right and specifically wrong** in
almost every particular, which is the most useful thing about it.

---

## 1. What was actually found, in order

| # | Hypothesis | Verdict | Decisive number |
|---|---|---|---|
| 1 | The repair pass is disabled in every render harness | **CONFIRMED**, but a discrepancy not the defect | `error_ratio` p99 `1.04e10 -> 35.6` |
| 2 | A substep cap saturates and advances anyway | **REFUTED** — no such mechanism in this port | `ab_floored` 0.000000, `n_cap_hits` 1.000000 |
| 3 | The failure is a cliff in `eta` | **REFUTED** — it is a slope | `2.13e5 -> 1.000` over four decades, 0 of 128 stuck |
| 4 | A single step overshoots its own interval | **CONFIRMED** | one step advanced the clock by **2.209e128** against an interval of 0.4 |
| 5 | Per-step limits fix it | **CONFIRMED** | `err>10` `0.1114 -> 0.0001`, overshoot `634 -> 0`, **+1.9% steps** |
| 6 | The wedges are the max-over-copies order statistic | **REFUTED** | roughness ratio **1.018** |
| 7 | They are constraint-switching creases in `min()` | **REFUTED** | soft-min: −4% roughness at +13% cost; harmonic form *worse* |
| 8 | They are the reference-body argmax | **REFUTED four ways** | switches *depleted* (0.648); transverse; drift moves `7.5e-5` decades under hysteresis |
| 9 | Softmax over the argmax will smooth them | **REFUTED** | roughness **up 3×**, drift up 324×, cost 3 arms/boundary |
| 10 | They are a `t_end` exposure artefact | **REFUTED** | hot set survives equal exposure at lift **3.40 of 4** |
| 11 | They are early dynamical sensitivity | **REFUTED** | density vs FTLE **−0.114**; AZ drift vs FTLE **−0.082** |

And the finding that reframed everything:

> **An independent integrator on the same initial conditions produces an unrelated drift field.**
> Leapfrog drift tracks FTLE at **+0.305** — what physics looks like. AZ drift tracks FTLE at
> **−0.082** and leapfrog drift at **−0.096**.

And then, with the argmax and the LC branch both excluded, the piece that does it:

> **The sync-boundary re-registration count.** Doubling the number of Cartesian → regularised →
> Cartesian round-trips **at fixed step size** moves the drift field by **0.444 decades** at the
> median and drops the wedge set's agreement with the baseline to **2.22 of a possible 4.0**.

| arm | n_sync | eta | steps p50 | hot lift | chord p50 |
|---|---|---|---|---|---|
| baseline | 125 | 0.0100 | 1.061e5 | — | — |
| **n=250 controlled** | 250 | 0.0200 | 1.014e5 | **2.222** | **4.442e-1** |
| **n=500 controlled** | 500 | 0.0400 | 9.999e4 | **1.885** | **6.555e-1** |
| n=250 confounded | 250 | 0.0100 | 1.999e5 | 2.053 | 1.069e0 |
| LC branch off | 125 | 0.0100 | 1.061e5 | 3.961 | 2.528e-6 |

`steps p50` is flat within 6% across the controlled rows, so `eta` held the step size and this is
not a step-size result wearing a re-registration label. The confounded row doubles it, as it must.

**The effect sizes separate the two factors cleanly**, because the hysteresis experiment varies
*which* reference is chosen while holding the round-trip count, and this one varies the round-trip
count while leaving the selection rule alone:

```
  LC branch unconditioned          2.5e-6 decades
  hysteresis, switches 17.8 -> 6.7 7.5e-5 decades
  re-registration x2, fixed step   4.4e-1 decades   <- 6000x the first, 175000x the second
```

**It is not which chart is chosen. It is how often the state is passed through one.**

---

## 2. The methodological failures, which are the point

### 2.1 A saturated mask has a lift of 1 by arithmetic

The first four correlation tests used masks covering 25–99.97% of the frame. `step count differs`
had a **base rate of 0.9997** — it could not discriminate anything, and its lift of 1.000 was
arithmetic, not evidence. The habit that eventually caught this was **printing the base rate above
every lift table**, so a saturated candidate is visible before its lift is read.

The converse also bit: in one region the mask was *empty* (`n_hot == 0`), making every
mask-derived statistic take a single value. A full mask and an empty mask are the same threshold
landing on either side of a distribution.

### 2.2 Measuring the wrong feature entirely

`edge_anatomy` spent 725 seconds correlating candidates against the top decile of **|∇ drift|** —
the *edges* of the drift field. The feature under discussion was regions of high **drift
magnitude**. Those are different objects, and every ratio quoted from that run was about the wrong
one.

It was caught by a three-word objection — *"the mask is wrong"* — after several confident
paragraphs built on it. **No amount of internal rigour catches this**; only someone looking at the
picture does.

### 2.3 A metric confounded by the thing it is measuring

Straightness was scored as structure-tensor anisotropy over a 9×9 window, giving a clean-looking
result: lines through `k ≤ 4`, "folding" by `k ≈ 5`. **The metric is confounded by density.** A
window on a *lattice* of straight lines contains several lines in different directions and reads
isotropic — indistinguishable from a tangle. The apparent folding may have been densification.

The repair — per-connected-component total-least-squares fitting — then produced its own null: the
differing set is **one 8-connected blob** covering 67.8% of the frame, so components said nothing.
Only after isolating thin lines (differing pixels with >60% coherent neighbours, 1.76% of the
frame) did the measurement work: `rms/extent` **0.0795 at median k = 3**, rising to 0.213 at
median k = 8.

Two failed metrics before one that worked, on a question that looked trivial.

### 2.4 A confound worth 6000×

Comparing two coordinate charts by integrating one interval under each and differencing the
endpoints seems obviously right. It is not: the landing residual is `O(h²)` — which is exactly why
the measured convergence order is 2.08 and not 3.06 — and `A·B` differs between charts, so **the
two arms stop at different physical times.**

Measured: the displacement explained by that time mismatch alone was **4968×** the corrected
signal, and `raw/corrected` was **6234**. Without the correction the experiment would have
measured the time transformations and reported a spectacular false positive.

This one came from the reviewer, not from the code.

### 2.5 The shifted control

Every spatial correlation carries a control in which one field is displaced by half a frame. It
preserves the marginals and destroys the alignment, so **a correlation that survives the shift is
about the two distributions rather than about where they are.**

It earned its place immediately: early-line density against AZ drift reads **+0.274 straight and
−0.184 shifted**. Without the control, +0.274 alone is unreadable.

### 2.6 Membership is not density, and the eye computes density

Five separate masks asked *is this pixel on a switching line*, returning lifts of 1.0–1.8. The
observation that broke the deadlock — from looking at the pictures, not the numbers — was that the
wedges sit where the lines are **crowded**. Density against drift: **ρ = +0.31, rising with window
radius**, shifted control −0.17, late-surface control −0.05.

A statistic the eye computes natively and no mask can express. It was still wrong as a mechanism
(§2.9), but it was a real signal that five careful tests could not see.

### 2.7 Amplitude is not coherence

`branch_jump` measured that crossing a chart boundary displaces the state by ~1.25× the local step
error, against 1.07× at a non-switching boundary, and this was written up as proving the selector
harmless. **That was stronger than the measurement supports.** A chart jump is *systematic across
its surface* where step error is incoherent, and a coherent perturbation of the same size draws an
edge where an incoherent one draws noise. The conclusion happened to survive, but on different
evidence.

### 2.8 A test whose subject does not execute

Three tests failed when the per-step limit shipped, and **every one failed correctly**:
`refined_pixels_are_repaired` tripped its own `n_ref > 0` guard — *nothing was flagged, so this
test has no subject* — because the fix deleted the damaged population the test exists to
characterise. That the fix invalidated three characterisation tests was the strongest single
corroboration in the suite, stronger than any error digit.

Each was pinned to the unlimited kernel **with an added arm justifying the pin**, or the pin would
read as a tolerance loosened to stay green.

### 2.9 The measurement proposed to confirm a mechanism refuted it

The density result (§2.6) suggested the wedges were a map of early dynamical sensitivity. The
proposed confirmation was to correlate density against FTLE. Measured: **−0.114**. And AZ drift
against FTLE: **−0.082**, while leapfrog drift against FTLE is **+0.305**.

The hypothesis died by its own proposed test. This is the best possible outcome for a hypothesis
and it only happens if the test is specified *before* the result is known.

### 2.10 You cannot audit an instrument with itself

Ten of the eleven hypotheses were tested with AZ quantities against other AZ quantities. None of
them could distinguish "property of the trajectories" from "property of the integrator", because
both live inside the same instrument.

**One run of a completely different integrator settled it in twelve seconds of compute.** The
independent arm was available the whole time.

### 2.11 Work that cannot be resumed will not finish

The convergence experiment was an ~80-minute indivisible block. It was killed three times, losing
everything each time. That is an experiment-design fault, not bad luck: **work must fit inside the
shortest interruption you expect.** Checkpointing per row made a kill cost one row.

The checkpoint is keyed on the full config provenance and **refuses to resume under different
settings** — a stale checkpoint read as current is this project's most expensive recurring failure.

---

## 3. The engineering that came out of it

Independent of the wedges, the investigation produced real fixes:

- **A per-step predictive limit.** `dtau <= f·d_min/(|v_rel|_max·A·B)` — one divide, no trial step,
  no retry, no branch. `error_ratio` p99 `7.1e9 -> 1.109`, overshoot `634 -> 0`, for **+1.9% of the
  steps**. Shipped as the default.
- **A tripwire.** `n_overshoot` counts steps carrying the interval clock past twice its interval.
  `dt > dt_left` is a bug, not a condition to handle; it went undetected for six days because
  nothing asserted it.
- **Single-source configuration.** `EnsembleCfg::production()` plus a **derived** provenance diff,
  so any config declares its own departures however it was built. The original defect —
  `refine_flagged: false` propagating by copy through five commits — was not that someone chose
  wrongly but that **nothing recorded the choice**.
- **`ab_floored` and `ab_min` plumbed.** They were computed on every march and read by nothing. A
  sticky bit nothing reads is indistinguishable from one that never fires.

And three candidate step-control mechanisms were measured and **rejected**: reject-and-retry (every
warp contains a retrying lane), the `A*B` growth clamp (bitwise inert — already shipped under
another name), and a global `eta` cut (still 153 overshoots at 4× the cost).

---

## 4. What the reviewer contributed that the process could not

Three corrections changed the direction of the work:

1. **"The argmax discontinuity is bounded and does not compound"** — wrong. A finite discrepancy at
   a switching surface *is* amplified by chaos. The real reason to prefer a hard choice is that
   each branch still integrates the actual equations, whereas `Σ wₖΦₖ` is not a flow map of
   anything.
2. **"You cannot be in two charts at once"** — a fine implementation remark and poor theory. The
   obstruction is that derivatives cannot be averaged when the independent variables and canonical
   structures differ.
3. **"Compare the arms at the same physical time"** — worth a factor of 6000 (§2.4).

Plus the two observations that redirected everything: *the mask is wrong* (§2.2), and *my old
leapfrog implementation didn't have these wedges* — which turned out to be the same result the
independent-integrator control produced, arriving a day earlier and for free.

**The pattern: the model was better at executing and instrumenting tests; the reviewer was better
at noticing that the wrong thing was being tested.** Neither substitutes for the other.

---

## 5. Open

- **Why re-registration costs what it does.** The mechanism is located but not explained. Each
  round-trip is `to_cartesian` then `to_reg` with a KS square-root and a fresh energy freeze; the
  candidates are the transform's own round-off, the re-freezing of `E`, and the landing residual
  at each boundary. The `n_sync` result says the count matters and does not say which of those it
  is.
- **Heggie 1974 global regularisation** is now the *specific* cure rather than a general one:
  its construction has **no reference body to re-choose and therefore no re-registration at all** —
  three relative vectors on equal footing, one KS transform each, one symmetric time
  transformation `dτ = dt/(R₁R₂R₃)` for the whole integration. The mechanism just identified is
  exactly the thing it does not have. Heggie's own caveat stands and is on the wrong axis for
  this purpose: he reports ~1.6× per step and calls it "significantly the weaker" on smooth
  problems, but his comparisons are **accuracy per step on one trajectory**, and the defect here
  is a discontinuity across *neighbouring* initial conditions that nobody was rendering in 1974.
- **`tie_surface`'s null is not index-matched**, which is why it reads depleted where the matched
  control reads enriched. The matched comparison is the trustworthy one; the number should not be
  quoted until fixed.

---

## 6. Reproduce

Every harness prints its full configuration as a provenance line and writes a `.cfg.txt` sidecar
beside every panel. Every colour ramp in a comparison is a **fixed constant shared across arms** —
an auto-ranged ramp per panel manufactures or hides the difference it is meant to show.

```sh
cargo run --release --example saturation_mask     512 results
cargo run --release --example cliff_ladder        256 5 128 results
cargo run --release --example step_limit          sample 128 48 results
cargo run --release --example step_limit          frame  192 48 results
cargo run --release --example warp_divergence     128 results
cargo run --release --example step_limit_gallery  384 results
cargo run --release --example edge_anatomy        384 results
cargo run --release --example soft_reference      256 results
cargo run --release --example wedge_id            256 results
cargo run --release --example switch_depth        256 results
cargo run --release --example line_id             256 results
cargo run --release --example wedge_edge          256 results
cargo run --release --example drift_time          256 results
cargo run --release --example line_density        256 results
cargo run --release --example independent_check   192 results
cargo run --release --example az_machinery        256 results
```
