# Notes from reading the brief and the reference

Working notes, not spec. `BRIEF.md` is authoritative; where this file disagrees with it, that
disagreement is the point and wants settling.

---

## 1. The shared-reference question — a recommendation

**BRIEF.md §3.** Should all `E+1` copies of a pixel share the nominal copy's AZ reference body?

**Recommendation: default UNSHARED, implement SHARED as the flag, and add a per-pixel
`ref_disagree` field so the question can be answered from one run rather than two.**

The reasoning, in three parts.

### The case for sharing is real

`error_ratio` and `spread_shape` are *cross-copy* statistics. If copies are integrated in different
regularisation charts, their integration errors are drawn from different distributions with
different structure. The energy spread then inflates for a reason that is not physics, and
`error_ratio` — whose whole value is that it has no threshold and no tuned constant — reads damage
that isn't there. Each trajectory's own `energy_drift` stays healthy throughout, so nothing else
flags it. That is exactly the failure mode the brief describes, and it is plausible.

### The case against sharing is also real, and is about correctness rather than diagnostics

AZ's guarantee rests on choosing the reference body **not in the longest side**, which makes the
unregularised separation `R3` the longest, so `R3 → 0` only in a genuine triple collision
(`reference/tb_az.py:145`). Force a copy onto another copy's reference and that copy may be
integrated with its reference body *inside* its longest side. Then `R3` is no longer the longest,
and an unregularised pair can close. The AZ guarantee is void for that copy — not degraded,
**void**.

So sharing trades a diagnostic artifact for a genuine integration failure. In general that is the
wrong direction.

### Why it is worse than a simple trade

The jitter is small (`0.5 ×` cell width), so copies usually agree on the reference anyway.
Disagreement concentrates where the triangle has two near-equal longest sides — where the
`argmax` in `choose_reference` is near-degenerate. Those are precisely the pixels where sharing is
most dangerous *and* where its diagnostic benefit would be largest. The two effects are not
independent; they are the same pixels. Neither default is safe there, which is why I would rather
**measure** than choose.

### The concrete suggestion

Record per pixel, per sync point, how many copies would have chosen a reference different from the
nominal copy's — accumulate to a `ref_disagree` count, alongside the existing `switches`. Then:

- Run unshared. Condition `error_ratio` on `ref_disagree`.
- If inflated `error_ratio` is uncorrelated with `ref_disagree`, the hypothesis is dead and the
  flag can default unshared permanently.
- If it correlates, the effect is real and quantified, and the shared run is worth its cost.

This answers BRIEF.md §8 experiment 2 with one run plus a confirmation, rather than two runs and a
difference of aggregates that carries no attribution. It costs one integer per pixel.

**One clarification for the flag's scope.** The brief says the reference "can change mid-run" and
raises sharing in the same breath. These are separate knobs: sharing is *across copies*, switching
is *across time*. The flag should govern sharing only. Freezing the reference across time would
break AZ outright — the whole point of re-choosing is to keep the reference out of the longest
side as the triangle deforms.

---

## 2. Three things in the brief I think are wrong or incomplete

Raised because the brief asks for it.

### 2.1 `d_min` as specified is not what the reference measures

BRIEF.md §4 defines `d_min` as "closest approach over the whole trajectory". The reference tracks
the minimum of `|R1|` and `|R2|` only (`reference/tb_az.py:201-203`) — the two *regularised* pairs.
`|R3|`, the unregularised side, never enters `d_min`.

At the moment the reference body is chosen, `|R3| ≥ max(|R1|,|R2|)`, so this looks safe. But the
reference is re-chosen only at sync points, and between them the triangle deforms freely. A close
approach on the `R3` pair inside a sync interval is invisible to `d_min`.

Whether this matters depends on `n_sync` and how fast the geometry turns over. It is at minimum a
discrepancy between the spec sentence and the validated code, and a Rust port matching the
reference to `1e-10` will inherit it. **Decide which behaviour is wanted before porting**, because
"closest approach over the whole trajectory" and "what `tb_az.py` returns" are not the same
quantity and the acceptance test cannot distinguish them.

### 2.2 The `error_ratio` acceptance test is in tension with its own specification

**Settled after PR #3, and the estimator changed with it.** §4's MAD requirement was chosen so the
statistic survives a non-finite copy, which is sound; but robustness to outliers is the opposite of
what a *detector* needs, and with 8 copies a single wild value sits above the median of eight
deviations and is arithmetically invisible. Damaged/healthy separation measured 1.06 with MAD and
59.51 with the maximum deviation from the median. `error_ratio` is now built on the maximum
deviation — NaN-safe by construction, since a non-finite copy gives an infinite deviation, which is
the correct answer where a std gives NaN. `error_ratio_mad` is dumped alongside.

The tension below is resolved by gating the **healthy** population: a genuinely damaged pixel now
reads five or six orders of magnitude above 1, as it should, so a grid-wide max is a measurement of
the worst pixel rather than a correctness criterion. Measured healthy p99 1.0228, healthy max
5.2087, bound set at 10.0. Median is 1.000000 on both populations.

**One thing not to conflate.** §4's `max` *aggregation across* footprints (Spearman +0.956 against
+0.599) and the `max` *estimator within* a footprint are independent decisions arrived at
separately. A per-pixel correlation of the two within-footprint estimators (measured −0.035 for
MAD, +0.032 for max deviation) is a different measurement at a different level and is not evidence
about the +0.956.

The original observation, kept for the record:

BRIEF.md §4 says `error_ratio` should be aggregated by `max`, treated as a **boolean flag**, and
that **its magnitude is unstable**. BRIEF.md §5 then requires it to equal `1.0000` — four decimal
places — as an acceptance test.

A max over ~10^6 pixels of a statistic whose magnitude is declared unstable is the least likely
quantity in the system to reproduce to four decimals. Either the test means the *median* or some
other robust aggregate (in which case §5 should say so), or it means max and the tolerance needs
stating as a bound (`max error_ratio < 1.01`, say) rather than an equality. As written a port can
fail this test while being correct, or pass it by aggregating the way that happens to be quiet.

### 2.3 Non-finite copies poison `d_min` silently

`reference/tb_az.py` correctly treats non-finite state as `done` so the step budget is not burned.
But `d_min` is updated with `np.minimum` *before* that check bites, and `np.minimum(x, nan)` is
`nan`. A copy that diverges leaves `d_min = NaN` for its slot.

This is consistent with "never discard a copy" — NaN *is* the measurement outcome — but only if
downstream reduction treats it as such. Any `min`/`max` over copies will propagate the NaN to the
pixel. The Rust port should decide this deliberately rather than inherit it: I would carry a
per-copy `finite` flag and let the pixel record "one copy undetermined" explicitly, which is the
same information without a value that silently contaminates every aggregate it touches.

### 2.4 `deep interior` is not a triple collision, and it does not fail

BRIEF §2.6 says the `deep interior` pixel "drives all three bodies together", is not
regularisable, "will fail however well the integrator is built", should hit the
triple-collision outcome, and cites 190 s per probe still failing.

Measured, in both implementations, on the same initial condition (body 0 at the origin,
bodies 1 and 2 at their Burrau positions — initial separations 2.236, 1.414, 3.0):

| | Rust (reference LC branch) | numpy reference |
|---|---|---|
| `d_min` | 2.2837794877e-5 | 2.2976014100e-5 |
| `\|dE/E\|` | 1.395286245e-7 | 1.393632170e-7 |
| switches | 2 | 2 |
| reaches `t = 13` | yes | yes |
| wall clock | ~0.1 s | 1.5 s |

It is an **ordinary binary encounter between bodies 0 and 2**. Sweeping `r_coll` from
`1e-4 R` up to `R` itself, pairs (0,1) and (1,2) never register at any threshold; only (0,2)
does, and it does so already at `1e-4 R`. There is no near-triple here by any definition.

The 190 s failure §2.6 records is almost certainly the **unregularised** integrator. A close
binary approach with a distant third body is precisely the case AZ regularises — the warning
predates the method that removes it. `examples/deep_interior.rs` is the probe.

The `d_min` figures differ by 0.6%, which is the sampling-offset mechanism established in
PR #2: near a close approach `d_min` is set by where a step lands relative to the crossing.

### 2.5 The escape arm never fires at `t = 13`

Over a 32x32 near-field grid at the project's horizon, `tb.classify` and the ported escape
test both return "bound" for all 1024 pixels. Burrau's escape happens later than `t_max = 13`:
the first firings appear at `t_max = 16` (2 of 1024), and by `t_max = 20` it is 109 of 1024.

Consequences worth stating rather than discovering later:

- At the project's horizon, `spread_event` and the outcome image are driven **entirely by the
  collision arm**. The one arm with a reference contributes nothing at `t = 13`.
- Escape is sampled at sync boundaries, so `t_end` for an escape has the resolution of the
  sync grid, and *which* events are seen depends on `n_sync` and `t_max`. An event at
  `t = 9.5` is caught at `t_max = 16` (boundaries 0.5 apart) and missed at `t_max = 13`
  (0.40625 apart). This is the reference's cadence, transcribed.
- The arm latches on first firing. At `t = 20`, 109 pixels fire but only 87 are still labelled
  escaping at the horizon: a body can be unbound and receding at one boundary and recaptured
  by the next.

### 2.6 `r_coll` has no plateau to sit on

The sweep BRIEF asks for, on 64x64 near-field at `t = 13` (`examples/r_coll_sweep.rs`):

| `r_coll/R` | escape | bounded | collision | triples | median `t_end` |
|---|---|---|---|---|---|
| 1e-4 | 0.0000 | 1.0000 | 0.0000 | 0.0000 | 13.0000 |
| 1e-3 | 0.0000 | 0.9758 | 0.0242 | 0.0000 | 13.0000 |
| 1e-2 | 0.0000 | 0.0000 | **1.0000** | 0.0000 | 1.8713 |

The collision fraction goes from nothing to everything across three decades, because the whole
grid's `d_min_true/R` lands inside **less than one decade** — min 5.909e-4, median 3.769e-3,
max 4.931e-3. Cross-checked from the `d_min` distribution of an unterminated run: 0.0000 /
0.0281 / 1.0000, which bounds the nominal-copy fractions from the expected side.

So no `r_coll` in this range is a physical event threshold; it is a readout of the `d_min`
distribution. The honest position is that **`d_min_true` is the primary quantity and the
collision label is derived from it**, with `r_coll` recorded in the dump header as the
parameter it is. The default of `1e-3` is picked because it separates the tail from the bulk
on this slice, not because anything physical happens there.

`triples` is 0.0000 at every threshold: no pixel on this slice ever has two pairs below
`r_coll` simultaneously, so the >=2-pair rule is exercised only by construction in
`tests/outcome_encoding.rs`, never by this data.

**And the branch cut does not reach it.** Outcome labels were compared between `lc_stable` and
the reference branch at all three thresholds: **0 flips of 4096** each time. Unlike
`spread_shape`, which the unstable branch destroys at f32, the outcome encoding is insensitive
to it here — measured, not assumed. That is a bounded result: at f64, on this slice, at these
thresholds. Step 6 repeats it at f32, where the conditioning error is six orders larger.

### 2.7 §2.4 lists six states and gives conditions for five

`state` is 3 bits over `{escape, bounded, collision, running, sim_failed, decode_failed}`, but
§2.4's table gives `running` as "still bound at `t_max`" and gives `bounded` no condition at
all. Read literally, one of the six is unreachable.

Resolved in `src/outcome.rs` by giving each a distinct job: `Bounded` is reaching `t_max` with
nothing having fired, `Running` is *not* reaching it — the step budget ran out, so the
trajectory is genuinely still running and its final state is not a terminal answer. That keeps
all six reachable and keeps "we integrated to the horizon" distinguishable from "we stopped
early", which the dump otherwise has nowhere to record. Flagged rather than assumed: if the
literal reading was intended, it is one line.

### 2.8 `spread_event` was defined over the wrong quantity

**Corrected in Step 6's PR, at the user's direction — the error was in the brief, not the
physics.** BRIEF §4 defined `spread_event` over the modal terminal `(state, detail)`. The
refinement work it derives from had concluded on something different: the **event class** is
which pair is *currently the tightest binary*, evaluated at every sync boundary. Terminal
outcome had been explicitly rejected as a contributor, because it is terminal-grain and
inverts under lockstep — early in the march nothing has terminated, so every copy agrees, so
the statistic reports maximum confidence at exactly the playhead where least is known.

`spread_event` is now disagreement over the event class at the playhead, joined with the
terminal class for copies that have terminated, normalisation unchanged. The terminal
`(state, detail)` encoding stays as the **outcome**, for classification and rendering; BRIEF
§2.4 is correct, it is simply not the spread contributor.

Measured, near-field 32x32, `E+1 = 8`, `eta = 0.01`, f64 — nonzero pixels of 1024:

| `t_max` | event class | terminal |
|---|---|---|
| 4 | 0 | 0 |
| 8 | **110** | **0** |
| 13 | 35 | 22 |
| 20 | 352 | 323 |

At `t = 13` the two are strictly nested: every pixel the terminal statistic flags, the event
class flags too, and 143 more — 165 against 22, a factor of 7.5, with **zero** pixels flagged
by the terminal one alone.

**But the gain is coverage, not lead time, and the "~4 time units earlier" framing does not
reproduce here.** On the 22 pixels both flag, the lead time is exactly zero: they fire at the
same boundary. A collision *is* the tightest pair reaching `r_coll`, and that usually settles
the tightest-pair identity at the same boundary it terminates on. What does reproduce is
horizon-independence — at `t_max = 8` it is 110 against 0, because the terminal statistic
cannot fire before something terminates and nothing has.

**A second thing worth knowing: the playhead value is a snapshot and can *un*-fire.** The
tightest-pair identity fluctuates, so copies that disagreed at one boundary can agree again at
the next — which is why the `t_max = 8` row has *more* nonzero pixels than the `t_max = 13`
row. Non-monotone in the horizon is not what a confidence flag should be, so
`spread_event_max`, the running max over boundaries, is dumped alongside; it flags the same
165 pixels at `t = 13` and never un-fires. Which one `ensemble_spread` should use is a
judgement and the spec one is the default.

`t_spread_event` records the first boundary at which the copies disagree — **NaN** when they
never do, not `t_max`, which would be indistinguishable from disagreeing at the last boundary.
Over the 165 pixels that fire: min 2.4375, median 8.9375, max 9.7500.

`ensemble_spread` at `t = 13`: median 0.001910, max 0.571429, with `spread_event` setting it
on 35 of 1024 pixels (22 under the old definition). At this horizon it is still mostly
`spread_shape`, but no longer only that.

### 2.9 `d_min_true` is primary; `r_coll` is a recorded parameter

Adopted as spec after §2.6's plateau finding. The collision label is **derived** from
`d_min_true`, not the other way round; `r_coll` is a parameter carried in every output header,
not a physical constant. The default stays `1e-3` with the honest justification: it separates
tail from bulk *on this slice*, and nothing more.

---

## 2b. The Levi-Civita branch cut

Found while porting AZ; written up in full in [`docs/lc-branch-cut.md`](docs/lc-branch-cut.md).

The original inverse LC map computes `u0 = sqrt((|rho| + rho.x)/2)` first and derives `u1`
from it. That sum cancels catastrophically when `rho` points along negative x, and the
division amplifies it. **The Burrau default sits exactly on the cut**: bodies 1 and 2 start at
the same `y`, so their separation is `(3, 0)`, and with reference body 2 that registers at
exactly 180 degrees before anything moves.

Worst case over 3600 orientations: 6.206e-11 unstable, 4.108e-16 stable at f64; **2.2e-2** at
f32 unstable against 5.96e-8 stable.

The defect is correctness, not precision: the cut is fixed in the coordinate frame, so
accuracy depends on the absolute orientation of a configuration. The physics is rotationally
invariant; the unstable implementation is not.

It is the leading candidate for the unresolved f32 dispute — orientation-dependence means
copies of one pixel can straddle the cut differently, so a cross-copy spread partly measures
registration error rather than dynamics, which is the exact shape of "drift looks fine but
the ensemble diagnostic breaks early".

Conditioning both sides also retired a negative result from PR #2: the `t=13` cross-check
went from 1.930e-10 to **2.718e-13**, so BRIEF §9's `~1e-10` stands and the amendment
proposed there was premature.

**And a limit of the instrument that produced that wrong reading.** The divergence-vs-horizon
table distinguishes *wrong algebra* (wrong intercept, or growth that is not exponential) from
*amplified ulp noise* (exponential growth from an `O(1e-16)` intercept). It does **not**
distinguish amplified ulp noise from a small error injected **repeatedly**. Branch-cut error at
each of the 32 registrations produces the same signature: small intercept, exponential envelope.
I read the curve as conclusive and it was under-determined. Keep the table — it is still the right
diagnostic for the first distinction — but do not treat a healthy-looking curve as ruling out a
per-step or per-sync error source.

---

## 3a. Standing rules, each promoted after it caught something

### Never conclude "no effect" from an aggregate without the per-pixel distribution

Three instances, two of them in a single PR:

- The f64 `spread_shape` rows were **identical to five printed digits** between LC branches.
  Checked per pixel: the branch changes **all 1024** pixels' drift, `spread_shape` and
  `error_ratio`, with a worst per-pixel change of 6.7%. The distribution does not move; the
  pixels do.
- Shared references move the `spread_shape` **median by 1%**. Per pixel: 268 of 1024 change,
  the worst by **1.86x** — an individual pixel nearly tripling.
- NOTES §1 anticipated exactly this, which is why `ref_disagree` is dumped per pixel rather
  than compared as a difference of aggregates.

Both would have been reported as "no effect" from the summary line alone. An aggregate can
only ever say the distribution did not move; it cannot say the pixels did not. Dump the
per-pixel comparison before writing "inert".

### A test that cannot fail is indistinguishable from a test that passes

Three instances:

- The `r_coll = 1e-2` label-flip count came back **zero** — but at that threshold every pixel
  collides regardless, so the label is saturated and *cannot* flip. Zero there is not
  reassurance, and the 152 flips at `1e-3` are the real answer.
- The scale-invariance test at `t_max = 6` passed while measuring nothing: no pixel terminated,
  every `t_end` was exactly the horizon, and the invariance was the rescaling's own arithmetic.
  It now asserts that some pixel terminates early.
- The finite-difference test on `Gamma` would have passed a sign error present in **both**
  `Gamma` and `deriv`. That is why the chain anchors on
  `gamma(s,E) == A*B*(energy_phys(s) - E)` first, and why a deliberately sign-flipped variant
  is kept as a `#[should_panic]` case.

The question that catches all three is the same: **what would have to be true for this test to
fire?** If the answer is "nothing in this configuration", the test is decoration.

### A statistic can report maximum confidence precisely when it is least informed

Not "a statistic can be noisy" — the inversion is the point. These fail *towards* confidence,
and they do it hardest exactly where the answer matters most.

- **`drift max` scatter reads 0.000 at `n <= 256`** and 1.000 at `n = 4096`. The apparently most
  stable quantity in the whole convergence table is stable because the tail is one pixel of
  16384 and small samples never draw it. It is stable at the wrong answer, and its stability
  *is* the symptom.
- **Outcome purity reads pure under lockstep.** A `spread_event` built on the terminal
  `(state, detail)` reports every copy in agreement early in the march — because nothing has
  terminated yet. Maximum confidence at the playhead where least is known. Measured: 0 of 1024
  pixels firing at `t_max = 8` against 110 for the event class.

Same inversion, different costumes: a quantity that has not yet had the chance to disagree, read
as a quantity that agrees. The check is to ask what the statistic would look like on a system
about which nothing is known — if the answer is "confident", the statistic is wrong.

### `n_sync` fixed while `t_max` varies compares different discretisations

`dtau = eta * dt_left / (A0*B0)` and `dt_left` comes from the sync grid, so changing `t_max` at
fixed `n_sync` changes the step size, and the runs are not the same trajectory sampled at
different playheads. Two instances:

- **The escape arm.** An event at `t = 9.5` is caught at `t_max = 16`, where the boundaries are
  0.5 apart, and missed at `t_max = 13`, where they are 0.40625 apart.
- **The `t_max` sweep in PR #5**, which I offered as evidence that `spread_event` un-fires. The
  unguarded running max across it reads 109, 297, 110, 488, 165 at `t = 4, 6, 8, 10, 13` — a
  running max cannot fall, which proves the rows are not nested. The un-firing is real; that
  evidence for it was not.

**Any sweep over `t_max` must scale `n_sync` with it**, or the horizon and the discretisation
move together and neither can be attributed. To vary the playhead alone, run once and evaluate at
each boundary — which is what `examples/latching_decision.rs` does.

### All five of these rules were found the same way

The two above join:

- never conclude "no effect" from an aggregate without the per-pixel distribution;
- a test that cannot fail is indistinguishable from a test that passes;
- and, below, the provenance of a number that could not be reproduced.

Every one of them came from the same question: **what would have to be true for this measurement
to see the thing it is aimed at?** Not "is the number right" — "could this number have come out
differently". An aggregate that cannot move, a test that cannot fire, a statistic that cannot
disagree, a sweep whose rows are not comparable, a lead time whose baseline was a different
criterion. In each case the arithmetic was correct and the measurement was still empty.

### Where the "~4 time units earlier" figure came from

Recorded because it could not be reproduced here and the reason matters more than the number.

The figure was measured against an **escape-based** terminal criterion, where the terminal
label genuinely lagged the event class. Here the terminal arm is *collision at `r_coll`*, and a
collision **is** the tightest pair reaching threshold — so both fire at the same boundary by
construction, and the measured lead is exactly zero on all 22 pixels where both fire.

The lead time was never a property of the event class. It was a property of what it was being
compared against. The justification that holds independently of the comparison is **coverage
and horizon-independence**: 165 pixels against 22 at `t = 13`, strictly nested with none
flagged by the terminal statistic alone, and 110 against 0 at `t_max = 8`. That is a stronger
claim than the lead time ever was.

## 3. Smaller observations

- `reference/tb_az.py` uses `eta=0.02` as its default; `tb_all_az` and the smoke test pass
  `eta=0.01`, which is what BRIEF.md §2.3 specifies. The bare default is the odd one out.
- The reference's `dtau` sizing (`eta*dt_left/(A0*B0)`) fixes the *first* physical step as a
  fraction of the remaining interval, so effective resolution depends on `n_sync` through
  `dt_left`. Worth pinning `n_sync` explicitly in the cross-check test rather than leaving it at a
  default, or the `1e-10` agreement target is comparing two different discretisations.

---

## 4. Step 6 — is Aarseth–Zare usable at f32?

Near-field 32x32, `t = 13`, `E+1 = 8`, `eta = 0.01`, `r_coll = 1e-3 R`, conditioned LC branch
unless stated. Initial conditions are generated once in f64 and cast down, so no row differs by
its initial condition.

### The four combinations

| prec | shared | drift med | drift max | er med | er max | ens_sp med | ref_disagree |
|---|---|---|---|---|---|---|---|
| f64 | off | 2.755e-9 | 3.070e-1 | 1.0000 | 2.6103e3 | 1.9095e-3 | 1174 |
| f64 | on | 2.755e-9 | 3.070e-1 | 1.0000 | 2.6103e3 | 1.8884e-3 | 0 |
| f32 | off | 9.293e-6 | 1.676e1 | 1.0039 | 1.1978e5 | 1.8875e-3 | 1152 |
| f32 | on | 9.307e-6 | 2.613e1 | 1.0028 | 1.8671e5 | 1.8777e-3 | 0 |

The drift tail, which the medians hide:

| prec | shared | p50 | p90 | p99 | max | `>1e-3` | `>1` |
|---|---|---|---|---|---|---|---|
| f64 | off | 2.755e-9 | 4.319e-8 | 9.257e-3 | 3.070e-1 | 22 | 0 |
| f64 | on | 2.755e-9 | 4.314e-8 | 1.108e-3 | 3.070e-1 | 11 | 0 |
| f32 | off | 9.293e-6 | 1.937e-5 | 1.740e-2 | 1.676e1 | 51 | **2** |
| f32 | on | 9.307e-6 | 1.939e-5 | 1.187e-2 | 2.613e1 | 42 | **2** |

**Answer: yes, with a named caveat.** The f32 median drift is 9.3e-6 — three and a half orders
worse than f64, and exactly what `eps ~ 1.19e-7` compounded over ~5000 RK4 steps predicts, not a
sign of a broken port. Outcome labels agree with f64 on 1022 of 1024 pixels. The caveat is the
tail: 2 pixels of 1024 lose more than the total energy of the system. Those pixels are not
data. `error_ratio` flags them, which is what it is for, and "never discard a copy" means they
are reported rather than removed.

### 1. The conditioned branch fixes `spread_shape` at f32

| prec | LC branch | median | p99 | max |
|---|---|---|---|---|
| f64 | conditioned | 1.9095e-3 | 2.1457e-1 | 2.8315e-1 |
| f64 | reference | 1.9095e-3 | 2.1457e-1 | 2.8315e-1 |
| f32 | conditioned | 1.8875e-3 | 2.2874e-1 | 3.3302e-1 |
| f32 | **reference** | **6.1582e-2** | 2.8454e-1 | 3.6668e-1 |

The reference branch at f32 inflates the median by **32x** against the f64 truth of 1.9095e-3.
The conditioned branch tracks it to **1.2%** at the median, 6.6% high at p99, 17.6% high at max.
That settles the prior dispute: the symptom was the branch cut, and conditioning removes it.

Two honest qualifications. The f64 rows are identical to five digits, which is a real result
and not a flag failing to reach the kernel — checked per pixel in `examples/flag_effect.rs`,
where the branch changes **all 1024** pixels' drift, `spread_shape` and `error_ratio`, with a
worst per-pixel `spread_shape` change of 6.7%. The distribution does not move; the pixels do.
And this run produced **no** NaN pixels at f32 on either branch, where PR #3 saw them — so the
NaN observation is configuration-dependent and should not be quoted as a general property.

### 2. The shared-reference flag does not help, and hurts the tail

Ratios are shared/unshared, so above 1 means sharing made it worse:

| prec | drift med | drift max | er max | `spread_shape` med |
|---|---|---|---|---|
| f64 | x1.000 | x1.000 | x1.000 | x0.9889 |
| f32 | x1.001 | **x1.559** | **x1.559** | x0.9948 |

Sharing lowers `spread_shape` by about 1%, in the direction the original hypothesis predicted
— but by 1%, against the 18.8x that motivated the concern. At f32 it makes the worst pixel
**56% worse**. **Recommendation: keep the default unshared**, consistent with NOTES §1.
Demoted, not eliminated: the flag stays and both settings stay measurable.

**And the aggregate hides something.** Sharing changes 152 of 1024 pixels' drift and 268 of
1024 pixels' `spread_shape`, with a worst per-pixel `spread_shape` change of **1.86x** — an
individual pixel nearly tripling while the median moves 1%. That is precisely the failure mode
NOTES §1 anticipated, and why `ref_disagree` is dumped per pixel rather than compared as a
difference of aggregates.

### 3. The branch cut DOES reach the outcome encoding at f32

| prec | `r_coll/R` | label flips (of 1024) |
|---|---|---|
| f64 | 1e-4 | 0 |
| f64 | 1e-3 | 0 |
| f64 | 1e-2 | 0 |
| f32 | 1e-4 | **37** |
| f32 | 1e-3 | **152** |
| f32 | 1e-2 | 0 |

At f64 it was inert (§2.6). At f32 the reference branch flips **152 of 1024** outcome labels at
the default `r_coll` — 14.8% of the grid classified differently by a registration artefact.
This is the failure mode a continuous-field check cannot see: a discrete label flips rather
than a number shifting.

The zero at `1e-2` is not reassurance. At that threshold every pixel collides regardless
(NOTES §2.6), so the label is saturated and cannot flip. The flips concentrate exactly where
the classification is actually deciding something.

### The floor divergence, stated up front

| constant | f64 | naively cast to f32 | f32 in use |
|---|---|---|---|
| `TINY` | 1e-300 | **0** | 1e-37 |
| `SYNC_EPS` | 1e-15 | 1e-15 | 1e-6 |
| `DRIFT_FLOOR` | 1e-30 | 1e-30 | 1e-30 |
| `DIST_FLOOR` | 1e-12 | 1e-12 | 1e-12 |

`1e-300` is exactly zero at f32, so the reference's guard would stop guarding. `ulp(13)` at f32
is 9.537e-7, which is 9.5e8 times the reference's `1e-15` slack, so `t < t_target - 1e-15`
degenerates to `t < t_target`.

One thing the floor does *not* cover: `TINY * TINY` underflows to zero at f32, and `A*B` is a
product of two floored quantities, so a doubly-degenerate state gives `dtau = eta*dt_left/(A*B)
= inf` rather than a large finite step. That is caught by the explicit `is_finite` test in the
RK4 loop, not by the floor. Worth knowing which guard is doing the work.

Gate (b) at both precisions, tolerances set by the type:

```
  f64: d_min 1.2881e-11 (tol 1e-10)   |dE/E| 2.9473e-14 (tol 1e-12)   steps 24839
  f32: d_min 1.2218e-11 (tol 1e-9)    |dE/E| 2.8553e-6  (tol 1e-4)    steps 24839
```

The `d_min` half survives the cast outright — f32 meets BRIEF §5's `1e-10` as written, because
regularisation still carries the trajectory through collision. The energy half cannot: `1e-12`
is five orders below f32 eps.

### A bug the shared-reference path found

Since Step 5b the nominal copy can terminate early, so its `refs` record is shorter than
`n_sync` and the shared policy indexed off the end — an outright panic the first time the two
features met. Fixed by falling back to the per-copy choice past the nominal's last boundary:
sharing applies where the nominal has a choice to share. Regression test in
`tests/f32_precision.rs`.

---

## 5. Should `spread_event` latch? — measured, not judged

The numpy work observed ensemble spread **falling 6x between `t=6` and `t=8`**
(diverge-then-reconverge), which is why the divergence accumulator latches. The same shape
appears here. But a **discrete** label has a failure mode a continuous divergence measure does
not: if two pairs are near-equal in separation, copies can disagree about which is *tightest*
without their trajectories having diverged at all, and a running max would latch that
permanently.

So: at the boundary where the copies first disagree, how close were the two tightest pairs?
`tie_ratio_at_disagree` is the second-tightest over tightest separation, minimised over the 8
copies. 1.0 is an exact tie. Near-field 32x32, `t = 13`, f64:

| population | min | p10 | median | p90 | max | n |
|---|---|---|---|---|---|---|
| all that ever disagree | 1.0000 | 1.0006 | 1.0040 | 1.0884 | 2.3587 | 165 |
| **un-fired** (disagreed, then re-agreed) | 1.0000 | 1.0006 | **1.0030** | 1.0193 | 1.1098 | 130 |
| still disagreeing at the playhead | 1.0004 | 1.0066 | 1.0797 | 1.1636 | 2.3587 | 35 |

**129 of the 130 un-fired pixels were at a near-tie** (below 1.1), median 1.0030 — essentially
exact ties. An unguarded running max would light 165 of 1024 pixels where 35 have genuinely
diverged: **79% of the firing pixels lit permanently for a labelling artefact.**

**But the tie ratio cannot be the guard.** The two populations are shifted, not separated —
genuine disagreements also sit near 1 (median 1.0797, p10 1.0066). A threshold at 1.1 would
admit 14 pixels and drop most of the genuine ones. Reported rather than fitted.

**Persistence separates them, and cleanly:**

| population | `n_disagree` median | max | longest run median | max |
|---|---|---|---|---|
| un-fired | 1.0 | 3.0 | **1.0** | **2.0** |
| still disagreeing | 10.0 | 10.0 | **10.0** | **10.0** |

| latch requires a run of | genuine kept | artefacts admitted |
|---|---|---|
| >= 1 | 35/35 | **130/130** |
| >= 2 | 22/35 | 1/130 |
| >= 3 | 22/35 | **0/130** |
| >= 4 | 22/35 | 0/130 |

`LATCH_RUN = 3`, chosen from that table rather than by eye. The 13 genuine pixels a run guard
alone misses are ones that began disagreeing within three boundaries of the horizon — censoring,
not a false negative — so the field joins the guarded latch with the playhead value.
`spread_event_latched` is then lit on **35/35 genuine and 0/130 artefact**.

### The latch is a no-op on this slice, and that is the result

Evaluated at every boundary of one `t = 13`, `n_sync = 32` run — nonzero pixels of 1024:

| k | t | playhead | latched | unguarded max |
|---|---|---|---|---|
| 3 | 1.6250 | 0 | 0 | 0 |
| 7 | 3.2500 | 0 | 0 | 26 |
| 15 | 6.5000 | 0 | 0 | 26 |
| 23 | 9.7500 | 22 | 22 | 165 |
| 31 | 13.0000 | 35 | 35 | 165 |

The guarded latch tracks the playhead value exactly. Every genuine disagreement on this slice
persists to the horizon, so latching adds nothing here — it is cheap insurance for regions where
one does re-agree, and it costs nothing to carry. The unguarded version over-reports by 4.7x.

**Recommendation: keep `spread_event` (playhead) as the spec field and as what
`ensemble_spread` uses. Do not adopt `spread_event_max`.** All three stay dumped.

### A correction to my own evidence in PR #5

I cited the `t_max` sweep — 110 nonzero pixels at `t_max = 8` against 35 at `t_max = 13` — as
evidence that the playhead value un-fires. **That comparison is invalid.** `n_sync` is fixed
while `t_max` varies, so the sync grid changes with the horizon, `dtau` changes with it, and the
rows are different discretisations rather than one run truncated at different playheads. The
unguarded running max across that sweep reads 109, 297, 110, 488, 165 at `t = 4, 6, 8, 10, 13` —
a running max cannot fall, which proves the rows are not nested.

The un-firing is real; the evidence I gave for it was not. The correct demonstration is within a
single run: 130 of the 165 pixels that ever disagree have re-agreed by the horizon. Same
cadence-dependence as the escape arm in PR #4, in a new place.

---

## 6. Experiment A — the refinement criterion, tested without a scheduler

BRIEF §8 experiment 1. The criterion compares a parent quad against its children, and a fine
uniform grid already *contains* every coarser scale by aggregation: pool a 2x2 block's copies to
synthesise the parent, compare against the children, and the whole exponent machinery is
testable with no quadtree. The absence of one is a feature of the test — nothing here can be an
artefact of a scheduler, because there is no scheduler.

`alpha = log2(spread_parent / spread_child)`, child taken as the median over the four.

### The control does not return 1.0, and that is the result

`alpha` for `sigma_E(0)` has true value **exactly 1.0**: `sigma_E(0)` is proportional to the
jitter and therefore to the cell width, so doubling the cell doubles it. Measured, near-field
64x64 -> 32x32 parents:

| E+1 | estimator | p10 | median | p90 | **p90−p10** | subsampled median |
|---|---|---|---|---|---|---|
| 8 | rms | 0.8525 | 1.0762 | 1.3321 | **0.4796** | 1.0191 |
| 8 | max_dev | 0.8787 | 1.1248 | 1.3855 | 0.5067 | 0.9865 |
| 16 | rms | 0.8614 | 1.0137 | 1.1817 | 0.3203 | 0.9999 |
| 16 | max_dev | 0.8613 | 1.0391 | 1.2115 | 0.3502 | 0.9730 |
| 32 | rms | 0.8902 | 0.9887 | 1.1025 | 0.2123 | 0.9884 |
| 64 | rms | 0.9169 | 0.9862 | 1.0663 | 0.1493 | 0.9949 |

**At the project's `E+1 = 8` the per-quad noise floor is an interdecile width of 0.48 in
`alpha`** — a factor of 1.4 in the ratio — on a quantity whose true value is exactly 1. It falls
as `1/sqrt(E)`: 0.48, 0.32, 0.21, 0.15 against ratios of 0.667, 0.663, 0.703 versus the
predicted 0.707.

There is also a **+7.6% median bias at `E+1 = 8`**, and the subsampled column identifies its
cause: drawing `E+1` of the parent's pooled `4(E+1)` copies puts the same sample count on both
sides and the median goes to 1.0191, and to 0.9999 at `E+1 = 16`. A parent pools four times as
many copies as a child, and a spread estimator's expectation depends on sample size. **Match the
counts, or the exponent is biased before any physics enters.** This is why the experiment uses
an rms deviation rather than `error_ratio`'s max deviation: an order statistic's bias with
sample size is much larger, as the `max_dev` rows show.

### The criterion discriminates regions, not individual quads

Fine grid 64x64, `t = 13`, f64. `alpha` for `spread_shape`:

| region | min | p10 | median | p90 | max |
|---|---|---|---|---|---|
| near-field | −0.6678 | −0.0989 | **0.1722** | 0.5324 | 10.7910 |
| body2 core | −0.5254 | 0.0551 | **0.3390** | 0.7312 | 10.7508 |
| mid-field | 0.7427 | 0.9546 | **1.1781** | 1.3930 | 1.7517 |
| far | 0.6000 | 0.9276 | **1.1716** | 1.4161 | 1.8814 |

**This is the criterion working.** In the tame regions `alpha ~ 1.17`: the shape spread scales
with cell width, as a smooth field must, so refining halves it and refinement pays. In the
chaotic regions `alpha ~ 0.17-0.34`: the parent spread is barely larger than the child's, so
refining buys almost nothing — the pixels are genuinely *undetermined*, not under-resolved. That
distinction is what the criterion exists to make.

The region separation is about 1.0 in `alpha`, roughly **twice the `E+1 = 8` noise floor of
0.48**. So at the project's ensemble size the criterion resolves *regions* comfortably and
*individual quads* not at all. A scheduler thresholding per quad at `E+1 = 8` would be acting on
noise; one thresholding on a region aggregate would not.

### Two smaller things

`alpha` for `sigma_E(t)` is **identical to `alpha` for `sigma_E(0)` to every printed digit** in
mid-field and far — the energy spread is unchanged from `t = 0` there, which is what a
well-integrated tame region should look like. In near-field and body2 core it departs and its
max reaches 10.35 and 12.84.

Those large values are a fragility, not a signal: a child quad with near-zero spread makes the
ratio explode. Same `0/0` shape as the `Gamma` residual normalised by `A*B`. Read the median and
the interdecile range; the max of a log-ratio is not a measurement.

---

## 7. Experiment B — which conclusions survive at large `n`?

BRIEF §8's reason for the whole build. Near-field, `t = 13`, `E+1 = 8`, `eta = 0.01`, f64.

### The resolution sweep, and why it is the weaker half

| quantity | 8x8 | 16x16 | 32x32 | 64x64 | 128x128 |
|---|---|---|---|---|---|
| drift median | 8.4847e-9 | 3.6442e-9 | 2.7549e-9 | 2.2884e-9 | 2.1312e-9 |
| drift p99 | 3.2745e-5 | 9.1063e-4 | 9.2568e-3 | 3.5450e-3 | 5.0117e-3 |
| **drift max** | 7.8271e-2 | 1.1508e-1 | 3.0705e-1 | 5.1173e-1 | **1.4909e4** |
| error_ratio median | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| **error_ratio p99** | 1.0095 | 3.5314 | 6.0864e1 | 4.7205e1 | **1.3655e2** |
| spread_shape median | 3.9465e-3 | 2.5833e-3 | 1.9095e-3 | 1.5660e-3 | 1.3244e-3 |
| spread_shape p99 | 1.2417e-1 | 2.1329e-1 | 2.1457e-1 | 1.2669e-1 | 2.9705e-2 |
| d_min_true median | 7.9563e-3 | 8.1360e-3 | 8.3647e-3 | 8.4287e-3 | 8.4693e-3 |
| frac collision | 4.6875e-2 | 3.1250e-2 | 2.5391e-2 | 2.4170e-2 | 2.2095e-2 |
| frac drift > 1e-3 | 1.5625e-2 | 1.1719e-2 | 2.1484e-2 | 1.3428e-2 | 1.4099e-2 |
| frac spread_event > 0 | 9.3750e-2 | 4.6875e-2 | 3.4180e-2 | 2.8320e-2 | 2.3315e-2 |

**A resolution sweep cannot separate "the estimate converged" from "the thing being estimated
changed"** — the jitter scales with cell width, so a finer grid measures a different physical
ensemble. `spread_shape median` falling monotonically 3.9e-3 -> 1.3e-3 is mostly that, not
convergence. Use it for orientation; the subsampling below is the measurement.

### Subsampling one fixed grid — the same physical quantity throughout

Truth is the full 128x128 grid (16384 pixels). Each cell is the **interdecile spread over 200
random draws of `n` pixels, as a fraction of the truth**. Below ~0.1 a conclusion drawn from `n`
pixels is stable.

| quantity | truth | n=16 | n=64 | n=256 | n=1024 | n=4096 |
|---|---|---|---|---|---|---|
| drift median | 2.1312e-9 | 1.544 | 1.070 | 0.851 | 0.445 | 0.227 |
| drift p99 | 5.0117e-3 | 7.716 | 2.021 | 3.051 | 2.056 | 1.146 |
| drift max | 1.4909e4 | **0.000** | **0.000** | **0.000** | 0.003 | 1.000 |
| error_ratio median | 1.0000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| error_ratio p99 | 1.3655e2 | 0.807 | 1.870 | 2.882 | 2.235 | 1.061 |
| spread_shape median | 1.3244e-3 | 0.514 | 0.243 | 0.154 | 0.065 | 0.033 |
| spread_shape p99 | 2.9705e-2 | 4.084 | 4.001 | 6.768 | 5.750 | 3.345 |
| d_min_true median | 8.4693e-3 | 0.370 | 0.234 | 0.140 | 0.060 | 0.034 |
| frac collision | 2.2095e-2 | 2.829 | 2.122 | 1.061 | 0.575 | 0.243 |
| frac er > 10 | 1.6235e-2 | 3.850 | 2.887 | 1.203 | 0.602 | 0.331 |
| frac drift > 1e-3 | 1.4099e-2 | 4.433 | 2.216 | 1.385 | 0.693 | 0.329 |
| frac spread_event > 0 | 2.3315e-2 | 2.681 | 2.010 | 1.005 | 0.545 | 0.251 |

**At `n = 64` — the size of every prior measurement in this project — every fraction has an
interdecile scatter of 2 to 4.4 times the quantity itself.** Two independent studies at that `n`
would routinely disagree by more than a factor of two. That is the 1.2x-turns-out-to-be-18.8x
mechanism, quantified rather than anecdotal.

What is safe at what `n`:

- **medians of well-behaved quantities**: `spread_shape` and `d_min_true` medians reach 0.15 at
  `n = 256` and 0.03 at `n = 4096`. Usable from a few hundred.
- **`drift median`** is slower — 0.85 at `n = 256`, 0.23 at `n = 4096` — because the underlying
  distribution spans twelve orders.
- **fractions** need `n >= 1000` to get under 1.0 and `n ~ 4000` for 0.25.
- **p99 of anything heavy-tailed is not estimable at these `n`.** `spread_shape p99` sits at
  3.3-6.8 at every `n` tested and does not improve. `drift p99` and `error_ratio p99` are at
  1.0-2.9 at `n = 4096`.

**The `drift max` row is the important one and it reads backwards.** It shows 0.000 scatter at
`n <= 256` — apparently the most stable quantity in the table — and 1.000 at `n = 4096`. It is
not stable; the tail is a *single pixel* of 16384, and at small `n` it is essentially never
drawn, so the statistic is stable at the wrong answer. **A max statistic's apparent stability at
small `n` is the statistic never seeing the tail.** That is this project's "a test that cannot
fail is indistinguishable from a test that passes" appearing inside a statistic.

### What 128x128 found that 64x64 could not

Seven pixels of 16384 (0.043%) have `|dE/E| > 1`, one of them `1.4909e4`. They are all finite,
all clustered in one corner of the slice, all with `d_min_true ~ 2e-3` against
`r_coll = 1e-3 R = 2.214e-3`, and all with `gamma_max ~ 1` — the regularised Hamiltonian
residual is order unity, so the trajectory is not being integrated.

| pixel | (jx,jy) | drift_max | error_ratio | d_min_true | gamma_max |
|---|---|---|---|---|---|
| 16110 | (110,125) | 1.4909e4 | 3.9615e8 | 1.8381e-3 | 1.0911 |
| 15989 | (117,124) | 4.0856e1 | 7.3239e5 | 1.7260e-3 | 0.7969 |
| 16351 | (95,127) | 6.3484 | 1.9212e5 | 2.0439e-3 | 0.9829 |
| 16229 | (101,126) | 6.1218 | 1.7145e5 | 2.2125e-3 | 0.9041 |

**`error_ratio` flags 7 of 7.** The field does exactly its job, and this is the strongest
evidence yet for the max-deviation switch: under MAD these would have been invisible.

**And it is a step-size problem, not a wrong equation** — checked, because CLAUDE.md requires it:

| pixel | eta=1e-2 | eta=3e-3 | eta=1e-3 | eta=3e-4 |
|---|---|---|---|---|
| 16110 | 1.4909e4 | 5.9589e-9 | 3.1161e-11 | 1.4885e-11 |
| 15989 | 4.0856e1 | 5.3896e-9 | 3.6647e-11 | 1.3174e-11 |
| 16351 | 6.3484 | 7.9098e-9 | 3.9721e-11 | 1.7923e-11 |
| 16229 | 6.1218 | 9.8056e-9 | 3.5737e-11 | 1.5836e-11 |

Thirteen orders of magnitude for a 3.3x change in `eta`. That is a **cliff, not a slope**: at
`eta = 1e-2` the step lands wrong relative to a close approach the run is not allowed to
terminate on; at `3e-3` it does not.

**Consequence for the spec: `eta = 1e-2` is not sufficient at 128x128.** BRIEF §2.3's production
value was set against grids where those pixels are not grid points. At 64x64 the worst drift is
5.1e-1; at 128x128 it is 1.5e4. The kernel does not fail silently — `error_ratio` catches all
seven — but a production run at `10^6` pixels should expect this class of pixel and either run
at `eta ~ 3e-3` or re-integrate flagged pixels at finer `eta`. **Re-integrating only the flagged
pixels is the cheaper option and needs no scheduler**: the flag already exists, and the
measurement above shows one refinement step is enough.

---

## 8. The refinement remedy, and two limits on it

BRIEF §2.5's remedy, implemented: run the grid, then re-integrate pixels `error_ratio` flags at
`eta/4`, up to three passes, recording the coarse value, the refined value and the `eta` used for
every pixel. Bounded, one extra evaluation of a shrinking subset, no tree and no state carried
between pixels.

Measured on 128x128 near-field:

| | no refinement | with refinement |
|---|---|---|
| drift max | 1.4909e4 | 5.0876e-4 |
| drift p99 | 5.0117e-3 | 3.2498e-5 |
| pixels `\|dE/E\| > 1` | 7 | **0** |
| pixels `\|dE/E\| > 1e-3` | 231 | **0** |
| pixels re-integrated | 0 | 266 (1.62%) |
| wall clock | 11.63 s | 12.00 s (+3%) |

### Limit 1: one pass is not always enough

My claim in PR #6 that one refinement step suffices was measured on near-field only. It does not
generalise. `deep interior` at 256x256 flags 9228 of 65536 pixels — 14% of the region — and one
pass takes it from `1.102e12` to `1.985e1`, still unusable. Three passes reach `1.146e-1` with
**0 pixels still flagged**.

`refine_max_passes` is a parameter and a pixel still flagged after the last pass keeps its
`error_ratio`, so an unrepaired pixel is reported rather than silently accepted. The 256x256
per-region picture, drift max before and after:

| region | refined | before | after |
|---|---|---|---|
| far | 0 | 4.702e-11 | 4.702e-11 |
| mid-field | 0 | 4.129e-10 | 4.129e-10 |
| body1 slice | 0 | 9.542e-8 | 9.542e-8 |
| near-field | 890 | 1.377e7 | 3.120e-4 |
| body2 core | 736 | 3.110e10 | 6.139e-4 |
| deep interior | 9228 | 1.102e12 | 2.543e-1 |

Three of six regions need no refinement at all. The cost is concentrated exactly where the
physics is. At 128x128 `deep interior` goes `1.219e11 -> 1.146e-1`, also with 0 still flagged.

**The pass budget is calibrated on f64.** f32 needs more passes for the same grid: at 128x128
near-field, one pass leaves 421 of 16384 still flagged and six leaves none
(`4.246e1 -> 5.370e-4`); at 256x256 with the default three, f32 leaves **1578 of 65536** still
flagged. Finer `eta` means more steps, and at f32 roundoff accumulation eats into what truncation
error gives back — convergence is slower, not absent. Raise `refine_max_passes` for f32 and read
`n_still_flagged` rather than assuming the default cleared it.

### Limit 2: `error_ratio` detects spread, not drift

After three passes, `deep interior` at 128x128 has **zero** pixels above the flag threshold and a
worst drift of `1.146e-1` — 11% energy error. There is no contradiction: `error_ratio` is
`sigma_E(t)/sigma_E(0)`, and an ensemble whose eight copies drift *together* has a low ratio
however badly they drift. BRIEF §4 already names this in a different form ("it says nothing about
whether the ensemble has decorrelated"); this is the same limitation from the other side.

**So the remedy repairs what `error_ratio` can see.** `energy_drift_max` is dumped per pixel and
is the quantity to threshold on if absolute conservation is what matters. Using `error_ratio` as
the flag was the right choice — it is the field that catches the catastrophic cases, 7 of 7 at
128x128 — but it is not a completeness guarantee and should not be quoted as one.

### A note on what the experiments measure

The experiment examples and every precision-comparison test pin `refine_flagged: false`.
Experiments A and B characterise the kernel whose behaviour motivated the second pass, and
measuring the repaired kernel would hide the thing being measured. Precision comparisons need it
off for a different reason: the pass is threshold-triggered on `error_ratio`, f32 and f64 flag
different pixel sets, and with it on the comparison would be of pipelines rather than of
arithmetic.

### The spread image was hiding its own structure

Unrelated to the physics but worth recording. `ensemble_spread` spans several decades with a
median near `1e-3`, so a linear `[0,1]` ramp painted every grid flat blue and the filament
structure was invisible. Now log scaled between the grid's own p1 and p99, with the window
printed alongside — a false-colour image without its scale is decoration. Non-finite is painted at
full scale: undetermined is the loudest thing a pixel can be and must not be shown as quiet.

---

## 9. The ensemble offsets were the wrong scheme, and fixing it changes less than expected

The spec calls for copy offsets from a **fixed low-discrepancy Halton (2,3) prefix indexed by
copy index**. The port inherited the reference's per-pixel PCG stream: pseudo-random, and
different in every footprint. Two properties were missing — even coverage at small `E`, and the
common-random-numbers structure that should make sampling noise cancel in the parent/child ratio
the refinement exponent is built from.

Implemented as `jitter::Scheme::Halton`, now the default; `Scheme::Pcg` reproduces every result
measured before the switch. The Halton path needs no RNG, no per-pixel seed and no ordering:
`halton_offset(k)` is a pure function of `k`.

### The scheme is decisively better at what it controls

Near-field 64x64 -> 32x32 parents, `alpha` for `sigma_E(0)`, whose true value is exactly 1.0:

| E+1 | scheme | median | p90−p10 | parent/child r |
|---|---|---|---|---|
| 4 | Halton | 1.6678 | **0.0242** | 0.9134 |
| 4 | Pcg | 1.2798 | 0.8397 | 0.1750 |
| 8 | Halton | 1.3857 | **0.0010** | 0.9998 |
| 8 | Pcg | 1.0762 | 0.4796 | 0.1746 |
| 16 | Halton | 1.1668 | **0.0005** | 1.0000 |
| 16 | Pcg | 1.0137 | 0.3203 | 0.2236 |
| 32 | Halton | 1.0661 | 0.0029 | 0.9983 |
| 32 | Pcg | 0.9887 | 0.2123 | 0.2742 |
| 64 | Halton | 1.0308 | 0.0018 | 0.9994 |
| 64 | Pcg | 0.9862 | 0.1493 | 0.3018 |

At `E+1 = 8` the noise floor falls from **0.4796 to 0.0010** — a factor of 480 — and the
parent/child correlation goes from 0.175 to **0.9998**. That is exactly the common-random-numbers
structure predicted, and it is not subtle.

Independently, the L2 star discrepancy of the offset set is lower at every `E` tested
(`tests/halton_offsets.rs`): ratios 0.748, 0.624, 0.489, 0.395 at `E+1 = 4, 8, 16, 32`. **Note
the advantage grows with `E` rather than shrinking** — the opposite of the expectation that it
would be largest at small `E`.

### But it barely moves the exponent anyone cares about

`alpha` for `spread_shape`, same grid, `t = 13`, `E+1 = 8`:

| scheme | median | p90−p10 | var |
|---|---|---|---|
| Halton | 0.1386 | 0.6326 | 5.331e-1 |
| Pcg | 0.1722 | 0.6313 | 5.725e-1 |

The interdecile scatter is **unchanged**. The variance falls by 6.9%, and that number is not a
coincidence: `var(alpha_shape)` drops by `5.725e-1 − 5.331e-1 = 3.94e-2`, against
`var(alpha_E) = 3.75e-2` under PCG. **The sampling noise adds in, it was about 7% of the total,
and Halton removes essentially all of it** — `var(alpha_E)` goes to `1.40e-7`, a factor of
267,000.

The other 93% is not sampling noise. It is chaotic divergence: `sigma_E(0)` is a smooth function
of position, so its ensemble spread is pure geometry and a better-placed offset set fixes it;
`spread_shape` at `t = 13` is dominated by trajectories that have separated, which no offset
scheme touches.

### So RESULTS.md's "per-quad noise floor of 0.48" was over-generalised

I wrote 0.48 as the per-quad floor of the refinement criterion. It is the sampling floor of the
**control**, and the criterion's actual per-quad scatter in `alpha_shape` is **0.63 under both
schemes**, mostly physics. The two are not the same quantity and I presented the first as though
it bounded the second.

**The conclusion it supported survives, with a different number and a different reason.** The
measured region separation is about 1.0 in `alpha`, which is ~1.6x the 0.63 scatter rather than
~2x the 0.48 floor. The criterion still resolves regions and not individual quads — and now for a
reason that cannot be bought off with more copies, since 93% of the scatter is not sampling.

### The +7.6% bias is not a sample-size artefact

That attribution was wrong. Under Halton the bias at `E+1 = 8` is **+38.6%**, and nearly
noise-free (p90−p10 = 0.0010). It falls as `1/E`:

| E+1 | \|median−1\| | x (E+1) |
|---|---|---|
| 4 | 0.6678 | 2.671 |
| 8 | 0.3857 | 3.086 |
| 16 | 0.1668 | 2.669 |
| 32 | 0.0661 | 2.117 |
| 64 | 0.0308 | 1.974 |
| 128 | 0.0173 | 2.220 |

The last column is constant to within ±25% with no trend, so the excess goes as `1/E`.

**The mechanism is geometric, and it is a property of the 2x2 aggregation as a parent surrogate,
not of the estimator.** With fixed offsets, a pooled 2x2 block is four *exact repeats* of one
offset pattern attached to four different cell centres — not `4(E+1)` samples of a genuinely
wider footprint. At small `E` the pooled spread is set by the child-centre separation rather than
by the offset set, and the surrogate diverges from a true parent ensemble. As `E` grows the
offset set fills the footprint and the surrogate converges, which is the `1/E`.

PCG's per-footprint randomisation partially masked this, and matched-count subsampling removed
most of what remained (1.0762 -> 1.0191), which is why I read it as sample-size bias. The
subsampling was removing a *symptom*. Under a fixed scheme there is no randomisation left to hide
behind and the geometry is visible directly.

Two consequences. The bias is **deterministic**, so it is calibratable in a way a noisy bias is
not. And **aggregation is a weaker parent surrogate than I claimed** — "a fine uniform grid
already contains every coarser scale" is true of the *positions* and false of the *ensemble*,
which the fixed scheme exposes.

### The `alpha_E` control variate buys nothing — a clean null both ways

Per-quad correlation between `(alpha_E − 1)` and `alpha_shape`, across quads:

| scheme | rho | fitted beta | floor ratio, regression | floor ratio, additive |
|---|---|---|---|---|
| Halton | −0.0789 | −153.76 | 1.0396 | 1.0008 |
| Pcg | −0.0419 | −0.1637 | 1.0544 | 1.2574 |

`rho` is near zero under both. The regression form makes the floor slightly *worse* in both cases
— fitting a coefficient on noise costs a degree of freedom and returns nothing. **Drop it.**

This is a different correlation from the −0.035 measured earlier, as flagged: that compared two
within-footprint *estimators* per pixel, this compares two per-quad *exponents* across quads.
Different level, different quantity. Both happen to be null, which is a coincidence of this
system rather than the same measurement twice.

The Halton `beta = −153.76` is not a result, it is a division by nothing: `var(alpha_E) = 1.4e-7`
under Halton, so the regression is fitting against a control with no variance left. **The control
variate is degenerate under Halton precisely because there is nothing left to correct** — which
is the outcome the switch was for.

The additive form (`beta = 1`) leaves the Halton floor unchanged (1.0008) and makes the PCG one
26% worse. It does shift the Halton median from 0.1386 to −0.2469, removing `alpha_E`'s geometric
bias — but **that correction is not justified**, because the two exponents do not share the bias:
the geometry that biases `alpha_E` is a smooth-function-of-position effect, and by `t = 13` chaos
has washed it out of `spread_shape`. Measured, the `alpha_shape` median moves only 0.1722 ->
0.1386 between schemes, not by the 0.31 that `alpha_E` moves.

### What was not done

`alpha` was **not** smoothed over neighbouring quads. It is the obvious variance reduction and it
is wrong here: `alpha` varies smoothly except at boundaries, and boundaries are exactly what a
refinement decision is about. Smoothing would blur the signal being detected.

### Two tests the switch broke, and what each one meant

**`seeding_golden.rs`** asserted that distinct pixels get distinct offsets. Under the fixed
prefix that is false *by construction* — it is the property being bought. The test now names
`Scheme::Pcg` explicitly and asserts the contrast in both directions, so the difference is a test
rather than a comment.

**`spread_branch_cut.rs`** failed harder and was worth chasing: the f32 `spread_shape` went from
0.5% off the f64 answer to **30%**. Measured per pixel rather than per mean, near-field 5x5:

| scheme | median rel err | p90 | max | worst pixel |
|---|---|---|---|---|
| Halton | 0.0156 | 0.0825 | **586.2** | f64 2.131e-4 / f32 1.251e-1 |
| Pcg | 0.0042 | 0.0469 | 0.1302 | f64 3.865e-3 / f32 4.368e-3 |

**One pixel of 25.** Its f64 spread is `2.131e-4` — a near-zero denominator — and the two
precisions took different branches at a close approach, so the relative error is 586 and the mean
over 25 pixels moved 30%. The median is 1.6%, and the other three regions are at 0.7-1.2%.

The fix is to the statistic, not the tolerance: the test now gates on the **median** per-pixel
relative error. This project's own rule about aggregates, arriving in its own test suite — a mean
over 25 pixels that one of them controls is not a measurement of the other 24.

The f64 truths also differ between schemes by up to 1.8x (near-field mean `1.606e-2` Halton
against `2.897e-2` PCG). That is expected and is not an error: different offsets measure a
different ensemble. Only the f32-tracks-f64 comparison is scheme-independent, which is why it is
the one gated.

---

## 10. The pooled parent is not a parent

Following the fixed-offset switch, the `+38.6%` bias was traced to the aggregation rather than
the estimator (§9). This settles it by building the thing the aggregation was standing in for.

A **true** parent quad at 2x cell width carries offsets scaled to *its* width. A pooled 2x2 block
carries offsets scaled to the *child's* width, four times over. They are not the same ensemble.
Rendering at `N` and `N/2` makes both sides real. `alpha` for `sigma_E(0)`, true value exactly
1.0, near-field, fine 64x64 against coarse 32x32:

| E+1 | scheme | pooled median | true median | pooled err | **true err** |
|---|---|---|---|---|---|
| 4 | Halton | 1.6678 | 1.0231 | 0.6678 | **0.0231** |
| 8 | Halton | 1.3857 | 1.0227 | 0.3857 | **0.0227** |
| 16 | Halton | 1.1668 | 1.0227 | 0.1668 | **0.0227** |
| 32 | Halton | 1.0661 | 1.0227 | 0.0661 | **0.0227** |
| 4 | Pcg | 1.2798 | 0.8843 | 0.2798 | 0.1157 |
| 8 | Pcg | 1.0762 | 0.9447 | 0.0762 | 0.0553 |
| 16 | Pcg | 1.0137 | 0.9762 | 0.0137 | 0.0238 |
| 32 | Pcg | 0.9887 | 0.9872 | 0.0113 | 0.0128 |

The pooled error runs `0.67 -> 0.07` with `E`; the true error is **flat at 0.0227**. Confirmed:
the bias is the surrogate. The remedy is to render twice, not to calibrate.

The 0.0227 residual is the two-resolution method's own cost. It does not shrink with `E`, and it
reappears independently below as the tame-region `alpha_shape` median of 1.0229 — same constant,
two unrelated quantities, so it is a property of the method rather than of either measurement.

**Under PCG the two errors partly cancel.** Pooling gives 4x the samples and per-footprint
randomisation blurs the surrogate mismatch, so the pooled figure reads *better* than the true one
at small `E` (0.0762 against 0.0553 at `E+1 = 8` is close; at `E+1 = 4` pooled 0.2798 against
true 0.1157 is not). Two wrongs partially offsetting is why this survived until the fixed scheme
removed the randomisation.

### The criterion, re-measured properly — and it is sharper than reported

True two-resolution `alpha_shape`, fixed Halton prefix, `t = 13`:

| region | p10 | median | p90 | interdecile |
|---|---|---|---|---|
| near-field | −0.6568 | 0.0368 | 0.6696 | **1.3264** |
| body2 core | −0.3673 | 0.1844 | 0.7405 | **1.1078** |
| mid-field | 1.0224 | 1.0229 | 1.0235 | **0.0010** |
| far | 1.0228 | 1.0230 | 1.0232 | **0.0004** |

Separation between region medians: **0.9862**.

**"Regions not quads" was too blunt, in both directions.** In the tame regions the exponent is
essentially exact — interdecile `0.0004`–`0.001` — so per-quad decisions there are trivial, and
the ~0.63 scatter pooling reported for those regions was **entirely surrogate error**. In the
chaotic regions the scatter is `1.1`–`1.3`, twice what pooling suggested, and it is chaotic
divergence rather than sampling noise.

That is the criterion working, not failing: **"not resolvable per quad" is the answer for a
chaotic quad.** The scatter is the measurement, not an error bar around one.

### Read the interdecile, not the variance

`alpha_shape` under Halton: variance `5.331e-1`, sd `0.7302`, interdecile `0.6326` (pooled) —
`interdecile / sd = 0.866` where a normal distribution gives 2.563, and **excess kurtosis 110.0**.

The variance lives in the tails; the interdecile describes the bulk. That resolves the apparent
tension between "variance fell 6.9% under Halton" and "scatter unchanged" — both true, of
different parts of one distribution. **A scheduler decides per typical quad, so the interdecile
is the measure and it is the one that did not move.** Do not quote the 6.9% as the improvement.

### Two things left alone

The **Halton advantage growing with `E`** (discrepancy ratios 0.748, 0.624, 0.489, 0.395 at
`E+1 = 4..32`) is reported, not chased. Low-discrepancy sequences are usually most valuable at
small `N`, so growing benefit suggests the raw unscrambled prefix's early terms are less well
distributed than its later ones — a known property in low dimensions. Scrambling is the standard
remedy if it ever matters; nothing here depends on it.

The **control variate** stays dropped. Two variance reductions targeted the same component and
the cheaper one won: the offset scheme had already removed what the control variate would have
corrected, which is why `beta = -153.76` was a division by nothing.

---

## 11. The scheduler — the criterion in a loop

Every measurement before this ran the criterion on **one split in isolation**. `prinq` descends
from one quad, adaptively, at a fixed playhead, and dumps the tree with every decision and reason.

Three regions: `near-field` (mixed), `deep interior` (richest), `far` (tame control). `N = 8`
footprints per quad axis, `E+1 = 8` copies per footprint, so **one quad is 512 trajectories**.

### The sweep had to run first, and it is the result

`tau_display` cannot be chosen before measuring, because the quad spread spans **six orders across
regions**: `~4e-8` in `far`, `~2e-3` in near-field, median `7.5e-5` but p90 `1.4e-1` in
`deep interior`. Any `tau` picked in advance would be the arbitrary constant that has already
disqualified two candidate designs here.

### `alpha_hi` dominates the criterion; `tau` is inert over most of its range

near-field, budget 2000 quads:

| tau | alpha_hi | quads | leaves | floor | keep | budget | depth |
|---|---|---|---|---|---|---|---|
| 1e-8 | **0.20** | 1997 | 1498 | 428 | 201 | **869** | **8** |
| 1e-8 | 0.50 | 25 | 19 | 5 | 14 | 0 | 3 |
| 1e-6 | **0.20** | 1997 | 1498 | 428 | 201 | **869** | **8** |
| 1e-6 | 0.50 | 25 | 19 | 5 | 14 | 0 | 3 |
| 1e-4 | 0.20 | 1997 | 1498 | 423 | 732 | 343 | 8 |
| 1e-3 | 0.20 | 441 | 331 | 131 | 200 | 0 | 8 |
| 1e-2 | 0.20 | 21 | 16 | 0 | 16 | 0 | 2 |

`tau = 1e-8` and `1e-6` produce **identical** trees — the spread never falls that low, so `tau`
never binds. Meanwhile `alpha_hi` from 0.20 to 0.50 collapses the tree **80x**, from 1997 quads at
depth 8 to 25 at depth 3.

near-field's `alpha` median is **+0.389**, so the threshold sits *inside* the distribution and a
small move flips most decisions. That is §4 question 7's "a criterion whose output is dominated by
an arbitrary threshold is not a criterion" — and it is the case that the sibling-reliability policy
(§3.3) exists to sidestep, since that one does not need a trustworthy `alpha` value at all.

Which knob binds is **region-dependent**: at `tau = 1e-4`, deep interior stops at 29 quads (its
median spread `7.5e-5` is below `tau`) while near-field runs to the cap (`2.3e-3` is above it).

### Coarse `N` OVER-splits — the opposite of the stated concern

The concern on record was that too low an `N` makes a quad misclassify itself as **coherent** by
undersampling its own area. Measured, at a fixed budget of 4000 quads and `tau = 1e-4`:

| region | N | traj/quad | leaves | depth | floor | keep | median alpha |
|---|---|---|---|---|---|---|---|
| far | 4 | 128 | 16 | 2 | 0 | 16 | 1.0013 |
| far | 7 | 392 | 16 | 2 | 0 | 16 | 1.0010 |
| far | 8 | 512 | 16 | 2 | 0 | 16 | 1.0010 |
| far | 16 | 2048 | 16 | 2 | 0 | 16 | 1.0010 |
| near-field | 4 | 128 | **106** | **5** | 57 | 49 | 0.2781 |
| near-field | 7 | 392 | 31 | 5 | 14 | 17 | 0.3857 |
| near-field | 8 | 512 | 19 | 3 | 5 | 14 | 0.3451 |
| near-field | 16 | 2048 | **16** | **2** | 8 | 8 | 0.2108 |
| deep interior | 4 | 128 | **40** | **4** | 24 | 16 | 0.3320 |
| deep interior | 7 | 392 | 19 | 3 | 7 | 12 | 0.4941 |
| deep interior | 8 | 512 | 19 | 3 | 6 | 13 | 0.3742 |
| deep interior | 16 | 2048 | **16** | **2** | 5 | 11 | 0.2488 |

Leaf count and depth fall **monotonically with `N`** in both chaotic regions: near-field 106 → 31 →
19 → 16, deep interior 40 → 19 → 19 → 16. A coarse quad does not call itself coherent; it calls
itself **uncertain** and demands refinement. The mechanism is the one §5 warns about — a noisy,
undersampled spread estimate biases toward *refine*, the conservative failure direction — and under
a budget that is not benign: **`N = 4` spends four times as many quads as `N = 16` to cover the same
region.** Cheaper quads, more of them, and the saving is smaller than it looks.

`far` is flat at 16 leaves for every `N`, with `alpha = 1.001` — the tame control both terminates
immediately and reproduces the uniform kernel's tame-region exponent independently.

**The `N = 7` CRN probe is inconclusive, and honestly so.** It was included because the parent–child
common-random-numbers overlap is 32.65% at `N = 7` against exactly 25.00% at every even `N`, so it
is the only lever that varies CRN strength. Measured, `N = 7` sits between `N = 4` and `N = 8` on
leaf count in near-field (31, between 106 and 19) and is *identical* to `N = 8` in deep interior
(19 leaves, depth 3). The `N` trend dominates; no separable CRN effect. It carries the highest
median `alpha` in two of three regions (0.3857, 0.4941), which is the direction more CRN would
predict, but one region disagrees and the effect is inside the scatter. **Reported as not
separated, rather than claimed.**

### Thrash is real, and falls with `N` for a reason the confound cannot explain

Adjacent leaf pairs whose spreads are within a factor of 1.5, at `tau = 1e-4`, `alpha_hi = 0.2`,
budget 4000 quads. "Thrash" is the fraction of those that sit at **different levels**:

| region | N | leaves | depth | similar pairs | diff level | thrash | edge share |
|---|---|---|---|---|---|---|---|
| far | 4 | 16 | 2 | 24 | 0 | 0.0000 | 25.0% |
| far | 8 | 16 | 2 | 24 | 0 | 0.0000 | 12.5% |
| far | 16 | 16 | 2 | 24 | 0 | 0.0000 | 6.2% |
| near-field | 4 | 2998 | 9 | 3579 | 1214 | **0.3392** | 25.0% |
| near-field | 8 | 2998 | 10 | 2529 | 551 | **0.2179** | 12.5% |
| near-field | 16 | 2569 | 10 | 2688 | 197 | **0.0733** | 6.2% |
| deep interior | 4 | 1936 | 10 | 3453 | 31 | 0.0090 | 25.0% |
| deep interior | 8 | 22 | 4 | 20 | 6 | 0.3000 | 12.5% |
| deep interior | 16 | 16 | 2 | 13 | 0 | 0.0000 | 6.2% |

**In near-field thrash is substantial and falls with `N`: 34% → 22% → 7%.** More samples per quad
means a less noisy spread estimate, so neighbours agree more often — which is the direct
confirmation that the thrash is *per-quad noise* rather than real structure.

**The edge-sharing confound cannot explain that trend, and the direction is the argument.** Shared
footprints make neighbours *more alike*, so they *suppress* apparent thrash. The sharing is 25% at
`N = 4` and 6.25% at `N = 16` — so the most-suppressed row is the one showing the **most** thrash.
The true `N = 4` figure is higher than 0.3392, and the fall with `N` is if anything understated.

**`far` reads 0.0000 at every `N` and that is not evidence.** The tree stops at the bootstrap, so
every leaf is level 2 and "different level" is structurally impossible. A uniform tree cannot
thrash. Depth is printed beside thrash for exactly this reason — the same "a test that cannot fail
is indistinguishable from a test that passes" rule, arriving inside a diagnostic of my own.

**`deep interior` is not comparable across `N`.** At `tau = 1e-4` its median spread (`7.5e-5`) sits
*below* `tau` at `N = 8`, so the descent stops at 22 leaves; at `N = 4` the noisier estimate reads
above `tau` and it runs to 1936. Trees of 1936 and 22 leaves do not have comparable thrash
statistics, and the 0.3000 at `N = 8` rests on 20 similar pairs. Reported, not averaged.

### §4 q1 and q2: it terminates, and the floor engages

`alpha_hi = 0.2`, `tau = 1e-4`, **no `max_level`**, budget 50 000 quads:

| region | quads used | leaves | depth | terminated | floored | budget hit |
|---|---|---|---|---|---|---|
| far | **21** | 16 | 2 | 100% | 0.0% | no |
| near-field | **4617** | 3463 | 12 | 100% | 17.6% | no |
| deep interior | **29** | 22 | 4 | 100% | 40.9% | no |

**The descent terminates of its own accord in every region, well inside the budget.** The
Wada-dense-boundary fear — that spread stays high however far you refine, flagged at the outset of
this work and never tested — does not materialise at this playhead and these thresholds. Not one
leaf hit the cap, the depth cap, or the precision floor.

**The floor branch engages, and hardest where it should**: 40.9% of leaves in `deep interior`,
17.6% in near-field, 0% in `far`. A tame region has nothing to floor — its spread is below `tau`,
so it exits through *keep*, which is the correct branch for "already resolved".

near-field's leaf-count-against-iteration is a clean saturation:
`1, 4, 16, 55, 127, 223, 412, 844, 1579, 2536, 3088, 3433, 3463` — the last three iterations add
345, then 30, then nothing.

**What terminates it is `tau`, not the floor.** In near-field's deepest two levels the exponent has
median **+3.945** (p10 +2.256, p90 +5.825) — far above the `alpha = 1` that a halving represents —
while the spread there is `1.4e-5`, below `tau = 1e-4`. Those quads split because their *parents*
were above `tau`, and the split collapsed the spread by ~2^4 rather than 2^1. So the descent ends
by crossing `tau` from above, and 82.4% of near-field's leaves exit through *keep*. The floor is
real but it is the minority branch.

That also means **the termination result is a statement about `tau`**, and `tau` is the knob the
sweep shows to be inert over four orders in this region. Termination at `tau = 1e-4` is not
evidence of termination at `tau = 1e-8`, where the same sweep exhausted a 2000-quad budget with 869
leaves still wanting to split. **Reported as bounded, not as general.**

### §4 q5: the sibling policy is 9x cheaper for the same depth

Equal budget 10 000 quads, `tau = 1e-4`, `alpha_hi = 0.2`:

| region | policy | quads | leaves | floor | keep | depth | median quad spread |
|---|---|---|---|---|---|---|---|
| far | alpha | 21 | 16 | 0 | 16 | 2 | 4.268e-8 |
| far | sibling | 21 | 16 | 0 | 16 | 2 | 4.268e-8 |
| near-field | alpha | **4617** | 3463 | 609 | 2854 | 12 | 5.775e-5 |
| near-field | sibling | **497** | 373 | 235 | 138 | 11 | 7.970e-4 |
| deep interior | alpha | 29 | 22 | 9 | 13 | 4 | 9.446e-5 |
| deep interior | sibling | 21 | 16 | 6 | 10 | 2 | 7.533e-5 |

In near-field the sibling policy reaches depth **11 against 12** for **497 quads against 4617** — a
factor of **9.3**. It floors 63% of its leaves against the alpha policy's 18%, which is exactly the
intended behaviour: where the four sibling exponents scatter, the unreliability *is* the answer and
no trustworthy `alpha` is needed.

**It is cheaper, not obviously better, and the distinction matters.** Its median quad spread is
`7.97e-4` against `5.78e-5` — an order of magnitude more uncertainty left on the table. The alpha
policy spent 4617 quads driving the median spread down by 13x. Whether that was worth it is a
budget question, not a correctness one, and this measurement does not settle it. What it does
settle is that the reliability signal **works as a floor detector**: it identifies the chaotic sea
and declines to spend there.

Leaf overlap in near-field: 264 shared, 3199 alpha-only, 109 sibling-only. The two trees are not
nested — the sibling policy refines 109 leaves the alpha policy does not — so it is not simply a
truncation of the other.

**The caveat stands and is now the next thing to do.** `alpha_sibling_spread` is the **range** of
four samples, which is itself a noisy statistic. This result makes the policy worth pursuing, so
characterising that noise is the next step rather than the first thing to trust.

### §4 q6: priority matters, but the choice of priority does not

The first attempt at this question answered nothing, and the reason is worth keeping. At a budget
of 10 000 the descent terminated at 4617 quads in near-field and 21–29 elsewhere, so the cap never
bound — and **when the cap does not bind, every ordering computes the same set and a jaccard of
1.0000 is structural, not a result.** The run had to be repeated at a budget *below* the natural
termination point. The table now prints a `cap hit` column so the reading cannot be made without
it.

Budget 1500 quads, `tau = 1e-4`, `alpha_hi = 0.2`:

| region | order | quads | leaves | depth | cap hit | vs spread | jaccard |
|---|---|---|---|---|---|---|---|
| far | spread / spread_area / shuffled | 21 | 16 | 2 | **no** | 0 | 1.0000 |
| near-field | spread | 1497 | 1123 | 8 | **yes** | 0 | 1.0000 |
| near-field | spread_area | 1497 | 1123 | 8 | **yes** | **0** | **1.0000** |
| near-field | shuffled | 1497 | 1123 | 8 | **yes** | **600** | **0.5784** |
| deep interior | spread / spread_area / shuffled | 29 | 22 | 4 | **no** | 0 | 1.0000 |

Where the budget binds, **ordering is load-bearing**: shuffling changes 600 of 1123 leaves, 42% of
the tree, at identical cost. But **spread and spread × area produce byte-identical trees** — the
symmetric difference is exactly zero. So the priority function is doing real work and the choice
between these two candidates is free; area weighting buys nothing here, because within one
iteration the frontier is mostly at a single level and the area factor is then a constant.

`far` and `deep interior` show 1.0000 for the other reason — their descents finish in 21 and 29
quads, so a 1500-quad cap cannot bind. Same non-answer as the first attempt, correctly labelled
this time rather than read as a result.

### §4 q3: the tree is sensible in near-field and **not** in deep interior

The overlay was initially drawn over the **outcome** image, which was the wrong base and hid this.
near-field's outcome image is 97.7% one colour, so a tree refining "where nothing is happening"
looked fine. The tree does not track outcome labels — it tracks `ensemble_spread` — so the spread
image is the direct check, and both are now written
(`sched-*_tree_outcome.png`, `sched-*_tree_spread.png`).

Against the spread base:

- **near-field** — dense refinement in coherent bands, sparse elsewhere. Plausibly right in shape.
  But the **brightest, thinnest spread filaments (the lower-left diagonals) sit in coarse quads.**
- **deep interior** — the tree fails. It leaves the large high-spread wedge in the top-left and the
  bright diagonal bands in the lower-right at **level 2**, while spending its only fine refinement
  on an unremarkable patch in the middle. 22 leaves, depth 4, against structure spanning the whole
  box.
- **far** — uniform at level 2, which is correct: there is no structure to track.

**The cause is `tau` and the aggregation together, and both are measured elsewhere here.** Deep
interior's median quad spread is `7.5e-5`, below `tau = 1e-4`, so nearly every quad is kept at
once. And a **median is blind to a thin filament crossing a quad**: most of that quad's footprints
are still in the smooth sea, so the median reads low however bright the filament is. §3.4's warning
was not a formality — it is the mechanism behind the failure in the picture.

**This is the honest answer to "does the tree look right?": in one region yes, in another no, and
the picture is what surfaced it.** It is also a caution about the picture itself — the first
overlay, on the wrong base, would have passed inspection.

### §3.4: the aggregation changes half the decisions — three schedulers, not one

Budget 6000 quads, `tau = 1e-4`, `alpha_hi = 0.2`:

| region | agg | quads | leaves | depth | floor | keep | cap hit | leaf jaccard vs median |
|---|---|---|---|---|---|---|---|---|
| far | median / mean / p90 | 21 | 16 | 2 | 0 | 16 | no | 1.0000 |
| near-field | median | 4617 | 3463 | 12 | 609 | 2854 | no | — |
| near-field | mean | 5997 | 4498 | 9 | 847 | 1139 | **yes** | **0.0963** |
| near-field | p90 | 1141 | 856 | **14** | 472 | 384 | no | **0.0283** |
| deep interior | median | 29 | 22 | 4 | 9 | 13 | no | — |
| deep interior | mean | 161 | 121 | 7 | 38 | 83 | no | 0.1260 |
| deep interior | p90 | 81 | 61 | 7 | 30 | 31 | no | 0.2388 |

**Decision-level disagreement over quads present in both trees** — the figure the brief asks for:

| region | mean vs median | p90 vs median |
|---|---|---|
| near-field | **54.1%** (749 of 1385) | **49.1%** (136 of 277) |
| deep interior | **34.5%** (10 of 29) | **34.5%** (10 of 29) |

**Half the shared decisions flip, and the resulting trees overlap by 3–13%.** This is not a
detail to pick quietly; it is three different schedulers wearing one name. §3.4's instruction not
to choose silently was right, and the reason turns out to be much larger than the kurtosis argument
that motivated it.

Each behaves distinctly, and the behaviours are intelligible:

- **median under-refines structure.** Blind to a thin filament crossing a quad, because most of
  that quad's footprints remain in the smooth sea. This is the mechanism behind the `deep interior`
  failure in the overlay — at `median` it stops at 29 quads and leaves the largest high-spread
  regions at level 2; at `mean` it takes 161 and at `p90` 81.
- **mean over-refines and blows the budget.** 5997 of 6000 quads in near-field, the only
  configuration in this study to hit a cap. Hostage to a single footprint, exactly as expected with
  excess kurtosis 110.
- **p90 refines deepest and narrowest.** Fewest quads (1141) but the greatest depth (14), and it
  **floors 55% of its leaves** against median's 18%. It sees the filament, observes that refining
  does not reduce the extreme, and floors — which is arguably the *correct* answer, since a
  filament genuinely is unresolvable at that scale.

**No recommendation is made here.** The three encode different intentions — "resolve the typical
footprint", "resolve the total", "resolve the worst" — and which is wanted is a display question
this measurement cannot settle. What it settles is that the choice is load-bearing and must be
stated wherever a tree is quoted.

---

## 12. The vertical slice — what the isolation was hiding

Every build before this one was deliberately isolated, and each isolation hid something the next
one found. This build put the camera, the screen floor, the adaptive render, SSAA, and the
linearised decoder together for the first time. Three of the findings below are in the **seams**
rather than in any component, which is where the brief predicted they would be.

### 12.1 A standing rule: a small difference can mean both sides are dead

**Before reading any agreement number, assert that each side still resolves what it is supposed to
resolve.** Two things that have both collapsed agree perfectly, and their agreement carries no
information whatever.

This is not the same rule as "a test that cannot fail" (§ CLAUDE.md), though it is a cousin. That
rule is about a *configuration* in which nothing could have gone wrong. This one is about a
*statistic* that reports success from two simultaneous failures. The tell is different: ask not
"what would make this fire" but "what would make this quantity small", and check that "both sides
lost their data" is not on the list.

Three catches in planning this build alone:

- a **curvature term on an affine chart** — `Slice::decode_pos` is a linspace, so `J_D` is
  constant, `x = x0 + J_D.delta` is exact, and "where does the linearisation start to matter"
  answers "never" at every depth;
- a **linearised f32 sum whose samples all collapse to `x0`** — at depth 40 the increment is
  ~1e-13 against an O(1) centre, so every sample is the same initial condition, agreeing
  perfectly with a direct path that collapsed too;
- an **`E` null that a veto-capped tree would have produced whatever `E` did** — the screen floor
  caps a tree at `4^6 = 4096` leaves, so an over-refining low-`E` run saturates and the sweep
  reports a null the veto manufactured.

The remedy in each case was to measure the *resolving power* first and the agreement second.
`decode::distinct` exists for exactly this: it counts how many of a quad's `N²` sample initial
conditions are actually different, and a divergence figure is only admissible where both sides
are fully distinct.

### 12.2 A collapsed decode makes the criterion maximally confident

The consequence of §12.1 inside the scheduler, and the sharpest seam in this build.

When a decode path has collapsed, every footprint in a quad is the **same** initial condition and
every copy is the **same** trajectory. `ensemble_spread` is then exactly zero. The criterion reads
zero spread as *"this quad is perfectly resolved"* and stops — confidently, with a small tidy
tree, having integrated nothing distinguishable at all.

This is the project's own standing pattern arriving from a new direction: *a statistic can report
maximum confidence precisely when it is least informed*. It has now been caught four times, in
`drift max` scatter, in terminal-outcome purity under lockstep, in the `Gamma` residual, and here.

The guard is not a threshold. It is to count distinct initial conditions per quad and treat a
collapsed quad as **undetermined**, in the same way a non-finite copy is a measurement outcome
rather than missing data.

### 12.3 The deep-zoom floor is a property of *where* you zoom, not of the renderer

PR #11 recorded a plain-f64 cell-width floor at level **45.87**. That figure is conditional on the
chart coordinate being of order 1, and the condition was never stated.

Measured over four configurations at matched settings (`results/output/decode_ladder.txt`), with
64 samples per quad:

| chart | centre | `direct_f64` still resolves 64/64 to | `direct_f32` to |
|---|---|---|---|
| `body_plane` | \|c\| ~ 3 | depth 44 | depth 14 |
| `body_plane` | 0 | **depth 55+ (no floor at all)** | **depth 55+** |
| `shape` | \|c\| ~ 3 | depth 35 | depth 14 |
| `shape` | 0 | depth 45 | depth 45 |

A quad centred at the chart origin has **no O(1) neighbour for the increment to be absorbed
into**, and therefore no cell-width floor in the tested range — on either precision. So 45.87 is
not a universal limit on zoom depth; it is the limit at coordinates of order one, and moving the
same box to the origin removes it entirely. Quote the coordinate magnitude alongside any floor
depth, or the number means nothing.

### 12.4 The linearised decoder buys ~24 levels over f32 and none over f64

The contract's claim is that quad-local relative coordinates extend usable zoom from ~23 to ~50+.
Measured on `body_plane` at \|c\| ~ 3:

| path | all 64 samples distinct to | collapsed to 1 by |
|---|---|---|
| `direct_f32` | depth 14 | depth 22 |
| `L-naive_f32` (the literal formula) | depth 14 | depth 22 |
| `L-split_f32` | depth 44 | depth 50 |
| `direct_f64` | depth 44 | depth 50 |

Two results, and the second is the one to carry:

**The literal formula buys nothing.** `L-naive` — `x0`, `J_D.delta` and the sum all in f32 —
collapses on **exactly the same curve** as forming the chart coordinate in f32 in the first place
(56, 18, 2, 1 at depths 16, 18, 20, 22 for both). Adding a ~1e-13 term to an O(1) f32 quantity is
the same operation as never having the term.

**The split form reaches f64's floor and stops there.** `L-split` — `x0` in f64 on the CPU,
`delta` and `J_D.delta` in f32, promoted and summed in f64 — tracks `direct_f64` rung for rung.
So the gain is ~24 levels *for an f32 consumer*, and **exactly zero over f64**. The "~50+" in the
contract is f64's floor, not something the linearisation creates. That follows from the bound
stated before the run: the initial conditions must be formed as absolute O(1) numbers before
integration, because the three-body separations are O(1) and no nonlinear integrator can carry
`(x0, delta)` separately through the march.

One exception, and it is modest: on the **nonlinear** chart at \|c\| ~ 3, `L-split` holds 64/64 to
depth 45 where `direct_f64` has fallen to 10/64. A single fused affine step loses fewer bits than
a decode through `cos`, `sin` and a renormalisation. Worth ~6 levels, from conditioning rather
than from the design's stated mechanism.

### 12.5 Where the linearisation actually matters is the coarse end, not the deep end

On the shape chart the linearisation error relative to the sample spacing runs `0.39` at
`half = 0.05`, `1.5e-3` at depth 8, `3.6e-7` at depth 20 — it falls as the box shrinks, because
the discarded term is `O(h²)` against a spacing of `O(h)`. It exceeds one sample spacing only at
`half >= 0.5`, i.e. boxes larger than any this project renders.

So the approximation is worst exactly where it is least needed and best exactly where it is used.
That is the opposite of the intuition that a linearisation "breaks down at depth", and it is worth
saying plainly because the intuition is load-bearing in the caching design.

On `body_plane` the same column is **structurally zero**, and is reported as structural rather
than as a measurement. What it does show at depth 44 is `0.55` — that is not curvature but
accumulation rounding reaching half a sample spacing, arriving exactly where distinctness starts
to fail.

### 12.6 "Zero spread" is not zero, and a collapse detector written that way cannot fire

Caught in this build's own instrumentation, which makes it the fourth catch of the same family
and the first one that was mine rather than the brief's.

The first version of `deep_zoom` detected a collapsed decode by testing
`red.spread_median == 0.0`. It reported **no collapse anywhere**, including at depth 40 where
exactly **1 of 64** initial conditions was distinct and every trajectory in the quad was the same
trajectory.

Identical inputs do not give an identically zero spread. `spread_shape` is the mean distance of
the copies' `shape_vec` from their centroid; eight identical unit vectors summed and divided by
eight do not return the value bitwise, so the residual is **5.551115e-17**. Measured directly:

```
  direct_f64: copies distinct 8/8  spread 2.351651e-14  sigma_E(0) 5.329e-15
  direct_f32: copies distinct 1/8  spread 5.551115e-17  sigma_E(0) 0.000e0
```

That residual is **exactly `2^-54`** — one rounding step, `f64::EPSILON / 4` — which is what makes
it structurally unreachable by an equality test rather than merely small. It is twelve orders
below `tau_display = 1e-4`. **No threshold anyone would set can separate a fully collapsed quad
from a perfectly resolved one**, and a small tidy tree built out
of nothing is exactly what a collapsed decode produces.

The fix is not a smaller epsilon. It is to stop asking a *statistic* whether the data was there
and ask the *data*: `decode::distinct` counts distinct initial conditions by bitwise comparison of
all twelve state components, which is exact and cannot drift. The spread is reported beside it, as
the number the criterion would have believed.

**And it connects to a limitation already on record.** `deep interior`'s floored quads carry a
median `worst_energy_drift` of **3.256e-1** — 33% energy error — while `error_ratio` on the same
quads does not flag them. That is exactly the blind spot the design notes record: `error_ratio` is
a ratio of *spreads*, so a drift correlated across the copies cancels in it. Two findings that
corroborate: the flag that should have caught the bad quads is structurally blind to that failure
mode, and `worst_energy_drift` beside it is what caught them. Both fields are needed, as specced,
and neither is redundant.

Note `sigma_E(0) = 0` in the collapsed row. That is a second, independent signal of the same
failure and it *is* exactly zero, because it is a spread of energies rather than a distance from a
computed centroid. It would make a serviceable secondary guard — but it is a symptom too, and the
distinct count is the measurement.

### 12.7 The collapse arrives from the leaves upward, so a root-level check is not enough

Measured in `deep_zoom.txt` at camera depth 14 under `direct_f32`: the **root quad still resolves
all 64 of its samples** while **16 of its 21 descendants have collapsed**. Children sit at half the
parent's cell width, so the failure begins at the leaves — exactly where the scheduler is spending
its budget — and works upward as the zoom deepens.

The consequence for any distinctness guard: it has to be **per quad, at the moment the quad is
computed**, not a once-per-frame check on the view. A frame whose root is fine can be built almost
entirely from collapsed leaves.

And the dangerous case is not the fully collapsed one. At depth 14 the partly-collapsed quads
reported a spread of **1.811e-7** — small, but not absurd, and nothing a sanity check would flag.
Full collapse at least produces the recognisable `5.551e-17`.

### 12.8 Tree size is slice-conditional to a factor of 4; the exponent is not

`slice_variety.txt`, at **one** fixed centre configuration (near-field's, which is Burrau's own),
varying only the 2-plane through it, with bases orthonormal in the 6D position metric so a unit of
chart coordinate moves the system the same distance in every case.

Leaf count spans **226 to 970** — a factor of **4.3**. Among pure rotations within one body's
plane it is only 403 to 526 (a factor of 1.3), so most of the variation comes from **which bodies
the plane moves**, not from the angle. The nonlinear shape chart is the most structured of all at
970 leaves.

The `alpha` distribution, by contrast, barely moves: median 0.172 to 0.289, p10 between -0.09 and
-0.17, p90 between 0.51 and 1.26, across every case including the nonlinear chart.

So the exponent the criterion reads is far more stable than the tree it produces. Two consequences:

- **Every leaf count in this repository is conditional on the slice, to about a factor of 4.** Say
  so when quoting one. The comparisons *within* a slice (with and without the veto, across `E`,
  across aggregation) are unaffected, because they share a slice.
- A criterion tuned on one slice family is more likely to transfer than a *budget* tuned on one.

**The control is exact and the check is on the right quantity.** `plane 0deg` reproduces
`body_plane` with `max |dIC| = 0` and an identical tree. The check compares initial conditions,
not trees: a tree is downstream of a chaotic integration, so checking it would be testing chaos
rather than the charts.

**And a gauge check nobody planned.** The three `shape phase` rows — 0.0, 0.4, 1.3 — are bitwise
identical in every column. The fibre phase is a global rotation and the three-body problem is
rotation-invariant, so they must be; if the Hopf inverse or the AZ port had broken rotational
invariance, they would have separated. Kept, because it costs nothing and it is the only
rotational-invariance check in the suite that runs through the *new* chart code.

### 12.9 The screen floor and `MAX_REL_DEPTH` are different caps, and reporting one hides the other

Raised in the PR #12 review: `deep interior` under mean or p90 reaches depth 7, and the screen
floor at `N = 8` on a 512² viewport sits at level 6. Does the aggregation fix collide with the
veto in exactly the configuration production runs?

**Yes, and it separates the two aggregations.** Measured in `agg_vs_floor.txt`:

| viewport | `MAX_REL_DEPTH` | agg | leaves | depth | cost of the cap |
|---|---|---|---|---|---|
| 512² | 6 | mean | 79 | 6 | **−42 of 121, 34.7%** |
| 512² | 6 | p90 | 58 | 6 | **−3 of 61, 4.9%** |
| 1024² | 7 | mean | 121 | 7 | none |
| 1024² | 7 | p90 | 61 | 7 | 4 leaves screen-floored |

Both keep an identical tree at levels 2–5; the whole difference is what piles up at the cap. So
**p90's fix survives the production viewport nearly intact and mean's does not** — a reason to
prefer p90 that has nothing to do with the statistic's own properties, and one neither the design
docs nor PR #11 could have found, because neither had a camera.

The collision is a *resolution* limit rather than a design conflict: `4^7 x 64 = 1,048,576`, so
level 7 is displayable at 1024² and the aggregation the region needs is affordable one viewport
step up.

**The trap, and it is the point of this note.** At 1024² with `MAX_REL_DEPTH` left at its default
6, the tree is **identical to the 512² tree** and the `screen` column reads **zero**. That reads
as "the viewport made no difference". It is wrong: `MAX_REL_DEPTH` had taken over as the binding
cap. The two coincide at 512² by construction — the contract's `MAX_REL_DEPTH <= screen floor`,
with 6 chosen to match — and **diverge at every larger viewport**, where the default silently
becomes the tighter of the two.

The first version of this run reported only `screen`, showed it falling to zero while the tree did
not grow by one quad, and would have been written up as "the viewport is inert". Two caps, one
column, wrong conclusion.

**And the two regions answer the collision question oppositely.** `deep interior` is
**criterion-bound** — the veto touches 4 of p90's 61 leaves at 512² and none at 2048², so a
viewport step hands the region back to the criterion. near-field is **view-bound at every viewport
tested**: at 1024² with `MAX_REL_DEPTH = 7` it still floors 576 of median's 844 leaves, 756 of
mean's 988, 88 of p90's 271; at 2048² with `MAX_REL_DEPTH = 8`, 2172 of mean's 2617 and 148 of
p90's 382. Uncapped, p90 reaches **depth 14** there. Its structure is dense at every scale, so more pixels buy more tree
and never reach the point where the criterion decides. A question about "the" collision has no
single answer; it is a per-region property and has to be reported as one.

**So: a scheduler's depth cap is two numbers, and a tree quoted with only one of them is
underspecified.** `MAX_REL_DEPTH` is a policy default; the screen floor is arithmetic. State both
wherever a tree is quoted, the same way §11 requires the aggregation to be stated.

---

## 13. Improving the criterion — what the brief got right, and where it got there differently

### 13.1 §1's conclusion is right and its mechanism is not

The brief reads `ensemble_spread` as a category error: a *within*-footprint statistic where only
*between*-footprint variation is reducible. Two things in the repo say the premise does not
describe this implementation.

`jitter_frac` defaults to **0.5** and `halton_offset` returns `[-1, 1)^2` scaled by cell width
per axis. So the copies span **±0.5 cell widths — the whole cell, edge to edge**. They are a
quasi-random sample of exactly the area the footprint stands for, not a cloud around a point.
The design-record quote §1.2 leans on — *"the ICs there are identical up to perturbation"* — is
about a different construction.

And the corroboration was already measured before this build: the Halton control's true `alpha`
is **exactly 1.0**. `alpha = log2(spread_parent/spread_child)` can only be 1 because splitting
halves the jitter ball along with the cell; an irreducible within-point statistic would have
`alpha == 0` by construction.

Measured directly, with counts matched:

| region | rho all | rho mix | scale (matched count) | count (matched extent) | n mix |
|---|---|---|---|---|---|
| near-field | 0.7240 | 0.5818 | 1.172 | 1.013 | 192 |
| far | 1.0000 | — | 9.555 | 1.003 | 0 |
| deep interior | 0.6828 | 0.6361 | 2.062 | 1.012 | 27 |

`count` is `within_pooled / between_shape` at **equal extent**: 1.01 everywhere. Matched for
extent and sample count the two arms are the **same estimator**. `scale` is
`between_matched / within` at **equal sample count**: 1.17 in the chaotic region and **9.56** in
the tame one, which is what a smooth field gives when the window widens by `N-1 = 7`, and what a
saturated field gives when widening buys nothing.

But `rho mix` — quads whose hot set is a proper subset, i.e. containing a transition — is only
**0.58–0.64**. At their actual settings the arms rank quads materially differently. So §1's
practical conclusion survives; what changes is the mechanism, and with it the fix: the fault is
the **aggregation**, and §3.1/§3.2 address it, not a new arm alone.

Two guards were needed to get a readable number. An unstratified `rho` is dominated by tame
quads where both arms read near zero and agree trivially — it would read high whatever happened
at the boundaries. And the first stratification chosen (`looks_like_boundary`) had a population
of **2, 0 and 4**: in a chaotic region most quads are *uniformly* hot, so there is no internal
hot/cold edge and they are correctly not boundaries. `med hot` is 1.000 in near-field. The
`mixed` stratum is the weaker, better-populated one.

### 13.2 The metric, and the shipped criterion coming last

`error(B)`: reference = the fully-refined tree at one sample per pixel, error = mean per-pixel
OKLab distance, criteria entered as **orderings** with no threshold consulted.

That last choice is not cosmetic. §13.1's `scale` factor means a threshold comparison would
score the 1.17-vs-9.56 rescaling instead of the signal. A ranking is invariant to any monotone
rescaling, so the confound simply does not arise.

`deep interior`, the region that matters, oracle-to-random gap wide enough to discriminate:

```
               B=      191       767      1535
  greedy_oracle     0.01386   0.00240   0.00004
  frac_hot_between  0.01446   0.00366   0.00002
  running_max       0.01786   0.00425   0.00158
  between/median    0.01882   0.00560   0.00003
  within/median     0.01861   0.01509   0.01386   <- shipped default
  random lo         0.01814   0.01288   0.00932
```

`within/median` is beaten by **random** at every budget past 383, in both regions. In near-field
it is flat at 0.00394 to `B=767` while `within/mean`, `between/median` and `max_of_both` all
reach the oracle's zero at `B=191`.

`far` **cannot be measured at all**: `error(root) = 0.00000`. The outcome image is featureless
at 512², so every criterion reads zero and none of it is data. That is the metric's own guard
firing, and it is reported as undefined rather than as agreement. It also reframes every earlier
leaf-count comparison on `far`: there was never an image there to get right.

### 13.3 A flat curve has two causes and they need different fixes

I predicted `within/median`'s flatness was a degenerate ranking. **It is not.**

| signal | distinct of 5461 | modal% | reading |
|---|---|---|---|
| within/median | 5418 | 0.3% | fine-grained, and actively bad |
| within/mean | 5461 | 0.0% | fine-grained, and good |
| frac_hot_within | 58 | 40.8% | degenerate — no ordering |
| layout | 78 | 40.8% | degenerate — no ordering |
| frac_hot_between | 65 | 33.9% | degenerate, **and the best criterion in `deep interior`** |

A bad ordering and no ordering produce the same flat `error(B)`, and the curve cannot separate
them. Counting distinct values can, which is why it is printed **above** the curves.

The last row is the one that resists a tidy story: a 65-valued signal beats a 4994-valued one.
Resolution is not what makes a ranking good. Nor is coverage — `term_grad` is **NaN on 97.1%**
of near-field and still reaches the oracle's zero by `B = 383`, because the 2.9% it does score
are exactly the structured quads.

### 13.4 Four measurements that could not have failed, caught in one build

- **The `sigma_E(0)` control for `alpha_sibling_spread`.** True `alpha` exactly 1.0, true sibling
  range exactly 0, no integration — and it reads **0.003, flat in both `N` and `E+1`**. The
  flatness is the tell: under the fixed Halton prefix the offsets *and* the footprint positions
  are fixed, so the quantity is deterministic and there is no sampling noise in it to measure.
  Kept as a geometric floor, labelled as one; part 2 varies a `Pcg` seed for a real draw.
- **`between_shape == 0.0` as a collapse test.** A genuinely uniform region has zero spread over
  perfectly distinct ICs. Tested on `decode::distinct` instead, which cannot confuse the two.
- **An unstratified `rho` between the two arms** (§13.1).
- **The `escaped` fraction.** `t_end` is set by whichever terminating event came first, so
  `deep interior` reads **0.99 terminated with the escape arm silent** — collisions. Reporting
  that as an escape fraction would have contradicted the standing "zero of 1024 near-field
  pixels escape at `t = 13`" while appearing to agree with it. `terminated_fraction` and
  `escape_fraction` are now separate columns.

### 13.5 §5's accumulators: the event arm already existed, and the shape arm is not a null

The brief says the temporal accumulators are specced and missing. Half of that is wrong:
`spread_event_max` is a running max over boundaries, `t_spread_event` is a first-divergence time
that is NaN rather than `t_max` when it never fires, and `spread_event_latched` is the
persistence-guarded latch with `LATCH_RUN = 3`. What was absent is the **continuous** arm.

Built, it is **not** the clean null the short horizon made plausible. `running_max` reaches
0.00158 at `B=1535` in `deep interior` where `within/median` sits at 0.01509 — third best of
everything tested there. `first_divergence` is degenerate (63 distinct, 82.3% modal) and
middling.

The driver had no extension point: `cart` is overwritten in place every boundary and
`AzOut::state` is final-only. `boundary_shapes` is pushed beside `tight`, behind a flag, and
reduced inside one footprint's evaluation so peak cost is one footprint's worth. It is **ragged**
— copies terminate at different boundaries under `stop_on_event` — and the reducer carries each
copy's last recorded shape forward rather than truncating to the shortest, which would discard
exactly the boundaries where the survivors are diverging.

### 13.6 §7's questions are mostly identities in the current design

`Camera::veto` reads `tile_size_px`, which depends on the quad's width and the camera's
`half_world` and `viewport` — **and not on `cx`/`cy` at all**. There is no view culling and no
cache. So "does the tree persist across a pan" is an identity, and reporting it as a finding
would be reporting one.

What is left open, and is measured: what *would* be evictable (`Camera::covers`, a predicate the
scheduler never consults), and whether a quad recomputed after leaving view comes back
**bitwise** what it was — the property a future cache would need to be sound.

---

## 14. Render the diagnostic field, not the science field

Three bugs in this sequence began the same way: an image looked wrong, the reference did not have
the artefact, and the cause was numerical rather than physical. The LC branch cut, the
escape-termination patchwork, and the `dtau` step control. Each was dismissed as physics at least
once.

The general lesson from the third is about *which field to look at*.

The production colouring is bivariate — hue from the shape sphere, lightness from a scalar — and
both channels are **science** fields. A numerical defect reaches them only after it has propagated:
it has to corrupt a trajectory badly enough to move a spread or flip an outcome label before
anything appears, and what appears then is a scatter of pixels that looks exactly like fractal
mixing. That is why the `dtau` blow-up survived a whole corpus.

`energy_drift_max` is already in the payload per footprint. Mapped directly — inferno ramp,
magenta for no-value, auto-ranged over its own p2-p98 — it showed **coherent arcs** of high drift
with the non-finite pixels sitting *inside* them: clustering 1.7x over chance, and 6 of 6
non-finite pixels having a high-drift neighbour. The defect at source, in one picture, on a field
that costs nothing to render because the number was always there.

So: `Scalar::Drift` and `colour::drift_rgb` are part of the standard render set alongside outcome
and spread, and `_drift.png` is written for **both** arms of any before/after — the "before" map is
the artefact worth keeping, because it is what the signature looks like for the next person.

Two cautions that came out of using it:

- **A ratio and its base rate move together.** The clustering ratio *rose* under the fix (2.381 ->
  3.868) while the counts fell (819 -> 435 hot, 11 -> 2 non-finite). What the fix removes is the
  diffuse population; what survives is the genuinely clustered core. Read the count beside the
  ratio or the ratio reads backwards.
- **Auto-ranging hides magnitude.** The p2-p98 window is per panel, so a clean field and a blown-up
  one both fill the ramp. The window is printed beside every panel for that reason —
  `(9.885e-10, 3.624e5)` before against `(8.559e-10, 4.387e4)` after is the improvement, and neither
  image says it on its own.

- **AND A DIAGNOSTIC FIELD IS SPECIFIC TO A CLASS OF DEFECT.** The boundary-overshoot bug
  (`RESULTS §24`) is invisible to this panel. The clamp buys **24,000x** on the figure-eight
  closure error and raises the convergence order from 1.13 to 3.06, while moving `near-field`'s
  median energy drift **37x the wrong way** (1.5e-9 -> 5.6e-8) and the NumPy smoke median from
  3.197e-9 to 4.047e-9. The overshoot displaces the state in *time*, and the AZ energy is nearly
  stationary along the flow, so a time displacement barely registers in `|dE/E|`. Energy drift
  finds a step that grew too large; it does not find a step that ended in the wrong place. Ask
  what the diagnostic would say about the defect you are looking for *before* reading it as clean
  -- the same question as *what would have to be true for this test to fire*, asked of an
  instrument rather than an assertion.

  For that class, the instrument is a **convergence order** on a system with a known exact answer.
  The Chenciner-Montgomery figure-eight is exactly periodic, so `|state(T) - state(0)|` is a pure
  error with no reference trajectory to integrate and no chaos to contaminate it, and it runs in
  under a second. **Read the order, not the error** -- an error falls for many reasons; only the
  order says the leading term changed.


---

## PARKED: the logarithmic Hamiltonian (logH), and why it is the next thing

Queued behind the 1024^2 gallery and `far_hierarchy`. Written down rather than left in a
conversation because it changes what the Heggie result means.

### What the literature check found

The AZ-vs-Heggie comparison is **not new**. Mikkola's *A Comparison of Regularization Methods for
few-body interactions* puts all three side by side in its Figure 6 — Heggie's global method, CHAIN,
and what he calls **"Zare's cartwheel method"** — and gives Heggie's algorithm in full (§5,
Eqs. 16-22), matching what `src/integrate/heggie/` transcribes.

**CHAIN has the same re-selection problem, checked more often than AZ checks it.** §7: *"There is
also a check for need of a new chain after every integration step"*, a new chain forming whenever an
unchained vector is shorter than the smallest chained vector at either end. So "the chart is rebuilt
as the configuration changes" is structural to CHAIN too, not peculiar to AZ.

**What could not be found in two searches:** any quantified cost of chain-switching or reference
re-selection — the literature says *when* the chain is rebuilt, not what rebuilding costs in error —
and nothing evaluating any method on coherence across **neighbouring initial conditions**. That is
consistent with the axis being unexplored, but two searches is not a literature review. Aarseth's
book chapter and Mikkola & Aarseth (1993) are the next reads.

### The gap in OUR work, and it is the important part

**Mikkola's own recommendation is neither method.** It is the logarithmic Hamiltonian

```
    Lambda = ln(T + B) - ln(U),      B = U - T
```

integrated with a leapfrog, or its non-canonical generalisation TTL. These are **algorithmic**
regularisation: a time transformation and a good integrator, with **no coordinate transformation at
all**. His words: *"In more general cases one may recommend the logH-method."*

And logH has **no chart whatsoever** — no reference body, no chain, nothing to re-select. That is a
*stronger* form of the property the Heggie win has been attributed to. **If the re-registration
mechanism is the real story, logH should match or beat Heggie.** If it does not, the mechanism is
wrong and the Heggie result needs a different explanation.

### Why it is cheap to try

`src/integrate/leapfrog.rs` already exists. logH is that plus a time transformation and a `B`
accumulator; Mikkola's toy implementation is about thirty lines:

```
    X(s):  dt = s/(T + B);  r += dt * p;  t += dt
    V(s):  p += (s/U) dU/dr;  B += (s/U) dU/dt
    step:  X(h/2) V(h) X(h/2)
```

Leapfrog alone is not enough — Mikkola uses it as the base for Gragg-Bulirsch-Stoer extrapolation.
**That is a real scope item and not a footnote**: comparing a GBS-extrapolated logH against an RK4
Heggie would be scoring the integrator, not the regularisation, and this project's whole AZ/Heggie
comparison is fair *because* both run the same RK4 under the same step control. Either give logH the
same RK4 treatment and accept it is not how the method is meant to be used, or run all three under
GBS. Say which.

### What it would settle

  1. Whether the re-registration mechanism explains the Heggie win, or only correlates with it.
  2. Whether the whole KS/coordinate-transformation family is the wrong branch for this application.
  3. `far` — logH is chartless, so if it also loses `far` the cause is not the chart at all.

### BUILT. Predictions recorded before the 256^2 run, and what was already seen

logH is implemented: `src/integrate/logh/`, four `Integrator` variants (`LogHLeapfrog`,
`LogHRk4`, and the two `Plain*` controls, which are the same code path with `LhTime::None`), 22
tests across `logh_hamiltonian_fd`, `logh_march` and `logh_seam`. Phase D is
`examples/logh_arms.rs` at 256^2 over the six cases the 1024^2 gallery finished.

**The four predictions, in the state they were written before the run:**

  1. `logh_lf >= heggie` overall, if re-registration is the mechanism.
  2. logH should also win `far`, where AZ wins today on all 65536 pixels by a flat 0.7-0.9
     decades. That win is attributed to AZ never re-registering there; a method that never
     re-registers anywhere should not lose it.
  3. The drift field tracks FTLE positively, as the leapfrog's does and AZ's does not.
  4. **Added during planning:** `logh_rk4` behaves as a Sundman-transformed RK4 and loses most
     of the advantage, because `K + B == U` on shell and RK4 evaluates both denominators at the
     same point. If it matches `logh_lf`, this is wrong and so is Mikkola & Merritt's sentence
     it rests on.

**AND PREDICTION 3's TARGET NEEDED CORRECTING BEFORE IT COULD BE USED.** The number on record is
FTLE-vs-drift `+0.3048` for leapfrog and `-0.0820` for AZ. Read against the shifted controls in
the same table (`+0.0240` and `-0.1022`), **AZ's is a null, not a negative correlation** -- its
control is larger in magnitude. And the FTLE is computed *on the plain leapfrog*
(`src/physics/ftle.rs:26`), so the row that works has the **same integrator on both sides**. The
honest target is "leapfrog's drift tracks FTLE and AZ's does not", with the shared-arm confound
carried alongside any number quoted. Fixing it properly means an FTLE per occupant, which is
propagating a tangent vector through a regularised chart and is a separate build.

**WHAT HAS ALREADY BEEN SEEN, AND WHY IT IS NOT A RESULT.** Two things landed during
implementation and both bear on the predictions, so they are recorded here rather than
discovered later in a table:

  - **Prediction 4 is already supported, on the sharpest fixture available.** KDK traverses the
    two-body radial collision for all three pairs at identical step counts; RK4 does not
    complete at any step size tried. That is `tests/logh_march.rs`, not a field.
  - **Prediction 2 looks likely to FAIL.** A 24^2 smoke run of `logh_arms` on `far` put
    `logh_lf` **5.3 decades behind Heggie** and `logh_rk4` 3.2 behind. That is a 576-pixel grid
    at `max_steps = 60000`, run to check the harness prints, and this project has been burned
    twice by coarse grids -- one understated a maximum eightfold, another overstated a median
    twenty-sixfold. **It is not a result and is not quoted as one.** It is recorded because
    seeing it before the real run is exactly the circumstance in which a prediction quietly gets
    softened, and the prediction above is left as written.

**And a third thing that is a result, from the tests rather than the field:** logH does **not**
meet BRIEF §5's collision gate. `d_min < 1e-10` yes at `eta <= 3e-5`; `|dE/E| < 1e-12` never --
its drift *rises* with penetration depth, 1.1e-9 to 2.3e-6, where Heggie reads `5.422e-27` with
`drift_reg` flat at 4.4e-15. Heggie's KS map removes the `1/r` from the Hamiltonian; logH only
slows the clock, so the encounter has to be **resolved** rather than removed, and there is no
second, better-conditioned energy to report instead. **A chartless method is not a coordinate
regularisation with the coordinates left out**, and if prediction 1 fails this is the first
place to look for why.

### THE FALSIFICATION TEST IS CONFOUNDED, AND THE CONFOUND IS MEASURED

Two of six cases at 256^2, unmasked kernel, diagnostic pass:

```
                     drift p50   vs heggie   frac better    evals p50   err>10
  far      az         2.824e-13      +0.89        1.0000      3.938e5        0
           heggie     2.189e-12          -             -      2.354e5        0
           logh_rk4    3.271e-9      -3.17        0.0000      2.345e5        0
           logh_lf     4.709e-7      -5.33        0.0000      1.477e5        0
           plain_rk4    3.165e2     -14.16        0.0000      1.056e5    65536
           plain_lf     3.130e5     -17.15        0.0000      1.056e5    65536

  deep_    az          1.345e-9      -0.15        0.4287      9.479e5      438
  interior heggie     9.298e-10          -             -      4.642e5        0
           logh_rk4    2.650e-8      -1.49        0.0168      4.623e5      129
           logh_lf     4.541e-5      -4.68        0.0001      1.942e5      598
           plain_rk4    2.240e2     -11.36        0.0000      1.056e5    65536
           plain_lf     1.869e2     -11.31        0.0000      1.056e5    65536
```

**Prediction 1 fails on both, by 5.33 and 4.68 decades.** The plan says that means the
re-registration attribution is wrong. **It does not, and the reason is now measured rather than
suspected.**

The syllogism was: logH has *no* re-registration, which is stronger than Heggie's *none after
the first*, so if re-registration is the mechanism logH should match or beat it. That is valid
only if logH is otherwise at least as good as Heggie, and it is not. `tests/logh_march.rs`
already shows logH missing BRIEF §5's collision gate -- `d_min < 1e-10` yes, `|dE/E| < 1e-12`
never, drift *rising* with penetration depth where Heggie's `drift_reg` is flat at 4.4e-15 --
because a KS map removes the `1/r` from the Hamiltonian and a time transformation only slows the
clock.

And **every region tested is collision-dominated**: `results/output/integrator_gallery.txt` gives
`coll` of 1048576/1048576 on `far`, 1033184 on `deep_interior`, 850590 on `preset_shape`. So
logH is being graded almost entirely on the one thing it is known not to do, and its deficit
there says nothing about re-registration either way. **This is the same class of error the
project keeps catching -- two things compared that differ in more than one way -- arrived at from
the other side: not a knob held fixed whose effect goes unattributed, but a knob NOT held fixed
that was treated as if it were.**

**THE UNCONFOUNDED MEASUREMENT IS CADENCE SENSITIVITY, NOT ABSOLUTE ACCURACY.** Double `n_sync`
at fixed step size (scaling `eta`, and `closure_k` with it) and measure how far the drift field
moves. AZ moves **0.44-0.52 decades**, Heggie **0.048**. That is a *difference within one
integrator*, so each occupant's own singularity handling cancels out of it, and logH -- which
re-registers zero times -- must land at Heggie's level or below if the mechanism is real,
however poor its absolute drift is. `examples/heggie_machinery.rs` is the harness; logH arms go
in after the sweep.

### AND PREDICTION 4 IS REFUTED IN THE DIRECTION I DID NOT CONSIDER

`logh_rk4` was predicted to lose most of the method's advantage, degenerating to a Sundman
transformation because `K + B == U` on shell. It does not lose to `logh_lf` -- it **beats** it,
by 2.16 decades on `far` and 3.23 on `deep_interior`, and by ~1.8 and ~2.5 after correcting for
the evaluation shortfall. The prediction assumed `logh_lf` was the ceiling.

The arithmetic behind it may still hold; its predicted *consequence* is false, because at these
depths a fourth-order stepper's accuracy outweighs whatever the leapfrog's structure buys. The
one place the structure is decisive is the **exact** collision, where RK4 does not complete at
all and KDK traverses for all three pairs. So there is a crossover in encounter depth, measured
at both ends and not in between.

### THE STEPPER-ONLY CONTROL GIVES OPPOSITE ANSWERS ON THE TWO CASES, WHICH IS WHY IT IS RUN

At exactly matched evaluations (1.056e5 on both arms, both cases): `far` reads `3.165e2` against
`3.130e5`, **three decades**, and `deep_interior` reads `2.240e2` against `1.869e2`, **0.08
decades**. So the stepper contributes heavily where the field is dominated by one clean
two-body encounter and essentially nothing where it is chaotic. A single case would have
licensed either "the stepper is in every comparison" or "the comparison is clean", and both
would have been wrong somewhere.

### CORRECTION, ONE CASE LATER: THE COLLISION EXPLANATION IS NOT SUPPORTED, AND THE CONTROL IS SATURATED

`near-field` landed and refutes the mechanism I offered above for logH's deficit. Decades against
Heggie, with each region's collision fraction beside it:

```
             case    coll   hg/az   logh_lf   logh_rk4   lf-rk4    plain_rk4    plain_lf
              far   1.000   -0.89     -5.33      -3.17    -2.16    3.165e+02   3.130e+05
    deep_interior   0.985   +0.16     -4.69      -1.45    -3.23    2.240e+02   1.869e+02
       near-field   0.024   +0.65     -4.94      -2.09    -2.85    1.486e+02   9.576e+01
```

**`logh_lf`'s deficit is flat at 4.7-5.3 decades while the collision fraction runs from 100% to
2.4%.** So "logH is being graded almost entirely on the one thing it is known not to do" is
wrong: the deficit is the same where collisions are rare. Withdrawn. The `logh_rk4` gaps do not
order by collision fraction either -- 3.17, 1.45, 2.09 against 1.000, 0.985, 0.024.

**AND THE STEPPER-ONLY CONTROL IS SATURATED ON EVERY FIELD CASE, SO ITS DIFFERENCES MEASURE
NOTHING.** Both `plain_*` arms carry `err>10` on **all 65536 pixels** of all three regions --
this project's own flag for *this pixel is not data* -- with drift at `1e2` to `1e5`. The "three
decades on `far`, 0.08 on `deep_interior`" quoted above is a difference between two meaningless
numbers, and the conclusion drawn from it (that the stepper's contribution is case-dependent) is
withdrawn with it. A fixed step of `4e-3` cannot resolve any of these regions, which is exactly
what `the_control_is_a_fixed_step_integrator` says it is.

Where the control is **not** saturated it does measure the stepper: on the figure-eight,
`LhTime::None` reads order 4.52 under RK4 and 1.85 under KDK, and at matched evaluations RK4 is
ahead by about **2 decades**. That is the same size as the measured `logh_lf` against `logh_rk4`
gap (2.16, 3.23, 2.85), so most of the leapfrog arm's deficit is its **order**, not its
regularisation.

**Which makes the like-for-like comparison `logh_rk4` against `heggie`, both RK4 under the same
step control: -3.17, -1.45, -2.09 decades.** That is the number the regularisation question turns
on, and it is two to three decades rather than five.

**What survives from the entry above:** the general confound -- logH differs from Heggie in more
than one way, so an absolute-accuracy comparison cannot isolate re-registration -- and therefore
the case for the cadence measurement, which is a difference within one integrator and cancels
whatever each occupant's absolute accuracy is. What does not survive is my account of *which*
difference was doing the work. I named collisions, and the collision fraction says no.

### THE SWEEP LANDED, AND THE SIGNAL IS THE VARIANCE, NOT THE ACCURACY

Six cases, eight arms, 256^2, diagnostic pass, unmasked kernel. Decades against Heggie, sorted,
with the **spread across the six regions**:

```
          az:  -4.50  -3.02  -1.67  -0.65  -0.16  +0.89    spread 5.39 decades
    logh_rk4:  -3.17  -2.09  -1.77  -1.76  -1.57  -1.45    spread 1.72
     logh_lf:  -5.33  -4.99  -4.94  -4.87  -4.69  -4.21    spread 1.12
```

**AZ's deficit against Heggie swings 5.4 decades across the corpus; logH's swings 1.7, and the
bare leapfrog's 1.1.** Both AZ and logH are measured against the same no-re-registration
baseline, so the comparison is like for like -- and the arm *with* re-registration is the one
whose error depends on the field, by a factor of three to five.

That is the strongest evidence in the whole run for the re-registration story, and it is **not**
the axis the falsification test was written on. Absolute accuracy asks "is logH better"; the
answer is no. Variance asks "does the chart machinery make the error field-dependent"; the answer
looks like yes. The leapfrog arm being the *most* consistent of the three fits the same reading:
its error is dominated by its ORDER, which is a property of the method, where AZ's is dominated
by events that depend on where in the chart the trajectory goes.

**Suggestive, not conclusive.** Six points, and a spread is a weak statistic on six points; the
cadence measurement remains the pre-registered test, because it is a difference *within* one
integrator rather than a spread across regions.

### AND logH BEATS AZ ON THE CHART CASES, WHICH THE FIRST FOUR ROWS DID NOT SHOW

```
    config_stability   -0.10       preset_shape      +1.46
          near-field   -1.44       preset_shape_h1   +2.74
                 far   -4.06
       deep_interior   -1.29
```

`preset_shape_h1` is the chart where AZ carries **2222 not-data pixels**; `logh_rk4` carries
**9**, and beats it by 2.74 decades. `preset_shape`: AZ 703, logH 11. **The robustness column is
where logH's advantage over AZ is unambiguous and it is large on every case** -- `err>10` of
423 -> 7, 438 -> 129, 703 -> 11, 2222 -> 9 -- even on the cases where its drift is worse.

Reporting the first four cases as "logH is competitive with AZ" understated it: on the two
presets it is decisively ahead. A four-case read of a six-case corpus, and the two missing cases
were the ones that differed.

### AND THE `plain_*` SATURATION IS 86%, NOT 100%, ON `config_stability`

The claim "err>10 on all 65536 pixels, both arms" was measured on `far`, `deep_interior` and
`near-field`, where it holds. `config_stability` reads 56199 and 47012. No practical difference
-- a drift of `8.4e2` is not data either way -- but the claim was about three cases and does not
generalise to the fourth.

### THE LANDING FIELD TEST USED A DIAGNOSTIC THIS PROJECT ALREADY DOCUMENTS AS BLIND TO IT

Both guards passed and the answer was flat:

```
    near-field  Rk4  land off   drift 3.292e-7   land resid 1.488e-5   corr 0
    near-field  Rk4  land on    drift 3.291e-7   land resid 1.815e-14  corr 7460158
```

The residual falls **nine orders**, the correction fires **7.5 million times**, and the drift
moves by **0.03%**. Same on every case and stepper.

**That is not a negative result, it is the wrong instrument, and CLAUDE.md:1292 says so
already:** *"A DIAGNOSTIC FIELD IS SPECIFIC TO A CLASS OF DEFECT, AND ENERGY DRIFT IS BLIND TO
THIS ONE... The overshoot displaces the state in TIME and the AZ energy is nearly stationary
along the flow."* That entry is about the overshoot clamp, which is the same defect one revision
earlier, and it ends *"ask what the diagnostic would say about the defect you are hunting BEFORE
reading it as clean."*

The figure-eight measured **closure** -- a position and velocity error, which a displacement in
time changes directly. The field harness measured **energy drift**, which is nearly invariant
along the flow, so a state displaced in time carries almost the same energy. The two tests
disagree because they measure different things, and the landing correction improves the one the
field harness cannot see.

**So the port question is unanswered, not answered no.** What it needs is a trajectory
diagnostic on a field -- a shape chord against a much finer reference, the way `dtau_mode` pairs
were compared -- not `energy_drift_max`.

### AND GBS LOOKS LIKE IT BEATS HEGGIE, WHICH CONTRADICTS THE ANSWER GIVEN TWO MESSAGES EARLIER

Decades, GBS against the `logh_arms` table:

```
                       vs heggie      vs az   vs logh_rk4
    config_stability      +0.98      +2.65         +2.75
          near-field      +1.20      +1.84         +3.29
       deep_interior      -3.43      -3.27         -1.97
```

Two of three cases ahead of Heggie, and `nonfin` runs RK4 14 -> GBS 0 on `config_stability` and
140 -> 0 on `deep_interior`. Cost is 3-4x the evaluations.

**This is NOT a controlled comparison and must not be quoted as one.** The GBS rows run with the
predictive step limit **off** and the `logh_arms` rows run with it **on**, which is the knob
measured to be worth up to 13 decades on `far` and to be fatal at an exact collision. Different
resolutions too (96^2 against 256^2), which per-trajectory statistics tolerate but which should
be said. The controlled run is `logh_gbs` as an arm of `logh_arms` under the same limit as
everything else, and it is the obvious next thing.

"logH does not beat Heggie" was answered on the bare leapfrog and the RK4 arm. Under the
configuration Mikkola actually recommends it may well; that is now a live question rather than a
settled one.

### THE LANDING CORRECTION DOES TRANSFER — AND ITS VALUE ORDERS INVERSELY WITH CHAOS

Re-run against a **trajectory** diagnostic (shape chord against a converged GBS reference at
`eta/16`, with a second reference at `eta/32` as the convergence guard) instead of energy drift.
`gain = chord_off / chord_on`:

```
                          t=1     t=2     t=4     t=8    chord off @ t=8
    near-field       Rk4  0.997  41.118  41.180  22.454     2.458e-5
                     Kdk  1.000   0.855   0.889   0.979
    config_stability Rk4  1.000   6.893   1.966   1.587     2.190e-5
                     Kdk  1.000   1.002   1.006   1.002
    deep_interior    Rk4  2.396   1.168   0.993   0.983     1.417e0   <- SATURATED
                     Kdk  0.990   0.998   1.000   0.997
```

**Up to 41x, where `energy_drift_max` read 1.00.** The blind-instrument diagnosis is confirmed
from the other side: the same correction, the same runs, a diagnostic that can see a displacement
in time.

**And the ordering is the mechanism, not scatter.** `near-field` is the tamest region in the
corpus -- `d_min/R = 3.8e-3`, 2.4% of pixels collide -- and gains 22-41x. `config_stability` is
chaotic and collision-heavy: 1.6-6.9x. `deep_interior` is the most chaotic: **nothing** past
`t = 2`. The landing contributes a roughly fixed error per boundary, so it matters in proportion
to how little else is going wrong.

**`deep_interior` at `t = 8` is a dead row and says so on its face**: `chord off = 1.417` against
a maximum of 2.000 on the shape sphere. That is chaotic saturation, which the harness was
designed to expose rather than average over.

**The KDK control never improves, in twelve cells out of twelve** -- 0.855 to 1.006, and
*slightly worse* where it moves at all, paying for corrections a second-order method cannot use.
That is what licenses reading the RK4 column: a correction that improved every arm would not be
removing an order-two cap.

### WHAT THAT MEANS FOR THE PORT, WITH THE DECISION LEFT OPEN

For AZ and Heggie the correction would buy **a large trajectory improvement in tame regions,
nothing where chaos dominates, for 2-7% more force evaluations**. And it is not an academic
quantity here: the shipping science field is `spread_shape`, which is trajectory-derived, so the
rendered field in tame regions would move.

Against: it invalidates every committed number in `results/`, and the regions it helps most are
the ones already at `1e-11` drift where nothing was visibly wrong.

Both halves are now measured rather than argued. The decision is not taken here.

---

# Session: GBS sweep, the no-discard fix, the secant landing, and two overnight stages

## THE GBS VERDICT: ONE CLEAN WIN IN SIX, AND IT NEVER WINS A TAIL

Controlled sweep, 256^2, matched force-evaluation convention, same step control, both limit arms,
`refine_flagged` off. Gain is against Heggie, positive means GBS better:

```
  near-field         +1.02          win, err>10 = 0
  config_stability   +1.33 median   bulk win, 131 not-data, p10 -1.51
  preset_shape_h1    +0.66 median   bulk win, 128 not-data, p99 2.9 dec worse
  preset_shape       -0.60          loss, 222 not-data
  deep_interior      -1.44          loss, 1557 not-data
  far                -4.20          uniform loss, p10..p90 spread only 0.35 dec
```

At **~3x the force evaluations** throughout, and **Heggie reads `err>10 = 0` on all six**. Wherever
GBS does not win outright it produces a not-data population Heggie does not have. `far` is the
extreme: a 0.35-decade spread across 65536 pixels is a constant offset, not a bad tail.

**`logh_rk4` lost all six**, by 1.46 to 3.17 decades. That settles the pre-registered prediction 4
in the direction it was written: on shell `K + B == U`, so an integrator evaluating both
denominators at the same point sees only a Sundman transformation.

**Six reference arms reproduced the committed corpus bit-for-bit**, including `preset_shape_h1` at
`7.835e-11`, which is checkable against the figure recorded in CLAUDE.md rather than against a
file.

## THE STEP LIMIT FOLLOWS A RULE -- AND IT IS RESOLUTION-SENSITIVE, WHICH I FIRST OVERSTATED

At 256^2, six cases, no exceptions: the predictive limit is essential exactly where overshoots
occur without it, and a tax where they do not.

```
                     nolim overshoot events   what the limit does
  far                          0              7.6x WORSE, +73% evals
  near-field                   0              mildly worse, +23% evals
  preset_shape_h1             67              err>10   906 ->  128
  preset_shape               117              err>10  1608 ->  222
  config_stability           234              err>10   718 ->  131
  deep_interior              274              err>10 10514 -> 1557
```

**The caveat, which I stated as "no exceptions" before measuring it:** at **128^2** on
`config_stability` the limit *hurts* -- `err>10` 28 with it against 9 without -- where at 256^2 it
helps 718 -> 131. Same case, same statistic, opposite verdict at two resolutions. The rule holds at
the resolution that ships; the confidence does not transfer downward.

Not a contradiction and worth stating exactly: the sweep's `over` column **sums** `n_overshoot`
events across pixels while `gbs_tail` **counts pixels** with any. 234 events and 1 pixel are
different quantities.

## `gbs_unconverged` WAS COMPUTED ON EVERY MARCH AND READ BY NOTHING

Third instance of this pattern after `ab_floored` and `ab_min`. It stopped at `LhOut` and never
reached `MarchOut` or `PixelOut`, so no render, dump, criterion or test could see it. Now plumbed.

It fires on **83% of pixels** on `config_stability`, at a rate of **6.3e-4 per macro-step** -- so
GBS is not failing constantly, it is that a trajectory of ~5.3e5 macro-steps almost always
accumulates one. A per-pixel boolean over a long march saturates for the same reason
`n_cap_hits > 0` did.

**Which makes it useless as a discriminator, and the harness labels that rather than scoring it.**
It covers 1.0000 of the hot set at a lift of **1.203**; chance alone would cover 23.3 of 28.
`budget_exhausted` is the sharp one -- lift **146** -- but covers only 3.6% of the tail. Nothing in
either tail is unexplained, and no single site attributes it.

## `logh_arms`' `budget` COLUMN IS OVER TWICE THE DENOMINATOR OF EVERY OTHER COLUMN IN ITS ROW

Counted over `px.iter().chain(dpx.iter())` -- both passes -- while `err>10`, `drift`, `nonfin`,
`steps`, `evals` and `over` are over the diagnostic pass alone. A truncated run in the science pass
cannot contribute to `err>10` at all, so "35 of 131" was never a subset claim.

**And the obvious repair for reading it was wrong.** I argued `nonfin = 0` proves the diagnostic
pass carries no truncation, since `budget_exhausted` sets `finite = false`. It does not: the drift
reduction *filtered non-finite copies out*, so a dead copy left the pixel finite. Which is the next
finding.

## A FAILED COPY WAS BEING DISCARDED FROM `energy_drift_max`, AND IT FAILED CHAOS-SELECTIVELY

`if dr.is_finite() { drift_max = drift_max.max(dr) }` **dropped** a diverged or truncated copy and
left the pixel reporting a healthy-looking maximum over the survivors. That is the no-discard rule
broken, and it breaks in the worst direction: a copy goes non-finite because its integration was
hard, integration is hard at a close encounter, and close encounters are what the instrument
exists to measure.

**`finite` means different things in different occupants.** `az::driver` sets it false only when
the state itself goes non-finite and leaves it `true` for a run it truncated; `logh::driver` sets
it false for both. The first patch keyed on `finite` and `tests/no_discard.rs` failed **64 of 64**.
The reduction now tests `o.finite && !o.budget_exhausted`, which needs no knowledge of the driver.

`+inf` and not `NaN`, because these are max-reductions and infinity is the absorbing element.
**`d_min` is deliberately untouched**: it is a *min*-reduction with no absorbing element for
"undetermined", and `-inf` would read as a collision at every threshold and rewrite the outcome
labels.

The test carries the control arm -- the same pixels at a generous budget must come back finite --
because a reduction hard-wired to return `inf` passes the property arm exactly as well as a correct
one.

## THE SECANT LANDING IS PORTED TO AZ, HEGGIE AND THE NUMPY REFERENCE, AND xcheck IS 4/4

`clamp_final_step` sizes the landing step from `dt/dtau` read *before* it, a first-order predictor,
so the clock is corrected to the boundary and the state is not. That `O(h^2)` residual is what
holds the measured convergence order at 2.08 (AZ) and 2.40 (Heggie) rather than at RK4's four.

Porting it to Rust alone broke `short_horizons_match_to_1e13` and `long_horizon_t13` while
`algebra_matches_the_reference` stayed green -- the equations untouched, the trajectory changed.
**Fix both sides or the cross-check disagrees for the right reason**: `reference/tb_az.py` carries
`LAND_ITERATE_DEFAULT` and a vectorised secant that leaves unselected lanes bit-for-bit alone.

With both sides ported, **xcheck is 4/4** -- two independently written implementations of the
correction, Rust scalar and NumPy vectorised, agreeing to 1e-10 on a Burrau trajectory. That is a
stronger statement than either side passing alone.

`land_iterate` defaults **true in all three `*Opts` structs and in `EnsembleCfg`**, deliberately:
this session was bitten three separate times by two option structs disagreeing about one field
(`step_limit` between `AzOpts`/`HgOpts`, `finite` between drivers, and the stale hardcoded
`land_iterate: false` in `pixel.rs` whose justification had expired).

**A fixture moved under it for the fifth time.**
`the_relative_mask_desaturates_where_the_absolute_one_does_not` ran at `t_max = 2.0` -- the
anti-pattern this project has recorded three times -- and went vacuous, because the landing removes
a *spatially varying* first-order error at every boundary, so the hot mask became coherent and
stopped fragmenting (40.3% -> 0%). Repaired by pinning to `t = 13` **and** replacing the calibrated
`rel_multi > 0.2` with a **paired** comparison against the absolute mask, which is the claim the
test is actually making and needs no constant. Re-pinning the constant would have been a tolerance
loosened to keep green.

## OVERNIGHT STAGE 0: `far`'s AZ WIN INVERTS AT f32, AND MY MECHANISM WAS WRONG

Predictions were written to `results/overnight/PREDICTIONS.md` before the run. Saturation guard
reads **LIVE at both precisions** -- 16368 / 16347 distinct drift values, distribution overlap
0.0004 -- so the comparison has a subject.

```
      az  f64  drift p50  2.823e-13   err>10   0
  heggie  f64  drift p50  2.189e-12   err>10   0    AZ better 0.890 dec, 100% of pixels
      az  f32  drift p50   1.175e-4   err>10  25
  heggie  f32  drift p50   2.022e-5   err>10   0    HEGGIE better 0.768 dec, 100% of pixels
```

**I predicted the advantage would WIDEN and argued against the brief's rule. The brief was right.**
I reasoned that `Gamma*`'s degree-six algebra would lose the most digits at low precision. The
measurement says AZ is the precision-sensitive arm: f64 -> f32 costs AZ a factor of **4e8** against
Heggie's **9e6**, about 44x more. **No replacement mechanism is offered** -- the failed one was
mine and inventing a second in the same breath is the move the original brief ruled out.

Chain coordinates (Stage 2) are justified on the brief's rule, and are a round-off fix predicted to
show hardest at exactly the precision where the effect lives.

## OVERNIGHT STAGE 1: REVERSIBILITY IS STABLE AND IS NOT MEASURING CHAOS

March to `t_max`, negate velocities, march back, compare COM-centred positions. No reference
trajectory, one scalar per pixel, 2x cost.

**Two harness defects had to be fixed before any number meant anything.** The residual was first
computed on **absolute** positions: AZ and Heggie reconstruct from relative coordinates and place
the COM at the origin, logH integrates absolute positions and leaves it where the decode put it, so
the two families sat a constant translation apart -- `(-0.0125, +2.4875)` on `far`, identical for
all three bodies. Every other comparison in this project is translation-invariant, which is why it
had never mattered. Read raw it reported the COM offset as a flat, `eta`-independent error of
`4.016e-1`. And `AzOpts::default()`/`HgOpts::default()` disagreed on the step limit, so the arms ran
under different step control.

**PREDICTION 1 (reversibility tracks FTLE better than drift) is REFUTED across all four cases.**
Null everywhere against the shifted control, and most decisively on `deep_interior`, the case with
the widest FTLE range (2.29-4.98), where every arm reads -0.01 to -0.04 against controls of the
same size.

**PREDICTION 1b fired and killed the normalisation arm.** `rev/amp` correlates with FTLE at -0.92
to -0.999 -- not "much less", a near-deterministic anti-correlation. It is `exp(-lambda t)` wearing
another name. Fully characterised: the correlation is set by whether the amplification's variation
exceeds the residual's and nothing else, which is why it reads -0.98 on `near-field` (FTLE spans
1.22-2.35) and only -0.09 for `logh_rk4` on `far` (FTLE spans 0.01, but that arm's `rev` spread is
64%).

**The premise, not the arithmetic, is what fails.** `deep_interior` has `lambda t ~ 42.6`, so
round-off times `e^{lambda t}` predicts a residual of order **100** -- complete irreversibility.
Measured: **7.175e-6**, eight orders out. Reversibility is not round-off carried through the tangent
flow, and no claim about what it *is* is offered here.

**What survives:** the magnitudes are stable and resolution-independent (`near-field` az
`3.943e-10` at 64^2 against `3.950e-10` at 128^2), and it reports a failure mode energy drift
cannot see -- **GBS saturates rather than converging, and worsens under refinement**
(`7.556e-12 -> 1.801e-11` over a 16x cut). It made one successful prediction: it ranked
`logh_gbs` worst on `far`, saturated, **before** the sweep's `far` GBS row existed, and that row
came back at -4.20 decades.

**Two cautions against my own tables.** The AZ/Heggie reversibility ordering **flips by case and is
never at matched cost** -- AZ wins `far` and `deep_interior` at 1.7x and 2x Heggie's evaluations,
Heggie wins `config_stability` at 1.2x AZ's -- so it is not an ordering. And the half-frame shifted
control reads 0.32-0.34 on `near-field`, large enough that it may not be destroying the spatial
relationship there, in which case it licenses nothing on that case in either direction.

## TTL IS BUILT AND VALIDATED, AND IT LOSES — MONOTONICALLY WORSE WITH MASS RATIO

Mikkola & Aarseth, CMDA **84** (2002) 343. `Omega = sum w_ij / r_ij` replaces logH's mass-weighted
`U`, with `W` carried in the state and advanced by `dW = (dOmega/dt) dt`.

**The weight choice makes the control an identity.** `w_ij = mbar^2` gives `Omega === U` exactly at
equal masses, so the `q = 1` row is an algebraic identity rather than a near-agreement. `w_ij = 1`
would have been simpler and left `Omega` on a different scale from `U`, so a comparison at fixed
`eta` would have scored the step size instead of the transformation.

`W` is registered **once at `t = 0`** and carried across every boundary — it is the analogue of
`B`, and re-seeding it per interval would discard the off-shell information the transformation runs
on and make the march depend on `n_sync`. `LhState` went 13 -> 14 components; `to_array13` is now
`to_array14` so a length mismatch is a compile error rather than a silent index shift.

**Validation, `tests/logh_ttl.rs`, five tests each with a non-vacuity arm:**

- `Omega === U` at equal masses to 1e-15, paired with the arm that they differ >10% at 90:1.
- `omega_dot` finite-differenced against `omega`, with a sign-flip mutation arm asserted to fire.
  Step `1e-3`, not `1e-8`: `Omega` is smooth and `O(1)`, so a small step surrenders digits to
  cancellation for nothing.
- **The `W`-vs-`Omega` gap converges at SECOND ORDER: `1.500e-5 -> 3.748e-6 -> 9.368e-7`, ratios
  `4.00` and `4.00`.** The first cut of this test asserted an absolute `< 1e-6`, measured `3.7e-6`,
  and would have been "fixed" by halving the step. An absolute tolerance cannot separate
  second-order error working as designed from a wrong `dW`; a convergence ratio can — a wrong term
  does not converge and a first-order one converges at 2x.
- End to end: equal-mass `|dr| 4.113e-14`, 90:1 ratio `|dr| 2.498e-5`.

**The sweep, `examples/ttl_mass_ratio.rs`, near-field configurations at 48^2, KDK, masses varied
and nothing else:**

```
      q   drift p50 logh    drift p50 ttl   gain(logH/TTL)   TTL cost
      1        8.363e-6         8.363e-6          -0.000      1.000x
      2        1.124e-5         5.786e-6          +0.288      0.999x
      5        1.147e-5         1.268e-5          -0.044      0.997x
     20        7.862e-6         2.334e-5          -0.472      0.996x
    100        2.683e-6         1.854e-5          -0.840      0.997x
   1000        2.369e-7         1.671e-5          -1.848      0.997x
```

Control exact at `q = 1`. `err>10` and `nonfin` **zero on every arm at every rung**, so neither
side is dead. Cost at parity throughout, so this is not a quality-for-work trade.

**The prediction — TTL beats logH at high ratio, ties at equal mass — is REFUTED on its first
arm.** No mechanism is offered.

The fact worth carrying: **the gap widens because logH IMPROVES 35x across the ladder**
(`8.4e-6 -> 2.4e-7`) while TTL stays flat near `1.7e-5`. Untested and stated as such: at `q = 1000`
the masses are `(0.4998, 0.4998, 0.0005)`, close to a two-body problem plus a test particle, and
whether that is still a mass-ratio sweep or a near-integrable limit is a separate measurement.

**IAS15 is not built.** Gauss-Radau nodes, a predictor-corrector loop and its own step control; its
only role here is a reference arm, since its per-lane variable work is already measured as fatal on
GPU (`warps hit 1.0000`).

## STAGE 2: CHAIN COORDINATES BUY 3.6x AT f32 AND NOTHING AT f64 — THE PREDICTION HOLDS

Mikkola & Aarseth, CMDA **57** (1993) 439. `src/integrate/logh/chain.rs`. For three bodies the
chain is two vectors: `X1 = r_b - r_a`, `X2 = r_c - r_b`, and the third separation is `|X1 + X2|`
— a **sum of small quantities** rather than a **difference of large ones**.

**The ordering is FROZEN, selected once at registration.** Choosing a chain ordering is a
*re-registration*, the same class of act as AZ picking a reference body, which is the mechanism
this whole logH investigation exists to isolate. A chain that re-selects at sync boundaries
reintroduces exactly what logH was built to have none of. The re-selecting variant is a named arm
and never the default.

It is a **diagnostic** integrator: no events, no closure escape, no `t_end`, no outcome labels, so
its numbers are **not comparable with the committed corpus**. It is deliberately too small to be
mistaken for a gallery arm.

**Result, `far` at 64^2, 4000 fixed fictitious steps of 1e-3, LogH + RK4, same ICs both arms:**

```
  prec   coords    drift p50    drift p90    drift p99   nonfin
   f64   direct    1.036e-15    2.691e-15    4.145e-15        0
   f64    chain    1.451e-15    3.728e-15    5.796e-15        0
   f32   direct     2.781e-6     3.113e-6     3.336e-6        0
   f32    chain     7.788e-7     1.223e-6     1.668e-6        0

  gain = log10(direct/chain):   f64 -0.146      f32 +0.553
```

**PREDICTION 2 — "chain helps at f32 and is near-invisible at f64" — CONFIRMED**, and it was
written before Stage 0 ran. The f64 control is what licenses the f32 reading: both f64 arms sit at
the ~1e-15 round-off floor, so `-0.146` is noise between two floors and not a cost.

### Four defects found building it, three by tests that carried a negative control

- **A SIGN ERROR IN BOTH RELATIVE ACCELERATIONS.** `a1` and `a2` had the `g1` term inverted. The
  march "worked" and produced trajectories; energy drift read **77.5**. Caught by differencing
  against `newton::accel` with a **swapped-pair negative control** — an index assertion alone would
  have passed on a transposition, the standing `shape_pl` lesson at a new site.
- **THE MARCH FIXTURE WAS UNRESOLVABLE, NOT WRONG.** The tight pair sits at `2.2e-3`, orbital
  period `~2.5e-4`, so an *unregularised* step of `1e-3` is four periods long. That is what the
  time transformation exists for, and the first cut of the test used `LhTime::None`.
- **A NON-COM-CENTRED FIXTURE GAVE A FLAT DRIFT CURVE — `1.890e-6` at three step sizes.** Read
  literally that is *the* wrong-equation signature. It was the frame: `to_cart` returns a
  COM-centred configuration because a chain has no COM degree of freedom, while `e0` was computed
  on the uncentred input, so the two energies differed by the COM kinetic energy — a constant,
  independent of step size. *A construction that assumes a COM-centred input returns a drifting
  system without one*, now at a third site, and the flat curve is what diagnosed it.
- **THE SEPARATION TEST PASSES VACUOUSLY AND THE FILE SAYS SO.** Chain and direct give
  **identical** f32 separations (`1.213e-7` both), because `from_cart` performs the same
  subtraction once. A one-shot conversion is *structurally incapable* of showing the benefit — it
  is dynamical, accruing over a march as the vectors are carried instead of reconstructed. The
  test now asserts they are identical, which is the honest claim; asserting an improvement there
  would have been a test that cannot fail in the direction that matters.

**And the convergence check replaced an absolute bound, twice.** The march drift reads
`2.068e-8 -> 1.275e-9 -> 8.358e-11`, ratios **16.2 and 15.3** — fourth order, twice. An absolute
threshold cannot separate "fourth-order truncation working as designed" from "the equations are
wrong", and this file had already caught a sign error a magnitude test would have passed. The
rungs are pinned **above** the fixture's own `~4e-12` round-off floor: below it the same ladder
reads ratios 5.8 then 1.4, which is the floor and not a failure.

## STAGE 4: IAS15 IS BUILT — MACHINE PRECISION, AND TWO INDEX BUGS ON THE WAY

Rein & Spiegel, MNRAS **446** (2015) 1424. `src/integrate/ias15.rs`. Fifteenth-order Gauss-Radau,
adaptive, error following Brouwer's law.

**It is a REFERENCE ARM, not a production candidate**, and that verdict is measured rather than
assumed: its predictor-corrector iterates a variable number of times per step — measured **min 5,
max 12** over 200 steps — which is the per-lane variable work already recorded as fatal on GPU
(`warps hit 1.0000`, worst lane 5.2 million retries). The variability is *asserted*, because it is
the property and not a defect to hide.

**Why it is worth having:** this project has no reference. `eta/256` came back **saturated** —
chord 2.000, antipodal, scoring a correct mode and a broken one alike.

```
  fixed step dt = 0.1                     drift 3.469e-16
  adaptive, t = 200, 5697 steps           drift 5.204e-15   (244432 evals, mean 5.99 iters)
```

### The conversion matrix is COMPUTED, not transcribed — and that decision paid twice

The `g -> b` conversion is 21 magic constants in the reference implementations. Two sign errors in
this project's AZ algebra were invisible until someone finite-differenced the Hamiltonian. So the
matrix is built by expanding the Newton basis into the monomial basis, and the test checks it
against the **defining identity** at random points with a perturbed-matrix negative control — a
transcribed table can only be tested against a copy of itself.

**The first cut was off by one and it still integrated.** `g_k` multiplies
`prod_{i=0..=k}(t - h_i)` and `h_0 = 0`, so every basis polynomial carries a factor of `t`: the
interpolant has no constant term and IAS15's `b_k` is the coefficient of `h^(k+1)`, not `h^k`.
Recording the product *before* multiplying, and indexing from `t^0`, shifts the matrix one place.
It compiled, it ran, it produced trajectories — **energy drift 7.8e-1 at `t = 200`** against the
`1e-15` the method is for, and a halving that bought **3.6x** where fifteenth order gives ~32000x.
Exactly the "wrong algebra looks like physics" failure mode, in new code, caught by an accuracy
test rather than by inspection.

**And the second bug was in the harness, not the method.** With the matrix fixed, the *fixed-step*
drift was already `3.469e-16` while the *adaptive* run read `1.282e-7`. The step controller was
being handed a **hardcoded constant** where `max |b_6|` belongs, so it degraded silently to
"multiply `dt` by a fixed factor every step" and looked like an integrator fault. `step` now
returns `b_6` rather than letting a caller supply it.

`next_dt` clamps its growth factor to `[0.2, 4.0]`: an unbounded step control is how this project
recorded a single step advancing the clock by `2.209e128`.
