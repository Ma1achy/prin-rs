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
