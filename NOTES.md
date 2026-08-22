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

## 3. Smaller observations

- `reference/tb_az.py` uses `eta=0.02` as its default; `tb_all_az` and the smoke test pass
  `eta=0.01`, which is what BRIEF.md §2.3 specifies. The bare default is the odd one out.
- The reference's `dtau` sizing (`eta*dt_left/(A0*B0)`) fixes the *first* physical step as a
  fraction of the remaining interval, so effective resolution depends on `n_sync` through
  `dt_left`. Worth pinning `n_sync` explicitly in the cross-check test rather than leaving it at a
  default, or the `1e-10` agreement target is comparing two different discretisations.
