# Overnight run — predictions, recorded BEFORE any number came back

Written at the start of the run. A prediction recorded after the numbers is not a prediction.

---

## Stage 0 — `far` at f32 vs f64, AZ and Heggie

**I disagree with the brief's decision rule, and I am recording that before the run rather than
after.**

The brief says: *"If AZ's `far` win COLLAPSES at f32 -> conditioning confirmed."* I think the
implication runs the other way, and the reasoning is short enough to check.

The conditioning story is that Heggie's `Gamma*` is **degree six** in the coordinates while AZ's
`Gamma` is linear in `A` and `B`. On `far`, body positions reach ~13 units, so `Gamma*` forms
intermediate quantities of order `13^6 ~ 4.8e6` and cancels them down to an `O(1)` result. That
costs roughly `log10(4.8e6) ~ 6.7` decimal digits of whatever the format has.

- f64 has ~15.95 digits. Losing 6.7 leaves ~9.2 — degraded but working. Measured: Heggie loses
  `far` by a flat 0.7-0.9 decades, on 100% of pixels.
- f32 has ~7.22 digits. Losing 6.7 leaves ~0.5. **There is nothing left.**

So if conditioning is the mechanism, **AZ's win should WIDEN sharply at f32, not collapse.**
Collapse would mean the penalty is insensitive to how many digits are available, which is close
to saying it is not a precision effect at all.

**PREDICTION 0a: AZ's `far` advantage over Heggie grows at f32 — the gap exceeds the f64 gap of
0.7-0.9 decades.**

**PREDICTION 0b, and this is the one I expect to bite: the test may be SATURATED and say
nothing.** `far`'s f64 drifts are `2.8e-13` (AZ) and `2.2e-12` (Heggie). f32's documented median
drift on this project is `9.3e-6`, which is *seven orders above both*. If f32 round-off dominates
for both arms they will land on a common floor and the comparison will be a difference between
two meaningless numbers — the failure already on record for the `plain_*` arms, where I withdrew
a "three decades" claim for exactly this reason.

So the harness prints a **saturation guard before the comparison**: the count of distinct drift
values per arm, and whether the two f32 distributions overlap. *A difference can be small because
both sides are right or because both are dead*, and this run is built to be able to tell me which.

**If 0b fires, Stage 0 does not settle Stage 2** and I will say so rather than reading the ratio.

---

## Stage 1 — reversibility

**PREDICTION 1: reversibility correlates with FTLE more strongly than drift does.**

Drift-vs-FTLE currently reads `+0.3048` for the leapfrog against a shifted control of `+0.0240`,
and `-0.0820` for AZ against a shifted control of `-0.1022` — which is a null, not a negative.
Reversibility on a chaotic field *is* round-off integrated through the tangent flow, so it should
track the tangent-flow exponent directly.

**PREDICTION 1b: the amplification-normalised form `|dx| / e^(lambda t)` should correlate with
FTLE MUCH LESS than the raw form.** That is the point of computing it: if dividing out the
amplification removes the FTLE correlation, the raw number was measuring chaos, and the normalised
one is the integrator-quality residual. If the normalised form still tracks FTLE strongly, the
normalisation is not doing what it is meant to and I should not report it as integrator quality.

**Known confound, stated up front:** `src/physics/ftle.rs` computes the FTLE **on the plain
leapfrog**, so the leapfrog arm shares an integrator with the field it is being correlated against
and every other arm does not. The shifted control is on every row.

**Reversibility is NOT a clean integrator-error measure** and the harness header says so: stepper
time-symmetry and adaptive step-control asymmetry both enter. RK4 is not time-symmetric, so
comparing reversibility *across steppers* scores symmetry rather than accuracy. Read down a
column.

---

## Stage 2 — chain coordinates

**PREDICTION 2: chain helps at f32 and is near-invisible at f64**, because it is a round-off fix
and f64 has digits to spare on this problem.

**Conditional on Stage 0.** If 0a holds, `far` is precision-limited and chain is the published
repair for exactly that. If Stage 0 comes back saturated (0b), this stage runs but its result on
`far` cannot be attributed and the file will say so.

---

## Stage 3 — TTL

**PREDICTION 3: TTL beats logH at high mass ratio and ties at equal mass.**

The equal-mass tie is the **control**, not a throwaway: a TTL arm that differs on an equal-mass
slice is measuring something other than the mass ratio, and the result would be uninterpretable.
Mass ratio is swept deliberately through `z8`/`z9`; the default charts are equal-mass and would
exercise nothing.

---

## Stage 4 — IAS15

**Not a production candidate, and its GPU verdict is already measured** — per-lane variable work,
the same property that gave reject-and-retry `warps hit 1.0000`. It is here as a
**near-machine-precision reference arm**, which this project does not have: `eta/256` came back
saturated at chord 2.000 for a correct mode and a broken one alike.

**PREDICTION 4: IAS15 provides a reference that `eta/256` could not** — i.e. it separates arms
that `eta/256` scored identically.

---

## Discipline

Force **evaluations**, not steps, in every table. Shared colour ramp constants. Provenance line on
every panel. Each stage writes on completion; a later crash loses nothing earlier.
