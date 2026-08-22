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

**`error_ratio` uses MAD internally, aggregates by max.**
Internally: MAD (`1.4826 × median|x − median|`), not a standard deviation. A std returns NaN on
precisely the pathological pixel the statistic exists to flag. Aggregate over pixels by **max**,
not median — max tracks damage at Spearman +0.956 against +0.599 for median. Treat the result as a
**boolean flag**; its magnitude is unstable.

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
