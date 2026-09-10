# The comparison-only rule, audited against the kernel that exists

**What was asked:** *verify that rule is actually followed in the shipped kernel, not just
specified.* This is that audit. Every claim carries `file:line` and every line was read back
directly, not taken from a summary.

**Read this as distance-to-target, not as a defect list.** `prin-rs` is the CPU research
instrument. It was written without the GPU constraint and several of the findings below are
*correct* choices for that job — `d < r_coll` is the natural way to write a collision test, and
`T::infinity()` is idiomatic Rust for "no bound in force". The finding is the **size of the gap**,
and it is large.

---

## 0. The rule, and which kernel it is about

From `principia_gpu_determinism_note.md`:

> **Anything a control-flow decision depends on is frozen at build time or reduced to a comparison —
> never left hostage to a floating-point operation the backend is free to compute its own way.**

**The rule is followed in the contract and in no kernel that exists.** `prin-rs` has **zero GPU
code** — dependencies are `num-traits`, `rayon`, `png`, `plotters`, `gif`, `color_quant`, and there
are no `.wgsl`, `.spv`, `.metal` or `.glsl` files anywhere in the tree.

And the CPU kernel is a **different integrator family** from the one the rule was written for. The
contract specifies `substep_bucket` / `N_sub` / `N_max` — a frozen-threshold bucket lookup. prin-rs
uses AZ/Heggie regularisation with `dtau` and `StepLimit::Predictive`. The repo says so itself, at
`src/integrate/mod.rs:63-66`:

> `principia_integrator_contract.md` is the **GLSL app's** contract and is not in this repo; its
> `substep_bucket`/`N_sub`/`N_max`/descriptor bit 5 appear nowhere in `src/`.

So "is the rule followed in the shipped kernel" resolves to: **the shipped kernel is not the kernel
the rule is about**, and the mechanism the contract relies on to make Tier B achievable is not
present to be audited. What follows audits the kernel that is here.

---

## 1. The finding that outranks all six rules

**`src/real.rs:60-83` gives f32 and f64 *different* branch thresholds.**

| constant | f64 | f32 | where it decides a branch |
|---|---|---|---|
| `SYNC_EPS` | `1e-15` | `1e-6` | the sync-skip test, `src/integrate/az/driver.rs:964` |
| `LAND_EPS_REL` | `1e-14` | `1e-5` | the landing test, `heggie:345`, `az:1005`, `logh:356` |
| `TINY` | `1e-300` | `1e-37` | the degeneracy floor, `heggie:354` |

The module's own header says it, at `src/real.rs:13-14`:

> This is the one place where f32 and f64 are not running the same algorithm, and it lives exactly
> where degenerate geometry lives.

**Consequence: the parity contract's Tier B — branch words bit-exact between CPU-f64 and GPU-f32 —
is unachievable for this kernel by construction, before any backend enters.** This is not backend
latitude and not a compiler's freedom; it is a deliberate, documented, precision-dependent threshold
sitting directly in a control-flow decision. Any determinism contract over this kernel has to say
**which precision is normative** before R1 can even be stated.

The choice was right for what it was for — the f64 arm carries the reference's literal values so the
NumPy cross-check is a clean equality, and `1e-15` is below f32 epsilon so the f32 sync test would
degenerate to `t < t_target`. It is simply incompatible with cross-precision bit-exactness.

---

## 2. The six rules

| rule | status | evidence |
|---|---|---|
| **R1** comparison-only | **violated, severely** | below |
| **R2 / R6** explicit `fma` | **violated** | `grep -rn "mul_add\|fma(" src/` → **0 hits, crate-wide** |
| **R3** clamp before float→int | **N/A in drivers, one adjacent** | `src/physics/ftle.rs:154` |
| **R4** no non-finite literals | **violated, in-kernel** | **27** `T::infinity()` in `src/integrate/` |
| **R5** flag-tested `while` | **violated** | **19 / 9 / 7** `break`s in az / heggie / logh |

### R1 — transcendentals and divisions deciding control flow

**The collision test is sqrt-decided in every integrator.** `Vec2::norm` is
`self.norm_sq().sqrt()` (`src/vec2.rs:42`), and:

- `src/integrate/heggie/driver.rs:431` — `if qi.norm() < r_coll {`
- `src/outcome.rs:499` — `if d < r_coll {`
- `src/outcome.rs:179` — `if *dk < r_coll {`

There are **46** `.norm()` call sites across `src/integrate/`, `src/outcome.rs` and
`src/physics/newton.rs`. The squared form is available for free — `norm_sq` *is* the sqrt's input.

Worse, **the threshold is not a constant either**: `r_coll = r_coll_frac * r0` where `r0` is a
runtime `(inertia(r, m) / mtot).sqrt()` (`src/physics/energy.rs:78`). Squaring both sides fixes the
left; the right still needs every backend to agree on one sqrt of a quotient — though only once per
trajectory rather than per step.

**Every step-limit comparison consumes a division.** `heggie/driver.rs:248` is the predictive limit
`f * d_min / denom`; `:359` is `(eta * dt_left / d).min(dtau_entry)`; `:367` is `if lim < dtau`.
AZ's `SoftMin` blend branches on a **`powf`** result (`az/driver.rs:645-647`). LogH branches on a
**ratio of two divisions** with a near-one fudge factor (`logh/driver.rs:391`) — the most
backend-fragile predicate in the kernel.

**The escape test uses all three at once** (`src/outcome.rs:311`):

```rust
dr.norm() > r_esc && dr.dot(dv) > T::zero() && spec > T::zero()
```

a sqrt comparison, a contractable dot product, and a value produced by a division by a sqrt'd
floored distance (`:295`).

**The landing test is exactly as documented** — `az/driver.rs:1005`,
`if s.t >= dt_left - land_tol::<T>(...) || bad` — comparing two independently accumulated floats
against a threshold that is itself a runtime multiply.

**One thing the kernel gets right and should be preserved.** The `TINY` floor tests are pure
`norm_sq` against a compile-time constant — `heggie/driver.rs:354`,
`if s.r(0) < T::TINY || ...`, where `s.r(i)` is `self.u[i].norm_sq()`, no sqrt. **That is the shape
every other test should have**, and it is already in the file.

### R2 / R6 — no `fma` anywhere

`grep -rn "mul_add\|fma(" src/` returns **zero hits in the entire crate**. Every squared distance
and dot product is `self.x * o.x + self.y * o.y` (`src/vec2.rs:26`) — the exact contractable
`a*b + c` the rule forbids — and it feeds every collision, escape and step-limit comparison above.

### R4 — `+inf` and `NaN` are load-bearing, not incidental

`T::infinity()` is the kernel's idiom for "no bound in force" and seeds min-folds rather than the
first element (`heggie:244,415`; `az:630…1025`; `logh:253…441`; 27 sites in `src/integrate/`).
`T::nan()` is a **deliberate sentinel** (`heggie:498`) relied on at `outcome.rs:361`.

This collides directly with the reduction semantics: `+inf` is the absorbing element that makes the
no-discard fix correct (`src/ensemble/pixel.rs:1077-1080`), while naga **rejects infinity literals**
outright. Reseeding the folds is not mechanical — it changes what an all-undetermined fold returns,
which is the exact question the no-discard work settled.

### R5 — every march loop is a bare `loop {`

`az/driver.rs:997`, `heggie/driver.rs:338`, `logh/driver.rs:349` are all literally `loop {`, with
four to six exits each, including labelled `break 'outer` unwinding two levels
(`heggie:343,351`). AZ adds a nested value-carrying retry `loop` (`az:1047-1064`). And the landing
secant is the named anti-pattern verbatim — `for _ in 0..opts.land_max_iters {` with a mid-body exit
(`heggie:388`, `az:1084`, `logh:400`), a data-dependent trip count whose termination is decided by a
division result.

This is the shape that produced the silent 1-of-3-body-pairs truncation through naga→MSL.

---

## 3. The good news, and it is substantial

**The step/deriv layer is already shader-shaped.** Measured across `az/rk4.rs`, `heggie/rk4.rs`,
both hamiltonians, `physics/newton.rs` and `physics/energy.rs`: **0 `Vec`, 0 `Box`, 0 `dyn`, 0
`std::`** in every one. States are already fixed arrays (`to_array9`, `to_array13`, `to_array14`).
There is **no rayon and no dynamic dispatch anywhere in `src/integrate/`**.

**Every blocker is in the driver layer, and every one is telemetry rather than physics** — 13 `Vec`
fields in `az/driver.rs` and 7 in `heggie/driver.rs` (reference history, boundary shapes, drift
history, escape flags, closure history), a `VecDeque` closure window, per-boundary `push`es, a
`sort_by`. A no-`std` fixed-array lowering of the step path is close; the driver needs its history
moved out of band.

---

## 4. Two structural notes, and one ordinary bug

**The decode is f64-pinned.** `src/physics/decoder.rs` contains **no `T: Real` at all** — the whole
z → initial-condition path is `f64`, as is `src/decode.rs`. So "generic over f32/f64 from one
source" is true of *integration* and not of *decode* — and decode is where the deep-zoom
linearisation lives.

**The 512-wide dispatch is legal under the shipped default and illegal under AZ + `RefPolicy::Shared`.**
Copies 1..7 consume copy 0's reference-body choices — but `src/ensemble/pixel.rs:808-810` states:

> Under Heggie `refs` is empty, so `Shared` hands over an empty slice and the per-copy path is taken
> anyway — the same behaviour, reached without a special case.

Heggie has been the default since 2026-09-02, so the copies are independent as shipped. Under AZ the
dependency is real and a 512-wide dispatch would need copy 0 to complete first.

**And an ordinary CPU bug found on the way.** `src/outcome.rs:289`:

```rust
let o: Vec<usize> = (0..3).filter(|&k| k != b).collect();
```

A heap allocation of two `usize`s, **per escape test — per boundary, per body, per copy**. A fixed
`[usize; 2]` is a drop-in replacement. This is a performance defect on the CPU path today,
independent of any GPU question.

---

## 5. Ranked fixes — recorded, not applied

**Nothing in `prin-rs` was modified.** These are ordered by value per unit of risk, for whenever the
GPU-kernel decision is taken.

| # | fix | invalidates the corpus? |
|---|---|---|
| 1 | `outcome.rs:289` — `Vec` → `[usize; 2]` | **no** — behaviour-preserving, and a CPU speed win |
| 2 | reseed min-folds with the first element instead of `T::infinity()` | **no** for a non-empty fold; **yes** in semantics for an empty one — see the no-discard note above |
| 3 | square the collision test — `d² < r_coll²` | **yes** — not bitwise identical at the boundary; flips labels on boundary pixels |
| 4 | explicit `fma` in `Vec2::dot` / `norm_sq` | **yes** — changes every downstream value |
| 5 | flag-tested `while`, zero `break` | no numerical change, large diff |
| 6 | make `SYNC_EPS` / `LAND_EPS_REL` precision-independent | **yes**, and it breaks the NumPy cross-check's clean equality — the reason they differ |

Fix 3 is *free* in arithmetic (`norm_sq` is already computed) and removes three sqrts per step, so
it is a performance win as well as a determinism one. It is listed third rather than first only
because it is corpus-invalidating.

**The decision this evidence should inform, and which is not taken here:** whether the GPU kernel is
a *port of prin-rs* or a *fresh kernel written to the rules*. If the latter, fixes 3-6 are not worth
making at all — prin-rs stays the f64 research instrument it is, and the rules are honoured in the
new kernel where they are cheap rather than retrofitted where they are expensive.
