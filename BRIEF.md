# Build brief: Principia uniform-resolution kernel (Rust)

**Audience:** an agent with **no prior context on this project**. Everything needed is here.

**Deliverable:** a Rust binary that renders a uniform-resolution slice of three-body initial-condition
space to an image plus a raw dump, with no adaptive refinement and no interaction. It exists to make
measurements that a NumPy harness cannot, and to be the real kernel later.

**Scope discipline — read this first.** No quadtree, no scheduler, no GUI, no streaming, no
interaction. Uniform grid, one pass, write files, exit. Every one of those omissions is deliberate;
adding any of them is out of scope.

---

## 1. What the program does

Sample a 2D slice of initial conditions on a uniform `W × H` grid. Each grid point is one
three-body simulation. Integrate every one forward to a playhead `t`, classify what happened, colour
it, write a PNG and a raw binary dump.

**The physics is the whole product.** The image is a diagnostic, not the point.

```
for each pixel (i,j):
    ic        = decode(slice_params, i, j)          # grid position -> initial conditions
    copies    = ic + E jittered perturbations       # the "ensemble"
    for each copy: integrate to t_max, record fields
    reduce copies -> per-pixel summary
write PNG + raw dump
```

---

## 2. Physics

### 2.1 System

Planar three-body, Newtonian, `G = 1`. Three point masses in a plane; each has position `(x,y)` and
velocity `(vx,vy)`. Total state is 12 numbers.

```
a_i = sum over j != i of  G m_j (r_j - r_i) / |r_j - r_i|^3
```

**Default configuration — Burrau (Pythagorean):** `m = (3, 4, 5)` at rest at the vertices of a 3-4-5
triangle. Released from rest, so initial velocities are zero and total angular momentum `L_z = 0`.

**Reference constants for this configuration** (use as unit tests): total mass `M = 12`, mass-weighted
hyperradius `R = sqrt(I/M) = 2.2361`, total energy `E = -12.8167`, crossing time
`sqrt(R^3/M) = 0.9652`, measured Lyapunov exponent `lambda ~ 0.7`. **One time unit is roughly one
crossing time.**

### 2.2 The slice

The IC slice varies **one body's position** over a 2D box; the other two stay fixed. Parameters:
which body (`0`, `1` or `2`), box centre `(cx, cy)`, box half-width.

Useful regions, all Burrau, format `(name, cx, cy, body)`:

```
near-field    ( 1.0,   3.0, 0)     mid-field    ( 1.0,   6.0, 0)
far           ( 1.0,  13.0, 0)     body2 core   ( 1.0,  -1.0, 2)
body2 mid     ( 1.0,  -5.0, 2)     body1 slice  (-2.0,  -1.0, 1)
body1 far     (-2.0,  -7.0, 1)
deep interior ( 0.0,   0.0, 0)     <- pathological, see §2.6
```

### 2.3 Integration — this is the hard part

Naive integration **fails**: gravity diverges as two bodies approach, and the three-body problem has
close encounters constantly. Two mechanisms, both required.

**(a) Adaptive timestep, scale-covariant.**

```
dt = eta * min over pairs of ( r_ij^{3/2} / sqrt(G (m_i + m_j)) )
```

This is the local two-body free-fall time. Under a rescaling `r -> alpha·r` it scales as
`alpha^{3/2}`, exactly as time does — so it introduces **no fixed length or time scale**. That matters:
see §2.5. Use `eta ~ 0.01`.

**(b) Levi-Civita / Aarseth–Zare regularisation.**

Adaptive stepping alone is not enough — a genuine close approach drives `dt` toward zero and the step
budget is exhausted. **Regularisation is a coordinate change that removes the singularity.**

*Levi-Civita, the single-pair case.* For a pair with separation `rho` (as a complex number), substitute
`rho = u^2` and change the time variable `dt = |rho| dtau`. In the transformed Hamiltonian the singular
`-G m_i m_j / |rho|` term becomes a **constant**, so nothing blows up. The regularised two-body problem
is a harmonic oscillator, and a trajectory passes through *exact* collision at machine precision
(verified in the reference implementation: `d_min = 1.35e-11` with energy drift `6.2e-15`).

*Aarseth–Zare, two pairs at once.* Single-pair LC only fixes its own pair; a close approach of either
other pair still destroys the integration. AZ regularises **two** pairs simultaneously. Choose the
reference body `a` to be the one **not in the longest side** of the triangle, so both regularised pairs
share it:

```
R1 = r_b - r_a     R2 = r_c - r_a     R3 = R2 - R1   (unregularised)
```

Because `(b,c)` is the longest side, `|R3| >= max(|R1|,|R2|)`, so `R3 -> 0` only when all three
separations vanish — a genuine triple collision, which is provably non-regularisable anyway.

```
dt = |R1| |R2| dtau        (vanishes when EITHER pair closes)
```

**A Python reference implementation of AZ exists and is correct** — see §6. Port it; do not re-derive
it. The algebra is error-prone and its failures are silent (see §5).

### 2.4 Termination

A simulation stops when any of these fires. Record which.

| outcome | condition |
|---|---|
| **escape** | one body unbound from the other two (relative energy positive, separation growing) |
| **collision** | a pair's separation falls below `r_coll` |
| **triple collision** | **two or more pairs** below `r_coll` simultaneously |
| **triple ejection** | all three mutually unbound (requires total `E > 0`) |
| **running** | still bound at `t_max` |

**The ≥2-pair rule matters.** By the triangle inequality, if `|AB| < r_coll` and `|AC| < r_coll` then
`|BC| < 2 r_coll` — so "exactly two pairs below threshold" is reachable and is already a near-triple.
Requiring all three would silently misclassify it as an ordinary binary collision.

**Encoding.** `state` is 3 bits `{escape, bounded, collision, running, sim_failed, decode_failed}` and
`detail` is 2 bits.

*The table above gives conditions for five of the six states.* `running` is written as "still bound
at `t_max`" and `bounded` is given no condition, so read literally one of the six is unreachable.
Implemented instead as: **`bounded`** = reached `t_max` with nothing having fired, **`running`** =
did *not* reach it, the step budget ran out and the final state is not a terminal answer. That keeps
all six reachable and keeps "integrated to the horizon" distinguishable from "stopped early". For `escape`, detail is the escaping body `0–2`; for `collision`, the colliding
pair `0–2`. **`detail = 3` means "all three"** — triple collision or triple ejection respectively.
One rule, both arms.

### 2.5 Length scales — canonical units, fixed at `t=0`

`r_coll` and an optional softening `epsilon` must be expressed **in canonical units** (as a fraction of
the initial hyperradius `R`), **evaluated once at `t=0` and never updated**.

This is not stylistic. Newtonian gravity has a **scale invariance** (`r -> alpha·r`, `t ->
alpha^{3/2}·t` leaves the dynamics unchanged), which this project deliberately quotients out. A length
in *absolute* units breaks it — measured: the same physical system gave answers differing by 1.66×
purely from an arbitrary choice of overall size. In canonical units it is preserved to 10 decimal
places across a 16× rescaling.

**Fixed at `t=0`, never co-moving.** Scaling `epsilon` with the *instantaneous* system size is
tempting and **catastrophic**: it makes the Hamiltonian time-dependent and destroys energy
conservation. Measured `|dE/E| = 3.06e-02`, *identical* at `dt=1e-4` and `dt=2e-5` — insensitive to
step size, which is the signature of a wrong equation rather than an accuracy problem.

**Defaults:** `r_coll` nonzero (a small fraction of `R`), `epsilon = 0`. Softening is a *different
force law*, so tag it in the output and never mix `eps>0` with `eps=0` data.

### 2.6 Known-pathological region

`deep interior` (centre `(0,0)`, body 0) was expected to drive all three bodies together — a
near-triple collision, not regularisable, failing however well the integrator is built, at 190 s
per probe.

**Measured under Aarseth–Zare, that is not what happens.** It is an ordinary binary encounter
between bodies 0 and 2: `d_min = 2.28e-5` (Rust) against `2.30e-5` (numpy), `|dE/E| ~ 1.4e-7`,
two reference switches, reaching `t = 13` in about a second in both implementations. Sweeping
`r_coll` from `1e-4 R` to `R`, pairs (0,1) and (1,2) never register at any threshold.

The 190 s failure is the **unregularised** integrator. A close binary approach with a distant
third body is exactly the case AZ regularises — the warning predates the method that removes it.
See `examples/deep_interior.rs` and NOTES §2.4.

Expect a **binary collision**, not a triple. A genuine triple still exists in principle and is
still non-regularisable; this pixel is not it.

---

## 3. The ensemble

Each pixel carries **`E` extra copies** (use `E = 7`, so 8 total), each with its initial condition
jittered by a small random offset. Their disagreement measures whether the pixel's value is
well-defined at all.

**The jitter must scale with the cell size** — `jitter = jitter_frac × cell_width`, `jitter_frac ~ 0.5`.
A *fixed* perturbation would make measured spreads drift with resolution for a purely trivial reason.

**Every pixel always carries exactly `E+1` copies. Never discard one.** If a copy integrates badly it
is still a *measurement outcome* — "this could not be determined" — which is the strongest available
statement that the pixel is undetermined. Removing it biases the sample toward the tame trajectories,
which on a chaos instrument is exactly backwards.

**Open design question worth settling while building:** AZ's reference body is chosen per trajectory
and can change mid-run. Should all `E+1` copies of a pixel be **forced to share** the nominal copy's
reference? Unshared references may make copies accumulate differently-*structured* error, which would
corrupt the ensemble spread while leaving each trajectory's own energy drift healthy. **Implement it
as a flag and measure both.** See §4 experiment 2.

---

## 4. Per-pixel outputs, and why each exists

Write **all** of these to the raw dump.

| field | definition |
|---|---|
| `state`, `detail` | outcome, per §2.4 |
| `t_end` | time of termination (censor at `t_max`) |
| `d_min` | closest approach over the whole trajectory |
| `energy_drift` | `\|E(t) - E(0)\| / \|E(0)\|` per copy — **integration quality** |
| `shape_vec` | see below |
| `spread_shape` | mean distance of the copies' `shape_vec` from their centroid, ÷ 2 |
| `spread_event` | fraction of copies not sharing the modal **event class**, ÷ `(1 - 1/(E+1))` |
| `error_ratio` | `sigma_E(t) / sigma_E(0)` — see below |
| `ensemble_spread` | `max(spread_shape, spread_event)` |

**`shape_vec`** — map the configuration to a point on the unit sphere, which quotients out
translation, rotation and scale so only the *shape* of the triangle remains. Mass-weighted Jacobi
coordinates, then the Hopf map:

```
rho = r1 - r0                    lam = r2 - (m0 r0 + m1 r1)/(m0+m1)
mu_rho = m0 m1/(m0+m1)           mu_lam = m2 (m0+m1)/M
rt = sqrt(mu_rho) rho            lt = sqrt(mu_lam) lam
A = |rt|^2   B = |lt|^2   I = A + B
p = rt.x lt.x + rt.y lt.y        q = rt.y lt.x - rt.x lt.y
n = ( (A-B)/I , 2p/I , 2q/I )    then normalise
```

**`spread_event` is defined over the *event class*, not the terminal outcome.** The event class
is **which pair is currently the tightest binary**, evaluated at every sync boundary and joined
with the terminal `(state, detail)` for copies that have terminated.

The terminal outcome was explicitly *rejected* as the contributor and must not be reinstated. It
is terminal-grain and inverts under lockstep: early in the march nothing has terminated, every copy
agrees, and the field reports maximum confidence at exactly the playhead where least is known. The
event class is defined at every playhead and needs no gate. Measured, near-field 32×32, nonzero
pixels of 1024: at `t_max = 8`, **110 against 0**; at `t_max = 13`, 165 against 22, strictly nested,
with none flagged by the terminal statistic alone. On the pixels both flag the lead time is zero —
the gain is coverage and horizon-independence, not earliness.

The `(state, detail)` encoding of §2.4 stays as the **outcome**, for classification and rendering.
It is correct; it is simply not the spread contributor.

Note also that the playhead value can *un*-fire, since the tightest-pair identity fluctuates:
within one `t = 13` run, 130 of the 165 pixels that ever disagree have re-agreed by the horizon.

**Do not latch it unguarded.** Of those 130, **129 were at a near-tie** — second-tightest over
tightest separation below 1.1, median 1.0030 — so the copies disagreed about which pair is
*tightest* without having diverged. A running max lights 165 pixels where 35 have genuinely
diverged. The tie ratio cannot be the guard either: genuine disagreements also sit near 1 (median
1.0797). **Persistence is the guard** — artefacts last one boundary (median run 1, max 2), genuine
divergence persists (median run 10) — and a run of 3 admits 0 of 130 artefacts. Dump
`spread_event_max` (unguarded) and `spread_event_latched` (guarded, joined with the playhead value)
alongside, and keep the playhead value as the field `ensemble_spread` uses.

**`d_min` is primary; `r_coll` is a recorded parameter.** The collision label is *derived* from
`d_min`, not the reverse. Measured on near-field 64×64: the collision fraction runs 0.0000 → 0.0242
→ 1.0000 across `r_coll/R ∈ {1e-4, 1e-3, 1e-2}` while the whole grid's `d_min/R` spans less than one
decade (5.909e-4 to 4.931e-3). No threshold in that range is a physical event boundary — it is a
readout of the `d_min` distribution. `r_coll` must appear in every output header. The default
`1e-3` separates tail from bulk **on this slice** and claims nothing more.

**`error_ratio` — the free correctness check.** Each trajectory conserves its own energy exactly, so
the ensemble's *spread* of energies is fixed at `t=0` and must stay constant. Any growth is **pure
integration error**, with no threshold and no tuned constant:

```
error_ratio = sigma_E(t) / sigma_E(0)      -- exactly 1.0 under exact dynamics
```

Two requirements, both learned the hard way:

- **The statistic inside it must be NaN-safe *and* sensitive to a single wild copy.** Both halves
  are required and they pull against each other. A standard deviation returns NaN the moment one
  copy is non-finite — precisely the pathological pixel this field exists to flag. MAD
  (`1.4826 × median|x - median|`) fixes that but overshoots: with `E+1 = 8` copies, one bad value
  sits above the median of eight deviations and is arithmetically invisible, so an estimator that a
  single wild copy cannot move is one that cannot *see* one. **Use the maximum deviation from the
  median**, with non-finite treated as an infinite deviation — NaN-safe by construction, and the
  correct answer rather than "could not compute".

  Measured on near-field 32×32 at `t=13`, f64, over 23 damaged and 1001 healthy pixels — separation
  is the damaged median over the healthy p99:

  | estimator | damaged median | healthy p99 | separation |
  |---|---|---|---|
  | MAD | 1.1369 | 1.0756 | **1.06** |
  | max deviation | 60.864 | 1.0228 | **59.51** |

  A pixel whose worst copy drifted 120× the total energy read 1.1369 under MAD — inside the healthy
  p99. Dump the MAD-based ratio alongside as `error_ratio_mad`; do not gate on it.
- **Aggregate by `max` over pixels, not median.** Max-aggregation tracks damage at Spearman +0.956
  against +0.599 for median. Treat it as a **boolean flag**; its magnitude is unstable.

  **These two bullets are independent decisions that both land on the word "max".** The +0.956 figure
  compares `error_ratio_max` against `error_ratio_median` *across* footprints, correlated against
  per-quad exponent damage. It says nothing about which estimator belongs *within* a footprint, and
  a per-pixel correlation of the two within-footprint estimators is a different measurement
  entirely (measured: −0.035 for MAD, +0.032 for max deviation — neither in tension with +0.956).

**Do not build an `L_z` version.** Released from rest, `v = 0`, so `L_z = 0` for *every* copy and
`sigma_Lz(0) = 0` — the ratio is `0/0`. Structurally undefined for this entire configuration family.

**Note what `error_ratio` cannot see.** It measures *integration error only*. Because each trajectory
conserves its own energy, a pixel whose copies have flown to completely different outcomes still
carries exactly its starting energy spread. It says nothing about whether the ensemble has
decorrelated — that is what `ensemble_spread` is for. The two are independent and both are needed.

---

## 5. Verification — do this before trusting any output

Three bugs in the reference work failed **silently** and looked like physics. Each has a signature.

1. **Wrong equations of motion.** Symptom: energy drift that **does not fall when you reduce the step
   size**. Accuracy insensitive to step size means a wrong *equation*, never a step-size problem.
   **Test: finite-difference the Hamiltonian and compare against your analytic derivatives.** This
   caught two sign errors in the reference AZ implementation that were otherwise invisible.
2. **A step-size rule that duplicates the method.** AZ's time transformation *already* shrinks the
   physical step at close approach — that is what regularisation buys. Shrinking the fictitious step
   too drove `dt -> 1e-13` and exhausted the budget, producing a false "this region is intractable".
3. **Non-finite values never satisfying a loop exit.** `NaN >= x` is `false`, so a diverged trajectory
   burns its entire step budget: **354 s against 3 s nominal**. Test `is_finite` explicitly.

**Required acceptance tests:**

| test | expected |
|---|---|
| two-body radial collision (equal masses from rest, third body far away) | passes through `d_min < 1e-10` with `\|dE/E\| < 1e-12` |
| gauge invariance: rescale ICs by `alpha in {0.25, 1, 4}`, rescale `t` by `alpha^{3/2}` | `shape_vec` spread **identical to ~10 decimals** |
| energy control: `error_ratio` at `t=13`, near-field | median `1.0000`; max over **healthy** pixels bounded (measured max 5.21, bound 10.0) |
| Burrau constants | `M=12`, `R=2.2361`, `E=-12.8167` |
| **cross-check against the Python reference at f64** | agreement to `~1e-10` on a small grid |

**Thresholds in this table are f64 thresholds.** At f32, `|dE/E| < 1e-12` is five orders below
`eps ~ 1.19e-7` and is a statement about the type rather than about the port; asserting it there
would fail for the wrong reason. Measured on gate (b): f64 `d_min 1.2881e-11`, `|dE/E| 2.9473e-14`;
f32 `d_min 1.2218e-11`, `|dE/E| 2.8553e-6`. The `d_min` half survives the cast outright, because
regularisation still carries the trajectory through collision. Report f32 numbers; gate them at f32
tolerances.

---

## 6. Reference implementation

A working NumPy implementation exists and is validated. **Port it; do not re-derive the algebra.**

| file | contents |
|---|---|
| `tb.py` | core: leapfrog, energy, pair distances, outcome classification, grid + ensemble construction |
| `tb_lc.py` | Levi-Civita transform and single-pair regularisation |
| `tb_az.py` | **Aarseth–Zare** — the one to port. Regularised Hamiltonian, RK4 in fictitious time |
| `tb_all_az.py` | AZ plus all per-pixel fields in one pass |
| `refine_test.py` | `shape_vec` (Hopf map), dispersion measures |

`tb_az.py` uses RK4, which is **not** symplectic or time-symmetric. It was chosen to prove the physics,
not to ship. If a symplectic alternative is straightforward, prefer it — but **match the reference at
f64 first**, then change one thing at a time.

---

## 7. Rust specifics

- **One kernel source, two precisions.** Generic over `f32`/`f64` (or a feature flag). This is the
  actual production architecture: f64 on CPU, f32 on GPU. **Both must be exercised from day one** —
  a specific f32 question is the reason this exists (§8, experiment 2).
- **CPU first.** Rayon over pixels. GPU comes later; correctness now.
- **No GPU, no windowing, no async.** `ndarray` or plain `Vec`, `rayon`, a PNG writer. That is all.
- **Deterministic.** Seed the jitter per pixel from `(i, j, seed)`, never from a global RNG, so any
  pixel is reproducible in isolation.
- **Output:** a PNG (colour by `state ⊕ detail`, plus a second image coloured by `ensemble_spread`)
  and a raw dump — a simple header plus a packed struct per pixel is fine.
- **Config:** a small TOML or CLI — slice region, `W`, `H`, `t_max`, `E`, `eta`, `r_coll`, `epsilon`,
  precision, shared-reference flag.

---

## 8. What this unlocks — the reason for building it

All previous measurements used 16–64 trajectories across 4–8 regions. **Conclusions repeatedly flipped
when re-tested at larger `n`** — a field excluded on a 1.2× effect turned out to be 18.8× at more
regions. That instability is the single biggest weakness in the work so far, and a uniform kernel at
`10^6` pixels ends it.

Three experiments become possible. **Do not build machinery for them; just make sure the raw dump
carries the fields.**

1. **The refinement criterion, without a scheduler.** The criterion compares a parent quad against its
   children — and **a fine uniform grid already contains every coarser scale by aggregation.**
   Aggregate 2×2 blocks to synthesise the parent, compare against the children, and the whole
   exponent machinery is testable with no quadtree at all.
2. **The f32 question, natively.** There is an unresolved dispute about whether AZ is usable in f32.
   One measurement says raw energy drift is fine (f32 AZ *beats* softened leapfrog at some horizons);
   another says the ensemble-spread diagnostic breaks early. The hypothesis is §3's reference-switching
   across copies. **The real kernel in real f32 settles it**, with no emulation gap. Run with the
   shared-reference flag both ways.
3. **Statistical convergence.** Which conclusions survive at `10^6` samples that were measured at
   `10^2`.

---

## 9. Definition of done

- Renders a `1024 × 1024` slice of Burrau near-field at `t = 13` in reasonable time
- **All acceptance tests in §5 pass**
- Matches the Python reference to `~1e-10` at f64
- Runs in both f32 and f64 from the same source
- Writes PNG plus raw dump; `deep interior` terminates as a triple collision rather than hanging
- No quadtree, no scheduler, no interaction — **if it grew one, that is a bug**
