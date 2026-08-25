# Chart reference — the maths, for implementation

**Source:** `principia_spec_revised.tex` §§ decoder, views, Burrau family. Equation numbers below
are that document's. **Where this and the LaTeX disagree, the LaTeX wins** — this is a
transcription for implementation, not a new derivation.

**One-line summary of the architecture.** A chart is a map `Φ : [0,1]² → Y` into some intermediate
space, followed by a **shared decoder** `D` and a **canonicaliser** `C`. Every chart feeds the same
`D ∘ C`. **The integrator is chart-agnostic** — it receives `(mᵢ, rᵢ, pᵢ)` and never knows which
chart produced them. So adding a chart means adding a `Φ`, nothing else.

```
(u,v) ∈ [0,1]²  --Φ-->  chart space  --D-->  (m, r, p)  --C-->  canonical IC  --> integrator
```

**Indexing note.** The LaTeX uses 1-based body indices in the Burrau section and 0-based in the
decoder. This document is **0-based throughout**. Inner pair is `(0,1)`, outer body is `2`.

---

## 0. The shared decoder `D` — every chart ends here

### 0.1 Masses

```
μ_k = μ_max · tanh(z_μk)                      k = 1,2      (logit saturation)
(m0, m1, m2) = softmax(0, μ1, μ2)
M01 = m0 + m1        M = m0 + m1 + m2 = 1
```

If `M01 < ε` emit `DEGENERATE(M01_TINY)`. `μ_max = 4` is the recorded default (an open
verification item — check before relying on it).

### 0.2 Configuration — hyperspherical mass-weighted Jacobi

Reduced masses:

```
μ_ρ = m0·m1 / M01              μ_λ = m2·M01        (M = 1)
```

**Canonical-frame decode** (the recommended path — it fixes rotation and scale gauge by
construction, so the canonicaliser is a no-op away from the seam):

```
ρ̃ = R̃ cos α · (1, 0)                        R̃ = 1        (scale gauge)
λ̃ = R̃ sin α · (cos β, sin β)                β ∈ [0, π]   (mirror fixed)

α = α_min + (π/2 − 2α_min)·σ(z_α)            σ = logistic sigmoid
β = π·σ(z_β)
```

`α_min` is a buffer keeping `‖ρ‖` away from zero. **Note the orientation, which is easy to get
backwards:** `‖ρ̃‖ = cos α`, so **small α is a LARGE inner-pair separation**; `α → π/2` is a tight
inner pair with a distant third body (hierarchical).

**Unweight and reconstruct positions:**

```
ρ = ρ̃ / √μ_ρ                  λ = λ̃ / √μ_λ

r01 = −m2·λ
r2  =  M01·λ
r0  = r01 − (m1/M01)·ρ
r1  = r01 + (m0/M01)·ρ
```

COM is at the origin by construction. `I = ‖ρ̃‖² + ‖λ̃‖² = R̃² = 1`.

### 0.3 Momentum — free Jacobi momenta

```
q_k = q_max·(2σ(z_qk) − 1)        k = 0..3
p_ρ = (q0, q1)                    p_λ = (q2, q3)

p0 = −p_ρ − (m0/M01)·p_λ
p1 =  p_ρ − (m1/M01)·p_λ
p2 =  p_λ
```

> **Transcription hazard.** The `m0` and `m1` factors are **crossed** relative to the position
> reconstruction (positions take `−m1/M01` on `r0`; momenta take `−m0/M01` on `p0`). This is as
> written in the spec. Transcribe it, then verify by asserting `Σpᵢ = 0` to machine precision — the
> test that catches a swap.

### 0.4 Canonicalisation `C`

With the canonical-frame decode both steps are no-ops away from the seam. Implement anyway, for
charts that bypass it:

```
φ = atan2(ρ_y, ρ_x)     rotate all rᵢ, pᵢ by R(−φ)
if λ_y < −δ_λ:          mirror (δ_λ = 1e−12 deadband)
```

### 0.5 Scale gauge

```
ℓ = √I        rᵢ ← rᵢ/ℓ        pᵢ ← √ℓ·pᵢ
```

No-op when `R̃ = 1` was enforced. **Note the asymmetric powers** — positions divide by `ℓ`,
momenta multiply by `√ℓ`. That is what makes the transformation canonical.

### 0.6 Energy normalisation — optional, and forbidden on some charts

```
η_E = √((E* − U)/K0)          pᵢ ← η_E·pᵢ          feasibility: E* ≥ U
```

**Each chart carries a `forbids_energy_normalisation` flag.** It must be `true` for `(Lz,E)` and
`(Lz,K)`, where energy is a chart coordinate or is enforced by the momentum construction — applying
it there would collapse the energy axis. **Enforce in code, not prose:** the validation pass
refuses a chart with the flag set combined with a non-zero `E*`.

### 0.7 Degeneracy — every pixel gets a label

`DEGENERATE(reason)` for numerical failure; `COLLISION_T0(pair)` with `t_event = 0` if
`r_min(0) < r_coll`. **No pixel is ever rejected.**

---

## 1. The 8D latent chart — the reference coordinate system

```
z = ( z_α, z_β | z_q0, z_q1, z_q2, z_q3 | z_μ1, z_μ2 ) ∈ ℝ⁸
      config       momentum                 mass
```

**No latent coordinate is spent on gauge.** Rotation and scale are fixed by the canonical-frame
decode, which is why the chart is 8D and not 10D.

### 1.1 Affine slices — the "no tilt" case

Slice centre `z0` and two direction vectors `q_a`:

```
z(u, v) = z0 + (2u − 1)·s_u·q_1 + (2v − 1)·s_v·q_2
```

**Axis-aligned (no tilt):** `q_1 = ê_i`, `q_2 = ê_j` for basis vectors of ℝ⁸. There are `C(8,2) =
28` such planes. The interesting named ones:

| slice | `q_1, q_2` | what it varies |
|---|---|---|
| **shape** | `ê_α, ê_β` | triangle geometry at fixed mass and momentum |
| **inner momentum** | `ê_q0, ê_q1` | `p_ρ` — the inner pair's relative momentum |
| **outer momentum** | `ê_q2, ê_q3` | `p_λ` — the third body against the inner COM |
| **mass** | `ê_μ1, ê_μ2` | the mass simplex at fixed geometry |
| **mixed** | `ê_α, ê_q2` | one config against one momentum axis |

**Arbitrary / oblique:** `q_1, q_2` any orthonormal pair in ℝ⁸. Generate by taking two random
Gaussian vectors and Gram–Schmidt. **Report the pair used**, or the slice is not reproducible.

> **A tilt is a rotation of the 2-plane, not a re-centering.** A 2-plane in 8D has 12 tilt
> axes (6 hidden dimensions × 2 basis vectors). Raw tilts **replace** rather than compose — the
> chart constructor is the commit mechanism.

---

## 2. Invariant-momentum charts `(Lz, E)` and `(Lz, K)`

These fix geometry and mass, and use the two axes for **conserved quantities** rather than raw
controls. `forbids_energy_normalisation = true`.

### 2.1 Feasibility, and the warp that makes every pixel valid

`E = U + K` with `K ≥ 0`, so `E ≥ U`, and

```
|Lz| ≤ √(2·I·K) = √(2·I·(E − U))
```

The feasible region is the interior of a parabola with apex at `(Lz, E) = (0, U)` — the rest start.
**Domain warp** (choose `K_max > 0`, exponent `γ_K ≥ 1`):

```
K(t)     = K_max · t^γ_K
E(t)     = U + K(t)
L_max(t) = √(2·I·K(t))
Lz(s,t)  = (2s − 1)·L_max(t)
```

with `(s,t) = (u,v) ∈ [0,1]²`. **This maps the unit square onto the feasible interior, so no pixel
is infeasible by construction** — which is why the warp exists rather than clamping.

For `(Lz, K)`: identical, with `K* = K(t)` directly. Simpler, since `K ≥ 0` is the natural
constraint.

### 2.2 Deterministic momentum construction

Given target `Lz` and `K*`, construct `pᵢ` with `Σpᵢ = 0`, `K = K*`, `Lz` as specified. Work in
velocities `vᵢ = pᵢ/mᵢ`; let `J(x,y) = (−y, x)` and `⟨a,b⟩_m = Σ mᵢ aᵢ·bᵢ`.

**(i) Minimal-energy rigid rotation realising `Lz`:**

```
ω = Lz / I              vᵢ^(L) = ω·J rᵢ              K_min = Lz²/(2I)
```

**(ii) A direction field that adds energy without changing `Lz`.** Deterministic seed family, tried
in order:

```
primary:    (ρ̇, λ̇) = (ρ, 0)
fallbacks:  (0, λ) ; (Jρ, 0) ; (0, Jλ)
```

Convert each to particle velocities:

```
v0 = −(m2/M)·λ̇ − (m1/M01)·ρ̇
v1 = −(m2/M)·λ̇ + (m0/M01)·ρ̇
v2 =  (M01/M)·λ̇
```

Then project out COM drift and angular momentum:

```
w⁽⁰⁾ = v
c = Σ mᵢ w⁽⁰⁾ᵢ                    w⁽¹⁾ᵢ = w⁽⁰⁾ᵢ − c/M
β_L = L(w⁽¹⁾)/I                   w⁽²⁾ᵢ = w⁽¹⁾ᵢ − β_L·J rᵢ        where L(w) = Σ mᵢ(rᵢ × wᵢ)_z
wᵢ = w⁽²⁾ᵢ / √(‖w⁽²⁾‖²_m)
```

**Seed selection:** take the first seed with `‖w⁽²⁾‖²_m > ε_w` (default `1e−10`); if several
qualify, **choose the largest `‖w⁽²⁾‖_m` for conditioning**. Emit `DEGENERATE` only if all four
fail.

**(iii) Mix to the target kinetic energy:**

```
K* = E − U          (or K(t) directly on the (Lz,K) chart)
if K* < K_min:  terminal
a = √(2·(K* − Lz²/(2I)))
vᵢ = vᵢ^(L) + a·wᵢ                pᵢ = mᵢ·vᵢ
```

**Verification that can fail:** assert `Σpᵢ = 0`, `Lz(p) = Lz_target`, and `K(p) = K*`, all to
machine precision, over random `(u,v)`. Three independent constraints; the construction should hit
all three exactly.

---

## 3. Shape-sphere chart `(θ, φ)`

The quotient of configuration space by translation, rotation and scale. `Φ_S²(u,v) = (n(θ,φ),
m_fixed, p_fixed)` — masses and momenta held, only the **shape** varies.

### 3.1 Forward map (already implemented as `shape_vec`)

```
a = ‖ρ̃‖²    b = ‖λ̃‖²    I = a + b
p = ρ̃ₓ λ̃ₓ + ρ̃_y λ̃_y                 (dot)
q = ρ̃_y λ̃ₓ − ρ̃ₓ λ̃_y                 (NEGATIVE of the standard 2D cross — sign convention matters)
n = ( (a−b)/I , 2p/I , 2q/I ) , normalised
```

### 3.2 Inverse — closed form

Given `n` on the sphere, moment of inertia `I` and fibre phase `φ_f`:

```
a = I(1 + n₀)/2      b = I(1 − n₀)/2      p = I·n₁/2      q = I·n₂/2
                                          (p² + q² = ab holds identically)

ρ̃ = √a·(cos φ_f, sin φ_f)
λ̃ = √b·(cos ψ,   sin ψ)        with  ψ = φ_f + atan2(−q, p)
```

**The `atan2(−q, p)` sign is correct given §3.1's `q`.** With the standard cross convention it
would be `atan2(q, p)`. Verified to machine precision.

Then unweight (`ρ = ρ̃/√μ_ρ`, `λ = λ̃/√μ_λ`) and reconstruct positions per §0.2.

**Round-trip test that can fail:** `shape_vec(decode(u,v)) == n(u,v)` to ~1e−14.

### 3.3 The chart map

Two options; state which is used.

**Spherical coordinates** (simple, has poles):
```
θ = π·v        φ = 2π·u        n = (cos θ, sin θ cos φ, sin θ sin φ)
```

**Exponential map about a centre `n0`** (no poles in view; the trig is where curvature lives):
```
n(u,v) = cos(r)·n0 + sin(r)·(d/‖d‖)      d = (2u−1)·s·e1 + (2v−1)·s·e2 ,  r = ‖d‖
```
with `(n0, e1, e2)` an orthonormal frame. **This is the nonlinear chart** — use it wherever a
linearised decoder is being tested, since an affine chart makes the curvature term identically
zero.

### 3.4 Landmarks at known fixed coordinates

Useful as overlays and as tests. **Collision singularities** (two bodies coincident) are three
points on the equator; **Euler configurations** (collinear) lie on the equator between them;
**Lagrange configurations** (equilateral) are the two poles. Their exact coordinates depend on the
mass ratios — compute them from §3.1 rather than hard-coding.

---

## 4. The Burrau family chart `(ν, K)` and friends

### 4.1 The discrete family — Euclid's parametrisation

For coprime `m > n > 0` with `m − n` odd:

```
a = m² − n²          b = 2mn          c = m² + n²
```

Every primitive Pythagorean triple arises exactly once, up to leg swap.

| (m,n) | triple | angles | a/b |
|---|---|---|---|
| (2,1) | (3,4,5) | 36.9°, 53.1° | 0.75 |
| (3,2) | (5,12,13) | 22.6°, 67.4° | 0.42 |
| (4,1) | (15,8,17) | 61.9°, 28.1° | 1.88 |
| (4,3) | (7,24,25) | 16.3°, 73.7° | 0.29 |

### 4.2 Positions and masses — Burrau convention

Right angle at the origin, normalised by hypotenuse, **each mass equal to its opposite side**:

```
r0 = (0,     0    )      m0 = c/(a+b+c)
r1 = (a/c,   0    )      m1 = b/(a+b+c)
r2 = (0,     b/c  )      m2 = a/(a+b+c)
```

Jacobi vectors for inner pair `(0,1)`:

```
ρ = ( a/c , 0 )
λ = ( −ab/(c(b+c)) , b/c )
```

**Sanity check that can fail:** at `(m,n) = (2,1)` this must give the classical Burrau
configuration, `(m0,m1,m2) = (5,4,3)/12`. Note this is the **normalised** form; the historical
convention uses masses `(3,4,5)` with a differently-scaled triangle. **State which is in use** —
they are the same system up to the scale gauge, but they are not the same numbers.

### 4.3 Continuous interpolation — the `ν` axis

Euclid's formulae hold for all real `m > n > 0`, so with `ν := n/m ∈ (0,1)`:

```
a ∝ 1 − ν²          b ∝ 2ν          c ∝ 1 + ν²
acute angle:  θ(ν) = arctan( (1 − ν²) / (2ν) )
```

Everything — positions, Jacobi vectors, masses — varies **smoothly** with `ν`. Primitive triples
sit at a countable set of `ν` values; the chart sweeps between them.

**The `(ν, K)` chart:** `ν` on one axis (triangle shape, Burrau masses following it), kinetic
energy `K` on the other via §2's construction. The spec calls this the *bifurcation strip* — the
right-triangle Burrau configurations are a **1D curve** inside a 2D map, which answers directly
whether the right-angle constraint is dynamically special or merely convenient.

### 4.4 Relaxations — each gives another chart

- **Mass simplex.** Fix the triangle, sweep `Δ₂ = {(m0,m1,m2) : Σmᵢ = 1, mᵢ > 0}`. The Burrau point
  `(c,b,a)/(a+b+c)` is one distinguished location. 2D by construction — barycentric coordinates map
  straight onto `[0,1]²` with a shear.
- **Rest start relaxed.** Replace zero momenta with §0.3 or §2.2.
- **Right angle relaxed.** Sweep the apex angle away from `π/2` at fixed side ratio.

---

## 5. Implementation plan

### 5.1 One trait, one dispatch

```rust
pub trait Chart {
    fn map(&self, u: f64, v: f64) -> ChartOut;      // Φ : [0,1]² → chart space
    fn forbids_energy_normalisation(&self) -> bool { false }
    fn name(&self) -> &str;                          // goes in every dump header
}
```

`ChartOut` is whatever `D` consumes — masses, `(α,β)` or an explicit `(ρ,λ)`, and momenta or a
`(Lz,K)` request. **`D` and `C` are shared and written once.** The integrator never sees a chart.

**Charts to implement**, in this order — each is a `Φ`, and the shared decoder is unchanged:

| # | chart | why this order |
|---|---|---|
| 1 | `Latent { z0, q1, q2 }` | the general case; axis-aligned and oblique are both instances |
| 2 | `ShapeSphere { n0, e1, e2, I, phase }` | the nonlinear one — needed for the linearised-decoder test |
| 3 | `Burrau { nu_range, k_range }` | the project's namesake family |
| 4 | `InvariantLE / InvariantLK` | most machinery, most verification |
| 5 | `MassSimplex` | cheap once the decoder exists |

`BodyPlane` (today's slice) stays as a chart and must reproduce **bit-for-bit** — it is the
Python cross-check's anchor.

### 5.2 Tests that can fail

- **`Σpᵢ = 0`** to machine precision, every chart, random `(u,v)` — catches the crossed-mass
  transcription hazard in §0.3
- **`I = 1`** after the canonical-frame decode
- **Shape round-trip** `shape_vec(decode(u,v)) == n(u,v)` to ~1e−14
- **`(Lz,E)`:** the constructed momenta hit `Lz` **and** `K` to machine precision, over random
  `(u,v)`, including near the parabola boundary where `K* → K_min`
- **Feasibility:** no `(u,v)` in `[0,1]²` produces `K* < K_min` under the warp — the warp's whole
  purpose
- **Seed fallback:** force the primary seed to degenerate (e.g. `ρ = 0`) and assert a fallback is
  selected rather than `DEGENERATE`
- **Burrau at `ν = 1/2`** reproduces `(3,4,5)`, and the classical configuration to a stated
  tolerance
- **`forbids_energy_normalisation`** is enforced — a config combining `(Lz,E)` with `E* ≠ 0` is
  **refused**, and a test asserts the refusal
- **Axis-aligned latent slice with `q1 = ê_α, q2 = ê_β`** equals a direct `(α,β)` sweep
- **`BodyPlane` bitwise unchanged**, and the Python cross-check green

### 5.3 What to report per chart

Leaf count, `alpha` distribution, tree shape. **But leaf counts are slice-conditional to 4.3× —
compare within a chart, never across.** The `alpha` distribution is the safer cross-chart quantity.

**And the standing caution applies with extra force here:** a chart that produces a prettier
picture is not a better chart. The measurement is whether the criterion behaves consistently
across charts, not which chart looks best.
