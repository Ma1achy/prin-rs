# Principia — Colour Composition

*Status: canonical. Single source of truth for the **colour occupant** — the internally-compositional
system that produces the `colour` and `brightness` values consumed by the render pipeline's
`combine` stage. Supersedes the implementation sections (§5–§9) of the shape-sphere colour-map
PDF and the mode-enumeration in `principia_debug_tooling_plan.md` §B–§G. The outer 4-slot pipeline
framing in `principia_gui_state_contract.md` §4 and the colour drill-down `principia_dd_colouring.md`
are amended to defer here (§8).*

*Design thesis: the PDF enumerates **products** where the system has a few **factors**. Nearly every
named colour map is one primitive family under different parameters; the LUT-spheres, the vMF map,
Voronoi, soft-Voronoi, basin-blend, the physics overlay, and custom N-pole are **the same primitive**.
Enumerating them is a maintenance liability and a second colouring path that does a strict subset of
what composition does. This document replaces the catalogue with an algebra, and re-expresses the
catalogue — and the entire debug-view set — as a **preset table over that algebra**.*

---

## 0. Scope & membrane position

Colour lives entirely on the **runtime-authorable hand-WGSL side** of the GPU↔CPU membrane. The
compute kernel (monomorphised Rust → SPIR-V/f64) is **untouched** by anything in this document;
determinism-law constraints do not apply here. Everything specified is **render-key**: editing it
never invalidates the survey cache and never re-integrates. (The one exception — a kernel-side
bring-up pattern — is quarantined to Appendix A and is the only sim-key item.)

The subsystem has three data sources, one read interface (§3 `ctx`), one primitive algebra (§1),
one pipeline shape (§4), one codegen path (§5). Debug views are presets over exactly this (§6).

---

## 1. Primitive algebra

A colour occupant is an expression tree over two primitive **families** and a small set of
**combinators**. Each node's output type is `vec3` (a colour) or `f32` (a scalar/lightness); the
slot's required output signature (`vec3` for `colour`, `f32` for `brightness`) constrains only the
root.

### 1.1 Family A — site-blend  →  `vec3`

The workhorse. Regenerates ~20 of the PDF modes.

```
SiteBlend {
  sites  : SiteSet          // §2 — an ordered set of unit directions p_i on S²
  kernel : Kernel           // how per-site weights are formed from d_i = n̂·p_i
  colours: SiteColours       // one colour per site
  space  : BlendSpace        // where the weighted mean is taken
}
→ blend( colours , weights(kernel, {n̂·p_i}) , space )
```

**Kernel** = `support × temperature`. This is the continuum the redesign asserts: hard assignment,
soft-Voronoi, and vMF blend are one primitive at different temperatures.

| kernel        | support   | weight law                                             | temperature |
|---------------|-----------|--------------------------------------------------------|-------------|
| `vmf(κ)`      | all sites | `w_i = exp(κ·(d_i − d_max))`                            | κ ∈ [0.5,12] |
| `topk(k,ks)`  | nearest k | softmax over the k largest `d_i`, sharpness `ks`       | ks ∈ [1,20] |
| `nearest`     | 1         | `w_i = [i = argmax_i d_i]`                              | — (= κ→∞ = topk(1,∞)) |

`nearest` is the κ→∞ / ks→∞ limit but **must be a discrete code path** (`exp` overflows). Fork (c)
is settled: the UI presents temperature as one **blend-sharpness dial with a hard detent at the top
edge** that snaps to the `nearest` path; `support` (all-site vMF vs top-k Voronoi) is the discrete
choice that distinguishes the vMF family from the Voronoi family. `d_max = max_i d_i` is subtracted
for numerical stability (shifts weights, not their ratios).

**BlendSpace** — an **explicit parameter**, not an accident of implementation:
- `oklab` — blend in OKLab (a,b) at the site colours' (L,C). Weighted-mean chroma **shrinks toward
  the boundary between sites**, so uncertainty reads as desaturation. This is a *feature* and is the
  perceptual default. It is a property of the blend space, **not** of any "hue table."
- `rgb` — blend in linear-ish sRGB. Used where site colours come from a perceptual LUT already
  (the LUT-sphere construction, §7) and further OKLab mixing is undesirable.

**SiteColours** — fork (b) is settled: **hue tables are dissolved.** A site's colour is *always* a
**swatch** (any OKLab colour). There is exactly one colour-assignment concept. The Full-OKLAB and
Okabe–Ito "hue tables" become **preset swatch-sets** (six swatches at fixed L, C, and the tabulated
hues). Colours may be authored directly, drawn from a palette generator (golden-angle, OI-cycle,
gradient A→B), or **sampled from a LUT at `i/N`** — the last of which is precisely how a LUT-sphere
is built (§7). Fidelity to the PDF is pinned by golden-image tests (§7), not by a special type.

### 1.2 Family B — field-ramp  →  `vec3` or `f32`

For everything that is *not* a directional blend: a scalar field mapped through a ramp.

```
FieldRamp {
  field : ScalarField        // → f32 (+ a validity lane, §3)
  ramp  : Ramp               // f32 → vec3   (colour)   OR   Compaction : f32 → f32 (lightness)
}
```

**ScalarField** sources (all read through `ctx`, §3):
- **payload** — any per-pixel kernel output: `state`, `ftle`, `energy_drift`, `Lz_drift`,
  `diffusion`, `d_min`, `word`/hash, decoded-IC quantities (E₀, K₀, V₀, L_z, virial, mass ratios,
  Jacobi ρ magnitudes/ratio/angle, min-pair-dist, …).
- **geometry-of-n̂** — `n_z`, `‖n̂‖`, azimuth/polar, `stability` = ½(1 − max_j n̂·b̂ⱼ) (BC distance),
  spherical harmonic `Yℓm(n̂)`, the 3-fold Turing wave-triple, 3-D value-noise octaves.
- **ctx lanes** — `quad.depth`, `quad.impurity`, `quad.spread`, `quad.priority`, `quad.cache_age`,
  `screen.uv`, `quad.uv`, etc. (§3). This is what dissolves the structural debug views.
- **derived operators** — `gradient_magnitude(map)` (finite-difference of another colour node —
  three evaluations at ±ε), `topk_margin(sites)` = `d_1 − d_2` (drives soft-Voronoi/basin edges).

**Ramp** (scalar → colour): `lut(name)`, `lerp(c0,c1)`, `diverging(c−,c0,c+)` (through a neutral),
`bands(field, n, line_col, base)` (a **band-mask**: line colour where the field falls in a periodic
band, base colour elsewhere — this is `grid`, `contours`, and the lattice classifiers `checker`,
`lat/lon-stripes`, `truchet` on (θ,φ)). **Compaction** (scalar → lightness, for the `brightness`
slot or before a ramp): `lin`, `log`, `symlog` (signed, through the midpoint), `cyclic` (phase),
`flag`. Every ramp/compaction carries an explicit **invalid-pixel colour/value** (§3, §6).

**Default ramps by field role.** Signed fields (energy, L_z, drifts) default to
diverging-through-neutral so the zero-crossing is a legible contour; positive fields to sequential;
angles to cyclic. **Magnitude / diagnostic fields default to *greyscale*, not a sequential LUT** —
because a greyscale magnitude field *is* a lightness, so the same field doubles as a natural
**brightness occupant** (§4.1). Polarity is **per-field**, set so the *salient* end is bright, and
is not uniform:

- **FTLE, diffusion, ensemble spread → white = high** (the *magnitude* pops — chaos, phase-space
  spreading, or ensemble dispersion; high value = bright, the standard convention, so for spread
  white = low certainty and black = high certainty);
- **time-to-event (`t_end`) → white = low / early** (quick-resolving pops; late and bounded darken,
  keeping bounded = black consistent with §1.4).

The payoff is composition: `combine(colour = shape-sphere map, brightness = FTLE-greyscale,
Replace-L)` modulates the position map's lightness by chaos (unstable brightens, regular darkens) —
a clean bivariate encoding (§4.1). Every default here is customisable; the greyscale and its polarity
are only the defaults, chosen so these fields compose well as the brightness channel.

### 1.3 Combinators  →  `vec3`

Compose sub-results. All are `vec3(+ctx) → vec3`.

- `mix_const(a, b, t)` — constant-weight blend (the "blend with second map" control).
- `mix_field(a, b, field)` — blend weight from a scalar field.
- `bandmask(base, field, band, line)` — overlay lines/tiles where a field is in-band (grid/contours
  over a base map; **quad-boundary overlay** on the normal render is this with
  `field = distance-to-quad-edge`).
- `site_overlay(base, SiteBlend*)` — additive site blobs over a base. **The physics overlay is
  exactly this**: `site_overlay(base, SiteBlend{ sites = physics(m), kernel = vmf(κ), colours =
  per-site })`. It was never a distinct node; it is Family A used as a combinator, with a physics
  site generator (§2) and per-blob strength `s`.

*The whole PDF catalogue is: two primitive families + four combinators. Adding a new map is wiring,
not a new pixel function.*

### 1.4 Categorical colour-assignment — the outcome-state default palette

Categorical fields (outcome `state`, encounter-word hash, F₂ conjugacy class) map through a
**palette**: one colour per class, hue carrying *identity, not magnitude*. Most categorical palettes are
**generated** — golden-angle for the many-class word/basin field (adjacent basins stay hue-separated),
an Okabe–Ito cycle for small class counts. This is the categorical case of §1.1's dissolved colour
assignment: a class→colour map is a swatch-set.

The outcome **`state`** field is the one categorical field with a **canonical default palette**, because
its nine terminal classes carry structure worth encoding in the colours themselves:

| class | colour | sRGB | | class | colour | sRGB |
|-------|--------|------|-|-------|--------|------|
| collision 1–2 | red     | `#DE2D2D` | | body 1 escape   | yellow  | `#F0DE32` |
| collision 1–3 | green   | `#2EBC4E` | | body 2 escape   | magenta | `#E034C6` |
| collision 2–3 | blue    | `#3462E0` | | body 3 escape   | cyan    | `#30C8DC` |
| bounded        | black   | `#141418` | | collision @ t=0 | orange  | `#F29620` |
| degenerate     | white   | `#ECECF0` | |                 |         |           |

The assignment is a **mnemonic, not arbitrary**: **collisions are additive primaries keyed by the
colliding pair** (1–2 → R, 1–3 → G, 2–3 → B); **escapes are subtractive primaries keyed by the
escaping body** (1 → Y, 2 → M, 3 → C). The two event *families* (collision vs escape) are therefore
separable at a glance while the pair/body identity stays legible; the three non-generic outcomes are
bounded (black), the t=0 collision (orange), and degenerate (white). The nine classes read from the
`state` enum **plus the `detail` union** — collision → pair id, escape → body id (so the R/G/B/Y/M/C
assignment lands on `detail`, not on a second field). There is **no separate `escaper` field**: the
escaping body *is* `detail | state=escape`, so “which body escaped” is already carried by this map's
escape colours. A standalone escaper view is therefore this map **filtered to the escape classes** — a
**categorical filter** (`show class ∈ {…}, mute the rest`), which is a general operation any categorical
mode admits (“just collisions”, “just body-2 escape”), not a distinct render mode.

Like every colour assignment in the system, **this is a default, not a fixed mapping** — the
class→colour swatch-set is user-editable. It is the canonical default the render-mode catalogue's
outcome-state row inherits (that catalogue is out of scope here; this palette is the one piece of it
that is settled).

---

## 2. Site-set kinds

Sites feed Family A and `site_overlay`. Two kinds, distinguished **in the type**, because one is
pure geometry and one depends on the decoded IC.

**Static generators** — functions of parameters only; **uniform-hoistable** (computed once, bound as
a uniform array):
- `axes6` — the six ±axis poles.  `corner8` — sign-octant corners.  `ico12` — icosahedron vertices.
- `fib(N)` — golden-angle Fibonacci lattice of N points.
- `ring(N, tilt, rot)` — N points on a rotatable great circle (custom N-pole).

**Physics generators** — functions of the **decoded IC** (the mass point), evaluated per pixel from
`ctx.payload` masses; **not bakeable**:
- `BC(m)` — binary-collision loci.  `Euler(m)` — collinear configs.  `Lagrange(m)` — equilateral
  poles. On the mass-weighted shape sphere every one of these **moves with (m₁,m₂,m₃)**, and when a
  slice axis (or a tilt) touches a `z_μ` dimension the masses are **per-pixel state**, so there is no
  per-slice constant to bake even in principle.

**Hoist optimisation.** When neither basis axis nor any active tilt touches a mass dimension, masses
are constant across the slice; the codegen detects this and hoists physics sites to uniforms —
**semantics per-pixel, cost per-slice** when possible, per-pixel otherwise.

**Preview corollary.** A shape-sphere / equirect preview is a single sphere and therefore renders at
**one mass point**. Default: the slice-centre masses `z₀`. When a pixel is inspected, the inspected
pixel's masses. The preview must **state which mass point it is showing**, because physics landmarks
(and any physics-dependent colouring) are only meaningful relative to it.

---

## 3. The `ctx` contract

`ctx` is a **first-class read interface** — the single contract that production colouring, the
codegen, the egui panel, and every debug preset program against. Making it rich is what lets debug
views be presets rather than a parallel system (§6), and what makes a working debug render a **live
test of the production read path** (same buffers, same bindings, same codegen). All lanes are
render-key.

| lane          | fields |
|---------------|--------|
| **screen**    | `pixel` (ivec2), `uv` (screenspace, vec2), `target_dims` (ivec2) |
| **chart**     | `slice_uv` (vec2 in [0,1]²), `z` (the full 8-D latent at this pixel, chart triple applied), `chart_id` |
| **quad**      | `index`, `depth`, `tl` (slice coords), `centre` (slice coords), `uv` (within-quad vec2), `state` (enum), and summary stats: `impurity`, `spread`, `suspect_frac`, `priority`, `cache_age`, `sample_count` |
| **tile/sample** | `tile_index` (within quad), `sample_index`, `N` (samples/quad), `E` (ensemble) |
| **payload**   | every per-pixel field written by the kernel — `state`, `ftle`, `energy_drift`, `Lz_drift`, `diffusion`, `d_min`, `word`/hash, `t_end`, decoded-IC quantities, masses, … |
| **validity**  | the sentinel/predicate lane paired with **every** field: `ftle_valid`, the diffusion `−1` sentinel, `sd_is_failed`, out-of-chart / saturated flags, `ftle_valid` etc. |

**Validity is not optional.** Every `ScalarField` returns `(value, valid)`. Every `Ramp`/`Compaction`
has an explicit **invalid colour/value**. Without this, debug views silently lie at exactly the
pixels they exist to expose (a NaN FTLE would ramp to *some* colour and look like data). The default
invalid colour is a conspicuous out-of-gamut-adjacent tone (spec: a fixed magenta), overridable per
node.

**Fragment-side recompute.** Because `ctx.chart.z` is present and the decode/encode are portable
WGSL, the fragment stage can *recompute* cheap quantities (decode `z` → shape/energy; `encode(decode
(z))` residual). This is what dissolves most of the §A kernel modes (§6) and enables **agreement
presets** (fragment-decode vs kernel-payload) as live cross-implementation checks.

---

## 4. Pipeline shape

Two-tier mutability, and it is deliberate: **the backbone is a fixed typed topology with switchable
occupants; the post chain is a variable-length ordered list.** Different mutability because they are
different kinds of thing.

```
            ┌─────────────┐
 colour  →  │             │
(Option)    │   combine   │ → post[0] → post[1] → … → post[k]  ═╡ display stage ╞═→ canvas
brightness→ │ (Replace-L/ │      (ordered chain, ≤ 8,           (fixed terminal,
(Option)    │  Multiply)  │       each vec3→vec3 +ctx)           settings not nodes)
            └─────────────┘
```

### 4.1 Backbone & `Option` occupants

`colour` and `brightness` are **`Option<Occupant>`**. Topology is **fixed** — occupants are switched
off, never structurally deleted. Fork settled: **ship None occupants; node "removal" sets None and
renders the node ghosted, not vanished** (keeps the topology legible — you can see where to click to
restore — while offering the gesture of removal). In the Rust state contract this is literally
`Option<Occupant>`: a free enum variant for typed `SetField`, serialisation, and undo. Structural
deletion would make topology variable and break the "occupants are data" invariant for zero semantic
gain.

**`None` is the identity element of `combine`.** Truth table (combiner-dependent):

| colour | brightness | Replace-L result | Multiply result | meaning |
|--------|------------|------------------|-----------------|---------|
| C      | B          | `OKLab(L=B, a=Cₐ, b=C_b)` | `C · B` | normal |
| C      | **None**   | `C` (colour's own L kept) | `C · 1` | "just the colour map" (default) |
| **None** | B        | `OKLab(L=B, 0, 0)` | `white · B` | **greyscale of the brightness field** (also the most CVD-robust encoding possible) |
| **None** | **None** | flat mid-grey `OKLab(0.6,0,0)` | flat mid-grey | well-defined, harmless, instantly visible |

This gives, for free, exactly the "use anything as colour, anything as brightness, or neither"
requirement: any field can occupy either slot (channel is independent of source — the *only*
constraint is the output signature, `vec3` vs `f32`), and either slot can be empty.

**The canonical bivariate family — and the honest replacement for "Stability × Hue".** The most
useful bivariate encodings are `colour = hue-carrier × brightness = magnitude-field`. The archetype
is **n̂ × ⟨field⟩**: shape-sphere position in the hue, a real dynamical magnitude in the lightness —
`FTLE × n̂`, `diffusion × n̂`, `ensemble-spread × n̂`, `time-to-event (t_end) × n̂`. Because those
magnitude fields default to greyscale-white=high (§1.2), the field drops into the brightness slot and
composes directly under Replace-L: chaotic / high-diffusion regions brighten the position map and
regular regions sink to black (and for `t_end`, with its inverted polarity, early-event regions
brighten while late / bounded regions darken). This **replaces the deleted "Stability × Hue"** — and it
replaces it with several encodings rather than one, because **"stability" is not a single metric**:
each brightness partner here is a distinct, named, *real* quantity (Lyapunov divergence, phase-space
diffusion, ensemble dispersion, time to the terminal event). `t_end` is deliberately **time to the
terminal event** — escape *or* collision *or* any class-ending event — not "escape time"; bounded
orbits reach no event within the window (sentinel).

Generally, the bivariate space is the **product ⟨hue carriers⟩ × ⟨brightness fields⟩**, not a fixed
list: categorical hue carriers compose the same way (outcome-state × FTLE), and any scalar can take
the brightness role. It is therefore a *row-family* of the render-mode catalogue, generated by the
backbone, not a set of hand-authored modes.

### 4.2 Post chain

Replaces the single post slot. An **ordered array of post nodes** (data), each `vec3(+ctx) → vec3`,
drawn from the combinator set (§1.3) and the post families: `site_overlay` (physics blobs, IC-dependent
per §2), `bandmask` (contours / quad-boundaries), `mix_const` / `mix_field` (blend-with-map),
fuzziness (ensemble-spread modulation), invert, tone. Codegen composes them in declaration order; the
node editor renders them as a linear run you can **insert into and reorder**. Bounded at **≤ 8** to
keep compile time and per-pixel cost sane. This is what unlocks "overlay *and* fuzziness *and* invert
simultaneously," which the single slot could not express.

### 4.3 Display stage (terminal, outside the pipeline)

A **fixed terminal stage** applied to the finished output, **as settings, not nodes** — the codegen
never sees them and they are not part of any occupant or preset:

```
gamut_clamp  →  cvd_sim(mode)  →  render→display scale
```

**CVD is here, not in the pipeline, and this is a category correction, not a convenience.** CVD
simulation asks "what does this *finished encoding* look like to a deuteranope/protanope/…?" — it
models the **viewer**, not the visualisation. That question is only well-posed on the output, so
anything applied after it is meaningless and anything that lets it be misplaced makes a diagnostic
that can lie. It is an **accessibility/display setting** alongside `render_scale`, and it therefore
applies uniformly to *everything* — main render, sphere preview, equirect unwrap, node thumbnails —
which is exactly what an accessibility audit wants: check the whole instrument at once. The CVD
matrices and linear-sRGB path are the Viénot/Brettel forms already in the reference artefacts.

---

## 5. Codegen & eject

**Graph → readable WGSL.** The occupant tree compiles to deterministic, readable WGSL:
- **one named function per node**, with stable identity across recompiles;
- comments carrying the node name and its parameter values;
- parameters bound as **uniforms**, so slider tweaks rebind rather than recompile;
- a small shared WGSL library the codegen calls into: `vmf_weight`, `nearest`, `topk`, `oklab↔srgb`,
  `lut_sample`, the static site arrays, the field functions, the decode/encode port. (This library
  is what replaces the 33 bespoke pixel functions.)

**View code.** Every pipeline stage exposes its generated WGSL snippet for reading. The graph *is* a
legible derivation of the shader — which doubles as pedagogy.

**Eject (fork (a), settled: per-node primary, whole-slot also available).** "Edit the shader" ejects
generated code into the existing **custom-WGSL occupant**:
- **per-node** — eject a single node's function to editable text; the rest of the slot stays
  graph-composed and live. Requires the codegen to expose **stable per-node function boundaries** and
  to splice an edited node back among generated peers, with **dependency tracking** so an upstream
  edit does not silently orphan a downstream node's inputs.
- **whole-slot** — eject the entire occupant function as one blob (the simpler mechanism; identical
  to today's custom-WGSL occupant).
- Eject is **one-way** (no decompilation). An ejected node/slot flags itself custom; its parameter
  widgets **grey out** (the code no longer derives from them); **live-compile with last-valid
  fallback** applies (already specified in `dd_colouring`); **"revert to generated"** is the way
  back. Same peek-vs-commit discipline as chart tilt vs promotion.

This is precisely the Unreal material-graph → HLSL relationship.

---

## 6. Debug views as presets

There is **one** colouring system. Debug views are presets over it (§1) reading `ctx` (§3), across
three data sources. `NORMAL` disappears — it was only ever "the user's pipeline."

**§B–§E payload field views → presets.** A field view *is* `colour = FieldRamp{field, recommended
ramp}, brightness = None, chain = []`. `f_edrift` = `energy_drift · diverging · symlog`; `f_wordhash`
= `word-hash · categorical`; `f_state` = `state · categorical palette`. As presets they inherit
compaction override, palette swap, the post chain, and per-stage shader visibility for free. The
"debug mode" dropdown becomes a **section of the preset library**.

**§F structural views → `ctx.quad` presets.** `s_depth` = `ctx.quad.depth → Viridis`; `s_impurity`
= `ctx.quad.impurity → Magma`; `s_state` = `ctx.quad.state → categorical`. Quad **boundary lines**
become a **post-chain `bandmask` step** on `distance-to-quad-edge` — strictly better than a mode,
because you can now overlay quad boundaries on the *normal* render.

**§A kernel modes → mostly presets, via fragment-side recompute (§3).** With `ctx.chart.z` present
and the decode/encode ported to WGSL:
- **UV view** = `colour = ctx.screen.uv → RG` (fragment addressing) or `ctx.quad.uv → RG` (structural
  addressing).
- **DECODE view** = fragment-side decode of `ctx.chart.z`, coloured.
- **ROUNDTRIP** = fragment-side `encode(decode(z))` residual, ramped.
- **Agreement presets** (new, and better than the originals) = `|E(fragment-decode) − ctx.payload.E₀|`
  and friends: WGSL-decode vs Rust-decode. These simultaneously test **write-addressing** (a dispatch
  scramble shows as spatial disagreement) and are a **live cross-implementation check** between the
  two decode ports — the project's two-references discipline, running on every debugged frame.

**Discipline 1 — debug presets ship locked.** Their diagnostic value is that ROUNDTRIP-red means the
same thing every time. Editing a debug preset **forks it to custom** via the same one-way eject; it
never mutates the named preset.

**Discipline 2 — validity first (see §3).** Every debug field carries its validity lane and every
ramp an explicit invalid-pixel colour, or the views lie at the pixels they exist to expose.

**Net:** one colouring system · three data sources (sample payload · quad attributes · fragment
recompute) · presets all the way down. `debug_tooling_plan` §B–§G are re-expressed as a preset table
(§7); only Appendix A remains kernel-side.

---

## 7. Preset table & golden-image obligation

Every currently-specified map (the PDF's Artefact-1 colour maps, Artefact-2 patterns, special modes,
the physics overlay) and every debug view is **recreated as a composition preset**. Representative
rows (schematic — full table lives with the preset library):

| preset | family / expression |
|--------|---------------------|
| VMF OKLAB | `SiteBlend{ axes6, vmf(κ=3), swatches=FullOKLAB-set, oklab }` |
| VMF Okabe–Ito | `SiteBlend{ axes6, vmf(κ), swatches=OI-set, oklab }` |
| Viridis (LUT-sphere) | `SiteBlend{ ring(16)+2poles, vmf(κ=4), colours=lut(viridis, i/N), rgb }` |
| Voronoi 6 | `SiteBlend{ axes6, nearest, swatches }` |
| Soft-Voronoi / basin-blend | `SiteBlend{ axes6 or fib(N), topk(2,ks), swatches }` |
| Octant / hemispheres / icosa / fibonacci | `SiteBlend{ corner8 / axes6 / ico12 / fib(N), nearest, swatches }` |
| Dot lattice | `SiteBlend{ fib(N), nearest, swatches }` masked by `bands(topk_margin…)`, `bgCol` base |
| Custom N-pole | `SiteBlend{ ring(N,tilt,rot), vmf(κ), swatches }` |
| Direction cosines | `FieldRamp` per channel `½(1+n)` |
| Checker / stripes / truchet | `bandmask(base, lattice-classifier(θ,φ), band, tileCol)` |
| Grid / contours | `bandmask(vmf-base, θφ-grid / iso-hue, band, lineCol)` |
| Gradient magnitude | `FieldRamp{ gradient_magnitude(vmf-map), lerp }` |
| Perlin / harmonics / Turing | `FieldRamp{ noise / Yℓm / wave-triple, lerp or diverging }` |
| **Stability × Hue** | pipeline preset: `colour = SiteBlend{axes6,vmf,swatches}` · `brightness = FieldRamp{stability, lin}` · `combine = Replace-L` |
| Physics overlay | post step `site_overlay(base, SiteBlend{ physics(m), vmf(κ=11 BC / 9 EL), per-site colours }, s)` |
| `s_depth`, `f_ftle`, ROUNDTRIP, … | §6 presets over `ctx` |

**Verification obligation.** The **two React reference artefacts are the oracle** (`ColourSphere` =
Artefact 1, `PatternSphere` = Artefact 2, with `physicsOverlay`/`stability_hue`). Every recreated
preset ships with a **golden-image test** against the corresponding reference output. This is a
**cross-implementation check in the project's established style** (like the shared-kernel-vs-
independent-integrator convergence reference): the composition engine and the reference artefact are
two implementations of the same maps, and agreement to tolerance certifies the port. A preset is not
"done" until its golden image matches.

---

## 8. Amendments to other docs

- **`principia_gui_state_contract.md` §4–§5** — "fixed 4-slot pipeline; occupants are data" becomes
  **"fixed typed *backbone* (`colour`/`brightness` → `combine`) with `Option` occupants · variable-
  length ordered *post chain* · fixed terminal *display stage* (gamut → CVD → scale, settings not
  nodes)."** The occupant model (source+reduction / channel / mapping, channel independent of source)
  is retained and points here for the compositional interior. CVD moves out of the pipeline.
- **`principia_dd_colouring.md`** — the vMF engine, LUT-sphere, physics overlay (Eq. 11), and
  combiner/compaction forms are retained as the *primitives* of §1 and cross-referenced; the
  mode-by-mode presentation is superseded by the preset table (§7). The `combine` L-ownership rules
  (Replace-L / Multiply) are unchanged and referenced by §4.1.
- **The shape-sphere colour-map PDF §5–§9** (implementation) — superseded by this document. The PDF
  remains the design rationale and the perceptual/CVD reference; its parameter ranges
  (L∈[0.35,0.90], C∈[0.05,0.22], κ∈[0.5,12], f∈[2,14], N∈[12,96], ks∈[1,20], s∈[0,1]) are adopted.
- **`principia_debug_tooling_plan.md` §B–§G** — re-expressed as the debug preset table (§6). §A is
  reduced to Appendix A.

---

## Appendix A — kernel-side bring-up mode

The **only** debug item that is *not* a render-key preset and *does* touch the kernel. A single
minimal mode — **"kernel writes a known pattern (e.g. `ctx`-derived UV or a fixed ramp) instead of
physics"** — for the bring-up situation where payload writes are so broken that the agreement /
fragment-recompute presets (§6) cannot even run (they depend on a trustworthy payload). It is
sim-key (it changes what the kernel computes), monomorphised on the Rust side under the determinism
law, and demoted from the debug catalogue to this appendix precisely because it is the sole exception
to "debug is presets over a shared data layer." Once payload writes are trusted, everything else in
§6 supersedes it.
