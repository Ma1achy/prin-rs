# The three prerequisites, before any refinement work

`corpus_scope.txt` settled that the criterion's input field was **re-ordered** between the corpus
kernel and production (`rho` 0.59-0.84), so no criterion comparison survives. These three answer
what has to be true before rebuilding.

Reproduce: `cargo run --release --example far_control -- 128 2000`,
`cargo run --release --example floor_justify -- 2000`, and the escape pair with
`cargo run --release --example criterion_metric -- 5 8 1e-4 13 <scratch> {production|legacy}`.

---

## 0. Found first, and it gates the other two: the harnesses were not running production

**Zero of nine** harnesses feeding the refinement work printed any provenance —
`criterion_metric`, `chart_gallery`, `balanced_march`, `oracle_audit`, `equal_budget`,
`signal_audit`, `between_vs_within`, `structure_metric`, `bivariate_colour`. That is the
`refine_flagged` failure exactly: *the failure was never the choice, it is that nothing recorded
the choice.* All nine now print `ens.provenance()`, which **derives** the declaration by diffing
against `production()`, so a config declares itself however it was built.

Three of them — the whole `error(B)` machinery — pinned `escape_rule: Reference` with
`stop_on_escape: true`, justified in a comment as *"every result in this diagnostic predates both
the distance gate and the closure criterion, and is quoted against that form"*. That expires the
moment those results are superseded, and one half contradicts a standing decision outright:
`stop_on_escape` is **off** in production, because closure certifies what escaped and is silent
on whether the displayed shape has settled.

**Measured, at `levels = 5`, both arms one binary:**

```
                    error(root)      terminated              escaped
  near-field  prod      0.29576   0.0234,   50 quads   0.0000,   0 quads
              legacy    0.29576   0.0234,   50 quads   0.0000,   0 quads
  far         prod      0.60165   1.0000, 1365 quads   0.0000,   0 quads
              legacy    0.60165   1.0000, 1365 quads   0.0000,   0 quads
  deep int.   prod      0.31026   0.9854, 1365 quads   0.0000,   0 quads
              legacy    0.25089   1.0000, 1365 quads   0.0625, 202 quads
```

`near-field` and `far` are **identical to five decimals** — the escape arm is silent there at
`t = 13`, as the standing result says. `deep interior` moves: **202 quads carry escapes under the
legacy pair and zero under production**, `error(root)` moves 24%, and every signal moves with it
(`within/median` spread `1.447e-2 -> 1.429e-1`, tenfold).

The mechanism is on record: `Reference` is `spec > 0 && receding`, which is transiently true
during a close encounter, and **0 of 895 `deep interior` trajectories firing it are still unbound
one boundary later**. Closure rejects them. So every `deep interior` `error(B)` number in the
corpus — the region the record names as the one where a criterion is decisively ahead — was taken
on a field carrying escapes the shipped criterion rejects, with `stop_on_escape` freezing the
displayed shape at them.

**Not attributed, and stated rather than implied:** the two arms move `escape_rule` *and*
`stop_on_escape` together, so this measures the **pair**. Separating them needs a third arm and
has not been run. *A knob held fixed is a knob whose effect is unattributed*, and two knobs moved
together are two unattributed effects.

Production is now the default; `legacy` is a **named control**, so the pre-rebuild numbers stay
reachable and stop being what runs when nobody chose.

---

## 2. `far` is still flat in the bulk, the tail is 0.5% of the frame, and the criteria no longer agree

At 128^2 the tail resolves better than the 64^2 pass could show. **The earlier "5-10% of the
region" was wrong** — it is above **p99**, not p90:

```
  arm        p01        p50        p90        p99       p999        max    p99/p01
  pre   9.267e-9   9.387e-9   9.466e-9   9.509e-9   9.522e-9   9.528e-9    1.026e0
  now   7.442e-9   7.499e-9   7.546e-9   7.583e-9   4.131e-6   4.407e-6    1.019e0
```

The bulk is flat to **1.9-2.6%** under both kernels, so the standing description of `far` as
featureless **survives for the bulk**. The tail is **82 pixels of 16384 (0.50%)**, in **11
components of 8-9 pixels** — neither dust nor one region, and the component sizes are what say so.

**And nothing explains it.** Cross-tabulated, the tail and the bulk are indistinguishable in every
candidate:

```
  population       n  spread p50   drift p50   t_end p50   d_min p50   ev doms      term
        tail      82    3.305e-6    7.933e-8     1.922e0    6.120e-3     0.000     1.000
        bulk   16302    7.498e-9    7.958e-8     1.922e0    6.118e-3     0.000     1.000
```

Identical drift, identical `t_end`, identical `d_min`, identical termination, and `ev doms 0.000`
so it is the **shape** arm and not the five-valued event staircase. Only the spread differs, by
**440x**. A genuine `spread_shape` feature with none of the obvious causes, and it is recorded as
unexplained rather than attributed.

**The operative question, and the answer is no.** Running the descent at a `tau` below `far`'s
bulk, the eleven criteria produce:

```
  tau=1e-4  leaves 16   2 distinct leaf sets   <- the degeneracy control, all `keep`
  tau=1e-6  leaves 19   6 distinct leaf sets
  tau=1e-8  leaves 64   7 distinct leaf sets
  tau=1e-9  leaves 64   6 distinct leaf sets
```

**Read the `tau = 1e-4` row as the control it is** — 16 leaves, everything `keep`, nothing ever
ranked. Agreement there is agreement for want of running. At every rung below the bulk the
criteria give **six or seven distinct trees**.

**A scope limit that matters, and conflating the two would be the error this measures against.**
The standing result — *"thirteen non-greedy rows agree to five digits; `within/median` on 21845
distinct values and `frac_hot_between` on 1 produce the identical leaf set"* — was measured by
`metric::replay` over a **cache**, where criteria enter purely as orderings and no `tau` exists.
This is the **descent**, with a `tau` gate, an `alpha` gate and a camera veto. They are different
machines answering different questions, and this result does not by itself overturn the replay
one. The replay question is re-answered by the `levels = 6` `criterion_metric` run.

---

## 3. `Floor` was detecting the integration bug, and `alpha` moved to its ideal value

```
region          arm  n_alpha    p10    p25    p50    p75    p90  frac<0.2
near-field      pre       60 -0.154  0.042  0.259  0.401  0.513     0.383
                now       84  0.884  0.943  1.011  1.058  1.118     0.000
body2 core      pre       56 -0.117  0.093  0.451  0.627  1.078     0.321
                now       84  0.672  0.818  1.053  1.182  1.392     0.012
deep interior   pre       28 -0.356 -0.227  0.256  0.733  1.144     0.464
                now       20 -0.030  0.777  1.265  1.463  2.288     0.150
mid-field       pre/now   20  ~0.99  ~0.99  ~1.00  ~1.01  ~1.01     0.000
far             pre/now   20  ~1.00  ~1.00  ~1.00  ~1.00  ~1.00     0.000
```

**`alpha` moves to ~1.0** — the value the fixed-Halton control has *exactly*, i.e. spread scaling
linearly with cell width and refinement paying precisely as the ideal predicts. Under the old
kernel the three chaotic regions read 0.26-0.45 with **32-46% below `alpha_lo`**, because a spread
inflated by integration failure is not reduced by halving the cell: the failure belongs to the
trajectory. `mid-field` and `far` sit at 1.00 under **both** arms — they were never broken, and
that is the control that says the shift is not an artefact of the harness.

**And it is not a mis-placed threshold.** Sweeping `alpha_lo` over `{0.05, 0.1, 0.2, 0.35, 0.5,
0.8}`, `Floor` leaf counts:

```
  near-field     pre  13 13 17 17 17 17   |  now  0 0 0 0 0 0    (46 -> 64 leaves)
  body2 core     pre   9 14 16 16 16 16   |  now  1 1 1 1 1 1    (43 -> 64 leaves)
  deep interior  pre   9  9  9  9  9  9   |  now  1 1 1 1 1 1    (22 -> 16 leaves)
  mid-field      pre   0  0  0  0  0  0   |  now  0 0 0 0 0 0
  far            pre   0  0  0  0  0  0   |  now  0 0 0 0 0 0
```

Pushing the threshold to **0.8** — four times the shipped value — still fires nothing, because the
distribution now sits above it. So this is the distribution moving, not the gate being in the wrong
place.

**One confound, named.** The two arms produce different trees, so the `alpha` distributions are
over different quad populations — the same shared-quads issue as churn, and here it is partly
circular because the tree shape is *caused by* `alpha`. The `alpha_lo` sweep does not depend on it,
which is why it is the discriminating arm and the distributions are the illustration.

### And the chart control refutes the strong reading

Concluding "`Floor` has no subject" from five Burrau regions would have been a claim about the
regions. `preset_shape` is on record as *the only tree in the corpus exercising the alpha gate* —
8 `floor` + 8 `keep`, 16 leaves against a complete 4096 — so it is the control, and it fires:

```
  preset_shape      pre   8  8  8  8  8  8   (16 leaves)   frac<0.2  0.400
                    now   8  8  8  8  8  8   (16 leaves)   frac<0.2  0.400
  preset_shape_h1   pre  15 15 16 16 16 16   (43 leaves)   frac<0.2  0.357
                    now   9  9  9  9  9  9   (46 leaves)   frac<0.2  0.167
```

**`preset_shape` is bitwise unmoved — 8 of 16 at every threshold under both kernels**, and the
record's count reproduces exactly under production. So the answer is neither "Floor was the bug"
nor "Floor is fine":

> **`Floor`'s Burrau firings were the integration bug; its chart firings are the chart.** It keeps
> its subject and stays in the descent, and every `Floor` count ever quoted from `near-field`,
> `body2 core` or `deep interior` is void.

### The mechanism differs between the two, and the chart one is a lead

On Burrau, `Floor` fired because `alpha` sat near **0** — refining failed to reduce a spread that
integration failure had inflated, which halving the cell cannot fix. On `preset_shape` `alpha`'s
p10 is **-11.868 (pre) and -12.567 (now)**: refining *multiplies* the spread by ~4000. `Floor`
fires there because the spread **grew**, not because it failed to shrink.

That is worth stating as a lead and not a conclusion. On a fractal boundary a spread that grows
under refinement is evidence there **is** structure — the children resolve a filament the parent's
coarse grid averaged over — and `Floor` responds by stopping. The standing account of
*`preset_shape` is where the criterion fails outright, the fractal core sits unrefined inside a
level-2 quad* attributes that to `Agg::Median` reading the quad resolved. **There is now a second
candidate at the same site**, and the two have not been separated: `alpha_lo` is a one-sided gate
applied to a two-sided quantity, and `alpha << 0` is not the same event as `alpha ~ 0`.
