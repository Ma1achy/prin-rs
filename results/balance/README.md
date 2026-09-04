# `results/balance` — what the 2:1 constraint costs, and why its own metric will not carry it

The 2:1 balance constraint says no two adjacent leaves may differ by more than one level.
`scheduler::balance_pass` implements it, `SchedCfg::balance` defaults to **false**, and
`descend_live_with` did not call it at all until this measurement — so every tree in `results/`
predating this is unbalanced. Whether to turn it on is what this measures.

Measured 2026-09-04, `Policy::Tolerance` at its committed defaults (`alpha_lo = 0.005`, dimension
floor and agreement arm on, stationarity off), 512² camera, `t = 13`, `refine_flagged = false`,
**budget 20000 — deliberately non-binding**. Reproduce:
`cargo run --release --example balance_census -- <root> "<charts>" 512 <k_frac>`.

## The prior was wrong, and it was a real question

`tests/slippy.rs` records that near-field is gap 1 at **all twenty-four** cells of
`alpha_hi × tau × n` it sweeps, and that the test has to drop to `n = 4` on `deep interior` to
produce a 2:1 violation at all. That made the pass look inert at the production `N = 8`. **That
sweep is `Policy::Alpha`-era.** Under the tolerance policy **all six charts violate 2:1**, at gap
2 or 3.

## The cost

| chart | dvar | quads off → on | **quad ×** | steps × | st/qd | fr/split k=0.25 | k=1.0 | drop |
|---|---|---|---|---|---|---|---|---|
| `config_stability` | 0.047 | 4729 → 4869 | **1.03** | 1.04 | 1.01 | 0.1134 | 0.0074 | 93% |
| `preset_shape` | 0.257 | 1837 → 2121 | 1.15 | 1.08 | 0.94 | 0.2396 | 0.0792 | 67% |
| `preset_shape_h1` | 0.309 | 2169 → 3141 | 1.45 | 1.36 | 0.94 | 0.2268 | 0.1096 | 52% |
| `preset_prho` | 0.586 | 977 → 1105 | 1.13 | 1.13 | 1.00 | 0.2754 | 0.0725 | 74% |
| `deep interior` | 1.516 | 225 → 397 | **1.76** | 1.82 | 1.03 | 0.2929 | 0.1717 | 41% |
| `near-field` | 2.098 | 125 → 189 | 1.51 | 1.51 | 1.00 | 0.2979 | 0.2340 | 21% |

**The tax is 1.03–1.76× the tree**, and it is real work: every forced split integrates
`N² (E+1) = 512` trajectories.

## `fr/split` — the number §4.4 asks for — is not a robust quantity

§4.4 asks to *"report what fraction of splits are balance-forced rather than criterion-driven."*
Measured at the two throttle settings, **every tree is identical and only the attribution moves**,
by **21–93%**. On `config_stability` the fraction is 0.113 or 0.007 — a factor of fifteen — for the
same tree. The mechanism is plain: at `k_frac = 1.0` the criterion splits its whole want-list each
round and reaches those quads first; throttled to 0.25 the balance pass gets there first and takes
the credit.

And the ordering does not survive either: **Spearman(`fr/split` at 0.25, `fr/split` at 1.0) = +0.771**
— the same statistic on the same six trees does not agree with itself across throttles. So it cannot
be used to compare charts, not merely to state an absolute share. **Quote `quad ×`**, which is
identical under both arms on all six.

*Caveat on the control, stated because it bounds the claim:* the budget is non-binding, so `k_frac`
only reorders within a round — deferred quads are re-decided next round and everything the criterion
wants eventually happens. That is **why** the trees are identical, and it means throttle-invariance
of `quad ×` is established only at a non-binding budget. Under a **frame** budget `k_frac` binds and
the trees will differ; that case is unmeasured and is the one Phase B needs.

## The geometry tax is cost-neutral in trajectories; the criterion's is not

`steps ÷ quad` runs **0.94–1.03** across all six. Against the committed cost ledger
(`results/payload/README.md`), where the *criterion's* `cpu/mem` runs **0.94–2.10**. The reason is
structural: the criterion selects on physics, which correlates with trajectory cost; balance selects
on geometry, which does not. So the geometry tax may be quoted in quads with at most a few percent
of error, where the criterion's may not.

## The mechanism is suggestive and not established

*"Balance costs in proportion to depth contrast, and depth contrast is what the criterion is for"*
reads **+0.943** against `fr/split` and **+0.771** against `quad ×`, at n = 6 where the two-tailed
5% bar is 0.886. The better correlation is on the throttle-contaminated column and the robust column
falls short, so this is recorded as a lead rather than a result.

Two candidate mechanisms were proposed and failed before it. **Tree size** — refuted by
`preset_shape_h1` (2169 quads, 1.45×, criterion splitting *more*) against `config_stability` (4729,
1.03×, splitting fewer). **Depth variance read off tree growth** — apparently refuted on the first
two charts, then found to hold on the split share; growth is confounded because it is forced splits
*plus the criterion's response*, and that response flips sign chart to chart (**+2, +14, −44, −56,
−103, +65**). A forced split changes what its descendants and its parent subsequently decide, in
both directions, and no account of the sign is offered here.

`violating_adjacencies` — adjacent leaf pairs differing by ≥2 levels, which is what the pass
actually fires on and a perimeter rather than a spread — is computed by the harness and was added
after this run. It is the next candidate.

## What this does not settle

Whether 2:1 should be on. It is a **rendering** requirement, and in this renderer it buys nothing
visible: texels are nearest-neighbour clipped to the quad box and quadtree leaves tile the root
exactly at any depth difference (`adaptive::coverage`: zero gaps, zero overlaps). What an unbalanced
tree produces is a resolution step, not a hole — and §4.5 already accepts resolution steps as honest
(*"big texels during motion are a deliberate choice"*). It becomes load-bearing under interpolation
across leaves or when quads are drawn as GPU geometry, where the standard remedy is render-side
stitching, which costs geometry and no integration.

So the cost measured here is paid in **physics** for a **rendering** property. `balance` stays a
flag; the default is a judgement about which renderer ships, not about these numbers.
