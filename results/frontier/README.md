# `results/frontier` — does the priority bucketing earn its place?

`src/frontier.rs` was built for §4.6's argument: *camera relevance changes for every quad on every
frame the camera moves, so during a gesture the naive version is genuinely recomputing all of it.*
It is complete, tested, and **referenced by nothing outside its own tests**.

It also did not yet do the job it is named for. `Frontier::top_k` flattens every bucket and sorts
the lot, so the per-frame `O(n log n)` the structure exists to remove was still being paid; the
bucketing saved only `reprioritise`. `top_k_bounded` is the early-out — walk bands from the top,
stop once no lower band can contend — and it is **exact**, because the derived factor is in `[0,1]`
and can only demote.

Whether it *bites* is empirical, and the failure to look for is saturation: this project has met it
at the hot mask (`n_hot = N²` in 99–100% of quads) and at the linear ramp. Measured from the
committed `.qcache` files at zero trajectories; the stored term at production settings is
`signal(Within, Median)` = `spread_median`, which they carry per quad.

Reproduce: `cargo run --release --example frontier_bands`. Full output in `output/`.

## Band occupancy is regional, and two of four charts are saturated

| chart | bands of 24 | modal share |
|---|---|---|
| `far` | 4 | **93.8%** |
| `near-field` | 6 | **83.8%** |
| `deep_interior` | 8 | 42.6% |
| `preset_shape_h1` | 10 | 23.7% |

## The verdict is `k`-dependent, and it splits exactly where the slippy map does

`scan/n` on the **visible** frontier (410 of 5461 quads), against a full sort at 1.000:

| chart | `k/n` = 0.01 | 0.05 | 0.25 | 0.50 |
|---|---|---|---|---|
| `near-field` | **0.098** | 0.351 | 1.000 | 1.000 |
| `deep_interior` | **0.061** | 0.146 | 0.810 | 0.810 |
| `far` | **0.098** | 0.351 | 1.000 | 1.000 |
| `preset_shape_h1` | **0.168** | 0.168 | 0.329 | 0.661 |

**At small `k` it saves 83–94% on every chart, saturated ones included. At the batch `k_frac = 0.25`
it saves nothing on `near-field` and `far`.**

So the structure earns its place in the **frame loop** and not in `descend_with`. A frame budget is
the small-`k` regime by construction — a handful of quads per frame against a frontier of thousands
— while the batch descent runs at `k_frac = 0.25`, where on half these charts the walk scans
everything and `top_k`'s single sort is the honest implementation. **Wiring it into the batch
`order_queue` would be adding a mechanism where it measurably does not pay.**

## Two dead arms, one caught by its own guard and one by reading the population

**The camera arm was an identity.** `Camera::framing` sets `half_world` to the *root* half-width,
so every quad lies wholly inside the viewport and `relevance` returns 1.0 on all of them. The first
cut measured the uniform arm twice and reported it as a camera result — identical to three decimals
in every cell, which is the tell. The harness now zooms to a quarter width offset toward a corner
and **asserts `relevance` varies** before reading the column.

**Then the live arm measured the wrong population.** With the whole tree in the frontier, 5051 of
5461 quads are fully off-screen and carry `relevance = 0`; their priority is 0, the k-th score never
leaves band 0, the walk can never stop, and `scan/n` reads **1.000 in every cell on every chart**.
True, and about a frontier no real frame holds — the brief's priority weights have visibility
dominating and *never compute off-screen*. The table above is the visible frontier.

Both are the same failure at two removes: *a difference can be small because both sides are right or
because both are dead*, and then *the population has to be the one the mechanism runs on*.

## A latent defect in the staleness check, found on the way

`agrees_with_rebuild` compares `top_k` against `rebuild`, and the two **disagreed on tied data**.
`top_k` flattens buckets top-down; `rebuild` reads an id-sorted list. Under a merely *stable* sort,
two entries of equal priority sitting in different bands come out in opposite orders — measured,
`0.2 × 6/7` and `0.4 × 3/7` are bitwise equal and land in bands 18 and 19. So the check with teeth
would have reported a **false** disagreement on any tied field: firing on the wrong thing rather
than failing to fire.

The ordering is now **total** — descending priority, `NaN` last, ties by id ascending — shared by
all three paths, so a disagreement between them can only be real.
`the_orderings_agree_on_exact_ties_across_bands` pins it, and asserts the fixture is an exact tie
straddling two bands so it cannot go vacuous.

And the first `top_k_bounded` inverted `band_of` analytically to get each band's upper bound. The
log round-trip does not land on the boundary: `band_upper(0)` returned `4.0616e-12`, which `band_of`
places back in **band 0** — a bound too small, stopping the walk early with a contender unseen,
unsound in the silent direction. The stopping test is now expressed through `band_of` itself
(`band_of(kth) > b`) with no inverse at all.
