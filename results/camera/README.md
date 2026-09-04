# `results/camera` — does `camera_bias` change anything, and what does the margin buy?

§4.3 wants priority to be a **product** of camera relevance and structure. `SchedCfg::camera_bias`
implements it — `priority = signal × Camera::relevance(cx, cy, half, margin)` — and defaults to
`None`. §3 offers three destination models and calls the third, *"simply widen the viewport margin
and drop prediction entirely"*, the honest baseline the others must beat. The margin is its only
parameter.

Measured 2026-09-04, `Policy::Tolerance` at committed defaults, 512² camera **zoomed to a quarter
width toward a corner** so `relevance` varies — `Camera::framing` sets `half_world` to the root
half-width, which makes every quad fully visible and the bias an identity. `rel span` is printed as
the liveness arm. Reproduce: `cargo run --release --example camera_probe`.

## The knob has exactly one regime, and it is the truncation regime

| chart | budget | margin | quads | moved | of | exhausted | rel span |
|---|---|---|---|---|---|---|---|
| `near-field` | 20000 | 0.0 / 0.5 / 2.0 | 389 | **0 / 0 / 0** | 389 | false | 0.60 / 1.00 / 1.00 |
| `near-field` | 155 | 0.0 / 0.5 / 2.0 | 153 | **17 / 16 / 12** | 93 / 93 / 81 | true | 0.60 / 1.00 / 1.00 |
| `deep interior` | 20000 | 0.0 / 0.5 / 2.0 | 1017 | **0 / 0 / 0** | 1017 | false | 0.49 / 1.00 / 1.00 |
| `deep interior` | 406 | 0.0 / 0.5 / 2.0 | 405 | **29 / 25 / 25** | 185 / 189 / 189 | true | 0.49 / 1.00 / 1.00 |

**Zero decisions move at a non-binding budget, on both charts, at every margin, with the arm
demonstrably live.** `camera_bias` reaches only `priority`, which reaches only `order_queue`, which
changes a tree **only where the ranking is truncated** — and under `k_frac < 1` a deferred quad is
re-decided next round, so when nothing truncates, everything the criterion wants eventually happens
and the ordering washes out. When the budget binds it moves 18% and 16% of shared decisions.

## The margin is a weak knob

17 → 16 → 12 and 29 → 25 → 25 across margins 0.0, 0.5, 2.0. Read `moved` together with the shared
count: at margin 2.0 on `near-field` the shared set falls 93 → 81, so the tree diverged *more* from
the unbiased baseline while the quads common to both agree slightly better. `moved` alone would
report that backwards.

So §3's option-3 baseline works, and its parameter is not where the leverage is. The swept path
(option 2) has to beat *this*, and §18's cursor bias covers the zoom case for free, since
scroll-to-zoom zooms toward the cursor.

## Two fixture failures, both mine, both about a constant that did not bind

The binding arm was first a **chosen constant, 400**, picked against a remembered 125-quad
`near-field` tree. The zoomed camera makes that tree 389 quads, so 400 did not bind and the
"binding" arm reproduced the non-binding one exactly — the same regime measured twice, which is the
failure the probe exists to avoid. The budget is now **derived per chart**: measure the
unconstrained tree, take 40%.

That fix then failed on `preset_shape_h1`, whose unconstrained arm read **19997 quads against the
20000 cap** — so `free` was a floor rather than a measurement and the derived budget inherited it.
That chart is dropped rather than reported. The general point: **a fixture constant needs deriving,
not choosing**, and a derived one needs its own guard that the thing it derived from was not itself
capped.

Under a real frame quota this whole class disappears, because the quota binds by construction.

## What this says about Phase A as a whole

All three Phase A knobs act through `order_queue`, and all three are inert or misleading in the
batch descent:

- **`camera_bias`** — 0 decisions moved unconstrained, 17/93 and 29/185 when the budget binds.
- **The frontier's bucketing** (`results/frontier/`) — no saving at `k_frac = 0.25`, 83–94% at
  frame-budget `k`.
- **The balance-forced share** (`results/balance/`) — identical trees, attribution shifting 21–93%
  with the throttle.

One conclusion reached three times: **these are frame-budget mechanisms, and the frame loop is the
precondition for measuring them rather than the phase after.**
