# The timeline bisect — where the picture changed, and why this slice is fragile

Two questions, run in parallel. The first is **which commit**; the second is **why this slice
and not the presets**. They have different answers and neither subsumes the other.

---

## THE SLICE

The user's saved config, as literals, in `harness/template.rs` — `Chart::config_stability()` is
deliberately **not** called, because it did not exist before `e53223d` and a window that is
nearly right reads as a physics disagreement.

```text
z0   = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116]
dimH = 0, dimV = 1, mag = 1.0
zoom = 0.63763, panX = 0.17512, panY = 0.17727
horizon 50, r_coll 0.005, 1024^2, one sample per pixel
```

Verified at every commit before anything is integrated, and printed in the log:

```text
masses [0.32735, 0.42763, 0.24502]  max|dm| 3.52e-6
window z_beta [-0.7478, +0.5275]   z_alpha [-0.5355, +0.7398]
```

Both match the reference UI exactly. **Equal masses would mean the decode is overridden and
every number here void**; the harness asserts it rather than printing it.

### Every parameter is explicit, and an absent one is the signal

`harness/run.sh` regenerates `bisect_slice.rs` **per commit** from one template, emitting a
literal for each `EnsembleCfg` field that exists at that commit and logging the ones that do
not. Defaults are exactly what changed across this range, so a harness that read them would
have measured the defaults. The per-commit generated sources are in `harness/harness_*.rs` —
that file *is* the record of what was pinned.

```text
f4084de   22 fields   absent: escape_every escape_confirm r_esc_frac escape_all_bodies
                              escape_rule closure_k stop_on_escape dtau_mode clamp_final_step
077b092   24 fields   absent: r_esc_frac escape_all_bodies escape_rule closure_k
                              stop_on_escape dtau_mode clamp_final_step
e53223d   26 fields   absent: escape_rule closure_k stop_on_escape dtau_mode clamp_final_step
71de13f   27 fields   absent: r_esc_frac escape_all_bodies dtau_mode clamp_final_step
5cc8dec   28 fields   absent: r_esc_frac escape_all_bodies clamp_final_step
f7d2a31   29 fields   absent: r_esc_frac escape_all_bodies
```

The two colour windows are **fixed constants shared by the whole strip** (`spread 6.85e-5 ..
4.955e-1`, `drift 1e-8 .. 4e7`). An auto-ranged ramp per panel would stretch each commit's own
p1–p99 to full scale — on a question about *bleaching* that manufactures or hides the thing
being looked for. Each panel's own auto-range is printed beside it and moves by under 25%.

---

## FOUR NULL CONTROLS, ALL BITWISE

A bisect whose harness is not deterministic across a checkout measures the checkout.

| pair | why it must be identical | result |
|---|---|---|
| `f4084de -> 84f9cbd -> 2596830 -> 483b630` | no integrator change over 34 hours | **bitwise identical** |
| `4b26466 -> 71de13f` | `4b26466` touches no `src/` — the user's "true before" | **bitwise identical** |
| `5cc8dec -> 220d928` | `220d928` touches no `src/` | **bitwise identical** |
| `f7d2a31 -> HEAD (working tree)` | the switch instrumentation, default off | **bitwise identical** |

The second row settles the timeline worry directly: **`4b26466`, the last commit before the
`dtau` work, produces the same 1048576 pixels as `71de13f`.** And the `fixoff` arm of
`results/dtau_fix` reproduces it too —

```text
71de13f   nonfin 30109  hot 0.9285  escape 0.2618  collision 0.3632
fixoff    nonfin 30109  hot 0.9285  escape 0.2618  collision 0.3632
5cc8dec   nonfin   178  hot 0.8558  escape 0.2048  collision 0.3477
fixon     nonfin   178  hot 0.8558  escape 0.2048  collision 0.3477
```

— so `DtauMode::FixedPerInterval` is a faithful reconstruction of the pre-fix commit, and the
comparisons made against it were **not** post-fix against post-fix. What was never captured is
the pre-fix *render*; the pre-fix *behaviour* is reachable and now checked against the commit.

---

## THE WALK — 08-25 18:02 TO 08-28 01:10

`strip_uniform.png`, `strip_outcome.png`, `strip_drift.png` — eight rendered states in order,
each with its numbers underneath.

```text
   commit      when       nonfin      hot   escape  bounded   collis   drift p50  |n0|>0.99
  f4084de  08-25 18:02     17032   0.9206   0.6760   0.0541   0.2668   6.159e-01     0.1032
  84f9cbd  08-26 08:22     17032   0.9206   0.6760   0.0541   0.2668   6.159e-01     0.1032
  2596830  08-26 19:38     17032   0.9206   0.6760   0.0541   0.2668   6.159e-01     0.1032
  483b630  08-27 00:16     17032   0.9206   0.6760   0.0541   0.2668   6.159e-01     0.1032
  077b092  08-27 00:33     17032   0.9206   0.6763   0.0541   0.2665   6.159e-01     0.1032
  e53223d  08-27 11:57      5450   0.9257   0.4899   0.1488   0.3604   2.107e+00     0.2664
  4b26466  08-27 13:52     30109   0.9285   0.2618   0.3693   0.3632   2.290e+00     0.2874
  71de13f  08-27 13:02     30109   0.9285   0.2618   0.3693   0.3632   2.290e+00     0.2874
  5cc8dec  08-27 21:19       178   0.8558   0.2048   0.4475   0.3477   3.362e-03     0.3048
  220d928  08-27 21:36       178   0.8558   0.2048   0.4475   0.3477   3.362e-03     0.3048
  f7d2a31  08-28 01:10       178   0.8541   0.2017   0.4505   0.3477   3.021e-03     0.3051

                 pair     flips     frac     moved     frac   chord p50   chord max
     f4084de->84f9cbd         0   0.0000         0   0.0000   0.000e+00   0.000e+00
     84f9cbd->2596830         0   0.0000         0   0.0000   0.000e+00   0.000e+00
     2596830->483b630         0   0.0000         0   0.0000   0.000e+00   0.000e+00
     483b630->077b092       347   0.0003         0   0.0000   0.000e+00   0.000e+00
     077b092->e53223d    418868   0.3995    653458   0.6232   2.743e-01   2.000e+00
     e53223d->71de13f    299211   0.2853    458488   0.4372   0.000e+00   2.000e+00
     71de13f->5cc8dec    348314   0.3322    909089   0.8670   1.113e-01   2.000e+00
     220d928->f7d2a31    145651   0.1389    979642   0.9343   1.242e-02   2.000e+00
```

### THERE IS NO EARLIER CANDIDATE — AND BEFORE 08-25 13:39 THERE IS NO SLICE

**`030de1a` (08-25 13:39) through `483b630` (08-27 00:16) are BITWISE IDENTICAL on this slice** —
34½ hours, the whole scheduler/criterion run of 08-26, zero pixels moved. `0114be4` and
`a320f50` sit inside that bracket and were not rendered individually; they are bounded by
endpoints that agree bitwise, which is not the same as rendering them, and is said that way
rather than claimed as proof.

Going further back stops for two different reasons, and **the mass gate caught both**:

| commit | when | what happens |
|---|---|---|
| `961a313` | 08-24 12:48 | `Chart::Latent` and `decoder::Latent` **do not exist**. No slice to render. |
| `be478e1` | 08-24 13:27 | same — and this is the commit that *creates* the decoded-mass path. |
| `30d713f` | 08-24 13:43 | chart exists, decodes to **`(0.31628, 0.48444, 0.19928)`** — gate fires, run refused. |
| `45e7dcb` | 08-25 03:09 | same masses, same refusal. |
| `030de1a` | 08-25 13:39 | gate passes. Earliest renderable state, and identical to everything through `483b630`. |

`be478e1` cannot be the culprit for a further reason than timing: **it is the commit that
introduced the latent chart's masses at all.** There is no earlier state in which this slice
exists to be broken.

`30d713f` and `45e7dcb` are not "this slice, rendered wrong" — they are a **different physical
system**, decoded before the three corrections of `f4084de` landed. Rendering them would have
produced a plausible-looking panel of the wrong problem, which is exactly the failure the gate
exists to stop.

### THE MASS PATH IS CLEAN AT HEAD, OVER 8388608 SYSTEMS

`examples/mass_audit.rs`, output in `results/output/mass_audit.txt`. Built by
`jitter::copies_with_path` with `evaluate_at`'s own arguments — `evaluate_at`'s next lines are
`integrate_az_opts(c.s, &c.m, ..)` and `energy(.., &c.m, ..)`, so this is the integrator's input
and not a parallel reconstruction of it.

```text
              case  masses (px 0)     max|dm| max|sum m-1|  max|sum p| max|sum m r|    m spread
  config_stability 0.32735,0.42763    3.525e-6   1.110e-16   2.861e-17   2.124e-16     0.000e0
      preset_shape 0.33333,0.33333     0.000e0     0.000e0     0.000e0   1.777e-16     0.000e0
       preset_prho 0.33333,0.33333     0.000e0     0.000e0     0.000e0   5.914e-17     0.000e0
    preset_plambda 0.33333,0.33333     0.000e0     0.000e0     0.000e0   5.914e-17     0.000e0
   preset_shape_pl 0.33333,0.33333     0.000e0     0.000e0     0.000e0   1.777e-16     0.000e0
```

Every pixel, all 8 copies. `max|dm| = 3.5e-6` is the rounding of the expected constants as
typed to 5 dp, not a discrepancy. **Total momentum is zero to `2.9e-17` and the COM is at the
origin to `2.1e-16`** — the two are asserted separately, because zero momentum does not imply
zero first moment and a construction that assumes a COM-centred input returns a drifting system
without one. `m spread` is exactly zero across every footprint, which is what makes
`evaluate_at`'s `copies[0].m` shortcut exact on a configuration chart rather than an
approximation.

**So the mass path is not the bug.** The one thing worth saying about that result is that an
equal-mass control could not have produced it: on a preset every mass error is invisible by
construction, so `config_stability` is the only row in that table carrying information.

### THE HISTORY IS LINEAR — THERE IS NO BRANCH COMPARISON TO MAKE

`git log --all --graph` restricted to `driver.rs` and `outcome.rs` is a **single straight line**
from `b497fa2` to `c7bdece`. No merge commit in `f4084de..f7d2a31` touches `driver.rs`,
`outcome.rs` or `pixel.rs` at all.

* `dtau-step-control` and `overshoot-clamp` are **the same commit**, `f7d2a31` — two names for
  one tip, not two divergent lines.
* `closure-criterion` is `4b26466` and `escape-distance-gate` is `e53223d`; both are **ancestors**
  of `f7d2a31` and both are already **merged into `origin/main`** (`66639b2`, PR #25).
* Only `f7d2a31` itself is outstanding, as PR #26.
* The single branch tip that is *not* an ancestor is `criterion-sweep` (`0e100e5`), whose `src/`
  tree is **identical** to `84f9cbd` — the same work, landed under another name via
  `sweep-onto-main`.
* Local `main` reads `0c070e4` (08-25 18:10). That is a stale local ref, not a structural fact;
  `origin/main` is `66639b2`.

The screenshots, the committed `dtau_fix` renders and this bisect are all on that one line. **The
walk did not need to become a branch comparison.**

### TRIPLE COLLISION IS IMPLEMENTED AND ON `main`

`src/outcome.rs` carries `State::is_triple`, `triple_ejection`, and the `>=2-pair` rule with its
triangle-inequality argument — *"`detail` for a collision mask: the pair index for one pair, `3`
for two or more"*. It is present at `0114be4`, `023c4ce` and `030de1a`, so it long predates the
closure work, and it is in `origin/main`. Not a local prototype.

The first pixel to move at all is at `077b092` — **347 flips, 0.03%**, and `chord p50` exactly
zero, so the shape field does not move at all there. Everything before the escape work is one
state.

### THE LARGEST MOVE IN THE RANGE IS NOT THE `dtau` FIX

`077b092 -> e53223d` — the escape gate's missing distance condition — flips **39.95%** of the
labels against the `dtau` fix's 33.22%, and the two escape commits together take the escape
fraction from **0.676 to 0.262** and the bounded fraction from **0.054 to 0.369**. Those are
correctness changes with their own evidence on record. But the user's screenshots sit at
`4b26466`, *after* both of them — so the state they consider original already contains the
largest label move in the whole range, and any account beginning at the `dtau` work has missed
the bigger half.

---

## THE BLEACHING IS AT `5cc8dec`, AND IT IS A LOSS OF TEXTURE, NOT A GAIN IN LIGHTNESS

The obvious statistics say the opposite of the report, and they are the wrong statistics.
Median OKLab lightness **falls** across the walk (0.796 -> 0.843 -> **0.614**) and the strictly
white population (`L>0.80, chroma<0.030`) falls monotonically 0.0179 -> 0.0150. The panel gets
*darker and more saturated*. What the eye reads as "the pale regions grew" is those regions
going **flat**: the fine mottled striation inside them is replaced by a near-uniform wash.

The statistic that matches the eye is local contrast — the 5x5 standard deviation of OKLab
lightness and chroma, read off the panel on the fixed ramp:

```text
   commit  magenta   L p50  Lvar p50  Cvar p50   C coh-x  C coh-y
  f4084de   0.0266  0.7957   0.05040   0.01542    0.6127   0.6567
  e53223d   0.0287  0.8425   0.04557   0.01726    0.6244   0.6641
  4b26466   0.0288  0.8432   0.04503   0.01735    0.6208   0.6617
  5cc8dec   0.0002  0.6143   0.02986   0.00557    0.7042   0.7265   <- here
  f7d2a31   0.0002  0.6116   0.02759   0.00461    0.7049   0.7277
```

`Cvar` collapses **3.1x** at `5cc8dec` and nowhere else, with a further 1.2x at `f7d2a31`. So the
user's read of *which* change is right.

**And what left was the incoherent component.** *Amplitude cannot tell a small real signal from
noise; coherence can.* The lag-1 neighbour correlation of chroma **rises** 0.62 -> 0.70 in the
same step in which the local contrast falls 3.1x. A texture that vanishes while what remains
becomes *more* spatially coherent was noise. Measured on the 55.5% of the frame lying at least
8 px from any pre-fix magenta pixel, `Cvar` falls by the same 2.9x — so this is not the halos
being cleaned up, it is the whole field.

---

## THE 10:1 RATIO: 2.87% WAS *UNUSABLE*, ~93% WAS *INACCURATE*

The framing that makes it look wrong is `nonfin`, which counts only total failure.

```text
before (4b26466):  nonfin 30109 (2.87%)   hot 0.9285   drift p50 2.2905
after  (5cc8dec):  nonfin   178 (0.02%)   hot 0.8558   drift p50 0.0034
```

Median energy drift before the fix is **2.29 — 229% of the total energy of the system**, over the
median pixel of the frame. The pre-fix field was not 2.87% broken with a good remainder; it was
NaN on 2.87% and quantitatively meaningless on most of the rest. Against a **682x** fall in
median drift, 33% of labels changing is a small number, not a large one. `strip_drift.png` shows
it without arithmetic: panels 1–5 are bright over most of the frame, panel 6 onward dark.

**Where the flips are** — enrichment of the 348314 flips inside the pre-fix magenta set, dilated
by `r`:

```text
71de13f->5cc8dec   magenta 30109   flips 348314 (0.3322)
                   r0 2.22   r1 1.91   r2 1.77   r4 1.62   r8 1.48   r16 1.36
```

Enriched 2.2x on the magenta pixels themselves and still 1.36x at 16 px, but the magenta set is
2.87% of the frame, so **about 94% of the flips are outside it**. Spread, not concentrated — and
given the drift map, that is what a correct fix looks like here rather than evidence against one.

---

## DOES `hot` FALL WHEN `eta` IS HALVED? YES — BUT THE LABELS DO NOT CONVERGE

256^2, HEAD, `n_sync = 32` held fixed across the ladder, because scaling it with `eta` would
compare different discretisations.

```text
       eta   nonfin      hot  drift p50   escape  bounded   collis   chord p50   chord max
   1.00e-2        8   0.8750   9.027e-3   0.2016   0.4516   0.3468         --          --
   5.00e-3        5   0.7643   3.619e-4   0.1787   0.4653   0.3560   5.529e-3     2.000e0
   2.50e-3        3   0.6339   1.593e-5   0.1546   0.4754   0.3700   7.436e-4     2.000e0
   1.25e-3        1   0.4871   7.281e-7   0.1171   0.4862   0.3967   1.249e-4     2.000e0
```

`hot` falls 0.875 -> 0.487 and median drift falls **12400x** over an 8x refinement — an observed
order near 4.5, so this is truncation error and **not** a wrong equation. The `hot` fraction is
high because the threshold `1e-6` sits far below where this slice's bulk lives, not because the
equation is wrong.

**But the classification is not converged.** The escape fraction moves monotonically
`0.2016 -> 0.1171` — a 42% relative change over 8x refinement with no sign of settling — and
`chord max` is **2.000, antipodal, at every rung**. The shape field itself converges (`chord p50`
falls ~6x per halving, consistent with the measured order 2.08); the *labels* do not. At horizon
50, sitting at the f64 predictability limit of ~52 at `lambda = 0.7`, the outcome image on this
slice is a picture of the discretisation as much as of the physics, and no code change reaches
that.

---

## REPRODUCE

```sh
RES=1024 OUT=<dir> bash results/timeline/harness/run.sh          # the walk
python3 results/timeline/harness/analyse.py <dir> 1024 <commits> # the tables
python3 results/timeline/harness/strip.py   <dir> 1024 384 <commits>
cargo run --release --example switch_study 384 <dir>             # the switch study
```

`run.sh` builds in its own `git worktree`; nothing in the working tree is reverted or checked
out. **The harness writes to an argument, never to `results/`** — see the standing rule this
project already has for that.
