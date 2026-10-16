# Chazy final-motion classification, as of time t

`src/chazy.rs`, `examples/chazy_probe.rs`, `examples/chazy_map.rs`, `tests/chazy.rs`.

**The question.** The payload carries an *event* (escape of body k, collision of pair k, bounded at
the horizon), but not the *final motion* in Chazy's (1922) sense. Can the final motion be given at
every playhead, as an approximation that holds up to the current time, and how close does it get
to the literature's standard way of deciding that a three-body interaction is over?

**The answer, measured on 1200 trajectories over two seeds:**

- **It never committed to a wrong class on a bound (E < 0) system, at any playhead, in either seed.**
- **It commits as early as the gold standard at its looser tolerance, and earlier than at its
  usual one.** Median time to the first correct label on E < 0 systems is 51.6 / 40.8 time units
  (seeds 1000 / 50000). Fewbody, the standard (Fregeau et al. 2004 §3.3.5), takes 69.6 / 69.6 at
  its δ = 1e-5 and 52.0 / 42.8 at δ = 1e-3.
- **It leaves fewer trajectories undecided at T = 1000:** 16 / 18 against Fewbody's 30 / 27.
- **The shape-settling gate does not earn its place.** Requiring the closure chord below τ = 0.01
  bought **zero** precision and cost a median 35–36 time units of delay (86.4 / 77.2 against
  51.6 / 40.8). The shipped default has it **off**, so the classifier needs no previous-shape latch,
  only `I_max`.
- **The Lagrange–Jacobi certificate holds.** No E > 0 trajectory was unresolved by the end of the
  truth run in either seed (none can stay bounded), and no E < 0 trajectory with a resolved truth was
  ever given a committed label that disagreed with it, H included.

The negative and unresolved results are below, and they are as much the finding as the table.

---

## 1. What the classifier does

Every snapshot, from the live state plus one running latch (`I_max`, the largest moment of inertia
seen so far). Nothing is re-run and nothing after `t` is used, so it passes the CLAUDE.md
live-playhead test.

| label | condition at the playhead | committed? |
|---|---|---|
| **H** | all three pair separations > `h_far R0`, all receding, every `K/U > 1 + eps` | yes |
| **HP**, **P** | as H, with one or two / all three pairs inside the parabolic band `\|K/U - 1\| <= eps` | yes |
| **HE** | the body most separated from its pair is at hierarchy ratio > `h_far`, receding, `K/U > 1 + eps`; the pair is bound (`K < U`) | yes |
| **PE** | as HE with the escaper inside the band | yes |
| **B-stable** | bound hierarchy that is Mardling–Aarseth (2001) stable | yes |
| **excursion** | the escaper is far and receding but bound to its pair: it will fall back | no |
| **transient** | interacting with E > 0: Lagrange–Jacobi, `I'' = 4E - 2V > 4E`, so it cannot stay bounded | no |
| **OS-candidate** | interacting, E < 0, and `I_max / I > 25`: has been far out and come back | no |
| **B-so-far** | interacting, E < 0 | no |

Defaults (`ChazyCfg::default()`): `gate = Hierarchy`, `h_far = 5`, `tau = INFINITY` (settling off),
`eps = 0.05`, `os_ratio = 25`, `mardling_aarseth = true`.

Two choices are physics rather than tuning:

- **The pair's binding is a sign, not a band.** A bound pair at pericentre has
  `K/U = (1 + e)/2`, which for `e > 0.9` sits inside a parabolic band around 1. The first cut read
  such pairs as marginal and called the system a two-stage H. `tests/chazy.rs::an_eccentric_pair_at_pericentre_is_still_bound`
  holds it, and reverting the sign test makes it fail (checked).
- **Parabolic is a band, not a class.** PE, HP and P need an exact zero of a two-body energy and
  have zero measure. A finite run cannot certify them. They are reported as labels so the map can
  show where the margin is thin, and scored with their hyperbolic neighbour (`Motion::coarse`).

## 2. The gold standard

Fewbody (Fregeau, Cheung, Portegies Zwart & Rasio 2004, MNRAS 352, 1, §3.3) is the standard
small-N scattering code's answer to "is this interaction over". Its classification pairs the bound
top-level nodes of smallest semimajor axis until none is bound to another. Its termination test
(§3.3.5) is all four of:

1. every pair of top-level nodes recedes;
2. the tidal perturbation `F_tid/F_rel` on every binary from the other nodes is below `δ`
   (the paper's runs use `δ = 1e-5`);
3. every tree is dynamically stable by Mardling & Aarseth (2001);
4. the top-level system has positive energy.

`chazy::fewbody_resolved` is that test for three bodies, with the planar Mardling–Aarseth criterion
`R_p^out / a_in > 2.8 [(1 + q)(1 + e_out)/sqrt(1 - e_out)]^{2/5} (1 - 0.3 i/180)`, `i` in {0, 180}.
It is the reference, not a candidate, but its tidal isolation is also offered as a gate
(`Gate::Tidal`) so the two tests can be compared inside one classifier.

Two naive baselines sit beside them as floors: `energy` (the tightest pair is the binary; the sign
of the third body's two-body energy decides) and `distance` (the farthest body beyond `10 R0` has
escaped).

## 3. The truth

The same trajectory, continued after the evaluation horizon in 10-unit chunks until the motion is
unambiguously asymptotic (hierarchy ratio or minimum separation above 1e3) or `t = 10000`, then
classified by **Chazy's own definition**: the growth exponent `beta = d ln r / d ln t` of each pair
separation over the second half of the run. All three `> 0.5` at minimum separation `> 50` is
H/HP/P (`beta > 0.85` hyperbolic). At hierarchy `> 100`, a bound pair (Kepler energy, which the
tidal ratio `< 2e-6` makes rigorous) with an escaper exponent `> 0.85` is HE, `0.5–0.85` PE.
Otherwise unresolved, or B-stable if Mardling–Aarseth stable at the end.

**The first truth was wrong, and the probe said so.** It measured the pair's own separation
exponent and required it under 0.3. An eccentric binary's separation at two arbitrary phases
differs by up to `(1+e)/(1-e)`, which reads as growth, so 40 genuine HE trajectories were
"unresolved" — and every method, Fewbody included, committed on all 40. Replaced with the Kepler
energy test, the unresolved count fell 49 → 13 and no method except `distance` commits on any of
them.

**PE truth is weak and is not leaned on.** A slightly hyperbolic escape still reads `beta < 0.85` at
`t ~ 1e3`. The PE/HE disagreements (4–5 per seed, every method) are this, and they are scored
coarse.

## 4. Results

`probe_n600_seed1000.txt`, `probe_n600_seed50000.txt`: 600 trajectories each. Latent decode,
`z_alpha, z_beta ~ N(0, 1.5^2)`, `z_mu ~ N(0, 1)`, and momenta `z_q ~ N(0, s^2)` with `s` = 0 for
half the sample, 0.3 for a quarter and 1.0 for a quarter. Heggie, `eta = 1e-2`, `r_coll = 0`
(regularised through binary collisions), snapshots every 0.4 (the escape rule's closure window).
0 of 1200 failed.

```
truth                   seed 1000   seed 50000
  H           E>0            270          274
  HE          E<0            279          290
  HE          E>0             30           20
  PE          E<0              8            4
  unresolved  E<0             13           12
```

### 4.1 E < 0, correct% / wrong% at playhead T (seed 1000; seed 50000 agrees within 6 points)

```
method                         T=5          T=10          T=20          T=50         T=100         T=200         T=500        T=1000
chazy h5 tau.01         1.4/0.00      1.7/0.00      4.2/0.00     13.6/0.00     31.4/0.00     53.0/0.00     74.6/0.00     87.8/0.00
chazy h5 tau.05         2.1/0.00      3.1/0.00     15.3/0.00     31.0/0.00     53.0/0.00     68.3/0.00     81.9/0.00     89.9/0.00
chazy h5 no-settle      4.5/0.00     10.5/0.00     24.7/0.00     44.3/0.00     57.1/0.00     70.0/0.00     82.2/0.00     89.9/0.00   <- default
chazy h3 no-settle      5.9/0.00     11.8/0.00     26.1/0.00     44.9/0.00     57.5/0.00     70.4/0.00     82.2/0.00     89.9/0.00
chazy h10 tau.01        1.4/0.00      1.7/0.00      4.2/0.00     13.6/0.00     31.4/0.00     53.0/0.00     74.6/0.00     87.8/0.00
chazy-tidal d1e-5       1.0/0.00      1.7/0.00     10.1/0.00     33.8/0.00     52.6/0.00     68.6/0.00     81.9/0.00     89.5/0.00
chazy-tidal d1e-3       2.8/0.00      9.1/0.00     24.0/0.00     43.9/0.00     57.1/0.00     70.4/0.00     82.2/0.00     89.9/0.00
fewbody d1e-5           1.0/0.00      1.7/0.00     10.1/0.00     33.8/0.00     52.6/0.00     67.9/0.00     81.5/0.00     89.5/0.00
fewbody d1e-3           2.8/0.00      9.1/0.00     24.0/0.00     43.9/0.00     56.4/0.00     69.7/0.00     81.9/0.00     89.9/0.00
naive energy            7.0/0.00     12.9/0.00     27.2/0.00     45.6/0.00     56.8/0.00     69.7/0.00     81.9/0.00     89.9/0.00
naive distance          8.0/0.00      8.7/0.00     23.0/0.00     50.5/0.00     68.3/0.00     85.0/0.00     94.4/0.00     97.2/0.00
```

### 4.2 Time to the first correct commit, E < 0 (p50 / p90; seed 1000 | seed 50000)

```
                         p50     p90  never  flips            p50     p90  never  flips
chazy h5 tau.01          86.4   485.2    26     0     |       77.2   402.0    26     0
chazy h5 no-settle       51.6   426.0    16     1     |       40.8   382.0    18     0     <- default
chazy-tidal d1e-5        72.4   462.4    26   144     |       69.6   416.0    25    59
fewbody d1e-5            69.6   462.4    30     0     |       69.6   414.4    27     1
fewbody d1e-3            52.0   458.4    29     0     |       42.8   357.2    27     1
naive energy             49.6   457.6    24     6     |       39.6   338.0    21     2
naive distance           44.4   235.6     2   285     |       35.6   181.2     0   289
```

`never` is resolved-truth E < 0 trajectories with no correct label by T = 1000. `flips` is coarse
label changes after the first commit, summed over all trajectories. `wrong-first` (committing to a
wrong class before the right one) is **0 for every method in both seeds**.

### 4.3 What it reads as

- **The default matches Fewbody at δ = 1e-3 and beats it at δ = 1e-5** on every horizon of
  4.1, at zero wrong commits. That is the "as good or nearly as good" bar, met.
- **The settling gate is pure delay.** `h5 tau.01` and `h10 tau.01` are identical on every E < 0
  row, because on these runs it is the settling that binds, not the hierarchy ratio. Settling is the escape rule's closure chord, and a
  hierarchical system's shape wobbles at the binary period with amplitude `~ a/d`, so the chord
  falls below 0.01 only once the escaper is of order a hundred binary separations away. The energy and
  hierarchy tests had decided long before.
- **The tidal gate flickers and the hierarchy gate does not.** `chazy-tidal` flips 144 times on 14
  trajectories (59 on 7 in the second seed). It is Fewbody's own isolation test, but Fewbody flips
  0 / 1 times: it commits only when a binary exists *and* is isolated, while the tidal gate inside the
  classifier lets the two-stage-H branch fire when a wide, barely-bound pair dips in and out of
  binding. With the hierarchy gate that branch is reachable only at ratio > 5, and flips fall to
  0–1.
- **The naive energy test is nearly as good on this sample, and that is reported, not hidden.**
  0 wrong, p50 49.6 / 39.6. It flips 6 / 2 times where the default flips 1 / 0, and it cannot
  produce PE, HP, transient, OS-candidate or excursion. The gates cost almost nothing and buy the
  stability and the labels. On a sample with more long, close interplay the snapshot test would
  latch transients — the 0-of-895 persistence result on `deep interior` is that failure.
- **Distance is fastest and wrong.** 0.85–1.4% wrong at late horizons, 285 flips, and it commits
  on unresolved trajectories at every playhead.

### 4.4 Not perfect, and where

- **Seed 1000, E > 0, T = 10: one wrong commit** (0.33%, 1 of 300) under the no-settle defaults and
  `h3`, withdrawn by T = 20. **Seed 50000: one temporary commit on an unresolved trajectory** at
  T = 20 and T = 100 (0 at the other playheads). Fewbody did neither. That is the price of the
  faster gate, and it is one trajectory in 600 per seed.
- **14–18 resolved HE trajectories are still B-so-far at T = 1000**, and 6–10 more are
  "excursion". They escape later than the playhead. That is the correct finite-time answer, not an
  error.
- **OS-candidate is not validated.** 2–3 per seed later escape as HE: the label means "has been far
  out and come back", not a prediction of oscillatory motion. Chazy's OS class is non-empty and
  believed to have zero measure; no finite run can certify it, and this one does not try.
- **B-stable never fired** on 1200 random starts. No Mardling–Aarseth-stable triple formed. It is
  covered by a constructed fixture only (`tests/chazy.rs`).
- **The sample is a choice.** The E > 0 half is mostly prompt dispersal and every method scores
  95–100% on it, so the interesting numbers are the E < 0 rows. A slice through a chaotic chart
  would weight the long-lived population more heavily.

### 4.5 The class mix over time (seed 1000, default classifier, % of all 600)

```
     T         H        HP         P        HE        PE  B-stable transient OS-candidate excursion  B-so-far
     5      43.8       0.2       0.0       6.0       0.0       0.0       2.2       0.3       2.0      45.5
    10      44.8       0.2       0.0       9.2       0.0       0.0       0.8       0.5       5.8      38.7
    20      45.0       0.0       0.0      16.3       0.0       0.0       0.5       0.3       5.7      32.2
    50      45.0       0.0       0.0      25.2       0.5       0.0       0.5       0.2       6.8      21.8
   100      45.0       0.0       0.0      31.5       0.7       0.0       0.2       0.7       6.5      15.5
   200      45.0       0.0       0.0      38.0       0.5       0.0       0.0       0.5       4.8      11.2
   500      45.0       0.0       0.0      44.0       0.3       0.0       0.0       0.3       4.0       6.3
  1000      45.0       0.0       0.0      47.7       0.3       0.0       0.0       0.3       2.7       4.0
```

This is the series the time-resolved view exists to show: bound-so-far draining into HE, transient
draining by T = 200 as Lagrange–Jacobi says it must, and the excursions as the population in
flight.

## 5. Maps

`chazy_map.rs` marches every pixel once and classifies it at every snapshot with the shipped
default; the strip is the class at `t = 10, 50, 200, 800`, left to right. Same trajectories in
every panel, so a pixel changing colour is the label evolving.

Legend: **H** `#3462E0`, **HP** `#7FA0F0`, **P** `#C0D0FA`, **HE** `#DE2D2D`, **PE** `#F29620`,
**B-stable** `#2EBC4E`, **B-so-far** `#141418`, **transient** `#E034C6`, **OS-candidate** `#F0DE32`,
**excursion** `#30C8DC`, **unresolved** `#4E4E4E`.

![preset_shape_h1](chazy_map_preset_shape_h1_n48.png)

`preset_shape_h1`, 48², free fall (E < 0 everywhere), 0 of 2304 failed:

```
       class     t=10     t=50    t=200    t=800
          HE   11.89%   32.99%   67.14%   90.49%
          PE    0.69%    0.95%    0.52%    0.22%
    B-so-far   77.17%   46.18%   17.06%    3.39%
OS-candidate    2.17%    1.39%    2.26%    1.04%
   excursion    8.07%   18.49%   13.02%    4.86%
```

The sea chart: the regular island and the smooth basins decide first, the fractal core drains
from bound-so-far through excursions into HE, and what is still black or cyan at `t = 800`
(3.4% and 4.9%) is what the map has not decided. No H, no transient: free fall has E < 0, and the
certificate holds pixel by pixel.

![latent_shape](chazy_map_latent_shape_n48.png)

`latent_shape`, 48², 0 of 2304 failed: 81.25% H and 18.75% HE from `t = 50` on, and 2.99%
transient at `t = 10` that has resolved by `t = 50`. A tame chart, the control for the one above.

**What the maps do not show.** HE is one colour, so which body escaped is not drawn; the event
class already carries it, and composing the two is the obvious next panel. And 48² is a probe, not
a render: the per-pixel cost to `t = 800` is about 0.45 CPU-seconds on the sea chart.

## 6. For the Principia implementation

- **Derived, not stored.** Everything but `I_max` is a function of the live state, the masses and
  `E_0`, which the read side already has (`read_side.wgsl` carries `E_0` and the live `r, p`). So the
  class is a generated read-side accessor beside `energy_drift`, recomputed at every playhead.
- **One latch.** `I_max` as a running max, the mirror of `d_min`. A log-quantised value fits the
  reserved descriptor bits 10–14 (32 half-octave levels up to `2^16`), leaving bit 15 for the open
  departure flag. With settling off there is no previous-shape latch to store.
- **Parity.** The class flips on the sign of a two-body energy, so it is Tier N (numerical), not
  Tier B: CPU and GPU agree to tolerance, and pixels within the band of a sign change read as the
  marginal labels rather than as confident ones.

## Reproduce

```bash
cargo run --release --example chazy_probe -- 600 1000 0.05 25 0 0 10000 1000    # ~3 min on 2 cores
cargo run --release --example chazy_probe -- 600 1000 0.05 25 0 0 10000 50000
cargo run --release --example chazy_map -- preset_shape_h1 48 results/chazy
cargo run --release --example chazy_map -- latent_shape 48 results/chazy
cargo test --release --test chazy
# event outcome per pixel with its time, at the production r_coll and closure escape rule:
cargo run --release --example outcome_map -- preset_shape_h1 48 800 results/outcome
```

`chazy_map` also writes the per-pixel class, structure, energy and per-body dumps described in its
module header, so the stain previews can be recoloured without re-integrating. `outcome_map` writes
each pixel's packed `(state, detail)` and `t_end`; a playhead shows a pixel's outcome once
`t_end <= playhead`. **Its escape `detail` is the reference rule's lowest-firing body index, not
necessarily the body that leaves** (see CLAUDE.md, "the escaping-body label is the lowest firing
index"): on `latent_shape` it reads body 0 where the Chazy escaper is body 2.

## References

- J. Chazy (1922), Ann. Sci. École Norm. Sup. 39, 29 — the final-motion classes.
- J. M. Fregeau, P. Cheung, S. F. Portegies Zwart, F. A. Rasio (2004), MNRAS 352, 1, §3.3 —
  Fewbody's classification and termination test ([arXiv:astro-ph/0401004](https://arxiv.org/abs/astro-ph/0401004)).
- R. A. Mardling, S. J. Aarseth (2001), MNRAS 321, 398 — hierarchical triple stability.
- E. M. Standish (1972), Celest. Mech. 6, 352 — rigorous sufficient conditions for return; the
  certificate this classifier does **not** implement, and the natural next one to add.
