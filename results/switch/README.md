# Reference-body switching — why THIS slice and not the presets

`config_stability` carries far more sharp discontinuity in its drift map than `preset_prho`
does. A chaotic field gives smooth gradients and fractal boundaries; a **discrete decision
boundary** gives a step. AZ chooses its reference body as the one *not* in the longest side —
a discrete choice — and every change of it re-derives the Levi-Civita registration
mid-trajectory. This slice has unequal masses `(0.32735, 0.42763, 0.24502)` and asymmetric
geometry; the presets sit at `z0 = 0` with equal thirds.

**One correction to the framing before the numbers.** The reference is chosen **once per sync
boundary**, not per RK4 step (`driver.rs`, the `for kk in 0..n_sync` head). So "choose it once
per sync interval rather than per step" is already the shipped behaviour, and switch times are
quantised to `t_max / n_sync` by construction — no finer statement is available from this data.

Run at 384^2, `keep_drift_hist` on so the drift series and `AzOut::refs` share one cadence.
Full output in `results/output/switch_study.txt`.

---

## 1. THE DISCONTINUITY CENSUS DOES NOT SEPARATE THE SLICES THE WAY IT WAS QUOTED

Measured on the committed 1024^2 drift panels, luminance gradient over adjacent pixels, magenta
excluded:

```text
                                        pairs      >100      >200      >400
config_stability_fixon_drift.png      2094393    0.0479    0.0026    0.0000
config_stability_fixoff_drift.png     1994370    0.0544    0.0016    0.0000
preset_prho_fixon_drift.png           2095100    0.0010    0.0000    0.0000
preset_shape_fixon_drift.png          2094996    0.0144    0.0032    0.0000
preset_plambda_fixon_drift.png        2095104    0.0065    0.0017    0.0000
preset_shape_pl_fixon_drift.png       2095100    0.0071    0.0024    0.0000
```

The **ratio reproduces**: `config_stability` runs **48x** `preset_prho` at `>100`, against the
quoted 50x at `>200`. The absolute percentages do not — 4.79% and 0.10% here against 6.19% and
0.12% — so the threshold is defined differently, and **`>400` is empty on every committed
panel**, so the 0.41% figure does not survive a luminance gradient.

**And the contrast is not "this slice against the presets".** `preset_shape` carries *more*
gradients above 200 than `config_stability` does (0.0032 against 0.0026). It is `preset_prho`
that is unusual — and it is unusual for a structural reason already on record: every pixel of it
is the *same triangle* at a different initial velocity.

---

## 2. SWITCH COUNT SEPARATES THEM 21x — AND DOES NOT ORDER THEM

```text
             slice   masses         mean sw   p50   p90   max   never sw   |grad|>100
  config_stability   (.327,.428,.245)  5.717     3    14    31     0.1069       0.0537
      preset_shape   (1/3,1/3,1/3)     1.023     0     4    25     0.8117       0.0212
    preset_plambda   (1/3,1/3,1/3)     1.463     1     3    24     0.0023       0.0003
       preset_prho   (1/3,1/3,1/3)     0.270     0     2    21     0.8679       0.0006
```

`config_stability` switches **21x** as often as `preset_prho` — 5.7 switches over 32 boundaries
against 0.27, and only 10.7% of its pixels never switch against 86.8%.

**But the count alone does not predict the discontinuity density.** `preset_plambda` switches
*more* on average than `preset_shape` (1.463 against 1.023) and carries **70x fewer** sharp
gradients. Its switches are all at the same boundary — `t_first` p10 = p50 = p90 = 0.406 — so
neighbouring pixels switch *together*, and a switch every pixel makes at the same time draws no
edge.

**The quantity that does order them is spatial variability of the switch history** — the
fraction of neighbour pairs whose `(count, first-switch time)` differ:

```text
  config_stability   58.4%     |grad|>100  0.0537
      preset_shape   13.6%                 0.0212
    preset_plambda    5.1%                 0.0003
       preset_prho    1.6%                 0.0006
```

Monotone across the top three, with the bottom two both at the floor. **It is not how often the
reference switches; it is how often two neighbouring pixels switch differently.**

---

## 3. THE ALIGNMENT TEST — CONFIRMED, WITH A CONTROL THAT KILLS THE NULL

For each adjacent pixel pair: is a large drift step more likely when the two pixels have
different switch histories? Reported against the matched-history population, and against a
**shifted control** in which each pixel's history is taken from a pixel 37 rows and 53 columns
away — same spatial statistics, alignment destroyed. Both fields are spatially smooth, so
without that control any two maps agree somewhat.

```text
                    |grad|>200      P(step|differs)   P(step|matches)   enrichment   SHIFTED
  config_stability       0.0037              0.0052            0.0015        3.43x     0.65x
      preset_shape       0.0036              0.0169            0.0015       10.90x     2.62x
    preset_plambda       0.0000              0.0006            0.0000          n/a     0.00x
       preset_prho       0.0000              0.0011            0.0000          n/a     0.00x
```

`config_stability` reads **3.43x**, against the 1.69x the `t = 0` reference map gave — the
proxy understated it by half, as predicted, because a trajectory to horizon 50 crosses many
boundaries the snapshot cannot see. **The shifted control reads 0.65x, below 1**, so none of
that enrichment is smoothness. `preset_shape` reads 10.90x with a shifted control of 2.62x, so
roughly a fifth of its enrichment is smoothness and the rest is real.

---

## 4. DRIFT ARRIVES AT SWITCHES — THE PAIRED INCREMENT

`|drift[k] - drift[k-1]|` at the boundaries where the reference **changed**, against the
boundaries where it **held**, within the same trajectory. A correlation between two maps can be
produced by both tracking a third thing; a paired increment across the switch cannot.

```text
             slice     switch med     hold med      PAIRED ratio p50    frac ratio>1
  config_stability       3.988e-9     1.169e-10                14.357          0.8274
      preset_shape       7.686e-10     8.815e-11               416.619         0.9885
    preset_plambda       1.605e-14     1.700e-15                 9.273         0.9780
       preset_prho       9.061e-12     6.045e-16             12927.000         0.9895
```

Paired within each pixel, the switch increment exceeds the hold increment on **82.7% to 99.0%**
of pixels in all four slices. The mechanism is present everywhere; what differs between slices
is how often, and how *incoherently*, it fires.

**The confound, stated.** A switch happens when the longest side changes, which is also when the
triangle is deforming fastest — so "drift is larger at switches" could be the encounter rather
than the re-registration. The argument against reading it that way is `preset_prho`: the
smooth, near-quiescent slice, where switch increments still run **four orders** above hold
increments (`9.06e-12` against `6.05e-16`) with nothing dramatic happening. That is suggestive,
not decisive, and a clean separation would need the boundaries matched on `d_min`.

---

## 5. THE COLLISION-CADENCE CANDIDATE POINTS THE OTHER WAY

The concern on record is that AZ tests collision separation every RK4 step while the reference
tests once per macro step. Measured on this slice at the shipping settings, the median physical
RK4 step is **~9.6e-3** (`t_end` over per-copy substeps), against the reference's
`dtMacro = 0.002`. So the reference samples collision about **five times more often in physical
time** than prin-rs does, not less. Whatever else is different, this is not a case of prin-rs
catching dips the reference misses.

---

## PANELS

`<slice>_switches.png` — switch count, inferno, scaled to that slice's own max.
`<slice>_tfirstswitch.png` — first-switch time; **magenta means never switched**, which is the
large-scale structure worth looking at: the never-switch set is exactly the coherent ribbons.
`<slice>_drift.png`, `<slice>_outcome.png` — for comparison, same run.

## WHAT IS NOT DONE

Remedies are costed nowhere here, deliberately. Hysteresis on the switch, or carrying the
registration continuously rather than re-deriving it, are the two that survive the "already
per-sync-boundary" correction — but the numbers come first.
