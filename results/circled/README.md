# What the circled ICs have in common

The regions marked on `results/closure/config_stability_stop0_uniform.png` (1024^2) are pale,
low-chroma wedges with straight edges and magenta speckle, sitting on the boundaries between the
coloured ribbons.

**The initial conditions need no integration**, so every number here is exact and free:
`grid::decode_state` gives masses, positions and velocities per pixel in milliseconds.
`examples/circled_ics.rs`, output in `results/output/circled_ics.txt`.

## The populations, and why there are three

`circled_mask.png` — yellow outlines are the hand-digitised circles, green is the
**pale/low-chroma** mask selected by property. The digitisation is dumped so it can be checked
against the photo rather than trusted.

The hand circles cover **25.4% of the frame** — far larger than the features inside them — so
that row is mostly background and its numbers are diluted. Read the **magenta** (3019 px) and
**dense-pale wedge** (1935 px) rows. They agree with each other, which is what says the finding
is about the field rather than about where a circle was drawn:

```text
  91.0% of the magenta lies inside a circle, against 25.4% expected by area  -- 3.58x
  56.7% of the pale mask lies inside a circle                                -- 2.23x
```

## What they have in common

```text
                       frame       magenta    dense wedges     (p50 unless stated)
  d(0,1)              1.5798        1.3713          1.3606     the tight pair
  d(0,2)              1.9883        1.8043          1.8828
  d(1,2)              1.7871        2.0836          1.9404     the wide pair
  d(1,2) p10          1.4211        1.8384          1.8377     <- the bottom is CUT OFF
  aspect dmax/dmin     1.428         1.574           1.579     more hierarchical
  alpha               0.8232        0.9392          0.9448     tighter inner pair
  beta                1.4849        1.8547          1.8032
  |Lz|                0.1662        0.1020          0.1061     ~60% of the frame median
  K/|U|               0.3761        0.3671          0.3677
  tightest pair (1,2) 0.2889        0.0007          0.0005     <- 413x DEPLETED
  AZ reference body 0 0.2902        0.7671          0.8656     2.6x
```

**One geometric fact says most of it.** These are **hierarchical** configurations ordered
`d(0,1) < d(0,2) < d(1,2)`, in which **bodies 1 and 2 start far apart** — `d(1,2)`'s 10th
percentile rises from 1.42 to 1.84, so the small-separation tail is essentially absent. With
`m = (0.32735, 0.42763, 0.24502)` that is the **heaviest** and the **lightest** body as the wide
pair, and the middle-mass body 0 closest to the heaviest. AZ then takes body 0 as its reference
— the body not in the longest side — so it regularises the two *shortest* pairs and leaves the
widest unregularised, which is correct behaviour and not itself a fault.

The near-hard version of that statement is the argmax form: **`tightest == (1,2)` occurs on
0.07% of the magenta against 28.9% of the frame** — a 413x depletion, close to an exclusion
rather than a shift.

Alongside it, consistently in both feature populations: **lower angular momentum** (|Lz| ~60% of
the frame median), a **tighter inner pair** (larger `alpha`, `cos alpha = ||rho~||`), and a
slightly lower virial ratio.

## What it is NOT — and this was the leading candidate

**Not a near-degeneracy.** If these were the pixels where AZ's `argmax` is a coin flip, the tie
statistics would be enriched. They are **depleted about 2.5x**: `d[2nd]/d[longest] > 0.95` is
0.0825 in the magenta against 0.2132 elsewhere, and the tightest-pair tie is 0.0391 against
0.1701. These configurations sit *well away* from the reference-body switching boundary at
`t = 0`, not on it.

`ic_class.png` is the direct test: the AZ reference body and tightest pair over the whole frame,
hue by reference and lightness by tightest pair. It is a **smooth six-sector pinwheel meeting at
one point** — the straight edges the eye picks up in the chart plane are *this*, and it does
**not** draw the circled wedges. The sectors are the polar structure already on record for this
slice; the wedges cross them.

## The limit of the finding, stated

Every row above is a **shift of overlapping distributions**, not a separation — the p10–p90 bands
overlap the frame's everywhere except `d(1,2)`'s lower tail. And the class is nowhere near
sufficient: `P(magenta | reference=0 and tightest=(0,1))` is **0.0070** against a base rate of
0.0029, a 2.4x lift on a class holding 22% of the frame.

So the initial conditions **constrain** where the artefact can appear and do not **draw** it. The
fine structure inside those constraints comes from the dynamics, not from the ICs — which is the
opposite of what the straight edges first suggested, and is why the class map was rendered rather
than argued about.
