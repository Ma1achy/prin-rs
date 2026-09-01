# Stages 2, 3 and 4 are NOT BUILT, and this file says so rather than leaving a gap

Stage 0 and Stage 1 ran. Chain coordinates, TTL and IAS15 did not, because each is a new
integrator or stepper and none could be written and *validated* in the time available. An
unvalidated integrator run overnight produces numbers that have to be withdrawn, which is worse
than no numbers -- this project has withdrawn three results that way already.

**Stage 2's gate is now OPEN — Stage 0 ran and adjudicated.** The saturation guard reads LIVE at
both precisions (16368 / 16347 distinct drift values, distribution overlap 0.0004), so the
comparison has a subject. `far`'s AZ advantage does not merely collapse at f32, it **inverts**:

    f64: gain p50 +0.890  frac az better 1.0000     (AZ better, every pixel)
    f32: gain p50 -0.768  frac az better 0.0000     (HEGGIE better, every pixel)

By the brief's decision rule that is conditioning confirmed and chain coordinates justified.

**The mechanism proposed in PREDICTIONS.md is refuted, and it was mine.** I argued `Gamma*`'s
degree-six algebra would be the fragile one at low precision, so AZ's advantage should WIDEN.
The measurement says AZ is the precision-sensitive arm: f64 -> f32 costs AZ a factor of 4e8
(2.8e-13 -> 1.2e-4) against Heggie's 9e6 (2.2e-12 -> 2.0e-5), about 44x more, and AZ picks up
25 not-data pixels where Heggie has none. **No replacement mechanism is offered here** -- the
standing instruction is not to reach for a second explanation in the same breath as the first
one failing.

What each stage needs, so the next session starts from a spec and not from scratch:

- **Chain coordinates** (Mikkola & Aarseth 1993). Two inter-particle vectors for three bodies.
  A round-off fix, so it must be run at BOTH precisions or the result is uninterpretable, and
  `far` is the predicted target.
- **TTL** (Mikkola & Aarseth 2002). Needs a deliberate mass-ratio sweep through `z8`/`z9`; the
  default charts are equal-mass and would exercise nothing. The equal-mass tie is the control.
- **IAS15** (Rein & Spiegel 2015). A reference arm, not a production candidate -- its per-lane
  variable work is the property already measured as fatal on GPU (`warps hit 1.0000`).
