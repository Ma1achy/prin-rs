# The Levi-Civita branch cut

A conditioning defect in the inverse LC map. It is present in the original
`reference/tb_lc.py`, it affects the Burrau default configuration **at `t = 0`**, and it is
the leading candidate to explain an unresolved f32 dispute in the prior numpy work.

## The mechanism

The LC map sends `u -> rho = u^2` (complex square). Its inverse must undo that, and the
original form always computes the first component and derives the second:

```python
u0 = sqrt((|rho| + rho.x) / 2)
u1 = rho.y / (2 * u0)                 # guarded only by u0 > 1e-300
```

When `rho` points along **negative x**, `|rho| + rho.x` is a difference of near-equal
numbers and cancels catastrophically. The division by the damaged `u0` then amplifies it.
The `1e-300` guard prevents division by zero and does nothing whatever for conditioning.

The fix is standard: compute whichever component is **larger** directly and derive the other.
`u0` is larger when `rho.x >= 0`, `u1` when `rho.x < 0`. Both branches satisfy
`rho_of_u(u) == rho`; the choice is otherwise free because `rho(u) = rho(-u)`.

## The Burrau default sits exactly on the cut

Not "near" it. Exactly on it.

Bodies 1 and 2 start at `(-2, -1)` and `(1, -1)` — the **same `y`** — so their separation is
`(3, 0)`, purely along the x-axis. Whenever the reference body is 2, the regularised pair
`(2, 1)` registers `rho = (-3, 0)`, at exactly 180 degrees, before anything has moved.

Measured over the cross-check geometry (near-field 3x3, nominal copies):

```
  t_max  n_sync  registrations   min u0/|u|    differ   max |du|/|u|
    0.5       1             18      0.000e0      1/18      1.570e-16
    1.0       2             36      0.000e0     10/36      4.104e-12
    2.0       5             90      0.000e0     37/90      1.593e-11
   13.0      32            576      0.000e0    255/576     1.659e-11
```

By `t = 13`, **255 of 576 registrations differ** between the two branches, by up to 1.7e-11
relative — four orders above ulp. Registration happens at every sync boundary, so the loss is
injected `n_sync` times per trajectory rather than once.

## Conditioning

Relative round-trip error, `rho -> u -> rho`, as a function of the angle of `rho`:

```
 angle (deg)       f64 ref    f64 stable       f32 ref    f32 stable
           0       0.000e0       0.000e0       0.000e0       0.000e0
          45     1.110e-16     1.110e-16       0.000e0       0.000e0
          90     1.608e-16     1.608e-16      5.960e-8      5.960e-8
         135     2.220e-16       0.000e0      1.788e-7       0.000e0
         170     3.886e-15     1.110e-16      2.325e-6      5.960e-8
         179     1.912e-13       0.000e0      9.835e-5      5.960e-8
       179.9     3.975e-12       0.000e0      2.213e-2      5.960e-8
      179.99      3.471e-9       0.000e0      1.745e-4       0.000e0
         180     1.225e-16       0.000e0     1.225e-16       0.000e0
```

Worst case over 3600 orientations: **6.206e-11 unstable, 4.108e-16 stable** at f64.

At f32 the unstable form reaches **2.2e-2 — a 2% error in a coordinate**. The stable form
holds 5.96e-8, which is f32 epsilon.

## The consequence is correctness, not precision

The cut is fixed along negative x in the coordinate frame. So the accuracy of a measurement
depends on the **absolute orientation** of the configuration.

**The physics is rotationally invariant. The unstable implementation of it is not.**

Rotating a configuration is a symmetry of the problem and must not change the numerics. That
makes this a defect rather than a tolerance question, independent of how large the loss
happens to be in any particular run.

## It cannot be inert on a cross-check

A natural expectation is that reconditioning a value already accurate at f64 should leave a
cross-check unchanged. It does not, and cannot: **a more accurate registration is a different
trajectory**, so it necessarily stops agreeing bit-for-bit with a reference that registers
less accurately.

The instrument is therefore two tables on identical code, differing only in the branch — not
an unchanged table. Running the Rust with the stable branch against the *unstable* reference:

```
      case   t_max   pinned-unstable   rust-stable-vs-unstable-ref
   az_t0p5     0.5         0.000e+00                     9.992e-16
     az_t1     1.0         0.000e+00                     1.611e-11
    az_t13    13.0         1.930e-10                     4.303e-09
```

With **both sides** conditioned, the comparison is re-established and is far better than it
ever was:

```
      case   t_max   both unstable    both stable
   az_t0p5     0.5      0.000e+00      0.000e+00
     az_t1     1.0      0.000e+00      4.441e-16
     az_t2     2.0      8.882e-16      3.553e-15
     az_t4     4.0      6.847e-12      2.576e-14
     az_t8     8.0      1.463e-10      1.465e-13
    az_t13    13.0      1.930e-10      2.718e-13
```

**Three orders of magnitude at `t = 13`.** This retires a negative result reported in PR #2:
the cross-check missing BRIEF §9's `~1e-10` target was attributed there to chaotic
amplification of ulp-level differences, and an amendment to §9 was proposed on that basis.
That was wrong. The shortfall was branch-cut error injected 32 times per trajectory, and §9's
figure stands as written.

## Effect on integration quality

At f64 the gain is real at short horizons and masked by RK4 truncation later:

```
  t_max |     med ref  med stable    gain |     max ref  max stable    gain
      1 |   4.400e-12   4.787e-13    9.2x |   1.163e-11   5.070e-13   22.9x
      2 |   1.764e-11   1.668e-11    1.1x |   5.100e-11   1.786e-11    2.9x
      4 |   2.688e-10   2.735e-10    1.0x |    9.736e-8    9.737e-8    1.0x
     13 |    3.400e-9    3.404e-9    1.0x |    1.241e-5    1.241e-5    1.0x
```

At f32 it is not masked:

```
  t_max | f32 reference    f32 stable     gain
    0.5 |      2.223e-7      1.497e-7     1.5x
      1 |      3.151e-3      3.723e-7  8464.7x
      2 |      4.692e-3      1.043e-6  4497.9x
      4 |      4.751e-3      2.683e-6  1771.2x
     13 |      5.873e-3      3.705e-6  1584.9x
```

## Bearing on the f32 dispute

The prior numpy work left an unresolved dispute: raw energy drift at f32 looked acceptable —
f32 AZ reportedly beat softened leapfrog at some horizons — while an ensemble diagnostic
broke down early. The working hypothesis was reference-body switching across ensemble copies
(BRIEF §3, §8 experiment 2).

**This defect is a better candidate, and it is not about arithmetic at all.**

- It is orientation-dependent, so copies of one pixel — which differ in configuration — can
  straddle the cut differently. A spread *across copies* would then partly measure
  registration error rather than dynamics. That is exactly the shape of "drift looks fine but
  the ensemble diagnostic breaks".
- It is ~5e8 times worse at f32 than f64, so it appears at f32 and hides at f64.
- Single-trajectory energy drift is comparatively insensitive to it (1585x at t=13, but from
  5.9e-3 to 3.7e-6 — both of which a drift-only test might wave through depending on the
  threshold), while a cross-copy statistic is not.

Any f32 conclusion drawn from the numpy harness before this patch was measuring AZ **plus a
known conditioning defect**, and a verdict of "AZ is unusable at f32" would have been right
about the symptom and wrong about the cause.

## Where it lives

- `src/integrate/az/lc.rs` — `u_of_rho` (stable, kernel default) and `u_of_rho_reference`
- `reference/tb_lc.py` — `u_of_rho_stable`, `u_of_rho_unstable`, and the `USE_STABLE_LC`
  module switch; the unstable form is retained because prior results were produced with it
- `tests/lc_conditioning.rs` — the tables above, asserted
- `examples/lc_cut_proximity.rs` — the registration-divergence measurement
- `examples/lc_branch_effect.rs` — the integration-quality comparison, f64 and f32
- `tools/xcheck/horizon.py [--lc-unstable]` — both cross-check paths
