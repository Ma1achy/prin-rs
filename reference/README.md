# Reference implementation (NumPy)

Validated. **Port `tb_az.py`; do not re-derive the algebra.**

| file | contents |
|---|---|
| `tb.py` | core — leapfrog, energy, pair distances, outcome classification, `burrau_grid` (grid + ensemble construction) |
| `tb_lc.py` | Levi-Civita transform, single-pair regularisation |
| `tb_az.py` | **Aarseth–Zare** — the one to port. Regularised Hamiltonian, RK4 in fictitious time |
| `tb_all_az.py` | AZ plus every per-pixel field in one pass |
| `refine_test.py` | `shape_vec` (Hopf map), dispersion measures |
| `tb_ftle.py` | Benettin FTLE shadows, moment of inertia (`inertia()` is used by `tb_all_az`) |

Pure NumPy, no other dependencies. Verified to import cleanly as a set.

## Smoke test — the number to match

```bash
python3 -c "
import numpy as np, tb, tb_all_az as AA, warnings; warnings.filterwarnings('ignore')
r0,v0,gid,_,_ = tb.burrau_grid(3,3, 1.0,3.0, 0.05, ens=3, jitter_frac=0.5, seed=0)
for mode in ['fixed', 'per-step-interval']:
    for cf in [False, True]:
        res = AA.integrate_all_az(r0,v0, t_max=13.0, n_sync=32, eta=0.01,
                                  dtau_mode=mode, clamp_final=cf)
        print(mode, 'clamp', int(cf), 'median |dE/E| =', np.median(res['drift']))
"
```

Measured here — the four arms, in the order the write-up quotes them:

| arm | `dtau_mode` | `clamp_final` | median `|dE/E|` | max |
|---|---|---|---|---|
| A | `fixed` | off | `3.196673558482950e-09` | `2.9587e-06` |
| B | `per-step-interval` | off | `4.462793760861922e-10` | `3.7858e-08` |
| C | `fixed` | on | `4.046561394010526e-09` | `1.4390e-05` |
| D | `per-step-interval` | on | `4.153297393626583e-10` | `1.4161e-07` |

**A is what every committed NumPy number in the corpus used. D is the default.**

**Three things to know before quoting these.**

A reproduces the committed pre-fix `tb_az.py` **bitwise** — `3.1966735584829495e-09` from
`git show HEAD:reference/tb_az.py` — so the mode switch is faithful and the change in row B is
the fix, not a transcription slip.

**Energy drift is nearly blind to the overshoot, and this table is the evidence.** The clamp
buys **24,000x** on the figure-eight closure error at `eta = 0.02` and raises the convergence
order from 1.13 to 3.06, while moving the median drift here by a factor of 1.3 *the wrong way*
(A → C) and 1.07 the right way (B → D). The overshoot displaces the state in *time*, and the
regularised Hamiltonian's energy is nearly stationary along the flow — so the science field the
`dtau` bug showed up in cannot see this one. Read §1 of `examples/overshoot.rs`, not this table,
for whether the clamp works.

**The number this README used to quote, `3.892633125701676e-09`, does not reproduce on the
unmodified committed reference either.** It is not a casualty of either change; it was already
wrong. A documented reproduction command can be wrong, and only running it finds out.

## `dtau_mode` and `clamp_final`

`integrate_az` and `integrate_all_az` take `dtau_mode`, matching
`src/integrate/az/driver.rs::DtauMode`. `dt = A*B*dtau`, so sizing `dtau` once per sync interval
makes the physical step `eta*dt_left` only while `A*B` stays near its entry value — a trajectory at
a close encounter *at a boundary* has a tiny `A0*B0`, so `dtau` is enormous and `dt` grows by orders
as the bodies separate. `'per-step-interval'` recomputes `A*B` per step with `dt_left` held fixed
and caps at the entry value; `'per-step-remaining'` puts the remaining time in the numerator and is
**Zeno by arithmetic** — `rem_{n+1} = rem_n (1-eta)`, so the interval is never completed. It is kept
as a named measurement axis, not as a candidate.

`clamp_final` (default `True`) lands the **final** step of each interval *on* the boundary
instead of past it. Without it the march exits by overshooting and only the clock is clipped
(`t += min(s.t, dt_left)`) — the state written back is the overshot one, a **first-order** error
at every boundary inside an RK4 march. It is the **partner** of `'per-step-interval'` and not an
independent knob: under `'fixed'` the overshoot is a fixed slice of time and neighbouring
trajectories overshoot alike, so the error is large but spatially *smooth*; under
`'per-step-interval'` the last step's size depends on the local `A*B`, so the overshoot varies
from trajectory to trajectory. Ship the step control without the clamp and a smooth large error
becomes a structured one.

**The nested-arc banding this was first proposed to explain is not caused by it** — all four arms
carry it, including the one predating both changes, and under outcome-class colouring it vanishes.
See `../RESULTS.md` §24.8. The defect is real and independently measured; it is not the cause of
that appearance.

The landing tolerance `LAND_EPS_REL` is **relative to `dt_left`**. An absolute slack is a
different tolerance at every scale under the project's `alpha^{3/2}` time rescaling, and the
bitwise scale-invariance test caught exactly that at `4.24e-15`.

`tb_az.py` uses RK4 — **not** symplectic or time-symmetric. It was built to prove the physics, not to
ship. Match it at f64 first, then change one thing at a time.

---

## The GLSL port and the closure criterion

`glsl_port.py` and `escape_criterion.py` are the reference for **`EscapeRule::Closure`**, attached
rather than derived. They are pure NumPy and self-contained.

- `glsl_port.py` — `decodeIC` and the Hopf map from
  `Ma1achy/principia-ii`, `src/shaders/principia/frag.glsl`. `MU_MAX = 5.0`, `Q_MAX = 2.0`,
  `ALPHA_MIN = 0.05`, mass saturation `MU_MAX*(2*sigmoid(z) - 1)`. Its integrator is a fixed-step
  **leapfrog**, not Aarseth–Zare.
- `escape_criterion.py` — `render()` carries the criterion:
  `fire = (dn < tau) & (E > 0)`, with `dn = |n(t) - n(t - win)|` a chord between the two **ends** of
  the window, `tau = 1e-3`, `win = 0.4`.

**Do not compare trajectories against these.** The integrators differ and diverge legitimately.
What was ported is the *criterion*, and that is what
`tests/outcome_encoding.rs::the_criterion_transcribes_the_reference_including_its_body_ordering`
checks — the same `(m, r, v)` and the same closure value into both, over a 40-row golden table
generated from `_rel` plus the criterion line. 20 rows fire and **14 of those return a body that is
not the one outside the tightest pair**, which is what makes it a test of the ordering rather than
of the arithmetic.
