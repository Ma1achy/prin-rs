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
    res = AA.integrate_all_az(r0,v0, t_max=13.0, n_sync=32, eta=0.01, dtau_mode=mode)
    print(mode, 'median |dE/E| =', np.median(res['drift']))
"
```

Measured here:

| `dtau_mode` | median `|dE/E|` | max |
|---|---|---|
| `fixed` (the behaviour every committed number used) | `3.196673558482950e-09` | `2.9587e-06` |
| `per-step-interval` (the default) | `4.462793760861922e-10` | `3.7858e-08` |

**Two things to know before quoting these.**

`fixed` reproduces the committed pre-fix `tb_az.py` **bitwise** — `3.1966735584829495e-09` from
`git show HEAD:reference/tb_az.py` — so the mode switch is faithful and the change in the second
row is the fix, not a transcription slip.

**The number this README used to quote, `3.892633125701676e-09`, does not reproduce on the
unmodified committed reference either.** It is not a casualty of the `dtau` change; it was already
wrong. A documented reproduction command can be wrong, and only running it finds out.

## `dtau_mode`

`integrate_az` and `integrate_all_az` take `dtau_mode`, matching
`src/integrate/az/driver.rs::DtauMode`. `dt = A*B*dtau`, so sizing `dtau` once per sync interval
makes the physical step `eta*dt_left` only while `A*B` stays near its entry value — a trajectory at
a close encounter *at a boundary* has a tiny `A0*B0`, so `dtau` is enormous and `dt` grows by orders
as the bodies separate. `'per-step-interval'` recomputes `A*B` per step with `dt_left` held fixed
and caps at the entry value; `'per-step-remaining'` puts the remaining time in the numerator and is
**Zeno by arithmetic** — `rem_{n+1} = rem_n (1-eta)`, so the interval is never completed. It is kept
as a named measurement axis, not as a candidate.

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
