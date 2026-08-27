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
res = AA.integrate_all_az(r0,v0, t_max=13.0, n_sync=32, eta=0.01)
print('median |dE/E| =', np.median(res['drift']))   # expect ~3.9e-09
"
```

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
