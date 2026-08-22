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
