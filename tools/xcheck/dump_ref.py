#!/usr/bin/env python3
"""Python side of the cross-check: emits a TSV the Rust side must reproduce.

Run from the repo root:
    python3 tools/xcheck/dump_ref.py --case algebra --out xcheck_out/ref_algebra.tsv
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "reference"))

import numpy as np  # noqa: E402
import tb  # noqa: E402
import tb_az  # noqa: E402
import tb_lc  # noqa: E402
import tb_ftle  # noqa: E402
import refine_test  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import cases  # noqa: E402


def dump_algebra(name, out):
    c = cases.CASES[name]
    cfgs = cases.random_configs(c["n"], c["seed"])
    r = np.array([g[0] for g in cfgs], dtype=np.float64)
    v = np.array([g[1] for g in cfgs], dtype=np.float64)

    energy = tb.energy(r, v, 0.0)
    pd = tb.pair_dists(r)
    inertia = tb_ftle.inertia(r)
    hyper = np.sqrt(inertia / tb.M.sum())
    n = refine_test.shape_vec(r)

    with open(out, "w") as f:
        for line in cases.header_lines(name):
            f.write(line + "\n")
        for i in range(c["n"]):
            row = [
                energy[i], pd[i, 0], pd[i, 1], pd[i, 2],
                inertia[i], hyper[i], n[i, 0], n[i, 1], n[i, 2],
            ]
            f.write(str(i) + "\t" + "\t".join("%.17e" % x for x in row) + "\n")
    print(f"wrote {out}: {c['n']} rows")


def dump_az(name, out):
    """Nominal copies only — `ens=0`, so no RNG participates on either side.

    `burrau_grid` draws jitter from ONE global PCG64 stream over the whole array, while
    BRIEF §7 requires per-pixel seeding from (i,j,seed). Those are incompatible, so jittered
    copies could never be the basis of a 1e-10 comparison between two implementations. The
    reference already contains the way out: `mask[::reps] = False` leaves copy 0 of every
    cell un-jittered and seed-independent, which is exactly what is compared here.

    The sync loop is stepped one boundary at a time so the chosen reference body can be
    recorded per boundary. That is equivalent to one monolithic call — asserted below, not
    assumed — because with n_sync=1 the sub-interval is the whole span and `t` lands exactly
    on each boundary.
    """
    c = cases.CASES[name]
    r0, v0, gid, _, hx = tb.burrau_grid(
        c["nx"], c["ny"], c["cx"], c["cy"], c["half"], body=c["body"], ens=0
    )
    n_sync, t_max = c["n_sync"], c["t_max"]
    kw = dict(eta=c["eta"], max_steps=c["max_steps"])

    e0 = tb.energy(r0, v0, 0.0)
    r, v = r0.copy(), v0.copy()
    refs = np.zeros((len(r0), n_sync), dtype=np.int64)
    dmin = np.full(len(r0), np.inf)
    switches = np.zeros(len(r0), dtype=np.int64)
    prev = np.full(len(r0), -1)
    t = np.zeros(len(r0))

    # dt_left must be formed exactly as the monolithic loop forms it:
    #     t_target = (kk+1)*t_max/n_sync,  dt_left = t_target - t
    # Using a precomputed `step = t_max/n_sync` instead differs at the last ulp, and that
    # difference amplifies to ~3e-13 by t=2 through the same chaotic growth the cross-check
    # is trying to measure. Measured, not assumed: it is what the assertion below caught.
    t_prev = 0.0
    for kk in range(n_sync):
        ref = tb_az.choose_reference(r)
        refs[:, kk] = ref
        switches += (prev >= 0) & (ref != prev)
        prev = ref.copy()
        t_target = (kk + 1) * t_max / n_sync
        out_k = tb_az.integrate_az(r, v, t_target - t_prev, n_sync=1, **kw)
        r, v = out_k["r"], out_k["v"]
        dmin = np.minimum(dmin, out_k["dmin"])
        t = t + out_k["t"]
        t_prev = t_target

    # The stepped loop must agree with one monolithic call, or the per-sync reference log is
    # describing a different trajectory from the one being compared.
    mono = tb_az.integrate_az(r0, v0, t_max, n_sync=n_sync, **kw)
    dev = np.max(np.abs(r - mono["r"])) + np.max(np.abs(v - mono["v"]))
    assert dev == 0.0, (
        f"stepped sync loop is not bit-identical to the monolithic call (dev={dev:e}); "
        "the per-sync reference log would then describe a different trajectory"
    )

    drift = np.abs((tb.energy(r, v, 0.0) - e0) / np.maximum(np.abs(e0), 1e-30))

    with open(out, "w") as f:
        for line in cases.header_lines(name):
            f.write(line + "\n")
        for i in range(len(r0)):
            row = list(r[i].ravel()) + list(v[i].ravel())
            row += [t[i], dmin[i], drift[i], float(switches[i])]
            row += [float(x) for x in refs[i]]
            f.write(str(i) + "\t" + "\t".join("%.17e" % x for x in row) + "\n")
    print(f"wrote {out}: {len(r0)} rows, n_sync={n_sync}, stepped-vs-monolithic dev={dev:.2e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True, choices=sorted(cases.CASES))
    ap.add_argument("--out", required=True)
    ap.add_argument("--lc-unstable", action="store_true",
                    help="use the reference's original (unconditioned) inverse LC branch")
    a = ap.parse_args()
    tb_lc.USE_STABLE_LC = not a.lc_unstable
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    kind = cases.CASES[a.case]["kind"]
    if kind == "algebra":
        dump_algebra(a.case, a.out)
    elif kind == "az":
        dump_az(a.case, a.out)
    else:
        raise SystemExit(f"unhandled case kind: {kind}")


if __name__ == "__main__":
    main()
