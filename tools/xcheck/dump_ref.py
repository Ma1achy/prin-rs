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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True, choices=sorted(cases.CASES))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    kind = cases.CASES[a.case]["kind"]
    if kind == "algebra":
        dump_algebra(a.case, a.out)
    else:
        raise SystemExit(f"unhandled case kind: {kind}")


if __name__ == "__main__":
    main()
