#!/usr/bin/env python3
"""The divergence-vs-horizon table.

A single pass/fail at t=13 cannot distinguish a transcription error from ulp-level
divergence amplified by chaos. The growth curve can: a correct port shows divergence rising
like exp(lambda t) from an O(1e-16) intercept, with lambda near the measured Lyapunov
exponent ~0.7. A wrong port shows the wrong intercept, or the wrong slope.
"""
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
import cases  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
OUT = os.path.join(ROOT, "xcheck_out")
CARGO = os.path.expanduser("~/.cargo/bin/cargo")

STATE_COLS = ["r0x", "r0y", "r1x", "r1y", "r2x", "r2y",
              "v0x", "v0y", "v1x", "v1y", "v2x", "v2y"]


def load(path):
    cols, rows = None, []
    for line in open(path):
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith("#"):
            if line.startswith("# columns="):
                cols = line.split("=", 1)[1].split(",")
            continue
        rows.append([float(x) for x in line.split("\t")])
    return cols, rows


def main():
    os.makedirs(OUT, exist_ok=True)
    which = "UNSTABLE (original)" if "--lc-unstable" in sys.argv else "STABLE"
    print(f"inverse LC branch, BOTH sides: {which}")
    print(f"{'case':>10}{'t_max':>8}{'n_sync':>8}{'max |dr|':>13}{'max rel':>13}{'refs':>7}{'drift(rs)':>12}")
    print("-" * 71)
    prev = None
    for name in cases.AZ_HORIZONS:
        c = cases.CASES[name]
        ref_p = os.path.join(OUT, f"ref_{name}.tsv")
        rs_p = os.path.join(OUT, f"rs_{name}.tsv")
        pycmd = [sys.executable, os.path.join(os.path.dirname(__file__), "dump_ref.py"),
                 "--case", name, "--out", ref_p]
        if "--lc-unstable" in sys.argv:
            pycmd.append("--lc-unstable")
        subprocess.run(pycmd, check=True, stdout=subprocess.DEVNULL, cwd=ROOT)
        cmd = [CARGO, "run", "--release", "--quiet", "--bin", "xcheck", "--",
               "--case", name, "--out", rs_p]
        if "--lc-unstable" not in sys.argv:
            cmd.append("--lc-stable")
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, cwd=ROOT)

        cols, a = load(ref_p)
        _, b = load(rs_p)
        idx = [cols.index(x) for x in STATE_COLS]
        ridx = [i for i, cn in enumerate(cols) if cn.startswith("ref")]
        dcol = cols.index("drift")

        mabs = mrel = 0.0
        for ra, rb in zip(a, b):
            for k in idx:
                d = abs(ra[k] - rb[k])
                mabs = max(mabs, d)
                mrel = max(mrel, d / max(abs(ra[k]), abs(rb[k]), 1e-300))
        refs_equal = all(ra[k] == rb[k] for ra, rb in zip(a, b) for k in ridx)
        drift = max(rb[dcol] for rb in b)

        growth = ""
        if prev is not None and prev[1] > 0:
            dt = c["t_max"] - prev[0]
            if mabs > 0 and dt > 0:
                lam = math.log(mabs / prev[1]) / dt
                growth = f"   growth exp({lam:.2f} t)"
        prev = (c["t_max"], mabs)

        print(f"{name:>10}{c['t_max']:>8}{c['n_sync']:>8}{mabs:>13.3e}{mrel:>13.3e}"
              f"{'ok' if refs_equal else 'DIFFER':>7}{drift:>12.2e}{growth}")

    print("-" * 71)
    print("Measured Lyapunov exponent for this configuration is ~0.7 (BRIEF §2.1), so a")
    print("correct port should show growth near exp(0.7 t) from an O(1e-16) intercept.")


if __name__ == "__main__":
    main()
