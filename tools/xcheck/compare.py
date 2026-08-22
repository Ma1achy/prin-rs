#!/usr/bin/env python3
"""Compare two cross-check dumps.

Reports **per column** — max absolute, max relative, and the row index where the max
relative deviation occurs. An aggregate max hides which term is wrong: a transcription
error in a single g2 term shows up as position diverging while d_min does not.

Exits nonzero if any column exceeds the tolerance.
"""
import argparse
import sys


def load(path):
    hdr, cols, rows = [], None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                hdr.append(line)
                if line.startswith("# columns="):
                    cols = line.split("=", 1)[1].split(",")
                continue
            rows.append([float(x) for x in line.split("\t")])
    if cols is None:
        raise SystemExit(f"{path}: no '# columns=' header line")
    return hdr, cols, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("rs")
    ap.add_argument("--tol", type=float, required=True)
    a = ap.parse_args()

    h1, c1, r1 = load(a.ref)
    h2, c2, r2 = load(a.rs)

    if h1 != h2:
        print("HEADER MISMATCH — the two sides are not running the same case:")
        for x, y in zip(h1, h2):
            if x != y:
                print(f"  ref: {x}\n  rs : {y}")
        return 2
    if len(r1) != len(r2):
        print(f"ROW COUNT MISMATCH: ref {len(r1)}, rs {len(r2)}")
        return 2

    ncol = len(c1)
    worst = 0.0
    print(f"{'column':<14}{'max abs':>14}{'max rel':>14}{'argmax row':>12}")
    print("-" * 54)
    fail = []
    for k in range(1, ncol):  # column 0 is the index
        mabs = mrel = 0.0
        arg = -1
        for i in range(len(r1)):
            x, y = r1[i][k], r2[i][k]
            if x != x and y != y:  # both NaN: agreement
                continue
            d = abs(x - y)
            scale = max(abs(x), abs(y), 1e-300)
            rel = d / scale
            if d > mabs:
                mabs = d
            if rel > mrel:
                mrel, arg = rel, i
        worst = max(worst, mrel)
        flag = "" if mrel <= a.tol else "   <-- FAIL"
        if mrel > a.tol:
            fail.append(c1[k])
        print(f"{c1[k]:<14}{mabs:>14.3e}{mrel:>14.3e}{arg:>12}{flag}")

    print("-" * 54)
    print(f"worst relative deviation: {worst:.3e}   tolerance: {a.tol:.1e}")
    if fail:
        print(f"FAIL: {', '.join(fail)}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
