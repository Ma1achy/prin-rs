"""Lay the per-commit panels out as a strip, in commit order, with the stop-the-eye numbers
under each. A strip is the point: the question is WHICH COMMIT, and that is a comparison
between neighbours, not a property of any one panel."""
import sys, os
from PIL import Image, ImageDraw, ImageFont

OUT, RES, TILE = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
KINDS = ["uniform", "outcome", "drift"]
COMMITS = sys.argv[4:]

stats = {}
for line in open(f"{OUT}/stats.tsv"):
    f = line.rstrip("\n").split("\t")
    if f[0] == "STAT":
        stats[f[1]] = f

try:
    F = ImageFont.truetype("/System/Library/Fonts/Menlo.ttc", 15)
    FB = ImageFont.truetype("/System/Library/Fonts/Menlo.ttc", 19)
except OSError:
    F = FB = ImageFont.load_default()

SUB = {c: t for c, t in (l.split(" ", 1) for l in open(f"{OUT}/subjects.txt").read().splitlines())}

for kind in KINDS:
    have = [c for c in COMMITS if os.path.exists(f"{OUT}/{c}_{kind}.png")]
    if not have:
        continue
    PAD, HDR, FTR = 6, 30, 92
    W = len(have) * (TILE + PAD) + PAD
    H = HDR + TILE + FTR
    im = Image.new("RGB", (W, H), (18, 18, 22))
    d = ImageDraw.Draw(im)
    for i, c in enumerate(have):
        x = PAD + i * (TILE + PAD)
        t = Image.open(f"{OUT}/{c}_{kind}.png").convert("RGB").resize((TILE, TILE), Image.NEAREST)
        im.paste(t, (x, HDR))
        d.text((x, 6), f"{i+1}. {c}", font=FB, fill=(235, 235, 240))
        s = stats.get(c)
        y = HDR + TILE + 6
        if s:
            rows = [
                f"nonfin {int(s[4]):>7}  hot {float(s[7]):.4f}",
                f"esc {float(s[8]):.4f} col {float(s[9]):.4f} bnd {float(s[10]):.4f}",
                f"drift ramp {float(s[11]):.2e}..{float(s[12]):.2e}",
                f"drift p50 {float(s[13]):.2e}  {float(s[3]):.0f}s",
                f"own spread p1/p99 {float(s[14]):.2e}..{float(s[15]):.2e}",
            ]
            for r in rows:
                d.text((x, y), r, font=F, fill=(190, 190, 200))
                y += 17
    im.save(f"{OUT}/strip_{kind}.png")
    print(f"{OUT}/strip_{kind}.png  {im.size}")
