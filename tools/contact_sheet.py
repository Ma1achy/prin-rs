"""Tile `results/logh_arms` panels into one comparable sheet per (case, field).

**The panels are already directly comparable and this does not rescale them.** The drift ramp is
a fixed `1e-12..1e2` in the harness, the gain ramp a fixed symmetric +/-4 decades, and the shape
window is taken once from the `az` arm and shared -- so a sheet is a LAYOUT, not a rendering
decision. An auto-ranged montage would undo all three and manufacture or hide the thing being
compared, which is on this project's record from the bleaching strip.

Nearest-neighbour paste, no resampling: if a sheet looks soft, a viewer is upscaling it.

    python3 tools/contact_sheet.py far deep_interior near-field
    python3 tools/contact_sheet.py --root results_regen far deep_interior near-field

**The output root is an ARGUMENT, not a constant.** It was hardcoded to `results/`, which is the
third site of that defect on this project -- `criterion_metric` wrote to a hardcoded `results/`
so a reduced-`levels` validation pass would overwrite committed artefacts, and `pan_sequence`
took a `viewport` argument while `frame_res` stayed a literal. Here it meant that running this
against a regenerated tree would write the sheets straight into the committed one, silently
pairing new panels with an old montage.
"""
import sys, os
from PIL import Image, ImageDraw, ImageFont

argv = sys.argv[1:]
ROOT = "results"
if argv and argv[0] == "--root":
    ROOT = argv[1]
    argv = argv[2:]

D = f"{ROOT}/logh_arms"
OUT = f"{ROOT}/logh_arms/sheets"
if not os.path.isdir(D):
    sys.exit(f"no such panel directory: {D}  (pass --root <dir> before the case names)")
os.makedirs(OUT, exist_ok=True)

ARMS = ["az", "heggie", "logh_rk4", "logh_lf", "logh_rk4_nolim", "plain_rk4"]
PAD, LAB, TITLE = 8, 22, 34

def font(sz):
    for p in ("/System/Library/Fonts/Supplemental/Arial.ttf",
              "/System/Library/Fonts/Helvetica.ttc"):
        if os.path.exists(p):
            try: return ImageFont.truetype(p, sz)
            except Exception: pass
    return ImageFont.load_default()

def sheet(case, field, arms, title, note=""):
    have = [(a, f"{D}/{case}_{a}_{field}.png") for a in arms]
    have = [(a, p) for a, p in have if os.path.exists(p)]
    if not have: return None
    ims = [(a, Image.open(p).convert("RGB")) for a, p in have]
    w, h = ims[0][1].size
    n = len(ims)
    W = n * w + (n + 1) * PAD
    H = TITLE + LAB + h + 2 * PAD + (LAB if note else 0)
    sh = Image.new("RGB", (W, H), (250, 250, 250))
    d = ImageDraw.Draw(sh)
    d.text((PAD, 7), title, fill=(15, 15, 15), font=font(19))
    for i, (a, im) in enumerate(ims):
        x = PAD + i * (w + PAD)
        d.text((x, TITLE + 3), a, fill=(30, 30, 30), font=font(15))
        sh.paste(im, (x, TITLE + LAB))
    if note:
        d.text((PAD, TITLE + LAB + h + 4), note, fill=(70, 70, 70), font=font(13))
    out = f"{OUT}/{case}_{field}.png"
    sh.save(out)
    return out

for case in argv:
    print(sheet(case, "drift",
        ARMS,
        f"{case} — energy drift, inferno ramp FIXED at 1e-12..1e2 across every case and arm",
        "magenta = undetermined (non-finite / vetoed). Brighter = worse. Identical window everywhere, so panels are comparable across arms AND across cases."))
    print(sheet(case, "gain_vs_heggie",
        ["logh_rk4", "logh_lf", "logh_rk4_nolim", "plain_rk4", "plain_lf"],
        f"{case} — log10(drift_heggie / drift_arm), fixed symmetric +/-4 decades",
        "BLUE = this arm beats Heggie.  RED = Heggie beats it.  near-white = agree.  magenta = either undetermined."))
    print(sheet(case, "uniform", ARMS,
        f"{case} — spread_shape, the shipping science field (window from the az arm, SHARED)",
        "Science pass, termination ON. Not the field to read a numerical defect from — that is the drift sheet."))
