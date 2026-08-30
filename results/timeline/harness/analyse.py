"""Read the per-commit dumps and answer the two questions the strip cannot answer by eye:
WHERE the outcome flips are, and whether they sit in the former magenta halos."""
import sys, os
import numpy as np

OUT = sys.argv[1]
RES = int(sys.argv[2])
COMMITS = sys.argv[3:]

DT = np.dtype([("outcome", "u1"), ("state", "u1"), ("nonfin", "u1"), ("pad", "u1"),
               ("drift", "<f4"), ("n0", "<f4"), ("n1", "<f4"), ("n2", "<f4"), ("tend", "<f4")])

D = {}
for c in COMMITS:
    p = f"{OUT}/{c}.bin"
    if os.path.exists(p):
        D[c] = np.fromfile(p, dtype=DT).reshape(RES, RES)

def dilate(mask, r):
    if r == 0:
        return mask
    out = mask.copy()
    for _ in range(r):
        n = out.copy()
        n[1:, :] |= out[:-1, :]; n[:-1, :] |= out[1:, :]
        n[:, 1:] |= out[:, :-1]; n[:, :-1] |= out[:, 1:]
        out = n
    return out

STATE = {0: "escape", 1: "bounded", 2: "collision", 3: "running", 6: "simfail", 7: "decodefail"}

print("== PER-COMMIT STATE CENSUS ==")
print(f"{'commit':>9} {'nonfin':>8} {'hot':>8} {'escape':>8} {'bounded':>8} {'collis':>8} "
      f"{'other':>8} {'drift p50':>11} {'drift p99':>11} {'|n0|>0.9':>9}")
for c in COMMITS:
    if c not in D: continue
    a = D[c].ravel(); n = a.size
    st = a["state"]; dr = a["drift"].astype(np.float64)
    fin = np.isfinite(dr)
    pale = np.abs(a["n0"]) > 0.9
    print(f"{c:>9} {int((a['nonfin']>0).sum()):>8} {(~(dr<=1e-6)).mean():>8.4f} "
          f"{(st==0).mean():>8.4f} {(st==1).mean():>8.4f} {(st==2).mean():>8.4f} "
          f"{(~np.isin(st,[0,1,2])).mean():>8.4f} "
          f"{np.percentile(dr[fin],50):>11.3e} {np.percentile(dr[fin],99):>11.3e} {pale.mean():>9.4f}")

print()
print("== ADJACENT-PAIR MOVEMENT ==")
print(f"{'pair':>21} {'flips':>9} {'frac':>8} {'moved':>9} {'frac':>8} {'chord p50':>11} {'chord max':>11}")
pairs = list(zip(COMMITS, COMMITS[1:]))
for a, b in pairs:
    if a not in D or b not in D: continue
    A, B = D[a].ravel(), D[b].ravel()
    fl = A["outcome"] != B["outcome"]
    ch = np.sqrt((A["n0"].astype(np.float64)-B["n0"])**2 + (A["n1"].astype(np.float64)-B["n1"])**2
                 + (A["n2"].astype(np.float64)-B["n2"])**2)
    mv = np.isfinite(ch) & (ch != 0)
    f = ch[np.isfinite(ch)]
    print(f"{a[:7]+'->'+b[:7]:>21} {int(fl.sum()):>9} {fl.mean():>8.4f} {int(mv.sum()):>9} "
          f"{mv.mean():>8.4f} {np.percentile(f,50):>11.3e} {f.max():>11.3e}")


print()
print("== THE BLEACHING, MEASURED IN OKLab ON THE RENDERED PANEL ==")
print("The production map is `oklab_to_srgb([L(scalar), a, b])` with `(a,b)` a vMF-weighted")
print("blend of LANDMARK colours -- so chroma collapses wherever the shape vector sits between")
print("sites, and high L with near-zero chroma is exactly the white the eye reads. Inverted")
print("here from the PNG itself, on the FIXED ramp shared by the strip. L in [0.30, 0.92],")
print("C_MAX 0.13. `white` is L>0.80 and chroma<0.030, magenta excluded.")
try:
    from PIL import Image

    def oklab(rgb):
        c = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
        r, g, bl = c[..., 0], c[..., 1], c[..., 2]
        l = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * bl
        m = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * bl
        s2 = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * bl
        l, m, s2 = np.cbrt(l), np.cbrt(m), np.cbrt(s2)
        L = 0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s2
        A = 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s2
        B = 0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s2
        return L, A, B

    print(f"{'commit':>9} {'white':>8} {'L p50':>8} {'L>0.85':>8} {'chroma p50':>11} "
          f"{'C<0.03':>8} {'magenta':>9} {'|n0|>0.99':>10}")
    for c in COMMITS:
        f = f"{OUT}/{c}_uniform.png"
        if not os.path.exists(f) or c not in D:
            continue
        im = np.asarray(Image.open(f).convert("RGB")).astype(np.float64) / 255.0
        mag = (im[..., 0] > 0.99) & (im[..., 1] < 0.01) & (im[..., 2] > 0.99)
        L, A, B = oklab(im)
        ch = np.hypot(A, B)
        ok = ~mag
        white = (L > 0.80) & (ch < 0.030) & ok
        aa = D[c].ravel()
        print(f"{c:>9} {white.mean():>8.4f} {np.median(L[ok]):>8.4f} {(L[ok] > 0.85).mean():>8.4f} "
              f"{np.median(ch[ok]):>11.4f} {(ch[ok] < 0.030).mean():>8.4f} {mag.mean():>9.4f} "
              f"{(np.abs(aa['n0']) > 0.99).mean():>10.4f}")
except ImportError:
    print("  PIL missing")

print()
print("== WHERE ARE THE FLIPS? enrichment near the BEFORE-arm's non-finite pixels ==")
print("`base` is the frame fraction inside the dilated magenta set; `in` is the flipped")
print("fraction inside it. enrichment = in/base. 1.0 means the flips ignore the halos.")
for a, b in pairs:
    if a not in D or b not in D: continue
    A, B = D[a], D[b]
    fl = (A["outcome"] != B["outcome"])
    mag = A["nonfin"] > 0
    if mag.sum() == 0 or fl.sum() == 0:
        print(f"{a[:7]+'->'+b[:7]:>21}  magenta {int(mag.sum()):>7}  flips {int(fl.sum()):>8}  -- skipped")
        continue
    row = f"{a[:7]+'->'+b[:7]:>21}  magenta {int(mag.sum()):>7}  flips {int(fl.sum()):>8} ({fl.mean():.4f}) "
    for r in (0, 1, 2, 4, 8, 16):
        d = dilate(mag, r)
        base = d.mean()
        inside = fl[d].sum() / fl.sum()
        row += f"  r{r}:{inside/base:5.2f}"
    print(row)

print()
print("== THE BLEACHING, MEASURED AS TEXTURE ==")
print("The eye reads the pale areas going FLAT, not going brighter -- median lightness in fact")
print("FALLS across the walk. So the statistic is local contrast: the 5x5 standard deviation of")
print("OKLab L and of chroma, read off `_uniform.png` on the fixed ramp. `coh` is the lag-1")
print("neighbour correlation of chroma: AMPLITUDE CANNOT TELL A SMALL REAL SIGNAL FROM NOISE,")
print("COHERENCE CAN, so a texture that vanishes while coherence RISES was the noise.")
try:
    from PIL import Image
    from scipy import ndimage

    def oklab2(rgb):
        c = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
        r, g, bl = c[..., 0], c[..., 1], c[..., 2]
        l = np.cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * bl)
        m = np.cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * bl)
        s2 = np.cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * bl)
        return (0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s2,
                1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s2,
                0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s2)

    K = np.ones((5, 5)) / 25.0

    def locsd(X):
        return np.sqrt(np.maximum(
            ndimage.convolve(X * X, K, mode="reflect")
            - ndimage.convolve(X, K, mode="reflect") ** 2, 0))

    def coh(X, m):
        o = []
        for ax in (0, 1):
            mm = m & np.roll(m, 1, axis=ax)
            o.append(np.corrcoef(X[mm], np.roll(X, 1, axis=ax)[mm])[0, 1])
        return o

    print(f"{'commit':>9} {'magenta':>8} {'L p50':>7} {'Lvar p50':>9} {'Cvar p50':>9} "
          f"{'Cvar|pale':>10} {'C coh-x':>8} {'C coh-y':>8} {'L coh-x':>8}")
    for c in COMMITS:
        f = f"{OUT}/{c}_uniform.png"
        if not os.path.exists(f):
            continue
        im = np.asarray(Image.open(f).convert("RGB")).astype(np.float64) / 255.0
        mag = (im[..., 0] > 0.99) & (im[..., 1] < 0.01) & (im[..., 2] > 0.99)
        L, A, B = oklab2(im)
        C = np.hypot(A, B)
        ok = ~mag
        pale = (L > 0.78) & ok
        cx, cy = coh(C, ok)
        lx, _ = coh(L, ok)
        print(f"{c:>9} {mag.mean():>8.4f} {np.median(L[ok]):>7.4f} {np.median(locsd(L)[ok]):>9.5f} "
              f"{np.median(locsd(C)[ok]):>9.5f} "
              f"{(np.median(locsd(C)[pale]) if pale.any() else float('nan')):>10.5f} "
              f"{cx:>8.4f} {cy:>8.4f} {lx:>8.4f}")
except ImportError:
    print("  PIL/scipy missing")
