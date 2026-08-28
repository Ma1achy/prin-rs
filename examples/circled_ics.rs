//! **What do the ICs under the circled regions have in common?**
//!
//! The regions are marked on `results/closure/config_stability_stop0_uniform.png`. They are pale,
//! low-chroma, magenta-speckled *wedges* with straight edges, sitting on the boundaries between
//! the coloured ribbons.
//!
//! # The population is selected two ways, and both are reported
//!
//! 1. **By hand**, from ellipses digitised off the marked-up screenshot. Written as fractions of
//!    the frame so the digitisation is inspectable, and dumped as an overlay PNG so it can be
//!    checked against the photo rather than trusted.
//! 2. **By property** — pale and low-chroma on the rendered panel. If the two masks agree, the
//!    hand digitisation is not doing the work and the finding is about the field rather than
//!    about where a circle was drawn.
//!
//! **Straight edges are the clue that says where to look.** A fractal boundary is a property of
//! the dynamics; a straight edge in the chart plane is a property of the *initial conditions*, and
//! those are computable with no integration at all. `grid::decode_state` gives masses, positions
//! and velocities per pixel in milliseconds, so every IC statistic here is exact and free.
//!
//! # Writes
//!
//! `<out>/circled_mask.png` and stdout.
use rayon::prelude::*;

use prin_rs::grid::{self, Chart};
use prin_rs::physics::{energy, newton, THIRD};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);

/// `(cx, cy, rx, ry)` as fractions of the frame, origin top-left, digitised from the marked-up
/// screenshot. Six marks: one large top-right, three small around the centre, one mid-left,
/// one large lower-right.
const CIRCLES: [(f64, f64, f64, f64); 6] = [
    (0.900, 0.253, 0.105, 0.210),
    (0.362, 0.462, 0.070, 0.058),
    (0.487, 0.510, 0.068, 0.068),
    (0.625, 0.535, 0.072, 0.045),
    (0.435, 0.720, 0.180, 0.092),
    (0.810, 0.845, 0.195, 0.157),
];

/// Minimal RGB8 read via the `png` crate, which is already a dependency of the writer.
fn read_rgb(path: &str) -> (Vec<u8>, usize, usize) {
    let dec = png::Decoder::new(std::fs::File::open(path).expect("open panel"));
    let mut r = dec.read_info().expect("png header");
    let mut buf = vec![0u8; r.output_buffer_size()];
    let info = r.next_frame(&mut buf).expect("png frame");
    let (w, h) = (info.width as usize, info.height as usize);
    let ch = info.color_type.samples();
    let out = if ch == 3 {
        buf[..w * h * 3].to_vec()
    } else {
        (0..w * h).flat_map(|i| [buf[i * ch], buf[i * ch + 1], buf[i * ch + 2]]).collect()
    };
    (out, w, h)
}

fn oklab(c: [f64; 3]) -> (f64, f64, f64) {
    let f = |x: f64| if x <= 0.04045 { x / 12.92 } else { ((x + 0.055) / 1.055).powf(2.4) };
    let (r, g, b) = (f(c[0]), f(c[1]), f(c[2]));
    let l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b).cbrt();
    let m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b).cbrt();
    let s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b).cbrt();
    (
        0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
        1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
        0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s,
    )
}

fn q(v: &mut Vec<f64>, p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[(((v.len() - 1) as f64) * p).round() as usize]
}

/// Everything about one pixel's INITIAL condition. No integration.
struct Ics {
    alpha: f64,
    beta: f64,
    d: [f64; 3],
    longest: usize,
    tightest: usize,
    reference: usize,
    dmin0: f64,
    dmax0: f64,
    aspect: f64,
    e0: f64,
    kinetic: f64,
    virial: f64,
    lz: f64,
    vmax: f64,
    area: f64,
    /// `d[2nd longest] / d[longest]`. **1.0 is the reference-body switching boundary** -- the two
    /// longest sides are equal, so `argmax` is a coin flip and the LC registration on either side
    /// of the line is a different one. A near-tie is exactly what draws a STRAIGHT edge in the
    /// chart plane, where a fractal boundary would not.
    tie_long: f64,
    /// `d[shortest] / d[2nd shortest]`. 1.0 is the boundary where the *tightest pair* identity
    /// flips, which is what `spread_event` and the closure criterion read.
    tie_tight: f64,
}

fn main() {
    let png: String = std::env::args().nth(1).unwrap_or_else(|| {
        "results/closure/config_stability_stop0_uniform.png".into()
    });
    let out: String = std::env::args().nth(2).unwrap_or_else(|| "circled_out".into());
    let _ = std::fs::create_dir_all(&out);

    let (rgb, w, h) = read_rgb(&png);
    println!("panel {png} -- {w}x{h}\n");

    let z0 = prin_rs::physics::decoder::Latent {
        z_alpha: Z[1], z_beta: Z[0],
        z_q: [Z[4], Z[5], Z[6], Z[7]], z_mu: [Z[8], Z[9]],
    };
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[1] = 1.0;
    q2[0] = 1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (cx, cy, half) = (2.0 * PAN.0 - 1.0 + ZOOM, 2.0 * PAN.1 - 1.0 + ZOOM, ZOOM);
    let sl = grid::Slice::body_plane(w, h, cx, cy, half, 0).with_chart(chart);

    // ---- the two masks -------------------------------------------------------------------
    let mut hand = vec![false; w * h];
    for y in 0..h {
        for x in 0..w {
            let (u, v) = ((x as f64 + 0.5) / w as f64, (y as f64 + 0.5) / h as f64);
            hand[y * w + x] = CIRCLES.iter().any(|&(a, b, rx, ry)| {
                ((u - a) / rx).powi(2) + ((v - b) / ry).powi(2) <= 1.0
            });
        }
    }
    let mut pale = vec![false; w * h];
    let mut magenta = vec![false; w * h];
    for y in 0..h {
        for x in 0..w {
            let o = (y * w + x) * 3;
            let p = [rgb[o], rgb[o + 1], rgb[o + 2]];
            let c = [p[0] as f64 / 255.0, p[1] as f64 / 255.0, p[2] as f64 / 255.0];
            magenta[y * w + x] = p[0] > 250 && p[1] < 5 && p[2] > 250;
            let (l, a, b) = oklab(c);
            // The circled wedges are the PALE, LOW-CHROMA population: high lightness with the
            // vMF blend collapsed. Thresholds are stated, and the overlay says whether they
            // select what was circled.
            pale[y * w + x] = !magenta[y * w + x] && l > 0.86 && a.hypot(b) < 0.045;
        }
    }

    // Overlay, so the digitisation is checkable rather than trusted.
    let mut buf = Vec::with_capacity(w * h * 3);
    for i in 0..w * h {
        let (x, y) = (i % w, i / w);
        let p = [rgb[i * 3], rgb[i * 3 + 1], rgb[i * 3 + 2]];
        let on_edge = hand[i] && {
            let n = |dx: isize, dy: isize| {
                let (a, b) = (x as isize + dx, y as isize + dy);
                a >= 0 && b >= 0 && (a as usize) < w && (b as usize) < h && hand[b as usize * w + a as usize]
            };
            !(n(1, 0) && n(-1, 0) && n(0, 1) && n(0, -1))
        };
        if on_edge {
            buf.extend_from_slice(&[255, 240, 0]);
        } else if pale[i] {
            buf.extend_from_slice(&[0, 255, 90]);
        } else {
            buf.extend_from_slice(&[p[0] / 3, p[1] / 3, p[2] / 3]);
        }
    }
    let _ = prin_rs::output::adaptive::save_rect(&format!("{out}/circled_mask.png"), w, h, &buf);
    println!("wrote {out}/circled_mask.png -- yellow = hand circles, green = the pale/low-chroma mask\n");

    // ---- the initial conditions, no integration -------------------------------------------
    let ic: Vec<Ics> = (0..w * h)
        .into_par_iter()
        .map(|k| {
            let (x, y) = sl.decode_pos(k);
            let s = grid::decode_state(&chart, 0, x, y);
            let d = newton::pair_dists(&s.s.r);
            let mut longest = 0usize;
            let mut tightest = 0usize;
            for j in 1..3 {
                if d[j] > d[longest] { longest = j; }
                if d[j] < d[tightest] { tightest = j; }
            }
            // `alpha`/`beta` are recoverable from the chart coordinate: the decoder maps
            // `(z_alpha, z_beta)` through `angles()`, and the slice adds the cell offset.
            let (alpha, beta) = prin_rs::physics::decoder::angles(Z[1] + y, Z[0] + x);
            let m = s.m;
            let e0 = energy::energy(&s.s.r, &s.s.v, &m, 0.0);
            let kin: f64 = (0..3).map(|i| 0.5 * m[i] * s.s.v[i].norm_sq()).sum();
            let lz: f64 = (0..3).map(|i| m[i] * (s.s.r[i].x * s.s.v[i].y - s.s.r[i].y * s.s.v[i].x)).sum();
            let vmax = (0..3).map(|i| s.s.v[i].norm()).fold(0.0, f64::max);
            // Twice the signed triangle area -- zero at a collinear configuration.
            let (a, b, c) = (s.s.r[0], s.s.r[1], s.s.r[2]);
            let area = ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)).abs();
            let dmin0 = d.iter().cloned().fold(f64::INFINITY, f64::min);
            let dmax0 = d.iter().cloned().fold(0.0, f64::max);
            let mut ds = d;
            ds.sort_by(|a, b| a.partial_cmp(b).unwrap());
            Ics {
                tie_long: ds[1] / ds[2],
                tie_tight: ds[0] / ds[1],
                alpha, beta, d, longest, tightest,
                reference: THIRD[longest],
                dmin0, dmax0, aspect: dmax0 / dmin0,
                e0, kinetic: kin, virial: kin / (e0 - kin).abs(),
                lz, vmax, area,
            }
        })
        .collect();

    let report = |name: &str, mask: &[bool]| {
        let idx: Vec<usize> = (0..w * h).filter(|&i| mask[i]).collect();
        if idx.is_empty() {
            println!("{name:>22}  EMPTY");
            return;
        }
        let n = idx.len();
        let col = |f: &dyn Fn(&Ics) -> f64| {
            let mut v: Vec<f64> = idx.iter().map(|&i| f(&ic[i])).filter(|x| x.is_finite()).collect();
            (q(&mut v.clone(), 0.1), q(&mut v.clone(), 0.5), q(&mut v, 0.9))
        };
        let frac = |f: &dyn Fn(&Ics) -> bool| idx.iter().filter(|&&i| f(&ic[i])).count() as f64 / n as f64;
        let (a1, a5, a9) = col(&|c| c.alpha);
        let (b1, b5, b9) = col(&|c| c.beta);
        let (m1, m5, m9) = col(&|c| c.dmin0);
        let (r1, r5, r9) = col(&|c| c.aspect);
        let (v1, v5, v9) = col(&|c| c.virial);
        let (e1, e5, e9) = col(&|c| c.e0);
        let (l1, l5, l9) = col(&|c| c.lz);
        let (ar1, ar5, ar9) = col(&|c| c.area);
        let (t1, t5, t9) = col(&|c| c.tie_long);
        let (g1, g5, g9) = col(&|c| c.tie_tight);
        // The three separations directly, rather than through an argmax class. An argmax is a
        // discretisation of these, and a class enrichment can be produced by a shift in any of
        // them -- so the continuous form is what says WHICH.
        let (p01a, p01b, p01c) = col(&|c| c.d[0]);
        let (p02a, p02b, p02c) = col(&|c| c.d[1]);
        let (p12a, p12b, p12c) = col(&|c| c.d[2]);
        println!(
            "{name:>22} n={n:>7} ({:.4} of frame)\n\
             {:>22}   alpha    p10 {a1:.4} p50 {a5:.4} p90 {a9:.4}\n\
             {:>22}   beta     p10 {b1:.4} p50 {b5:.4} p90 {b9:.4}\n\
             {:>22}   d_min0   p10 {m1:.4} p50 {m5:.4} p90 {m9:.4}\n\
             {:>22}   aspect   p10 {r1:.3} p50 {r5:.3} p90 {r9:.3}\n\
             {:>22}   |area|   p10 {ar1:.4} p50 {ar5:.4} p90 {ar9:.4}\n\
             {:>22}   E0       p10 {e1:.4} p50 {e5:.4} p90 {e9:.4}\n\
             {:>22}   K/|U|    p10 {v1:.4} p50 {v5:.4} p90 {v9:.4}\n\
             {:>22}   Lz       p10 {l1:.4} p50 {l5:.4} p90 {l9:.4}\n\
             {:>22}   d(0,1)   p10 {p01a:.4} p50 {p01b:.4} p90 {p01c:.4}\n\
             {:>22}   d(0,2)   p10 {p02a:.4} p50 {p02b:.4} p90 {p02c:.4}\n\
             {:>22}   d(1,2)   p10 {p12a:.4} p50 {p12b:.4} p90 {p12c:.4}\n\
             {:>22}   tie_long p10 {t1:.4} p50 {t5:.4} p90 {t9:.4}   frac>0.95 {:.4}  >0.99 {:.4}\n\
             {:>22}   tie_tght p10 {g1:.4} p50 {g5:.4} p90 {g9:.4}   frac>0.95 {:.4}  >0.99 {:.4}\n\
             {:>22}   ref body 0:{:.4} 1:{:.4} 2:{:.4}   tightest pair 0:{:.4} 1:{:.4} 2:{:.4}",
            n as f64 / (w * h) as f64, "", "", "", "", "", "", "", "", "", "", "", "",
            frac(&|c| c.tie_long > 0.95), frac(&|c| c.tie_long > 0.99), "",
            frac(&|c| c.tie_tight > 0.95), frac(&|c| c.tie_tight > 0.99), "",
            frac(&|c| c.reference == 0), frac(&|c| c.reference == 1), frac(&|c| c.reference == 2),
            frac(&|c| c.tightest == 0), frac(&|c| c.tightest == 1), frac(&|c| c.tightest == 2),
        );
    };

    // The pale mask includes the fine striation everywhere, which is not what was circled. The
    // circled features are the DENSE, contiguous patches of it: pale with a pale majority in a
    // 9x9 neighbourhood. Reported beside the raw mask so the thresholding is visible.
    let mut dense = vec![false; w * h];
    for y in 4..h - 4 {
        for x in 4..w - 4 {
            if !pale[y * w + x] {
                continue;
            }
            let mut c = 0;
            for dy in 0..9 {
                for dx in 0..9 {
                    if pale[(y + dy - 4) * w + (x + dx - 4)] {
                        c += 1;
                    }
                }
            }
            dense[y * w + x] = c * 4 > 81; // >=25% pale in a 9x9
        }
    }

    // ---- the joint class, mapped over the WHOLE frame -----------------------------------
    // `PAIRS = [(0,1),(0,2),(1,2)]`, `THIRD = [2,1,0]`, so `reference = THIRD[longest]`. The
    // class that matters below is `reference 0` (longest side is (1,2)) with `tightest (0,1)`.
    // Mapping it over the whole frame and rendering it is the test: if the class map reproduces
    // the circled wedges, the population IS the class and nothing further needs describing.
    let mut cbuf = Vec::with_capacity(w * h * 3);
    for c in &ic {
        // hue by reference body, lightness by tightest pair
        let base: [u8; 3] = match c.reference {
            0 => [235, 70, 70],
            1 => [70, 110, 235],
            _ => [70, 200, 110],
        };
        let k = 0.45 + 0.275 * c.tightest as f64;
        cbuf.extend_from_slice(&[
            (base[0] as f64 * k) as u8,
            (base[1] as f64 * k) as u8,
            (base[2] as f64 * k) as u8,
        ]);
    }
    let _ = prin_rs::output::adaptive::save_rect(&format!("{out}/ic_class.png"), w, h, &cbuf);
    println!("wrote {out}/ic_class.png -- hue = AZ reference body at t=0, lightness = tightest pair\n");

    let joint: Vec<bool> = ic.iter().map(|c| c.reference == 0 && c.tightest == 0).collect();
    let jn = joint.iter().filter(|x| **x).count();
    let inside = |m: &[bool]| {
        let t = m.iter().filter(|x| **x).count().max(1);
        (0..w * h).filter(|&i| m[i] && joint[i]).count() as f64 / t as f64
    };
    println!(
        "JOINT CLASS  reference=0 (longest side is (1,2))  AND  tightest=(0,1)\n\
         {:>22} {:.4} of the frame\n\
         {:>22} {:.4} of the hand-circled set        enrichment {:.2}x\n\
         {:>22} {:.4} of the pale/low-chroma set     enrichment {:.2}x\n\
         {:>22} {:.4} of the DENSE pale wedges       enrichment {:.2}x\n\
         {:>22} {:.4} of the magenta set             enrichment {:.2}x\n",
        "", jn as f64 / (w * h) as f64,
        "", inside(&hand), inside(&hand) / (jn as f64 / (w * h) as f64),
        "", inside(&pale), inside(&pale) / (jn as f64 / (w * h) as f64),
        "", inside(&dense), inside(&dense) / (jn as f64 / (w * h) as f64),
        "", inside(&magenta), inside(&magenta) / (jn as f64 / (w * h) as f64),
    );

    let all = vec![true; w * h];
    let outside: Vec<bool> = (0..w * h).map(|i| !hand[i] && !pale[i]).collect();
    println!("== INITIAL CONDITIONS, NO INTEGRATION ==\n");
    report("whole frame", &all);
    println!();
    report("hand-circled", &hand);
    println!();
    report("pale/low-chroma", &pale);
    println!();
    report("DENSE pale (the wedges)", &dense);
    println!();
    report("magenta (non-finite)", &magenta);
    println!();
    report("everything else", &outside);
    println!();

    // Agreement between the two selections -- if the property mask lands inside the circles,
    // the hand digitisation is not doing the work.
    let both = (0..w * h).filter(|&i| hand[i] && pale[i]).count();
    let hn = hand.iter().filter(|x| **x).count();
    let pn = pale.iter().filter(|x| **x).count();
    println!(
        "mask agreement: hand {hn}, pale {pn}, both {both}  --  \
         {:.4} of the pale mask lies inside a circle ({:.4} expected by area), enrichment {:.2}x",
        both as f64 / pn.max(1) as f64,
        hn as f64 / (w * h) as f64,
        (both as f64 / pn.max(1) as f64) / (hn as f64 / (w * h) as f64),
    );
    println!(
        "magenta: {} of {} ({:.4}); inside a circle {:.4}, enrichment {:.2}x",
        magenta.iter().filter(|x| **x).count(), w * h,
        magenta.iter().filter(|x| **x).count() as f64 / (w * h) as f64,
        (0..w * h).filter(|&i| magenta[i] && hand[i]).count() as f64
            / magenta.iter().filter(|x| **x).count().max(1) as f64,
        ((0..w * h).filter(|&i| magenta[i] && hand[i]).count() as f64
            / magenta.iter().filter(|x| **x).count().max(1) as f64)
            / (hn as f64 / (w * h) as f64),
    );
}
