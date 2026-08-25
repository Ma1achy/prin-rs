//! OKLab, for measuring how far two renders are apart.
//!
//! The repo's only image comparison is `ssaa::resolve_effect`, a Chebyshev distance over raw
//! 8-bit sRGB channels. That is fine for "did this pixel move at all", which is what it was
//! written for, and wrong for "how different do these two images look", which is what the §2
//! metric needs: sRGB is perceptually non-uniform, so an equal channel step is a much larger
//! visible change in the dark end than the light end, and a criterion would be scored partly on
//! where in the ramp its errors happened to land.
//!
//! Transcribed from Bjorn Ottosson's published constants. **Ported, not re-derived** — the
//! matrices are the kind of algebra that fails silently, and the round-trip test in
//! `tests/criterion.rs` is what makes a transcription error visible.

/// sRGB 8-bit to linear-light, the standard piecewise curve (not a plain 2.2 power).
fn to_linear(c: u8) -> f64 {
    let x = c as f64 / 255.0;
    if x <= 0.04045 {
        x / 12.92
    } else {
        ((x + 0.055) / 1.055).powf(2.4)
    }
}

fn from_linear(x: f64) -> u8 {
    let x = x.clamp(0.0, 1.0);
    let s = if x <= 0.0031308 {
        x * 12.92
    } else {
        1.055 * x.powf(1.0 / 2.4) - 0.055
    };
    (s * 255.0).round().clamp(0.0, 255.0) as u8
}

/// Linear sRGB to OKLab.
pub fn linear_to_oklab(r: f64, g: f64, b: f64) -> [f64; 3] {
    let l = 0.412_221_470_8 * r + 0.536_332_536_3 * g + 0.051_445_992_9 * b;
    let m = 0.211_903_498_2 * r + 0.680_699_545_1 * g + 0.107_396_956_6 * b;
    let s = 0.088_302_461_9 * r + 0.281_718_837_6 * g + 0.629_978_700_5 * b;

    let (l_, m_, s_) = (l.cbrt(), m.cbrt(), s.cbrt());

    [
        0.210_454_255_3 * l_ + 0.793_617_785_0 * m_ - 0.004_072_046_8 * s_,
        1.977_998_495_1 * l_ - 2.428_592_205_0 * m_ + 0.450_593_709_9 * s_,
        0.025_904_037_1 * l_ + 0.782_771_766_2 * m_ - 0.808_675_766_0 * s_,
    ]
}

/// OKLab back to linear sRGB.
pub fn oklab_to_linear(lab: [f64; 3]) -> (f64, f64, f64) {
    let l_ = lab[0] + 0.396_337_777_4 * lab[1] + 0.215_803_757_3 * lab[2];
    let m_ = lab[0] - 0.105_561_345_8 * lab[1] - 0.063_854_172_8 * lab[2];
    let s_ = lab[0] - 0.089_484_177_5 * lab[1] - 1.291_485_548_0 * lab[2];

    let (l, m, s) = (l_ * l_ * l_, m_ * m_ * m_, s_ * s_ * s_);

    (
        4.076_741_662_1 * l - 3.307_711_591_3 * m + 0.230_969_929_2 * s,
        -1.268_438_004_6 * l + 2.609_757_401_1 * m - 0.341_319_396_5 * s,
        -0.004_196_086_3 * l - 0.703_418_614_7 * m + 1.707_614_701_0 * s,
    )
}

/// sRGB 8-bit triple to OKLab.
pub fn srgb_to_oklab(c: [u8; 3]) -> [f64; 3] {
    linear_to_oklab(to_linear(c[0]), to_linear(c[1]), to_linear(c[2]))
}

/// OKLab back to an sRGB 8-bit triple, clamped into gamut.
pub fn oklab_to_srgb(lab: [f64; 3]) -> [u8; 3] {
    let (r, g, b) = oklab_to_linear(lab);
    [from_linear(r), from_linear(g), from_linear(b)]
}

/// Euclidean distance in OKLab — the space is built so this is a perceptual difference.
pub fn delta(a: [u8; 3], b: [u8; 3]) -> f64 {
    let (x, y) = (srgb_to_oklab(a), srgb_to_oklab(b));
    let d = [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// Per-pixel OKLab distance between two RGB8 buffers: `(mean, p99, max, fraction moved)`.
///
/// The fraction is reported alongside the aggregates on purpose. *"Never conclude 'no effect'
/// from an aggregate without the per-pixel distribution"* — a mean can sit still while every
/// pixel moves, which has been measured twice in this project in a single PR.
pub fn image_error(a: &[u8], b: &[u8]) -> (f64, f64, f64, f64) {
    assert_eq!(a.len(), b.len(), "images must match in size");
    let n = a.len() / 3;
    if n == 0 {
        return (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
    }
    let mut d: Vec<f64> = Vec::with_capacity(n);
    for k in 0..n {
        d.push(delta(
            [a[3 * k], a[3 * k + 1], a[3 * k + 2]],
            [b[3 * k], b[3 * k + 1], b[3 * k + 2]],
        ));
    }
    let moved = d.iter().filter(|&&x| x > 0.0).count() as f64 / n as f64;
    let mean = d.iter().sum::<f64>() / n as f64;
    let max = d.iter().cloned().fold(0.0f64, f64::max);
    let p99 = crate::quad::quantile(&mut d, 0.99);
    (mean, p99, max, moved)
}
