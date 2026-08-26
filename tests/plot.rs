//! Tests for the figure writer.
//!
//! Two of these fire on faults the previous figures actually had: a NaN drawn as a valid datum
//! at the top edge, and every series in `far` drawn on top of one another so fifteen of
//! seventeen were invisible.

use prin_rs::output::plot::{all_zero, palette, Figure, Series};

fn s(label: &str, ys: &[f64]) -> Series {
    Series::new(
        label,
        ys.iter().enumerate().map(|(i, &y)| ((i + 1) as f64 * 10.0, y)).collect(),
        (200, 100, 100),
    )
}

#[test]
fn a_non_finite_point_is_counted_as_dropped_and_never_as_live() {
    // The old `ly()` sent NaN through `f64::clamp` into `as isize`, which saturates to 0, so a
    // NaN was drawn at the TOP of the frame as though it were the largest value in the series.
    let a = s("term_grad", &[1e-2, f64::NAN, f64::NAN, 1e-4, 0.0]);
    assert_eq!(a.census(), (2, 1, 2), "(finite positive, exact zero, dropped)");

    // The control: a series with no NaN reports none, so the arm above is not passing by
    // counting everything as dropped.
    let b = s("between", &[1e-2, 1e-3, 1e-4]);
    assert_eq!(b.census(), (3, 0, 0));
}

#[test]
fn an_exact_zero_is_located_and_is_not_confused_with_a_small_value() {
    let a = s("first_div", &[1e-2, 1e-9, 0.0, 0.0]);
    assert_eq!(a.first_zero(), Some(30.0));

    // 1e-9 is small but is NOT zero. Under the old figure both snapped to the same floor pixel
    // and `error(B) = 0` -- the result on several of these curves -- could not be read off.
    let b = s("within", &[1e-2, 1e-9, 1e-12]);
    assert_eq!(b.first_zero(), None);
}

#[test]
fn a_set_of_series_that_are_all_zero_is_detected() {
    // `far` had error(root) = 0.00000, so all 17 series were exactly zero at every budget. The
    // figure drew one flat line; 15 series were completely overdrawn, including the white
    // greedy_lookahead_1 control.
    let z: Vec<Series> = (0..5).map(|i| s(&format!("c{i}"), &[0.0, 0.0, 0.0])).collect();
    assert!(all_zero(&z));

    // The controls, both directions: one live series is enough to make it a real figure, and an
    // empty set is not "all zero" (there is nothing to be zero).
    let mut mixed = z.clone();
    mixed.push(s("live", &[1e-3, 1e-4, 0.0]));
    assert!(!all_zero(&mixed));
    assert!(!all_zero(&[]));

    // A series that is entirely NaN is not zero either -- it has no value at all.
    assert!(!all_zero(&[s("nan", &[f64::NAN, f64::NAN])]));
}

#[test]
fn the_palette_gives_distinct_colours_for_as_many_series_as_asked() {
    // The old table had 8 entries for 17 series and wrapped, so two criteria could not be told
    // apart by colour.
    for n in [3usize, 8, 17, 24] {
        let p = palette(n);
        assert_eq!(p.len(), n);
        let mut worst = u32::MAX;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = (p[i].0 as i32 - p[j].0 as i32).unsigned_abs()
                    + (p[i].1 as i32 - p[j].1 as i32).unsigned_abs()
                    + (p[i].2 as i32 - p[j].2 as i32).unsigned_abs();
                worst = worst.min(d);
            }
        }
        assert!(worst > 20, "palette({n}) has two entries only {worst} apart in sRGB");
    }
}

#[test]
fn a_figure_writes_both_a_png_and_an_svg_and_the_svg_carries_its_labels_as_text() {
    // The SVG is the readable artefact: it stores text as text, so it does not depend on this
    // machine's font resolution. If the labels were rasterised the figure would be as
    // unreadable elsewhere as the 3x5 glyph set was here.
    let dir = std::env::temp_dir().join("prin_plot_test");
    let _ = std::fs::create_dir_all(&dir);
    let stem = dir.join("fig");
    let stem = stem.to_str().unwrap();

    let fig = Figure {
        title: "near-field — error(B)".into(),
        x_label: "budget B (quads computed)".into(),
        y_label: "mean per-pixel OKLab distance".into(),
        series: vec![
            s("between/median", &[1e-2, 3e-3, 1e-3, 0.0]),
            s("term_grad/median", &[1e-2, f64::NAN, f64::NAN, 0.0]),
        ],
        y_lo: 1e-6,
        y_hi: 0.2,
        notes: vec!["a note that must survive to the reader".into()],
    };
    fig.save(stem).unwrap();

    let png = std::fs::metadata(format!("{stem}.png")).unwrap();
    assert!(png.len() > 2000, "png is {} bytes -- suspiciously empty", png.len());

    let svg = std::fs::read_to_string(format!("{stem}.svg")).unwrap();
    for needle in [
        "near-field",
        "budget B (quads computed)",
        "a note that must survive to the reader",
        "between/median",
    ] {
        assert!(svg.contains(needle), "the svg lost the text {needle:?}");
    }
    // The dropped-point count reaches the reader rather than being silently absorbed.
    assert!(
        svg.contains("term_grad/median (2/4)"),
        "the svg does not report the dropped points in the label"
    );
}

#[test]
fn a_degenerate_figure_says_so_instead_of_drawing_a_line() {
    let dir = std::env::temp_dir().join("prin_plot_test");
    let _ = std::fs::create_dir_all(&dir);
    let stem = dir.join("degenerate");
    let stem = stem.to_str().unwrap();

    let fig = Figure {
        title: "far — error(B)".into(),
        x_label: "budget B".into(),
        y_label: "OKLab".into(),
        series: (0..17).map(|i| s(&format!("c{i}"), &[0.0, 0.0, 0.0])).collect(),
        y_lo: 1e-6,
        y_hi: 0.2,
        notes: vec![],
    };
    fig.save(stem).unwrap();

    let svg = std::fs::read_to_string(format!("{stem}.svg")).unwrap();
    assert!(
        svg.contains("every series is exactly 0 at every budget"),
        "the degenerate case was drawn as a curve"
    );
    assert!(
        !svg.contains("mean per-pixel"),
        "the degenerate panel should not carry a y axis, which would imply a measurement"
    );
}
