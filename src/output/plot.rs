//! Line figures for the `error(B)` curves, on `plotters`.
//!
//! The previous version was hand-rolled: a 3x5 dot-matrix glyph set at 1:1 with no case
//! distinction and blanks for unsupported characters, so `random[1]` rendered as `random`. That
//! is what the dependency buys. The dependency does **not** buy the three things that were
//! actually wrong with the figures, and those are handled here.
//!
//! # 1. A NaN is dropped and the drop is reported, never drawn
//!
//! The old `ly()` fed NaN through `f64::clamp`, which propagates it, into `as isize`, which
//! saturates to `0` — so **a NaN point was silently drawn at the top edge as a valid datum**.
//! `term_grad` is NaN on 97.1% of `near-field` quads, so those curves had spurious lines pinned
//! to the top of the frame. Here a non-finite point is dropped and the series label carries its
//! live count (`term_grad (4/11)`), because a criterion that declines to answer must not look
//! like one answering confidently. A high NaN fraction is a property to read — `term_grad` is
//! NaN on 97.1% of near-field and still reaches the oracle's zero by `B = 383`, so the 2.9% it
//! scores are the right quads.
//!
//! # 2. An exact zero gets its own band, not the floor
//!
//! Zero is not representable on a log axis. Snapping it to the bottom of the log panel makes it
//! indistinguishable from a small finite value, and `error(B) = 0` is the most important point
//! on several of these curves. The figure is split: a log panel above, and below a rule, a
//! **zero band** in which each series that reaches zero gets its own row. Rows rather than one
//! shared line because `curve_far_t13.png` had all 17 series at exactly zero at every budget,
//! **15 of them completely overdrawn** — including the white `greedy_oracle` — so the figure
//! showed one flat line and looked like a plot of one criterion.
//!
//! # 3. A figure that cannot distinguish its inputs says so
//!
//! If every series is identically zero the figure carries a stated panel rather than a picture.
//! `far` has `error(root) = 0.00000`: there is no measurable image there at `512^2`, and a
//! figure that renders it as a line is the plotting equivalent of a test that cannot fail.
//!
//! # Output
//!
//! Both PNG and SVG, from one description. The SVG is the readable artefact — it stores text as
//! text, so it renders in whatever face the viewer has and does not depend on this machine's
//! font resolution. The PNG is the raster twin, for folding into an APNG.

use std::error::Error;

use plotters::coord::Shift;
use plotters::prelude::*;

use crate::output::oklab;

/// A named series. `points` may contain non-finite `y`; they are dropped at draw time and
/// counted into the label.
#[derive(Clone, Debug)]
pub struct Series {
    pub label: String,
    pub points: Vec<(f64, f64)>,
    pub rgb: (u8, u8, u8),
    /// Drawn dashed. Used for the controls, so a reader can tell a reference from a candidate.
    pub dashed: bool,
}

impl Series {
    pub fn new(label: impl Into<String>, points: Vec<(f64, f64)>, rgb: (u8, u8, u8)) -> Self {
        Series { label: label.into(), points, rgb, dashed: false }
    }
    pub fn dashed(mut self) -> Self {
        self.dashed = true;
        self
    }
    /// `(finite_positive, exact_zero, dropped_non_finite)`.
    pub fn census(&self) -> (usize, usize, usize) {
        let mut pos = 0;
        let mut zero = 0;
        let mut nan = 0;
        for &(_, y) in &self.points {
            if !y.is_finite() {
                nan += 1;
            } else if y > 0.0 {
                pos += 1;
            } else {
                zero += 1;
            }
        }
        (pos, zero, nan)
    }
    /// The smallest `x` at which this series is exactly zero, if any.
    pub fn first_zero(&self) -> Option<f64> {
        self.points
            .iter()
            .filter(|&&(_, y)| y.is_finite() && y <= 0.0)
            .map(|&(x, _)| x)
            .fold(None, |a: Option<f64>, x| Some(a.map_or(x, |b| b.min(x))))
    }
    /// Label with the live-point count appended when any point was dropped.
    fn display_label(&self) -> String {
        let (pos, zero, nan) = self.census();
        if nan == 0 {
            self.label.clone()
        } else {
            format!("{} ({}/{})", self.label, pos + zero, pos + zero + nan)
        }
    }
}

/// One figure: a title, axis names, and the series to draw.
pub struct Figure {
    pub title: String,
    pub x_label: String,
    pub y_label: String,
    pub series: Vec<Series>,
    /// Log-panel y range. Exact zeros go to the zero band regardless.
    pub y_lo: f64,
    pub y_hi: f64,
    /// Extra lines under the title — the parameters a reader needs to not misread the figure.
    pub notes: Vec<String>,
}

const BG: RGBColor = RGBColor(18, 18, 22);
const FG: RGBColor = RGBColor(215, 215, 225);
const AXIS: RGBColor = RGBColor(120, 120, 132);
const GRID: RGBColor = RGBColor(46, 46, 54);

/// `n` perceptually distinct colours, from OKLCh rather than an eight-entry table.
///
/// The old palette had 8 entries for 17 series and repeated, so two criteria could not be told
/// apart by colour at all. Hue is spread evenly and lightness alternates, which separates
/// adjacent entries on a second axis when `n` is large enough that hue alone is not enough.
pub fn palette(n: usize) -> Vec<(u8, u8, u8)> {
    (0..n)
        .map(|i| {
            let h = 360.0 * (i as f64) / (n.max(1) as f64) + 20.0;
            let l = if i % 2 == 0 { 0.78 } else { 0.62 };
            let c = if i % 3 == 0 { 0.15 } else { 0.12 };
            let r = h.to_radians();
            let [x, y, z] = [l, c * r.cos(), c * r.sin()];
            let v = oklab::oklab_to_srgb([x, y, z]);
            (v[0], v[1], v[2])
        })
        .collect()
}

/// Whether every series is exactly zero wherever it has a value.
///
/// A figure in this state is not a picture of a comparison; the caller is told rather than
/// handed a line. `far` is the case on record.
pub fn all_zero(series: &[Series]) -> bool {
    !series.is_empty()
        && series.iter().all(|s| {
            let (pos, zero, _) = s.census();
            pos == 0 && zero > 0
        })
}

impl Figure {
    /// Write `{stem}.png` and `{stem}.svg`.
    pub fn save(&self, stem: &str) -> Result<(), Box<dyn Error>> {
        let (w, h) = (1400u32, 800u32);
        if let Some(dir) = std::path::Path::new(stem).parent() {
            let _ = std::fs::create_dir_all(dir);
        }
        {
            // `with_buffer` rather than `BitMapBackend::new`: writing a file directly needs
            // plotters' `image` feature, which pulls the whole `image` crate to encode a PNG
            // this repo already has an encoder for. Draw into our own buffer and hand it to
            // `adaptive::save_rect`, which is the writer every other image here goes through.
            let mut buf = vec![0u8; (w as usize) * (h as usize) * 3];
            {
                let root = BitMapBackend::with_buffer(&mut buf, (w, h)).into_drawing_area();
                self.render(root)?;
            }
            crate::output::adaptive::save_rect(
                &format!("{stem}.png"),
                w as usize,
                h as usize,
                &buf,
            )?;
        }
        {
            let svg_path = format!("{stem}.svg");
            let root = SVGBackend::new(&svg_path, (w, h)).into_drawing_area();
            self.render(root)?;
        }
        Ok(())
    }

    fn render<DB: DrawingBackend>(&self, root: DrawingArea<DB, Shift>) -> Result<(), Box<dyn Error>>
    where
        DB::ErrorType: 'static,
    {
        root.fill(&BG)?;

        // Header: title plus the parameters, so the figure is not read without them.
        let head_h = 26 + 18 * self.notes.len() as i32;
        let (head, body) = root.split_vertically(head_h as u32);
        head.draw_text(&self.title, &("sans-serif", 20).into_font().color(&FG), (14, 4))?;
        for (i, note) in self.notes.iter().enumerate() {
            head.draw_text(
                note,
                &("sans-serif", 13).into_font().color(&AXIS),
                (14, 26 + 18 * i as i32),
            )?;
        }

        if all_zero(&self.series) {
            // Not a picture. Every series is exactly zero at every budget, so there is nothing
            // to compare and a line would be read as a result.
            body.draw_text(
                "every series is exactly 0 at every budget",
                &("sans-serif", 30).into_font().color(&FG),
                (60, 120),
            )?;
            for (i, line) in [
                "There is no measurable image here at this resolution: error(root) = 0.",
                "The reference tree and every budgeted tree render the same pixels, so the",
                "criteria are not distinguishable and no ordering read from this region means",
                "anything. This panel is drawn instead of a curve because a flat line at the",
                "floor looks like a result and is not one.",
            ]
            .iter()
            .enumerate()
            {
                body.draw_text(
                    line,
                    &("sans-serif", 15).into_font().color(&FG.mix(0.75)),
                    (60, 170 + 22 * i as i32),
                )?;
            }
            root.present()?;
            return Ok(());
        }

        // Zero band height: one row per series that reaches zero, plus the rule.
        let zero_rows: Vec<&Series> = self.series.iter().filter(|s| s.first_zero().is_some()).collect();
        let band_h = if zero_rows.is_empty() { 0 } else { 58 + 14 * zero_rows.len() as u32 };
        let body_h = body.dim_in_pixel().1;
        let (top, band) = body.split_vertically(body_h.saturating_sub(band_h));

        // x runs from the smallest budget present, not from 1. The old figure anchored the axis
        // at v = 1 and left 165 px permanently blank.
        let xs: Vec<f64> = self.series.iter().flat_map(|s| s.points.iter().map(|p| p.0)).collect();
        let x_lo = xs.iter().cloned().fold(f64::INFINITY, f64::min).max(1.0);
        let x_hi = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max).max(x_lo * 2.0);

        let mut chart = ChartBuilder::on(&top)
            .margin(10)
            .margin_right(260)
            .x_label_area_size(if band_h == 0 { 42 } else { 8 })
            .y_label_area_size(72)
            .build_cartesian_2d(
                (x_lo..x_hi).log_scale(),
                (self.y_lo..self.y_hi).log_scale(),
            )?;

        {
            // When the zero band is present it owns the x axis. Drawing it twice puts two rows
            // of identical tick labels a few pixels apart and reads as two different axes.
            let mut mesh = chart.configure_mesh();
            if band_h == 0 {
                mesh.x_desc(self.x_label.clone());
            } else {
                mesh.disable_x_axis();
            }
            // Budgets are quad counts. A log axis defaults to one decimal, so `383` printed as
            // `383.0` and the ladder read as a continuous quantity rather than a count.
            mesh.x_label_formatter(&|v: &f64| format!("{}", v.round() as i64))
                .y_desc(&self.y_label)
                .label_style(("sans-serif", 13).into_font().color(&AXIS))
                .axis_style(AXIS)
                .bold_line_style(GRID)
                .light_line_style(BG)
                .draw()?;
        }

        for s in &self.series {
            let c = RGBColor(s.rgb.0, s.rgb.1, s.rgb.2);
            let style = ShapeStyle::from(&c).stroke_width(2);
            // Only finite positive points reach the log panel. Non-finite are DROPPED, not
            // clamped; exact zeros go to the band below.
            let pts: Vec<(f64, f64)> = s
                .points
                .iter()
                .filter(|&&(_, y)| y.is_finite() && y > 0.0)
                .cloned()
                .collect();
            if pts.is_empty() {
                // Still needs a legend entry, or a criterion that only ever reached zero would
                // silently vanish from the figure.
                chart
                    .draw_series(std::iter::empty::<Circle<(f64, f64), i32>>())?
                    .label(s.display_label())
                    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], c.stroke_width(2)));
                continue;
            }
            if s.dashed {
                chart
                    .draw_series(DashedLineSeries::new(pts.iter().cloned(), 8, 6, style))?
                    .label(s.display_label())
                    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], c.stroke_width(2)));
            } else {
                chart
                    .draw_series(LineSeries::new(pts.iter().cloned(), style))?
                    .label(s.display_label())
                    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], c.stroke_width(2)));
            }
        }

        chart
            .configure_series_labels()
            .position(SeriesLabelPosition::UpperRight)
            .margin(6)
            .legend_area_size(24)
            .label_font(("sans-serif", 13).into_font().color(&FG))
            .background_style(BG.mix(0.85))
            .border_style(AXIS)
            .draw()?;

        // ---- the zero band ----
        if band_h > 0 {
            let mut zc = ChartBuilder::on(&band)
                .margin_left(10)
                .margin_right(260)
                .margin_top(20)
                .x_label_area_size(34)
                .y_label_area_size(72)
                .build_cartesian_2d((x_lo..x_hi).log_scale(), 0f64..(zero_rows.len() as f64))?;
            zc.configure_mesh()
                .disable_y_mesh()
                .disable_y_axis()
                .x_label_formatter(&|v: &f64| format!("{}", v.round() as i64))
                .x_desc(&self.x_label)
                .label_style(("sans-serif", 13).into_font().color(&AXIS))
                .axis_style(AXIS)
                .bold_line_style(GRID)
                .light_line_style(BG)
                .draw()?;
            band.draw_text(
                "error(B) = 0 exactly",
                &("sans-serif", 12).into_font().color(&FG),
                (14, 2),
            )?;

            // One row per series, so overlapping zeros are all visible. This is the fix for
            // `far`, where 15 of 17 series were drawn on top of each other.
            for (row, s) in zero_rows.iter().enumerate() {
                let y = zero_rows.len() as f64 - 0.5 - row as f64;
                let x0 = s.first_zero().unwrap().max(x_lo);
                let c = RGBColor(s.rgb.0, s.rgb.1, s.rgb.2);
                zc.draw_series(LineSeries::new(
                    vec![(x0, y), (x_hi, y)],
                    ShapeStyle::from(&c).stroke_width(3),
                ))?;
                zc.draw_series(std::iter::once(Circle::new(
                    (x0, y),
                    4,
                    ShapeStyle::from(&c).filled(),
                )))?;
            }
        }

        root.present()?;
        Ok(())
    }
}
