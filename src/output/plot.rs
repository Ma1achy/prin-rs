//! A minimal line plot, written to RGB8.
//!
//! Deliberately small: there is no plotting crate here and adding one to draw a dozen curves
//! would be a dependency for a figure. This does axes, a log-y grid, and labelled series, which
//! is what an `error(B)` curve needs and nothing more.
//!
//! **Log y, always.** The curves in this project span four decades and the interesting part is
//! the bottom one — a linear axis puts every good criterion on the same flat line at the floor
//! and makes them indistinguishable, which is the opposite of what the figure is for. Zero is
//! not representable on a log axis, so an exact zero is drawn **at the floor and marked**, never
//! silently dropped: `error(B) = 0` is the most important point on several of these curves.

/// A named series of `(x, y)` points.
pub struct Series<'a> {
    pub label: &'a str,
    pub points: Vec<(f64, f64)>,
    pub rgb: [u8; 3],
    /// Drawn as a dashed line. Used for the controls, so a reader can tell at a glance which
    /// curves are references rather than candidates.
    pub dashed: bool,
}

pub struct Plot {
    pub w: usize,
    pub h: usize,
    px: Vec<u8>,
    bg: [u8; 3],
}

const PAD_L: usize = 64;
const PAD_R: usize = 150;
const PAD_T: usize = 28;
const PAD_B: usize = 40;

impl Plot {
    pub fn new(w: usize, h: usize) -> Self {
        let bg = [18, 18, 22];
        Plot { w, h, px: bg.iter().cycle().take(w * h * 3).cloned().collect(), bg }
    }

    fn set(&mut self, x: isize, y: isize, rgb: [u8; 3]) {
        if x < 0 || y < 0 || x as usize >= self.w || y as usize >= self.h {
            return;
        }
        let o = (y as usize * self.w + x as usize) * 3;
        self.px[o] = rgb[0];
        self.px[o + 1] = rgb[1];
        self.px[o + 2] = rgb[2];
    }

    fn line(&mut self, x0: isize, y0: isize, x1: isize, y1: isize, rgb: [u8; 3], dashed: bool) {
        let (dx, dy) = ((x1 - x0).abs(), -(y1 - y0).abs());
        let (sx, sy) = (if x0 < x1 { 1 } else { -1 }, if y0 < y1 { 1 } else { -1 });
        let (mut x, mut y, mut err) = (x0, y0, dx + dy);
        let mut n = 0usize;
        loop {
            if !dashed || (n / 4) % 2 == 0 {
                // 2px wide, so a curve is legible against the grid at a glance.
                self.set(x, y, rgb);
                self.set(x, y + 1, rgb);
            }
            n += 1;
            if x == x1 && y == y1 {
                break;
            }
            let e2 = 2 * err;
            if e2 >= dy {
                err += dy;
                x += sx;
            }
            if e2 <= dx {
                err += dx;
                y += sy;
            }
        }
    }

    /// A 3x5 dot-matrix glyph set — enough for labels and axis numbers.
    fn glyph(c: char) -> [u8; 5] {
        match c.to_ascii_lowercase() {
            '0' => [0b111, 0b101, 0b101, 0b101, 0b111],
            '1' => [0b010, 0b110, 0b010, 0b010, 0b111],
            '2' => [0b111, 0b001, 0b111, 0b100, 0b111],
            '3' => [0b111, 0b001, 0b111, 0b001, 0b111],
            '4' => [0b101, 0b101, 0b111, 0b001, 0b001],
            '5' => [0b111, 0b100, 0b111, 0b001, 0b111],
            '6' => [0b111, 0b100, 0b111, 0b101, 0b111],
            '7' => [0b111, 0b001, 0b010, 0b010, 0b010],
            '8' => [0b111, 0b101, 0b111, 0b101, 0b111],
            '9' => [0b111, 0b101, 0b111, 0b001, 0b111],
            'a' => [0b010, 0b101, 0b111, 0b101, 0b101],
            'b' => [0b110, 0b101, 0b110, 0b101, 0b110],
            'c' => [0b011, 0b100, 0b100, 0b100, 0b011],
            'd' => [0b110, 0b101, 0b101, 0b101, 0b110],
            'e' => [0b111, 0b100, 0b110, 0b100, 0b111],
            'f' => [0b111, 0b100, 0b110, 0b100, 0b100],
            'g' => [0b011, 0b100, 0b101, 0b101, 0b011],
            'h' => [0b101, 0b101, 0b111, 0b101, 0b101],
            'i' => [0b111, 0b010, 0b010, 0b010, 0b111],
            'j' => [0b001, 0b001, 0b001, 0b101, 0b010],
            'k' => [0b101, 0b110, 0b100, 0b110, 0b101],
            'l' => [0b100, 0b100, 0b100, 0b100, 0b111],
            'm' => [0b101, 0b111, 0b111, 0b101, 0b101],
            'n' => [0b101, 0b111, 0b111, 0b111, 0b101],
            'o' => [0b010, 0b101, 0b101, 0b101, 0b010],
            'p' => [0b110, 0b101, 0b110, 0b100, 0b100],
            'q' => [0b010, 0b101, 0b101, 0b111, 0b011],
            'r' => [0b110, 0b101, 0b110, 0b101, 0b101],
            's' => [0b011, 0b100, 0b010, 0b001, 0b110],
            't' => [0b111, 0b010, 0b010, 0b010, 0b010],
            'u' => [0b101, 0b101, 0b101, 0b101, 0b111],
            'v' => [0b101, 0b101, 0b101, 0b101, 0b010],
            'w' => [0b101, 0b101, 0b111, 0b111, 0b101],
            'x' => [0b101, 0b101, 0b010, 0b101, 0b101],
            'y' => [0b101, 0b101, 0b010, 0b010, 0b010],
            'z' => [0b111, 0b001, 0b010, 0b100, 0b111],
            '-' => [0b000, 0b000, 0b111, 0b000, 0b000],
            '_' => [0b000, 0b000, 0b000, 0b000, 0b111],
            '/' => [0b001, 0b001, 0b010, 0b100, 0b100],
            '.' => [0b000, 0b000, 0b000, 0b000, 0b010],
            ':' => [0b000, 0b010, 0b000, 0b010, 0b000],
            '=' => [0b000, 0b111, 0b000, 0b111, 0b000],
            '^' => [0b010, 0b101, 0b000, 0b000, 0b000],
            '(' => [0b001, 0b010, 0b010, 0b010, 0b001],
            ')' => [0b100, 0b010, 0b010, 0b010, 0b100],
            _ => [0; 5],
        }
    }

    pub fn text(&mut self, x: usize, y: usize, s: &str, rgb: [u8; 3]) {
        for (i, c) in s.chars().enumerate() {
            let g = Self::glyph(c);
            for (r, row) in g.iter().enumerate() {
                for col in 0..3 {
                    if row & (1 << (2 - col)) != 0 {
                        self.set((x + i * 4 + col) as isize, (y + r) as isize, rgb);
                    }
                }
            }
        }
    }

    /// Draw the axes and every series. `y_floor` is where an exact zero is drawn.
    pub fn draw(&mut self, title: &str, series: &[Series], x_max: f64, y_lo: f64, y_hi: f64) {
        let axis = [90, 90, 100];
        let grid = [40, 40, 48];
        let fg = [210, 210, 220];

        let plot_w = self.w - PAD_L - PAD_R;
        let plot_h = self.h - PAD_T - PAD_B;

        // log10 x (budget spans 5 -> 5461) and log10 y.
        let lx = |v: f64| -> isize {
            let t = (v.max(1.0).log10()) / (x_max.max(10.0).log10());
            (PAD_L as f64 + t * plot_w as f64) as isize
        };
        let ly = |v: f64| -> isize {
            // An exact zero is not representable on a log axis. It is drawn at the floor and
            // marked, never dropped: it is the most important point on several of these curves.
            let vv = if v <= y_lo { y_lo } else { v };
            let t = (vv.log10() - y_lo.log10()) / (y_hi.log10() - y_lo.log10());
            (PAD_T as f64 + (1.0 - t.clamp(0.0, 1.0)) * plot_h as f64) as isize
        };

        // decade grid
        let mut d = y_lo.log10().ceil() as i32;
        while (d as f64) <= y_hi.log10() {
            let y = ly(10f64.powi(d));
            for x in PAD_L..(PAD_L + plot_w) {
                self.set(x as isize, y, grid);
            }
            self.text(6, (y - 2).max(0) as usize, &format!("1e{d}"), axis);
            d += 1;
        }
        for &b in &[10.0f64, 100.0, 1000.0] {
            if b <= x_max {
                let x = lx(b);
                for y in PAD_T..(PAD_T + plot_h) {
                    self.set(x, y as isize, grid);
                }
                self.text(x as usize - 6, self.h - PAD_B + 8, &format!("{}", b as u64), axis);
            }
        }
        self.text(8, 8, title, fg);
        self.text(PAD_L, self.h - PAD_B + 20, "budget b (quads computed)", axis);
        self.text(6, PAD_T - 14, "error(b) oklab", axis);

        for (i, s) in series.iter().enumerate() {
            let pts: Vec<(isize, isize)> =
                s.points.iter().map(|&(x, y)| (lx(x), ly(y))).collect();
            for w in pts.windows(2) {
                self.line(w[0].0, w[0].1, w[1].0, w[1].1, s.rgb, s.dashed);
            }
            // An exact zero gets a tick at the floor so it reads as reached, not missing.
            for (j, &(_, y)) in s.points.iter().enumerate() {
                if y == 0.0 {
                    let x = pts[j].0;
                    for k in 0..5 {
                        self.set(x, ly(y_lo) - k, s.rgb);
                    }
                    break;
                }
            }
            let ly0 = PAD_T + 6 + i * 10;
            if ly0 + 6 < self.h {
                for k in 0..10 {
                    self.set((self.w - PAD_R + 4 + k) as isize, (ly0 + 2) as isize, s.rgb);
                }
                self.text(self.w - PAD_R + 18, ly0, s.label, fg);
            }
        }
    }

    pub fn save(&self, path: &str) -> std::io::Result<()> {
        let _ = self.bg;
        crate::output::adaptive::save_rect(path, self.w, self.h, &self.px)
    }

    pub fn pixels(&self) -> &[u8] {
        &self.px
    }
}
