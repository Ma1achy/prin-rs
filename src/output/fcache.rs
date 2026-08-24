//! `PRQF` — the **per-footprint** cache, so a colouring change never costs an integration again.
//!
//! # Why this exists
//!
//! `PRQC` (`qcache.rs`) stores per-**quad** reductions and a baked `err_sum`. `err_sum` is a
//! quad's summed OKLab distance to the reference *were it drawn as a leaf* — which makes the
//! greedy replay a static priority queue, and makes the whole cache a function of the colouring.
//! So when the production colouring changed, every `error(B)` number in PR #13 became a
//! measurement against a target that no longer ships, and there was **no way to recover them
//! short of re-integrating 2.8 million trajectories per region.**
//!
//! This stores what a colouring reads, per footprint. With it, `error(B)` under any colouring —
//! present or future — is a replay of a file rather than a march. The integration is paid once.
//!
//! # What is stored, and what is deliberately not
//!
//! Only the fields a colour map or the reserved-null path reads. Not the ensemble, not the
//! per-copy outcomes, not the boundary shape vectors: those feed the *criterion*, and the
//! criterion's scalars are already in `PRQC`. Splitting it this way keeps each file the size of
//! the question it answers — at 5461 quads x 64 footprints x 14 `f64` this is about 39 MB per
//! region, the same budget `PRQC` was sized against.
//!
//! # The one thing to be careful of
//!
//! A `Row` reconstitutes into a `PixelOut` with everything else at [`Default`]. That is correct
//! for colouring and **wrong for anything else**: a reconstituted footprint has `state = 0`
//! (`Escape`) unless the packed byte says otherwise, empty `copy_outcomes`, and zero drift. It
//! is not a footprint; it is the colour-relevant projection of one. [`Row::to_pixel`] says so at
//! its own call site, and nothing here hands one to `scheduler::reduce`.

use std::collections::HashMap;
use std::io::{self, Read, Write};

use crate::ensemble::pixel::PixelOut;
use crate::metric::Key;

pub const MAGIC: &[u8; 4] = b"PRQF";
pub const VERSION: u32 = 1;

pub const FIELDS: &[&str] = &[
    "shape_x", "shape_y", "shape_z",
    "packed_outcome", "n_nonfinite",
    "ensemble_spread", "spread_shape", "spread_event",
    "ftle", "diffusion",
    "t_end", "error_ratio", "d_min_true", "energy_drift_max",
];

/// The colour-relevant projection of one footprint.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Row {
    pub shape: [f64; 3],
    /// `state << 2 | detail`, the same packing `ssaa::packed_rgb` reads.
    pub packed: u8,
    pub n_nonfinite: u8,
    pub ensemble_spread: f64,
    pub spread_shape: f64,
    pub spread_event: f64,
    pub ftle: f64,
    pub diffusion: f64,
    pub t_end: f64,
    pub error_ratio: f64,
    pub d_min_true: f64,
    pub energy_drift_max: f64,
}

impl Row {
    pub fn of(p: &PixelOut) -> Row {
        Row {
            shape: p.shape_vec,
            packed: (p.state << 2) | (p.detail & 0b11),
            n_nonfinite: p.n_nonfinite,
            ensemble_spread: p.ensemble_spread,
            spread_shape: p.spread_shape,
            spread_event: p.spread_event,
            ftle: p.ftle,
            diffusion: p.diffusion,
            t_end: p.t_end,
            error_ratio: p.error_ratio,
            d_min_true: p.d_min_true,
            energy_drift_max: p.energy_drift_max,
        }
    }

    /// Back to a `PixelOut` **for colouring only**.
    ///
    /// Every field not listed in [`FIELDS`] comes back at its default, so this is a projection
    /// and not a round trip. Handing one of these to anything that reads the ensemble — the
    /// scheduler, `error_ratio`'s neighbourhood, the temporal accumulators — would read zeros
    /// as measurements. `tests/fcache.rs` asserts the colouring round-trips bitwise and that
    /// the projection is *not* claimed to be more than that.
    pub fn to_pixel(&self) -> PixelOut {
        PixelOut {
            shape_vec: self.shape,
            state: self.packed >> 2,
            detail: self.packed & 0b11,
            outcome: self.packed,
            n_nonfinite: self.n_nonfinite,
            ensemble_spread: self.ensemble_spread,
            spread_shape: self.spread_shape,
            spread_event: self.spread_event,
            ftle: self.ftle,
            diffusion: self.diffusion,
            t_end: self.t_end,
            error_ratio: self.error_ratio,
            d_min_true: self.d_min_true,
            energy_drift_max: self.energy_drift_max,
            ..Default::default()
        }
    }

    fn to_f64s(self) -> [f64; 14] {
        [
            self.shape[0],
            self.shape[1],
            self.shape[2],
            self.packed as f64,
            self.n_nonfinite as f64,
            self.ensemble_spread,
            self.spread_shape,
            self.spread_event,
            self.ftle,
            self.diffusion,
            self.t_end,
            self.error_ratio,
            self.d_min_true,
            self.energy_drift_max,
        ]
    }

    fn from_f64s(v: &[f64]) -> Row {
        Row {
            shape: [v[0], v[1], v[2]],
            packed: v[3] as u8,
            n_nonfinite: v[4] as u8,
            ensemble_spread: v[5],
            spread_shape: v[6],
            spread_event: v[7],
            ftle: v[8],
            diffusion: v[9],
            t_end: v[10],
            error_ratio: v[11],
            d_min_true: v[12],
            energy_drift_max: v[13],
        }
    }
}

/// Every footprint of a complete tree, addressed by quad.
#[derive(Clone, Debug)]
pub struct Footprints {
    pub region: String,
    pub chart: String,
    pub cx: f64,
    pub cy: f64,
    pub half: f64,
    pub body: usize,
    pub levels: u32,
    pub n: usize,
    pub res: usize,
    pub t_max: f64,
    pub quads: HashMap<Key, Vec<Row>>,
}

impl Footprints {
    /// The geometry this file describes matches the cache it is about to recolour.
    ///
    /// Checked rather than assumed: a footprint file from a different region or resolution would
    /// recolour without complaint and produce an `error(B)` curve for a tree that was never
    /// integrated. That is the failure mode a self-describing header exists to prevent, and a
    /// header only prevents it if something reads it.
    pub fn agrees_with(&self, c: &crate::metric::Cache) -> Result<(), String> {
        let mismatch = |what: &str, a: String, b: String| Err(format!("{what}: file {a}, cache {b}"));
        if self.region != c.region {
            return mismatch("region", self.region.clone(), c.region.clone());
        }
        if self.levels != c.levels || self.n != c.n || self.res != c.res {
            return mismatch(
                "geometry",
                format!("levels={} n={} res={}", self.levels, self.n, self.res),
                format!("levels={} n={} res={}", c.levels, c.n, c.res),
            );
        }
        for (k, v) in [("cx", (self.cx, c.cx)), ("cy", (self.cy, c.cy)), ("half", (self.half, c.half))] {
            if v.0.to_bits() != v.1.to_bits() {
                return mismatch(k, format!("{:?}", v.0), format!("{:?}", v.1));
            }
        }
        if self.quads.len() != c.quads.len() {
            return mismatch(
                "quad count",
                self.quads.len().to_string(),
                c.quads.len().to_string(),
            );
        }
        Ok(())
    }
}

pub fn write<W: Write>(w: &mut W, f: &Footprints) -> io::Result<()> {
    w.write_all(MAGIC)?;
    w.write_all(&VERSION.to_le_bytes())?;

    let header = format!(
        "region={} body={} cx={:?} cy={:?} half={:?} levels={} n={} res={} t_max={}\n\
         chart={}\n\
         quads={} footprints_per_quad={}\n\
         note=this is the COLOUR-RELEVANT PROJECTION of each footprint, not the footprint. Every \
field not listed below comes back at its default, so a reconstituted PixelOut must never be \
handed to the scheduler or to anything that reads the ensemble.\n\
         note=PRQC stores per-quad reductions and a BAKED err_sum, which is a function of the \
colouring. This file is what makes error(B) under a new colouring a replay rather than a \
re-integration.\n\
         fields={}\n",
        f.region,
        f.body,
        f.cx,
        f.cy,
        f.half,
        f.levels,
        f.n,
        f.res,
        f.t_max,
        f.chart,
        f.quads.len(),
        f.n * f.n,
        FIELDS.join(","),
    );
    let hb = header.as_bytes();
    w.write_all(&(hb.len() as u32).to_le_bytes())?;
    w.write_all(hb)?;

    // Sorted by (level, iy, ix), matching PRQC, so the two files line up and diff.
    let mut keys: Vec<Key> = f.quads.keys().cloned().collect();
    keys.sort_by_key(|&(l, ix, iy)| (l, iy, ix));

    w.write_all(&(keys.len() as u64).to_le_bytes())?;
    w.write_all(&(FIELDS.len() as u32).to_le_bytes())?;
    for k in keys {
        let rows = &f.quads[&k];
        w.write_all(&(k.0).to_le_bytes())?;
        w.write_all(&(k.1).to_le_bytes())?;
        w.write_all(&(k.2).to_le_bytes())?;
        w.write_all(&(rows.len() as u32).to_le_bytes())?;
        for r in rows {
            for v in r.to_f64s() {
                w.write_all(&v.to_le_bytes())?;
            }
        }
    }
    Ok(())
}

pub fn read<R: Read>(r: &mut R) -> io::Result<Footprints> {
    let bad = |m: &str| io::Error::new(io::ErrorKind::InvalidData, m.to_string());

    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if &magic != MAGIC {
        return Err(bad("not a PRQF file"));
    }
    let mut u32b = [0u8; 4];
    r.read_exact(&mut u32b)?;
    let v = u32::from_le_bytes(u32b);
    if v != VERSION {
        return Err(bad(&format!("PRQF version {v}, expected {VERSION}")));
    }
    r.read_exact(&mut u32b)?;
    let hlen = u32::from_le_bytes(u32b) as usize;
    let mut hb = vec![0u8; hlen];
    r.read_exact(&mut hb)?;
    let header = String::from_utf8_lossy(&hb).into_owned();

    let field = |name: &str| -> Option<String> {
        header
            .split_whitespace()
            .find_map(|t| t.strip_prefix(&format!("{name}=")).map(str::to_string))
    };
    // `chart` carries the chart's full parameters and contains spaces, so it gets its own line
    // and is read line-wise. A space-containing value on a shared line silently truncates at the
    // first space and the header stops being self-describing -- caught by the round-trip test,
    // which is why the round trip compares the chart string rather than only the geometry.
    let line_field = |name: &str| -> Option<String> {
        header
            .lines()
            .find_map(|l| l.trim().strip_prefix(&format!("{name}=")).map(str::to_string))
    };
    let num = |name: &str| -> f64 {
        field(name).and_then(|s| s.parse().ok()).unwrap_or(f64::NAN)
    };

    let mut u64b = [0u8; 8];
    r.read_exact(&mut u64b)?;
    let nq = u64::from_le_bytes(u64b) as usize;
    r.read_exact(&mut u32b)?;
    let nf = u32::from_le_bytes(u32b) as usize;
    if nf != FIELDS.len() {
        return Err(bad(&format!("PRQF has {nf} fields, this build expects {}", FIELDS.len())));
    }

    let mut quads: HashMap<Key, Vec<Row>> = HashMap::with_capacity(nq);
    for _ in 0..nq {
        r.read_exact(&mut u32b)?;
        let l = u32::from_le_bytes(u32b);
        r.read_exact(&mut u32b)?;
        let ix = u32::from_le_bytes(u32b);
        r.read_exact(&mut u32b)?;
        let iy = u32::from_le_bytes(u32b);
        r.read_exact(&mut u32b)?;
        let nrows = u32::from_le_bytes(u32b) as usize;
        let mut rows = Vec::with_capacity(nrows);
        let mut buf = vec![0f64; nf];
        for _ in 0..nrows {
            for slot in buf.iter_mut() {
                r.read_exact(&mut u64b)?;
                *slot = f64::from_le_bytes(u64b);
            }
            rows.push(Row::from_f64s(&buf));
        }
        quads.insert((l, ix, iy), rows);
    }

    Ok(Footprints {
        region: field("region").unwrap_or_default(),
        chart: line_field("chart").unwrap_or_default(),
        cx: num("cx"),
        cy: num("cy"),
        half: num("half"),
        body: num("body") as usize,
        levels: num("levels") as u32,
        n: num("n") as usize,
        res: num("res") as usize,
        t_max: num("t_max"),
        quads,
    })
}
