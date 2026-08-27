//! Run configuration. A small CLI; no TOML unless it earns its place.

use crate::ensemble::pixel::EnsembleCfg;
use crate::grid::{self, Slice};
use crate::integrate::az::RefPolicy;
use crate::render::Precision;

#[derive(Clone, Debug)]
pub struct Config {
    pub slice: Slice,
    pub ens: EnsembleCfg,
    pub precision: Precision,
    pub out: String,
    pub region: String,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            slice: grid::region("near-field", 64, 64, grid::Chart::BodyPlane.default_half())
                .unwrap(),
            ens: EnsembleCfg::default(),
            precision: Precision::F64,
            out: "out".into(),
            region: "near-field".into(),
        }
    }
}

pub const USAGE: &str = "\
prin — uniform-resolution three-body initial-condition kernel

  --region <name>     one of BRIEF §2.2's regions (default near-field)
  --size <n>          grid is n x n (default 64)
  --half <x>          box half-width (default: the chart's own, 0.05 for the body plane)
  --t-max <t>         playhead (default 13)
  --n-sync <n>        sync boundaries (default 32)
  --eta <x>           timestep parameter (default 0.01)
  --copies <n>        E+1 copies per pixel (default 8)
  --jitter-frac <x>   jitter as a fraction of cell width (default 0.5)
  --seed <n>          per-pixel seed, PCG scheme only (default 0)
  --jitter-pcg        use the reference's per-pixel PCG stream instead of fixed Halton
  --precision <p>     f32 | f64 (default f64)
  --shared-reference  force all copies onto the nominal copy's reference body
  --lc-unstable       use the reference's unconditioned inverse LC branch
  --r-coll <x>        collision radius as a FRACTION of R, fixed at t=0 (default 1e-3)
  --r-esc <x>         escape distance gate as a FRACTION of R, fixed at t=0 (default 5);
                      0 restores the numpy reference's ungated test
  --escape-one-body   test only the body outside the tightest pair (the numpy reference)
  --no-stop           record events but integrate every copy to t_max anyway
  --no-refine         skip the second pass over error_ratio-flagged pixels
  --refine-threshold <x>  error_ratio above which a pixel is re-integrated (default 10)
  --refine-eta <x>    factor applied to eta on each refinement pass (default 0.25)
  --refine-passes <n> maximum refinement passes (default 3)
  --out <stem>        output stem (default out)
";

pub fn parse(args: &[String]) -> Result<Config, String> {
    let mut c = Config::default();
    let get = |i: usize| -> Result<String, String> {
        args.get(i + 1).cloned().ok_or_else(|| format!("{} needs a value", args[i]))
    };
    let mut size = 64usize;
    // No literal here. `half` means something different in every chart family -- a body
    // position in Burrau units on the body plane, a sigmoid pre-image on the latent chart -- and
    // a single shared default silently meant both. See `Chart::default_half`.
    let mut half: Option<f64> = None;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--region" => { c.region = get(i)?; i += 1; }
            "--size" => { size = get(i)?.parse().map_err(|e| format!("{e}"))?; i += 1; }
            "--half" => { half = Some(get(i)?.parse().map_err(|e| format!("{e}"))?); i += 1; }
            "--t-max" => { c.ens.t_max = get(i)?.parse().map_err(|e| format!("{e}"))?; i += 1; }
            "--n-sync" => { c.ens.n_sync = get(i)?.parse().map_err(|e| format!("{e}"))?; i += 1; }
            "--eta" => { c.ens.eta = get(i)?.parse().map_err(|e| format!("{e}"))?; i += 1; }
            "--copies" => {
                let n: usize = get(i)?.parse().map_err(|e| format!("{e}"))?;
                if n < 1 { return Err("--copies must be at least 1".into()); }
                c.ens.n_extra = n - 1;
                i += 1;
            }
            "--jitter-frac" => { c.ens.jitter_frac = get(i)?.parse().map_err(|e| format!("{e}"))?; i += 1; }
            "--seed" => { c.ens.seed = get(i)?.parse().map_err(|e| format!("{e}"))?; i += 1; }
            "--precision" => {
                c.precision = match get(i)?.as_str() {
                    "f32" => Precision::F32,
                    "f64" => Precision::F64,
                    other => return Err(format!("unknown precision: {other}")),
                };
                i += 1;
            }
            "--jitter-pcg" => c.ens.jitter_scheme = crate::ensemble::jitter::Scheme::Pcg,
            "--shared-reference" => c.ens.ref_policy = RefPolicy::Shared,
            "--lc-unstable" => c.ens.lc_stable = false,
            "--r-coll" => { c.ens.r_coll_frac = get(i)?.parse().map_err(|e| format!("{e}"))?; i += 1; }
            "--r-esc" => { c.ens.r_esc_frac = get(i)?.parse().map_err(|e| format!("{e}"))?; i += 1; }
            "--escape-one-body" => c.ens.escape_all_bodies = false,
            "--no-stop" => c.ens.stop_on_event = false,
            "--no-refine" => c.ens.refine_flagged = false,
            "--refine-threshold" => { c.ens.refine_threshold = get(i)?.parse().map_err(|e| format!("{e}"))?; i += 1; }
            "--refine-passes" => { c.ens.refine_max_passes = get(i)?.parse().map_err(|e| format!("{e}"))?; i += 1; }
            "--refine-eta" => { c.ens.refine_eta_factor = get(i)?.parse().map_err(|e| format!("{e}"))?; i += 1; }
            "--out" => { c.out = get(i)?; i += 1; }
            "-h" | "--help" => return Err(USAGE.into()),
            other => return Err(format!("unknown argument: {other}\n\n{USAGE}")),
        }
        i += 1;
    }
    // `grid::region` only ever builds a `BodyPlane` slice, so the fallback is 0.05 today and
    // this changes no `prin` behaviour. The point is that the number now has one home.
    let half = half.unwrap_or_else(|| grid::Chart::BodyPlane.default_half());
    c.slice = grid::region(&c.region, size, size, half)
        .ok_or_else(|| format!("unknown region: {}\n\n{USAGE}", c.region))?;
    Ok(c)
}
