//! Outputs: a PNG pair for looking at, and a raw dump for measuring.

pub mod adaptive;
pub mod png;
pub mod plot;
pub mod qcache;
pub mod raw;
pub mod apng;
pub mod ckpt;
pub mod colour;
pub mod fcache;
pub mod gifout;
pub mod oklab;
pub mod compose;
pub mod palette;
pub mod siteblend;
pub mod ssaa;
pub mod tree;
pub mod viridis;
pub mod wire;

/// Write `<stem>.cfg.txt` beside a rendered panel, naming every departure from
/// [`EnsembleCfg::production`](crate::ensemble::pixel::EnsembleCfg::production).
///
/// **This is where the six-day failure lived.** The `.raw` and `.prnq` dumps have carried a full
/// settings header since they were written; the PNGs carry nothing, and the harnesses that make
/// them printed nothing either. So `refine_flagged: false` propagated by copy through five
/// commits invisibly, and `results/README.md` went on asserting the opposite. A convention cannot
/// fail silently if the value is in the log.
///
/// `extra` is for whatever the harness varies that is not in the config — the window, the arm
/// label, the resolution. Errors are returned rather than swallowed: a sidecar that silently
/// failed to write would reproduce the defect exactly.
pub fn provenance_sidecar(
    png_path: &str,
    cfg: &crate::ensemble::pixel::EnsembleCfg,
    extra: &str,
) -> std::io::Result<()> {
    let stem = png_path.strip_suffix(".png").unwrap_or(png_path);
    let body = format!(
        "image={png_path}\nconfig={}\n{}{}",
        cfg.provenance(),
        extra,
        if extra.ends_with('\n') || extra.is_empty() { "" } else { "\n" }
    );
    std::fs::write(format!("{stem}.cfg.txt"), body)
}
