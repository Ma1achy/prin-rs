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

/// A line sink that writes to stdout **and** to a log file, so a run's table is committed
/// beside the panels it describes.
///
/// `results/output/chart_gallery.txt` was the 25 August run while the panels beside it were
/// from 3 September: the regeneration committed pictures and sidecars but not its stdout,
/// and the txt kept describing a corpus that no longer existed. A harness that tees its own
/// log cannot leave the log behind. Interior mutability so the handle can be shared into the
/// closures a harness builds its rows in.
pub struct Log {
    file: std::cell::RefCell<Option<std::io::BufWriter<std::fs::File>>>,
}

impl Log {
    /// Tee to `path`, creating its parent. A path that cannot be opened logs to stdout only,
    /// and says so once, rather than silently writing nothing.
    pub fn tee(path: &str) -> Log {
        if let Some(dir) = std::path::Path::new(path).parent() {
            let _ = std::fs::create_dir_all(dir);
        }
        let file = match std::fs::File::create(path) {
            Ok(f) => Some(std::io::BufWriter::new(f)),
            Err(e) => {
                println!("  (log file `{path}` could not be opened: {e}; stdout only)");
                None
            }
        };
        Log { file: std::cell::RefCell::new(file) }
    }

    /// stdout only.
    pub fn stdout() -> Log {
        Log { file: std::cell::RefCell::new(None) }
    }

    pub fn line(&self, s: &str) {
        use std::io::Write;
        println!("{s}");
        if let Some(f) = self.file.borrow_mut().as_mut() {
            let _ = writeln!(f, "{s}");
            let _ = f.flush();
        }
    }
}

/// `println!` that also lands in a [`Log`]: `logln!(log, "...", args)`.
#[macro_export]
macro_rules! logln {
    ($log:expr) => {
        $log.line("")
    };
    ($log:expr, $($arg:tt)*) => {
        $log.line(&format!($($arg)*))
    };
}
