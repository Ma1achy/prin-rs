//! **Does the priority bucketing earn its place? Band occupancy, from the committed caches.**
//!
//! `src/frontier.rs` exists because *"camera relevance changes for every quad on every frame the
//! camera moves, so during a gesture the naive version is genuinely recomputing all of it"*
//! (§4.6). But `Frontier::top_k` flattens every bucket and sorts the lot — the per-frame
//! `O(n log n)` the structure was built to remove is **still paid**. `top_k_bounded` walks bands
//! from the top and stops once no lower band can contend, which is exact; whether it *bites* is an
//! empirical question about how the signal distributes across bands.
//!
//! **The failure to look for is saturation**, which this project has already met twice — at the
//! hot mask (`n_hot = N^2` in 99-100% of quads) and at the linear ramp. If the stored priorities
//! pile into two or three bands, the walk cannot stop early, the scan fraction goes to 1, and the
//! frontier is a `HashMap` with extra steps. That would make it the third mechanism here to
//! compute, sort and change nothing — after `k_frac = 1.0` and the disabled repair pass.
//!
//! **Zero trajectories.** The stored term at production settings is `signal(Within, Median)` =
//! `spread_median`, which every `.qcache` carries per quad. So this reads committed files.
//!
//! The bound to read it against: at `k = f*n` the walk must collect `f*n` entries before it can
//! stop, so **`scan/n >= k/n` by construction**. The column that matters is the *overshoot*,
//! `scan / k` — 1.0 is perfect, and `n/k` is the full sort.
//!
//! Run: `cargo run --release --example frontier_bands -- [dir=results/payload] [charts]`

use prin_rs::camera::Camera;
use prin_rs::frontier::{band_of, Frontier, BANDS};
use prin_rs::output::qcache;

fn tok(header: &str, name: &str) -> Option<f64> {
    header
        .split_whitespace()
        .find_map(|t| t.strip_prefix(&format!("{name}=")))
        .and_then(|v| v.parse().ok())
}

fn main() {
    let dir: String = std::env::args().nth(1).unwrap_or_else(|| "results/payload".into());
    let charts: Vec<String> = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "near-field,deep_interior,far,preset_shape_h1".into())
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();

    println!("frontier_bands: {dir}, {BANDS} log bands over [1e-12, 1e2]");
    println!("  stored term = spread_median (signal(Within, Median), the production criterion)");
    println!("  `overshoot` = scanned / k: 1.00 is perfect, n/k is the full sort `top_k` pays.\n");

    for name in &charts {
        let path = format!("{dir}/{name}_t13_L6.qcache");
        let Ok(file) = std::fs::File::open(&path) else {
            println!("{name:>18}  no cache at {path}, skipped");
            continue;
        };
        let qr = match qcache::read(&mut std::io::BufReader::new(file)) {
            Ok(q) => q,
            Err(e) => {
                println!("{name:>18}  {e}");
                continue;
            }
        };

        // Every quad of the complete tree, as the frontier would hold them.
        let mut f = Frontier::new();
        let mut hist = vec![0usize; BANDS];
        for i in 0..qr.rows.len() {
            let stored = qr.get(i, "spread_median");
            f.insert(i, stored);
            hist[band_of(stored)] += 1;
        }
        let n = f.len();
        let occupied = hist.iter().filter(|&&c| c > 0).count();
        let modal = *hist.iter().max().unwrap_or(&0);

        println!("{name:>18}  n = {n}, bands occupied {occupied} of {BANDS}, modal share {:.1}%",
                 100.0 * modal as f64 / n.max(1) as f64);
        let occ: Vec<String> = hist
            .iter()
            .enumerate()
            .filter(|(_, &c)| c > 0)
            .map(|(b, c)| format!("{b}:{c}"))
            .collect();
        println!("{:>18}  {}", "", occ.join(" "));

        // Two derive arms. The uniform one is the best case for the early-out; the camera one is
        // what a real frame pays, and it can only demote, so it can only make the walk longer.
        let (cx, cy, half) = (
            tok(&qr.header, "cx").unwrap_or(0.0),
            tok(&qr.header, "cy").unwrap_or(0.0),
            tok(&qr.header, "half").unwrap_or(0.05),
        );
        let res = tok(&qr.header, "res").unwrap_or(512.0) as usize;
        // **The camera arm has to be zoomed in or it is an identity.** `Camera::framing` sets
        // `half_world` to the ROOT half-width, so every quad lies wholly inside the viewport and
        // `relevance` returns 1.0 on all of them — the first cut of this measured the uniform arm
        // twice and reported it as a camera result, identical to three decimals in every cell.
        // A quarter-width viewport offset toward a corner makes the term vary, which is what a
        // real frame pays.
        let cam = Camera::framing(cx + half * 0.4, cy + half * 0.4, half * 0.25, res);
        let boxes: Vec<(f64, f64, f64)> =
            (0..qr.rows.len()).map(|i| (qr.get(i, "cx"), qr.get(i, "cy"), qr.get(i, "half"))).collect();

        // State the arm is live before reading it.
        let rels: Vec<f64> = boxes.iter().map(|&(qx, qy, qh)| cam.relevance(qx, qy, qh, 0.0)).collect();
        let rel_spread = rels.iter().cloned().fold(f64::MIN, f64::max)
            - rels.iter().cloned().fold(f64::MAX, f64::min);
        let off = rels.iter().filter(|&&r| r == 0.0).count();
        println!("{:>18}  camera arm: relevance spans {rel_spread:.3}, {off} of {n} quads fully off-screen{}",
                 "", if rel_spread <= 0.0 { "  <- DEAD ARM, this column is the uniform one" } else { "" });
        // `assert!`, not `debug_assert!`: this only ever runs in release, so a debug assertion
        // here is a guard that cannot fire — written three times in one session before it stuck.
        assert!(rel_spread > 0.0, "the camera arm is dead: relevance is constant, so it is the uniform arm");

        println!("{:>18}  {:>8} {:>9} {:>9} {:>10} {:>9} {:>10}",
                 "", "k/n", "scan/n u", "overshoot", "scan/n vis", "overshoot", "full sort");
        // **The frontier a real frame holds is the VISIBLE quads, not the whole tree.** With the
        // off-screen ones in, 92% of entries carry `relevance = 0`, their priority is 0, the k-th
        // score never leaves band 0 and the walk can never stop — `scan/n = 1.000` everywhere.
        // That is a true statement about the wrong population: the brief's priority weights have
        // visibility dominating and *never compute off-screen*, so those quads are not ranked.
        let mut vis = Frontier::new();
        let mut vis_n = 0usize;
        for i in 0..qr.rows.len() {
            let (qx, qy, qh) = boxes[i];
            if cam.relevance(qx, qy, qh, 0.0) > 0.0 {
                vis.insert(i, qr.get(i, "spread_median"));
                vis_n += 1;
            }
        }
        println!("{:>18}  visible frontier: {vis_n} of {n} quads", "");

        for frac in [0.01, 0.05, 0.25, 0.5] {
            let k = ((n as f64 * frac).ceil() as usize).max(1);
            let kv = ((vis_n as f64 * frac).ceil() as usize).max(1);
            let (_, su) = f.top_k_bounded(k, |_| 1.0);
            let (_, sc) = vis.top_k_bounded(kv, |i| {
                let (qx, qy, qh) = boxes[i];
                cam.relevance(qx, qy, qh, 0.0)
            });
            println!("{:>18}  {frac:>8.2} {:>9.3} {:>9.2}x {:>10.3} {:>9.2}x {:>10.2}x",
                     "", su as f64 / n as f64, su as f64 / k as f64,
                     sc as f64 / vis_n.max(1) as f64, sc as f64 / kv as f64, n as f64 / k as f64);
        }
        println!();
    }
    println!("If `scan/n` sits near 1.00 the signal is band-saturated and the buckets buy nothing;");
    println!("`top_k` would then be the honest implementation and the structure should be dropped.");
}
