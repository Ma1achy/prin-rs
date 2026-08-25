//! `prinq` — the minimal refinement scheduler.
//!
//! Starts from one quad, refines adaptively using the criterion at a fixed playhead, dumps the
//! tree with every decision and its reason, exits. **Batch job.**
//!
//! It exists to answer what only appears when the criterion runs **in a loop** — nothing before
//! this ran it more than once. See `SCHEDULER_BRIEF.md`.
//!
//! **No eviction, no caching, no async, no promotion, no camera, no interaction.** Everything cut
//! is deliberate; if one of them appears here, it is a bug.

use std::fs::File;
use std::io::BufWriter;

use prin_rs::ensemble::pixel::EnsembleCfg;
use prin_rs::grid;
use prin_rs::output::{png, tree as treeout};
use prin_rs::quad::Agg;
use prin_rs::render::{self, Precision};
use prin_rs::scheduler::{self, Order, Policy, SchedCfg};
use prin_rs::spatial::HotRule;

const USAGE: &str = "\
prinq — minimal refinement scheduler for prin-rs

  --region <name>       one of BRIEF §2.2's regions (default near-field)
  --half <x>            root quad half-width (default 0.05)
  --samples <n>         N, footprints per quad axis (default 8; must be >= 2)
  --budget <n>          cap on QUADS computed (default 2000). One quad is N^2*(E+1)
                        trajectories: 512 at N=8, E+1=8, about 47 ms
  --bootstrap <n>       levels split unconditionally before any decision (default 2)
  --max-level <n>       optional depth cap; omit to run to the budget (§4 q1 needs it omitted)
  --tau <x>             tau_display, the spread below which a quad is kept (default 1e-2)
  --hot-rule <r>        abs | q:<fraction> (default q:0.5). Which mask the shape criteria
                        read: `abs` cuts at --tau, `q:0.5` at the quad's own median. Both
                        masks are always computed and dumped; this picks the one that feeds
                        `layout_rel` and `grad_rms`. Note n_hot is fixed by the rule under
                        any quantile, so frac_hot is not a signal there
  --alpha-hi <x>        split at or above this exponent (default 0.5)
  --alpha-lo <x>        floor below this exponent (default 0.2)
  --alpha-band <f>      set alpha_lo = f * alpha_hi; give it AFTER --alpha-hi. Without it,
                        raising alpha_hi alone leaves a zero-width keep band
  --sib-tau <x>         floor above this sibling range, sibling policy (default 0.5)
  --policy <p>          alpha | sibling (default alpha)
  --order <o>           spread | spread-area | shuffled (default spread)
  --agg <a>             mean | median | p90 (default median)
  --copies <n>          E+1 copies per footprint (default 8)
  --t-max <t>           playhead (default 13)
  --eta <x>             timestep parameter (default 0.01)
  --precision <p>       f32 | f64 (default f64)
  --overlay <n>         also render the region uniformly at n x n and draw leaf boundaries
  --seed <n>            shuffle seed
  --out <stem>          output stem (default tree)
";

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.iter().any(|a| a == "--help" || a == "-h") {
        print!("{USAGE}");
        return;
    }

    let mut region = "near-field".to_string();
    let mut half = 0.05f64;
    let mut cfg = SchedCfg::default();
    let mut ens = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let mut precision = Precision::F64;
    let mut overlay: Option<usize> = None;
    let mut out = "tree".to_string();

    let get = |i: usize| -> String {
        args.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value", args[i]))
    };
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--region" => { region = get(i); i += 1; }
            "--half" => { half = get(i).parse().unwrap(); i += 1; }
            "--samples" => { cfg.n = get(i).parse().unwrap(); i += 1; }
            "--budget" => { cfg.budget = get(i).parse().unwrap(); i += 1; }
            "--bootstrap" => { cfg.bootstrap_levels = get(i).parse().unwrap(); i += 1; }
            "--max-level" => { cfg.max_level = Some(get(i).parse().unwrap()); i += 1; }
            "--tau" => { cfg.tau_display = get(i).parse().unwrap(); i += 1; }
            "--hot-rule" => {
                // `abs` uses --tau; `q:<f>` uses the quad's own quantile. Both masks are always
                // computed and dumped -- this picks which one the mask-derived CRITERIA read.
                let v = get(i);
                cfg.hot_rule = match v.split_once(':') {
                    Some(("q", f)) => HotRule::Quantile(f.parse().expect("quantile")),
                    _ if v == "abs" => HotRule::AbsTau(cfg.tau_display),
                    _ => panic!("--hot-rule takes `abs` or `q:<fraction>`, got {v}"),
                };
                i += 1;
            }
            "--alpha-hi" => { cfg.alpha_hi = get(i).parse().unwrap(); i += 1; }
            "--alpha-band" => {
                // Set alpha_lo as a fraction of alpha_hi, so the keep band is not silently
                // zero-width when only --alpha-hi is given.
                let f: f64 = get(i).parse().unwrap();
                cfg.alpha_lo = cfg.alpha_hi * f;
                i += 1;
            }
            "--alpha-lo" => { cfg.alpha_lo = get(i).parse().unwrap(); i += 1; }
            "--sib-tau" => { cfg.sib_tau = get(i).parse().unwrap(); i += 1; }
            "--policy" => { cfg.policy = Policy::parse(&get(i)).expect("policy"); i += 1; }
            "--order" => { cfg.order = Order::parse(&get(i)).expect("order"); i += 1; }
            "--agg" => { cfg.agg = Agg::parse(&get(i)).expect("agg"); i += 1; }
            "--copies" => { ens.n_extra = get(i).parse::<usize>().unwrap() - 1; i += 1; }
            "--t-max" => { ens.t_max = get(i).parse().unwrap(); i += 1; }
            "--eta" => { ens.eta = get(i).parse().unwrap(); i += 1; }
            "--precision" => {
                precision = match get(i).as_str() {
                    "f32" => Precision::F32,
                    "f64" => Precision::F64,
                    o => panic!("unknown precision {o}"),
                };
                i += 1;
            }
            "--overlay" => { overlay = Some(get(i).parse().unwrap()); i += 1; }
            "--seed" => { cfg.seed = get(i).parse().unwrap(); i += 1; }
            "--out" => { out = get(i); i += 1; }
            other => panic!("unknown flag {other}\n\n{USAGE}"),
        }
        i += 1;
    }

    let root = grid::region(&region, 2, 2, half).unwrap_or_else(|| panic!("unknown region {region}"));
    let (tree, st) = scheduler::descend(root.cx, root.cy, half, root.body, &cfg, &ens, precision);

    let hist = tree.depth_histogram();
    let leaves: Vec<usize> = tree.leaves().collect();
    let count = |d: prin_rs::quad::Decision| {
        leaves.iter().filter(|&&i| tree.nodes[i].decision == d).count()
    };
    use prin_rs::quad::Decision as D;

    println!("region {region} half={half} body={}  N={} copies={}  precision={}",
             root.body, cfg.n, ens.n_extra + 1, precision.name());
    println!("policy={} order={} agg={} tau={} alpha_hi={} alpha_lo={} sib_tau={} bootstrap={}",
             cfg.policy.name(), cfg.order.name(), cfg.agg.name(),
             cfg.tau_display, cfg.alpha_hi, cfg.alpha_lo, cfg.sib_tau, cfg.bootstrap_levels);
    println!();
    println!("  {} quads computed of a {} budget, {} footprints, {} trajectories",
             st.quads_computed, cfg.budget, st.footprints,
             st.footprints * (ens.n_extra + 1));
    println!("  {:.2} s wall, {:.1} ms/quad", st.wall_seconds,
             1e3 * st.wall_seconds / st.quads_computed.max(1) as f64);
    println!("  {} iterations, {} leaves", st.iterations, leaves.len());
    println!();
    println!("  leaf decisions:  split {}  floor {}  keep {}  precision_floor {}  max_level {}  budget {}",
             count(D::Split), count(D::Floor), count(D::Keep),
             count(D::PrecisionFloor), count(D::MaxLevel), count(D::BudgetExhausted));
    println!("  budget exhausted: {}", st.budget_exhausted);
    println!();
    print!("  depth histogram:");
    for (l, n) in hist.iter().enumerate() {
        if *n > 0 {
            print!(" {l}:{n}");
        }
    }
    println!();
    print!("  leaves per iteration:");
    for n in &st.leaves_per_iteration {
        print!(" {n}");
    }
    println!();

    let dump = format!("{out}.tree");
    let mut f = BufWriter::new(File::create(&dump).expect("create dump"));
    treeout::write(&mut f, &tree, &cfg, &ens, &st, &region, precision.name()).expect("write dump");
    println!();
    println!("wrote {dump}");

    if let Some(res) = overlay {
        let slice = grid::Slice::body_plane(res, res, root.cx, root.cy, half, root.body).with_chart(tree.chart);
        let base = render::render(&slice, &ens, precision);
        treeout::overlay(&out, "tree_outcome", &tree, &base, res, png::outcome_rgb)
            .expect("overlay");

        // The spread base is the direct check: the tree tracks ensemble_spread, not outcome
        // labels, and near-field's outcome image is 97.7% one colour. Log scaled over the
        // grid's own p1..p99, as the uniform kernel's spread image is.
        let mut fin: Vec<f64> = base
            .iter()
            .map(|p| p.ensemble_spread)
            .filter(|x| x.is_finite() && *x > 0.0)
            .collect();
        fin.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let q = |f: f64| fin[(((fin.len() - 1) as f64) * f).round() as usize];
        let (lo, hi) = (q(0.01).max(1e-12), q(0.99));
        let (ll, lh) = (lo.ln(), hi.ln());
        treeout::overlay(&out, "tree_spread", &tree, &base, res, move |p| {
            let v = p.ensemble_spread;
            let t = if !v.is_finite() {
                1.0
            } else if v <= 0.0 {
                0.0
            } else {
                ((v.ln() - ll) / (lh - ll)).clamp(0.0, 1.0)
            };
            [(255.0 * t.powf(0.6)) as u8,
             (255.0 * (t * (1.0 - t) * 4.0).powf(0.8)) as u8,
             (255.0 * (1.0 - t).powf(0.6)) as u8]
        })
        .expect("overlay");
        println!("wrote {out}_tree_outcome.png and {out}_tree_spread.png ({res}x{res} base)");
        println!("spread base log scaled over [{lo:.3e}, {hi:.3e}]; the tree tracks spread, not");
        println!("outcome labels, so that is the base the tree should be checked against.");
    }

    println!();
    println!("The tree and the decisions are the output; the image is a diagnostic. A threshold");
    println!("chosen because the picture looked right is an arbitrary constant — sweep, do not tune.");
}
