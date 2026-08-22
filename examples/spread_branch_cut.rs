//! Does `spread_shape` inherit the LC branch-cut sensitivity?
//!
//! The concern is specific: copies of one pixel differ in configuration, so they can straddle
//! the cut differently. A spread *across copies* would then partly measure registration error
//! rather than dynamics — which is the exact shape of "single-trajectory drift looks fine but
//! the ensemble diagnostic breaks".
//!
//! Measured rather than assumed, given how the last "this should be inert" went.
use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::grid;

fn main() {
    println!("spread_shape with the reference vs the conditioned inverse LC branch");
    println!("(near-field is the region whose pair (1,2) separation sits exactly on the cut)\n");
    println!(
        "{:>14}{:>8}{:>14}{:>14}{:>12}{:>12}",
        "region", "prec", "unstable", "stable", "max |d|", "max rel"
    );

    for region in ["near-field", "body2 core", "body1 slice", "mid-field"] {
        let s = grid::region(region, 5, 5, 0.05).unwrap();
        for prec in ["f64", "f32"] {
            let mut worst_abs = 0.0f64;
            let mut worst_rel = 0.0f64;
            let mut mean_u = 0.0f64;
            let mut mean_s = 0.0f64;
            for i in 0..s.npix() {
                let cu = EnsembleCfg { t_max: 13.0, lc_stable: false, ..Default::default() };
                let cs = EnsembleCfg { t_max: 13.0, lc_stable: true, ..Default::default() };
                let (u, st) = if prec == "f64" {
                    (evaluate::<f64>(&s, i, &cu), evaluate::<f64>(&s, i, &cs))
                } else {
                    (evaluate::<f32>(&s, i, &cu), evaluate::<f32>(&s, i, &cs))
                };
                mean_u += u.spread_shape;
                mean_s += st.spread_shape;
                let d = (u.spread_shape - st.spread_shape).abs();
                worst_abs = worst_abs.max(d);
                worst_rel = worst_rel.max(d / st.spread_shape.max(1e-300));
            }
            let n = s.npix() as f64;
            println!(
                "{region:>14}{prec:>8}{:>14.4e}{:>14.4e}{worst_abs:>12.2e}{worst_rel:>12.2e}",
                mean_u / n,
                mean_s / n
            );
        }
    }
    println!();
    println!("Columns 3 and 4 are the grid-mean spread_shape under each branch. If the");
    println!("conditioning defect leaks into the ensemble diagnostic, the two differ by more");
    println!("than the dynamics justify - and it should show at f32 first.");
}
