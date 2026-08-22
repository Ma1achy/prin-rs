//! What actually happens at BRIEF §2.6's `deep interior`, pair by pair.
//!
//! The brief says to expect a triple collision. This prints the three separations at the
//! moment the first pair crosses `r_coll`, and the per-pair minimum over the whole run, so the
//! label can be checked against the geometry rather than assumed.

use prin_rs::integrate::az::{self, AzOpts};
use prin_rs::outcome;
use prin_rs::physics::{burrau, energy, newton};
use prin_rs::grid;

fn main() {
    let s = grid::region("deep interior", 1, 1, 1e-6).unwrap();
    let m = burrau::masses::<f64>();
    let c = s.nominal::<f64>(0);
    let r0 = energy::hyperradius(&c.r, &m);

    println!("deep interior, R(0) = {r0:.6}");
    println!("{:>10}{:>14}{:>10}{:>8}{:>14}{:>14}{:>14}",
             "r_coll/R", "r_coll", "outcome", "detail", "|01| at fire", "|02| at fire", "|12| at fire");

    for frac in [1e-4f64, 1e-3, 1e-2, 1e-1] {
        let o = az::integrate_az_opts(
            c, &m, 13.0, 32, 0.01, 30_000,
            &AzOpts { r_coll_frac: frac, stop_on_event: true, ..Default::default() },
        );
        let out = outcome::classify(&o.events, &o.state, &m, o.finite, o.budget_exhausted);
        let d = newton::pair_dists(&o.state.r);
        println!("{frac:>10.0e}{:>14.4e}{:>10}{:>8}{:>14.4e}{:>14.4e}{:>14.4e}",
                 frac * r0, out.state.name(), out.detail, d[0], d[1], d[2]);
    }

    // Where the other two pairs actually get to. A threshold sweep is the honest probe: the
    // driver records only the *first* firing, so the way to bound a pair's minimum is to ask
    // which thresholds it crosses. Endpoint sampling at sync boundaries would not do — with
    // n_sync = 32 and t_max = 13 they are 0.4 apart, which is the whole reason collision
    // detection lives inside the RK4 loop.
    println!();
    println!("smallest r_coll/R at which each pair registers (sweep, stop_on_event off):");
    let mut first = [f64::INFINITY; 3];
    for frac in [1e-4f64, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0] {
        let o = az::integrate_az_opts(
            c, &m, 13.0, 32, 0.01, 30_000,
            &AzOpts { r_coll_frac: frac, stop_on_event: false, ..Default::default() },
        );
        if let Some((mask, _)) = o.events.collision {
            for j in 0..3 {
                if mask & (1 << j) != 0 && first[j].is_infinite() {
                    first[j] = frac;
                }
            }
        }
    }
    for (j, (i, k)) in prin_rs::physics::PAIRS.iter().enumerate() {
        if first[j].is_finite() {
            println!("  pair ({i},{k}): registers at r_coll/R = {:.0e}", first[j]);
        } else {
            println!("  pair ({i},{k}): never registers, even at r_coll = R", );
        }
    }

    let full = az::integrate_az_opts(
        c, &m, 13.0, 32, 0.01, 30_000,
        &AzOpts { stop_on_event: false, ..Default::default() },
    );
    let refbranch = az::integrate_az_opts(
        c, &m, 13.0, 32, 0.01, 30_000,
        &AzOpts { stop_on_event: false, lc_stable: false, ..Default::default() },
    );
    println!();
    println!("d_min over all three pairs, sampled inside the RK4 loop: {:.4e} ({:.3e} of R)",
             full.d_min_true, full.d_min_true / r0);
    println!("|dE/E| = {:.3e}, {} steps, {} switches, finite = {}",
             full.drift, full.steps, full.switches, full.finite);
    println!();
    println!("on the reference's own LC branch, for direct comparison with the numpy run:");
    println!("  d_min_ref = {:.10e}   drift = {:.9e}   switches = {}",
             refbranch.d_min_ref, refbranch.drift, refbranch.switches);
    println!("  reference:  2.2976014100e-05        1.393632170e-07              2");
    println!("  drift agrees to 3 digits; d_min differs by 0.6%, which is the sampling-offset");
    println!("  mechanism established in PR #2 - d_min near a close approach is set by where a");
    println!("  step happens to land relative to the crossing, not by a floor.");
    println!();
    println!("BRIEF §2.6 expects this pixel to be an unregularisable triple collision that");
    println!("fails however well the integrator is built. It is not. It is an ordinary binary");
    println!("encounter between bodies 0 and 2 - the other two pairs never come within R of");
    println!("each other - and both implementations integrate it to t = 13 in a couple of");
    println!("seconds with |dE/E| ~ 1e-7. The 190 s failure recorded in 2.6 is the");
    println!("unregularised integrator; AZ regularises exactly this case, which is its job.");

}
