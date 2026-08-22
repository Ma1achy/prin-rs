//! Does the escape arm ever fire on the near-field slice, and does it agree with the
//! reference's own event detection?
use rayon::prelude::*;
use prin_rs::integrate::az::{self, AzOpts};
use prin_rs::physics::burrau;
use prin_rs::{grid, outcome};

fn main() {
    for t_max in [13.0f64, 16.0, 20.0, 30.0] {
        run(t_max);
    }
}

fn run(t_max: f64) {
    let s = grid::region("near-field", 32, 32, 0.05).unwrap();
    let m = burrau::masses::<f64>();
    let out: Vec<(u8, bool, f64)> = (0..s.npix())
        .into_par_iter()
        .map(|i| {
            let o = az::integrate_az_opts(
                s.nominal::<f64>(i), &m, t_max, 32, 0.01, 30_000,
                &AzOpts { stop_on_event: false, ..Default::default() },
            );
            let legacy = outcome::classify_legacy(&o.state, &m);
            (legacy, o.events.escape.is_some(), o.events.escape.map(|(_, t)| t).unwrap_or(f64::NAN))
        })
        .collect();
    let mut byclass = [0usize; 4];
    for (c, _, _) in &out {
        byclass[*c as usize] += 1;
    }
    println!("t_max = {t_max}: classify_legacy {:?} (0,1,2 = escaping body, 3 = bound)", byclass);
    println!("  escape arm fired at some sync boundary on {} of {} pixels",
             out.iter().filter(|(_, e, _)| *e).count(), out.len());
    let mut ts: Vec<f64> = out.iter().filter_map(|(_, e, t)| e.then_some(*t)).collect();
    ts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if !ts.is_empty() {
        println!("  escape times: min {:.3} median {:.3} max {:.3}",
                 ts[0], ts[ts.len()/2], ts[ts.len()-1]);
    }
}
