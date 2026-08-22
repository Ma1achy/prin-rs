//! Render a uniform-resolution slice: one pass, PNG plus raw dump, exit.
//!
//! No quadtree, no scheduler, no GUI, no streaming, no interaction. Each omission is
//! deliberate. If this grows a scheduler, that is a bug.

use std::fs::File;
use std::io::BufWriter;
use std::time::Instant;

use prin_rs::{config, output, render};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cfg = match config::parse(&args) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2);
        }
    };

    let t0 = Instant::now();
    let pixels = render::render(&cfg.slice, &cfg.ens, cfg.precision);
    let elapsed = t0.elapsed().as_secs_f64();
    let s = render::summarise(&pixels);

    let dump = format!("{}.raw", cfg.out);
    let mut w = BufWriter::new(File::create(&dump).expect("cannot create raw dump"));
    output::raw::write(&mut w, &cfg.slice, &cfg.ens, cfg.precision.name(), &pixels)
        .expect("raw dump failed");
    drop(w);
    output::png::write_pair(&cfg.out, &cfg.slice, &pixels).expect("png failed");

    println!("region {} {}x{} half={} body={}  precision={}  t_max={} eta={} copies={}",
             cfg.region, cfg.slice.nx, cfg.slice.ny, cfg.slice.half, cfg.slice.body,
             cfg.precision.name(), cfg.ens.t_max, cfg.ens.eta, cfg.ens.n_extra + 1);
    println!("  {:.2} s wall, {:.2} ms/pixel", elapsed, 1e3 * elapsed / s.n as f64);
    println!();
    println!("  error_ratio   max {:.4e}  p99 {:.4e}  median {:.6}  argmax pixel {}",
             s.error_ratio_max, s.error_ratio_p99, s.error_ratio_median, s.error_ratio_argmax);
    println!("  sigma_E(0)    median {:.6e}   <- proportional to cell width; see below",
             s.sigma_e_0_median);
    println!("  |dE/E|        median {:.3e}  max {:.3e}", s.drift_median, s.drift_max);
    println!("  d_min_gap     median {:.3e}  max {:.3e}", s.d_min_gap_median, s.d_min_gap_max);
    println!("  ref_disagree  total {}", s.ref_disagree_total);
    println!("  t_end         median {:.4}", s.t_end_median);
    print!("  outcomes     ");
    for k in 0..6u8 {
        if s.state_fracs[k as usize] > 0.0 {
            print!(" {} {:.3}",
                   prin_rs::outcome::State::from_bits(k).unwrap().name(),
                   s.state_fracs[k as usize]);
        }
    }
    println!();
    println!("  pixels with a non-finite copy: {} of {}", s.n_pixels_with_nonfinite, s.n);
    println!();
    println!("wrote {dump}, {}_outcome.png, {}_spread.png", cfg.out, cfg.out);
    println!();
    println!("error_ratio uses the maximum deviation within a footprint and aggregates across");
    println!("pixels by max. It is a boolean flag: its magnitude is unstable.");
    println!("It also inflates with resolution for a trivial reason - sigma_E(0) is");
    println!("proportional to the jitter and so to cell width, while integration error is");
    println!("not. Use sigma_e_0 and sigma_e_t from the dump before comparing across sizes.");
}
