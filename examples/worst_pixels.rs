//! Where does the damage sit, and does error_ratio flag the same pixels drift does?
use prin_rs::ensemble::pixel::{evaluate, EnsembleCfg};
use prin_rs::grid;

fn main() {
    let s = grid::region("near-field", 32, 32, 0.05).unwrap();
    let cfg = EnsembleCfg { refine_flagged: false, ..Default::default() };
    let px: Vec<_> = (0..s.npix()).map(|i| evaluate::<f64>(&s, i, &cfg)).collect();

    let mut by_drift: Vec<usize> = (0..px.len()).collect();
    by_drift.sort_by(|&a, &b| px[b].energy_drift_max.partial_cmp(&px[a].energy_drift_max).unwrap());

    println!("worst 8 pixels by max |dE/E| over the ensemble:");
    println!("{:>7}{:>13}{:>13}{:>15}{:>13}", "pixel", "drift_max", "er (MAD)", "er (max dev)", "d_min_true");
    for &i in by_drift.iter().take(8) {
        let p = &px[i];
        println!("{i:>7}{:>13.3e}{:>13.4}{:>15.4e}{:>13.3e}",
                 p.energy_drift_max, p.error_ratio_mad, p.error_ratio, p.d_min_true);
    }

    let mut d: Vec<f64> = px.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |f: f64| d[((d.len() - 1) as f64 * f).round() as usize];
    println!("\ndrift_max distribution: median {:.3e}  p90 {:.3e}  p99 {:.3e}  max {:.3e}",
             q(0.5), q(0.9), q(0.99), q(1.0));
    println!("pixels with drift_max > 1e-3: {}", d.iter().filter(|x| **x > 1e-3).count());

    // Does error_ratio flag the damaged pixels? Rank correlation between the two.
    let mut rd: Vec<usize> = (0..px.len()).collect();
    rd.sort_by(|&a, &b| px[a].energy_drift_max.partial_cmp(&px[b].energy_drift_max).unwrap());
    let mut re: Vec<usize> = (0..px.len()).collect();
    re.sort_by(|&a, &b| px[a].error_ratio_mad.partial_cmp(&px[b].error_ratio_mad).unwrap());
    let mut rr: Vec<usize> = (0..px.len()).collect();
    rr.sort_by(|&a, &b| px[a].error_ratio.partial_cmp(&px[b].error_ratio).unwrap());
    let mut rank_d = vec![0.0f64; px.len()];
    let mut rank_e = vec![0.0f64; px.len()];
    for (r, &i) in rd.iter().enumerate() { rank_d[i] = r as f64; }
    for (r, &i) in re.iter().enumerate() { rank_e[i] = r as f64; }
    let n = px.len() as f64;
    let mean = (n - 1.0) / 2.0;
    let (mut num, mut da, mut db) = (0.0, 0.0, 0.0);
    for i in 0..px.len() {
        let (x, y) = (rank_d[i] - mean, rank_e[i] - mean);
        num += x * y; da += x * x; db += y * y;
    }
    let mad_rho = num / (da.sqrt() * db.sqrt());

    let mut rank_r = vec![0.0f64; px.len()];
    for (r, &i) in rr.iter().enumerate() { rank_r[i] = r as f64; }
    let (mut num2, mut da2, mut db2) = (0.0, 0.0, 0.0);
    for i in 0..px.len() {
        let (x, y) = (rank_d[i] - mean, rank_r[i] - mean);
        num2 += x * y; da2 += x * x; db2 += y * y;
    }
    let range_rho = num2 / (da2.sqrt() * db2.sqrt());

    println!("\nSpearman against drift_max over ALL pixels:");
    println!("  error_ratio built on MAD        {mad_rho:+.4}");
    println!("  error_ratio built on max dev    {range_rho:+.4}");
    println!();
    println!("Both near zero - but Spearman over every pixel is the wrong summary for a");
    println!("statistic BRIEF §4 itself calls a boolean flag. 1024 healthy pixels sit at");
    println!("ratio ~= 1 plus noise, and their ranks swamp the 23 damaged ones. Detection is");
    println!("the question that matters:");
    println!();

    let damaged: Vec<usize> = (0..px.len()).filter(|&i| px[i].energy_drift_max > 1e-3).collect();
    let healthy: Vec<usize> = (0..px.len()).filter(|&i| px[i].energy_drift_max <= 1e-3).collect();
    println!("  {} damaged pixels (drift_max > 1e-3), {} healthy", damaged.len(), healthy.len());
    println!();
    println!("{:>16}{:>14}{:>14}{:>16}", "statistic", "damaged med", "healthy p99", "separation");
    for (name, f) in [
        ("MAD", (|p: &prin_rs::ensemble::pixel::PixelOut| p.error_ratio_mad) as fn(&_) -> f64),
        ("max deviation", |p: &prin_rs::ensemble::pixel::PixelOut| p.error_ratio),
    ] {
        let mut dv: Vec<f64> = damaged.iter().map(|&i| f(&px[i])).filter(|x| x.is_finite()).collect();
        let mut hv: Vec<f64> = healthy.iter().map(|&i| f(&px[i])).filter(|x| x.is_finite()).collect();
        dv.sort_by(|a, b| a.partial_cmp(b).unwrap());
        hv.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let dmed = dv[dv.len() / 2];
        let hp99 = hv[((hv.len() - 1) as f64 * 0.99).round() as usize];
        println!("{name:>16}{dmed:>14.4e}{hp99:>14.4e}{:>16.2}", dmed / hp99);
    }
    println!();
    println!("Separation is the damaged median over the healthy p99: how far a threshold");
    println!("could be placed between them. Below ~1 the flag cannot separate the two");
    println!("populations at all.");
}
