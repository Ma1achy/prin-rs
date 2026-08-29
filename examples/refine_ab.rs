//! **The refinement pass, on and off, on the same window.** `EnsembleCfg::default()` has
//! `refine_flagged: true`; every render harness in `examples/` sets it to `false`, including the
//! one that made the committed panel. This renders both arms so the difference is a picture.
use rayon::prelude::*;
use prin_rs::ensemble::pixel::{EnsembleCfg, PixelOut};
use prin_rs::grid::{self, Chart};
use prin_rs::outcome::{EscapeRule, CLOSURE_TAU};
use prin_rs::output::colour::Scalar;
use prin_rs::output::{adaptive, colour, png};

const Z: [f64; 10] = [-0.098, 0.11, -0.034, 0.067, -0.093, -0.066, 0.027, -0.114, 0.107, -0.116];
const ZOOM: f64 = 0.637_63;
const PAN: (f64, f64) = (0.175_12, 0.177_27);
const WINDOW: f64 = 0.4;
const DLO: f64 = 1e-12;
const DHI: f64 = 1e2;

fn ramp(x: f64) -> [u8; 3] {
    const S: [[f64; 3]; 5] = [[0.0,0.0,0.015],[0.34,0.06,0.43],[0.72,0.21,0.33],[0.98,0.55,0.04],[0.99,1.0,0.64]];
    let t = x.clamp(0.0,1.0)*4.0; let i=(t.floor() as usize).min(3); let f=t-i as f64;
    let mut o=[0u8;3]; for k in 0..3 { o[k]=(255.0*(S[i][k]*(1.0-f)+S[i+1][k]*f)).clamp(0.0,255.0) as u8; } o
}

fn main() {
    let res: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(192);
    let u: f64 = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(0.942);
    let v: f64 = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(0.789);
    let h: f64 = std::env::args().nth(4).and_then(|s| s.parse().ok()).unwrap_or(0.0522);
    let tag: String = std::env::args().nth(5).unwrap_or_else(|| "B10".into());
    let out: String = std::env::args().nth(6).unwrap_or_else(|| "refine_ab".into());
    let _ = std::fs::create_dir_all(&out);

    let z0 = prin_rs::physics::decoder::Latent { z_alpha: Z[1], z_beta: Z[0],
        z_q: [Z[4],Z[5],Z[6],Z[7]], z_mu: [Z[8],Z[9]] };
    let (mut q1, mut q2) = ([0.0f64;8],[0.0f64;8]); q1[1]=1.0; q2[0]=1.0;
    let chart = Chart::Latent { z0, q1, q2 };
    let (cx0,cy0,half0)=(2.0*PAN.0-1.0+ZOOM, 2.0*PAN.1-1.0+ZOOM, ZOOM);
    let cx = cx0 + (2.0*u-1.0)*half0;
    let cy = cy0 + (2.0*v-1.0)*half0;
    let half = 2.0*h*half0;
    let (t_max, r_coll) = (50.0f64, 0.005f64);
    let n_sync = (t_max/WINDOW).round().max(4.0) as usize;
    let m = grid::decode_state(&chart, 0, cx, cy).m;
    let sites = colour::landmarks(&m);
    let sl = grid::Slice::body_plane(res,res,cx,cy,half,0).with_chart(chart);
    println!("{tag} at {res}^2, cell {:.3e}\n", 2.0*half/res as f64);

    for refine in [false, true] {
        let ens = EnsembleCfg { refine_flagged: refine, t_max, n_sync, r_coll_frac: r_coll,
            escape_rule: EscapeRule::Closure(CLOSURE_TAU), closure_k: 1,
            stop_on_escape: false, ..Default::default() };
        let t0 = std::time::Instant::now();
        let px: Vec<PixelOut> = (0..sl.npix()).into_par_iter()
            .map(|k| prin_rs::ensemble::pixel::evaluate::<f64>(&sl,k,&ens)).collect();
        let (lo,hi)=colour::range(&px,Scalar::ShapeSpread);
        let l=(DLO.log10(),DHI.log10());
        let (mut a,mut b,mut c)=(Vec::new(),Vec::new(),Vec::new());
        for p in &px {
            a.extend_from_slice(&colour::rgb(p,Scalar::ShapeSpread,&sites,lo,hi));
            b.extend_from_slice(&png::outcome_rgb(p));
            c.extend_from_slice(&if p.n_nonfinite>0 || !p.energy_drift_max.is_finite() { [255,0,255] }
                else { ramp((p.energy_drift_max.max(1e-300).log10()-l.0)/(l.1-l.0)) });
        }
        let s = if refine {"on"} else {"off"};
        for (k,buf) in [("uniform",&a),("outcome",&b),("drift",&c)] {
            let _ = adaptive::save_rect(&format!("{out}/{tag}_refine{s}_{k}.png"),res,res,buf);
        }
        let mut d: Vec<f64> = px.iter().map(|p| p.energy_drift_max).filter(|x| x.is_finite()).collect();
        d.sort_by(|x,y| x.partial_cmp(y).unwrap());
        let mut e: Vec<f64> = px.iter().map(|p| p.error_ratio).filter(|x| x.is_finite()).collect();
        e.sort_by(|x,y| x.partial_cmp(y).unwrap());
        let mut sp: Vec<f64> = px.iter().map(|p| p.spread_shape).filter(|x| x.is_finite()).collect();
        sp.sort_by(|x,y| x.partial_cmp(y).unwrap());
        println!("  refine {s:>3}  {:.0}s  drift p50 {:.3e}  max {:.3e}   error_ratio p50 {:.4e}   \
                  spread_shape p50 {:.3e}   nonfin {}   esc {:.4} col {:.4} bnd {:.4}   ramp ({lo:.2e},{hi:.2e})",
            t0.elapsed().as_secs_f64(), d[d.len()/2], d[d.len()-1], e[e.len()/2], sp[sp.len()/2],
            px.iter().filter(|p| p.n_nonfinite>0).count(),
            px.iter().filter(|p| p.state==0).count() as f64/px.len() as f64,
            px.iter().filter(|p| p.state==2).count() as f64/px.len() as f64,
            px.iter().filter(|p| p.state==1).count() as f64/px.len() as f64);
    }
}
