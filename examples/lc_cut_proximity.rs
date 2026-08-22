//! Does the cross-check geometry actually visit the LC branch cut?
//!
//! The branch cut is at rho pointing along -x, i.e. u0/|u| -> 0. If the near-field
//! trajectories never go near it, the stable inverse is inert at f64 and the cross-check
//! cannot move. If they do go near it, the two branches register different u and the
//! trajectories separate.
use prin_rs::grid;
use prin_rs::integrate::az::lc;
use prin_rs::integrate::az::{reference_body::choose_reference, reference_body::triple, AzSystem};
use prin_rs::physics::burrau;
use prin_rs::Real;

fn main() {
    let m = burrau::masses::<f64>();
    let s = grid::region("near-field", 3, 3, 0.05).unwrap();

    for &(t_max, n_sync) in [(0.5f64, 1usize), (1.0, 2), (2.0, 5), (13.0, 32)].iter() {
        let mut worst_cos = 1.0f64;   // smallest u0/|u| seen at any registration
        let mut worst_ulp = 0.0f64;   // largest |u_ref - u_stable| / |u|
        let mut n_reg = 0usize;
        let mut n_differ = 0usize;

        for i in 0..s.npix() {
            // Walk the sync boundaries exactly as the driver does, using the reference
            // branch for the trajectory so this measures the geometry the gate actually runs.
            let mut cart = s.nominal::<f64>(i);
            let mut t = 0.0f64;
            for kk in 0..n_sync {
                let t_target = ((kk + 1) as f64) * t_max / (n_sync as f64);
                let a = choose_reference(&cart.r);
                let (ra, rb, rc) = triple(a);
                let sys = AzSystem::new(ra, rb, rc, m).with_reference_lc();

                let r1 = cart.r[rb] - cart.r[ra];
                let r2 = cart.r[rc] - cart.r[ra];
                for rho in [r1, r2] {
                    n_reg += 1;
                    let ur = lc::u_of_rho_reference(rho);
                    let us = lc::u_of_rho(rho);
                    let cos = ur.x.abs() / ur.norm().max(f64::TINY);
                    worst_cos = worst_cos.min(cos);
                    let d = (ur - us).norm() / us.norm().max(f64::TINY);
                    if d > 0.0 {
                        n_differ += 1;
                    }
                    worst_ulp = worst_ulp.max(d);
                }

                let _ = sys;
                let out = prin_rs::integrate::az::integrate_az_lc(
                    cart, &m, t_target - t, 1, 0.01, 30_000, None, false,
                );
                cart = out.state;
                t = t_target;
            }
        }
        println!(
            "t_max={t_max:>5}  n_sync={n_sync:>3}  registrations={n_reg:>5}  \
             min u0/|u| = {worst_cos:.3e}  differ = {n_differ}/{n_reg}  max |du|/|u| = {worst_ulp:.3e}"
        );
    }
    println!();
    println!("u0/|u| near 1 means rho points along +x, far from the cut. Near 0 means it sits");
    println!("on the cut, where the reference branch cancels. 'differ' counts registrations");
    println!("where the two branches produce a different u at all, even by one ulp.");
}
