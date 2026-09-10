//! Stage 1 measurement: what THIS machine's backends actually report.
//!
//! The determinism note's standing caveat ("the spike's evidence is Metal-only") exists because
//! the backend actually exercised was never printed. So this enumerates every adapter wgpu can
//! see and prints the numbers the fit question turns on -- MEASURED. Everything it cannot see
//! (WebGPU's guaranteed floor, Vulkan's required minimums, D3D12) is cited in the write-up and
//! labelled there as cited, never mixed into this output.

const N: u32 = 8; // src/scheduler.rs:419
const COPIES: u32 = 8; // src/ensemble/pixel.rs:204 (n_extra = 7, "the pixel always carries E + 1")
const SIMSTATE_B: u32 = 136; // memory-tiers: hot SimState, FTLE-on

fn main() {
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
    let adapters = pollster::block_on(instance.enumerate_adapters(wgpu::Backends::all()));

    println!("wgpu 30.0.1 -- adapters visible: {}", adapters.len());
    if adapters.is_empty() {
        println!("NONE reachable. That is itself the finding.");
        return;
    }

    for a in &adapters {
        let i = a.get_info();
        let l = a.limits();
        println!("\n=== {:?} / {} ({:?}) ===", i.backend, i.name, i.device_type);
        println!("  driver : {} {}", i.driver, i.driver_info);

        println!("  max_compute_invocations_per_workgroup : {}", l.max_compute_invocations_per_workgroup);
        println!("  max_compute_workgroup_size_x/y/z      : {} / {} / {}",
            l.max_compute_workgroup_size_x, l.max_compute_workgroup_size_y, l.max_compute_workgroup_size_z);
        println!("  max_compute_workgroup_storage_size    : {} B", l.max_compute_workgroup_storage_size);
        println!("  max_compute_workgroups_per_dimension  : {}", l.max_compute_workgroups_per_dimension);
        println!("  max_storage_buffers_per_shader_stage  : {}", l.max_storage_buffers_per_shader_stage);
        println!("  max_storage_buffer_binding_size       : {} B", l.max_storage_buffer_binding_size);
        println!("  max_buffer_size                       : {} B", l.max_buffer_size);
        println!("  max_bind_groups                       : {}", l.max_bind_groups);
        println!("  max_uniform_buffer_binding_size       : {} B", l.max_uniform_buffer_binding_size);
        // wgpu 30 renamed push constants to "immediate data", tracking the WebGPU proposal.
        println!("  max_immediate_size (push constants)   : {} B", l.max_immediate_size);

        // Print the whole feature set rather than probing named constants: the set is what is
        // on the record, and a name that moved between wgpu versions must not silently read
        // as "unsupported".
        println!("  features : {:?}", a.features());

        let per_traj = N * N * COPIES;
        let per_fp = N * N;
        let cap = l.max_compute_invocations_per_workgroup;
        let staged = per_traj * SIMSTATE_B;
        let verdict = |v: u32, cap: u32| {
            if v <= cap { "fits".to_string() } else { format!("OVER by {:.2}x", v as f64 / cap as f64) }
        };
        println!("  -- fit at N={N}, E+1={COPIES} --");
        println!("     1 thread / trajectory ({per_traj:>4}) vs cap {cap:<5} : {}", verdict(per_traj, cap));
        println!("     1 thread / footprint  ({per_fp:>4}) vs cap {cap:<5} : {}", verdict(per_fp, cap));
        println!("     quad SimState in workgroup mem ({staged} B)  : {}",
            verdict(staged, l.max_compute_workgroup_storage_size));
    }
}
