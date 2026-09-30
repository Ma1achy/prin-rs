//! Dispatch the WGSL mirror and read the discrete words back.

use crate::{Discrete, N_MAX};
use wgpu::util::DeviceExt;

pub struct GpuRun {
    pub backend: String,
    pub adapter: String,
    pub driver: String,
    pub out: Vec<Discrete>,
}

pub fn run(d2s: &[f32], thr: &[f32; N_MAX as usize]) -> Option<GpuRun> {
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
    let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::None,
        force_fallback_adapter: false,
        compatible_surface: None,
        ..Default::default()
    }))
    .ok()?;
    let info = adapter.get_info();
    // Print what was actually exercised. The Metal-only caveat on the earlier spike exists because
    // this was not printed -- an assumed backend is not a measured one.
    let (backend, name, driver) = (
        format!("{:?}", info.backend),
        info.name.clone(),
        format!("{} {}", info.driver, info.driver_info),
    );

    let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
        label: None,
        required_features: wgpu::Features::empty(),
        required_limits: wgpu::Limits::downlevel_defaults(),
        memory_hints: wgpu::MemoryHints::default(),
        trace: wgpu::Trace::Off,
        ..Default::default()
    }))
    .ok()?;

    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("discrete"),
        source: wgpu::ShaderSource::Wgsl(include_str!("kernel.wgsl").into()),
    });

    let n = d2s.len();
    let mk = |data: &[u8], usage| {
        device.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: None, contents: data, usage })
    };
    let b_in = mk(bytemuck_cast(d2s), wgpu::BufferUsages::STORAGE);
    let b_thr = mk(bytemuck_cast(thr), wgpu::BufferUsages::STORAGE);
    let out_bytes = (n * 10 * 4) as u64;
    let b_out = device.create_buffer(&wgpu::BufferDescriptor {
        label: None,
        size: out_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    let b_read = device.create_buffer(&wgpu::BufferDescriptor {
        label: None,
        size: out_bytes,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: None,
        layout: None,
        module: &shader,
        entry_point: Some("main"),
        compilation_options: Default::default(),
        cache: None,
    });
    let bind = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: None,
        layout: &pipeline.get_bind_group_layout(0),
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: b_in.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: b_thr.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: b_out.as_entire_binding() },
        ],
    });

    let mut enc = device.create_command_encoder(&Default::default());
    {
        let mut pass = enc.begin_compute_pass(&Default::default());
        pass.set_pipeline(&pipeline);
        pass.set_bind_group(0, &bind, &[]);
        pass.dispatch_workgroups(((n + 63) / 64) as u32, 1, 1);
    }
    enc.copy_buffer_to_buffer(&b_out, 0, &b_read, 0, out_bytes);
    queue.submit([enc.finish()]);

    let slice = b_read.slice(..);
    slice.map_async(wgpu::MapMode::Read, |_| {});
    device.poll(wgpu::PollType::Wait { submission_index: None, timeout: None }).ok()?;
    let mapped = slice.get_mapped_range().ok()?;
    let words: Vec<u32> = mapped
        .chunks_exact(4)
        .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect();

    let out = words
        .chunks_exact(10)
        .map(|w| Discrete {
            bucket_frozen: w[0],
            bucket_trans: w[1],
            coll_sq: w[2],
            coll_sqrt: w[3],
            packed: w[4],
            roundtrip: w[5],
            pairsum_flag: w[6],
            pairsum_brk: w[7],
            packed_ctl: w[8],
            roundtrip_ctl: w[9],
        })
        .collect();

    Some(GpuRun { backend, adapter: name, driver, out })
}

fn bytemuck_cast(f: &[f32]) -> &[u8] {
    // f32 and u8 have no alignment conflict in this direction.
    unsafe { std::slice::from_raw_parts(f.as_ptr() as *const u8, std::mem::size_of_val(f)) }
}
