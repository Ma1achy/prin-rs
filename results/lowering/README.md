# Does the kernel actually lower? — backend limits, and the one answer

**Question asked, before the implementation plan is built on the assumption.** Not precision — f64
on CPU and f32 on GPU is the stated design and the linearised decoder already handles the f32 floor.
The question is the **discrete parity surface**: packed descriptors, integer step counts, enum
values, word arithmetic. Integers and comparisons, so they *should* lower identically.

**The answer to the one question — does the current dispatch shape fit inside the intersection of
all backends — is NO, in four independent ways, and three of them are invisible on Metal.**

---

## 0. Two corrections to the premise, both from the design corpus itself

**A spike already ran, and the comparison-only rule was earned rather than specified.**
`principia_gpu_determinism_note.md` records each of its six rules as a *measured* failure. The
f32-quantise hypothesis **forked 5 of 272** boundary states; a controlled test ruled out **both**
fast-math and FMA contraction as the cause and identified inherent transcendental latitude — GPU
`pow(4.000000477, 1.5)` → exactly `8.0` against libm's `8.000001907`, a 2-ulp gap landing on a
`ceil` boundary. And rule 5 cost a **silent miscompile**: a triangular nested loop through
naga→MSL summed one body-pair of three (`H = −0.577` against `−1.732`), no error, no crash.

So of the five risks named in the question, **#1 is answered by the contract** (the rule exists
because it was needed) and **#2 is answered and negative** — contraction was excluded as the cause,
though rules 2 and 6 still mandate explicit `fma` for a different reason.

**What is untested is what that note already names as the standing pre-Paper-2 action item:** *"the
spike's bit-identity evidence is Metal-only… branch-word bit-identity across a different naga
backend + driver is the one thing the spike could not witness on an M-series machine."* That is the
real gap, and §2 shows this machine cannot close it.

---

## 1. The limits table

**Provenance is per-row, and it earned its keep immediately.** A first pass cited WebGPU's
`maxComputeWorkgroupStorageSize` as 49152 from a spec fetch. The **measured** Metal value is 32768 —
and Dawn's WebGPU runs *on* Metal, so a WebGPU default above 32768 is impossible. The measurement
falsified the citation. Re-sourced from `wgpu-types-30.0.1/src/limits.rs:421-462`, the real baseline
is **16384**; the same fetch was also wrong about storage-buffer and buffer size (it returned Metal's
numbers, not the baseline's). **Every number below says where it came from.**

| limit | **Metal, measured here** | **WebGPU baseline, cited** | Vulkan required min, cited |
|---|---|---|---|
| `maxComputeInvocationsPerWorkgroup` | **1024** | **256** | 128 |
| `maxComputeWorkgroupSize` x/y/z | 1024 / 1024 / 1024 | 256 / 256 / 64 | 128 / 128 / 64 |
| `maxComputeWorkgroupStorageSize` | **32768 B** | **16384 B** | 16384 B |
| `maxComputeWorkgroupsPerDimension` | 65535 | 65535 | 65535 |
| `maxStorageBuffersPerShaderStage` | 29 | **8** | 4 |
| `maxStorageBufferBindingSize` | 4294967292 B | 128 MiB | 128 MiB |
| `maxBufferSize` | 4294967295 B | 256 MiB | — |
| `maxBindGroups` | 8 | **4** | 4 (descriptor sets) |
| `maxUniformBufferBindingSize` | 4294967295 B | 64 KiB | 16 KiB |
| push constants (`maxImmediateSize`) | **256 B** | **0 B — absent** | 128 B |
| f64 in shader | **absent** | absent | optional |
| f16 in shader | present | optional extension | optional |
| **64-bit integers** | **`SHADER_INT64` present** | **absent — WGSL has i32/u32 only** | optional |
| **float32 atomics** | **`SHADER_FLOAT32_ATOMIC` present** | **absent — `atomic<T>` is u32/i32 only** | extension |

- **Metal column**: measured on this machine (Apple M3 Pro, macOS 26.2, Metal 4) by
  `spikes/lowering/src/bin/limits.rs`; raw output in `metal_limits.txt`.
- **WebGPU column**: `wgpu-types-30.0.1/src/limits.rs:421-462` (`Limits::defaults()`, the WebGPU
  baseline every conformant implementation must meet) and the WGSL spec for the type system.
- **Vulkan column**: the spec's required-minimums table. **Theoretical** — real desktop drivers
  exceed it universally. WebGPU's is the floor that actually binds a shipped web artefact.

---

## 2. Only one backend is reachable from this machine, and that is the finding

`enumerate_adapters(Backends::all())` returns **1 adapter: Metal**. No Vulkan SDK, no MoltenVK, no
D3D12 (impossible on macOS). Chrome is installed, but its WebGPU goes Dawn→Metal — **the same MSL
compiler**, which is precisely why the original spike's evidence is Metal-only and why running it
again here would reproduce the limitation rather than close it.

So the cross-backend leg belongs in **CI on Linux with Mesa lavapipe** (software Vulkan): it
exercises naga's SPIR-V backend and a non-Apple compiler, needs no hardware, and CI is where the
determinism note says the gate lives permanently.

---

## 3. The answer

Dispatch shape, from the repo rather than assumed: `n: 8` (`src/scheduler.rs:419`), `n_extra: 7` → 8
copies (`src/ensemble/pixel.rs:204`), and the code's own comment at `src/scheduler.rs:252` — *"At
`N=8`, `E+1=8` one quad is 512 trajectories"*. The lowering contract's intended mapping is
`workgroups: samplesPerQuad(tier)`, N×N per quad. Hot `SimState` is 136 B (memory-tiers).

| what | needed | WebGPU | Metal |
|---|---|---|---|
| one invocation per **trajectory** | 512 | **✗ 2.00× over** | ✓ fits (1024) |
| one invocation per **footprint**, looping 8 copies | 64 | ✓ | ✓ |
| a quad's `SimState` staged in workgroup memory | 69632 B | **✗ 4.25× over** | **✗ 2.12× over** |
| per-quad params as push constants | any | **✗ feature does not exist** | ✓ 256 B |
| `u64` counters (`total_substeps`, `n_cap_hits`, …) | required | **✗ no 64-bit integers** | ✓ present |
| float atomics for the `max` reductions | wanted | **✗ integer atomics only** | ✓ present |

**The natural mapping fits on Metal and fails on WebGPU — and the spike ran on Metal.** That is the
same shape as the standing caveat, one level down: the machine that is easiest to test on is the one
whose limits are loosest, so it certifies nothing about the target.

**Note the inversion, because it means there is no single binding backend.** WebGPU binds on
invocations (256 against 1024) and on capability; **Metal binds on workgroup storage** (32768 against
WebGPU's… also-insufficient 16384 — both fail, Metal by less). Quoting "WebGPU is the constraint"
alone would be wrong.

### What each failure costs

1. **512 invocations** — not fatal. One invocation per *footprint* (64) looping the 8 copies fits
   everywhere with room. But the natural "one workgroup = one quad, one thread = one trajectory"
   mapping is unavailable, and the copy loop becomes serial per thread.
2. **Workgroup storage** — not fatal; the quad's payload lives in a storage buffer, not workgroup
   memory. It does mean no whole-quad staging for the reduction.
3. **No push constants in WebGPU** — the lowering contract's per-quad `quadUniforms` must be a
   uniform buffer with dynamic offsets. A mechanism change, not a design change.
4. **No 64-bit integers in WGSL, and no float atomics** — these are the real ones, and they are
   **design changes rather than ports**:
   - Seven payload counters are `u64` (`src/ensemble/pixel.rs:487-607`) and `MarchOut::steps` /
     `force_evals` are **`usize`** (`src/integrate/mod.rs:324,328`). Every one is a Tier-B
     bit-exact quantity. They become `u32`, or hi/lo pairs, or they do not exist on the GPU.
   - The reductions are `max` over `f64` (`src/scheduler.rs:544-570`), so they become workgroup
     tree reductions rather than atomic max.
   - **And both work natively on Metal.** A spike written here would use `SHADER_INT64` and
     `SHADER_FLOAT32_ATOMIC` without noticing, and break in a browser.

### One more collision, between the repo and the validator

`+inf` is **load-bearing** in those reductions as the absorbing element — `src/ensemble/pixel.rs:1077-1080`,
the no-discard fix, where a max-fold must saturate to `+inf` rather than drop an undetermined
footprint. The determinism note's rule 4 records that **naga rejects infinity literals** ("Float
literal is infinite" — WGSL cannot express it). The semantics the repo deliberately adopted and the
validator's constraint meet head-on, and the fix is not mechanical: seeding with the first element
changes what an all-undetermined fold returns.

---

## 4. Reproduction

```sh
cd spikes/lowering && $HOME/.cargo/bin/cargo run --release --bin limits
```

`spikes/lowering/` is a throwaway crate, deliberately **not** a member of prin-rs — prin-rs has no
`[workspace]` and auto-discovers only `src/`, `examples/`, `tests/`, `benches/`, so the directory is
invisible to it and **nothing in prin-rs is modified**.

The probe prints the whole feature set rather than probing named constants: `wgpu` 30 renamed push
constants to "immediates" and split `Features` into WGPU/WebGPU halves, and a constant that moved
between versions must not silently read as "unsupported".

## 5. What this does not answer

- **Real-driver behaviour.** lavapipe (the CI leg) is a software rasteriser. It exercises naga's
  SPIR-V backend and a non-Apple compiler; it is not an AMD/NVIDIA/Qualcomm driver. It narrows the
  standing action item and does not close it.
- **D3D12** — unreachable from macOS or Linux. Open.
- **Whether the branch words actually agree.** That is the spike (`SPIKE.md`), not the limits
  table. This document says the shape does not fit; it does not say the arithmetic forks.
