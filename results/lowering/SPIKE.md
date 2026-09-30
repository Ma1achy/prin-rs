# The spike: does the discrete surface survive lowering?

**Answer: yes, on two independent backends — 0 forks across every comparison-only arm over 705
boundary states, on Metal (naga→MSL, Apple) and on Mesa lavapipe (naga→SPIR-V, LLVM 20.1.2).** The
measurement is worth something because the control forked in both runs: 105 on Metal, 82 on lavapipe.

That narrows the standing action item — *"the spike's bit-identity evidence is Metal-only"* — to two
backends and two compilers. It does not close it: lavapipe is a software rasteriser, not an
AMD/NVIDIA/Qualcomm driver, and D3D12 stays unreachable.

**And the two backends fork on DIFFERENT STATES, which is the sharpest result here.** See §"No
strictest backend exists" below.

---

## What was lowered, and what was deliberately not

Four discrete things, each already real, mirrored in `src/lib.rs` (Rust) and `src/kernel.wgsl` (WGSL):

| thing | shape | where it comes from |
|---|---|---|
| `pack` | `((state & 7) << 2) \| (detail & 3)` | copied from `src/outcome.rs:157`, not imported |
| `from_bits` | the six-variant code table | `src/outcome.rs:128` |
| collision | `d² < r_coll²`, squared form, build-time constant | the **contract's** shape, not prin-rs's |
| `substep_bucket` | `d²` against 64 frozen f32 thresholds | `principia_integrator_contract.md` — **absent from prin-rs**, and the mechanism the contract relies on to make Tier B achievable, so the thing most worth testing |

**Scope stated rather than implied.** The WGSL is hand-written, not rust-gpu. Both routes go through
naga, so this exercises the standing gap — *a different naga backend and driver* — and does **not**
re-test the shared-source toolchain, which the original spike already cleared. A rust-gpu leg would
answer a different question and is not what the action item names.

---

## The table

`spikes/lowering/src/bin/parity.rs`, 705 states: ±4 ulp around each of the 64 `ceil` boundaries,
plus ±64 ulp around the collision threshold. Raw output in `parity_metal.txt`.

```
field           rule       f32-vs-f64   gpu-vs-f32
bucket_frozen   clean               0            0
bucket_trans    CONTROL            22           83
coll_sq         clean               0            0
coll_sqrt       CONTROL             0            0
packed          clean               0            0
roundtrip       clean               0            0
pairsum_flag    clean               0            0
pairsum_brk     CONTROL             0            0
packed_ctl      propagate          22           83
roundtrip_ctl   propagate           0            0
```

`gpu` is Metal on an Apple M3 Pro; the lavapipe column reads the same except `bucket_trans` 60 and
`packed_ctl` 60 (`parity_lavapipe.txt`, adapter `llvmpipe (LLVM 20.1.2, 256 bits)`, Mesa 25.2.8).
`adapter.get_info()` is printed by the binary rather than assumed, because the Metal-only caveat on
the earlier spike exists precisely because that was not printed.

**The clean arms read 0 in both columns.** The frozen-threshold bucket, the squared collision test,
the packing, the round-trip and the flag-tested loop all agree bitwise across f64, f32 and Metal.

**The `bucket_trans` control forks 22 times CPU-to-CPU and 83 times CPU-to-GPU.** That is the
original spike's failure reproduced independently, on this machine, in this harness — a runtime
`ceil((r_sub/√d²)^1.5)`. So the harness *can* see a fork. Without that column the zeros above would
be a test that cannot fail.

---

## The propagation arm is the part worth quoting

`packed` and `packed_ctl` are the **same packing expression** — identical shifts and masks — fed by
different buckets. `packed` reads 0 forks; `packed_ctl` reads 83.

**A packed descriptor is exactly as deterministic as the branches feeding it, and no more.** The
bit-packing is not where parity is won or lost; it never was. Auditing the packing would find
nothing, on either arm.

And `roundtrip_ctl` reads **0 forks on a contaminated input** — `from_bits` masks to bits 2-4, and
the fork lands in bits 0-1. **A decoder can round-trip perfectly on a value that is already wrong.**
Any parity check placed at the decode would pass here.

---

## No strictest backend exists, so validating against one GPU cannot certify another

The control's fork count differs between backends — 83 states on Metal, 60 on lavapipe, from the
**identical** 705 inputs. A count alone would admit a comfortable reading: that one backend is
simply looser and validating against it covers the other. The sets say otherwise.

```
lavapipe 60   metal 83
intersection 26   lavapipe-only 34   metal-only 57
lavapipe subset of metal?  False
```

**Only 26 states fork on both.** 57 fork on Metal alone and 34 on lavapipe alone — the sets cross,
neither nests. So there is no strictest backend to test against and no single GPU whose agreement
implies another's. **Every backend has states the others get right**, and a validation suite pinned
to one of them would pass while a second backend disagreed on 34 inputs it never sampled.

This is the same shape as the limits table one level up: *there is no single binding backend* —
WebGPU binds on invocations, Metal on workgroup storage. Here it is the same conclusion measured on
arithmetic rather than read off a spec, and it is the argument for the comparison-only rule as a
**construction constraint** rather than a thing to check for after the fact. A rule you can verify
per-backend would be a weaker rule; this one has to hold by construction because no finite set of
backends can stand in for the rest.

Fork sets are printed by the harness and are in `parity_metal.txt` and `parity_lavapipe.txt`.

## Two controls read zero, and they read zero for different reasons

A silent control is a failed measurement until its cause is established. Both were chased.

**`coll_sqrt` — a control with a subject, that genuinely did not fork.** The first cut of the sweep
never approached `r_coll` at all: it targeted `ceil` boundaries, so the sqrt-decided collision test
had nothing to decide and its zero meant nothing. With ±64 ulp around the threshold added, the sweep
puts **3 of 705 states within 1 ULP of `r_coll`, nearest 0 ULP** — measured, not argued — and it
still reads 0. WGSL specifies `sqrt` to **1 ULP, not correctly-rounded**, so it *may* fork on another
backend; on Metal it did not, at a state where one ulp would have flipped the branch.

**That refines `AUDIT.md`'s severity ranking with a measurement.** The audit lists 46 sqrt-decided
branch sites as R1 violations. They are still violations — the specification permits the latitude —
but the measured fork came from `pow`, at 83 of 705, and not from `sqrt`, at 0 of 705 with the
threshold straddled. **`pow` is the one that bites; `sqrt` is the one that is merely unguaranteed.**

**`pairsum_brk` — a control whose historical bug appears to be fixed.** This is the triangular nested
loop with a mid-body `break`, the shape that summed one body-pair of three through naga→MSL
(`H = −0.577` against `−1.732`), silently. It reads 0 forks. Agreement alone would not settle it:
**two backends agreeing on a wrong answer is the miscompile happening identically on both.** So both
loop shapes are checked against the closed form `3·seed + 21`, and both are **0 wrong of 705**. The
break-form is correct here, not merely consistent. The bug is absent at `wgpu` 30.0.1 / naga on
Metal; rule 5 remains cheap insurance, and this is evidence about one compiler version on one
backend, not a repeal.

---

## Verdict semantics, and why the CI gate reads the text

The binary prints one of three verdicts and **exits 0 for all three, deliberately**:

- `PASS` — clean arms 0 forks **and** controls non-zero.
- `INCONCLUSIVE` — clean arms 0 forks and controls **also** 0. Not a pass. The harness is
  insensitive and the zeros are not evidence.
- `FAIL` — a comparison-only arm forked.

A conventional exit-code gate cannot tell `PASS` from `INCONCLUSIVE`, so a blind control would ship
as a green tick. The workflow greps for `^verdict: PASS`.

---

## What this does not answer

- **Real drivers.** lavapipe is a software rasteriser. Two compilers is better than one and is not
  a driver survey; the crossing fork sets are the direct evidence that a third would disagree with
  both. D3D12 is unreachable from macOS or Linux and remains open.
- **prin-rs's own integrator.** The spike lowers the *contract's* step control. prin-rs's
  `dtau`/`StepLimit::Predictive` is division-decided throughout (`AUDIT.md` §2), and `src/real.rs`
  gives f32 and f64 different branch thresholds, so **Tier B over prin-rs's kernel is unachievable by
  construction whatever this table says**. The surface lowers; the shipped kernel is not on it.
- **The dispatch shape still does not fit** (`README.md` §3). This spike says the arithmetic survives;
  it says nothing about the four ways the shape exceeds the WebGPU baseline.

## Reproduction

```sh
cd spikes/lowering && $HOME/.cargo/bin/cargo run --release --bin parity
```

`spikes/lowering/` is a throwaway crate, deliberately **not** a member of prin-rs, which has no
`[workspace]` and auto-discovers only `src/`, `examples/`, `tests/`, `benches/`. **Nothing in prin-rs
is modified** — `git status` over `src/`, `tests/`, `examples/` and `Cargo.toml` is empty.
