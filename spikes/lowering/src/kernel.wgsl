// The discrete surface, mirrored from src/lib.rs. Hand-written WGSL rather than rust-gpu: both
// routes go through naga, so this exercises the standing gap (a different naga backend and driver)
// and does NOT re-test the shared-source toolchain, which the original spike already cleared.

const N_MAX: u32 = 64u;
const R_SUB: f32 = 0.05;
const GAMMA: f32 = 1.5;

@group(0) @binding(0) var<storage, read>       d2s: array<f32>;
@group(0) @binding(1) var<storage, read>       thr: array<f32>;   // 64 frozen thresholds
@group(0) @binding(2) var<storage, read_write> out: array<u32>;   // 10 words per input

// --- comparison-only arm -----------------------------------------------------------------------
// No sqrt, no divide, no pow. Every branch input is a value that arrived as bits.
// Loop shape is the contract's canonical wrapper: `done` flag tested in the condition, zero break.
fn bucket_frozen(d2: f32) -> u32 {
    var n: u32 = 1u;
    var i: u32 = 0u;
    var done: bool = false;
    while (!done) {
        if (d2 < thr[i]) { n = i + 2u; }
        i = i + 1u;
        done = i >= N_MAX;
    }
    if (n > N_MAX) { return N_MAX; }
    return n;
}

// --- NEGATIVE CONTROL: R1-violating twin ---------------------------------------------------------
// sqrt, divide and pow all feeding a ceil whose integer output is a branch decision. This is the
// shape that forked 5 of 272 in the original spike. If this agrees too, the harness is blind.
fn bucket_trans(d2: f32) -> u32 {
    let r_min = sqrt(d2);
    let q = R_SUB / r_min;
    let p = pow(q, GAMMA);
    var c = ceil(p);
    if (c < 1.0) { c = 1.0; }
    if (c > f32(N_MAX)) { c = f32(N_MAX); }
    return u32(c);
}

fn pack(state: u32, detail: u32) -> u32 { return ((state & 7u) << 2u) | (detail & 3u); }

fn from_bits(code: u32) -> u32 {
    let s = (code >> 2u) & 7u;
    if (s <= 5u) { return s; }
    return 255u;
}

// --- loop-shape arms -----------------------------------------------------------------------------
// The triangular nested loop over body pairs, in both shapes. Integer arithmetic, so the answer is
// exact and any disagreement is a miscompile rather than rounding. The `_brk` twin is the shape
// that summed one pair of three through naga->MSL with no error and no crash.
fn pairsum_flag(seed: u32) -> u32 {
    var acc: u32 = 0u;
    var i: u32 = 0u;
    var di: bool = false;
    while (!di) {
        var j: u32 = i + 1u;
        var dj: bool = false;
        while (!dj) {
            acc = acc + (seed + i * 16u + j);
            j = j + 1u;
            dj = j >= 3u;
        }
        i = i + 1u;
        di = i >= 2u;
    }
    return acc;
}

// NEGATIVE CONTROL: counted `for` with a mid-body `break`, two levels.
fn pairsum_brk(seed: u32) -> u32 {
    var acc: u32 = 0u;
    for (var i: u32 = 0u; i < 3u; i = i + 1u) {
        if (i >= 2u) { break; }
        for (var j: u32 = 0u; j < 3u; j = j + 1u) {
            if (j <= i) { continue; }
            if (j >= 3u) { break; }
            acc = acc + (seed + i * 16u + j);
        }
    }
    return acc;
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= arrayLength(&d2s)) { return; }
    let d2 = d2s[k];

    let r_coll: f32 = 0.02;
    let bf = bucket_frozen(d2);
    let bt = bucket_trans(d2);
    let state = bf % 6u;
    let packed = pack(state, bf & 3u);        // fed only by comparison-only values
    let packed_ctl = pack(state, bt & 3u);    // fed by the R1-violating bucket: propagation arm

    let b = k * 10u;
    out[b + 0u] = bf;
    out[b + 1u] = bt;
    out[b + 2u] = select(0u, 1u, d2 < r_coll * r_coll);
    out[b + 3u] = select(0u, 1u, sqrt(d2) < r_coll);   // NEGATIVE CONTROL: sqrt-decided collision
    out[b + 4u] = packed;
    out[b + 5u] = from_bits(packed);
    out[b + 6u] = pairsum_flag(bf);
    out[b + 7u] = pairsum_brk(bf);
    out[b + 8u] = packed_ctl;
    out[b + 9u] = from_bits(packed_ctl);
}
