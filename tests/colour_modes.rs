//! **The two specced default colour modes.** Fidelity pinned by values, not by types.
//!
//! §7 asks for golden images over a preset table. `principia_colour_composition.md` is **not in
//! this repo**, so the preset table is unavailable and these pin the properties the spec states
//! in prose instead: the exact palette bytes, the mnemonic structure, the kernel limit, the
//! desaturation feature, the hoist, the invalid colour, and the inverted `t_end` polarity. When
//! the document lands, §7's table replaces the palette test and the rest stand.
use prin_rs::output::colour::{self, Scalar};
use prin_rs::output::palette::{self, EventClass, Filter, Swatches};
use prin_rs::output::siteblend::{
    BlendSpace, Brightness, Kernel, PhysicsGen, SiteBlend, Sites, StaticGen,
};

// -------------------------------------------------------------------------------------------
// MODE 1
// -------------------------------------------------------------------------------------------

/// The palette bytes, exactly. **A generated palette cannot encode the mnemonic**, so this is
/// pinned rather than described.
#[test]
fn the_canonical_palette_is_exact_and_the_mnemonic_holds() {
    let sw = Swatches::default();
    let expect = [
        (EventClass::Collision12, [0xDE, 0x2D, 0x2D], "red"),
        (EventClass::Collision13, [0x2E, 0xBC, 0x4E], "green"),
        (EventClass::Collision23, [0x34, 0x62, 0xE0], "blue"),
        (EventClass::Escape1, [0xF0, 0xDE, 0x32], "yellow"),
        (EventClass::Escape2, [0xE0, 0x34, 0xC6], "magenta"),
        (EventClass::Escape3, [0x30, 0xC8, 0xDC], "cyan"),
        (EventClass::Bounded, [0x14, 0x14, 0x18], "black"),
        (EventClass::CollisionAtZero, [0xF2, 0x96, 0x20], "orange"),
        (EventClass::Degenerate, [0xEC, 0xEC, 0xF0], "white"),
    ];
    for (c, rgb, name) in expect {
        assert_eq!(sw.get(c), rgb, "{} should be {name}", c.name());
    }

    // The mnemonic itself: collisions are ADDITIVE primaries, escapes SUBTRACTIVE. An additive
    // primary has ONE dominant channel; a subtractive primary has TWO. That is the structure a
    // golden-angle cycle cannot reproduce, so it is asserted rather than trusted.
    let dominant = |c: EventClass| {
        let v = sw.get(c);
        let mx = *v.iter().max().unwrap() as i32;
        v.iter().filter(|&&x| (x as i32) > mx / 2).count()
    };
    for c in [EventClass::Collision12, EventClass::Collision13, EventClass::Collision23] {
        assert_eq!(dominant(c), 1, "{} is additive: one dominant channel", c.name());
    }
    for c in [EventClass::Escape1, EventClass::Escape2, EventClass::Escape3] {
        assert_eq!(dominant(c), 2, "{} is subtractive: two dominant channels", c.name());
    }
}

/// The swatch table is a **default, not a fixture**.
#[test]
fn the_swatch_set_is_user_editable() {
    let mut sw = Swatches::default();
    sw.set(EventClass::Bounded, [1, 2, 3]);
    assert_eq!(sw.get(EventClass::Bounded), [1, 2, 3]);
    assert_eq!(sw.get(EventClass::Collision12), [0xDE, 0x2D, 0x2D], "editing one must not move others");
}

/// **Muted and invalid must be distinguishable.** Not selected and not known are different
/// things, and a reader has to be able to tell them apart.
#[test]
fn a_muted_class_and_an_invalid_pixel_do_not_render_alike() {
    let f = Filter::collisions();
    let sw = Swatches::default();
    assert!(f.passes(EventClass::Collision12));
    assert!(!f.passes(EventClass::Escape1), "the collisions filter must mute escapes");
    assert_ne!(f.muted, sw.get(EventClass::Invalid), "muted must differ from invalid");
    // And the filter is a general operation: escapes-only is the same machinery.
    assert!(Filter::escapes().passes(EventClass::Escape2));
    assert!(!Filter::escapes().passes(EventClass::Collision12));
}

/// **The inverted polarity, and it is what keeps bounded black.** `t_end` white = LOW / EARLY.
#[test]
fn t_end_brightness_is_early_is_white_and_ftle_is_the_opposite() {
    let (lo, hi) = (0.0, 50.0);
    let early = colour::range_norm(Scalar::TEnd, 1.0, lo, hi).unwrap();
    let late = colour::range_norm(Scalar::TEnd, 49.0, lo, hi).unwrap();
    println!("t_end: early {early:.4}, late {late:.4}");
    assert!(early > late, "t_end must render EARLY as bright: {early} vs {late}");
    assert!(colour::lightness(early) > colour::lightness(late));

    // A footprint that never resolves sits at t_max and takes the darkest lightness, which is
    // what makes `bounded = black` fall out of the ramp rather than need a special case.
    let never = colour::range_norm(Scalar::TEnd, hi, lo, hi).unwrap();
    assert!(never <= 1e-12, "a never-resolving footprint must be darkest, got {never}");

    // The opposite convention, on the same machinery. The inconsistency is deliberate.
    let f_lo = colour::range_norm(Scalar::Ftle, 0.1, 0.0, 5.0).unwrap();
    let f_hi = colour::range_norm(Scalar::Ftle, 4.9, 0.0, 5.0).unwrap();
    assert!(f_hi > f_lo, "FTLE must render HIGH as bright: {f_hi} vs {f_lo}");
}

/// The legend states the polarity. A reader must not have to infer it.
#[test]
fn the_legend_states_the_inverted_polarity() {
    let t = palette::legend(&Swatches::default(), &Filter::default());
    assert!(t.contains("WHITE = LOW / EARLY"), "legend must state the polarity:\n{t}");
    assert!(t.contains("OPPOSITE"), "legend must flag the inconsistency with FTLE:\n{t}");
    assert!(t.contains("collision 1-2") && t.contains("body 1 escape"));
}

// -------------------------------------------------------------------------------------------
// MODE 2
// -------------------------------------------------------------------------------------------

/// **`Nearest` is a discrete path, not a large `kappa`.** The point of the variant is that the
/// limit is unreachable numerically: `exp` overflows first.
#[test]
fn nearest_is_a_discrete_path_and_a_large_kappa_is_not_a_substitute() {
    let d = [0.9f64, 0.2, -0.5];
    let w = Kernel::Nearest.weights(&d).unwrap();
    assert_eq!(w, vec![1.0, 0.0, 0.0], "nearest is an indicator on the argmax");

    // The dial's hard detent.
    assert_eq!(Kernel::from_dial(1.0, None), Kernel::Nearest);
    assert!(matches!(Kernel::from_dial(0.99, None), Kernel::Vmf(_)));

    // vMF approaches it and stays finite because `d_max` is subtracted -- without that shift
    // `exp(kappa*d)` overflows and the weights come back NaN rather than Voronoi.
    let w12 = Kernel::Vmf(12.0).weights(&d).unwrap();
    assert!(w12.iter().all(|x| x.is_finite()), "vMF at the top of its range must stay finite");
    assert!(w12[0] > w12[1] && w12[1] > w12[2]);
    let huge = Kernel::Vmf(1e4).weights(&d).unwrap();
    assert!(huge.iter().all(|x| x.is_finite()), "the shift must hold even far past the range");
}

/// `TopK` truncates support; `Vmf` does not. That is the `support` axis, separate from
/// temperature.
#[test]
fn topk_truncates_support_and_vmf_does_not() {
    let d = [0.9f64, 0.5, 0.1, -0.3];
    let v = Kernel::Vmf(3.0).weights(&d).unwrap();
    assert!(v.iter().all(|&x| x > 0.0), "vMF gives every site weight");
    let t = Kernel::TopK(2, 3.0).weights(&d).unwrap();
    assert_eq!(t.iter().filter(|&&x| x > 0.0).count(), 2, "top-2 keeps exactly two");
    assert!(t[0] > 0.0 && t[1] > 0.0 && t[2] == 0.0 && t[3] == 0.0);
}

/// **Uncertainty reads as desaturation, and it is a feature.** A direction equidistant from two
/// oppositely-coloured sites must blend to lower chroma than one sitting on a site.
#[test]
fn oklab_blending_desaturates_toward_a_site_boundary() {
    let blend = SiteBlend {
        sites: Sites::Static(StaticGen::Axes6),
        kernel: Kernel::Vmf(4.0),
        // +x and -x get opposite (a,b); the rest are neutral.
        colours: vec![
            [0.7, 0.2, 0.0], [0.7, -0.2, 0.0],
            [0.7, 0.0, 0.0], [0.7, 0.0, 0.0],
            [0.7, 0.0, 0.0], [0.7, 0.0, 0.0],
        ],
        space: BlendSpace::Oklab,
    };
    let m = [1.0 / 3.0; 3];
    let on_site = blend.blend([1.0, 0.0, 0.0], &m).unwrap();
    let between = blend.blend([0.0, 0.0, 1.0], &m).unwrap();
    let chroma = |c: [f64; 3]| c[1].hypot(c[2]);
    println!("chroma on-site {:.5}, between {:.5}", chroma(on_site), chroma(between));
    assert!(
        chroma(between) < chroma(on_site),
        "a boundary direction must desaturate: {} vs {}",
        chroma(between),
        chroma(on_site)
    );
}

/// **The hoist is exact.** `decoder::Latent` puts `z_mu` at indices 6 and 7.
#[test]
fn physics_sites_hoist_only_when_no_axis_touches_a_mass_dimension() {
    let phys = Sites::Physics(PhysicsGen::Bc);
    let stat = Sites::Static(StaticGen::Ico12);
    let (mut q1, mut q2) = ([0.0f64; 8], [0.0f64; 8]);
    q1[0] = 1.0;
    q2[1] = 1.0;
    assert!(phys.hoistable_over(&q1, &q2), "a shape-only slice holds the masses constant");
    assert!(stat.hoistable_over(&q1, &q2), "static sites always hoist");

    q2[6] = 0.3; // z_mu1
    assert!(!phys.hoistable_over(&q1, &q2), "an axis touching z_mu makes masses per-pixel");
    assert!(stat.hoistable_over(&q1, &q2), "static sites are unaffected by a mass axis");

    // And the semantics are per-pixel regardless: the sites really do move with the masses, so
    // hoisting is an optimisation and not a different answer.
    let a = phys.points(&[0.4, 0.35, 0.25]);
    let b = phys.points(&[0.2, 0.3, 0.5]);
    assert!(a != b, "physics sites must move with the masses or the hoist is vacuous");
}

/// Every static generator returns unit vectors of the advertised count.
#[test]
fn static_generators_are_unit_and_the_right_size() {
    for (g, n) in [
        (StaticGen::Axes6, 6),
        (StaticGen::Corner8, 8),
        (StaticGen::Ico12, 12),
        (StaticGen::Fib(37), 37),
        (StaticGen::Ring(9, 0.4, 0.1), 9),
    ] {
        let p = Sites::Static(g).points(&[1.0 / 3.0; 3]);
        assert_eq!(p.len(), n, "{g:?}");
        for v in &p {
            let r = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            assert!((r - 1.0).abs() < 1e-12, "{g:?} produced a non-unit direction: {r}");
        }
    }
}

/// **The window is fixed and shareable, and auto-range announces itself.**
#[test]
fn the_brightness_window_is_fixed_by_default_and_auto_range_declares_itself() {
    let b = Brightness::fixed(Scalar::Ftle, 0.0, 5.0);
    assert!(!b.auto_ranged);
    let p = b.provenance();
    assert!(p.contains("Ftle") && p.contains("window="), "{p}");
    assert!(!p.contains("AUTO-RANGED"), "a fixed window must not claim to be auto-ranged");

    let auto = Brightness { auto_ranged: true, ..b };
    assert!(
        auto.provenance().contains("AUTO-RANGED"),
        "auto-range must announce itself: {}",
        auto.provenance()
    );
}

/// Both modes emit a provenance line carrying the axes the spec names.
#[test]
fn mode_two_declares_kernel_temperature_space_and_sites() {
    let blend = SiteBlend {
        sites: Sites::Physics(PhysicsGen::Lagrange),
        kernel: Kernel::TopK(3, 7.5),
        colours: vec![[0.7, 0.1, 0.0], [0.7, -0.1, 0.0]],
        space: BlendSpace::Oklab,
    };
    let p = blend.provenance();
    println!("{p}");
    for want in ["mode=siteblend", "topk", "k=3", "ks=7.5", "space=Oklab", "physics(Lagrange)"] {
        assert!(p.contains(want), "provenance is missing {want}: {p}");
    }
}

/// An undetermined direction yields `None`, never a valid-looking blend.
#[test]
fn an_undetermined_direction_is_none_and_never_a_uniform_fallback() {
    let blend = SiteBlend {
        sites: Sites::Static(StaticGen::Axes6),
        kernel: Kernel::Vmf(4.0),
        colours: vec![[0.7, 0.0, 0.0]; 6],
        space: BlendSpace::Oklab,
    };
    let m = [1.0 / 3.0; 3];
    assert!(blend.blend([f64::NAN, 0.0, 0.0], &m).is_none());
    assert!(Kernel::Vmf(4.0).weights(&[f64::NAN, 0.1]).is_none());
    assert!(Kernel::Nearest.weights(&[]).is_none());
}

// -------------------------------------------------------------------------------------------
// THE BACKBONE — Option occupants, the truth table, and one supersampler
// -------------------------------------------------------------------------------------------

use prin_rs::output::compose::{self, Combiner, NEUTRAL};

/// §4.1's truth table, all four cells, both combiners. **`None` is the identity element.**
#[test]
fn none_is_the_identity_element_of_combine() {
    let c = [0.42, 0.10, -0.05];
    let b = 0.80;

    // C + B
    assert_eq!(compose::combine(Some(c), Some(b), Combiner::ReplaceL), [b, c[1], c[2]]);
    // C + None -- the colour keeps its OWN lightness. "Just the colour map", the default.
    assert_eq!(compose::combine(Some(c), None, Combiner::ReplaceL), c);
    assert_eq!(compose::combine(Some(c), None, Combiner::Multiply), c);
    // None + B -- greyscale of the brightness field, and the most CVD-robust encoding possible.
    assert_eq!(compose::combine(None, Some(b), Combiner::ReplaceL), [b, 0.0, 0.0]);
    assert_eq!(compose::combine(None, Some(b), Combiner::Multiply), [b, 0.0, 0.0]);
    // None + None -- well-defined, harmless, instantly visible as "nothing is wired here".
    assert_eq!(compose::combine(None, None, Combiner::ReplaceL), NEUTRAL);
    assert_eq!(compose::combine(None, None, Combiner::Multiply), NEUTRAL);

    // Replace-L keeps the channels independent; Multiply does not, which is why it is not the
    // default: under Multiply the site's own L scales the scalar.
    let m = compose::combine(Some(c), Some(b), Combiner::Multiply);
    assert!((m[0] - c[0] * b).abs() < 1e-15, "multiply scales L");
    assert!(m[1] != c[1], "multiply also scales chroma, so the channels are NOT independent");
}

/// **The supersampler is independent of the colouring.** It is handed a closure and knows nothing
/// else — so a map it has never heard of gets supersampling by passing one.
#[test]
fn resolve_is_independent_of_what_produced_the_samples() {
    // A colour function with no relation to any map in this project.
    let made_up = |i: usize| Some([0.5, 0.1 * i as f64, -0.2 * i as f64]);
    let r = compose::resolve(4, made_up).unwrap();
    // mean of i = 0..3 is 1.5
    assert!((r[0] - 0.5).abs() < 1e-15);
    assert!((r[1] - 0.15).abs() < 1e-15, "a = 0.1 * mean(i) = 0.15, got {}", r[1]);
    assert!((r[2] + 0.30).abs() < 1e-15, "b = -0.2 * mean(i) = -0.30, got {}", r[2]);
}

/// **One undetermined sub-sample makes the pixel undetermined.** The no-discard rule, at the
/// resolve path — averaging the survivors is what biases a chaos instrument toward the tame.
#[test]
fn one_undetermined_sub_sample_makes_the_pixel_undetermined() {
    let with_hole = |i: usize| if i == 2 { None } else { Some([0.5, 0.0, 0.0]) };
    assert!(compose::resolve(4, with_hole).is_none(), "a hole must not be averaged away");
    assert!(compose::resolve(0, |_| Some([0.5, 0.0, 0.0])).is_none(), "no samples is undetermined");
    // And the control: the same function without the hole resolves.
    assert!(compose::resolve(4, |_| Some([0.5, 0.0, 0.0])).is_some());
}

/// **The mean is taken after the map, never before it** — and on a nonlinear map the two differ.
#[test]
fn averaging_after_the_map_is_not_the_same_as_averaging_directions() {
    let blend = SiteBlend {
        sites: Sites::Static(StaticGen::Axes6),
        kernel: Kernel::Vmf(6.0),
        colours: vec![
            [0.7, 0.2, 0.0], [0.7, -0.2, 0.0],
            [0.7, 0.0, 0.2], [0.7, 0.0, -0.2],
            [0.7, 0.1, 0.1], [0.7, -0.1, -0.1],
        ],
        space: BlendSpace::Oklab,
    };
    let m = [1.0 / 3.0; 3];
    let dirs = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    let after = compose::resolve(3, |i| blend.blend(dirs[i], &m)).unwrap();
    let mean_dir = {
        let mut v = [0.0f64; 3];
        for d in &dirs {
            for k in 0..3 {
                v[k] += d[k] / 3.0;
            }
        }
        let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        [v[0] / n, v[1] / n, v[2] / n]
    };
    let before = blend.blend(mean_dir, &m).unwrap();
    let d = (0..3).fold(0.0f64, |w, k| w.max((after[k] - before[k]).abs()));
    println!("resolve-after {after:?}\nblend-before  {before:?}\nmax diff {d:.6}");
    assert!(d > 1e-6, "the orders must differ on a nonlinear map, else the caveat is vacuous");
}
