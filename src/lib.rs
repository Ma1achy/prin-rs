//! prin-rs — uniform-resolution three-body initial-condition kernel.
//!
//! `BRIEF.md` is the spec; `CLAUDE.md` is the working agreement. The physics is the
//! product — the image is a diagnostic.

pub mod real;
pub mod vec2;
pub mod config;
pub mod ensemble;
pub mod grid;
pub mod integrate;
pub mod outcome;
pub mod physics;
pub mod output;
pub mod render;
pub mod rng;

pub use real::Real;
pub use vec2::Vec2;
