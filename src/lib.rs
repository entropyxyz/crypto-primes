#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg, doc_auto_cfg))]
#![doc = include_str!("../README.md")]
#![deny(unsafe_code)]
#![warn(
    clippy::mod_module_files,
    missing_docs,
    missing_debug_implementations,
    missing_copy_implementations,
    rust_2018_idioms,
    trivial_casts,
    trivial_numeric_casts,
    unused_qualifications
)]

extern crate alloc;

pub mod hazmat;
mod presets;
mod traits;
mod uint_traits;

pub use presets::{
    generate_prime_with_rng, generate_safe_prime_with_rng, is_prime_with_rng,
    is_safe_prime_with_rng,
};
pub use traits::RandomPrimeWithRng;

#[cfg(feature = "default-rng")]
pub use presets::{generate_prime, generate_safe_prime, is_prime, is_safe_prime};

pub use uint_traits::{UintLike, UintModLike};
