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

pub use presets::{is_prime_with_rng, is_safe_prime_with_rng, prime_with_rng, safe_prime_with_rng};

#[cfg(feature = "default-rng")]
pub use presets::{is_prime, is_safe_prime, prime, safe_prime};
