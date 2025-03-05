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
    unused_qualifications,
    clippy::unwrap_used
)]

extern crate alloc;

mod generic;
pub mod hazmat;
mod presets;
mod traits;

pub use generic::{sieve_and_find, SieveIterator};
pub use presets::{
    fips_is_prime_with_rng, fips_is_safe_prime_with_rng, generate_prime_with_rng, generate_safe_prime_with_rng,
    is_prime, is_safe_prime,
};
pub use traits::{RandomPrimeWithRng, SieveFactory};

#[cfg(feature = "default-rng")]
pub use presets::{generate_prime, generate_safe_prime};
#[cfg(all(feature = "default-rng", feature = "multicore"))]
pub use presets::{par_generate_prime, par_generate_safe_prime};
#[cfg(feature = "multicore")]
pub use presets::{par_generate_prime_with_rng, par_generate_safe_prime_with_rng};

#[cfg(feature = "multicore")]
pub use generic::par_sieve_and_find;
