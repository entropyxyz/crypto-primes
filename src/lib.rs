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

pub use generic::sieve_and_find;
pub use presets::{fips_is_prime, fips_is_safe_prime, is_prime, is_safe_prime, random_prime, random_safe_prime};
pub use traits::SieveFactory;

#[cfg(feature = "multicore")]
pub use presets::{par_random_prime, par_random_safe_prime};

#[cfg(feature = "multicore")]
pub use generic::par_sieve_and_find;
