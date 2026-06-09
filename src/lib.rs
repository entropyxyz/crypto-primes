#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg))]
#![doc = include_str!("../README.md")]

#[cfg(any(feature = "alloc", test))]
extern crate alloc;

mod error;
pub mod fips;
mod generic;
pub mod hazmat;
mod presets;

#[cfg(feature = "multicore")]
pub mod multicore;

pub use error::Error;
pub use generic::sieve_and_find;
pub use presets::{Flavor, is_prime, random_prime};
