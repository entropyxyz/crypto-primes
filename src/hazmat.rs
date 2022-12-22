//! Components to build your own primality test.
//! Handle with care.

mod jacobi;
mod lucas;
mod miller_rabin;
mod precomputed;
#[cfg(test)]
pub(crate) mod primes;
#[cfg(test)]
pub(crate) mod pseudoprimes;
mod sieve;

pub use lucas::{is_strong_lucas_prime, BruteForceBase, LucasBase, SelfridgeBase};
pub use miller_rabin::MillerRabin;
pub use sieve::{random_odd_uint, sieve_once, Sieve};
