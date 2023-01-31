use crypto_bigint::Uint;
use rand_core::CryptoRngCore;

use crate::{is_prime_with_rng, is_safe_prime_with_rng, prime_with_rng, safe_prime_with_rng};

/// Provides a generic way to access methods for random prime number generation
/// and primality checking, wrapping the standalone functions ([`is_prime_with_rng`] etc).
pub trait RandomPrimeWithRng {
    /// Returns a random prime of size `bit_length` using the provided RNG.
    ///
    /// Panics if `bit_length` is less than 2, or greater than the bit size of the target `Uint`.
    ///
    /// See [`is_prime_with_rng`] for details about the performed checks.
    fn prime_with_rng(rng: &mut impl CryptoRngCore, bit_length: usize) -> Self;

    /// Returns a random safe prime (that is, such that `(n - 1) / 2` is also prime)
    /// of size `bit_length` using the provided RNG.
    ///
    /// Panics if `bit_length` is less than 3, or greater than the bit size of the target `Uint`.
    ///
    /// See [`is_prime_with_rng`] for details about the performed checks.
    fn safe_prime_with_rng(rng: &mut impl CryptoRngCore, bit_length: usize) -> Self;

    /// Checks probabilistically if the given number is prime using the provided RNG.
    ///
    /// See [`is_prime_with_rng`] for details about the performed checks.
    fn is_prime_with_rng(rng: &mut impl CryptoRngCore, num: &Self) -> bool;

    /// Checks probabilistically if the given number is a safe prime using the provided RNG.
    ///
    /// See [`is_prime_with_rng`] for details about the performed checks.
    fn is_safe_prime_with_rng(rng: &mut impl CryptoRngCore, num: &Self) -> bool;
}

impl<const L: usize> RandomPrimeWithRng for Uint<L> {
    fn prime_with_rng(rng: &mut impl CryptoRngCore, bit_length: usize) -> Self {
        prime_with_rng(rng, bit_length)
    }
    fn safe_prime_with_rng(rng: &mut impl CryptoRngCore, bit_length: usize) -> Self {
        safe_prime_with_rng(rng, bit_length)
    }
    fn is_prime_with_rng(rng: &mut impl CryptoRngCore, num: &Self) -> bool {
        is_prime_with_rng(rng, num)
    }
    fn is_safe_prime_with_rng(rng: &mut impl CryptoRngCore, num: &Self) -> bool {
        is_safe_prime_with_rng(rng, num)
    }
}
