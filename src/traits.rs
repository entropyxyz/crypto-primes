use crypto_bigint::{Integer, RandomBits, RandomMod};
use rand_core::CryptoRng;

use crate::{generate_prime_with_rng, generate_safe_prime_with_rng, is_prime, is_safe_prime};

/// A type producing sieves for random prime generation.
pub trait SieveFactory {
    /// The type of items returning by the sieves.
    type Item;

    /// The resulting sieve.
    type Sieve: Iterator<Item = Self::Item>;

    /// Makes a sieve given an RNG and the previous exhausted sieve (if any).
    ///
    /// Returning `None` signals that the prime generation should stop.
    fn make_sieve<R: CryptoRng + ?Sized>(
        &mut self,
        rng: &mut R,
        previous_sieve: Option<&Self::Sieve>,
    ) -> Option<Self::Sieve>;
}

/// Provides a generic way to access methods for random prime number generation
/// and primality checking, wrapping the standalone functions ([`is_prime_with_rng`] etc).
pub trait RandomPrimeWithRng {
    /// Returns a random prime of size `bit_length` using the provided RNG.
    ///
    /// Panics if `bit_length` is less than 2, or greater than the bit size of the target `Uint`.
    ///
    /// See [`is_prime_with_rng`] for details about the performed checks.
    fn generate_prime_with_rng<R: CryptoRng + ?Sized>(rng: &mut R, bit_length: u32) -> Self;

    /// Returns a random safe prime (that is, such that `(n - 1) / 2` is also prime)
    /// of size `bit_length` using the provided RNG.
    ///
    /// Panics if `bit_length` is less than 3, or greater than the bit size of the target `Uint`.
    ///
    /// See [`is_prime_with_rng`] for details about the performed checks.
    fn generate_safe_prime_with_rng<R: CryptoRng + ?Sized>(rng: &mut R, bit_length: u32) -> Self;

    /// Probabilistically checks if the given number is prime using the provided RNG.
    ///
    /// See [`is_prime_with_rng`] for details about the performed checks.
    fn is_prime_with_rng<R: CryptoRng + ?Sized>(&self, rng: &mut R) -> bool;

    /// Probabilistically checks if the given number is a safe prime using the provided RNG.
    ///
    /// See [`is_prime_with_rng`] for details about the performed checks.
    fn is_safe_prime_with_rng<R: CryptoRng + ?Sized>(&self, rng: &mut R) -> bool;
}

impl<T> RandomPrimeWithRng for T
where
    T: Integer + RandomBits + RandomMod,
{
    fn generate_prime_with_rng<R: CryptoRng + ?Sized>(rng: &mut R, bit_length: u32) -> Self {
        generate_prime_with_rng(rng, bit_length)
    }
    fn generate_safe_prime_with_rng<R: CryptoRng + ?Sized>(rng: &mut R, bit_length: u32) -> Self {
        generate_safe_prime_with_rng(rng, bit_length)
    }
    fn is_prime_with_rng<R: CryptoRng + ?Sized>(&self, _rng: &mut R) -> bool {
        is_prime(self)
    }
    fn is_safe_prime_with_rng<R: CryptoRng + ?Sized>(&self, _rng: &mut R) -> bool {
        is_safe_prime(self)
    }
}

#[cfg(test)]
mod tests {
    use crypto_bigint::{BoxedUint, U64};
    use rand_core::{OsRng, TryRngCore};

    use super::RandomPrimeWithRng;

    #[test]
    fn uint_impl() {
        assert!(!U64::from(15u32).is_prime_with_rng(&mut OsRng.unwrap_err()));
        assert!(U64::from(19u32).is_prime_with_rng(&mut OsRng.unwrap_err()));

        assert!(!U64::from(13u32).is_safe_prime_with_rng(&mut OsRng.unwrap_err()));
        assert!(U64::from(11u32).is_safe_prime_with_rng(&mut OsRng.unwrap_err()));

        assert!(U64::generate_prime_with_rng(&mut OsRng.unwrap_err(), 10).is_prime_with_rng(&mut OsRng.unwrap_err()));
        assert!(U64::generate_safe_prime_with_rng(&mut OsRng.unwrap_err(), 10)
            .is_safe_prime_with_rng(&mut OsRng.unwrap_err()));
    }

    #[test]
    fn boxed_uint_impl() {
        assert!(!BoxedUint::from(15u32).is_prime_with_rng(&mut OsRng.unwrap_err()));
        assert!(BoxedUint::from(19u32).is_prime_with_rng(&mut OsRng.unwrap_err()));

        assert!(!BoxedUint::from(13u32).is_safe_prime_with_rng(&mut OsRng.unwrap_err()));
        assert!(BoxedUint::from(11u32).is_safe_prime_with_rng(&mut OsRng.unwrap_err()));

        assert!(
            BoxedUint::generate_prime_with_rng(&mut OsRng.unwrap_err(), 10).is_prime_with_rng(&mut OsRng.unwrap_err())
        );
        assert!(BoxedUint::generate_safe_prime_with_rng(&mut OsRng.unwrap_err(), 10)
            .is_safe_prime_with_rng(&mut OsRng.unwrap_err()));
    }
}
