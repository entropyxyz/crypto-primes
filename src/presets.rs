use crypto_bigint::{Integer, Uint};
use rand_core::{CryptoRng, OsRng, RngCore};

use crate::hazmat::{
    is_strong_lucas_prime, random_odd_uint, sieve_once, MillerRabin, SelfridgeBase, Sieve,
};

/// Returns a random prime of size `bit_length` using [`OsRng`] as the RNG.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
pub fn prime<const L: usize>(bit_length: usize) -> Uint<L> {
    prime_with_rng(&mut OsRng, bit_length)
}

/// Returns a random safe prime (that is, such that `(n - 1) / 2` is also prime)
/// of size `bit_length` using [`OsRng`] as the RNG.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
pub fn safe_prime<const L: usize>(bit_length: usize) -> Uint<L> {
    safe_prime_with_rng(&mut OsRng, bit_length)
}

/// Checks probabilistically if the given number is prime using [`OsRng`] as the RNG.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
pub fn is_prime<const L: usize>(num: &Uint<L>) -> bool {
    is_prime_with_rng(&mut OsRng, num)
}

/// Checks probabilistically if the given number is a safe prime
/// (that is, such that `(n - 1) / 2` is also prime)
/// using [`OsRng`] as the RNG.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
pub fn is_safe_prime<const L: usize>(num: &Uint<L>) -> bool {
    is_safe_prime_with_rng(&mut OsRng, num)
}

/// Returns a random prime of size `bit_length` using the provided RNG.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
pub fn prime_with_rng<const L: usize>(
    rng: &mut (impl CryptoRng + RngCore),
    bit_length: usize,
) -> Uint<L> {
    loop {
        let start = random_odd_uint::<L, _>(rng, bit_length);
        let sieve = Sieve::new(&start, bit_length);
        for num in sieve {
            if _is_prime_with_rng(rng, &num) {
                return num;
            }
        }
    }
}

/// Returns a random safe prime (that is, such that `(n - 1) / 2` is also prime)
/// of size `bit_length` using the provided RNG.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
pub fn safe_prime_with_rng<const L: usize>(
    rng: &mut (impl CryptoRng + RngCore),
    bit_length: usize,
) -> Uint<L> {
    loop {
        let start = random_odd_uint::<L, _>(rng, bit_length);
        let sieve = Sieve::new(&start, bit_length);
        for num in sieve {
            if _is_safe_prime_with_rng(rng, &num) {
                return num;
            }
        }
    }
}

/// Checks probabilistically if the given number is prime using the provided RNG.
///
/// Performed checks:
/// - Trial division by a number of small primes;
/// - Miller-Rabin check with base 2;
/// - Strong Lucas check with Selfridge base (a.k.a. Baillie method A);
/// - Miller-Rabin check with a random base.
///
/// See [`MillerRabin`] and [`is_strong_lucas_prime`] for more details about the checks.
///
/// The second and the third checks constitute the Baillie-PSW primality test[^Baillie1980];
/// the third one is a precaution that follows the approach of GMP (as of v6.2.1).
/// At the moment of the writing there are no known composites passing
/// the Baillie-PSW test[^Baillie2021];
/// it is conjectured that they may exist, but their size is larger than the numbers
/// that are used in practice.
///
/// [^Baillie1980]: R. Baillie, S. S. Wagstaff, "Lucas pseudoprimes",
///       Math. Comp. 35 1391-1417 (1980),
///       DOI: [10.2307/2006406](https://dx.doi.org/10.2307/2006406),
///       <http://mpqs.free.fr/LucasPseudoprimes.pdf>
///
/// [^Baillie2021]: R. Baillie, A. Fiori, S. S. Wagstaff,
///       "Strengthening the Baillie-PSW primality test",
///       Math. Comp. 90 1931-1955 (2021),
///       DOI: [10.1090/mcom/3616](https://doi.org/10.1090/mcom/3616)
pub fn is_prime_with_rng<const L: usize>(
    rng: &mut (impl CryptoRng + RngCore),
    num: &Uint<L>,
) -> bool {
    if let Some(is_prime) = sieve_once(num) {
        return is_prime;
    }
    _is_prime_with_rng(rng, num)
}

/// Checks probabilistically if the given number is prime using the provided RNG.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
pub fn is_safe_prime_with_rng<const L: usize>(
    rng: &mut (impl CryptoRng + RngCore),
    num: &Uint<L>,
) -> bool {
    if let Some(is_prime) = sieve_once(num) {
        return is_prime;
    }
    _is_safe_prime_with_rng(rng, num)
}

/// Checks for safe prime assuming that `num` was already pre-sieved.
fn _is_safe_prime_with_rng<const L: usize>(
    rng: &mut (impl CryptoRng + RngCore),
    num: &Uint<L>,
) -> bool {
    debug_assert!(bool::from(num.is_odd()));
    if !_is_prime_with_rng(rng, num) {
        return false;
    }
    if !is_prime_with_rng(rng, &(num >> 1)) {
        return false;
    }
    true
}

/// Checks for primality assuming that `num` was already pre-sieved.
fn _is_prime_with_rng<const L: usize>(rng: &mut (impl CryptoRng + RngCore), num: &Uint<L>) -> bool {
    debug_assert!(bool::from(num.is_odd()));
    let mr = MillerRabin::new(num);
    if !mr.check_base_two() {
        return false;
    }
    if !is_strong_lucas_prime(num, SelfridgeBase, true) {
        return false;
    }
    if !mr.check_random_base(rng) {
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use crypto_bigint::Uint;
    use rand_core::OsRng;

    use super::is_prime_with_rng;
    use crate::hazmat::primes;

    fn test_large_primes<const L: usize>(nums: &[Uint<L>]) {
        for num in nums {
            assert!(is_prime_with_rng(&mut OsRng, num));
        }
    }

    #[test]
    fn large_primes() {
        test_large_primes(primes::PRIMES_128);
        test_large_primes(primes::PRIMES_256);
        test_large_primes(primes::PRIMES_384);
        test_large_primes(primes::PRIMES_512);
        test_large_primes(primes::PRIMES_1024);
    }
}
