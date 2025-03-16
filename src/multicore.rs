//! Prime-finding functions that can parallelize across multiple cores.

use crypto_bigint::{Integer, RandomBits, RandomMod};
use rand_core::CryptoRng;
use rayon::iter::{ParallelBridge, ParallelIterator};

use crate::{
    hazmat::{SetBits, SmallPrimesSieveFactory},
    presets::{is_prime, is_safe_prime},
    traits::SieveFactory,
};

/// Sieves through the results of `sieve_factory` using a thread pool with `threadcount` threads,
/// and returns the first item for which `predicate` is `true`.
///
/// If `sieve_factory` signals that no more results can be created, returns `None`.
pub fn sieve_and_find<R, S, F>(rng: &mut R, sieve_factory: S, predicate: F, threadcount: usize) -> Option<S::Item>
where
    R: CryptoRng + Clone + Send + Sync,
    S: Send + Sync + SieveFactory,
    S::Sieve: Send,
    S::Item: Send,
    F: Sync + Fn(&mut R, &S::Item) -> bool,
{
    let threadpool = rayon::ThreadPoolBuilder::new()
        .num_threads(threadcount)
        .build()
        .expect("If the platform can spawn threads, then this call will work.");

    let mut iter_rng = rng.clone();
    let iter = SieveIterator::new(&mut iter_rng, sieve_factory)?;

    threadpool.install(|| {
        iter.par_bridge().find_any(|c| {
            let mut rng = rng.clone();
            predicate(&mut rng, c)
        })
    })
}

/// A structure that chains the creation of sieves, returning the results from one until it is exhausted,
/// and then creating a new one.
#[derive(Debug)]
struct SieveIterator<'a, R: ?Sized, S: SieveFactory> {
    sieve_factory: S,
    sieve: S::Sieve,
    rng: &'a mut R,
}

impl<'a, R, S> SieveIterator<'a, R, S>
where
    R: CryptoRng + ?Sized,
    S: SieveFactory,
{
    /// Creates a new chained iterator producing results from sieves returned from `sieve_factory`.
    pub fn new(rng: &'a mut R, sieve_factory: S) -> Option<Self> {
        let mut sieve_factory = sieve_factory;
        let sieve = sieve_factory.make_sieve(rng, None)?;
        Some(Self {
            sieve_factory,
            rng,
            sieve,
        })
    }
}

impl<R, S> Iterator for SieveIterator<'_, R, S>
where
    R: CryptoRng + ?Sized,
    S: SieveFactory,
{
    type Item = S::Item;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(result) = self.sieve.next() {
                return Some(result);
            }

            self.sieve = self.sieve_factory.make_sieve(self.rng, Some(&self.sieve))?;
        }
    }
}

/// Returns a random prime of size `bit_length` using the provided RNG.
///
/// Uses `threadcount` cores to parallelize the prime search.
///
/// Panics if `bit_length` is less than 2, or greater than the bit size of the target `Uint`.
///
/// Panics if the platform is unable to spawn threads.
pub fn random_prime<T, R>(rng: &mut R, bit_length: u32, threadcount: usize) -> T
where
    T: Integer + RandomBits + RandomMod,
    R: CryptoRng + Send + Sync + Clone,
{
    sieve_and_find(
        rng,
        SmallPrimesSieveFactory::new(bit_length, SetBits::Msb),
        |_rng, candidate| is_prime(candidate),
        threadcount,
    )
    .expect("will produce a result eventually")
}

/// Returns a random safe prime (that is, such that `(n - 1) / 2` is also prime)
/// of size `bit_length` using the provided RNG.
///
/// Uses `threadcount` cores to parallelize the prime search.
///
/// Panics if `bit_length` is less than 3, or greater than the bit size of the target `Uint`.
/// Panics if the platform is unable to spawn threads.
///
/// See [`is_prime`] for details about the performed checks.
pub fn random_safe_prime<T, R>(rng: &mut R, bit_length: u32, threadcount: usize) -> T
where
    T: Integer + RandomBits + RandomMod,
    R: CryptoRng + Send + Sync + Clone,
{
    sieve_and_find(
        rng,
        SmallPrimesSieveFactory::new_safe_primes(bit_length, SetBits::Msb),
        |_rng, candidate| is_safe_prime(candidate),
        threadcount,
    )
    .expect("will produce a result eventually")
}

#[cfg(test)]
mod tests {
    use crypto_bigint::{nlimbs, BoxedUint, U128};
    use rand_core::{OsRng, TryRngCore};

    use super::{is_prime, random_prime, random_safe_prime};

    #[test]
    fn parallel_prime_generation() {
        for bit_length in (28..=128).step_by(10) {
            let p: U128 = random_prime(&mut OsRng.unwrap_err(), bit_length, 4);
            assert!(p.bits_vartime() == bit_length);
            assert!(is_prime(&p));
        }
    }

    #[test]
    fn parallel_prime_generation_boxed() {
        for bit_length in (28..=128).step_by(10) {
            let p: BoxedUint = random_prime(&mut OsRng.unwrap_err(), bit_length, 2);
            assert!(p.bits_vartime() == bit_length);
            assert!(p.to_words().len() == nlimbs!(bit_length));
            assert!(is_prime(&p));
        }
    }

    #[test]
    fn parallel_safe_prime_generation() {
        for bit_length in (28..=128).step_by(10) {
            let p: U128 = random_safe_prime(&mut OsRng.unwrap_err(), bit_length, 8);
            assert!(p.bits_vartime() == bit_length);
            assert!(is_prime(&p));
        }
    }

    #[test]
    fn parallel_safe_prime_generation_boxed() {
        for bit_length in (28..=128).step_by(10) {
            let p: BoxedUint = random_safe_prime(&mut OsRng.unwrap_err(), bit_length, 4);
            assert!(p.bits_vartime() == bit_length);
            assert!(p.to_words().len() == nlimbs!(bit_length));
            assert!(is_prime(&p));
        }
    }
}
