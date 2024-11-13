use rand_core::CryptoRngCore;

#[cfg(feature = "multicore")]
use rayon::iter::{ParallelBridge, ParallelIterator};

use crate::SieveFactory;

/// Sieves through the results of `sieve_factory` and returns the first item for which `predicate` is `true`.
///
/// If `sieve_factory` signals that no more results can be created, returns `None`.
pub fn sieve_and_find<R, S, T>(rng: &mut R, sieve_factory: S, predicate: impl Fn(&mut R, &T) -> bool) -> Option<T>
where
    S: SieveFactory<T>,
    R: CryptoRngCore,
{
    // We could use `SieveIterator` here, but it requires cloning the `rng`.
    // Unlike the parallel version, it is avoidable here.

    let mut sieve = sieve_factory.make_sieve(rng, None)?;

    loop {
        if let Some(value) = sieve.find(|num| predicate(rng, num)) {
            return Some(value);
        }
        if let Some(new_sieve) = sieve_factory.make_sieve(rng, Some(&sieve)) {
            sieve = new_sieve;
        } else {
            return None;
        }
    }
}

/// Sieves through the results of `sieve_factory` using a thread pool with `threadcount` threads,
/// and returns the first item for which `predicate` is `true`.
///
/// If `sieve_factory` signals that no more results can be created, returns `None`.
#[cfg(feature = "multicore")]
pub fn par_sieve_and_find<R, S, T, F>(rng: &mut R, sieve_factory: S, predicate: F, threadcount: usize) -> Option<T>
where
    R: CryptoRngCore + Clone + Send + Sync,
    T: Send,
    S: Send + Sync + SieveFactory<T>,
    S::Sieve: Send,
    F: Sync + Fn(&mut R, &T) -> bool,
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
pub struct SieveIterator<'a, R: CryptoRngCore, T, S: SieveFactory<T>> {
    sieve_factory: S,
    sieve: S::Sieve,
    rng: &'a mut R,
}

impl<'a, R: CryptoRngCore, T, S: SieveFactory<T>> SieveIterator<'a, R, T, S> {
    /// Creates a new chained iterator producing results from sieves returned from `sieve_factory`.
    pub fn new(rng: &'a mut R, sieve_factory: S) -> Option<Self> {
        let sieve = sieve_factory.make_sieve(rng, None)?;
        Some(Self {
            sieve_factory,
            rng,
            sieve,
        })
    }
}

impl<'a, R: CryptoRngCore, T, S: SieveFactory<T>> Iterator for SieveIterator<'a, R, T, S> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(result) = self.sieve.next() {
                return Some(result);
            }

            self.sieve = self.sieve_factory.make_sieve(self.rng, Some(&self.sieve))?;
        }
    }
}
