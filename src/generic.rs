use rand_core::CryptoRng;

use crate::SieveFactory;

/// Sieves through the results of `sieve_factory` and returns the first item for which `predicate` is `true`.
///
/// If `sieve_factory` signals that no more results can be created, returns `None`.
pub fn sieve_and_find<R, S>(
    rng: &mut R,
    sieve_factory: S,
    predicate: impl Fn(&mut R, &S::Item) -> bool,
) -> Option<S::Item>
where
    S: SieveFactory,
    R: CryptoRng + ?Sized,
{
    // We could use `SieveIterator` here, but it requires cloning the `rng`.
    // Unlike the parallel version, it is avoidable here.

    let mut sieve_factory = sieve_factory;
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

#[cfg(test)]
mod tests {
    use rand_core::{CryptoRng, OsRng, TryRngCore};

    use super::sieve_and_find;
    use crate::SieveFactory;

    #[test]
    fn test_exhaustable_sieve_factory() {
        // Test the logic handling the case of the sieve factory not being able to produce new sieves.
        struct TestSieveFactory {
            count: usize,
        }

        impl SieveFactory for TestSieveFactory {
            type Item = usize;
            type Sieve = core::ops::Range<usize>;

            fn make_sieve<R: CryptoRng + ?Sized>(
                &mut self,
                _rng: &mut R,
                previous_sieve: Option<&Self::Sieve>,
            ) -> Option<Self::Sieve> {
                self.count += 1;
                if previous_sieve.is_none() {
                    Some(self.count * 10..(self.count * 10 + 2))
                } else {
                    None
                }
            }
        }

        let factory = TestSieveFactory { count: 0 };
        let result = sieve_and_find(&mut OsRng.unwrap_err(), factory, |_rng, num| *num == 11);
        assert!(result.is_some());

        let factory = TestSieveFactory { count: 0 };
        let result = sieve_and_find(&mut OsRng.unwrap_err(), factory, |_rng, num| *num == 20);
        assert!(result.is_none());
    }
}
