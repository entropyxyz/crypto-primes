//! Miller-Rabin primality test.

use crypto_bigint::{Integer, Limb, Monty, NonZero, Odd, PowBoundedExp, RandomMod, Square};
use rand_core::CryptoRngCore;

use super::Primality;

/// Precomputed data used to perform Miller-Rabin primality test[^Pomerance1980].
///
/// The numbers that pass it are commonly called "strong probable primes"
/// (or "strong pseudoprimes" if they are, in fact, composite).
///
/// [^Pomerance1980]:
///   C. Pomerance, J. L. Selfridge, S. S. Wagstaff "The Pseudoprimes to 25*10^9",
///   Math. Comp. 35 1003-1026 (1980),
///   DOI: [10.2307/2006210](https://dx.doi.org/10.2307/2006210)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MillerRabin<T: Integer> {
    candidate: T,
    bit_length: u32,
    montgomery_params: <<T as Integer>::Monty as Monty>::Params,
    one: <T as Integer>::Monty,
    minus_one: <T as Integer>::Monty,
    s: u32,
    d: T,
}

impl<T: Integer + RandomMod> MillerRabin<T> {
    /// Initializes a Miller-Rabin test for `candidate`.
    pub fn new(candidate: &Odd<T>) -> Self {
        let params = <T as Integer>::Monty::new_params_vartime(candidate.clone());
        let m_one = <T as Integer>::Monty::one(params.clone());
        let m_minus_one = -m_one.clone();

        let one = T::one_like(candidate.as_ref());

        // Find `s` and odd `d` such that `candidate - 1 == 2^s * d`.
        let (s, d) = if candidate.as_ref() == &one {
            (0, one)
        } else {
            let candidate_minus_one = candidate.wrapping_sub(&one);
            let s = candidate_minus_one.trailing_zeros_vartime();
            // Will not overflow because `candidate` is odd and greater than 1.
            let d = candidate_minus_one
                .overflowing_shr_vartime(s)
                .expect("shift should be within range by construction");
            (s, d)
        };

        Self {
            candidate: candidate.as_ref().clone(),
            bit_length: candidate.bits_vartime(),
            montgomery_params: params,
            one: m_one,
            minus_one: m_minus_one,
            s,
            d,
        }
    }

    /// Perform a Miller-Rabin check with a given base.
    pub fn test(&self, base: &T) -> Primality {
        // TODO: it may be faster to first check that gcd(base, candidate) == 1,
        // otherwise we can return `Composite` right away.

        let base = <T as Integer>::Monty::new(base.clone(), self.montgomery_params.clone());

        // Implementation detail: bounded exp gets faster every time we decrease the bound
        // by the window length it uses, which is currently 4 bits.
        // So even when the bound isn't low enough that the number can fit
        // in a smaller number of limbs, there is still a performance gain
        // from specifying the bound.
        let mut test = base.pow_bounded_exp(&self.d, self.bit_length);

        if test == self.one || test == self.minus_one {
            return Primality::ProbablyPrime;
        }
        for _ in 1..self.s {
            test = test.square();
            if test == self.one {
                return Primality::Composite;
            } else if test == self.minus_one {
                return Primality::ProbablyPrime;
            }
        }
        Primality::Composite
    }

    /// Perform a Miller-Rabin check with base 2.
    pub fn test_base_two(&self) -> Primality {
        self.test(&T::from_limb_like(Limb::from(2u32), &self.candidate))
    }

    /// Perform a Miller-Rabin check with a random base (in the range `[3, candidate-2]`)
    /// drawn using the provided RNG.
    ///
    /// Note: panics if `candidate == 3` (so the range above is empty).
    pub fn test_random_base(&self, rng: &mut impl CryptoRngCore) -> Primality {
        // We sample a random base from the range `[3, candidate-2]`:
        // - we have a separate method for base 2;
        // - the test holds trivially for bases 1 or `candidate-1`.
        if self.candidate.bits() < 3 {
            panic!("No suitable random base possible when `candidate == 3`; use the base 2 test.")
        }

        let range = self.candidate.wrapping_sub(&T::from(4u32));
        // Can unwrap here since `candidate` is odd, and `candidate >= 4` (as checked above)
        let range_nonzero =
            NonZero::new(range).expect("the range should be non-zero by construction");
        // This should not overflow as long as `random_mod()` behaves according to the contract
        // (that is, returns a number within the given range).
        let random = T::random_mod(rng, &range_nonzero)
            .checked_add(&T::from(3u32))
            .expect("addition should not overflow by construction");
        self.test(&random)
    }
}

#[cfg(test)]
mod tests {

    use alloc::format;
    use core::num::NonZeroU32;

    use crypto_bigint::{Integer, Odd, RandomMod, Uint, U1024, U128, U1536, U256, U64};
    use rand_chacha::ChaCha8Rng;
    use rand_core::{CryptoRngCore, OsRng, SeedableRng};

    #[cfg(feature = "tests-exhaustive")]
    use num_prime::nt_funcs::is_prime64;

    use super::MillerRabin;
    use crate::hazmat::{primes, pseudoprimes, random_odd_integer, Sieve};

    #[test]
    fn miller_rabin_derived_traits() {
        let mr = MillerRabin::new(&Odd::new(U64::ONE).unwrap());
        assert!(format!("{mr:?}").starts_with("MillerRabin"));
        assert_eq!(mr.clone(), mr);
    }

    #[test]
    #[should_panic(
        expected = "No suitable random base possible when `candidate == 3`; use the base 2 test."
    )]
    fn random_base_range_check() {
        let mr = MillerRabin::new(&Odd::new(U64::from(3u32)).unwrap());
        mr.test_random_base(&mut OsRng);
    }

    fn is_spsp(num: u32) -> bool {
        pseudoprimes::STRONG_BASE_2.iter().any(|x| *x == num)
    }

    fn random_checks<T: Integer + RandomMod>(
        rng: &mut impl CryptoRngCore,
        mr: &MillerRabin<T>,
        count: usize,
    ) -> usize {
        (0..count)
            .map(|_| -> usize { mr.test_random_base(rng).is_probably_prime().into() })
            .sum()
    }

    fn test_composites(numbers: &[u32], expected_result: bool) {
        let mut rng = ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");
        for num in numbers.iter() {
            let base_2_false_positive = is_spsp(*num);
            let actual_expected_result = if base_2_false_positive {
                true
            } else {
                expected_result
            };

            // A random base MR test is expected to report a composite as a prime
            // with about 1/4 probability. So we're expecting less than
            // 35 out of 100 false positives, seems to work.

            let mr = MillerRabin::new(&Odd::new(U64::from(*num)).unwrap());
            assert_eq!(
                mr.test_base_two().is_probably_prime(),
                actual_expected_result
            );
            let reported_prime = random_checks(&mut rng, &mr, 100);
            assert!(
                reported_prime < 35,
                "reported as prime in {reported_prime} out of 100 tests",
            );
        }
    }

    #[test]
    fn trivial() {
        let mut rng = ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");
        let start =
            random_odd_integer::<U1024>(&mut rng, NonZeroU32::new(1024).unwrap(), U1024::BITS);
        for num in Sieve::new(start.as_ref(), NonZeroU32::new(1024).unwrap(), false).take(10) {
            let mr = MillerRabin::new(&Odd::new(num).unwrap());

            // Trivial tests, must always be true.
            assert!(mr.test(&1u32.into()).is_probably_prime());
            assert!(mr.test(&num.wrapping_sub(&1u32.into())).is_probably_prime());
        }
    }

    #[test]
    fn mersenne_prime() {
        let mut rng = ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");

        // Mersenne prime 2^127-1
        let num = Odd::new(U128::from_be_hex("7fffffffffffffffffffffffffffffff")).unwrap();

        let mr = MillerRabin::new(&num);
        assert!(mr.test_base_two().is_probably_prime());
        for _ in 0..10 {
            assert!(mr.test_random_base(&mut rng).is_probably_prime());
        }
    }

    #[test]
    fn strong_fibonacci_pseudoprimes() {
        let mut rng = ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");

        for num in pseudoprimes::STRONG_FIBONACCI.iter() {
            let mr = MillerRabin::new(&Odd::new(*num).unwrap());
            assert!(!mr.test_base_two().is_probably_prime());
            for _ in 0..1000 {
                assert!(!mr.test_random_base(&mut rng).is_probably_prime());
            }
        }
    }

    #[test]
    fn strong_pseudoprimes_base_2() {
        // These known exceptions for the base 2 MR test.
        test_composites(pseudoprimes::STRONG_BASE_2, true);
    }

    #[test]
    fn lucas_pseudoprimes() {
        // Cross-test against the pseudoprimes that circumvent the Lucas test.
        // We expect the MR test to correctly classify them as composites (most of the time).
        test_composites(pseudoprimes::EXTRA_STRONG_LUCAS, false);
        test_composites(pseudoprimes::STRONG_LUCAS, false);
        test_composites(pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS, false);
        test_composites(pseudoprimes::FIBONACCI, false);
        test_composites(pseudoprimes::BRUCKMAN_LUCAS, false);
        test_composites(pseudoprimes::LUCAS, false);
    }

    #[test]
    fn large_carmichael_number() {
        let mr = MillerRabin::new(&Odd::new(pseudoprimes::LARGE_CARMICHAEL_NUMBER).unwrap());

        // It is known to pass MR tests for all prime bases <307
        assert!(mr.test_base_two().is_probably_prime());
        assert!(mr.test(&U1536::from(293u64)).is_probably_prime());

        // A test with base 307 correctly reports the number as composite.
        assert!(!mr.test(&U1536::from(307u64)).is_probably_prime());
    }

    fn test_large_primes<const L: usize>(nums: &[Uint<L>]) {
        let mut rng = ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");
        for num in nums {
            let mr = MillerRabin::new(&Odd::new(*num).unwrap());
            assert!(mr.test_base_two().is_probably_prime());
            for _ in 0..10 {
                assert!(mr.test_random_base(&mut rng).is_probably_prime());
            }
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

    #[cfg(feature = "tests-exhaustive")]
    #[test]
    fn exhaustive() {
        // Test all the odd numbers up to the limit where we know the false positives,
        // and compare the results with the reference.
        for num in (3..pseudoprimes::EXHAUSTIVE_TEST_LIMIT).step_by(2) {
            let res_ref = is_prime64(num.into());

            let spsp = is_spsp(num);

            let mr = MillerRabin::new(&Odd::new(U64::from(num)).unwrap());
            let res = mr.test_base_two().is_probably_prime();
            let expected = spsp || res_ref;
            assert_eq!(
                res, expected,
                "Miller-Rabin: n={num}, expected={expected}, actual={res}",
            );
        }
    }
}
