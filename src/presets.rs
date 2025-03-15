use crypto_bigint::{Integer, Odd, RandomBits, RandomMod, Word};
use rand_core::CryptoRng;

#[cfg(feature = "default-rng")]
use rand_core::{OsRng, TryRngCore};

use crate::{
    generic::sieve_and_find,
    hazmat::{
        lucas_test, AStarBase, LucasCheck, MillerRabin, Primality, SelfridgeBase, SetBits, SmallPrimesSieveFactory,
    },
};

#[cfg(feature = "multicore")]
use crate::generic::par_sieve_and_find;

/// Returns a random prime of size `bit_length` using [`OsRng`] as the RNG.
///
/// See [`is_prime`] for details about the performed checks.
#[cfg(feature = "default-rng")]
pub fn generate_prime<T: Integer + RandomBits + RandomMod>(bit_length: u32) -> T {
    generate_prime_with_rng(&mut OsRng.unwrap_err(), bit_length)
}

/// Returns a random prime of size `bit_length` using [`OsRng`] as the RNG.
///
/// See [`is_prime`] for details about the performed checks.
///
/// Uses `threadcount` cores to parallelize the prime search.
///
/// Panics if `bit_length` is less than 2, or greater than the bit size of the target `Uint`.
///
/// Panics if the platform is unable to spawn threads.
#[cfg(all(feature = "default-rng", feature = "multicore"))]
pub fn par_generate_prime<T: Integer + RandomBits + RandomMod>(bit_length: u32, threadcount: usize) -> T {
    par_generate_prime_with_rng(&mut OsRng.unwrap_err(), bit_length, threadcount)
}

/// Returns a random safe prime (that is, such that `(n - 1) / 2` is also prime) of size
/// `bit_length` using [`OsRng`] as the RNG.
///
/// See [`is_prime`] for details about the performed checks.
#[cfg(feature = "default-rng")]
pub fn generate_safe_prime<T: Integer + RandomBits + RandomMod>(bit_length: u32) -> T {
    generate_safe_prime_with_rng(&mut OsRng.unwrap_err(), bit_length)
}

/// Returns a random safe prime (that is, such that `(n - 1) / 2` is also prime) of size
/// `bit_length` using [`OsRng`] as the RNG.
///
/// See [`is_prime`] for details about the performed checks.
///
/// Uses `threadcount` cores to parallelize the prime search.
///
/// Panics if `bit_length` is less than 3, or greater than the bit size of the target `Uint`.
///
/// Panics if the platform is unable to spawn threads.
#[cfg(all(feature = "default-rng", feature = "multicore"))]
pub fn par_generate_safe_prime<T: Integer + RandomBits + RandomMod>(bit_length: u32, threadcount: usize) -> T {
    par_generate_safe_prime_with_rng(&mut OsRng.unwrap_err(), bit_length, threadcount)
}

/// Returns a random prime of size `bit_length` using the provided RNG.
///
/// Panics if `bit_length` is less than 2, or greater than the bit size of the target `Uint`.
///
/// See [`is_prime`] for details about the performed checks.
pub fn generate_prime_with_rng<T: Integer + RandomBits + RandomMod, R: CryptoRng + ?Sized>(
    rng: &mut R,
    bit_length: u32,
) -> T {
    sieve_and_find(
        rng,
        SmallPrimesSieveFactory::new(bit_length, SetBits::Msb),
        |_rng, candidate| is_prime(candidate),
    )
    .expect("will produce a result eventually")
}

/// Returns a random safe prime (that is, such that `(n - 1) / 2` is also prime)
/// of size `bit_length` using the provided RNG.
///
/// Panics if `bit_length` is less than 3, or greater than the bit size of the target `Uint`.
///
/// See [`is_prime`] for details about the performed checks.
pub fn generate_safe_prime_with_rng<T: Integer + RandomBits + RandomMod, R: CryptoRng + ?Sized>(
    rng: &mut R,
    bit_length: u32,
) -> T {
    sieve_and_find(
        rng,
        SmallPrimesSieveFactory::new_safe_primes(bit_length, SetBits::Msb),
        |_rng, candidate| is_safe_prime(candidate),
    )
    .expect("will produce a result eventually")
}

/// Returns a random prime of size `bit_length` using the provided RNG.
///
/// Uses `threadcount` cores to parallelize the prime search.
///
/// Panics if `bit_length` is less than 2, or greater than the bit size of the target `Uint`.
///
/// Panics if the platform is unable to spawn threads.
#[cfg(feature = "multicore")]
pub fn par_generate_prime_with_rng<T, R>(rng: &mut R, bit_length: u32, threadcount: usize) -> T
where
    T: Integer + RandomBits + RandomMod,
    R: CryptoRng + Send + Sync + Clone,
{
    par_sieve_and_find(
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
#[cfg(feature = "multicore")]
pub fn par_generate_safe_prime_with_rng<T, R>(rng: &mut R, bit_length: u32, threadcount: usize) -> T
where
    T: Integer + RandomBits + RandomMod,
    R: CryptoRng + Send + Sync + Clone,
{
    par_sieve_and_find(
        rng,
        SmallPrimesSieveFactory::new_safe_primes(bit_length, SetBits::Msb),
        |_rng, candidate| is_safe_prime(candidate),
        threadcount,
    )
    .expect("will produce a result eventually")
}

fn equals_primitive<T: Integer>(num: &T, primitive: Word) -> bool {
    num.bits_vartime() <= u16::BITS && num.as_ref()[0].0 == primitive
}

/// Checks if the given number is prime.
///
/// Performed tests:
/// - Miller-Rabin test with base 2;
/// - [`LucasCheck::Bpsw21`] test with [`AStarBase`].
///
/// See [`MillerRabin`] and [`lucas_test`] for more details about the tests.
///
/// This is the recommended approach by Baillie et al[^Baillie2021],
/// improving on the BPSW'80 test[^Baillie1980].
/// At the moment of the writing there are no known composites that are reported to be primes
/// by either of the approaches[^Baillie2021];
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
pub fn is_prime<T: Integer + RandomMod>(candidate: &T) -> bool {
    if equals_primitive(candidate, 1) {
        return false;
    }

    if equals_primitive(candidate, 2) {
        return true;
    }

    let odd_candidate: Odd<T> = match Odd::new(candidate.clone()).into() {
        Some(x) => x,
        None => return false,
    };

    let mr = MillerRabin::new(odd_candidate.clone());

    if !mr.test_base_two().is_probably_prime() {
        return false;
    }

    match lucas_test(odd_candidate, AStarBase, LucasCheck::Bpsw21) {
        Primality::Composite => false,
        Primality::Prime => true,
        Primality::ProbablyPrime => true,
    }
}

/// Checks if the given number is a safe prime.
///
/// See [`is_prime`] for details about the performed checks.
pub fn is_safe_prime<T: Integer + RandomMod>(candidate: &T) -> bool {
    // Since, by the definition of safe prime, `(candidate - 1) / 2` must also be prime,
    // and therefore odd, `candidate` has to be equal to 3 modulo 4.
    // 5 is the only exception, so we check for it.
    if equals_primitive(candidate, 5) {
        return true;
    }

    // Safe primes are always of the form 4k + 3 (i.e. n ≡ 3 mod 4)
    // The last two digits of a binary number give you its value modulo 4.
    // Primes p=4n+3 will always end in 11 in binary because p ≡ 3 mod 4.
    if candidate.as_ref()[0].0 & 3 != 3 {
        return false;
    }

    is_prime(candidate) && is_prime(&candidate.wrapping_shr_vartime(1))
}

/// Probabilistically checks if the given number is prime using the provided RNG
/// according to FIPS-186.5[^FIPS] standard.
///
/// Performed checks:
/// - `mr_iterations` of Miller-Rabin check with random bases;
/// - Regular Lucas check with Selfridge base (see [`SelfridgeBase`] for details), if `add_lucas_test` is `true`.
///
/// See [`MillerRabin`] and [`lucas_test`] for more details about the checks;
/// use [`minimum_mr_iterations`](`crate::hazmat::minimum_mr_iterations`)
/// to calculate the number of required iterations.
///
/// [^FIPS]: FIPS-186.5 standard, <https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf>
pub fn fips_is_prime_with_rng<T: Integer + RandomMod>(
    rng: &mut (impl CryptoRng + ?Sized),
    candidate: &T,
    mr_iterations: usize,
    add_lucas_test: bool,
) -> bool {
    if equals_primitive(candidate, 1) {
        return false;
    }

    if equals_primitive(candidate, 2) {
        return true;
    }

    let odd_candidate: Odd<T> = match Odd::new(candidate.clone()).into() {
        Some(x) => x,
        None => return false,
    };

    // The random base test only makes sense when `candidate > 3`.
    if !equals_primitive(candidate, 3) {
        let mr = MillerRabin::new(odd_candidate.clone());
        for _ in 0..mr_iterations {
            if !mr.test_random_base(rng).is_probably_prime() {
                return false;
            }
        }
    }

    if add_lucas_test {
        match lucas_test(odd_candidate, SelfridgeBase, LucasCheck::Strong) {
            Primality::Composite => return false,
            Primality::Prime => return true,
            _ => {}
        }
    }

    true
}

/// Probabilistically checks if the given number is a safe prime using the provided RNG
/// according to FIPS-186.5[^FIPS] standard.
///
/// See [`fips_is_prime_with_rng`] for details about the performed checks.
///
/// [^FIPS]: FIPS-186.5 standard, <https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf>
pub fn fips_is_safe_prime_with_rng<T: Integer + RandomMod>(
    rng: &mut (impl CryptoRng + ?Sized),
    candidate: &T,
    mr_iterations: usize,
    add_lucas_test: bool,
) -> bool {
    // Since, by the definition of safe prime, `(candidate - 1) / 2` must also be prime,
    // and therefore odd, `candidate` has to be equal to 3 modulo 4.
    // 5 is the only exception, so we check for it.
    if equals_primitive(candidate, 5) {
        return true;
    }

    // Safe primes are always of the form 4k + 3 (i.e. n ≡ 3 mod 4)
    // The last two digits of a binary number give you its value modulo 4.
    // Primes p=4n+3 will always end in 11 in binary because p ≡ 3 mod 4.
    if candidate.as_ref()[0].0 & 3 != 3 {
        return false;
    }

    fips_is_prime_with_rng(rng, candidate, mr_iterations, add_lucas_test)
        && fips_is_prime_with_rng(rng, &candidate.wrapping_shr_vartime(1), mr_iterations, add_lucas_test)
}

#[cfg(test)]
mod tests {
    use crypto_bigint::{nlimbs, BoxedUint, CheckedAdd, Integer, RandomMod, Uint, Word, U128, U64};
    use num_prime::nt_funcs::is_prime64;
    use rand_core::{OsRng, TryRngCore};

    use super::{
        fips_is_prime_with_rng, fips_is_safe_prime_with_rng, generate_prime, generate_prime_with_rng,
        generate_safe_prime, generate_safe_prime_with_rng, is_prime, is_safe_prime,
    };
    use crate::hazmat::{minimum_mr_iterations, primes, pseudoprimes};

    fn fips_is_prime<T: Integer + RandomMod>(num: &T) -> bool {
        let mr_iterations = minimum_mr_iterations(128, 100).unwrap();
        fips_is_prime_with_rng(&mut OsRng.unwrap_err(), num, mr_iterations, false)
    }

    fn fips_is_safe_prime<T: Integer + RandomMod>(num: &T) -> bool {
        let mr_iterations = minimum_mr_iterations(128, 100).unwrap();
        fips_is_safe_prime_with_rng(&mut OsRng.unwrap_err(), num, mr_iterations, false)
    }

    fn test_large_primes<const L: usize>(nums: &[Uint<L>]) {
        for num in nums {
            assert!(is_prime(num));
            assert!(fips_is_prime(num));
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

    fn test_pseudoprimes(nums: &[u32]) {
        for num in nums {
            assert!(!is_prime(&U64::from(*num)));
            assert!(!fips_is_prime(&U64::from(*num)));
        }
    }

    #[test]
    fn pseudoprimes() {
        test_pseudoprimes(pseudoprimes::EXTRA_STRONG_LUCAS);
        test_pseudoprimes(pseudoprimes::STRONG_LUCAS);
        test_pseudoprimes(pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS);
        test_pseudoprimes(pseudoprimes::STRONG_BASE_2);
        test_pseudoprimes(pseudoprimes::FIBONACCI);
        test_pseudoprimes(pseudoprimes::BRUCKMAN_LUCAS);
        test_pseudoprimes(pseudoprimes::LUCAS);

        for num in pseudoprimes::STRONG_FIBONACCI {
            assert!(!is_prime(num));
            assert!(!fips_is_prime(num));
        }

        assert!(!is_prime(&pseudoprimes::LARGE_CARMICHAEL_NUMBER));
        assert!(!fips_is_prime(&pseudoprimes::LARGE_CARMICHAEL_NUMBER));
    }

    fn test_cunningham_chain<const L: usize>(length: usize, num: &Uint<L>) {
        let mut next = *num;
        for i in 0..length {
            assert!(is_prime(&next));
            assert!(fips_is_prime(&next));

            // The start of the chain isn't a safe prime by definition
            if i > 0 {
                assert!(is_safe_prime(&next));
                assert!(fips_is_safe_prime(&next));
            }

            next = next.wrapping_shl_vartime(1).checked_add(&Uint::<L>::ONE).unwrap();
        }

        // The chain ended.
        assert!(!is_prime(&next));
        assert!(!fips_is_prime(&next));
    }

    #[test]
    fn cunningham_chains() {
        for (length, num) in primes::CUNNINGHAM_CHAINS_128 {
            test_cunningham_chain(*length, num);
        }
        for (length, num) in primes::CUNNINGHAM_CHAINS_512 {
            test_cunningham_chain(*length, num);
        }
    }

    #[test]
    fn prime_generation() {
        for bit_length in (28..=128).step_by(10) {
            let p: U128 = generate_prime(bit_length);
            assert!(p.bits_vartime() == bit_length);
            assert!(is_prime(&p));
            assert!(fips_is_prime(&p));
        }
    }

    #[test]
    fn prime_generation_boxed() {
        for bit_length in (28..=128).step_by(10) {
            let p: BoxedUint = generate_prime(bit_length);
            assert!(p.bits_vartime() == bit_length);
            assert!(p.to_words().len() == nlimbs!(bit_length));
            assert!(is_prime(&p));
            assert!(fips_is_prime(&p));
        }
    }

    #[test]
    fn safe_prime_generation() {
        for bit_length in (28..=128).step_by(10) {
            let p: U128 = generate_safe_prime(bit_length);
            assert!(p.bits_vartime() == bit_length);
            assert!(is_safe_prime(&p));
            assert!(fips_is_safe_prime(&p));
        }
    }

    #[test]
    fn safe_prime_generation_boxed() {
        for bit_length in (28..=189).step_by(10) {
            let p: BoxedUint = generate_safe_prime(bit_length);
            assert!(p.bits_vartime() == bit_length);
            assert!(p.to_words().len() == nlimbs!(bit_length));
            assert!(is_safe_prime(&p));
            assert!(fips_is_safe_prime(&p));
        }
    }

    #[test]
    fn corner_cases_is_prime() {
        for num in 0u64..30 {
            let is_prime_ref = is_prime64(num);
            let is_safe_prime_ref = is_prime_ref && is_prime64(num / 2);

            let num_uint = U64::from(num);

            let is_prime_test = is_prime(&num_uint);
            assert_eq!(
                is_prime_ref, is_prime_test,
                "num={num}, expected={is_prime_ref}, actual={is_prime_test}"
            );

            let is_safe_prime_test = is_safe_prime(&num_uint);
            assert_eq!(
                is_safe_prime_ref, is_safe_prime_test,
                "num={num}, expected={is_safe_prime_ref}, actual={is_safe_prime_test}"
            );

            let is_prime_test = fips_is_prime(&num_uint);
            assert_eq!(
                is_prime_ref, is_prime_test,
                "num={num}, expected={is_prime_ref}, actual={is_prime_test}"
            );

            let is_safe_prime_test = fips_is_safe_prime(&num_uint);
            assert_eq!(
                is_safe_prime_ref, is_safe_prime_test,
                "num={num}, expected={is_safe_prime_ref}, actual={is_safe_prime_test}"
            );
        }
    }

    #[test]
    fn inconclusive_sieving_result() {
        // Coverage test.
        // This number is a product of two primes larger than the maximum prime in `SMALL_PRIMES`,
        // so the initial sieving cannot tell if it is prime or not,
        // and a full primality test is run.
        assert!(!is_safe_prime(&U64::from(17881u32 * 17891u32)));
        assert!(!fips_is_safe_prime(&U64::from(17881u32 * 17891u32)));
    }

    #[test]
    #[should_panic(expected = "random_odd_integer() failed: BitLengthTooLarge { bit_length: 65, bits_precision: 64 }")]
    fn generate_prime_too_many_bits() {
        let _p: U64 = generate_prime_with_rng(&mut OsRng.unwrap_err(), 65);
    }

    #[test]
    #[should_panic(expected = "random_odd_integer() failed: BitLengthTooLarge { bit_length: 65, bits_precision: 64 }")]
    fn generate_safe_prime_too_many_bits() {
        let _p: U64 = generate_safe_prime_with_rng(&mut OsRng.unwrap_err(), 65);
    }

    fn is_prime_ref(num: Word) -> bool {
        num_prime::nt_funcs::is_prime(&num, None).probably()
    }

    #[test]
    fn corner_cases_generate_prime() {
        for bits in 2..5 {
            for _ in 0..100 {
                let p: U64 = generate_prime(bits);
                let p_word = p.as_words()[0];
                assert!(is_prime_ref(p_word));
            }
        }
    }

    #[test]
    fn corner_cases_generate_safe_prime() {
        for bits in 3..5 {
            for _ in 0..100 {
                let p: U64 = generate_safe_prime(bits);
                let p_word = p.as_words()[0];
                assert!(is_prime_ref(p_word) && is_prime_ref(p_word / 2));
            }
        }
    }
}

#[cfg(all(test, feature = "multicore"))]
mod multicore_tests {
    use super::{is_prime, par_generate_prime, par_generate_safe_prime};
    use crypto_bigint::{nlimbs, BoxedUint, U128};

    #[test]
    fn parallel_prime_generation() {
        for bit_length in (28..=128).step_by(10) {
            let p: U128 = par_generate_prime(bit_length, 4);
            assert!(p.bits_vartime() == bit_length);
            assert!(is_prime(&p));
        }
    }

    #[test]
    fn parallel_prime_generation_boxed() {
        for bit_length in (28..=128).step_by(10) {
            let p: BoxedUint = par_generate_prime(bit_length, 2);
            assert!(p.bits_vartime() == bit_length);
            assert!(p.to_words().len() == nlimbs!(bit_length));
            assert!(is_prime(&p));
        }
    }

    #[test]
    fn parallel_safe_prime_generation() {
        for bit_length in (28..=128).step_by(10) {
            let p: U128 = par_generate_safe_prime(bit_length, 8);
            assert!(p.bits_vartime() == bit_length);
            assert!(is_prime(&p));
        }
    }

    #[test]
    fn parallel_safe_prime_generation_boxed() {
        for bit_length in (28..=128).step_by(10) {
            let p: BoxedUint = par_generate_safe_prime(bit_length, 4);
            assert!(p.bits_vartime() == bit_length);
            assert!(p.to_words().len() == nlimbs!(bit_length));
            assert!(is_prime(&p));
        }
    }
}

#[cfg(test)]
#[cfg(feature = "tests-openssl")]
mod tests_openssl {
    use alloc::format;
    use core::num::NonZero;

    use crypto_bigint::U128;
    use openssl::bn::{BigNum, BigNumContext};
    use rand_core::{OsRng, TryRngCore};

    use super::{fips_is_prime_with_rng, generate_prime, is_prime};
    use crate::hazmat::{minimum_mr_iterations, random_odd_integer, SetBits};

    fn openssl_is_prime(num: &BigNum, ctx: &mut BigNumContext) -> bool {
        num.is_prime(64, ctx).unwrap()
    }

    fn to_openssl(num: &U128) -> BigNum {
        BigNum::from_hex_str(&format!("{num:x}")).unwrap()
    }

    fn from_openssl(num: &BigNum) -> U128 {
        U128::from_be_hex(&num.to_hex_str().unwrap())
    }

    #[test]
    fn openssl_cross_check() {
        let mut ctx = BigNumContext::new().unwrap();

        // Generate primes, let OpenSSL check them
        for _ in 0..100 {
            let p: U128 = generate_prime(128);
            let p_bn = to_openssl(&p);
            assert!(openssl_is_prime(&p_bn, &mut ctx), "OpenSSL reports {p} as composite",);
        }

        let mr_iterations = minimum_mr_iterations(U128::BITS, 100).unwrap();

        // Generate primes with OpenSSL, check them
        let mut p_bn = BigNum::new().unwrap();
        for _ in 0..100 {
            p_bn.generate_prime(128, false, None, None).unwrap();
            let p = from_openssl(&p_bn);
            assert!(is_prime(&p), "we report {p} as composite");
            assert!(
                fips_is_prime_with_rng(&mut OsRng.unwrap_err(), &p, mr_iterations, false),
                "we report {p} as composite"
            );
        }

        // Generate random numbers, check if our test agrees with OpenSSL
        for _ in 0..100 {
            let p = random_odd_integer::<U128, _>(&mut OsRng, NonZero::new(128).unwrap(), SetBits::Msb).unwrap();
            let p_bn = to_openssl(&p);
            let expected = openssl_is_prime(&p_bn, &mut ctx);

            let actual = is_prime(p.as_ref());
            assert_eq!(
                actual, expected,
                "difference between OpenSSL and us: OpenSSL reports {expected}, we report {actual}",
            );

            let actual = fips_is_prime_with_rng(&mut OsRng.unwrap_err(), p.as_ref(), mr_iterations, false);
            assert_eq!(
                actual, expected,
                "difference between OpenSSL and us: OpenSSL reports {expected}, we report {actual}",
            );
        }
    }
}

#[cfg(test)]
#[cfg(feature = "tests-gmp")]
mod tests_gmp {
    use core::num::NonZero;

    use crypto_bigint::U128;
    use rand_core::{OsRng, TryRngCore};
    use rug::{
        integer::{IsPrime, Order},
        Integer,
    };

    use super::{fips_is_prime_with_rng, generate_prime, is_prime};
    use crate::hazmat::{minimum_mr_iterations, random_odd_integer, SetBits};

    fn gmp_is_prime(num: &Integer) -> bool {
        matches!(num.is_probably_prime(25), IsPrime::Yes | IsPrime::Probably)
    }

    fn to_gmp(num: &U128) -> Integer {
        Integer::from_digits(num.as_words(), Order::Lsf)
    }

    fn from_gmp(num: &Integer) -> U128 {
        U128::from_words(num.to_digits(Order::Lsf).try_into().unwrap())
    }

    #[test]
    fn gmp_cross_check() {
        // Generate primes, let GMP check them
        for _ in 0..100 {
            let p: U128 = generate_prime(128);
            let p_bn = to_gmp(&p);
            assert!(gmp_is_prime(&p_bn), "GMP reports {p} as composite");
        }

        let mr_iterations = minimum_mr_iterations(U128::BITS, 100).unwrap();

        // Generate primes with GMP, check them
        for _ in 0..100 {
            let start =
                random_odd_integer::<U128, _>(&mut OsRng.unwrap_err(), NonZero::new(128).unwrap(), SetBits::Msb)
                    .unwrap();
            let start_bn = to_gmp(&start);
            let p_bn = start_bn.next_prime();
            let p = from_gmp(&p_bn);
            assert!(is_prime(&p), "we report {p} as composite");
            assert!(
                fips_is_prime_with_rng(&mut OsRng.unwrap_err(), &p, mr_iterations, false),
                "we report {p} as composite"
            );
        }

        // Generate random numbers, check if our test agrees with GMP
        for _ in 0..100 {
            let p = random_odd_integer::<U128, _>(&mut OsRng.unwrap_err(), NonZero::new(128).unwrap(), SetBits::Msb)
                .unwrap();
            let p_bn = to_gmp(&p);
            let expected = gmp_is_prime(&p_bn);

            let actual = is_prime(p.as_ref());
            assert_eq!(
                actual, expected,
                "difference between GMP and us: GMP reports {expected}, we report {actual}",
            );

            let actual = fips_is_prime_with_rng(&mut OsRng.unwrap_err(), p.as_ref(), mr_iterations, false);
            assert_eq!(
                actual, expected,
                "difference between GMP and us: GMP reports {expected}, we report {actual}",
            );
        }
    }
}
