use core::num::NonZeroU32;

use crypto_bigint::{Integer, Limb, Odd, RandomBits, RandomMod};
use rand_core::CryptoRngCore;

#[cfg(feature = "default-rng")]
use rand_core::OsRng;

use crate::hazmat::{
    lucas_test, random_odd_integer, AStarBase, LucasCheck, MillerRabin, Primality, Sieve,
};

/// Returns a random prime of size `bit_length` using [`OsRng`] as the RNG.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
#[cfg(feature = "default-rng")]
pub fn generate_prime<T: Integer + RandomBits + RandomMod>(bit_length: u32) -> T {
    generate_prime_with_rng(&mut OsRng, bit_length)
}

/// Returns a random safe prime (that is, such that `(n - 1) / 2` is also prime) of size
/// `bit_length` using [`OsRng`] as the RNG.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
#[cfg(feature = "default-rng")]
pub fn generate_safe_prime<T: Integer + RandomBits + RandomMod>(bit_length: u32) -> T {
    generate_safe_prime_with_rng(&mut OsRng, bit_length)
}

/// Probabilistically checks if the given number is prime using [`OsRng`] as the RNG.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
#[cfg(feature = "default-rng")]
pub fn is_prime<T: Integer + RandomMod>(num: &T) -> bool {
    is_prime_with_rng(&mut OsRng, num)
}

/// Probabilistically checks if the given number is a safe prime (that is, such that `(n - 1) / 2` is
/// also prime) using [`OsRng`] as the RNG.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
#[cfg(feature = "default-rng")]
pub fn is_safe_prime<T: Integer + RandomMod>(num: &T) -> bool {
    is_safe_prime_with_rng(&mut OsRng, num)
}

/// Returns a random prime of size `bit_length` using the provided RNG.
///
/// Panics if `bit_length` is less than 2, or greater than the bit size of the target `Uint`.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
pub fn generate_prime_with_rng<T: Integer + RandomBits + RandomMod>(
    rng: &mut impl CryptoRngCore,
    bit_length: u32,
) -> T {
    if bit_length < 2 {
        panic!("`bit_length` must be 2 or greater.");
    }
    let bit_length = NonZeroU32::new(bit_length).expect("`bit_length` should be non-zero");
    loop {
        let start = random_odd_integer::<T>(rng, bit_length);
        let sieve = Sieve::new(start.get(), bit_length, false);
        for num in sieve {
            if is_prime_with_rng(rng, &num) {
                return num;
            }
        }
    }
}

/// Returns a random safe prime (that is, such that `(n - 1) / 2` is also prime)
/// of size `bit_length` using the provided RNG.
///
/// Panics if `bit_length` is less than 3, or greater than the bit size of the target `Uint`.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
pub fn generate_safe_prime_with_rng<T: Integer + RandomBits + RandomMod>(
    rng: &mut impl CryptoRngCore,
    bit_length: u32,
) -> T {
    if bit_length < 3 {
        panic!("`bit_length` must be 3 or greater.");
    }
    let bit_length = NonZeroU32::new(bit_length).expect("`bit_length` should be non-zero");
    loop {
        let start = random_odd_integer::<T>(rng, bit_length);
        let sieve = Sieve::new(start.get(), bit_length, true);
        for num in sieve {
            if is_safe_prime_with_rng(rng, &num) {
                return num;
            }
        }
    }
}

/// Probabilistically checks if the given number is prime using the provided RNG.
///
/// Performed checks:
/// - Miller-Rabin check with base 2;
/// - Strong Lucas check with A* base (see [`AStarBase`] for details);
/// - Miller-Rabin check with a random base.
///
/// See [`MillerRabin`] and [`lucas_test`] for more details about the checks.
///
/// The first two checks constitute the Baillie-PSW primality test[^Baillie1980];
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
pub fn is_prime_with_rng<T: Integer + RandomMod>(rng: &mut impl CryptoRngCore, num: &T) -> bool {
    if num == &T::from_limb_like(Limb::from(2u32), num) {
        return true;
    }

    match Odd::new(num.clone()).into() {
        Some(x) => _is_prime_with_rng(rng, x),
        None => false,
    }
}

/// Probabilistically checks if the given number is a safe prime using the provided RNG.
///
/// See [`is_prime_with_rng`] for details about the performed checks.
pub fn is_safe_prime_with_rng<T: Integer + RandomMod>(
    rng: &mut impl CryptoRngCore,
    num: &T,
) -> bool {
    // Since, by the definition of safe prime, `(num - 1) / 2` must also be prime,
    // and therefore odd, `num` has to be equal to 3 modulo 4.
    // 5 is the only exception, so we check for it.
    if num == &T::from_limb_like(Limb::from(5u32), num) {
        return true;
    }

    // Safe primes are always of the form 4k + 3 (i.e. n ≡ 3 mod 4)
    // The last two digits of a binary number give you its value modulo 4.
    // Primes p=4n+3 will always end in 11​ in binary because p ≡ 3 mod 4.
    if num.as_ref()[0].0 & 3 != 3 {
        return false;
    }

    // These are ensured to be odd by the check above.
    let odd_num = Odd::new(num.clone()).expect("`num` is be odd here given the checks above");
    let odd_half_num = Odd::new(num.wrapping_shr_vartime(1)).expect("The binary rep of `num` is `11` so shifting right by one is guaranteed to end in one, so it's odd");

    _is_prime_with_rng(rng, odd_num) && _is_prime_with_rng(rng, odd_half_num)
}

/// Checks for primality.
fn _is_prime_with_rng<T: Integer + RandomMod>(rng: &mut impl CryptoRngCore, num: Odd<T>) -> bool {
    let candidate = num.clone();
    let mr = MillerRabin::new(num);

    if !mr.test_base_two().is_probably_prime() {
        return false;
    }

    match lucas_test(candidate, AStarBase, LucasCheck::Strong) {
        Primality::Composite => return false,
        Primality::Prime => return true,
        _ => {}
    }

    // The random base test only makes sense when `num > 3`.
    if mr.bit_length() > 2 && !mr.test_random_base(rng).is_probably_prime() {
        return false;
    }

    true
}

#[cfg(test)]
mod tests {
    use crypto_bigint::{nlimbs, BoxedUint, CheckedAdd, Uint, Word, U128, U64};
    use num_prime::nt_funcs::is_prime64;
    use rand_core::OsRng;

    use super::{
        generate_prime, generate_prime_with_rng, generate_safe_prime, generate_safe_prime_with_rng,
        is_prime, is_safe_prime,
    };
    use crate::hazmat::{primes, pseudoprimes};

    fn test_large_primes<const L: usize>(nums: &[Uint<L>]) {
        for num in nums {
            assert!(is_prime(num));
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
        }

        assert!(!is_prime(&pseudoprimes::LARGE_CARMICHAEL_NUMBER));
    }

    fn test_cunningham_chain<const L: usize>(length: usize, num: &Uint<L>) {
        let mut next = *num;
        for i in 0..length {
            assert!(is_prime(&next));

            // The start of the chain isn't a safe prime by definition
            if i > 0 {
                assert!(is_safe_prime(&next));
            }

            next = next
                .wrapping_shl_vartime(1)
                .checked_add(&Uint::<L>::ONE)
                .unwrap();
        }

        // The chain ended.
        assert!(!is_prime(&next));
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
        }
    }

    #[test]
    fn prime_generation_boxed() {
        for bit_length in (28..=128).step_by(10) {
            let p: BoxedUint = generate_prime(bit_length);
            assert!(p.bits_vartime() == bit_length);
            assert!(p.to_words().len() == nlimbs!(bit_length));
            assert!(is_prime(&p));
        }
    }

    #[test]
    fn safe_prime_generation() {
        for bit_length in (28..=128).step_by(10) {
            let p: U128 = generate_safe_prime(bit_length);
            assert!(p.bits_vartime() == bit_length);
            assert!(is_safe_prime(&p));
        }
    }

    #[test]
    fn safe_prime_generation_boxed() {
        for bit_length in (28..=189).step_by(10) {
            let p: BoxedUint = generate_safe_prime(bit_length);
            assert!(p.bits_vartime() == bit_length);
            assert!(p.to_words().len() == nlimbs!(bit_length));
            assert!(is_safe_prime(&p));
        }
    }

    #[test]
    fn corner_cases_is_prime() {
        for num in 0u64..30 {
            let is_prime_ref = is_prime64(num);
            let is_safe_prime_ref = is_prime_ref && is_prime64(num / 2);

            let num_uint = U64::from(num);
            let is_prime_test = is_prime(&num_uint);
            let is_safe_prime_test = is_safe_prime(&num_uint);

            assert_eq!(is_prime_ref, is_prime_test);
            assert_eq!(is_safe_prime_ref, is_safe_prime_test);
        }
    }

    #[test]
    fn inconclusive_sieving_result() {
        // Coverage test.
        // This number is a product of two primes larger than the maximum prime in `SMALL_PRIMES`,
        // so the initial sieving cannot tell if it is prime or not,
        // and a full primality test is run.
        assert!(!is_safe_prime(&U64::from(17881u32 * 17891u32)));
    }

    #[test]
    #[should_panic(expected = "`bit_length` must be 2 or greater")]
    fn generate_prime_too_few_bits() {
        let _p: U64 = generate_prime_with_rng(&mut OsRng, 1);
    }

    #[test]
    #[should_panic(expected = "`bit_length` must be 3 or greater")]
    fn generate_safe_prime_too_few_bits() {
        let _p: U64 = generate_safe_prime_with_rng(&mut OsRng, 2);
    }

    #[test]
    #[should_panic(
        expected = "try_random_bits() failed: BitLengthTooLarge { bit_length: 65, bits_precision: 64 }"
    )]
    fn generate_prime_too_many_bits() {
        let _p: U64 = generate_prime_with_rng(&mut OsRng, 65);
    }

    #[test]
    #[should_panic(
        expected = "try_random_bits() failed: BitLengthTooLarge { bit_length: 65, bits_precision: 64 }"
    )]
    fn generate_safe_prime_too_many_bits() {
        let _p: U64 = generate_safe_prime_with_rng(&mut OsRng, 65);
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

#[cfg(test)]
#[cfg(feature = "tests-openssl")]
mod tests_openssl {
    use alloc::format;
    use core::num::NonZeroU32;

    use crypto_bigint::U128;
    use openssl::bn::{BigNum, BigNumContext};
    use rand_core::OsRng;

    use super::{generate_prime, is_prime};
    use crate::hazmat::random_odd_integer;

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
            assert!(
                openssl_is_prime(&p_bn, &mut ctx),
                "OpenSSL reports {p} as composite",
            );
        }

        // Generate primes with OpenSSL, check them
        let mut p_bn = BigNum::new().unwrap();
        for _ in 0..100 {
            p_bn.generate_prime(128, false, None, None).unwrap();
            let p = from_openssl(&p_bn);
            assert!(is_prime(&p), "we report {p} as composite");
        }

        // Generate random numbers, check if our test agrees with OpenSSL
        for _ in 0..100 {
            let p = random_odd_integer::<U128>(&mut OsRng, NonZeroU32::new(128).unwrap());
            let actual = is_prime(p.as_ref());
            let p_bn = to_openssl(&p);
            let expected = openssl_is_prime(&p_bn, &mut ctx);
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
    use core::num::NonZeroU32;

    use crypto_bigint::U128;
    use rand_core::OsRng;
    use rug::{
        integer::{IsPrime, Order},
        Integer,
    };

    use super::{generate_prime, is_prime};
    use crate::hazmat::random_odd_integer;

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

        // Generate primes with GMP, check them
        for _ in 0..100 {
            let start = random_odd_integer::<U128>(&mut OsRng, NonZeroU32::new(128).unwrap());
            let start_bn = to_gmp(&start);
            let p_bn = start_bn.next_prime();
            let p = from_gmp(&p_bn);
            assert!(is_prime(&p), "we report {p} as composite");
        }

        // Generate random numbers, check if our test agrees with GMP
        for _ in 0..100 {
            let p = random_odd_integer::<U128>(&mut OsRng, NonZeroU32::new(128).unwrap());
            let actual = is_prime(p.as_ref());
            let p_bn = to_gmp(&p);
            let expected = gmp_is_prime(&p_bn);
            assert_eq!(
                actual, expected,
                "difference between GMP and us: GMP reports {expected}, we report {actual}",
            );
        }
    }
}
