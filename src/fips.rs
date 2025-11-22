//! Functions implementing the functionality prescribed by FIPS-186.5 standard[^FIPS].
//!
//! [^FIPS]: FIPS-186.5 standard, <https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf>

use core::num::NonZero;

use crypto_bigint::{Odd, RandomMod, Unsigned};
use rand_core::CryptoRng;

use crate::{
    hazmat::{LucasCheck, MillerRabin, Primality, SelfridgeBase, SmallFactorsSieve, equals_primitive, lucas_test},
    presets::Flavor,
};

/// Options for FIPS primality testing.
#[derive(Copy, Clone, Debug, Default)]
pub struct FipsOptions {
    /// Use an additional regular Lucas check with Selfridge base (see [`SelfridgeBase`] for details).
    ///
    /// `false` by default.
    pub add_lucas_test: bool,
    /// Use a trial division as explained in Appendix B.3.
    ///
    /// `false` by default.
    ///
    /// Note that it is a performance optimization for the cases when you expect the candidates to be random
    /// (and thus likely to have small factors).
    /// It does not affect the failure probability of the primality check.
    pub add_trial_division_test: bool,
}

/// Probabilistically checks if the given number is prime using the provided RNG
/// according to FIPS-186.5[^FIPS] standard.
///
/// By default, performs `mr_iterations` of Miller-Rabin check with random bases.
/// See [`MillerRabin`] and [`lucas_test`] for more details about the checks;
/// use [`minimum_mr_iterations`](`crate::hazmat::minimum_mr_iterations`)
/// to calculate the number of required iterations.
///
/// Additional checks can be specified in the [`FipsOptions`] structure.
///
/// [^FIPS]: FIPS-186.5 standard, <https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf>
pub fn is_prime<T>(
    rng: &mut (impl CryptoRng + ?Sized),
    flavor: Flavor,
    candidate: &T,
    mr_iterations: usize,
    options: FipsOptions,
) -> bool
where
    T: Unsigned + RandomMod,
{
    match flavor {
        Flavor::Any => {}
        Flavor::Safe => return is_safe_prime(rng, candidate, mr_iterations, options),
    }

    if equals_primitive(candidate, 0) || equals_primitive(candidate, 1) {
        return false;
    }

    if equals_primitive(candidate, 2) {
        return true;
    }

    let canditate_bits = match NonZero::new(candidate.bits()) {
        Some(x) => x,
        None => return false,
    };

    if options.add_trial_division_test {
        let mut sieve =
            SmallFactorsSieve::new(candidate.clone(), canditate_bits, flavor.eq(&Flavor::Safe)).expect(concat![
                "canditate_bits is always less than or equal to candidate.bits_precision(), ",
                "which is the only way to get an error from this function"
            ]);
        sieve.update_residues();
        if sieve.current_is_composite() {
            return false;
        }
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

    if options.add_lucas_test {
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
/// See [`fips_is_prime`] for details about the performed checks.
///
/// [^FIPS]: FIPS-186.5 standard, <https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf>
fn is_safe_prime<T>(
    rng: &mut (impl CryptoRng + ?Sized),
    candidate: &T,
    mr_iterations: usize,
    options: FipsOptions,
) -> bool
where
    T: Unsigned + RandomMod,
{
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

    is_prime(rng, Flavor::Any, candidate, mr_iterations, options)
        && is_prime(
            rng,
            Flavor::Any,
            &candidate.wrapping_shr_vartime(1),
            mr_iterations,
            options,
        )
}
