//! An iterator for weeding out multiples of small primes,
//! before proceeding with slower tests.

use alloc::{vec, vec::Vec};

use crypto_bigint::{Integer, Limb, NonZero, Random, Uint, Zero};
use rand_core::CryptoRngCore;

use crate::hazmat::{
    precomputed::{SmallPrime, RECIPROCALS, SMALL_PRIMES},
    Primality,
};

/// Returns a random odd integer with given bit length
/// (that is, with both `0` and `bit_length-1` bits set).
///
/// Panics if `bit_length` is 0 or is greater than the bit size of the target `Uint`.
pub fn random_odd_uint<const L: usize>(rng: &mut impl CryptoRngCore, bit_length: usize) -> Uint<L> {
    if bit_length == 0 {
        panic!("Bit length must be non-zero");
    }

    if bit_length > Uint::<L>::BITS {
        panic!(
            "The requested bit length ({}) is larger than the chosen Uint size",
            bit_length
        );
    }

    // TODO: not particularly efficient, can be improved by zeroing high bits instead of shifting
    let mut random = Uint::<L>::random(rng);
    if bit_length != Uint::<L>::BITS {
        random >>= Uint::<L>::BITS - bit_length;
    }

    // Make it odd
    random |= Uint::<L>::ONE;

    // Make sure it's the correct bit size
    random |= Uint::<L>::ONE << (bit_length - 1);

    random
}

/// An iterator returning numbers with up to and including given bit length,
/// starting from a given number, that are not multiples of the first 2048 small primes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sieve<const L: usize> {
    // Instead of dividing `Uint` by small primes every time (which is slow),
    // we keep a "base" and a small increment separately,
    // so that we can only calculate the residues of the increment.
    base: Uint<L>,
    incr: Residue,
    incr_limit: Residue,
    safe_primes: bool,
    residues: Vec<SmallPrime>,
    max_bit_length: usize,
    produces_nothing: bool,
    starts_from_exception: bool,
    last_round: bool,
}

// The type we use to calculate incremental residues.
// Should be >= `SmallPrime` in size.
type Residue = u32;

// The maximum increment that won't overflow the type we use to calculate residues of increments:
// we need `(max_prime - 1) + max_incr <= Type::MAX`.
const INCR_LIMIT: Residue = Residue::MAX - SMALL_PRIMES[SMALL_PRIMES.len() - 1] as Residue + 1;

impl<const L: usize> Sieve<L> {
    /// Creates a new sieve, iterating from `start` and
    /// until the last number with `max_bit_length` bits,
    /// producing numbers that are not non-trivial multiples
    /// of a list of small primes in the range `[2, start)`.
    ///
    /// Note that `start` is adjusted to `2`, or the next `1 mod 2` number (`safe_prime = false`);
    /// and `5`, or `3 mod 4` number (`safe_prime = true`).
    ///
    /// If `safe_primes` is `true`, it only produces the candidates
    /// that can possibly be safe primes (that is, 5, and those equal to 3 modulo 4).
    pub fn new(start: &Uint<L>, max_bit_length: usize, safe_primes: bool) -> Self {
        let mut base = *start;

        // This is easier than making all the methods generic enough to handle these corner cases.
        let produces_nothing = max_bit_length < base.bits()
            || (!safe_primes && max_bit_length < 2)
            || (safe_primes && max_bit_length < 3);

        // Add exceptions to the produced candidates - the only ones that don't fit
        // the general pattern of incrementing the base by 2 or by 4.
        let mut starts_from_exception = false;
        if !safe_primes {
            if base <= Uint::<L>::from(2u32) {
                starts_from_exception = true;
                base = Uint::<L>::from(3u32);
            } else {
                // Adjust the base so that we hit correct numbers when incrementing it by 2.
                base |= Uint::<L>::ONE;
            }
        } else if base <= Uint::<L>::from(5u32) {
            starts_from_exception = true;
            base = Uint::<L>::from(7u32);
        } else {
            // Adjust the base so that we hit correct numbers when incrementing it by 4.
            // If we are looking for safe primes, we are only interested
            // in the numbers == 3 mod 4.
            base |= Uint::<L>::from(3u32);
        }

        // Only calculate residues by primes up to and not including `base`,
        // because when we only have the resiude,
        // we cannot distinguish between a prime itself and a multiple of that prime.
        let residues_len = if Uint::<L>::from(SMALL_PRIMES[SMALL_PRIMES.len() - 1]) >= base {
            SMALL_PRIMES
                .iter()
                .enumerate()
                .find(|(_i, p)| Uint::<L>::from(**p) >= base)
                .map(|(i, _p)| i)
                .unwrap_or(SMALL_PRIMES.len())
        } else {
            // This will be the majority of use cases
            SMALL_PRIMES.len()
        };

        Self {
            base,
            incr: 0, // This will ensure that `update_residues()` is called right away.
            incr_limit: 0,
            safe_primes,
            residues: vec![0; residues_len],
            max_bit_length,
            produces_nothing,
            starts_from_exception,
            last_round: false,
        }
    }

    fn update_residues(&mut self) -> bool {
        if self.incr_limit != 0 && self.incr <= self.incr_limit {
            return true;
        }

        if self.last_round {
            return false;
        }

        // Set the new base.
        self.base = self.base.wrapping_add(&self.incr.into());

        self.incr = 0;

        // Re-calculate residues.
        for (i, rec) in RECIPROCALS.iter().enumerate().take(self.residues.len()) {
            let (_quo, rem) = self.base.ct_div_rem_limb_with_reciprocal(rec);
            self.residues[i] = rem.0 as SmallPrime;
        }

        // Find the increment limit.
        // Note that the max value is the same regardless of the value of `self.safe_prime`,
        // since `(2^n - 1) = 3 mod 4` (for n > 1).
        let max_value = (Uint::<L>::ONE << self.max_bit_length).wrapping_sub(&Uint::<L>::ONE);
        let incr_limit = max_value.wrapping_sub(&self.base);
        self.incr_limit = if incr_limit > INCR_LIMIT.into() {
            INCR_LIMIT
        } else {
            // We are close to `2^max_bit_length - 1`.
            // Mark this round as the last.
            self.last_round = true;
            // Can unwrap here since we just checked above that `incr_limit <= INCR_LIMIT`,
            // and `INCR_LIMIT` fits into `Residue`.
            let incr_limit_small: Residue = incr_limit.as_words()[0].try_into().unwrap();
            incr_limit_small
        };

        true
    }

    // Returns `true` if the current `base + incr` is divisible by any of the small primes.
    fn current_is_composite(&self) -> bool {
        for (i, m) in self.residues.iter().enumerate() {
            let r = (*m as Residue + self.incr) % (SMALL_PRIMES[i] as Residue);
            if r == 0 {
                return true;
            }
        }
        false
    }

    // Returns the restored `base + incr` if it is not composite (wrt the small primes),
    // and bumps the increment unconditionally.
    fn maybe_next(&mut self) -> Option<Uint<L>> {
        let result = if self.current_is_composite() {
            None
        } else {
            Some(self.base.wrapping_add(&self.incr.into()))
        };

        // If we are looking for safe primes, we are only interested in the numbers == 3 mod 4.
        self.incr += if self.safe_primes { 4 } else { 2 };
        result
    }

    fn next(&mut self) -> Option<Uint<L>> {
        // Corner cases handled here

        if self.produces_nothing {
            return None;
        }

        if self.starts_from_exception {
            self.starts_from_exception = false;
            return Some(Uint::<L>::from(if self.safe_primes { 5u32 } else { 2u32 }));
        }

        // Main loop

        while self.update_residues() {
            match self.maybe_next() {
                Some(x) => return Some(x),
                None => continue,
            };
        }
        None
    }
}

impl<const L: usize> Iterator for Sieve<L> {
    type Item = Uint<L>;

    fn next(&mut self) -> Option<Self::Item> {
        Sieve::<L>::next(self)
    }
}

/// Performs trial division of the given number by a list of small primes starting from 2.
/// Returns `Some(primality)` if there was a definitive conclusion about `num`'s primality,
/// and `None` otherwise.
pub fn sieve_once<const L: usize>(num: &Uint<L>) -> Option<Primality> {
    // Our reciprocals start from 3, so we check for 2 separately
    if num == &Uint::<L>::from(2u32) {
        return Some(Primality::Prime);
    }
    if num.is_even().into() {
        return Some(Primality::Composite);
    }
    for small_prime in SMALL_PRIMES.iter() {
        let (quo, rem) = num.div_rem_limb(NonZero::new(Limb::from(*small_prime)).unwrap());
        if rem.is_zero().into() {
            let primality = if quo == Uint::<L>::ONE {
                Primality::Prime
            } else {
                Primality::Composite
            };
            return Some(primality);
        }
    }
    None
}

#[cfg(test)]
mod tests {

    use alloc::format;
    use alloc::vec::Vec;

    use crypto_bigint::U64;
    use num_prime::nt_funcs::factorize64;
    use rand_chacha::ChaCha8Rng;
    use rand_core::{OsRng, SeedableRng};

    use super::{random_odd_uint, Sieve};
    use crate::hazmat::precomputed::SMALL_PRIMES;

    #[test]
    fn random() {
        let max_prime = SMALL_PRIMES[SMALL_PRIMES.len() - 1];

        let mut rng = ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");
        let start: U64 = random_odd_uint(&mut rng, 32);
        for num in Sieve::new(&start, 32, false).take(100) {
            let num_u64: u64 = num.into();
            assert!(num_u64.leading_zeros() == 32);

            let factors_and_powers = factorize64(num_u64);
            let factors = factors_and_powers.into_keys().collect::<Vec<_>>();

            assert!(factors[0] > max_prime as u64);
        }
    }

    fn check_sieve(start: u32, bit_length: usize, safe_prime: bool, reference: &[u32]) {
        let test = Sieve::new(&U64::from(start), bit_length, safe_prime).collect::<Vec<_>>();
        assert_eq!(test.len(), reference.len());
        for (x, y) in test.iter().zip(reference.iter()) {
            assert_eq!(x, &U64::from(*y));
        }
    }

    #[test]
    fn empty_sequence() {
        check_sieve(1, 1, false, &[]); // no primes of 1 bits
        check_sieve(1, 2, true, &[]); // no safe primes of 2 bits
        check_sieve(64, 6, true, &[]); // 64 is 7 bits long
    }

    #[test]
    fn small_range() {
        check_sieve(1, 2, false, &[2, 3]);
        check_sieve(2, 2, false, &[2, 3]);
        check_sieve(3, 2, false, &[3]);

        check_sieve(1, 3, false, &[2, 3, 5, 7]);
        check_sieve(3, 3, false, &[3, 5, 7]);
        check_sieve(5, 3, false, &[5, 7]);
        check_sieve(7, 3, false, &[7]);

        check_sieve(1, 4, false, &[2, 3, 5, 7, 9, 11, 13, 15]);
        check_sieve(3, 4, false, &[3, 5, 7, 9, 11, 13, 15]);
        check_sieve(5, 4, false, &[5, 7, 11, 13]);
        check_sieve(7, 4, false, &[7, 11, 13]);
        check_sieve(9, 4, false, &[11, 13]);
        check_sieve(13, 4, false, &[13]);
        check_sieve(15, 4, false, &[]);

        check_sieve(1, 3, true, &[5, 7]);
        check_sieve(3, 3, true, &[5, 7]);
        check_sieve(5, 3, true, &[5, 7]);
        check_sieve(7, 3, true, &[7]);

        // start is adjusted to 5 here, so multiples of 3 can be excluded
        check_sieve(1, 4, true, &[5, 7, 11]);

        check_sieve(5, 4, true, &[5, 7, 11]);
        check_sieve(7, 4, true, &[7, 11]);
        check_sieve(9, 4, true, &[11]);
        check_sieve(13, 4, true, &[]);
    }

    #[test]
    fn random_below_max_length() {
        for _ in 0..10 {
            let r: U64 = random_odd_uint(&mut OsRng, 50);
            assert_eq!(r.bits(), 50);
        }
    }

    #[test]
    #[should_panic(expected = "Bit length must be non-zero")]
    fn random_odd_uint_0bits() {
        let _p: U64 = random_odd_uint(&mut OsRng, 0);
    }

    #[test]
    #[should_panic(expected = "The requested bit length (65) is larger than the chosen Uint size")]
    fn random_odd_uint_too_many_bits() {
        let _p: U64 = random_odd_uint(&mut OsRng, 65);
    }

    #[test]
    fn sieve_derived_traits() {
        let s = Sieve::new(&U64::ONE, 10, false);
        assert!(format!("{s:?}").starts_with("Sieve"));
        assert_eq!(s.clone(), s);
    }
}
