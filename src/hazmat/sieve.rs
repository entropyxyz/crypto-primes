//! An iterator for weeding out multiples of small primes,
//! before proceeding with slower tests.

use alloc::{vec, vec::Vec};

use crypto_bigint::{CheckedAdd, Odd, Random, Uint};
use rand_core::CryptoRngCore;

use crate::hazmat::precomputed::{SmallPrime, RECIPROCALS, SMALL_PRIMES};

/// Returns a random odd integer with given bit length
/// (that is, with both `0` and `bit_length-1` bits set).
///
/// Panics if `bit_length` is 0 or is greater than the bit size of the target `Uint`.
pub fn random_odd_uint<const L: usize>(
    rng: &mut impl CryptoRngCore,
    bit_length: u32,
) -> Odd<Uint<L>> {
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
    random |= Uint::<L>::ONE.wrapping_shl_vartime(bit_length - 1);

    Odd::new(random).expect("ensured to be odd")
}

// The type we use to calculate incremental residues.
// Should be >= `SmallPrime` in size.
type Residue = u32;

// The maximum increment that won't overflow the type we use to calculate residues of increments:
// we need `(max_prime - 1) + max_incr <= Type::MAX`.
const INCR_LIMIT: Residue = Residue::MAX - SMALL_PRIMES[SMALL_PRIMES.len() - 1] as Residue + 1;

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
    max_bit_length: u32,
    produces_nothing: bool,
    starts_from_exception: bool,
    last_round: bool,
}

impl<const L: usize> Sieve<L> {
    /// Creates a new sieve, iterating from `start` and
    /// until the last number with `max_bit_length` bits,
    /// producing numbers that are not non-trivial multiples
    /// of a list of small primes in the range `[2, start)` (`safe_primes = false`)
    /// or `[2, start/2)` (`safe_primes = true`).
    ///
    /// Note that `start` is adjusted to `2`, or the next `1 mod 2` number (`safe_primes = false`);
    /// and `5`, or `3 mod 4` number (`safe_primes = true`).
    ///
    /// Panics if `max_bit_length` is zero or greater than the size of the target `Uint`.
    ///
    /// If `safe_primes` is `true`, both the returned `n` and `n/2` are sieved.
    pub fn new(start: &Uint<L>, max_bit_length: u32, safe_primes: bool) -> Self {
        if max_bit_length == 0 {
            panic!("The requested bit length cannot be zero");
        }

        if max_bit_length > Uint::<L>::BITS {
            panic!(
                "The requested bit length ({}) is larger than the chosen Uint size",
                max_bit_length
            );
        }

        // If we are targeting safe primes, iterate over the corresponding
        // possible Germain primes (`n/2`), reducing the task to that with `safe_primes = false`.
        let (max_bit_length, base) = if safe_primes {
            (max_bit_length - 1, start.wrapping_shr_vartime(1))
        } else {
            (max_bit_length, *start)
        };

        let mut base = base;

        // This is easier than making all the methods generic enough to handle these corner cases.
        let produces_nothing = max_bit_length < base.bits() || max_bit_length < 2;

        // Add the exception to the produced candidates - the only one that doesn't fit
        // the general pattern of incrementing the base by 2.
        let mut starts_from_exception = false;
        if base <= Uint::<L>::from(2u32) {
            starts_from_exception = true;
            base = Uint::<L>::from(3u32);
        } else {
            // Adjust the base so that we hit odd numbers when incrementing it by 2.
            base |= Uint::<L>::ONE;
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
        // Should not overflow since `incr` is never greater than `incr_limit`,
        // and the latter is chosen such that it doesn't overflow when added to `base`
        // (see the rest of this method).
        self.base = self
            .base
            .checked_add(&self.incr.into())
            .expect("Integer overflow");

        self.incr = 0;

        // Re-calculate residues.
        for (i, rec) in RECIPROCALS.iter().enumerate().take(self.residues.len()) {
            let (_quo, rem) = self.base.div_rem_limb_with_reciprocal(rec);
            self.residues[i] = rem.0 as SmallPrime;
        }

        // Find the increment limit.
        let max_value = Uint::<L>::ONE
            .overflowing_shl(self.max_bit_length)
            .unwrap_or(Uint::ZERO)
            .wrapping_sub(&Uint::<L>::ONE);
        let incr_limit = max_value.wrapping_sub(&self.base);
        self.incr_limit = if incr_limit > Uint::<L>::from(INCR_LIMIT) {
            INCR_LIMIT
        } else {
            // We are close to `2^max_bit_length - 1`.
            // Mark this round as the last.
            self.last_round = true;
            // Can unwrap here since we just checked above that `incr_limit <= INCR_LIMIT`,
            // and `INCR_LIMIT` fits into `Residue`.
            let incr_limit_small: Residue = incr_limit.as_words()[0]
                .try_into()
                .expect("ensured to fit within `Residue`");
            incr_limit_small
        };

        true
    }

    // Returns `true` if the current `base + incr` is divisible by any of the small primes.
    fn current_is_composite(&self) -> bool {
        for (i, m) in self.residues.iter().enumerate() {
            let d = SMALL_PRIMES[i] as Residue;
            let r = (*m as Residue + self.incr) % d;

            if r == 0 {
                return true;
            }

            // A trick from "Safe Prime Generation with a Combined Sieve" by Michael J. Wiener
            // (https://eprint.iacr.org/2003/186).
            // Remember that the check above was for the `(n - 1)/2`;
            // If `(n - 1)/2 mod d == (d - 1)/2`, it means that `n mod d == 0`.
            // In other words, we are checking the remainder of `n mod d`
            // for virtually no additional cost.
            if self.safe_primes && r == (d - 1) >> 1 {
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
            // The overflow should never happen here since `incr`
            // is never greater than `incr_limit`, and the latter is chosen such that
            // it does not overflow when added to `base` (see `update_residues()`).
            let mut num: Uint<L> = self
                .base
                .checked_add(&self.incr.into())
                .expect("Integer overflow");
            if self.safe_primes {
                num = num.wrapping_shl_vartime(1) | Uint::<L>::ONE;
            }
            Some(num)
        };

        self.incr += 2;
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
        Self::next(self)
    }
}

#[cfg(test)]
mod tests {

    use alloc::format;
    use alloc::vec::Vec;

    use crypto_bigint::{Odd, U64};
    use num_prime::nt_funcs::factorize64;
    use rand_chacha::ChaCha8Rng;
    use rand_core::{OsRng, SeedableRng};

    use super::{random_odd_uint, Sieve};
    use crate::hazmat::precomputed::SMALL_PRIMES;

    #[test]
    fn random() {
        let max_prime = SMALL_PRIMES[SMALL_PRIMES.len() - 1];

        let mut rng = ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");
        let start: Odd<U64> = random_odd_uint(&mut rng, 32);
        for num in Sieve::new(&start, 32, false).take(100) {
            let num_u64: u64 = num.into();
            assert!(num_u64.leading_zeros() == 32);

            let factors_and_powers = factorize64(num_u64);
            let factors = factors_and_powers.into_keys().collect::<Vec<_>>();

            assert!(factors[0] > max_prime as u64);
        }
    }

    fn check_sieve(start: u32, bit_length: u32, safe_prime: bool, reference: &[u32]) {
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

        // In the following three cases, the "half-start" would be set to 3,
        // and since every small divisor equal or greater than the start is not tested
        // (because we can't distinguish between the remainder being 0
        // and the number being actually equal to the divisor),
        // no divisors will actually be tested at all, so 15 (a composite)
        // is included in the output.
        check_sieve(1, 4, true, &[5, 7, 11, 15]);
        check_sieve(5, 4, true, &[5, 7, 11, 15]);
        check_sieve(7, 4, true, &[7, 11, 15]);

        check_sieve(9, 4, true, &[11]);
        check_sieve(13, 4, true, &[]);
    }

    #[test]
    #[should_panic(expected = "The requested bit length cannot be zero")]
    fn sieve_zero_bits() {
        let _sieve = Sieve::new(&U64::ONE, 0, false);
    }

    #[test]
    #[should_panic(expected = "The requested bit length (65) is larger than the chosen Uint size")]
    fn sieve_too_many_bits() {
        let _sieve = Sieve::new(&U64::ONE, 65, false);
    }

    #[test]
    fn random_below_max_length() {
        for _ in 0..10 {
            let r: Odd<U64> = random_odd_uint(&mut OsRng, 50);
            assert_eq!(r.bits(), 50);
        }
    }

    #[test]
    #[should_panic(expected = "Bit length must be non-zero")]
    fn random_odd_uint_0bits() {
        let _p: Odd<U64> = random_odd_uint(&mut OsRng, 0);
    }

    #[test]
    #[should_panic(expected = "The requested bit length (65) is larger than the chosen Uint size")]
    fn random_odd_uint_too_many_bits() {
        let _p: Odd<U64> = random_odd_uint(&mut OsRng, 65);
    }

    #[test]
    fn sieve_derived_traits() {
        let s = Sieve::new(&U64::ONE, 10, false);
        assert!(format!("{s:?}").starts_with("Sieve"));
        assert_eq!(s.clone(), s);
    }
}
