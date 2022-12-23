//! An iterator for weeding out multiples of small primes,
//! before proceeding with slower tests.

use alloc::{vec, vec::Vec};

use crypto_bigint::{Integer, Limb, Random, Uint, Zero};
use rand_core::{CryptoRng, RngCore};

use crate::hazmat::precomputed::{SmallPrime, SMALL_PRIMES, SMALL_PRIMES_RECIPROCALS};

/// Returns a random odd integer with given bit length
/// (that is, with both `0` and `bit_length-1` bits set).
pub fn random_odd_uint<const L: usize, R: CryptoRng + RngCore + ?Sized>(
    rng: &mut R,
    bit_length: usize,
) -> Uint<L> {
    if bit_length > Limb::BIT_SIZE * L {
        panic!(
            "The requested bit length ({}) is larger than the chosen Uint size",
            bit_length
        );
    }

    // TODO: not particularly efficient, can be improved by zeroing high bits instead of shifting
    let mut random = Uint::<L>::random(rng);
    if bit_length != Limb::BIT_SIZE * L {
        random >>= Limb::BIT_SIZE * L - bit_length;
    }

    // Make it odd
    random |= Uint::<L>::ONE;

    // Make sure it's the correct bit size
    random |= Uint::<L>::ONE << (bit_length - 1);

    random
}

/// An iterator returning numbers with up to and including given bit length,
/// starting from a given number, that are not multiples of the first 2048 small primes.
#[derive(Clone, Debug)]
pub struct Sieve<const L: usize> {
    // Instead of dividing `Uint` by small primes every time (which is slow),
    // we keep a "base" and a small increment separately,
    // so that we can only calculate the residues of the increment.
    base: Uint<L>,
    incr: Residue,
    incr_limit: Residue,
    residues: Vec<SmallPrime>,
    max_bit_length: usize,
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
    /// until the last number with `max_bit_length` bits.
    pub fn new(start: &Uint<L>, max_bit_length: usize) -> Self {
        // Since 2 is not in `SMALL_PRIMES`, we need the starting number to be pre-selected
        // to not be divisible by it.
        if start.is_even().into() {
            panic!("Sieving must be started from an odd number");
        }

        Self {
            base: *start,
            incr: 0, // This will ensure that `update_residues()` is called right away.
            incr_limit: 0,
            residues: vec![0; SMALL_PRIMES.len()],
            max_bit_length,
            last_round: false,
        }
    }

    fn update_residues(&mut self) -> bool {
        if self.incr >= self.incr_limit {
            if self.last_round {
                return false;
            }

            // Set the new base.
            self.base = self.base.wrapping_add(&self.incr.into());
            self.incr = 0;

            // Re-calculate residues.
            for (i, reciprocal) in SMALL_PRIMES_RECIPROCALS
                .iter()
                .enumerate()
                .take(self.residues.len())
            {
                let (_quo, rem) = self.base.ct_div_rem_limb_with_reciprocal(reciprocal);
                self.residues[i] = rem.0 as SmallPrime;
            }

            // Find the increment limit.
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
        }
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
        self.incr += 2;
        result
    }
}

impl<const L: usize> Iterator for Sieve<L> {
    type Item = Uint<L>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.update_residues() {
            match self.maybe_next() {
                Some(x) => return Some(x),
                None => continue,
            };
        }
        None
    }
}

/// Performs trial division of the given number by a list of small primes.
/// Returns `Some(is_prime)` if there was a definitive conclusion about `num`'s primality,
/// and `None` otherwise.
pub fn sieve_once<const L: usize>(num: &Uint<L>) -> Option<bool> {
    // Our reciprocals start from 3, so we check for 2 separately
    if num == &Uint::<L>::from(2u32) {
        return Some(true);
    }
    if num.is_even().into() {
        return Some(false);
    }
    for reciprocal in SMALL_PRIMES_RECIPROCALS.iter() {
        let (quo, rem) = num.ct_div_rem_limb_with_reciprocal(reciprocal);
        if rem.is_zero().into() {
            return Some(quo == Uint::<L>::ONE);
        }
    }
    None
}

#[cfg(test)]
mod tests {

    use alloc::vec::Vec;

    use number_theory::NumberTheory;
    use rand_chacha::ChaCha8Rng;
    use rand_core::SeedableRng;

    use super::{random_odd_uint, Sieve};
    use crate::hazmat::precomputed::SMALL_PRIMES;

    #[test]
    fn random() {
        let max_prime = SMALL_PRIMES[SMALL_PRIMES.len() - 1];

        let mut rng = ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");
        let start = random_odd_uint::<1, _>(&mut rng, 32);
        for num in Sieve::new(&start, 32).take(100) {
            let num_u64: u64 = num.into();
            assert!(num_u64.leading_zeros() == 32);

            let factors_and_powers = num_u64.factor();
            let mut factors = factors_and_powers
                .chunks(2)
                .map(|pair| pair[0])
                .collect::<Vec<_>>();
            factors.sort_unstable();

            assert!(factors[0] > max_prime as u64);
        }
    }
}
