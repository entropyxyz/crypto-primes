//! Miller-Rabin primality test; the numbers that pass it
//! are also called "strong probable-/pseudoprimes".
//!
//! See C. Pomerance, J. L. Selfridge, S. S. Wagstaff "The Pseudoprimes to 25*10^9",
//! Math. Comp. 35 1003-1026 (1980) (DOI: [10.2307/2006210](https://dx.doi.org/10.2307/2006210))

use rand_core::{CryptoRng, RngCore};

use crypto_bigint::{
    modular::{
        runtime_mod::{DynResidue, DynResidueParams},
        MulResidue, PowResidue,
    },
    Integer, NonZero, RandomMod, Uint,
};

pub struct MillerRabin<const L: usize> {
    candidate: Uint<L>,
    montgomery_params: DynResidueParams<L>,
    one_m: DynResidue<L>,
    cand_minus_one_m: DynResidue<L>,
    s: u32,
    d: Uint<L>,
}

impl<const L: usize> MillerRabin<L> {
    pub fn new(candidate: &Uint<L>) -> Self {
        let params = DynResidueParams::<L>::new(*candidate);
        let one_m = DynResidue::<L>::new(Uint::<L>::ONE, params);
        let cand_minus_one = candidate.wrapping_sub(&Uint::<L>::ONE);
        let cand_minus_one_m = DynResidue::<L>::new(cand_minus_one, params);
        let (s, d) = decompose(candidate);
        Self {
            candidate: *candidate,
            montgomery_params: params,
            one_m,
            cand_minus_one_m,
            s,
            d,
        }
    }

    pub fn check(&self, basis: &Uint<L>) -> bool {
        // TODO: it may be faster to first check that gcd(basis, candidate) == 1,
        // otherwise we can return `false` right away.

        let basis_m = DynResidue::<L>::new(*basis, self.montgomery_params);
        let mut test_m = basis_m.pow(&self.d);

        if test_m == self.one_m || test_m == self.cand_minus_one_m {
            return true;
        }
        for _ in 1..self.s {
            test_m = test_m.square();
            if test_m == self.one_m {
                return false;
            } else if test_m == self.cand_minus_one_m {
                return true;
            }
        }
        false
    }

    pub fn check_basis_two(&self) -> bool {
        self.check(&Uint::<L>::from(2u32))
    }

    pub fn check_random_basis<R: CryptoRng + RngCore>(&self, rng: &mut R) -> bool {
        // We sample a random basis from the range `[3, candidate-2]`:
        // - we have a separate method for basis 2;
        // - the test holds trivially for bases 1 or `candidate-1`.
        let range = self.candidate.wrapping_sub(&Uint::<L>::from(4u32));
        let range_nonzero = NonZero::new(range).unwrap();
        let random =
            Uint::<L>::random_mod(rng, &range_nonzero).wrapping_add(&Uint::<L>::from(3u32));
        self.check(&random)
    }
}

/// For the given odd `n`, finds `s` and odd `d` such that `n - 1 == 2^s * d`.
fn decompose<const L: usize>(n: &Uint<L>) -> (u32, Uint<L>) {
    let mut d = n.wrapping_sub(&Uint::<L>::ONE);
    let mut s = 0;

    while d.is_even().into() {
        d >>= 1;
        s += 1;
    }

    (s, d)
}

#[cfg(test)]
mod tests {
    use crypto_bigint::{Uint, Wrapping};
    use rand_chacha::ChaCha8Rng;
    use rand_core::SeedableRng;

    use super::MillerRabin;
    use crate::hazmat::{pseudoprimes, random_odd_uint, Sieve};

    fn test_sequence(numbers: &[u32]) {
        let mut rng = ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");
        for num in numbers.iter() {
            let mr = MillerRabin::new(&Uint::<1>::from(*num));
            assert!(!mr.check_basis_two(), "{}", num);
            for _ in 0..10 {
                assert!(!mr.check_random_basis(&mut rng));
            }
        }
    }

    #[test]
    fn trivial() {
        let mut rng = ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");
        let start = random_odd_uint::<16, _>(&mut rng, 1024);
        for num in Sieve::new(&start, 1024).take(10) {
            let mr = MillerRabin::new(&num);

            // Trivial tests, must always be true.
            assert!(mr.check(&1u32.into()));
            assert!(mr.check(&num.wrapping_sub(&1u32.into())));
        }
    }

    #[test]
    fn mersenne_prime() {
        let mut rng = ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");

        // Mersenne prime 2^127-1
        let num = Uint::<2>::from_be_hex("7fffffffffffffffffffffffffffffff");

        let mr = MillerRabin::new(&num);
        assert!(mr.check_basis_two());
        for _ in 0..10 {
            assert!(mr.check_random_basis(&mut rng));
        }
    }

    #[test]
    fn strong_fibonacci_pseudoprimes() {
        let mut rng = ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");

        for num in pseudoprimes::STRONG_FIBONACCI.iter() {
            let num_hi = num >> 32;
            let num_lo = num & ((1 << 32) - 1);
            let uint2 = Uint::<2>::from(num_lo) | (Uint::<2>::from(num_hi) << 32);
            let mr = MillerRabin::new(&uint2);
            assert!(!mr.check_basis_two());
            for _ in 0..1000 {
                assert!(!mr.check_random_basis(&mut rng));
            }
        }
    }

    #[test]
    fn strong_pseudoprimes_base_2() {
        // These are composites, but the basis 2 MR test will pass for them.
        for psp in pseudoprimes::STRONG_BASE_2 {
            let num = Uint::<2>::from(*psp);
            let mr = MillerRabin::new(&num);
            assert!(mr.check_basis_two());
        }
    }

    #[test]
    fn extra_strong_lucas_pseudoprimes() {
        // Cross-test against the pseudoprimes that circumvent the Lucas test.
        // We expect the MR test to correctly classify them as composites.
        test_sequence(pseudoprimes::EXTRA_STRONG_LUCAS);
    }

    #[test]
    fn carmichael_number() {
        let p = Wrapping(Uint::<21>::from_be_hex(concat!(
            "0000000000000000",
            "00000000000000000000000000000000",
            "00000000000000000000000000000000",
            "00000000000000000000000000000000",
            "00000000000000000000000000000000",
            "00000000000000000000000000000000",
            "00000000000000000000000000000000",
            "00000000000000000002acf5bff060d5",
            "cb98a808aef25df948beef3f7c06c5c5",
            "c9efed1e53d19363db74537f89d3cadc",
            "c9be464335664ed21d0b6dbc736e0e23",
        )));
        let c1 = Wrapping(Uint::<21>::ONE);
        let c313 = Wrapping(Uint::<21>::from(313u64));
        let c353 = Wrapping(Uint::<21>::from(353u64));

        // This is a known Carmichael number - a composite that passes MR test with some bases.
        let n = p * ((p - c1) * c313 + c1) * ((p - c1) * c353 + c1);

        let mr = MillerRabin::new(&n.0);

        // It is known to pass MR tests for all prime bases <307
        assert!(mr.check_basis_two());
        assert!(mr.check(&Uint::<21>::from(293u64)));

        // A test with basis 307 correctly reports the number as composite.
        assert!(!mr.check(&Uint::<21>::from(307u64)));
    }
}
