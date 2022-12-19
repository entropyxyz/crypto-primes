use crypto_bigint::Uint;
use rand_core::{CryptoRng, OsRng, RngCore};

use crate::hazmat::{is_strong_lucas_prime, random_odd_uint, MillerRabin, SelfridgeBase, Sieve};

pub fn prime<const L: usize>(bit_length: usize) -> Uint<L> {
    prime_with_rng(&mut OsRng, bit_length)
}

pub fn is_prime<const L: usize>(rng: &mut (impl CryptoRng + RngCore), num: &Uint<L>) -> bool {
    let mr = MillerRabin::new(num);
    if !mr.check_basis_two() {
        return false;
    }
    if !is_strong_lucas_prime::<SelfridgeBase, L>(num, true) {
        return false;
    }
    if !mr.check_random_basis(rng) {
        return false;
    }
    return true;
}

pub fn prime_with_rng<const L: usize>(
    rng: &mut (impl CryptoRng + RngCore),
    bit_length: usize,
) -> Uint<L> {
    loop {
        let start = random_odd_uint::<L, _>(rng, bit_length);
        let sieve = Sieve::new(&start, bit_length);
        for num in sieve {
            if is_prime(rng, &num) {
                return num;
            }
        }
    }
}
