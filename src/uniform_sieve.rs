// Actors in this play are:
// n:	the number of bits in the prime we're looking for, e.g. 512
// l:	the number of top bits that are re-sampled on every iteration, e.g. 64
// m:	a product of all small odd primes up to a bound ß chosen such that the bit size of m is n - l, e.g. ß = 512 - 64
//		For n=64, we can fit the product of the first 15 odd primes and get 16294579238595022365.
//		For n=128, we can fit the product of the first 24 odd primes and get 1152783981972759212376551073665878035.
// b:	picked uniformely at random among integers less than m and coprime to m (use unit generation algorithm from JP06)
// λ:	Carmichael's function: the LCM of the λ(p)s for each prime used to compute m. Each λ(p) is simply p-1 because each prime appears once and we exclude 2.

// fn lambda(m: u64) -> u64 {}
use crate::hazmat::precomputed::SMALL_PRIMES;
use crate::{hazmat::binary_gcd, is_prime};

use crypto_bigint::modular::Retrieve;
use crypto_bigint::{
    modular::{MontyForm, MontyParams},
    Integer, Monty, Odd, U1024, U256,
};
use crypto_bigint::{BitOps, Bounded, Constants, FixedInteger, NonZero, PowBoundedExp, RandomBits, RandomMod, U128};
use rand::{thread_rng, Rng};

fn jp09_unitgen_draft(m: u64, lambda_m: u32) -> u64 {
    let mut rng = thread_rng();
    // 1. sample k in [1, m[
    let mut k = rng.gen_range(1..m);
    // 2. set `u = 1-(k^lambda_m) mod m`; rewrite as `(1 - k^lambda_m + m) mod m` if k^lambda_m < m, else compute (k^lambda_m) mod m
    let mut k_lambda_m = (k as u128).checked_pow(lambda_m).expect("k^lambda_m exploded");
    if k_lambda_m > m as u128 {
        k_lambda_m = k_lambda_m % m as u128;
    }
    let mut k_lambda_m = k_lambda_m as u64;
    let mut u = (1 + m - k_lambda_m) % m;

    // if u != 0 {
    // 	pick random r from [1,m[
    // 	k = k + r*u mod m
    // 	go to step 2
    // }
    while u != 0 {
        let r = rng.gen_range(1..m);
        k = (k + r * u) % m;
        k_lambda_m = k.checked_pow(lambda_m).expect("k^lmbda_m exploded in loop");
        if k_lambda_m > m {
            k_lambda_m = k_lambda_m % m;
        }
        u = (1 + m - k_lambda_m) % m;
    }
    // return k
    k
}

fn jp09_unitgen_draft2(m: u64, lambda_m: u32) -> U1024 {
    let mut rng = thread_rng();
    // 1. sample k in [1, m[
    let k = rng.gen_range(1..m);
    let k = U1024::from(k);
    // 2. set `u = 1-(k^lambda_m) mod m`; rewrite as `(1 - k^lambda_m + m) mod m` if k^lambda_m < m, else compute (k^lambda_m) mod m
    let m_big = U1024::from(m);
    let lambda_m_big = U1024::from(lambda_m);
    let prms = <U1024 as Integer>::Monty::new_params_vartime(Odd::new(m_big).expect("m is odd"));
    let one = MontyForm::one(prms);
    let zero = MontyForm::zero(prms);

    let mut k = MontyForm::new(&k, prms);
    let mut u = one - k.pow_bounded_exp(&lambda_m_big, 32);

    println!("[jp09] Entering loop with u={:?}", u);
    // if u != 0 {
    // 	pick random r from [1,m[
    // 	k = k + r*u mod m
    // 	go to step 2
    // }
    while u != zero {
        let r = U1024::from(rng.gen_range(1..m));
        let r = MontyForm::new(&r, prms);
        k = k + r * u;
        u = one - k.pow_bounded_exp(&lambda_m_big, 32);
    }
    k.retrieve()
}

fn jp09_unitgen<T>(m: u64, lambda_m: u32) -> T
where
    T: Integer,
    T::Monty: Retrieve<Output = T> + Copy,
    <<T as Integer>::Monty as Monty>::Params: Copy,
{
    let mut rng = thread_rng();
    // 1. sample k in [1, m[
    let k = rng.gen_range(1..m);
    let k = T::from(k);
    // 2. set `u = 1-(k^lambda_m) mod m`; rewrite as `(1 - k^lambda_m + m) mod m` if k^lambda_m < m, else compute (k^lambda_m) mod m
    let m_big = T::from(m);
    let lambda_m_big = T::from(lambda_m);
    let prms = <T as Integer>::Monty::new_params_vartime(Odd::new(m_big).expect("m is odd"));
    let one = T::Monty::one(prms);
    let zero = T::Monty::zero(prms);

    let mut k = T::Monty::new(k, prms);
    let mut u = one - k.pow_bounded_exp(&lambda_m_big, 32);

    // if u != 0 {
    // 	pick random r from [1,m[
    // 	k = k + r*u mod m
    // 	go to step 2
    // }
    while u != zero {
        let r = T::from(rng.gen_range(1..m));
        let r = T::Monty::new(r, prms);
        k = k + r * u;
        u = one - k.pow_bounded_exp(&lambda_m_big, 32);
    }
    k.retrieve()
}

fn algo2_draft(b: U1024, m: u64) -> U1024 {
    let mut rng = thread_rng();
    let m_big = U1024::from(m);
    let a_max = U1024::BITS - m.leading_zeros();
    println!("[algo2]\n\tb:{:?}\n\tm:{:?}\n\ta_max:{:?}", b, m_big, a_max);
    let a = U1024::random_bits(&mut rng, a_max);
    let mut p = a * m_big + b;
    while !is_prime(&p) {
        let a = U1024::random_bits(&mut rng, a_max);
        println!("[algo2, loop] a: {:?}", a);
        p = a * m_big + b;
    }
    p
}

fn algo2<T>(b: T, m: u64) -> T
where
    T: Integer + Constants + Bounded + RandomBits + RandomMod + Copy,
{
    let mut rng = thread_rng();
    let m_big = T::from(m);
    let a_max = T::BITS - u64::BITS + m.leading_zeros();
    // println!(
    //     "[algo2]\n\tb:{:?}\n\tm:{:?}\n\ta_max:{:?}\n\tm.leading_zeros:{:?}",
    //     b,
    //     m_big,
    //     a_max,
    //     m.leading_zeros()
    // );
    let a = T::random_bits(&mut rng, a_max);
    let mut p = a * m_big + b;
    while !is_prime(&p) {
        let a = T::random_bits(&mut rng, a_max);
        // println!("[algo2, loop] a: {:?}", a);
        p = a * m_big + b;
    }
    p
}

/// Bla bla TODO(dp)
pub fn generate_prime<T>() -> T
where
    T: FixedInteger + RandomBits + RandomMod,
    // T: Integer + Constants + Bounded + RandomBits + RandomMod + Copy,
    T::Monty: Retrieve<Output = T> + Copy,
    <<T as Integer>::Monty as Monty>::Params: Copy,
{
    let (m, lambda_m) = right_consts();
    let unit = jp09_unitgen(m, lambda_m);
    // let unit = jp09_unitgen(M_32 as u64, LAMBDA_M_32 as u32);
    // let unit = jp09_unitgen(M_64, LAMBDA_M_64);
    algo2(unit, M_64)
}

fn bla128() -> (u128, u64) {
    // The product of the 24 first primes\2 fits in a u128
    // Primes: 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97
    let (m, lambda_m) = SMALL_PRIMES[..24].iter().fold((1u128, 1u128), |mut acc, p| {
        acc.0 = acc.0.checked_mul(*p as u128).unwrap();
        let lcm = {
            let p_minus_one: u128 = *p as u128 - 1;
            let gcd = binary_gcd128(acc.1, p_minus_one);
            acc.1 * (*p as u128 / gcd)
        };
        acc.1 = lcm;
        acc
    });
    (m, lambda_m.try_into().expect("lambda_m does not fit into a u64"))
}

fn bla64() -> (u64, u32) {
    // The product of the 15 first primes\2 fits in a u64
    // Primes: 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
    let (m, lambda_m) = SMALL_PRIMES[..15].iter().fold((1u64, 1u64), |mut acc, p| {
        acc.0 = acc.0.checked_mul(*p as u64).unwrap();
        let lcm = {
            let p_minus_one: u64 = *p as u64 - 1;
            let gcd = binary_gcd(acc.1, p_minus_one);
            acc.1 * (*p as u64 / gcd)
        };

        acc.1 = lcm;
        acc
    });
    (m, lambda_m.try_into().expect("lambda_m does not fit into a u32"))
}

fn bla32() -> (u32, u16) {
    // The product of the 9 first odd primes fits in a u32
    // Primes: 3, 5, 7, 11, 13, 17, 19, 23, 29
    let (m, lambda_m) = SMALL_PRIMES[..9].iter().fold((1u32, 1u64), |mut acc, p| {
        acc.0 = acc.0.checked_mul(*p as u32).unwrap();
        let lambda_p = {
            let p_minus_one: u64 = *p as u64 - 1;
            let gcd = binary_gcd(acc.1, p_minus_one);
            acc.1 * (*p as u64 / gcd)
        };

        acc.1 = lambda_p;
        acc
    });
    (m, lambda_m.try_into().expect("lambda_m does not fit into a u16"))
}

fn binary_gcd128(mut n: u128, mut m: u128) -> u128 {
    // Using identities 2 and 3:
    // gcd(2ⁱn, 2ʲm) = 2ᵏ gcd(n, m) with n, m odd and k = min(i, j)
    // 2ᵏ is the greatest power of two that divides both 2ⁱn and 2ʲm
    let i = n.trailing_zeros();
    n >>= i;
    let j = m.trailing_zeros();
    m >>= j;
    let k = core::cmp::min(i, j);

    loop {
        // Swap if necessary so n ≤ m
        if n > m {
            core::mem::swap(&mut n, &mut m);
        }

        // Identity 4: gcd(n, m) = gcd(n, m-n) as n ≤ m and n, m are both odd
        m -= n;
        // m is now even

        if m == 0 {
            // Identity 1: gcd(n, 0) = n
            // The shift by k is necessary to add back the 2ᵏ factor that was removed before the loop
            return n << k;
        }

        // Identity 3: gcd(n, 2ʲ m) = gcd(n, m) as n is odd
        m >>= m.trailing_zeros();
    }
}

const M_128: u128 = 1152783981972759212376551073665878035;
// const LAMBDA_M_128: u64 = 39419059680u64;
const LAMBDA_M_128: u128 = 39419059680u128;
const M_64: u64 = 16294579238595022365;
const LAMBDA_M_64: u32 = 16576560;
const M_32: u32 = 3234846615;
const LAMBDA_M_32: u16 = 55440;

// const ARY: [(u128, u64); 3] = [
//     (M_128, LAMBDA_M_128),
//     (M_64 as u128, LAMBDA_M_64 as u64),
//     (M_32 as u128, LAMBDA_M_32 as u64),
// ];

trait Bla<MT> {
    const M: MT;
    const LAMBDA_M: MT;
}

impl Bla<U128> for U256 {
    const M: U128 = U128::from_u128(M_128);
    const LAMBDA_M: U128 = U128::from_u128(LAMBDA_M_128);
}

struct Beef<T> {
    m: T,
    lambda_m: T,
}

impl<T: From<u128>> Beef<T> {
    fn new() -> Self {
        Self {
            m: T::from(M_128),
            lambda_m: T::from(LAMBDA_M_128),
        }
    }
}

impl From<U128> for Beef<U128> {
    fn from(value: U128) -> Self {
        Self {
            m: U128::from_u128(M_128),
            lambda_m: U128::from_u128(LAMBDA_M_128),
        }
    }
}

fn right_consts<T: Integer>() -> (T, T) {
    match T::LIMBS {}
}

const A_MAX_M_32: u32 = 992;
#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::U1024;
    #[test]
    fn unitgen64() {
        let unit = jp09_unitgen(M_64, LAMBDA_M_64);
        println!("Unit={:?}", unit);
        let p = algo2::<U1024>(unit, M_64);
        println!("prime={:?}", p);
    }
    #[test]
    fn unitgen() {
        let unit = jp09_unitgen(M_32 as u64, LAMBDA_M_32 as u32);
        println!("Unit={:?}", unit);
        let p = algo2::<U1024>(unit, M_32 as u64);
        println!("prime={:?}", p);
    }

    #[test]
    fn lcm_works() {
        fn lcm(a: u64, b: u64) -> u64 {
            a * (b / binary_gcd(a, b))
        }
        assert_eq!(lcm(2, 4), 4);
        let numbers = vec![2u64, 4, 6, 10, 12];
        let lcm = numbers.into_iter().reduce(|a, b| lcm(a, b)).unwrap();
        assert_eq!(lcm, 60);
    }

    #[test]
    fn primes_in_a_u32() {
        let o = bla32();
        println!("o={:?}", o);
        assert!(o.0 == 3234846615); // This is `m`
        assert_eq!(o.1, 55440u16); // This is `lambda(m)`, aka Carmichael's function
    }
    #[test]
    fn primes_in_a_u64() {
        let o = bla64();
        println!("o={:?}", o);
        assert!(o.0 == 16294579238595022365); // This is `m`
        assert_eq!(o.1, 16576560u32); // This is `lambda(m)`, aka Carmichael's function
    }
    #[test]
    fn primes_in_a_u128() {
        let o = bla128();
        println!("o={:?}", o);
        assert_eq!(o.0, 1152783981972759212376551073665878035); // m
        assert_eq!(o.1, 39419059680u64); // This is `lambda(m)`, aka Carmichael's function
    }
}
