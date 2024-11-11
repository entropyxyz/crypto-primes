use crypto_bigint::modular::Retrieve;
use crypto_bigint::{
    Bounded, Constants, FixedInteger, Integer, Monty, NonZero, Odd, PowBoundedExp, RandomBits, RandomMod, U1024, U128,
    U2048, U256, U4096, U512, U64,
};
use rand::thread_rng;

use crate::hazmat::precomputed::SMALL_PRIMES;
use crate::{hazmat::binary_gcd, is_prime};
/// Generate primes using a pseudo-uniform distribution.
/// Actors in this play are:
///     n:	the number of bits in the prime we're looking for, e.g. 512
///     l:	the number of top bits that are re-sampled on every iteration, e.g. 64
///     m:	a product of all small odd primes up to a bound ß chosen such that the bit size of m is n - l, e.g. ß = 512 - 64
///     b:	picked uniformely at random among integers less than m and coprime to m (use unit generation algorithm from JP06)
///     λ:	Carmichael's function: the LCM of the λ(p)s for each prime used to compute m. Each λ(p) is simply p-1 because each prime appears once and we exclude 2.

// TODO(dp): Proper docs here, explaining when this is useful and why, discussion about the distribution quality etc.
// TODO(dp): The `l` value isn't used anywhere and instead I'm assuming we want to re-sample
// `T::BITS/2` bits on every iteration, so for a 128-bit prime, resample 64 bits. I think this is
// way too much, but need more benchmarking + statistics to know what a "good" value actually means.
pub trait UniformGeneratePrime<T>
where
    T: Integer + Constants + Bounded + RandomBits + RandomMod + Copy,
{
    /// Generate a prime.
    fn generate_prime() -> T;
    // TODO(dp): missing
    // generate_prime_with_rng
    // generate_safe_prime
    // generate_safe_prime_with_rng
}

macro_rules! impl_generate_prime {
    ($(($name:ident, $bits:expr, $m:expr, $lambda_m:expr)),+) => {
        $(
            impl UniformGeneratePrime<$name> for $name {
                fn generate_prime() -> $name {
                    debug_assert!($m.len() == (2*$name::BITS/8) as usize, "expected m to be {} long, instead it's {}", 2*$name::BITS/8, $m.len());
                    const M: $name = $name::from_be_hex($m);
                    const LAMBDA_M: $name = $name::from_be_hex($lambda_m);
                    let unit = jp06_unitgen(M, LAMBDA_M);
                    algorithm2(unit, M)
                }
            }
        )+

    };
}

impl_generate_prime! {
    (U64, 64, "00000000C0CFD797", "000000000000D890")
    ,
    (U128, 128, "0000000000000000E221F97C30E94E1D", "00000000000000000000000000FCF030")
    ,
    (U256, 256, "000000000000000000000000000000005797D47C51681549D734E4FC4C3EAF7F", "0000000000000000000000000000000000000000000000000000002DE3CB9560")
    ,
    (U512, 512, "0000000000000000000000000000000000000000000000000000000000000000DBF05B6F5654B3C0F5243551439586889F155887819AED2AC05B93352BE98677", "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000097262FAE1826B96DBE40")
    ,
    (U1024, 1024, "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000106AA9FB7646FA6EB0813C28C5D5F09F077EC3BA238BFB99C1B631A203E81187233DB117CBC384056EF04659A4A11DE49F7ECB29BADA8F980DECECE92E30C48F", "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000468A4F19CEA99B05F593F7F03EBA946B4AF700")
    ,
    (U2048, 2048, "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002C85FF870F24BE80F62B1BA6C20BD72B837EFDF121206D87DB56B7D69FA4C021C107C3CA206FE8FA7080EF576EFFC82F9B10F5750656B7794B16AFD70996E91AEF6E0AD15E91B071AC9B24D98B233AD86EE055518E58E56638EF18BAC5C74CB35BBB6E5DAE2783DD1C0CE7DEC4FC70E5186D411DF36368F061AA36011F30179", "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000009220E37A82CBB6007A9DEE07DE852B1FD11D7594688264F7F40E71355F33B7EBFC100")
    ,
    (U4096, 4096, "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002465A7BD85011E1C9E0527929FFF268C82EF7EFA416863BAA5ACDB0971DBA0CCAC3EE4999345029F2CF810B99E406AAC5FCE5DD69D1C717DAEA5D18AB913F456505679BC91C57D46D9888857862B36E2EDE2E473C1F0AB359DA25271AFFE15FF240E299D0B04F4CD0E4D7C0E47B1A7BA007DE89AAE848FD5BDCD7F9815564EB060AE14F19CB50C291F0BBD8ED1C4C7F8FC5FBA51662001939B532D92DAC844A8431D400C832D039F5F900B278A75219C2986140C79045D7759540854C31504DC56F1DF5EEBE7BEE447658B917BF696D6927F2E2428FBEB340E515CB9835D63871BE8BBE09CF13445799F2E67788151571A93B4C1EEE55D1B9072E0B2F5C4607F", "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000099F002D0502D70DF7B64150D19B477A781987C56EBCDB7637842794D22D0D022010DC9BA1DB141FC9D07A9587E7A130D9DF6F9B67812800D00")
}

// From https://marcjoye.github.io/papers/JP06pgen.pdf, page 6, fig. 4
// Parameters: m and λ(m)
// Output: a random unit k ∈ (Z/mZ)∗
// 1. Randomly choose k ∈ [1, m[
// 2. Set U ← (1 − kλ(m)) mod m
// 3. If (U ! = 0) then
//  (a) Choose a random r ∈ [1, m[
//  (b) Set k ← k + rU (mod m)
//  (c) Go to Step 2
// 4. Output k
// TODO(dp): probably unify this with "algorithm2" yeah?
#[inline(always)]
fn jp06_unitgen<T>(m: T, lambda_m: T) -> T
where
    T: FixedInteger + RandomBits,
    T::Monty: Retrieve<Output = T> + Copy,
    <<T as Integer>::Monty as Monty>::Params: Copy,
{
    let mut rng = thread_rng();
    // 1. sample k in [1, m[
    // TODO(dp): removing one bit is sketchy, do I have to loop here to get a k smaller than m?
    let k = T::random_bits(&mut rng, m.bits_vartime() - 1);
    debug_assert!(k < m, "k must be sampled in [1, m[, but {:?}>={:?}", k, m);

    // 2. set `u = 1-(k^lambda_m) mod m`; rewrite as `(1 - k^lambda_m + m) mod m` if k^lambda_m < m, else compute (k^lambda_m) mod m
    let prms = <T as Integer>::Monty::new_params_vartime(Odd::new(m).expect("m is odd"));
    let one = T::Monty::one(prms);
    let zero = T::Monty::zero(prms);

    let mut k = T::Monty::new(k, prms);
    let mut u = one - k.pow_bounded_exp(&lambda_m, 32);

    // if u != 0 {
    // 	pick random r from [1,m[
    // 	k = k + r*u mod m
    // 	go to step 2
    // }
    while u != zero {
        // TODO(dp): removing one bit is sketchy, do I have to loop here to get an r smaller than m?
        let r = T::random_bits(&mut rng, m.bits_vartime() - 1);
        debug_assert!(r < m, "r must be sampled in [1, m[, but {:?}>={:?}", r, m);
        let r = T::Monty::new(r, prms);
        k = k + r * u;
        u = one - k.pow_bounded_exp(&lambda_m, 32);
    }
    k.retrieve()
}

// From https://eprint.iacr.org/2011/481.pdf, page 4:
// Algorithm 2 More eﬃcient method
// 1:   b $← {1,…, m−1}
// 2:   u ← (1−b^λ(m)) mod m
// 3:   if u != 0 then
// 4:       r $← {1,…, m−1}
// 5:       b ← b + ru mod m
// 6:       goto step 2
// 7:   end if
// 8:   repeat
// 9:       a $← {0,…, ⌊2n/m⌋−1}
// 10:      p ← am + b
// 11:  until p is prime
// 12:  return p
// NOTE: Steps 1-7 are implemented in `jp06_unitgen`, where `b` is referred to as `k`.
// TODO(dp): probably should unify the two functions.
#[inline(always)]
fn algorithm2<T>(b: T, m: T) -> T
where
    T: Integer + Constants + Bounded + RandomBits + RandomMod + Copy,
{
    let mut rng = thread_rng();
    let a_max = T::BITS - m.leading_zeros();
    let a = T::random_bits(&mut rng, a_max);
    let mut p = a * m + b;
    while !is_prime(&p) {
        let a = T::random_bits(&mut rng, a_max);
        p = a * m + b;
    }
    p
}

// Calculates two constants, `m` and `λ(m)`, used in the uniform sieving algorithms.
// `m` is the product of all odd primes such that the product fits in a `T`.
// `λ(m)`, aka Carmichael's function is the LCM of the λ(p)s for each prime used to compute m. In
// our case, each λ(p) is simply p-1 because each prime appears once (and we exclude 2).
// Constant values for common `T`s sized from 32 to 2028:
// For u32, 9 primes,       m = 0xC0CFD797
//                          λ(m) = 0x0000D890
// For U64, 15 primes,      m = 0xE221F97C30E94E1D
//                          λ(m) = 0x0000000000FCF030
// For U128, 25 primes,     m = 0x5797D47C51681549D734E4FC4C3EAF7F
//                          λ(m) = 0x00000000000000000000002DE3CB9560
// For U256, 43 primes,     m = 0xDBF05B6F5654B3C0F5243551439586889F155887819AED2AC05B93352BE98677
//                          λ(m) = 0x0000000000000000000000000000000000000000000097262FAE1826B96DBE40
// For U512, 74 primes,     m = 0x106AA9FB7646FA6EB0813C28C5D5F09F077EC3BA238BFB99C1B631A203E81187233DB117CBC384056EF04659A4A11DE49F7ECB29BADA8F980DECECE92E30C48F
//                          λ(m) = 0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000468A4F19CEA99B05F593F7F03EBA946B4AF700
// For U1024, 130 primes,   m = 0x02C85FF870F24BE80F62B1BA6C20BD72B837EFDF121206D87DB56B7D69FA4C021C107C3CA206FE8FA7080EF576EFFC82F9B10F5750656B7794B16AFD70996E91AEF6E0AD15E91B071AC9B24D98B233AD86EE055518E58E56638EF18BAC5C74CB35BBB6E5DAE2783DD1C0CE7DEC4FC70E5186D411DF36368F061AA36011F30179
//                          λ(m) = 0x00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000009220E37A82CBB6007A9DEE07DE852B1FD11D7594688264F7F40E71355F33B7EBFC100
// For U2048, 232 primes,   m = 0x2465A7BD85011E1C9E0527929FFF268C82EF7EFA416863BAA5ACDB0971DBA0CCAC3EE4999345029F2CF810B99E406AAC5FCE5DD69D1C717DAEA5D18AB913F456505679BC91C57D46D9888857862B36E2EDE2E473C1F0AB359DA25271AFFE15FF240E299D0B04F4CD0E4D7C0E47B1A7BA007DE89AAE848FD5BDCD7F9815564EB060AE14F19CB50C291F0BBD8ED1C4C7F8FC5FBA51662001939B532D92DAC844A8431D400C832D039F5F900B278A75219C2986140C79045D7759540854C31504DC56F1DF5EEBE7BEE447658B917BF696D6927F2E2428FBEB340E515CB9835D63871BE8BBE09CF13445799F2E67788151571A93B4C1EEE55D1B9072E0B2F5C4607F
//                          λ(m) = 0x0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000099F002D0502D70DF7B64150D19B477A781987C56EBCDB7637842794D22D0D022010DC9BA1DB141FC9D07A9587E7A130D9DF6F9B67812800D00
// TODO(dp): Would love for this to be `const`. Or not here at all. Does it have any value, in some form?
#[allow(unused)]
fn calculate_m_and_lambda_m<T>() -> (T, T)
where
    T: FixedInteger + crypto_bigint::Gcd<Output = T>,
{
    let mut m = T::ONE;
    let mut lambda_m = T::ONE;
    for (i, prime) in SMALL_PRIMES.into_iter().enumerate() {
        let prime_t = T::from(prime);
        let prod_overflowed = m.checked_mul(&prime_t);
        if prod_overflowed.is_none().into() {
            println!("Breaking at i={i}");
            break;
        } else {
            m = prod_overflowed.unwrap();
            let p_minus_one = prime_t - T::ONE;
            let gcd = lambda_m.gcd(&p_minus_one);
            let gcd_nz = NonZero::new(gcd).unwrap();
            lambda_m = lambda_m * prime_t.div(gcd_nz);
        }
    }

    (m, lambda_m)
}

// Special case for u32s, given crypto-bigint does not provide such small uints.
#[allow(unused)]
fn calculate_m_and_lambda_m_32() -> (u32, u16) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::{U1024, U128, U2048, U256, U512, U64};
    #[test]
    fn generate_m_and_lambda_m() {
        let (m, lambda_m) = calculate_m_and_lambda_m_32();
        println!("32 bits: m={m:X?}, lambda_m={lambda_m:X?}");

        let (m, lambda_m) = calculate_m_and_lambda_m::<U64>();
        println!("64 bits: m={m:?}, lambda_m={lambda_m:?}");

        let (m, lambda_m) = calculate_m_and_lambda_m::<U128>();
        println!("128 bit: m={m:?}, lambda_m={lambda_m:?}");

        let (m, lambda_m) = calculate_m_and_lambda_m::<U256>();
        println!("256 bits: m={m:?}, lambda_m={lambda_m:?}");

        let (m, lambda_m) = calculate_m_and_lambda_m::<U512>();
        println!("512 bit: m={m:?}, lambda_m={lambda_m:?}");

        let (m, lambda_m) = calculate_m_and_lambda_m::<U1024>();
        println!("1024 bits: m={m:?}, lambda_m={lambda_m:?}");

        let (m, lambda_m) = calculate_m_and_lambda_m::<U2048>();
        println!("2048 bit: m={m:?}, lambda_m={lambda_m:?}");
    }
    #[test]
    fn genprime() {
        // This isn't actually a test, just here to generate constants
        let p: U64 = U64::generate_prime();
        println!("64 bit prime={p:?}");
        let p = U128::generate_prime();
        println!("128 bit prime={p:?}");
        let p = U256::generate_prime();
        println!("256 bit prime={p:?}");
        let p = U512::generate_prime();
        println!("512 bit prime={p:?}");
        let p = U1024::generate_prime();
        println!("1024 bit prime={p:?}");
        let p = U2048::generate_prime();
        println!("2048 bit prime={p:?}");
        // This is very slow
        // let p = U4096::generate_prime();
        // println!("4096 bit prime={p:?}");
    }

    // TODO(dp): test for statisical properties
}
