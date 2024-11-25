/// Generate primes using a pseudo-uniform distribution.
///
/// The algorithm used in this module is described in two papers:
/// 1. ["Fast Generation of Prime Numbers on Portable Devices: An
///    Update"](https://marcjoye.github.io/papers/JP06pgen.pdf), aka JP06
/// 2. ["Close to Uniform Prime Number Generation With Fewer Random
///    Bits"](https://eprint.iacr.org/2011/481.pdf)
///
/// The second paper builds and improves on the first.
///
/// The motivation for the papers and algorithms is to provide a close to uniform prime generation
/// algorithm that is at the same time performant and consumes as little randomness as possible. The
/// yardstick here is searching for primes by random sampling (which is trivially uniform, given a
/// uniform CSPRNG). Compared to the main algorithm provided by this crate (known in literature as
/// "PRIMEINC"), this algorithm is much slower. Depending on the bitsize and parameter selection,
/// it's between 1.3 to 5x slower.
///
/// The main use case for this prime generator is when a large number of primes need to be generated
/// using the same sieve and/or on platforms with limited amounts of available entropy (e.g. IoT
/// devices, TEEs etc).
///
/// The available API is not complete; notably this generator samples primes in the whole bitspace
/// of the Uint and does not support searching for safe primes. PRs welcome!
///
/// Actors in this play are:
///     n: the number of bits in the prime we're looking for, e.g. 512
///     l: the number of top bits that are re-sampled on every iteration, e.g. 64
///     m: a product of all small odd primes up to a bound ß chosen such that the bit size of m is
///     `n - l`, e.g. ß = 512 - 64
///     b: picked uniformly at random among integers less than m and coprime to m (use unit
///     generation algorithm from JP06)
///     λ: Carmichael's function: the LCM of the λ(p)s for each prime used to compute m. Each λ(p)
///     is simply p-1 because each prime appears once and we exclude 2.
use crypto_bigint::modular::Retrieve;
use crypto_bigint::{
    Bounded, Constants, FixedInteger, Integer, Monty, NonZero, Odd, PowBoundedExp, RandomBits, RandomMod, U1024, U128,
    U2048, U256, U4096, U512, U64,
};
use rand_core::CryptoRngCore;
#[cfg(feature = "default-rng")]
use rand_core::OsRng;

use crate::is_prime_with_rng;

/// Prime search using a uniform sieve.
pub trait UniformSieve<T>
where
    T: Integer + Constants + Bounded + RandomBits + RandomMod + Copy,
{
    /// Returns a random prime using the provided CSPRNG.
    ///
    /// See [`is_prime_with_rng`][crate::presets::is_prime_with_rng] for details about the performed checks.
    fn generate_prime_with_rng(rng: &mut impl CryptoRngCore) -> T;

    /// Returns a random prime using [`OsRng`] as the CSPRNG.
    ///
    /// See [`is_prime_with_rng`][crate::presets::is_prime_with_rng] for details about the performed checks.
    #[cfg(feature = "default-rng")]
    fn generate_prime() -> T;

    // TODO(dp): missing
    // bit size of desired prime
    // generate_safe_prime
    // generate_safe_prime_with_rng
}

macro_rules! impl_generate_prime {
    ($(($name:ident, $m:expr, $lambda_m:expr, $a_max:expr)),+) => {
        $(
            impl UniformSieve<$name> for $name {
                fn generate_prime_with_rng(rng: &mut impl CryptoRngCore) -> Self {
                    const M: $name = $name::from_be_hex($m);
                    const LAMBDA_M: $name = $name::from_be_hex($lambda_m);
                    const A_MAX: $name = $name::from_be_hex($a_max);
                    let unit = jp06_unitgen(rng, M, LAMBDA_M);
                    let a_max = NonZero::new(A_MAX).expect("A_MAX is pre-calculated and known-good");
                    algorithm2(rng, unit, M, &a_max)
                }

                #[cfg(feature = "default-rng")]
                fn generate_prime() -> Self {
                    Self::generate_prime_with_rng(&mut OsRng)
                }
            }
        )+

    };
}

// Macro arguments: `uint type`, `m`, `λ(m)`, `a_max`
impl_generate_prime! {
    // Here `l` is 32, i.e. 64/2, i.e. using a 32 bit `m` and a 16 bit `λ(m)`
    (U64,   "00000017592B1B49", "000000000000D890", "000000000AF6E233")
    ,
    // Here `l` is 64, i.e. 128/2, i.e. using a 64 bit `m`, a 32 bit `λ(m)`, resulting in an `a_max` that fits in 64 bits
    (U128,  "00000000000000341DD47F9F45C500AF", "0000000000000000000000001CA73570", "000000000000000004E97D6751832168")
    ,
    // Here `l` is 128, i.e. 256/2, i.e. using a 128 bit `m` and a 64 bit `λ(m)`
    (U256,  "0000000000000011F29EDACF6E2E9A7C551B02D2F73F5BB75C429A6CC2E79B83",
            "000000000000000000000000000000000000000000000000004607AF4E8D0720",
            "0000000000000000000000000000000000000000000000000E437DB94E575A97")
    ,
    // Here `l` is 64, i.e. 512/8, i.e. using a 448 bit `m`, a 128 bit `λ(m)`, resulting in an `a_max` that fits in 64 bits
    (U512,  "00000000000000017C55C0F0A4177201114F0EE33F80274A3D96ACA2550BECC44E0A7D91FE964F9EB063A1988371CE58DBCE32B9F167C9555FA02A47B545338F",
            "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000272829291267EAF47536965B3B0AD1500",
            "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000AC4FAEF0A02ACA2D")
    ,
    // Here `l` is 64, i.e. 1024/16, i.e. using a 960 bit `m` and a ~256 bit `λ(m)`, resulting in an `a_max` that fits in 64 bits
    (U1024, "00000000000000016CC3AC9DC18F442E3F73D34147E253920667D800DA63CE9FEA167C22335097E1B207A8BAE729F4AB07F1BA5062481BC1166E2E5A42FF2393A9D7A7F5A195FB3FA7318A51407B138E9C8C54557AD65B9080FF8DB0F7672097CBC0DBC6EDA3CF451C20ACF1D003D909D87DA7BDA10D2C4772F7EEF762520D17",
            "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000291A6B42D1C7D2A7184D13E36F65773BBEFB4FA7996101300D49F09962A361F00",
            "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000B3AAAB7FFC52941B")
    ,
    // The performance of this algorithm sort of drops off a cliff for 2048 bit primes. I don't know why this is.
    // Trying different values for `l` (i.e. the number of bits that are re-randomized at each iteration), it seems like 128 bits (1/16) is the
    // most performant for U2048s. Oddly enough, doing the same for U1024 has next to no impact: using a 128 bit l is just as fast as a 64 bit l.
    // 128 bit l
    (U2048, "0000000000000000000000000000036C6E3E13331EC5D79DF6E1F12D4902E52EBD74678CAD34836377EA3F6A7D98535DCC65584C96430DE3D89EEA5FE88736A988E75874E2611F2AA3C99150555A3E1182E47F11BB9EB2D91C889DBDA6FEA36073EEBEB093E2F516994336516084203898201A73FD1637632DC4A98D066988261CC95023087D7E3D42EF0690252173F1A4146B660A6A019379DF5262CA76A1A3693F5DAC4443009AFA505F992AC852CE0D87FE456E94052755FDDBE7603D9DD469FDDFDE07938684E62FE2DF0F03A7F587C5F664895882508170D56D37E47F0F28B4475707DB8239C1507721CD10E3B3BE6073144C01738A4F17E7216D4847E5",
            "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000013682DE4C81E1A88DA9A9C674E3B34E1FEDD57D8A8F9F52AA39390FFE682B454EC328F9A1C455977E59489DE32175FEED2696B12E43100",
            "000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000004AC6A9A539EF07EDCED266DC673757")
    ,
    // Here `l` is 1024. Again, this is a rather surprising result, but the speedup is real and consistent compared to 64, 128, 256, 512 and 1024 bits
    (U4096, "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000029D24EFCB456540B7F56DF456A5AE69A265A164BCA43EE254756A73D5602520588A9C2E7B352D4FE5E15AB0103EA4E9CD3819FE76A97031810A118EBF9B15F349A8DE319B7863465A2BC14B24EF6E81D5C8F2540C0DC75C4E5C65CBB49CC01AEA791972EEBB02CB82BE2AF30AACDD7A5615FE5690448CC4E1E7DF448D1DBCC308E22EF51E01BAFD7F98B9175542931645CCE165681D0164466A24EEA419693B9CDADBE07BC3E11881ECE8FD6E60E96B1C5E026DFD241882CC74C66FEDC34509A0E8BB3184FFD4AB5C5650652FC1AE4D3F09F145F8A0CA17A6DC945C63F68A2B84EC79CB99AB75C88A6AEA2B60452C0DDDDA4FEA54AEEF0465DEC92DD581819365FC29378F05AFCCEC8D6378BA416D79C94B6E548C06801700AFAD22ED70F297F9F962B2A888798DEF444DAE08C7C6C1451AA2915B794471B71725005C7ABAE35563F05D5400497C64219C7F84BE819AE3EF3392227572A7C86569B057E1CCC0E5D175E0880C50D621EED19E68ED72228624E1D7F33A060C5B987DCF18ACAD38D9",
            "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001E16878440322FD1A0A8949F895ADA13E884B2EA9EB0DB3320A87B4A9A5653D89D36C428CDABF68B6F86A69BDBCF709FAD8207D78607C99FBC822BF7E8941CAE4F15ECB065A7623FD58C3162F00CDA0C5583E0152D00",
            "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000061F0A4B463AE4066954ADFEF304499A3A0CD686F8DD5B809EC68A0C774374E6356DBA3A66D71E6511E5AEDC4E16FD63695812261AE39EA1813B3E1FF274B4B5DA00D8439FEF0C49D5582537CDBC876A8701358FFD5C24B1ABB27ABF926A97EEE22707ABE31B2EA950DD23732723322D95657CF7192F37D2F851E1211334783E8")
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
#[inline(always)]
fn jp06_unitgen<T>(rng: &mut impl CryptoRngCore, m: T, lambda_m: T) -> T
where
    T: FixedInteger + RandomBits,
    T::Monty: Retrieve<Output = T> + Copy,
    <<T as Integer>::Monty as Monty>::Params: Copy,
{
    // 1. sample k in [1, m[
    let m_bits = m.bits_vartime();
    let mut k = T::random_bits(rng, m_bits);
    while k >= m {
        k = T::random_bits(rng, m_bits);
    }
    // 2. set `u = 1-(k^lambda_m) mod m`
    let prms = <T as Integer>::Monty::new_params_vartime(
        Odd::new(m).expect("m is a known, pre-calculated, non-zero, odd value"),
    );
    let one = T::Monty::one(prms);
    let zero = T::Monty::zero(prms);
    let m_monty_plus_one = T::Monty::new(m, prms) + one;

    let mut k = T::Monty::new(k, prms);
    let mut u = m_monty_plus_one - k.pow_bounded_exp(&lambda_m, lambda_m.bits_vartime());

    // if u != 0 {
    // 	pick random r from [1,m[
    // 	k = k + r*u mod m
    // 	go to step 2
    // }
    while u != zero {
        let mut r = T::random_bits(rng, m_bits);
        while r >= m {
            r = T::random_bits(rng, m_bits);
        }
        let r = T::Monty::new(r, prms);
        k += r * u;
        u = m_monty_plus_one - k.pow_bounded_exp(&lambda_m, lambda_m.bits_vartime());
    }
    k.retrieve()
}

// From https://eprint.iacr.org/2011/481.pdf, page 4:
// "Algorithm 2 More eﬃcient method"
// 1:   b $← {1,…, m−1}
// 2:   u ← (1−b^λ(m)) mod m
// 3:   if u != 0 then
// 4:       r $← {1,…, m−1}
// 5:       b ← b + ru mod m
// 6:       goto step 2
// 7:   end if
// 8:   repeat
// 9:       a $← {0,…, ⌊2^n/m⌋−1}
// 10:      p ← am + b
// 11:  until p is prime
// 12:  return p
// NOTE: Steps 1-7 are implemented in `jp06_unitgen`, where `b` is referred to as `k`.
#[allow(unused)]
#[inline(always)]
fn algorithm2<T>(rng: &mut impl CryptoRngCore, b: T, m: T, a_max: &NonZero<T>) -> T
where
    T: Integer + Bounded + RandomBits + RandomMod + Copy,
{
    let a_max = a_max.get();
    let a_max_bits = a_max.bits_vartime();
    let mut a = T::random_bits(rng, a_max_bits);
    while a >= a_max {
        a = T::random_bits(rng, a_max_bits);
    }

    let mut p = a * m + b;
    while !is_prime_with_rng(rng, &p) {
        a = T::random_bits(rng, a_max_bits);
        while a >= a_max {
            a = T::random_bits(rng, a_max_bits);
        }
        p = a * m + b;
    }
    p
}

// Calculates three constants, `m`, `λ(m)` and `a_max`, used in the uniform sieving algorithms.
// `m` is the product of all odd primes such that the product fits in a `T`.
// `λ(m)`, aka Carmichael's function is the LCM of the λ(p)s for each prime used to compute m. In
// our case, each λ(p) is simply p-1 because each prime appears once (and we exclude 2).
// `a_max` is the upper bound on the part of a prime candidate that is re-randomized on each iteration and set to 2^n/m -1
// Find constant values for common `T`s sized from 64 to 4096 in the test `generate_constants`.
#[cfg(test)]
fn calculate_constants<T>(a_max_bits: u32) -> (T, T, T)
where
    T: FixedInteger + crypto_bigint::Gcd<Output = T>,
{
    use crate::hazmat::precomputed::SMALL_PRIMES;
    use tracing::trace;

    let mut m = T::ONE;
    let mut lambda_m = T::ONE;
    let a_max_bits = core::cmp::max(32, a_max_bits);
    for (i, prime) in SMALL_PRIMES.iter().enumerate() {
        let prime_t = T::from(*prime);
        let prod_overflowed = m.checked_mul(&prime_t);
        let a_max = T::MAX.div(NonZero::new(m).unwrap()) - T::ONE;

        // Stop on one of two conditions: 1. the product of small primes does not fit in a `T`; 2. `a_max` would be too small to be useful
        if prod_overflowed.is_none().into() {
            trace!(
                "{} bits: Stopping after {} iterations, because prime={} overflows m. Last prime in m is {}",
                T::BITS,
                i - 1,
                prime,
                SMALL_PRIMES[i - 1]
            );
            break;
        } else if a_max.bits_vartime() <= a_max_bits {
            trace!(
                "{} bits: Stopping after {} iterations, because prime={} would make a_max smaller than {a_max_bits} bits. Last prime in m is {}",
                T::BITS, i-1, prime, SMALL_PRIMES[i - 1]
            );
            break;
        } else {
            m = prod_overflowed.unwrap();
            // Update λ(m) with λ(p) for this prime, which is just `p-1`
            let p_minus_one = prime_t - T::ONE;
            // LCM
            let gcd = lambda_m.gcd(&p_minus_one);
            let gcd_nz = NonZero::new(gcd).unwrap();
            lambda_m = lambda_m * prime_t.div(gcd_nz);
        }
    }
    trace!(
        "{} bits: dividing T::MAX with m={m:?} that has {} bits and {} leading zeros.",
        T::BITS,
        m.bits_vartime(),
        m.leading_zeros()
    );

    let a_max = T::MAX.div(NonZero::new(m).unwrap()) - T::ONE;
    trace!(
        "{} bits: a_max={a_max:?}, has {} bits and {} leading zeros.",
        T::BITS,
        a_max.bits_vartime(),
        a_max.leading_zeros()
    );
    (m, lambda_m, a_max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stat_utils::check_distribution_quality;

    use crypto_bigint::{U1024, U128, U2048, U256, U512, U64};
    use tracing::info;

    #[test_log::test]
    fn generate_constants() {
        let (m, lambda_m, a_max) = calculate_constants::<U64>(32);
        info!("64 bits: \nm={m:?}, \nlambda_m={lambda_m:?}, \na_max={a_max:?}");
        assert_eq!(m, U64::from_be_hex("00000017592B1B49"));
        assert_eq!(lambda_m, U64::from_be_hex("000000000000D890"));
        assert_eq!(a_max, U64::from_be_hex("000000000AF6E233"));

        let (m, lambda_m, a_max) = calculate_constants::<U128>(64);
        info!("128 bits: \nm={m:?}, \nlambda_m={lambda_m:?}, \na_max={a_max:?}");
        assert_eq!(m, U128::from_be_hex("00000000000000341DD47F9F45C500AF"));
        assert_eq!(lambda_m, U128::from_be_hex("0000000000000000000000001CA73570"));
        assert_eq!(a_max, U128::from_be_hex("000000000000000004E97D6751832168"));

        let (m, lambda_m, a_max) = calculate_constants::<U256>(64);
        info!("256 bits: \nm={m:?}, \nlambda_m={lambda_m:?}, \na_max={a_max:?}",);
        #[rustfmt::skip]
        assert_eq!(m, U256::from_be_hex("0000000000000011F29EDACF6E2E9A7C551B02D2F73F5BB75C429A6CC2E79B83"));
        #[rustfmt::skip]
        assert_eq!(lambda_m, U256::from_be_hex("000000000000000000000000000000000000000000000000004607AF4E8D0720"));
        #[rustfmt::skip]
        assert_eq!(a_max, U256::from_be_hex("0000000000000000000000000000000000000000000000000E437DB94E575A97"));

        let (m, lambda_m, a_max) = calculate_constants::<U512>(64);
        info!("512 bits: \nm={m:?}, \nlambda_m={lambda_m:?}, \na_max={a_max:?}");
        assert_eq!(m, U512::from_be_hex("00000000000000017C55C0F0A4177201114F0EE33F80274A3D96ACA2550BECC44E0A7D91FE964F9EB063A1988371CE58DBCE32B9F167C9555FA02A47B545338F"));
        assert_eq!(lambda_m, U512::from_be_hex("00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000272829291267EAF47536965B3B0AD1500"));
        assert_eq!(a_max, U512::from_be_hex("0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000AC4FAEF0A02ACA2D"));

        let (m, lambda_m, a_max) = calculate_constants::<U1024>(64);
        info!("1024 bits: \nm={m:?}, \nlambda_m={lambda_m:?}, \na_max={a_max:?}");
        assert_eq!(m, U1024::from_be_hex("00000000000000016CC3AC9DC18F442E3F73D34147E253920667D800DA63CE9FEA167C22335097E1B207A8BAE729F4AB07F1BA5062481BC1166E2E5A42FF2393A9D7A7F5A195FB3FA7318A51407B138E9C8C54557AD65B9080FF8DB0F7672097CBC0DBC6EDA3CF451C20ACF1D003D909D87DA7BDA10D2C4772F7EEF762520D17"));
        assert_eq!(lambda_m, U1024::from_be_hex("00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000291A6B42D1C7D2A7184D13E36F65773BBEFB4FA7996101300D49F09962A361F00"));
        assert_eq!(a_max, U1024::from_be_hex("000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000B3AAAB7FFC52941B"));

        let (m, lambda_m, a_max) = calculate_constants::<U2048>(128);
        info!("2048 bits: \nm={m:?}, \nlambda_m={lambda_m:?}, \na_max={a_max:?}");
        assert_eq!(m, U2048::from_be_hex("0000000000000000000000000000036C6E3E13331EC5D79DF6E1F12D4902E52EBD74678CAD34836377EA3F6A7D98535DCC65584C96430DE3D89EEA5FE88736A988E75874E2611F2AA3C99150555A3E1182E47F11BB9EB2D91C889DBDA6FEA36073EEBEB093E2F516994336516084203898201A73FD1637632DC4A98D066988261CC95023087D7E3D42EF0690252173F1A4146B660A6A019379DF5262CA76A1A3693F5DAC4443009AFA505F992AC852CE0D87FE456E94052755FDDBE7603D9DD469FDDFDE07938684E62FE2DF0F03A7F587C5F664895882508170D56D37E47F0F28B4475707DB8239C1507721CD10E3B3BE6073144C01738A4F17E7216D4847E5"));
        assert_eq!(lambda_m, U2048::from_be_hex("00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000013682DE4C81E1A88DA9A9C674E3B34E1FEDD57D8A8F9F52AA39390FFE682B454EC328F9A1C455977E59489DE32175FEED2696B12E43100"));
        assert_eq!(a_max, U2048::from_be_hex("000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000004AC6A9A539EF07EDCED266DC673757"));

        let (m, lambda_m, a_max) = calculate_constants::<U4096>(1024);
        info!("4096 bits: \nm={m:?}, \nlambda_m={lambda_m:?}, \na_max={a_max:?}");
        assert_eq!(m, U4096::from_be_hex("00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000029D24EFCB456540B7F56DF456A5AE69A265A164BCA43EE254756A73D5602520588A9C2E7B352D4FE5E15AB0103EA4E9CD3819FE76A97031810A118EBF9B15F349A8DE319B7863465A2BC14B24EF6E81D5C8F2540C0DC75C4E5C65CBB49CC01AEA791972EEBB02CB82BE2AF30AACDD7A5615FE5690448CC4E1E7DF448D1DBCC308E22EF51E01BAFD7F98B9175542931645CCE165681D0164466A24EEA419693B9CDADBE07BC3E11881ECE8FD6E60E96B1C5E026DFD241882CC74C66FEDC34509A0E8BB3184FFD4AB5C5650652FC1AE4D3F09F145F8A0CA17A6DC945C63F68A2B84EC79CB99AB75C88A6AEA2B60452C0DDDDA4FEA54AEEF0465DEC92DD581819365FC29378F05AFCCEC8D6378BA416D79C94B6E548C06801700AFAD22ED70F297F9F962B2A888798DEF444DAE08C7C6C1451AA2915B794471B71725005C7ABAE35563F05D5400497C64219C7F84BE819AE3EF3392227572A7C86569B057E1CCC0E5D175E0880C50D621EED19E68ED72228624E1D7F33A060C5B987DCF18ACAD38D9"));
        assert_eq!(lambda_m, U4096::from_be_hex("0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001E16878440322FD1A0A8949F895ADA13E884B2EA9EB0DB3320A87B4A9A5653D89D36C428CDABF68B6F86A69BDBCF709FAD8207D78607C99FBC822BF7E8941CAE4F15ECB065A7623FD58C3162F00CDA0C5583E0152D00"));
        assert_eq!(a_max, U4096::from_be_hex("00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000061F0A4B463AE4066954ADFEF304499A3A0CD686F8DD5B809EC68A0C774374E6356DBA3A66D71E6511E5AEDC4E16FD63695812261AE39EA1813B3E1FF274B4B5DA00D8439FEF0C49D5582537CDBC876A8701358FFD5C24B1ABB27ABF926A97EEE22707ABE31B2EA950DD23732723322D95657CF7192F37D2F851E1211334783E8"));
    }

    #[test_log::test]
    fn uniform_sieve_primes() {
        let p: U64 = U64::generate_prime();
        info!("64 bit prime={p:?}");
        assert!(is_prime(&p));

        let p = U128::generate_prime();
        info!("128 bit prime={p:?}");
        assert!(is_prime(&p));

        let p = U256::generate_prime();
        info!("256 bit prime={p:?}");
        assert!(is_prime(&p));

        let p = U512::generate_prime();
        info!("512 bit prime={p:?}");
        assert!(is_prime(&p));

        let p = U1024::generate_prime();
        info!("1024 bit prime={p:?}");
        assert!(is_prime(&p));

        let p = U2048::generate_prime();
        info!("2048 bit prime={p:?}");
        assert!(is_prime(&p));
        // // This is very slow
        // let p = U4096::generate_prime();
        // info!("4096 bit prime={p:?}");
        // assert!(is_prime(&p));
    }

    #[test_log::test]
    fn check_distribution_quality_u64() {
        check_distribution_quality::<U64>();
    }
    #[test_log::test]
    fn check_distribution_quality_u128() {
        check_distribution_quality::<U128>();
    }
    #[test_log::test]
    fn check_distribution_quality_u256() {
        check_distribution_quality::<U256>();
    }
    #[test_log::test]
    fn check_distribution_quality_u512() {
        check_distribution_quality::<U512>();
    }
    #[cfg(feature = "tests-exhaustive")]
    #[test_log::test]
    fn check_distribution_quality_u1024() {
        check_distribution_quality::<U1024>();
    }
    #[cfg(feature = "tests-exhaustive")]
    #[test_log::test]
    fn check_distribution_quality_u2048() {
        check_distribution_quality::<U2048>();
    }
}
