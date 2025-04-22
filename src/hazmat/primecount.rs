use core::f64;
use crypto_bigint::{Concat, NonZero, Split, Uint};
#[allow(unused_imports)]
use num_traits::float::FloatCore as _;

use super::log2::log2;

/// Estimate the number of primes smaller than x using the asymptotic expansion of Li(x) with 4 terms, i.e.:
///
///   `𝜋(𝑥) ~ x/ln x * (1 + 1!/ln x + 2!/ln^2 x + 3!/ln^3 x).`
///
/// # Error considerations
///
/// There are two components to the error:
///
/// 1. The truncation error from using a limited number of terms (4) of the asymptotic expansion of Li(x): `|Li(x) -
///    Li_approx(x)|`
/// 2. The theoretical error, given by the absolute difference between π(x) and Li(x): `|π(x) - Li(x)|`
///
/// ## Truncation Error
///
/// The size of the truncation error for the 4-term approximation `Li_approx(x)` is roughly the size of the first
/// omitted term, which is `24x / (ln x)^5`.
///
/// | Bits | Truncation error (abs)     | Truncation error (rel) |
/// | :--- | :------------------------- | :----------------------|
/// | 1024 | ~2^981                     | ~2^-33                 |
/// | 2048 | ~2^2000                    | ~2^-37                 |
/// | 4096 | ~2^4043                    | ~2^-41                 |
///
/// *(Relative error is calculated as the first omitted term divided by / Li_approx(x)`)*
///
/// ## Theoretical Error
///
/// Assuming RH, we can use Schoenfeld’s bound `|π(x) - Li(x)| < (1 / 8π) * √x * ln x`. This gives a very small
/// theoretical error bound compared to the best known unconditional bounds (by Trudgian).
///
/// | Bits | Error Bound (abs) | Error Bound (rel) |
/// | :--- | :---------------- | :-----------------|
/// | 1024 | < 2^517           | < 2^-497          |
/// | 2048 | < 2^1030          | < 2^-1007         |
/// | 4096 | < 2^2055          | < 2^-2029         |
///
/// *(Relative error bound is calculated as `Schoenfeld’s Bound (abs) / Li_approx(x)`)*
///
/// ## Discussion
/// Assuming RH, the dominant error term in estimating `π(x)` with `Li_approx(x)` comes from truncating the asymptotic
/// expansion of Li(x). While this truncation error is extremely small in relative terms, it is large in absolute terms.
/// Improving the Li(x) approximation to be the same order of magnitude as the Schoenfeld bound would require using
/// hundreds of terms from the asymptotic series (or employing more complex methods like continued fractions for the
/// remainder, which is beyond the scope of this simple approximation). In relative terms the relative error bound is
/// approximately ~2^-33 (for 1024 bits) down to ~2^-41 (for 4096 bits). This corresponds to an extremely small
/// percentage error (significantly less than 10^-8 %).
///
/// It should be noted that while Li(x) is generally smaller than π(x) for 'small' x, it is known that the sign of
/// `π(x) - Li(x)` changes infinitely often. It has been proven that there must be a crossing below ~10^316 (~2^1051),
/// which is well within the ranges used in this library. Thus, users should be aware that the estimate provided here
/// can be both greater than and smaller than the actual value of π(x).
pub fn estimate_pi_x<const LIMBS: usize, const RHS_LIMBS: usize>(x: &Uint<LIMBS>) -> Uint<LIMBS>
where
    Uint<LIMBS>: Concat<Output = Uint<RHS_LIMBS>>,
    Uint<RHS_LIMBS>: Split<Output = Uint<LIMBS>>,
{
    // Number of bits to scale by (effectively fractional bits during division)
    const SCALE_BITS: u32 = 64;

    let ln_x = ln(x);

    if !ln_x.is_finite() || ln_x <= 1.0 {
        return Uint::ZERO;
    }

    // Calculate f64 denominators L, L^2, L^3, L^4 and round to u64
    let ln_x_2 = ln_x * ln_x;
    let ln_x_3 = ln_x_2 * ln_x;
    let ln_x_4 = ln_x_3 * ln_x;

    // Round to nearest u64. Ensure results are non-zero for division.
    let d1 = ln_x.round().max(1.0) as u64;
    let d2 = ln_x_2.round().max(1.0) as u64;
    let d3 = ln_x_3.round().max(1.0) as u64;
    let d4 = ln_x_4.round().max(1.0) as u64;

    let d1 = NonZero::new(Uint::<RHS_LIMBS>::from(d1)).expect("at least 1 by construction");
    let d2 = NonZero::new(Uint::<RHS_LIMBS>::from(d2)).expect("at least 1 by construction");
    let d3 = NonZero::new(Uint::<RHS_LIMBS>::from(d3)).expect("at least 1 by construction");
    let d4 = NonZero::new(Uint::<RHS_LIMBS>::from(d4)).expect("at least 1 by construction");

    let x_wide: Uint<RHS_LIMBS> = x.concat(&Uint::ZERO);
    let term1_scaled = (x_wide << SCALE_BITS).wrapping_div(&d1);
    let term2_scaled = (x_wide << SCALE_BITS).wrapping_div(&d2);

    // Term 3: (2x << SCALE_BITS) / d3  (Factorial 2!)
    let x_times_2: Uint<RHS_LIMBS> = x_wide << 1;
    let term3_scaled = (x_times_2 << SCALE_BITS).wrapping_div(&d3);

    // Term 4: (6x << SCALE_BITS) / d4  (Factorial 3!)
    let six = Uint::<LIMBS>::from(6u64);
    let x_times_6 = x_wide.saturating_mul(&six);
    let term4_scaled = (x_times_6 << SCALE_BITS).wrapping_div(&d4);

    let sum_scaled = term1_scaled
        .wrapping_add(&term2_scaled)
        .wrapping_add(&term3_scaled)
        .wrapping_add(&term4_scaled);

    // Descale by right-shifting
    let li_x = sum_scaled >> SCALE_BITS;
    let (lo, hi) = li_x.split();
    assert_eq!(hi, Uint::ZERO, "De-scaling should leave the high half zero");
    lo
}

// Calculate the natural logarithm of a big integer using the relation ln(x) = log₂(x) / log₂(e).
// Uses fixed-point arithmetic for large values of x (> 2^53).
#[inline]
fn ln<const LIMBS: usize>(x: &Uint<LIMBS>) -> f64 {
    if x <= &Uint::ONE {
        return 0.0;
    }
    let ilog2_x = x.bits_vartime().saturating_sub(1);
    // if x is small enough to be cast losslessly to an f64 we use the normal 64-bit log2() from `libm`.
    if ilog2_x < f64::MANTISSA_DIGITS {
        let x = {
            #[cfg(target_pointer_width = "64")]
            {
                x.as_limbs()[0].0
            }
            #[cfg(target_pointer_width = "32")]
            {
                // Small enough to fit in 32 bits
                if ilog2_x < 32 {
                    x.as_limbs()[0].0 as u64
                } else {
                    let lo = x.as_limbs()[0].0;
                    let hi = x.as_limbs()[1].0;
                    (hi as u64) << 32 | lo as u64
                }
            }
        };
        return log2(x as f64) * f64::consts::LN_2;
    }
    // x can be approximated by M*2^shift, where shift is `ilog2(x) - 52` and M is the integer value represented by the
    // top 53 bits of x.
    // log2(x) ~ log2(M*2^shift) ~ log2(M) + shift ~ log2(M) + ilog2(x) - 52
    let shift = ilog2_x.saturating_sub(f64::MANTISSA_DIGITS - 1);
    let shifted_x = x.wrapping_shr_vartime(shift);
    #[cfg(target_pointer_width = "64")]
    let fraction = shifted_x.as_limbs()[0].0 as f64;
    #[cfg(target_pointer_width = "32")]
    let fraction = {
        let lo = shifted_x.as_limbs()[0].0;
        let hi = shifted_x.as_limbs()[1].0;
        let fraction = (hi as u64) << 32 | lo as u64;
        fraction as f64
    };

    // Fraction is now m * 2^52, where m is the top 53 bits of x. Take log2(m) and subtract 52 to scale the result back
    // to the expected range.
    let fraction = log2(fraction) - (f64::MANTISSA_DIGITS - 1) as f64;
    let log2_x = ilog2_x as f64 + fraction;
    log2_x * f64::consts::LN_2
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use alloc::vec::Vec;
    use crypto_bigint::{Random, U1024, U128, U256, U512};
    use rand_core::SeedableRng;

    #[test]
    fn test_ln_x_small() {
        let test_cases: Vec<(u128, f64)> = vec![
            (1, 0.0),
            (3, 1.09861228866811),
            (4, 1.3862943611198906),
            (5, 1.6094379124341004),
            (10u128.pow(2), 4.60517018598809),
            (10u128.pow(3), 6.90775527898214),
            (10u128.pow(4), 9.21034037197618),
            (10u128.pow(5), 11.51292546497023),
            (10u128.pow(6), 13.81551055796427),
            (10u128.pow(7), 16.11809565095832),
            (10u128.pow(8), 18.42068074395237),
            (10u128.pow(9), 20.72326583694641),
            (10u128.pow(10), 23.02585092994046),
            (10u128.pow(11), 25.3284360229345),
            (10u128.pow(12), 27.63102111592855),
            (10u128.pow(13), 29.93360620892259),
            (10u128.pow(14), 32.23619130191664),
            (10u128.pow(15), 34.53877639491069),
            (10u128.pow(16), 36.84136148790473),
            (10u128.pow(17), 39.14394658089878),
            (10u128.pow(18), 41.44653167389282),
            (10u128.pow(19), 43.74911676688687),
            (10u128.pow(20), 46.05170185988091),
            (10u128.pow(21), 48.35428695287496),
            (10u128.pow(22), 50.65687204586901),
            (10u128.pow(23), 52.95945713886305),
            (10u128.pow(24), 55.262042231857096),
            (10u128.pow(25), 57.56462732485114),
            (10u128.pow(26), 59.86721241784519),
            (10u128.pow(27), 62.16979751083923),
            (10u128.pow(28), 64.47238260383328),
            (10u128.pow(29), 66.77496769682733),
            (10u128.pow(30), 69.07755278982137),
            (10u128.pow(31), 71.38013788281542),
            (10u128.pow(32), 73.68272297580946),
            (10u128.pow(33), 75.9853080688035),
            (10u128.pow(34), 78.28789316179756),
            (10u128.pow(35), 80.5904782547916),
            (10u128.pow(36), 82.89306334778564),
            (10u128.pow(37), 85.19564844077969),
            (10u128.pow(38), 87.49823353377374),
        ];
        for (x, expected) in test_cases.iter() {
            let x_big = U128::from_u128(*x);
            let result = ln(&x_big);
            assert!(
                (result - *expected).abs() < f64::EPSILON * 100.0,
                "x: {x}, mine: {result}, expected: {expected}"
            );
        }
    }

    #[test]
    fn test_ln_x_large() {
        let mut rng = rand_chacha::ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");
        let expected_256 = 172.29242258375361;

        let x = U256::random(&mut rng);
        let ln_x = ln(&x);
        assert!(
            (ln_x - expected_256).abs() < f64::EPSILON,
            "x: {x}, ln x: {ln_x}, expected: {expected_256}"
        );

        let expected_512 = 354.2665608707877;
        let mut rng = rand_chacha::ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");
        let x = U512::random(&mut rng);
        let ln_x = ln(&x);
        assert!(
            (ln_x - expected_512).abs() < f64::EPSILON,
            "x: {x}, ln x: {ln_x}, expected: {expected_512}"
        );

        let expected_1024 = 708.828942204781;
        let mut rng = rand_chacha::ChaCha8Rng::from_seed(*b"01234567890123456789012345678901");
        let x = U1024::random(&mut rng);
        let ln_x = ln(&x);
        assert!(
            (ln_x - expected_1024).abs() < f64::EPSILON,
            "x: {x}, ln x: {ln_x}, expected: {expected_1024}"
        );
    }

    #[test]
    fn pi_x_2_500() {
        let x = Uint::ONE << 500;
        let sage_est = U1024::from_be_hex("0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000bda54dd744907290defac4f74bec507fdafd96e123c49bea56826f73702a469b67453a13a6abc40e81b760a0a5fd95870dbb8bbe99973c246c49561e101");
        let estimate = estimate_pi_x(&x);
        let delta = if sage_est > estimate {
            sage_est - estimate
        } else {
            estimate - sage_est
        };
        const MIN_BIT_DIFFERENCE: u32 = 10;
        assert!(
            sage_est.bits_vartime() - delta.bits_vartime() >= MIN_BIT_DIFFERENCE,
            "Estimate not close enough: delta has {} bits, expected {} bits. Difference should be >= {}",
            delta.bits_vartime(),
            sage_est.bits_vartime(),
            MIN_BIT_DIFFERENCE
        );
    }

    #[test]
    fn pi_x_max() {
        let x = U1024::MAX;
        let estimate = estimate_pi_x(&x);
        let sage_est = Uint::from_be_hex("005c7682fe13533e630c22e716b35439b3dc61f1d4898d78a36dd9c9afc0745a06d3a0deb93b77423f6d11c107283fcfdb8ae17de22b5197972f37cb480a2737fe8d0f15202bb43bc1863b05f6d3849f865b95242eaec9789dcf3b40e92504d98258f80b394ebec1c63d1186f9552689076f709c2fd8497b5f78d82cea2c2137");
        let delta = if sage_est > estimate {
            sage_est - estimate
        } else {
            estimate - sage_est
        };
        const MIN_BIT_DIFFERENCE: u32 = 12;
        assert!(
            sage_est.bits_vartime() - delta.bits_vartime() >= MIN_BIT_DIFFERENCE,
            "Estimate not close enough: delta has {} bits, expected {} bits. Difference should be >= {}",
            delta.bits_vartime(),
            sage_est.bits_vartime(),
            MIN_BIT_DIFFERENCE
        );
    }
    #[test]
    fn pi_x_random() {
        // This is a ChaCha8Rng with seed b"01234567890123456789012345678901"
        let x = Uint::from_be_hex("62A211E0907141403FD3EB60A82EAB701524710BDB024EB68DFF309389258B632EB9975D29F028F5137AC9DE870EB622D2D45A0D3A9C5801E8A3109BED220F82890E108F1778E5523E3E89CCD5DEDB667E6C17E940E9D4C3F58575C86CB76403017AD59D33AC084D2E58D81F8BB87A61B44677037A7DBDE04814256570DCBD7A");
        let sage_est = U1024::from_be_hex("0023ac3184a0c4c8e9025e0ae9b44d7980cee1baacf69032bb898677841fac0e516fa6bc8c1d1d3bb282622aa62c49f2d8e622d2f9aa80af3140c8c2251363017c99621943c90ab55a6dd69a678110233254a1a3c50ceb1cdb516e7220a7514a17b20114c7bef6f316e94cf7c9181187d70e751bda2e18695fa71e8015b8cf1c");
        let estimate = estimate_pi_x(&x);
        let delta = if sage_est > estimate {
            sage_est - estimate
        } else {
            estimate - sage_est
        };
        const MIN_BIT_DIFFERENCE: u32 = 12;
        assert!(
            sage_est.bits_vartime() - delta.bits_vartime() >= MIN_BIT_DIFFERENCE,
            "Estimate not close enough: delta has {} bits, expected {} bits. Difference should be >= {}",
            delta.bits_vartime(),
            sage_est.bits_vartime(),
            MIN_BIT_DIFFERENCE
        );
    }

    #[test]
    fn pi_estimates_for_known_values() {
        // List of known values for π(x), expressed as a tuple of `(π(x), exponent)`, where the `exponent` is used with base 10.
        let pi_xs: Vec<(u128, u32)> = vec![
            // The error is large for small x, so we skip them.
            // (4, 1),
            // (25, 2),
            // (168, 3),
            // (1229, 4),
            // (9592, 5),
            (78498, 6),
            (664579, 7),
            (5761455, 8),
            (50847534, 9),
            (455052511, 10),
            (4118054813, 11),
            (37607912018, 12),
            (346065536839, 13),
            (3204941750802, 14),
            (29844570422669, 15),
            (279238341033925, 16),
            (2623557157654233, 17),
            (24739954287740860, 18),
            (234057667276344607, 19),
            (2220819602560918840, 20),
            (21127269486018731928, 21),
            (201467286689315906290, 22),
            (1925320391606803968923, 23),
            (18435599767349200867866, 24),
            (176846309399143769411680, 25),
            (1699246750872437141327603, 26),
            (16352460426841680446427399, 27),
            (157589269275973410412739598, 28),
            (1520698109714272166094258063, 29),
        ];
        for (pi_x, exponent) in pi_xs.iter() {
            let pi_x_wide = U256::from_u128(*pi_x);
            let n = U256::from_u128(10u128.pow(*exponent));
            let estimate = estimate_pi_x(&n);
            let delta = if pi_x_wide > estimate {
                pi_x_wide - estimate
            } else {
                estimate - pi_x_wide
            };
            let delta = uint_to_u128(&delta);
            let estimate_128 = uint_to_u128(&estimate);
            let error = (delta as f64 / *pi_x as f64) * 100.0;
            assert!(
                error < 2.2,
                "10^{exponent}:\t{pi_x} - {estimate_128} = {delta}, err: {error:.2}"
            );
        }
    }

    fn uint_to_u128<const LIMBS: usize>(x: &Uint<LIMBS>) -> u128 {
        let limbs = x.as_limbs();
        #[cfg(target_pointer_width = "32")]
        return (limbs[3].0 as u128) << 96
            | (limbs[2].0 as u128) << 64
            | (limbs[1].0 as u128) << 32
            | limbs[0].0 as u128;
        #[cfg(target_pointer_width = "64")]
        return ((limbs[1].0 as u128) << 64) | limbs[0].0 as u128;
    }
}
