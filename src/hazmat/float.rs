//! Const-context floating point functions that are currently not present in `core`.

use core::f64;
use crypto_bigint::Uint;
#[allow(unused_imports)]
use num_traits::float::FloatCore as _;

/// Calculates `base^exp`.
const fn pow(mut base: f64, mut exp: u32) -> f64 {
    let mut result = 1.;
    while exp > 0 {
        if exp & 1 == 1 {
            result *= base;
        }
        base *= base;
        exp >>= 1;
    }
    result
}

/// Calculates `2^exp`.
pub(crate) const fn two_powi(exp: u32) -> f64 {
    pow(2f64, exp.abs_diff(0))
}

/// Calculates `floor(x)`.
// Taken from `libm` crate.
const fn floor(x: f64) -> f64 {
    const TOINT: f64 = 1. / f64::EPSILON;

    let ui = x.to_bits();
    let e = ((ui >> 52) & 0x7ff) as i32;

    if (e >= 0x3ff + 52) || (x == 0.) {
        return x;
    }
    /* y = int(x) - x, where int(x) is an integer neighbor of x */
    let y = if (ui >> 63) != 0 {
        x - TOINT + TOINT - x
    } else {
        x + TOINT - TOINT - x
    };
    /* special case because of non-nearest rounding modes */
    if e < 0x3ff {
        return if (ui >> 63) != 0 { -1. } else { 0. };
    }
    if y > 0. {
        x + y - 1.
    } else {
        x + y
    }
}

/// Calculates a lower bound approximation of `2^exp` where `0 <= exp <= 1`.
const fn two_powf_normalized_lower_bound(exp: f64) -> f64 {
    debug_assert!(exp >= 0.);
    debug_assert!(exp <= 1.);

    // Use the first four terms of the Taylor expansion to calculate `2^exp`.
    // The error is under 2%.
    //
    // Since it is a monotonous function, `res <= 2^exp`.

    let exp_2 = exp * exp;
    let exp_3 = exp_2 * exp;

    // The coefficients are `ln(2)^n / n!`, where `n` is the power of the corresponding term.
    const LN_2: f64 = f64::consts::LN_2;
    const C1: f64 = LN_2;
    const C2: f64 = LN_2 * LN_2 / 2.;
    const C3: f64 = LN_2 * LN_2 * LN_2 / 6.;

    1. + C1 * exp + C2 * exp_2 + C3 * exp_3
}

/// Calculates an approximation of `2^exp` where `exp < 0`.
/// The approximation is guaranteed to always be greater than `2^exp`.
pub(crate) const fn two_powf_upper_bound(exp: f64) -> f64 {
    debug_assert!(exp < 0.);

    let positive_exp = -exp;

    let int_part = floor(positive_exp);
    let frac_part = positive_exp - int_part;

    let int_res = two_powi(int_part as u32);
    let frac_res = two_powf_normalized_lower_bound(frac_part);

    // `int_res * frac_res <= 2^(int_part + frac_part)`,
    // so when we invert it, we get the upper bound approximation instead.
    1. / (int_res * frac_res)
}

/// Calculates `floor(sqrt(x))`.
pub(crate) const fn floor_sqrt(x: u32) -> u32 {
    if x < 2 {
        return x;
    }

    // Initialize the binary search bounds.
    let mut low = 1;
    let mut high = x / 2;

    while low <= high {
        let mid = (low + high) / 2;
        let mid_squared = mid * mid;

        // Check if `mid` is the floor of the square root of `x`.
        if mid_squared <= x && (mid + 1) * (mid + 1) > x {
            break;
        } else if mid_squared < x {
            low = mid + 1;
        } else {
            high = mid - 1
        }
    }

    (low + high) / 2
}

// Copied from the `libm` crate, which is licensed under the MIT license.

/* origin: FreeBSD /usr/src/lib/msun/src/e_log2.c */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */
/*
 * Return the base 2 logarithm of x.  See log.c for most comments.
 *
 * Reduce x to 2^k (1+f) and calculate r = log(1+f) - f + f*f/2
 * as in log.c, then combine and scale in extra precision:
 *    log2(x) = (f - f*f/2 + r)/log(2) + k
 */

/// The base 2 logarithm of `x` (f64).
#[allow(clippy::all)]
pub(crate) const fn log2(mut x: f64) -> f64 {
    const IVLN2HI: f64 = 1.44269504072144627571e+00; /* 0x3ff71547, 0x65200000 */
    const IVLN2LO: f64 = 1.67517131648865118353e-10; /* 0x3de705fc, 0x2eefa200 */
    const LG1: f64 = 6.666666666666735130e-01; /* 3FE55555 55555593 */
    const LG2: f64 = 3.999999999940941908e-01; /* 3FD99999 9997FA04 */
    const LG3: f64 = 2.857142874366239149e-01; /* 3FD24924 94229359 */
    const LG4: f64 = 2.222219843214978396e-01; /* 3FCC71C5 1D8E78AF */
    const LG5: f64 = 1.818357216161805012e-01; /* 3FC74664 96CB03DE */
    const LG6: f64 = 1.531383769920937332e-01; /* 3FC39A09 D078C69F */
    const LG7: f64 = 1.479819860511658591e-01; /* 3FC2F112 DF3E5244 */

    let x1p54 = f64::from_bits(0x4350000000000000); // 0x1p54 === 2 ^ 54

    let mut ui: u64 = x.to_bits();
    let hfsq: f64;
    let f: f64;
    let s: f64;
    let z: f64;
    let r: f64;
    let mut w: f64;
    let t1: f64;
    let t2: f64;
    let y: f64;
    let mut hi: f64;
    let lo: f64;
    let mut val_hi: f64;
    let mut val_lo: f64;
    let mut hx: u32;
    let mut k: i32;

    hx = (ui >> 32) as u32;
    k = 0;
    if hx < 0x00100000 || (hx >> 31) > 0 {
        if ui << 1 == 0 {
            return -1. / (x * x); /* log(+-0)=-inf */
        }
        if (hx >> 31) > 0 {
            return (x - x) / 0.0; /* log(-#) = NaN */
        }
        /* subnormal number, scale x up */
        k -= 54;
        x *= x1p54;
        ui = x.to_bits();
        hx = (ui >> 32) as u32;
    } else if hx >= 0x7ff00000 {
        return x;
    } else if hx == 0x3ff00000 && ui << 32 == 0 {
        return 0.;
    }

    /* reduce x into [sqrt(2)/2, sqrt(2)] */
    hx += 0x3ff00000 - 0x3fe6a09e;
    k += (hx >> 20) as i32 - 0x3ff;
    hx = (hx & 0x000fffff) + 0x3fe6a09e;
    ui = ((hx as u64) << 32) | (ui & 0xffffffff);
    x = f64::from_bits(ui);

    f = x - 1.0;
    hfsq = 0.5 * f * f;
    s = f / (2.0 + f);
    z = s * s;
    w = z * z;
    t1 = w * (LG2 + w * (LG4 + w * LG6));
    t2 = z * (LG1 + w * (LG3 + w * (LG5 + w * LG7)));
    r = t2 + t1;

    /* hi+lo = f - hfsq + s*(hfsq+R) ~ log(1+f) */
    hi = f - hfsq;
    ui = hi.to_bits();
    ui &= (-1i64 as u64) << 32;
    hi = f64::from_bits(ui);
    lo = f - hi - hfsq + s * (hfsq + r);

    val_hi = hi * IVLN2HI;
    val_lo = (lo + hi) * IVLN2LO + lo * IVLN2HI;

    /* spadd(val_hi, val_lo, y), except for not using double_t: */
    y = k as f64;
    w = y + val_hi;
    val_lo += (y - w) + val_hi;
    val_hi = w;

    val_lo + val_hi
}

// Calculate the natural logarithm of a big integer using the relation ln(x) = log₂(x) / log₂(e).
// Uses fixed-point arithmetic for large values of x (> 2^53).
#[inline]
pub(crate) fn ln<const LIMBS: usize>(x: &Uint<LIMBS>) -> f64 {
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
    use alloc::vec;
    use alloc::vec::Vec;

    use crypto_bigint::{Random, U1024, U128, U256, U512};
    use float_cmp::assert_approx_eq;
    use proptest::prelude::*;
    use rand_core::SeedableRng;

    use super::{floor, floor_sqrt, ln, log2, pow, two_powf_normalized_lower_bound};

    #[test]
    fn sqrt_corner_cases() {
        assert_eq!(floor_sqrt(0), 0);
        assert_eq!(floor_sqrt(1), 1);
        assert_eq!(floor_sqrt(2), 1);
    }

    proptest! {
        #[test]
        fn fuzzy_pow(base in 0..100u32, exp in 0..30u32) {
            let base_f = base as f64 / 100.;
            let test = pow(base_f, exp);
            let reference = base_f.powf(exp as f64);
            assert_approx_eq!(f64, test, reference, ulps = 20);
        }

        #[test]
        fn fuzzy_floor(x in proptest::num::f64::NORMAL) {
            let test = floor(x);
            let reference = x.floor();
            assert_approx_eq!(f64, test, reference);
        }

        #[test]
        fn fuzzy_two_powf_upper_bound(exp in 0..1000) {
            let exp_f = exp as f64 / 1000.;
            let test = two_powf_normalized_lower_bound(exp_f);
            let reference = 2f64.powf(exp_f);
            assert!(test <= reference);
            assert!((reference - test) / reference <= 0.02);
        }

        #[test]
        fn fuzzy_floor_sqrt(x in 0..100000u32) {
            let x_f = x as f64;
            let test = floor_sqrt(x);
            let reference = x_f.sqrt().floor() as u32;
            assert_eq!(test, reference);
        }
    }

    #[test]
    fn test_log2_special_cases() {
        // Test lines 56-57: log(±0) = -inf
        let log_zero = log2(0.0);
        assert!(log_zero.is_infinite() && log_zero.is_sign_negative(), "log2(0.0) is −∞");
        let log_neg_zero = log2(-0.0);
        assert!(
            log_neg_zero.is_infinite() && log_neg_zero.is_sign_negative(),
            "log2(-0.0) is −∞"
        );

        // Test lines 59-60: log(-#) = NaN
        let log_neg_one = log2(-1.0);
        assert!(log_neg_one.is_nan(), "log2(-1.0) is NaN");
        let log_neg_inf = log2(f64::NEG_INFINITY);
        assert!(log_neg_inf.is_nan(), "log2(-Inf) is NaN");

        // Test lines 67-68: log(Inf) = Inf, log(NaN) = NaN
        let log_inf = log2(f64::INFINITY);
        assert!(log_inf.is_infinite() && log_inf.is_sign_positive(), "log2(Inf) is ∞");
        assert!(log2(f64::NAN).is_nan(), "log2(NaN) is NaN");

        // Test lines 69-70: log(1.0) = 0.0
        assert_eq!(log2(1.0), 0.0, "log2(1.0) is zero");

        // Test lines 55 & 62-66: Subnormal input
        let log_subnormal = log2(f64::MIN_POSITIVE); // Smallest positive f64 > 0
        assert!(log_subnormal.is_finite(), "log2(subnormal) is finite");
        // log2(MIN_POSITIVE) is log2(2^-1074) = -1074, but the log2 code approximates it as -1022 because floating
        // point maths is weird.
        assert!(log_subnormal < -1000.0, "log2(subnormal) is smaller than 1000");
    }

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
}
