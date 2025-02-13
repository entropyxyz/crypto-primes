//! Const-context floating point functions that are currently not present in `core`.

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
    const LN_2: f64 = core::f64::consts::LN_2;
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

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;
    use proptest::prelude::*;

    use super::{floor, floor_sqrt, pow, two_powf_normalized_lower_bound};

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
}
