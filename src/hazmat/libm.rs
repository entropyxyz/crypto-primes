//! Functions extracted from `libm` to avoid a dependency on it
//! (which sometimes breaks builds due to linking to C libraries).
//!
//! Some reformatting and clippy exceptions were added.

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

#[expect(clippy::excessive_precision)]
const IVLN2HI: f64 = 1.44269504072144627571e+00; /* 0x3ff71547, 0x65200000 */

#[expect(clippy::excessive_precision)]
const IVLN2LO: f64 = 1.67517131648865118353e-10; /* 0x3de705fc, 0x2eefa200 */

#[expect(clippy::excessive_precision)]
const LG1: f64 = 6.666666666666735130e-01; /* 3FE55555 55555593 */

#[expect(clippy::excessive_precision)]
const LG2: f64 = 3.999999999940941908e-01; /* 3FD99999 9997FA04 */

#[expect(clippy::excessive_precision)]
const LG3: f64 = 2.857142874366239149e-01; /* 3FD24924 94229359 */

#[expect(clippy::excessive_precision)]
const LG4: f64 = 2.222219843214978396e-01; /* 3FCC71C5 1D8E78AF */

#[expect(clippy::excessive_precision)]
const LG5: f64 = 1.818357216161805012e-01; /* 3FC74664 96CB03DE */

#[expect(clippy::excessive_precision)]
const LG6: f64 = 1.531383769920937332e-01; /* 3FC39A09 D078C69F */

#[expect(clippy::excessive_precision)]
const LG7: f64 = 1.479819860511658591e-01; /* 3FC2F112 DF3E5244 */

/// The base 2 logarithm of `x` (f64).
pub(crate) fn log2(mut x: f64) -> f64 {
    let x1p54 = f64::from_bits(0x4350000000000000); // 0x1p54 === 2 ^ 54

    let mut ui = x.to_bits();
    let mut k = 0i32;

    #[expect(clippy::as_conversions)] // pre-shifted, so no overflow can occur
    let mut hx = (ui >> 32) as u32;

    if hx < 0x00100000 || (hx >> 31) > 0 {
        if ui << 1 == 0 {
            return -1. / (x * x); /* log(+-0)=-inf */
        }
        if (hx >> 31) > 0 {
            // We could return `f64::NAN` here, but then we would lose bit-to-bit compatibility with `libm::log2`.
            // The difference occurs e.g. for `x = NaN` (it still returns a NaN, but a different one).
            #[expect(clippy::eq_op)]
            return (x - x) / 0.0; /* log(-#) = NaN */
        }
        /* subnormal number, scale x up */
        k -= 54;
        x *= x1p54;
        ui = x.to_bits();

        #[expect(clippy::as_conversions)] // pre-shifted, so no overflow can occur
        let t = (ui >> 32) as u32;
        hx = t;
    } else if hx >= 0x7ff00000 {
        return x;
    } else if hx == 0x3ff00000 && ui << 32 == 0 {
        return 0.;
    }

    /* reduce x into [sqrt(2)/2, sqrt(2)] */
    hx += 0x3ff00000 - 0x3fe6a09e;

    #[expect(clippy::as_conversions, clippy::cast_possible_wrap)] // pre-shifted, so no overflow can occur
    let k_mod = (hx >> 20) as i32 - 0x3ff;
    k += k_mod;
    hx = (hx & 0x000fffff) + 0x3fe6a09e;
    ui = (u64::from(hx) << 32) | (ui & 0xffffffff);
    x = f64::from_bits(ui);

    let f = x - 1.0;
    let hfsq = 0.5 * f * f;
    let s = f / (2.0 + f);
    let z = s * s;
    let w = z * z;
    let t1 = w * (LG2 + w * (LG4 + w * LG6));
    let t2 = z * (LG1 + w * (LG3 + w * (LG5 + w * LG7)));
    let r = t2 + t1;

    /* hi+lo = f - hfsq + s*(hfsq+R) ~ log(1+f) */
    let mut hi = f - hfsq;
    ui = hi.to_bits();
    ui &= u64::MAX << 32;
    hi = f64::from_bits(ui);
    let lo = f - hi - hfsq + s * (hfsq + r);

    let mut val_hi = hi * IVLN2HI;
    let mut val_lo = (lo + hi) * IVLN2LO + lo * IVLN2HI;

    /* spadd(val_hi, val_lo, y), except for not using double_t: */
    let y: f64 = k.into();
    let w = y + val_hi;
    val_lo += (y - w) + val_hi;
    val_hi = w;

    val_lo + val_hi
}

/// A uint of the same width as the float
type F64Int = u64;

/// The bitwidth of the float type.
const F64_BITS: u32 = 64;

/// The bitwidth of the significand (does not include the implicit bit).
const F64_SIG_BITS: u32 = 52;

/// A mask for the sign bit
const F64_SIGN_MASK: F64Int = 1 << (F64_BITS - 1);

/// A mask for the significand
const F64_SIG_MASK: F64Int = (1 << F64_SIG_BITS) - 1;

/// The bitwidth of the exponent.
const F64_EXP_BITS: u32 = F64_BITS - F64_SIG_BITS - 1;

/// The saturated (maximum bitpattern) value of the exponent, i.e. the infinite
/// representation.
///
/// This shifted fully right, use `EXP_MASK` for the shifted value.
const F64_EXP_SAT: u32 = (1 << F64_EXP_BITS) - 1;

/// The exponent bias value.
const F64_EXP_BIAS: u32 = F64_EXP_SAT >> 1;

/// Constructs an `f64` from its parts. Inputs are treated as bits and shifted into position.
#[expect(clippy::as_conversions)] // converting `u32` to `u64` is safe
const fn f64_from_parts(negative: bool, exponent: u32, significand: F64Int) -> f64 {
    let sign: F64Int = if negative { 1 } else { 0 };
    f64::from_bits(
        (sign << (F64_BITS - 1))
            | (((exponent & F64_EXP_SAT) as F64Int) << F64_SIG_BITS)
            | (significand & F64_SIG_MASK),
    )
}

/// Returns the exponent, not adjusting for bias, not accounting for subnormals or zero.
#[expect(clippy::as_conversions)] // pre-shifted, so the conversion is safe
const fn f64_ex(x: f64) -> u32 {
    (x.to_bits() >> F64_SIG_BITS) as u32 & F64_EXP_SAT
}

/// Extract the exponent and adjust it for bias, not accounting for subnormals or zero.
#[expect(clippy::as_conversions, clippy::cast_possible_wrap)] // the constant fits into i32
const fn f64_exp_unbiased(x: f64) -> i32 {
    f64_ex(x) as i32 - (F64_EXP_BIAS as i32)
}

/* SPDX-License-Identifier: MIT
 * origin: musl src/math/trunc.c */
const fn trunc(x: f64) -> f64 {
    let xi: F64Int = x.to_bits();
    let e = f64_exp_unbiased(x);

    // The represented value has no fractional part, so no truncation is needed
    #[expect(clippy::as_conversions, clippy::cast_possible_wrap)] // the constant fits into i32
    if e >= F64_SIG_BITS as i32 {
        return x;
    }

    #[expect(clippy::as_conversions, clippy::cast_sign_loss)] // in the branch where `as` is used, `e` is non-negative
    let clear_mask = if e < 0 {
        // If the exponent is negative, the result will be zero so we clear everything
        // except the sign.
        !F64_SIGN_MASK
    } else {
        // Otherwise, we keep `e` fractional bits and clear the rest.
        F64_SIG_MASK >> (e as u32)
    };

    let cleared = xi & clear_mask;

    // Now zero the bits we need to truncate and return.
    f64::from_bits(xi ^ cleared)
}

pub(crate) const fn round(x: f64) -> f64 {
    let f0p5 = f64_from_parts(false, F64_EXP_BIAS - 1, 0); // 0.5
    let f0p25 = f64_from_parts(false, F64_EXP_BIAS - 2, 0); // 0.25

    trunc(x + (f0p5 - f0p25 * f64::EPSILON).copysign(x))
}

/* SPDX-License-Identifier: MIT
 * origin: musl src/math/floor.c */
/// Generic `floor` algorithm.
///
/// Note that this uses the algorithm from musl's `floorf` rather than `floor` or `floorl` because
/// performance seems to be better (based on icount) and it does not seem to experience rounding
/// errors on i386.
pub(crate) const fn floor(x: f64) -> f64 {
    let zero = 0;

    let mut ix = x.to_bits();
    let e = f64_exp_unbiased(x);

    // If the represented value has no fractional part, no truncation is needed.
    #[expect(clippy::as_conversions, clippy::cast_possible_wrap)] // constant fits into i32
    if e >= F64_SIG_BITS as i32 {
        return x;
    }

    if e >= 0 {
        // |x| >= 1.0
        #[expect(clippy::as_conversions, clippy::cast_sign_loss)] // checked that e >= 0
        let m = F64_SIG_MASK >> (e as u32);
        if ix & m == zero {
            // Portion to be masked is already zero; no adjustment needed.
            return x;
        }

        if x.is_sign_negative() {
            ix += m;
        }

        ix &= !m;
        f64::from_bits(ix)
    } else {
        // |x| < 1.0, zero or inexact with truncation

        if (ix & !F64_SIGN_MASK) == 0 {
            return x;
        }

        if x.is_sign_positive() {
            // 0.0 <= x < 1.0; rounding down goes toward +0.0.
            0.
        } else {
            // -1.0 < x < 0.0; rounding down goes toward -1.0.
            -1.
        }
    }
}

#[cfg(test)]
mod tests {
    use libm::{floor as libm_floor, log2 as libm_log2, round as libm_round};
    use proptest::prelude::*;

    use super::{floor, log2, round};

    proptest! {
        #[test]
        fn fuzzy_round(x in proptest::num::f64::ANY) {
            let test = round(x);
            let reference = libm_round(x);
            assert_eq!(test.to_bits(), reference.to_bits());
        }

        #[test]
        fn fuzzy_floor(x in proptest::num::f64::ANY) {
            let test = floor(x);
            let reference = libm_floor(x);
            assert_eq!(test.to_bits(), reference.to_bits());
        }

        #[test]
        fn fuzzy_log2(x in proptest::num::f64::ANY) {
            let test = log2(x);
            let reference = libm_log2(x);
            assert_eq!(test.to_bits(), reference.to_bits());
        }
    }
}
