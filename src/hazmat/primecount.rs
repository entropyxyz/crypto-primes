use crypto_bigint::{Uint, U256};

/// Estimate the number of primes smaller than x using the asymptotic expansion of Li(x) with three terms. 𝜋(𝑥)/Li(x)
/// converges to 1 from below, which means that the estimation of number of primes for the ranges relevant here is a
/// slight overestimate.
/// Uses the formula 𝜋(𝑥) ~ (x/ln(x)) * (1 + 1!/ln(x) + 2!/(ln(x)^2) + 3!/(ln(x)^3)).
pub fn estimate_pi_x<const LIMBS: usize>(x: &Uint<LIMBS>) -> Uint<LIMBS> {
    let ln_x = ln(x);
    let term1 = x / Uint::from_u64(ln_x);
    let ln_x = ln_x as f64;
    let ln_x_sq = ln_x * ln_x;
    let term2 = (1.0 + 1.0 / ln_x + 2.0 / (ln_x_sq) + 6.0 / (ln_x_sq * ln_x)) as u64;
    term1 * Uint::<LIMBS>::from_u64(term2)
}

// Calculate the natural logarithm of a big integer using the relation ln(x) = log₂(x) / log₂(e).
// Returns a truncated integer representation of ln(x), which means this is a "decent" approximation for large integers
// in the range of 2^64 to 2^8192. Smaller integers will be less accurate.
fn ln<const LIMBS: usize>(x: &Uint<LIMBS>) -> u64 {
    // ln(2) multiplied by 2^128
    const SCALED_LN2: U256 = U256::from_be_hex("00000000000000000000000000000000B17217F7D1CF79ABC9E3B39803F2F6AF");

    let log2_x = U256::from(x.bits_vartime().saturating_sub(1));
    // The error in this approximation is "lower 128 bits of the product"/2^128, which is on the order of 10^-40.
    let ln_x: U256 = (log2_x * SCALED_LN2) >> 128;
    #[cfg(target_pointer_width = "32")]
    let ln_x = (ln_x.as_limbs()[1].0 as u64) << 32 | ln_x.as_limbs()[0].0 as u64;
    #[cfg(target_pointer_width = "64")]
    let ln_x = ln_x.as_limbs()[0].0;
    ln_x & !((ln_x == 0) as u64) // Ensure we return 0 for ln(1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use alloc::vec::Vec;
    use crypto_bigint::U256;

    #[test]
    fn natural_logarithm() {
        for exponent in 1..=35 {
            let x = U256::from_u128(10u128.pow(exponent));
            let ln_x = ln(&x);
            let expected_ln_x = (exponent as f64 * core::f64::consts::LN_10) as u64; // ln(10^exponent) = exponent * ln(10)
            let delta = ln_x.abs_diff(expected_ln_x);
            assert!(
                delta < 2,
                "ln(10^{exponent}) = {ln_x}, expected {expected_ln_x}, delta {delta}"
            );
        }
    }

    // TODO(dp): Need a variant for 32 bit, or at least a helper function to extract u64 from limbs.
    #[cfg(target_pointer_width = "64")]
    #[test]
    fn pi_estimates() {
        let pi_xs: Vec<(u128, u32)> = vec![
            // The error is very large for small x, so we skip them.
            // (4, 1),
            // (25, 2),
            // (168, 3),
            // (1229, 4),
            // (9592, 5),
            // (78498, 6),
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
            let approx_pi_x = estimate_pi_x(&n);
            let delta = if pi_x_wide > approx_pi_x {
                pi_x_wide - approx_pi_x
            } else {
                approx_pi_x - pi_x_wide
            };
            let delta = (delta.as_limbs()[1].0 as u128) << 64 | delta.as_limbs()[0].0 as u128;
            let approx_pi_x_128 = (approx_pi_x.as_limbs()[1].0 as u128) << 64 | approx_pi_x.as_limbs()[0].0 as u128;
            let error = (delta as f64 / *pi_x as f64) * 100.0;
            assert!(
                error < 5.0, // For large x, this error is much better, well below 1.
                "10^{exponent}:\t{pi_x} - {approx_pi_x_128} = {delta}, err: {error:.2}"
            );
        }
    }
}
