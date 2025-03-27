use crypto_bigint::Uint;

/// Estimate the number of primes smaller than x using the asymptotic expansion of Li(x) with three terms.
/// 𝜋(𝑥)/Li(x) converges to 1 from below, which means that the estimation of number of primes for the ranges relevant
/// here is a slight overestimate.
/// The numerical approximations involved – in particular log₂ x ~ floor(log₂ x) – means that the approximation for ln x
/// is an underestimate, which means the approximated Li(x) is >= actual Li(x).
/// Uses the formula 𝜋(𝑥) ~ x/ln x * (1 + 1!/ln x + 2!/ln^2 x + 3!/ln^3 x).
pub fn estimate_pi_x<const LIMBS: usize>(x: &Uint<LIMBS>) -> Uint<LIMBS> {
    let ln_x = ln(x);
    let term1 = x / Uint::from_u64(ln_x.round() as u64);
    let ln_x_sq = ln_x * ln_x;
    let term2 = (1.0 + 1.0 / ln_x + 2.0 / (ln_x_sq) + 6.0 / (ln_x_sq * ln_x)).round() as u64;
    term1 * Uint::<LIMBS>::from_u64(term2)
}

// Calculate the natural logarithm of a big integer using the relation ln(x) = log₂(x) / log₂(e).
#[inline]
fn ln<const LIMBS: usize>(x: &Uint<LIMBS>) -> f64 {
    if x == &Uint::ONE {
        return 0.0;
    }
    let log2_x = x.bits_vartime().saturating_sub(1) as f64;
    log2_x * core::f64::consts::LN_2
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use alloc::vec::Vec;
    use crypto_bigint::U256;

    #[test]
    fn pi_estimates() {
        let pi_xs: Vec<(u128, u32)> = vec![
            // The error is large for small x, so we skip them.
            // (4, 1),
            // (25, 2),
            // (168, 3),
            // (1229, 4),
            // (9592, 5),
            // (78498, 6),
            // (664579, 7),
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
            let delta = uint_to_u128(&delta);
            let approx_pi_x_128 = uint_to_u128(&approx_pi_x);
            let error = (delta as f64 / *pi_x as f64) * 100.0;
            assert!(
                error < 5.0, // For large x, this error is much better, well below 1.
                "10^{exponent}:\t{pi_x} - {approx_pi_x_128} = {delta}, err: {error:.2}"
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
        return (limbs[1].0 as u128) << 64 | limbs[0].0 as u128;
    }
}
