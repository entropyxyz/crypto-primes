use crypto_bigint::Uint;
#[allow(unused_imports)]
use num_traits::float::FloatCore as _;

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
/// | Bits | Error Bound (abs) | Error Bound (rel) 		|
/// | :--- | :---------------- | :--------------------- |
/// | 1024 | < 2^517           | < 2^-497               |
/// | 2048 | < 2^1030          | < 2^-1007              |
/// | 4096 | < 2^2055          | < 2^-2029              |
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
/// It should be noted that while Li(x) is generally smaller than π(x) for 'small' x, it is known that the sign of `π(x)
/// - Li(x)` changes infinitely often. It has been proven that there must be a crossing below ~10^316 (~2^1051), which
/// is well within the ranges used in this library. Thus, users should be aware that the estimate provided here can be
/// both greater than and smaller than the actual value of π(x).
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
