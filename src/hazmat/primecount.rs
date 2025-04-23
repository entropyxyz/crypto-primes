use crate::hazmat::float::ln;
use core::f64;
use crypto_bigint::{Concat, NonZero, Split, Uint};
#[allow(unused_imports)]
use num_traits::float::FloatCore;

/// Estimate the number of primes smaller than x using the asymptotic expansion of Li(x) with 4 terms, i.e.:
///
///   `𝜋(𝑥) ~ x/ln x * (1 + 1!/ln x + 2!/ln^2 x + 3!/ln^3 x).`
///
/// # Panics
///
/// The number of limbs must be at least 2. For smaller values, use precalculated values for π(x) from
/// e.g. https://sweet.ua.pt/tos/primes.html.
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
///
/// # Sources
///
/// Wikipedia, [Prime Counting Function][wikipedia-pcf]
/// Pierre Dusart (2010), [ESTIMATES OF SOME FUNCTIONS OVER PRIMES WITHOUT R.H.][dusart2010].
/// Pierre Dusart (2018), [Explicit estimates of some functions over primes][dusart2018].
/// Lowell Schoenfeld (1976), [Sharper Bounds for the Chebyshev Functions θ(x) and ψ(x)](schoenfeld).
/// Trudgian, T. S. (2014). [Updating the error term in the prime number theorem](trudgian)
///
/// [dusart2010]: https://arxiv.org/pdf/1002.0442
/// [dusart2018]: https://www.researchgate.net/publication/309522478_Explicit_estimates_of_some_functions_over_primes
/// [wikipedia-pcf]: https://en.wikipedia.org/wiki/Prime-counting_function
/// [schoenfeld]: https://www.jstor.org/stable/2005976
/// [trudgian]: https://arxiv.org/abs/1401.2689
///
pub fn estimate_pi_x<const LIMBS: usize, const RHS_LIMBS: usize>(x: &Uint<LIMBS>) -> Uint<LIMBS>
where
    Uint<LIMBS>: Concat<Output = Uint<RHS_LIMBS>>,
    Uint<RHS_LIMBS>: Split<Output = Uint<LIMBS>>,
{
    #[cfg(target_pointer_width = "64")]
    assert!(
        LIMBS >= 2,
        "LIMBS must be at least 2; for smaller values, use precalculated values for π(x)"
    );
    #[cfg(target_pointer_width = "32")]
    assert!(
        LIMBS >= 4,
        "LIMBS must be at least 4; for smaller values, use precalculated values for π(x)"
    );
    // Number of bits to scale by (fractional bits during division)
    const SCALE_BITS: u32 = 64;
    const TOTAL_SCALE_BITS: u32 = 2 * SCALE_BITS;
    // Scaling factor for the denominators of the expansion terms, 2^64.
    const DENOM_SCALE_FACTOR: f64 = (1u128 << SCALE_BITS) as f64;

    let ln_x = ln(x);

    if !ln_x.is_finite() || ln_x <= 1.0 {
        return Uint::ZERO;
    }

    // Calculate f64 denominators L, L^2, L^3, L^4 and round to u64
    let ln_x_2 = ln_x * ln_x;
    let ln_x_3 = ln_x_2 * ln_x;
    let ln_x_4 = ln_x_3 * ln_x;

    // Scale up a float by 2^64, round, cast to u128 and then to a wide `Uint`.
    let f64_to_scaled_uint = |value: f64| -> NonZero<Uint<RHS_LIMBS>> {
        let scaled = value * DENOM_SCALE_FACTOR;
        let scaled = (scaled.round()).max(1.0) as u128;
        let denom = Uint::<RHS_LIMBS>::from_u128(scaled);
        NonZero::new(denom).expect("max(1.0) ensures value is at least 1")
    };

    let d1 = f64_to_scaled_uint(ln_x);
    let d2 = f64_to_scaled_uint(ln_x_2);
    let d3 = f64_to_scaled_uint(ln_x_3);
    let d4 = f64_to_scaled_uint(ln_x_4);

    let x_wide: Uint<RHS_LIMBS> = x.concat(&Uint::ZERO);

    let term1_scaled = (x_wide << TOTAL_SCALE_BITS).wrapping_div(&d1);
    let term2_scaled = (x_wide << TOTAL_SCALE_BITS).wrapping_div(&d2);

    // Term 3: (2x << TOTAL_SCALE_BITS) / d3  (Factorial 2!)
    let term3_scaled = (x_wide << (TOTAL_SCALE_BITS + 1)).wrapping_div(&d3);

    // Term 4: (6x << TOTAL_SCALE_BITS) / d4  (Factorial 3!)
    let six = Uint::<LIMBS>::from(6u64);
    let x_times_6 = x_wide.saturating_mul(&six);
    let term4_scaled = (x_times_6 << TOTAL_SCALE_BITS).wrapping_div(&d4);

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

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use alloc::vec::Vec;
    use crypto_bigint::{U1024, U256, U64};

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
        // Want the delta to be "far away" from Sage's estimate.
        const MIN_BIT_DIFFERENCE: u32 = 29;
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

        // Want the delta to be "far away" from Sage's estimate.
        const MIN_BIT_DIFFERENCE: u32 = 33;
        assert!(
            sage_est.bits_vartime() - delta.bits_vartime() >= MIN_BIT_DIFFERENCE,
            "Estimate not close enough: delta has {} bits, expected {} bits. Difference should be >= {}; delta: {delta}, sage_est: {sage_est}, estimate: {estimate}",
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
        // Want the delta to be "far away" from Sage's estimate.
        const MIN_BIT_DIFFERENCE: u32 = 34;
        assert!(
            sage_est.bits_vartime() - delta.bits_vartime() >= MIN_BIT_DIFFERENCE,
            "Estimate not close enough: delta has {} bits, expected {} bits. Difference should be >= {}",
            delta.bits_vartime(),
            sage_est.bits_vartime(),
            MIN_BIT_DIFFERENCE
        );
    }

    #[test]
    #[should_panic(expected = "LIMBS must be at least ")]
    fn pi_x_tiny() {
        let x = U64::from_u64(10000000);
        estimate_pi_x(&x);
    }

    #[test]
    fn pi_x_estimates_for_known_values() {
        // List of known values for π(x), expressed as a tuple of `(π(x), exponent)`, where the `exponent` is used with base 10.
        let pi_xs: Vec<(u128, u32)> = vec![
            // The error is large for small x, so we skip them.
            // (4, 1),
            // (25, 2),
            // (168, 3),
            (1229, 4),
            (9592, 5),
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
