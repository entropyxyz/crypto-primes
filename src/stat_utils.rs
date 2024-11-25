use alloc::{vec, vec::Vec};
use crypto_bigint::{Bounded, Constants, Integer, NonZero, RandomBits, RandomMod};
use statrs::distribution::{ChiSquared, ContinuousCDF};
use tracing::{debug, info};

use crate::hazmat::UniformSieve;

/// Calculate the distribution quality for a prime generator. In this example the we work with the
/// [`UniformSieve`].
pub fn check_distribution_quality<T>()
where
    T: Integer + Copy + Bounded + Constants + RandomBits + RandomMod + UniformSieve<T>,
{
    let sample_count = 1000;
    let num_intervals = 10;
    let expected_count = sample_count as f64 / num_intervals as f64;

    let primes = collect_primes::<T>(sample_count);
    let counts = distribution(primes, num_intervals);
    debug!("Found {counts:?} primes in each interval. Expected count per interval={expected_count}");

    // Chi-Square test for uniformity
    let (chi_square_stat, p_value) = chi_square_test(&counts, expected_count);
    info!(
        "Sampled {sample_count} {} bit primes. p={p_value}, chi_squared={chi_square_stat}",
        T::BITS
    );
    assert!(p_value > 0.05, "Expected a p_value larger than 0.05, instead found {p_value}, failing to prove that the distribution is uniform.");
    assert!(chi_square_stat < 16.92, "For p 0.05, the critical chi-square value is 16.92, instead we got {chi_square_stat}, failing the hypothesis that the distribution is uniform.");
}

// Count primes in each subinterval. Returns stats.
fn distribution<T>(primes: Vec<T>, intervals: usize) -> Vec<u64>
where
    T: Integer + Constants,
{
    let interval_width = T::MAX.div(NonZero::new(T::from(intervals as u64)).unwrap());
    let interval_width_nz = NonZero::new(interval_width).unwrap();
    let mut counts = vec![0; intervals];

    for prime in &primes {
        let interval = prime.clone().div(interval_width_nz.clone()).as_ref()[0].0 as usize;
        if interval < intervals {
            counts[interval] += 1;
        }
    }

    counts
}

// Collect primes
fn collect_primes<T>(num_primes: usize) -> Vec<T>
where
    T: Integer + Copy + Bounded + Constants + RandomBits + RandomMod + UniformSieve<T>,
{
    (0..num_primes).fold(Vec::with_capacity(num_primes), |mut acc, _| {
        acc.push(T::generate_prime());
        acc
    })
}
// Calculate Chi-Square statistic and p-value
fn chi_square_test(counts: &[u64], expected: f64) -> (f64, f64) {
    let chi_square_stat: f64 = counts
        .iter()
        .map(|&count| {
            let diff = count as f64 - expected;
            diff * diff / expected
        })
        .sum();

    // Calculate p-value
    let chi_square_dist = ChiSquared::new((counts.len() - 1) as f64).unwrap();
    let p_value = 1.0 - chi_square_dist.cdf(chi_square_stat);

    (chi_square_stat, p_value)
}
