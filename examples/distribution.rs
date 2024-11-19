use crypto_bigint::{FixedInteger, NonZero, RandomBits, RandomMod};
use crypto_primes::uniform_sieve::UniformGeneratePrime;
use statrs::distribution::{ChiSquared, ContinuousCDF};

/// Chi-squared statistical test for some collection of primes.

fn main() {
    let sample_count = 5000;
    let num_intervals = 20;
    type T = crypto_bigint::U512;

    let end = T::MAX;
    let (counts, elapsed) = distribution::<T>(sample_count, num_intervals, end);
    let expected_count = sample_count as f64 / num_intervals as f64;

    // Chi-Square test for uniformity
    let (chi_square_stat, p_value) = chi_square_test(&counts, expected_count);

    println!("Prime Counts in Each Interval: {:?}", counts);
    println!("Collected {sample_count} primes in {elapsed:?}");
    println!("Expected Count per Interval: {:.2}", expected_count);
    println!("Chi-Square Statistic: {:.2}", chi_square_stat);
    println!("P-Value: {:.4}", p_value);

    if p_value < 0.05 {
        println!("The prime distribution significantly deviates from uniformity.");
    } else {
        println!("The prime distribution is close to uniform.");
    }
}

// Divide interval and count primes in each subinterval. Returns stats and the time elapsed.
fn distribution<T>(num_primes: usize, intervals: usize, end: T) -> (Vec<u64>, std::time::Duration)
where
    T: FixedInteger + RandomBits + RandomMod,
    T: UniformGeneratePrime<T>,
{
    let now = std::time::Instant::now();
    let primes: Vec<T> = collect_primes(num_primes);
    let elapsed = now.elapsed();
    println!("primes={:?}", &primes[0..5]);
    let start = T::ZERO;
    assert!(end > start);
    let interval_width = (end - start).div(NonZero::new(T::from(intervals as u64)).unwrap());
    let interval_width_nz = NonZero::new(interval_width).unwrap();
    println!("intervals_width={interval_width:?}");
    // Initialize counters for each interval
    let mut counts = vec![0; intervals];

    // Count primes in each interval
    for &prime in &primes {
        let interval = (prime - start).div(interval_width_nz).as_ref()[0].0 as usize;
        if interval < intervals {
            counts[interval] += 1;
        } else {
            println!("interval={interval} is outside the intervals count {intervals}");
        }
    }

    (counts, elapsed)
}

// Collect primes
fn collect_primes<T>(num_primes: usize) -> Vec<T>
where
    T: FixedInteger + RandomBits + RandomMod,
    T: UniformGeneratePrime<T>,
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
