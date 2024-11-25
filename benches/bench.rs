use core::num::NonZero;

use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use crypto_bigint::{nlimbs, BoxedUint, Integer, Odd, RandomBits, Uint, U1024, U128, U2048, U256, U512};
use rand_chacha::ChaCha8Rng;
use rand_core::{CryptoRngCore, OsRng, SeedableRng};

#[cfg(feature = "tests-gmp")]
use rug::{integer::Order, Integer as GmpInteger};

#[cfg(feature = "tests-openssl")]
use openssl::bn::BigNum;

use crypto_primes::{
    generate_prime_with_rng, generate_safe_prime_with_rng,
    hazmat::{
        lucas_test, random_odd_integer, AStarBase, BruteForceBase, LucasCheck, MillerRabin, SelfridgeBase, Sieve,
    },
    is_prime_with_rng, is_safe_prime_with_rng,
};
#[cfg(feature = "multicore")]
use crypto_primes::{par_generate_prime_with_rng, par_generate_safe_prime_with_rng};

fn make_rng() -> ChaCha8Rng {
    ChaCha8Rng::from_seed(*b"01234567890123456789012345678901")
}

fn random_odd_uint<T: RandomBits + Integer>(rng: &mut impl CryptoRngCore, bit_length: u32) -> Odd<T> {
    random_odd_integer::<T>(rng, NonZero::new(bit_length).unwrap())
}

fn make_sieve<const L: usize>(rng: &mut impl CryptoRngCore) -> Sieve<Uint<L>> {
    let start = random_odd_uint::<Uint<L>>(rng, Uint::<L>::BITS);
    Sieve::new(start.get(), NonZero::new(Uint::<L>::BITS).unwrap(), false)
}

fn make_presieved_num<const L: usize>(rng: &mut impl CryptoRngCore) -> Odd<Uint<L>> {
    let mut sieve = make_sieve(rng);
    Odd::new(sieve.next().unwrap()).unwrap()
}

fn bench_uniform_sieve(c: &mut Criterion) {
    use crypto_primes::uniform_sieve::UniformSieve;
    let mut group = c.benchmark_group("Uniform sieve");

    let mut rng = make_rng();
    group.bench_function("(U128) Random prime", |b| {
        b.iter(|| U128::generate_prime_with_rng(&mut rng));
    });

    let mut rng = make_rng();
    group.bench_function("(U256) Random prime", |b| {
        b.iter(|| U256::generate_prime_with_rng(&mut rng));
    });

    let mut rng = make_rng();
    group.bench_function("(U512) Random prime", |b| {
        b.iter(|| U256::generate_prime_with_rng(&mut rng));
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) Random prime", |b| {
        b.iter(|| U1024::generate_prime_with_rng(&mut rng));
    });

    let mut rng = make_rng();
    group.bench_function("(U2048) Random prime", |b| {
        b.iter(|| U2048::generate_prime_with_rng(&mut rng));
    });
}
fn bench_sieve(c: &mut Criterion) {
    let mut group = c.benchmark_group("Sieve");

    group.bench_function("(U128) random start", |b| {
        b.iter(|| random_odd_uint::<U128>(&mut OsRng, 128))
    });

    group.bench_function("(U128) creation", |b| {
        b.iter_batched(
            || random_odd_uint::<U128>(&mut OsRng, 128),
            |start| Sieve::new(start.get(), NonZero::new(128).unwrap(), false),
            BatchSize::SmallInput,
        )
    });

    // 5 is the average number of pre-sieved samples we need to take before we encounter a prime
    group.bench_function("(U128) average sieve samples for a prime (5)", |b| {
        b.iter_batched(
            || make_sieve::<{ nlimbs!(128) }>(&mut OsRng),
            |sieve| sieve.take(5).for_each(drop),
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U1024) random start", |b| {
        b.iter(|| random_odd_uint::<U1024>(&mut OsRng, 1024))
    });

    group.bench_function("(U1024) creation", |b| {
        b.iter_batched(
            || random_odd_uint::<U1024>(&mut OsRng, 1024),
            |start| Sieve::new(start.get(), NonZero::new(1024).unwrap(), false),
            BatchSize::SmallInput,
        )
    });

    // 42 is the average number of pre-sieved samples we need to take before we encounter a prime
    group.bench_function("(U1024) average sieve samples for a prime (42)", |b| {
        b.iter_batched(
            || make_sieve::<{ nlimbs!(1024) }>(&mut OsRng),
            |sieve| sieve.take(42).for_each(drop),
            BatchSize::SmallInput,
        )
    });

    // 42^2 is the average number of pre-sieved samples we need to take
    // before we encounter a safe prime
    group.bench_function("(U1024) average sieve samples for a safe prime (42^2)", |b| {
        b.iter_batched(
            || make_sieve::<{ nlimbs!(1024) }>(&mut OsRng),
            |sieve| sieve.take(42 * 42).for_each(drop),
            BatchSize::SmallInput,
        )
    });

    group.finish()
}

fn bench_miller_rabin(c: &mut Criterion) {
    let mut group = c.benchmark_group("Miller-Rabin");

    group.bench_function("(U128) creation", |b| {
        b.iter_batched(
            || random_odd_uint::<U128>(&mut OsRng, 128),
            MillerRabin::<U128>::new,
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U128) random base test (pre-sieved)", |b| {
        b.iter_batched(
            || MillerRabin::new(make_presieved_num::<{ nlimbs!(128) }>(&mut OsRng)),
            |mr| mr.test_random_base(&mut OsRng),
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U1024) creation", |b| {
        b.iter_batched(
            || random_odd_uint::<U1024>(&mut OsRng, 1024),
            MillerRabin::<U1024>::new,
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U1024) random base test (pre-sieved)", |b| {
        b.iter_batched(
            || MillerRabin::new(make_presieved_num::<{ nlimbs!(1024) }>(&mut OsRng)),
            |mr| mr.test_random_base(&mut OsRng),
            BatchSize::SmallInput,
        )
    });
}

fn bench_lucas(c: &mut Criterion) {
    let mut group = c.benchmark_group("Lucas");

    let mut rng = make_rng();
    group.bench_function("(U128) Selfridge base, strong check (pre-sieved)", |b| {
        b.iter_batched(
            || make_presieved_num::<{ nlimbs!(128) }>(&mut rng),
            |n| lucas_test(n, SelfridgeBase, LucasCheck::Strong),
            BatchSize::SmallInput,
        )
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) Selfridge base, strong check (pre-sieved)", |b| {
        b.iter_batched(
            || make_presieved_num::<{ nlimbs!(1024) }>(&mut rng),
            |n| lucas_test(n, SelfridgeBase, LucasCheck::Strong),
            BatchSize::SmallInput,
        )
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) A* base, Lucas-V check (pre-sieved)", |b| {
        b.iter_batched(
            || make_presieved_num::<{ nlimbs!(1024) }>(&mut rng),
            |n| lucas_test(n, AStarBase, LucasCheck::LucasV),
            BatchSize::SmallInput,
        )
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) brute force base, almost extra strong (pre-sieved)", |b| {
        b.iter_batched(
            || make_presieved_num::<{ nlimbs!(1024) }>(&mut rng),
            |n| lucas_test(n, BruteForceBase, LucasCheck::AlmostExtraStrong),
            BatchSize::SmallInput,
        )
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) brute force base, extra strong (pre-sieved)", |b| {
        b.iter_batched(
            || make_presieved_num::<{ nlimbs!(1024) }>(&mut rng),
            |n| lucas_test(n, BruteForceBase, LucasCheck::ExtraStrong),
            BatchSize::SmallInput,
        )
    });

    // A number triggering a slow path through the Lucas test, invoking all the available checks:
    // - U_{d} checked, but is not 0
    // - V_{d * 2^t} checked for t == 0..s-1, but no V = 0 found
    // - s = 5, so the previous step has multiple checks
    // - Q != 1 (since we're using Selfridge base)
    let slow_path = Odd::new(U1024::from_be_hex(concat![
        "D1CB9F1B6F3414A4B40A7E51C53C6AE4689DFCDC49FF875E7066A229D704EA8E",
        "6B674231D8C5974001673C3CE7FF9D377C8564E5182165A23434BC7B7E6C0419",
        "FD25C9921B0E9C90AF2570DB0772E1A9C82ACABBC8FC0F0864CE8A12124FA29B",
        "7F870924041DFA13EE5F5541C1BF96CA679EFAE2C96F5F4E9DF6007185198F5F"
    ]))
    .unwrap();

    group.bench_function("(U1024) Selfridge base, strong check, slow path", |b| {
        b.iter(|| {
            lucas_test(slow_path, SelfridgeBase, LucasCheck::Strong);
        })
    });

    group.finish();
}

fn bench_presets(c: &mut Criterion) {
    let mut group = c.benchmark_group("Presets");

    group.bench_function("(U128) Prime test", |b| {
        b.iter_batched(
            || random_odd_uint::<U128>(&mut OsRng, 128),
            |num| is_prime_with_rng(&mut OsRng, num.as_ref()),
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U128) Safe prime test", |b| {
        b.iter_batched(
            || random_odd_uint::<U128>(&mut OsRng, 128),
            |num| is_safe_prime_with_rng(&mut OsRng, num.as_ref()),
            BatchSize::SmallInput,
        )
    });

    let mut rng = make_rng();
    group.bench_function("(U128) Random prime", |b| {
        b.iter(|| generate_prime_with_rng::<U128>(&mut rng, 128))
    });

    let mut rng = make_rng();
    group.bench_function("(U256) Random prime", |b| {
        b.iter(|| generate_prime_with_rng::<U256>(&mut rng, 128))
    });

    let mut rng = make_rng();
    group.bench_function("(U512) Random prime", |b| {
        b.iter(|| generate_prime_with_rng::<U512>(&mut rng, 128))
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) Random prime", |b| {
        b.iter(|| generate_prime_with_rng::<U1024>(&mut rng, 1024))
    });

    let mut rng = make_rng();
    group.bench_function("(U2048) Random prime", |b| {
        b.iter(|| generate_prime_with_rng::<U2048>(&mut rng, 1024))
    });

    let mut rng = make_rng();
    group.bench_function("(U128) Random safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<U128>(&mut rng, 128))
    });

    group.sample_size(20);
    let mut rng = make_rng();
    group.bench_function("(U1024) Random safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<U1024>(&mut rng, 1024))
    });

    let mut rng = make_rng();
    group.bench_function("(Boxed128) Random safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<BoxedUint>(&mut rng, 128))
    });

    group.sample_size(20);
    let mut rng = make_rng();
    group.bench_function("(Boxed1024) Random safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<BoxedUint>(&mut rng, 1024))
    });

    group.finish();

    // A separate group for bounded tests, to make it easier to run them separately.
    let mut group = c.benchmark_group("Presets (bounded)");

    let mut rng = make_rng();
    group.bench_function("(U128) Random safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<U128>(&mut rng, 128))
    });

    // The performance should scale with the prime size, not with the Uint size.
    // So we should strive for this test's result to be as close as possible
    // to that of the previous one and as far away as possible from the next one.
    group.bench_function("(U256) Random 128 bit safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<U256>(&mut rng, 128))
    });

    // The upper bound for the previous test.
    group.bench_function("(U256) Random 256 bit safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<U256>(&mut rng, 256))
    });

    group.finish();
}

#[cfg(feature = "multicore")]
fn bench_multicore_presets(c: &mut Criterion) {
    let mut group = c.benchmark_group("Presets (multicore)");
    let mut rng = make_rng();
    group.bench_function("(U128) Random prime", |b| {
        b.iter(|| par_generate_prime_with_rng::<U128>(&mut rng, 128, num_cpus::get()))
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) Random prime", |b| {
        b.iter(|| par_generate_prime_with_rng::<U1024>(&mut rng, 1024, num_cpus::get()))
    });

    let mut rng = make_rng();
    group.bench_function("(U128) Random safe prime", |b| {
        b.iter(|| par_generate_safe_prime_with_rng::<U128>(&mut rng, 128, num_cpus::get()))
    });

    group.sample_size(20);
    let mut rng = make_rng();
    group.bench_function("(U1024) Random safe prime", |b| {
        b.iter(|| par_generate_safe_prime_with_rng::<U1024>(&mut rng, 1024, num_cpus::get()))
    });

    let mut rng = make_rng();
    group.bench_function("(Boxed128) Random safe prime", |b| {
        b.iter(|| par_generate_safe_prime_with_rng::<BoxedUint>(&mut rng, 128, num_cpus::get()))
    });

    group.sample_size(20);
    let mut rng = make_rng();
    group.bench_function("(Boxed1024) Random safe prime", |b| {
        b.iter(|| par_generate_safe_prime_with_rng::<BoxedUint>(&mut rng, 1024, num_cpus::get()))
    });
}
#[cfg(not(feature = "multicore"))]
fn bench_multicore_presets(_c: &mut Criterion) {}

#[cfg(feature = "tests-gmp")]
fn bench_gmp(c: &mut Criterion) {
    let mut group = c.benchmark_group("GMP");

    fn random<const L: usize>(rng: &mut impl CryptoRngCore) -> GmpInteger {
        let num = random_odd_uint::<Uint<L>>(rng, Uint::<L>::BITS).get();
        GmpInteger::from_digits(num.as_words(), Order::Lsf)
    }

    group.bench_function("(U128) Random prime", |b| {
        b.iter_batched(
            || random::<{ nlimbs!(128) }>(&mut OsRng),
            |n| n.next_prime(),
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U1024) Random prime", |b| {
        b.iter_batched(
            || random::<{ nlimbs!(1024) }>(&mut OsRng),
            |n| n.next_prime(),
            BatchSize::SmallInput,
        )
    });

    group.finish();
}

#[cfg(not(feature = "tests-gmp"))]
fn bench_gmp(_c: &mut Criterion) {}

#[cfg(feature = "tests-openssl")]
fn bench_openssl(c: &mut Criterion) {
    let mut group = c.benchmark_group("OpenSSL");

    group.bench_function("(U128) Random prime", |b| {
        b.iter_batched(
            || BigNum::new().unwrap(),
            |p| {
                let mut p = p;
                p.generate_prime(128, false, None, None).unwrap()
            },
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U1024) Random prime", |b| {
        b.iter_batched(
            || BigNum::new().unwrap(),
            |p| {
                let mut p = p;
                p.generate_prime(1024, false, None, None).unwrap()
            },
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U128) Random safe prime", |b| {
        b.iter_batched(
            || BigNum::new().unwrap(),
            |p| {
                let mut p = p;
                p.generate_prime(128, true, None, None).unwrap()
            },
            BatchSize::SmallInput,
        )
    });

    group.sample_size(20);
    group.bench_function("(U1024) Random safe prime", |b| {
        b.iter_batched(
            || BigNum::new().unwrap(),
            |p| {
                let mut p = p;
                p.generate_prime(1024, true, None, None).unwrap()
            },
            BatchSize::SmallInput,
        )
    });

    group.finish();
}

#[cfg(not(feature = "tests-openssl"))]
fn bench_openssl(_c: &mut Criterion) {}

#[cfg(feature = "tests-glass-pumpkin")]
fn bench_glass_pumpkin(c: &mut Criterion) {
    use crypto_bigint::Limb;
    use crypto_primes::hazmat::{lucas_test, AStarBase, LucasCheck, MillerRabin, Primality};

    // The `glass-pumpkin` implementation is doing a different number of M-R checks than this crate.
    // For a fair comparison we make a custom implementation here,
    // using the same number of checks that `glass-pumpkin` does.
    fn required_checks(bits: u32) -> usize {
        ((bits as f64).log2() as usize) + 5
    }

    // Mimics the sequence of checks `glass-pumpkin` does to find a prime.
    fn prime_like_gp(bit_length: u32, rng: &mut impl CryptoRngCore) -> BoxedUint {
        loop {
            let start = random_odd_integer::<BoxedUint>(rng, NonZero::new(bit_length).unwrap()).get();
            let sieve = Sieve::new(start, NonZero::new(bit_length).unwrap(), false);
            for num in sieve {
                let odd_num = Odd::new(num.clone()).unwrap();

                let mr = MillerRabin::new(odd_num.clone());
                if (0..required_checks(bit_length)).any(|_| !mr.test_random_base(rng).is_probably_prime()) {
                    continue;
                }

                match lucas_test(odd_num, AStarBase, LucasCheck::Strong) {
                    Primality::Composite => continue,
                    Primality::Prime => return num,
                    _ => {}
                }

                return num;
            }
        }
    }

    // Mimics the sequence of checks `glass-pumpkin` does to find a safe prime.
    fn safe_prime_like_gp(bit_length: u32, rng: &mut impl CryptoRngCore) -> BoxedUint {
        loop {
            let start = random_odd_integer::<BoxedUint>(rng, NonZero::new(bit_length).unwrap()).get();
            let sieve = Sieve::new(start, NonZero::new(bit_length).unwrap(), true);
            for num in sieve {
                let odd_num = Odd::new(num.clone()).unwrap();

                let limbs: &[Limb] = num.as_ref();
                if limbs[0].0 & 3 != 3 {
                    continue;
                }

                let half = num.wrapping_shr_vartime(1);
                let odd_half = Odd::new(half.clone()).unwrap();

                let checks = required_checks(bit_length) - 5;

                let mr = MillerRabin::new(odd_num.clone());
                if (0..checks).any(|_| !mr.test_random_base(rng).is_probably_prime()) {
                    continue;
                }

                if lucas_test(odd_num, AStarBase, LucasCheck::Strong) == Primality::Composite {
                    continue;
                }

                let mr = MillerRabin::new(odd_half.clone());
                if (0..checks).any(|_| !mr.test_random_base(rng).is_probably_prime()) {
                    continue;
                }

                match lucas_test(odd_half, AStarBase, LucasCheck::Strong) {
                    Primality::Composite => continue,
                    Primality::Prime => return num,
                    _ => {}
                }

                return num;
            }
        }
    }

    let mut group = c.benchmark_group("glass-pumpkin");

    let mut rng = make_rng();
    group.bench_function("(U1024) Random prime (crypto-primes default)", |b| {
        b.iter(|| generate_prime_with_rng::<BoxedUint>(&mut rng, 1024))
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) Random prime (crypto-primes mimicking glass-pumpkin)", |b| {
        b.iter(|| prime_like_gp(1024, &mut rng))
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) Random prime", |b| {
        b.iter(|| glass_pumpkin::prime::from_rng(1024, &mut rng))
    });

    group.sample_size(20);

    let mut rng = make_rng();
    group.bench_function("(U1024) Random safe prime (crypto-primes default)", |b| {
        b.iter(|| generate_safe_prime_with_rng::<BoxedUint>(&mut rng, 1024))
    });

    let mut rng = make_rng();
    group.bench_function(
        "(U1024) Random safe prime (crypto-primes mimicking glass-pumpkin)",
        |b| b.iter(|| safe_prime_like_gp(1024, &mut rng)),
    );

    let mut rng = make_rng();
    group.bench_function("(U1024) Random safe prime", |b| {
        b.iter(|| glass_pumpkin::safe_prime::from_rng(1024, &mut rng))
    });
}

#[cfg(not(feature = "tests-glass-pumpkin"))]
fn bench_glass_pumpkin(_c: &mut Criterion) {}

criterion_group!(
    benches,
    bench_uniform_sieve,
    bench_sieve,
    bench_miller_rabin,
    bench_lucas,
    bench_presets,
    bench_multicore_presets,
    bench_gmp,
    bench_openssl,
    bench_glass_pumpkin,
);
criterion_main!(benches);
