use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use crypto_bigint::{nlimbs, Uint, U1024};
use rand_chacha::ChaCha8Rng;
use rand_core::{CryptoRngCore, OsRng, SeedableRng};

#[cfg(feature = "tests-gmp")]
use rug::{integer::Order, Integer};

#[cfg(feature = "tests-openssl")]
use openssl::bn::BigNum;

use crypto_primes::{
    hazmat::{
        lucas_test, random_odd_uint, AStarBase, BruteForceBase, LucasCheck, MillerRabin,
        SelfridgeBase, Sieve,
    },
    is_prime_with_rng, is_safe_prime_with_rng, prime_with_rng, safe_prime_with_rng,
};

fn make_rng() -> ChaCha8Rng {
    ChaCha8Rng::from_seed(*b"01234567890123456789012345678901")
}

fn make_sieve<const L: usize>(rng: &mut impl CryptoRngCore) -> Sieve<L> {
    let start: Uint<L> = random_odd_uint(rng, Uint::<L>::BITS);
    Sieve::new(&start, Uint::<L>::BITS, false)
}

fn make_presieved_num<const L: usize>(rng: &mut impl CryptoRngCore) -> Uint<L> {
    let mut sieve = make_sieve(rng);
    sieve.next().unwrap()
}

fn bench_sieve(c: &mut Criterion) {
    let mut group = c.benchmark_group("Sieve");

    group.bench_function("(U128) random start", |b| {
        b.iter(|| random_odd_uint::<{ nlimbs!(128) }>(&mut OsRng, 128))
    });

    group.bench_function("(U128) creation", |b| {
        b.iter_batched(
            || random_odd_uint::<{ nlimbs!(128) }>(&mut OsRng, 128),
            |start| Sieve::new(&start, 128, false),
            BatchSize::SmallInput,
        )
    });

    // 5 is the average number of pre-sieved samples we need to take before we encounter a prime
    group.bench_function("(U128) 5 samples", |b| {
        b.iter_batched(
            || make_sieve::<{ nlimbs!(128) }>(&mut OsRng),
            |sieve| sieve.take(5).for_each(drop),
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U1024) random start", |b| {
        b.iter(|| random_odd_uint::<{ nlimbs!(1024) }>(&mut OsRng, 1024))
    });

    group.bench_function("(U1024) creation", |b| {
        b.iter_batched(
            || random_odd_uint::<{ nlimbs!(1024) }>(&mut OsRng, 1024),
            |start| Sieve::new(&start, 1024, false),
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U1024) 5 samples", |b| {
        b.iter_batched(
            || make_sieve::<{ nlimbs!(1024) }>(&mut OsRng),
            |sieve| sieve.take(5).for_each(drop),
            BatchSize::SmallInput,
        )
    });

    group.finish()
}

fn bench_miller_rabin(c: &mut Criterion) {
    let mut group = c.benchmark_group("Miller-Rabin");

    group.bench_function("(U128) creation", |b| {
        b.iter_batched(
            || random_odd_uint::<{ nlimbs!(128) }>(&mut OsRng, 128),
            |start| MillerRabin::new(&start),
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U128) random base test (pre-sieved)", |b| {
        b.iter_batched(
            || MillerRabin::new(&make_presieved_num::<{ nlimbs!(128) }>(&mut OsRng)),
            |mr| mr.test_random_base(&mut OsRng),
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U1024) creation", |b| {
        b.iter_batched(
            || random_odd_uint::<{ nlimbs!(1024) }>(&mut OsRng, 1024),
            |start| MillerRabin::new(&start),
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U1024) random base test (pre-sieved)", |b| {
        b.iter_batched(
            || MillerRabin::new(&make_presieved_num::<{ nlimbs!(1024) }>(&mut OsRng)),
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
            |n| lucas_test(&n, SelfridgeBase, LucasCheck::Strong),
            BatchSize::SmallInput,
        )
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) Selfridge base, strong check (pre-sieved)", |b| {
        b.iter_batched(
            || make_presieved_num::<{ nlimbs!(1024) }>(&mut rng),
            |n| lucas_test(&n, SelfridgeBase, LucasCheck::Strong),
            BatchSize::SmallInput,
        )
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) A* base, Lucas-V check (pre-sieved)", |b| {
        b.iter_batched(
            || make_presieved_num::<{ nlimbs!(1024) }>(&mut rng),
            |n| lucas_test(&n, AStarBase, LucasCheck::LucasV),
            BatchSize::SmallInput,
        )
    });

    let mut rng = make_rng();
    group.bench_function(
        "(U1024) brute force base, almost extra strong (pre-sieved)",
        |b| {
            b.iter_batched(
                || make_presieved_num::<{ nlimbs!(1024) }>(&mut rng),
                |n| lucas_test(&n, BruteForceBase, LucasCheck::AlmostExtraStrong),
                BatchSize::SmallInput,
            )
        },
    );

    let mut rng = make_rng();
    group.bench_function("(U1024) brute force base, extra strong (pre-sieved)", |b| {
        b.iter_batched(
            || make_presieved_num::<{ nlimbs!(1024) }>(&mut rng),
            |n| lucas_test(&n, BruteForceBase, LucasCheck::ExtraStrong),
            BatchSize::SmallInput,
        )
    });

    // A number triggering a slow path through the Lucas test, invoking all the available checks:
    // - U_{d} checked, but is not 0
    // - V_{d * 2^t} checked for t == 0..s-1, but no V = 0 found
    // - s = 5, so the previous step has multiple checks
    // - Q != 1 (since we're using Selfridge base)
    let slow_path = U1024::from_be_hex(concat![
        "D1CB9F1B6F3414A4B40A7E51C53C6AE4689DFCDC49FF875E7066A229D704EA8E",
        "6B674231D8C5974001673C3CE7FF9D377C8564E5182165A23434BC7B7E6C0419",
        "FD25C9921B0E9C90AF2570DB0772E1A9C82ACABBC8FC0F0864CE8A12124FA29B",
        "7F870924041DFA13EE5F5541C1BF96CA679EFAE2C96F5F4E9DF6007185198F5F"
    ]);

    group.bench_function("(U1024) Selfridge base, strong check, slow path", |b| {
        b.iter(|| {
            lucas_test(&slow_path, SelfridgeBase, LucasCheck::Strong);
        })
    });

    group.finish();
}

fn bench_presets(c: &mut Criterion) {
    let mut group = c.benchmark_group("Presets");

    group.bench_function("(U128) Prime test", |b| {
        b.iter_batched(
            || random_odd_uint::<{ nlimbs!(128) }>(&mut OsRng, 128),
            |num| is_prime_with_rng(&mut OsRng, &num),
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U128) Safe prime test", |b| {
        b.iter_batched(
            || random_odd_uint::<{ nlimbs!(128) }>(&mut OsRng, 128),
            |num| is_safe_prime_with_rng(&mut OsRng, &num),
            BatchSize::SmallInput,
        )
    });

    let mut rng = make_rng();
    group.bench_function("(U128) Random prime", |b| {
        b.iter(|| prime_with_rng::<{ nlimbs!(128) }>(&mut rng, 128))
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) Random prime", |b| {
        b.iter(|| prime_with_rng::<{ nlimbs!(1024) }>(&mut rng, 1024))
    });

    let mut rng = make_rng();
    group.bench_function("(U128) Random safe prime", |b| {
        b.iter(|| safe_prime_with_rng::<{ nlimbs!(128) }>(&mut rng, 128))
    });

    group.sample_size(20);
    let mut rng = make_rng();
    group.bench_function("(U1024) Random safe prime", |b| {
        b.iter(|| safe_prime_with_rng::<{ nlimbs!(1024) }>(&mut rng, 1024))
    });

    group.finish();
}

#[cfg(feature = "tests-gmp")]
fn bench_gmp(c: &mut Criterion) {
    let mut group = c.benchmark_group("GMP");

    fn random<const L: usize>(rng: &mut impl CryptoRngCore) -> Integer {
        let num: Uint<L> = random_odd_uint(rng, Uint::<L>::BITS);
        Integer::from_digits(num.as_words(), Order::Lsf)
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

criterion_group!(
    benches,
    bench_sieve,
    bench_miller_rabin,
    bench_lucas,
    bench_presets,
    bench_gmp,
    bench_openssl,
);
criterion_main!(benches);
