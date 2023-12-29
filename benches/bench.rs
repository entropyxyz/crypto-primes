use core::num::NonZeroU32;

use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use crypto_bigint::{nlimbs, BoxedUint, Integer, Odd, RandomBits, Uint, U1024, U128, U256};
use rand_chacha::ChaCha8Rng;
use rand_core::{CryptoRngCore, OsRng, SeedableRng};

#[cfg(feature = "tests-gmp")]
use rug::{integer::Order, Integer as GmpInteger};

#[cfg(feature = "tests-openssl")]
use openssl::bn::BigNum;

use crypto_primes::{
    generate_prime_with_rng, generate_safe_prime_with_rng,
    hazmat::{
        lucas_test, random_odd_integer, AStarBase, BruteForceBase, LucasCheck, MillerRabin,
        SelfridgeBase, Sieve,
    },
    is_prime_with_rng, is_safe_prime_with_rng,
};

fn make_rng() -> ChaCha8Rng {
    ChaCha8Rng::from_seed(*b"01234567890123456789012345678901")
}

fn random_odd_uint<T: RandomBits + Integer>(
    rng: &mut impl CryptoRngCore,
    bit_length: u32,
    bits_precision: u32,
) -> Odd<T> {
    random_odd_integer::<T>(rng, NonZeroU32::new(bit_length).unwrap(), bits_precision)
}

fn make_sieve<const L: usize>(rng: &mut impl CryptoRngCore) -> Sieve<Uint<L>> {
    let start = random_odd_uint::<Uint<L>>(rng, Uint::<L>::BITS, Uint::<L>::BITS);
    Sieve::new(&start, NonZeroU32::new(Uint::<L>::BITS).unwrap(), false)
}

fn make_presieved_num<const L: usize>(rng: &mut impl CryptoRngCore) -> Odd<Uint<L>> {
    let mut sieve = make_sieve(rng);
    Odd::new(sieve.next().unwrap()).unwrap()
}

fn bench_sieve(c: &mut Criterion) {
    let mut group = c.benchmark_group("Sieve");

    group.bench_function("(U128) random start", |b| {
        b.iter(|| random_odd_uint::<U128>(&mut OsRng, 128, 128))
    });

    group.bench_function("(U128) creation", |b| {
        b.iter_batched(
            || random_odd_uint::<U128>(&mut OsRng, 128, 128),
            |start| Sieve::new(start.as_ref(), NonZeroU32::new(128).unwrap(), false),
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
        b.iter(|| random_odd_uint::<U1024>(&mut OsRng, 1024, 1024))
    });

    group.bench_function("(U1024) creation", |b| {
        b.iter_batched(
            || random_odd_uint::<U1024>(&mut OsRng, 1024, 1024),
            |start| Sieve::new(start.as_ref(), NonZeroU32::new(1024).unwrap(), false),
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
            || random_odd_uint::<U128>(&mut OsRng, 128, 128),
            |n| MillerRabin::new(&n),
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
            || random_odd_uint::<U1024>(&mut OsRng, 1024, 1024),
            |n| MillerRabin::new(&n),
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
    let slow_path = Odd::new(U1024::from_be_hex(concat![
        "D1CB9F1B6F3414A4B40A7E51C53C6AE4689DFCDC49FF875E7066A229D704EA8E",
        "6B674231D8C5974001673C3CE7FF9D377C8564E5182165A23434BC7B7E6C0419",
        "FD25C9921B0E9C90AF2570DB0772E1A9C82ACABBC8FC0F0864CE8A12124FA29B",
        "7F870924041DFA13EE5F5541C1BF96CA679EFAE2C96F5F4E9DF6007185198F5F"
    ]))
    .unwrap();

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
            || random_odd_uint::<U128>(&mut OsRng, 128, 128),
            |num| is_prime_with_rng(&mut OsRng, num.as_ref()),
            BatchSize::SmallInput,
        )
    });

    group.bench_function("(U128) Safe prime test", |b| {
        b.iter_batched(
            || random_odd_uint::<U128>(&mut OsRng, 128, 128),
            |num| is_safe_prime_with_rng(&mut OsRng, num.as_ref()),
            BatchSize::SmallInput,
        )
    });

    let mut rng = make_rng();
    group.bench_function("(U128) Random prime", |b| {
        b.iter(|| generate_prime_with_rng::<U128>(&mut rng, 128, 128))
    });

    let mut rng = make_rng();
    group.bench_function("(U1024) Random prime", |b| {
        b.iter(|| generate_prime_with_rng::<U1024>(&mut rng, 1024, 1024))
    });

    let mut rng = make_rng();
    group.bench_function("(U128) Random safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<U128>(&mut rng, 128, 128))
    });

    group.sample_size(20);
    let mut rng = make_rng();
    group.bench_function("(U1024) Random safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<U1024>(&mut rng, 1024, 1024))
    });

    let mut rng = make_rng();
    group.bench_function("(Boxed128) Random safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<BoxedUint>(&mut rng, 128, 128))
    });

    group.sample_size(20);
    let mut rng = make_rng();
    group.bench_function("(Boxed1024) Random safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<BoxedUint>(&mut rng, 1024, 1024))
    });

    group.finish();

    // A separate group for bounded tests, to make it easier to run them separately.
    let mut group = c.benchmark_group("Presets (bounded)");

    let mut rng = make_rng();
    group.bench_function("(U128) Random safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<U128>(&mut rng, 128, 128))
    });

    // The performance should scale with the prime size, not with the Uint size.
    // So we should strive for this test's result to be as close as possible
    // to that of the previous one and as far away as possible from the next one.
    group.bench_function("(U256) Random 128 bit safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<U256>(&mut rng, 128, 256))
    });

    // The upper bound for the previous test.
    group.bench_function("(U256) Random 256 bit safe prime", |b| {
        b.iter(|| generate_safe_prime_with_rng::<U256>(&mut rng, 256, 256))
    });

    group.finish();
}

#[cfg(feature = "tests-gmp")]
fn bench_gmp(c: &mut Criterion) {
    let mut group = c.benchmark_group("GMP");

    fn random<const L: usize>(rng: &mut impl CryptoRngCore) -> GmpInteger {
        let num = random_odd_uint::<Uint<L>>(rng, Uint::<L>::BITS, Uint::<L>::BITS).get();
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
