use rand_core::OsRng;

use criterion::{
    criterion_group, criterion_main, measurement::Measurement, BenchmarkGroup, Criterion,
};
use crypto_bigint::U1024;

use crypto_primes::hazmat::{
    is_lucas_prime, random_odd_uint, BruteForceBase, LucasCheck, MillerRabin, SelfridgeBase, Sieve,
};

fn bench_sieve<'a, M: Measurement>(group: &mut BenchmarkGroup<'a, M>) {
    let start: U1024 = random_odd_uint(&mut OsRng, 1024);
    group.bench_function("(U1024) Sieve, 1000 samples", |b| {
        b.iter(|| Sieve::new(&start, 1024).take(1000).for_each(drop))
    });
}

fn bench_miller_rabin<'a, M: Measurement>(group: &mut BenchmarkGroup<'a, M>) {
    let start: U1024 = random_odd_uint(&mut OsRng, 1024);
    group.bench_function("(U1024) Miller-Rabin creation", |b| {
        b.iter(|| {
            MillerRabin::new(&start);
        })
    });

    let start: U1024 = random_odd_uint(&mut OsRng, 1024);
    let mut sieve = Sieve::new(&start, 1024);
    group.bench_function(
        "(U1024) Sieve + Miller-Rabin creation + random base check",
        |b| {
            b.iter(|| {
                let mr = MillerRabin::new(&sieve.next().unwrap());
                mr.check_random_base(&mut OsRng);
            })
        },
    );
}

fn bench_lucas<'a, M: Measurement>(group: &mut BenchmarkGroup<'a, M>) {
    let start: U1024 = random_odd_uint(&mut OsRng, 1024);
    let mut sieve = Sieve::new(&start, 1024);
    group.bench_function("(U1024) Sieve + Lucas test (Selfridge base)", |b| {
        b.iter(|| {
            is_lucas_prime(&sieve.next().unwrap(), SelfridgeBase, LucasCheck::Strong);
        })
    });

    group.bench_function(
        "(U1024) Sieve + Lucas test (brute force base, almost extra strong)",
        |b| {
            b.iter(|| {
                is_lucas_prime(
                    &sieve.next().unwrap(),
                    BruteForceBase,
                    LucasCheck::AlmostExtraStrong,
                );
            })
        },
    );

    group.bench_function(
        "(U1024) Sieve + Lucas test (brute force base, extra strong)",
        |b| {
            b.iter(|| {
                is_lucas_prime(
                    &sieve.next().unwrap(),
                    BruteForceBase,
                    LucasCheck::ExtraStrong,
                );
            })
        },
    );

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

    group.bench_function("(U1024) Lucas test (Selfridge base), slow path", |b| {
        b.iter(|| {
            is_lucas_prime(&slow_path, SelfridgeBase, LucasCheck::Strong);
        })
    });
}

fn bench_primality_tests(c: &mut Criterion) {
    let mut group = c.benchmark_group("primality tests");
    bench_sieve(&mut group);
    bench_miller_rabin(&mut group);
    bench_lucas(&mut group);
    group.finish();
}

criterion_group!(benches, bench_primality_tests);
criterion_main!(benches);
