use rand_core::OsRng;

use criterion::{
    criterion_group, criterion_main, measurement::Measurement, BenchmarkGroup, Criterion,
};

use crypto_primes::hazmat::{lucas, lucas2, lucas_new, random_odd_uint, MillerRabin, Sieve};

fn bench_sieve<'a, M: Measurement>(group: &mut BenchmarkGroup<'a, M>) {
    let start = random_odd_uint::<16, _>(&mut OsRng, 1024);
    group.bench_function("(U1024) Sieve, 1000 samples", |b| {
        b.iter(|| Sieve::new(&start, 1024).take(1000).for_each(drop))
    });
}

fn bench_miller_rabin<'a, M: Measurement>(group: &mut BenchmarkGroup<'a, M>) {
    let start = random_odd_uint::<16, _>(&mut OsRng, 1024);
    group.bench_function("(U1024) Miller-Rabin creation", |b| {
        b.iter(|| {
            MillerRabin::new(&start);
        })
    });

    let start = random_odd_uint::<16, _>(&mut OsRng, 1024);
    let mut sieve = Sieve::new(&start, 1024);
    group.bench_function(
        "(U1024) Sieve + Miller-Rabin creation + random basis check",
        |b| {
            b.iter(|| {
                let mr = MillerRabin::new(&sieve.next().unwrap());
                mr.check_random_basis(&mut OsRng);
            })
        },
    );
}

fn bench_lucas<'a, M: Measurement>(group: &mut BenchmarkGroup<'a, M>) {
    let start = random_odd_uint::<16, _>(&mut OsRng, 1024);
    let mut sieve = Sieve::new(&start, 1024);
    group.bench_function("(U1024) Sieve + Lucas test (Selfridge base)", |b| {
        b.iter(|| {
            lucas_new::is_strong_lucas_prime::<lucas_new::SelfridgeBase, 16>(
                &sieve.next().unwrap(),
                true,
            );
        })
    });

    group.bench_function(
        "(U1024) Sieve + Lucas test (brute force base, almost extra strong)",
        |b| {
            b.iter(|| {
                lucas_new::is_strong_lucas_prime::<lucas_new::BruteForceBase, 16>(
                    &sieve.next().unwrap(),
                    false,
                );
            })
        },
    );

    group.bench_function(
        "(U1024) Sieve + Lucas test (brute force base, extra strong)",
        |b| {
            b.iter(|| {
                lucas_new::is_strong_lucas_prime::<lucas_new::BruteForceBase, 16>(
                    &sieve.next().unwrap(),
                    true,
                );
            })
        },
    );
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
