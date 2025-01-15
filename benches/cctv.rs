use std::io::BufRead;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use crypto_bigint::U1024;
use rand_chacha::ChaCha8Rng;

use crypto_primes::is_prime_with_rng;
use rand_core::SeedableRng;

/// CCTV stands for Community Cryptography Test Vectors[1]. This benchmark uses the
/// "rsa.bench.2048.txt" test vector, which is a file of 708 1024-bit long candidates for prime
/// testing. The series of candidates in the test vecotr is an average representative sequence of
/// candidates that can be tested across different implementations. There are two primes in the
/// file, the first at line 354 and the other on the last line. Unless there's a bug, the second
/// half of the vector is not traversed in this benchmark.
///
/// [1]: https://github.com/C2SP/CCTV
fn bench_cctv(c: &mut Criterion) {
    let mut group = c.benchmark_group("CCTV RSA 1024-bit candidates");
    group.sample_size(10);
    let mut rng = ChaCha8Rng::from_seed([123; 32]);
    let candidates: Vec<U1024> = std::fs::read("./benches/rsa.bench.2048.txt")
        .expect("file present")
        .lines()
        .map(|candidate_hex| U1024::from_be_hex(&candidate_hex.unwrap()))
        .collect();

    assert!(
        is_prime_with_rng(&mut rng, &candidates[353]),
        "Line 354 is a prime. This is a bug."
    );
    assert!(
        is_prime_with_rng(&mut rng, &candidates[707]),
        "Line 708 is a prime. This is a bug."
    );

    group.bench_function("all", |b| {
        b.iter(|| {
            for candidate in &candidates {
                black_box(is_prime_with_rng(&mut rng, candidate));
            }
        });
    });
}

criterion_group!(benches, bench_cctv,);
criterion_main!(benches);
