use std::io::BufRead;

use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use crypto_bigint::U1024;
use rand_chacha::ChaCha8Rng;

use crypto_primes::{is_prime_with_rng, sieve_and_find, SieveFactory};
use rand_core::SeedableRng;

struct HardcodedSieve {
    candidates: Vec<String>,
    idx: usize,
}

impl HardcodedSieve {
    fn new() -> Self {
        Self {
            candidates: std::fs::read("./benches/rsa.bench.2048.txt")
                .expect("file present")
                .lines()
                .map(|l| l.unwrap().to_owned())
                .collect(),
            idx: 0,
        }
    }
}

impl Iterator for HardcodedSieve {
    type Item = U1024;
    fn next(&mut self) -> Option<Self::Item> {
        self.idx += 1;
        self.candidates.get(self.idx).map(|str| U1024::from_be_hex(str))
    }
}

impl SieveFactory for HardcodedSieve {
    type Item = U1024;
    type Sieve = Self;
    fn make_sieve(
        &mut self,
        _rng: &mut impl rand_core::CryptoRngCore,
        _previous_sieve: Option<&Self::Sieve>,
    ) -> Option<Self::Sieve> {
        Some(Self {
            candidates: self.candidates.clone(),
            idx: 0,
        })
    }
}

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
    group
        .sample_size(10)
        .measurement_time(std::time::Duration::from_secs(10));
    let mut rng = ChaCha8Rng::from_seed([123; 32]);
    let candidates: Vec<U1024> = std::fs::read("./benches/rsa.bench.2048.txt")
        .expect("file present")
        .lines()
        .map(|candidate_hex| U1024::from_be_hex(&candidate_hex.unwrap()))
        .collect();

    group.bench_function("manual", |b| {
        b.iter(|| {
            for candidate in &candidates {
                if is_prime_with_rng(&mut rng, candidate) {
                    break;
                }
            }
        });
    });

    group.bench_function("factory API", |b| {
        b.iter_batched(
            HardcodedSieve::new,
            |sieve| sieve_and_find(&mut rng, sieve, is_prime_with_rng),
            BatchSize::SmallInput,
        );
    });
}

criterion_group!(benches, bench_cctv,);
criterion_main!(benches);
