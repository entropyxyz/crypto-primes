# Prime number tools for `crypto-bigint`

[![crate][crate-image]][crate-link]
[![Docs][docs-image]][docs-link]
![License][license-image]
[![Build Status][build-image]][build-link]
[![Coverage][coverage-image]][coverage-link]

This library implements prime number generation and primality checking for [`crypto-bigint`](https://crates.io/crates/crypto-bigint) integers.
In particular:

- Generating random primes and safe primes of given bit size;
- Sieving iterator;
- Miller-Rabin test;
- Strong and extra strong Lucas tests, and Lucas-V test.


The library is no-std compatible and contains no unsafe code.

Most users will be using the small set of functions exported from the top level, providing "pre-packaged" prime finding functionality with sane defaults.

## Example

Find a 196 bit prime returned in a 256-bit long `crypto_bigint::U256`:

```rust
use crypto_bigint::U256;
let prime = crypto_primes::generate_prime::<U256>(196);
assert!(crypto_primes::is_prime(&prime));
```

Find a 64 bit safe prime returned in a `crypto_bigint::U1024`:

```rust
use crypto_bigint::U1024;
let prime = crypto_primes::generate_safe_prime::<U1024>(64);
assert!(crypto_primes::is_safe_prime(&prime));
```

## Advanced

Advanced users can use the [`hazmat`][hazmat-lnk] module in the library to build a custom prime finding solution that best fit their needs, e.g. by picking different Lucas bases or running Miller-Rabin tests with particular bases.

## Features

The following features are available:

- `default-rng`: Use the OS default CSPRNG, `OsRng`. Enabled by default.
- `multicore`: Enables additional parallel prime finding functions. Disabled by default.


[crate-image]: https://img.shields.io/crates/v/crypto-primes.svg
[crate-link]: https://crates.io/crates/crypto-primes
[docs-image]: https://docs.rs/crypto-primes/badge.svg
[docs-link]: https://docs.rs/crypto-primes/
[license-image]: https://img.shields.io/crates/l/crypto-primes
[build-image]: https://github.com/entropyxyz/crypto-primes/workflows/crypto-primes/badge.svg?branch=master&event=push
[build-link]: https://github.com/entropyxyz/crypto-primes/actions?query=workflow%3Acrypto-primes
[coverage-image]: https://codecov.io/gh/entropyxyz/crypto-primes/branch/master/graph/badge.svg
[coverage-link]: https://codecov.io/gh/entropyxyz/crypto-primes
[hazmat-lnk]: https://docs.rs/crypto-primes/latest/crypto_primes/hazmat
