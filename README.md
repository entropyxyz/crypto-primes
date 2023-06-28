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

See the documentation for the specific tests for more information and references.


[crate-image]: https://img.shields.io/crates/v/crypto-primes.svg
[crate-link]: https://crates.io/crates/crypto-primes
[docs-image]: https://docs.rs/crypto-primes/badge.svg
[docs-link]: https://docs.rs/crypto-primes/
[license-image]: https://img.shields.io/crates/l/crypto-primes
[build-image]: https://github.com/entropyxyz/crypto-primes/workflows/crypto-primes/badge.svg?branch=master&event=push
[build-link]: https://github.com/entropyxyz/crypto-primes/actions?query=workflow%3Acrypto-primes
[coverage-image]: https://codecov.io/gh/entropyxyz/crypto-primes/branch/master/graph/badge.svg
[coverage-link]: https://codecov.io/gh/entropyxyz/crypto-primes
