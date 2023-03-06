# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [0.2.0] - 2023-03-06

### Changed

- Bumped MSRV to 1.65 (`crypto-bigint` requirement). ([#18])
- Bumped `crypto-bigint` to 0.5 with the corresponding changes in the API. ([#13], [#18])
- `Sieve::new()` now does not panic; if the parameters imply an empty output, that's what the iterator generates. ([#14])


### Added

- `RandomPrimeWithRng` trait for use in generic code. ([#12])


### Fixed

- Added a `gcd(Q, n) == 1` check to the Lucas test. ([#11])
- Added a simple 3 mod 4 check for safe primes, speeding up their generation. ([#14])
- Fixed multiple corner cases in `Sieve` for small ranges. ([#14])


[#11]: https://github.com/nucypher/rust-umbral/pull/11
[#12]: https://github.com/nucypher/rust-umbral/pull/12
[#13]: https://github.com/nucypher/rust-umbral/pull/13
[#14]: https://github.com/nucypher/rust-umbral/pull/14
[#18]: https://github.com/nucypher/rust-umbral/pull/18


## [0.1.0] - 2023-01-20

Initial release.


[Unreleased]: https://github.com/entropyxyz/crypto-primes/compare/v0.2.0...HEAD
[0.1.0]: https://github.com/nucypher/rust-umbral/releases/tag/v0.1.0
[0.2.0]: https://github.com/nucypher/rust-umbral/releases/tag/v0.2.0
