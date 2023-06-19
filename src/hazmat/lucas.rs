//! Lucas primality test.
use crypto_bigint::{
    modular::runtime_mod::{DynResidue, DynResidueParams},
    Integer, Limb, Uint, Word,
};

use super::{
    gcd::gcd,
    jacobi::{jacobi_symbol, JacobiSymbol},
    Primality,
};

/// The maximum number of attempts to find `D` such that `(D/n) == -1`.
// This is widely believed to be impossible.
// So if we exceed it, we will panic reporting the value of `n`.
const MAX_ATTEMPTS: u32 = 10000;

/// The number of attempts to find `D` such that `(D/n) == -1`
/// before checking that `n` is a square (in which case such `D` does not exist).
// This check is relatively expensive compared to calculating the Jacobi symbol
// (~30x for 1024-bit numbers, ~100x for 2048 bit).
// On the other hand, if `n` is a non-square we expect to find a `D`
// in just a few attempts on average (an estimate for the Selfridge method
// can be found in [^Baillie1980], section 7; for the brute force method
// it seems to be about the same).
const ATTEMPTS_BEFORE_SQRT: u32 = 30;

/// A method for selecting the base `(P, Q)` for the Lucas primality test.
pub trait LucasBase {
    /// Given an odd integer, returns `Ok((P, Q))` on success,
    /// or `Err(Primality)` if the primality for the given integer was discovered
    /// during the search for a base.
    fn generate<const L: usize>(&self, n: &Uint<L>) -> Result<(u32, i32), Primality>;
}

/// "Method A" for selecting the base given in Baillie & Wagstaff[^Baillie1980],
/// attributed to Selfridge.
///
/// Try `D = 1 - 4Q = 5, -7, 9, -11, 13, ...` until `Jacobi(D, n) = -1`.
/// Return `P = 1, Q = (1 - D) / 4)`.
///
/// [^Baillie1980]:
///   R. Baillie, S. S. Wagstaff, "Lucas pseudoprimes",
///   Math. Comp. 35 1391-1417 (1980),
///   DOI: [10.2307/2006406](https://dx.doi.org/10.2307/2006406),
///   <http://mpqs.free.fr/LucasPseudoprimes.pdf>
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct SelfridgeBase;

impl LucasBase for SelfridgeBase {
    fn generate<const L: usize>(&self, n: &Uint<L>) -> Result<(u32, i32), Primality> {
        let mut d = 5_i32;
        let n_is_small = n.bits_vartime() < (Limb::BITS - 1);
        // Can unwrap here since it won't overflow after `&`
        let small_n: u32 = (n.as_words()[0] & Word::from(u32::MAX)).try_into().unwrap();
        let mut attempts = 0;
        loop {
            if attempts >= MAX_ATTEMPTS {
                panic!("internal error: cannot find (D/n) = -1 for {:?}", n)
            }

            if attempts >= ATTEMPTS_BEFORE_SQRT {
                let sqrt_n = n.sqrt();
                if &sqrt_n.wrapping_mul(&sqrt_n) == n {
                    return Err(Primality::Composite);
                }
            }

            let j = jacobi_symbol(d, n);

            if j == JacobiSymbol::MinusOne {
                break;
            }
            if j == JacobiSymbol::Zero {
                // Modification of Method A by Baillie, in an example to OEIS:A217120
                // (https://oeis.org/A217120/a217120_1.txt):
                // If `d == (+,-)n`, (e.g., `n` = 5 or 11) try the next `d` instead of quitting;
                // this small modification of Selfridge's method A
                // enables 5 and 11 to be classified as Lucas probable primes.
                // Otherwise GCD(D, n) > 1, and therefore n is not prime.
                if !(n_is_small && small_n == d.abs_diff(0)) {
                    return Err(Primality::Composite);
                }
            }

            attempts += 1;
            d = -d;
            d += d.signum() * 2;
        }

        // No remainder by construction of `d`.
        let q = (1 - d) / 4;

        Ok((1, q))
    }
}

/// "Method A*" for selecting the base given in Baillie & Wagstaff[^Baillie1980].
///
/// Same as [`SelfridgeBase`], but returns `(P = 5, Q = 5)` if the Selfridge base set `Q = -1`.
///
/// [^Baillie1980]:
///   R. Baillie, S. S. Wagstaff, "Lucas pseudoprimes",
///   Math. Comp. 35 1391-1417 (1980),
///   DOI: [10.2307/2006406](https://dx.doi.org/10.2307/2006406),
///   <http://mpqs.free.fr/LucasPseudoprimes.pdf>
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct AStarBase;

impl LucasBase for AStarBase {
    fn generate<const L: usize>(&self, n: &Uint<L>) -> Result<(u32, i32), Primality> {
        SelfridgeBase
            .generate(n)
            .map(|(p, q)| if q == -1 { (5, 5) } else { (p, q) })
    }
}

/// "Method C" for selecting the base given by Baillie[^Baillie].
///
/// Try `P = 3, 4, 5, ...` until `Jacobi(D, n) = -1`, where `D = P^2 - 4Q`.
/// Returns the found `P`, and `Q = 1`.
///
/// [^Baillie]:
///   R. Baillie, Mathematica code for extra strong Lucas pseudoprimes,
///   <https://oeis.org/A217719/a217719.txt>
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct BruteForceBase;

impl LucasBase for BruteForceBase {
    fn generate<const L: usize>(&self, n: &Uint<L>) -> Result<(u32, i32), Primality> {
        let mut p = 3_u32;
        let mut attempts = 0;

        loop {
            if attempts >= MAX_ATTEMPTS {
                panic!("internal error: cannot find (D/n) = -1 for {:?}", n)
            }

            if attempts >= ATTEMPTS_BEFORE_SQRT {
                let sqrt_n = n.sqrt();
                if &sqrt_n.wrapping_mul(&sqrt_n) == n {
                    return Err(Primality::Composite);
                }
            }

            // Can unwrap here since `p` is always small (see the condition above).
            let j = jacobi_symbol((p * p - 4).try_into().unwrap(), n);

            if j == JacobiSymbol::MinusOne {
                break;
            }
            if j == JacobiSymbol::Zero {
                // D = P^2 - 4 = (P - 2)(P + 2).
                // If (D/n) == 0 then D shares a prime factor with n.
                // Since the loop proceeds in increasing P and starts with P - 2 == 1,
                // the shared prime factor must be P + 2.
                // If P + 2 == n, then n is prime; otherwise P + 2 is a proper factor of n.
                let primality = if n == &Uint::<L>::from(p + 2) {
                    Primality::Prime
                } else {
                    Primality::Composite
                };
                return Err(primality);
            }

            attempts += 1;
            p += 1;
        }

        Ok((p, 1))
    }
}

/// For the given odd `n`, finds `s` and odd `d` such that `n + 1 == 2^s * d`.
fn decompose<const L: usize>(n: &Uint<L>) -> (u32, Uint<L>) {
    debug_assert!(bool::from(n.is_odd()));

    // Need to be careful here since `n + 1` can overflow.
    // Instead of adding 1 and counting trailing 0s, we count trailing ones on the original `n`.

    let mut n = *n;
    let mut s = 0;

    while n.is_odd().into() {
        n >>= 1;
        s += 1;
    }

    // This won't overflow since the original `n` was odd, so we right-shifted at least once.
    (s, n.wrapping_add(&Uint::<L>::ONE))
}

/// The checks to perform in the Lucas test.
///
/// Given the Lucas sequence built from some base `(P, Q)` (see [`LucasBase`])
/// up to the elements `V(d)`, `U(d)`, where `d * 2^s == n - (D/n)`, `d` odd, and `D = P^2 - 4Q`,
/// the checks are defined as follows:
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum LucasCheck {
    /// Introduced by Baillie & Wagstaff[^Baillie1980].
    /// If either of the following is true:
    /// - any of `V(d*2^r) == 0` for `0 <= r < s`,
    /// - `U(d) == 0`,
    /// report the number as prime.
    ///
    /// If the base is [`SelfridgeBase`], known false positives constitute OEIS:A217255[^A217255].
    ///
    /// [^Baillie1980]:
    ///   R. Baillie, S. S. Wagstaff, "Lucas pseudoprimes",
    ///   Math. Comp. 35 1391-1417 (1980),
    ///   DOI: [10.2307/2006406](https://dx.doi.org/10.2307/2006406),
    ///   <http://mpqs.free.fr/LucasPseudoprimes.pdf>
    ///
    /// [^A217255]: <https://oeis.org/A217255>
    Strong,

    /// A [`LucasCheck::ExtraStrong`] without checking for `U(d)`.
    /// That is, if either of the following is true:
    /// - any of `V(d*2^r) == 0` for `0 <= r < s`,
    /// - `V(d) == ±2`,
    /// report the number as prime.
    ///
    /// Note: the second condition is only checked if `Q == 1`,
    /// otherwise it is considered to be true.
    ///
    /// If the base is [`BruteForceBase`], some known false positives
    /// are listed by Jacobsen[^Jacobsen].
    ///
    /// Note: this option is intended for testing against known pseudoprimes;
    /// do not use unless you know what you are doing.
    ///
    /// [^Jacobsen]:
    ///   D. Jacobsen, "Pseudoprime Statistics, Tables, and Data",
    ///   <http://ntheory.org/pseudoprimes.html>
    AlmostExtraStrong,

    /// Introduced by Mo[^Mo1993], and also described by Grantham[^Grantham2001].
    /// If either of the following is true:
    /// - any of `V(d*2^r) == 0` for `0 <= r < s`,
    /// - `U(d) == 0` and `V(d) == ±2`,
    /// report the number as prime.
    ///
    /// Note that this check only differs from [`LucasCheck::Strong`] if `Q == 1`.
    ///
    /// If the base is [`BruteForceBase`], known false positives constitute OEIS:A217719[^A217719].
    ///
    /// [^Mo1993]:
    ///   Zhaiyu Mo, "Diophantine equations, Lucas sequences and pseudoprimes",
    ///   graduate thesis, University of Calgary, Calgary, AB (1993)
    ///   DOI: [10.11575/PRISM/10820](https://dx.doi.org/10.11575/PRISM/10820)
    ///
    /// [^Grantham2001]:
    ///   J. Grantham, "Frobenius pseudoprimes",
    ///   Math. Comp. 70 873-891 (2001),
    ///   DOI: [10.1090/S0025-5718-00-01197-2](https://dx.doi.org/10.1090/S0025-5718-00-01197-2)
    ///
    /// [^A217719]: <https://oeis.org/A217719>
    ExtraStrong,

    /// Introduced by Baillie et al[^Baillie2021].
    /// If `V(n+1) == 2 Q`, report the number as prime.
    ///
    /// If the base is [`AStarBase`], some known false positives
    /// are provided by Baillie et al[^Baillie2021].
    ///
    /// [^Baillie2021]: R. Baillie, A. Fiori, S. S. Wagstaff,
    ///   "Strengthening the Baillie-PSW primality test",
    ///   Math. Comp. 90 1931-1955 (2021),
    ///   DOI: [10.1090/mcom/3616](https://doi.org/10.1090/mcom/3616)
    LucasV,
}

/// Performs the primality test based on Lucas sequence.
/// See [`LucasCheck`] for possible checks, and the implementors of [`LucasBase`]
/// for the corresponding bases.
pub fn lucas_test<const L: usize>(
    candidate: &Uint<L>,
    base: impl LucasBase,
    check: LucasCheck,
) -> Primality {
    // The comments in this function use references in `LucasCheck`, plus this one:
    //
    // [^Crandall2005]:
    //   R. Crandall, C. Pomerance, "Prime numbers: a computational perspective",
    //   2nd ed., Springer (2005) (ISBN: 0-387-25282-7, 978-0387-25282-7)

    if candidate.is_even().into() {
        return Primality::Composite;
    }

    // Find the base for the Lucas sequence.
    let (p, q) = match base.generate(candidate) {
        Ok((p, q)) => (p, q),
        Err(primality) => return primality,
    };
    let discriminant = (p * p) as i32 - 4 * q;

    // If either is true, it allows us to optimize certain parts of the calculations.
    let p_is_one = p == 1;
    let q_is_one = q == 1;

    // See the references for the specific checks in the docstrings for [`LucasCheck`].

    // All of the definitions require gcd(n, 2QD) == 1.
    // We know gcd(n, D) = 1 by construction of the base (D is chosen such that (D/n) != 0).
    // We know gcd(n, 2) = 1 (we checked for it earlier).
    // In practice, gcd(n, Q) = 1 is always true, because the Lucas test is preceded by a sieve,
    // and since `Q` is always small, division by it would have been already checked.
    // But in order to avoid an implicit assumption that a sieve has been run,
    // we check that gcd(n, Q) = 1 anyway - again, since `Q` is small,
    // it does not noticeably affect the performance.
    let abs_q = q.abs_diff(0);
    if abs_q != 1 && gcd(candidate, abs_q) != 1 && candidate > &Uint::<L>::from(abs_q) {
        return Primality::Composite;
    }

    // Find d and s, such that d is odd and d * 2^s = (n - (D/n)).
    let (s, d) = decompose(candidate);

    // Some constants in Montgomery form

    let params = DynResidueParams::<L>::new(candidate);

    let zero = DynResidue::<L>::zero(params);
    let one = DynResidue::<L>::one(params);
    let two = one + one;
    let minus_two = -two;

    // Convert Q to Montgomery form

    let q = if q_is_one {
        one
    } else {
        let abs_q = DynResidue::<L>::new(&Uint::<L>::from(q.abs_diff(0)), params);
        if q < 0 {
            -abs_q
        } else {
            abs_q
        }
    };

    // Convert P to Montgomery form

    let p = if p_is_one {
        one
    } else {
        DynResidue::<L>::new(&Uint::<L>::from(p), params)
    };

    // Compute d-th element of Lucas sequence (U_d(P, Q), V_d(P, Q)), where:
    //
    // V_0 = 2
    // U_0 = 0
    //
    // U_{2k} = U_k V_k
    // V_{2k} = V_k^2 - 2 Q^k
    //
    // U_{k+1} = (P U_k + V_k) / 2
    // V_{k+1} = (D U_k + P V_k) / 2
    //
    // (The propagation method is due to [^Baillie2021], Eqs. 13, 14, 16, 17)
    // We can therefore start with k=0 and build up to k=d in log2(d) steps.

    // Starting with k = 0
    let mut vk = two; // keeps V_k
    let mut uk = DynResidue::<L>::zero(params); // keeps U_k
    let mut qk = one; // keeps Q^k

    // D in Montgomery representation - note that it can be negative.
    let abs_d = DynResidue::<L>::new(&Uint::<L>::from(discriminant.abs_diff(0)), params);
    let d_m = if discriminant < 0 { -abs_d } else { abs_d };

    for i in (0..d.bits_vartime()).rev() {
        // k' = k * 2

        let u_2k = uk * vk;
        let v_2k = vk.square() - (qk + qk);
        let q_2k = qk.square();

        uk = u_2k;
        vk = v_2k;
        qk = q_2k;

        if d.bit_vartime(i) {
            // k' = k + 1

            let (p_uk, p_vk) = if p_is_one { (uk, vk) } else { (p * uk, p * vk) };

            let u_k1 = (p_uk + vk).div_by_2();
            let v_k1 = (d_m * uk + p_vk).div_by_2();
            let q_k1 = qk * q;

            uk = u_k1;
            vk = v_k1;
            qk = q_k1;
        }
    }

    // Now k=d, so vk = V_d and uk = U_d.

    // Check for the first sufficient condition in various strong checks.

    if check == LucasCheck::Strong && uk == zero {
        // Strong check: `U_d == 0 mod n`.
        return Primality::ProbablyPrime;
    } else if check == LucasCheck::ExtraStrong || check == LucasCheck::AlmostExtraStrong {
        // Extra strong check (from [^Mo1993]): `V_d == ±2 mod n` and `U_d == 0 mod n`.
        //
        // Note that the first identity only applies if `Q = 1`, since it is a consequence
        // of a property of Lucas series: `V_k^2 - 4 Q^k = D U_k^2 mod n`.
        // If `Q = 1` we can easily decompose the left side of the equation
        // leading to the check above.
        //
        // If `Q != 1` we just consider it passed (we don't have a corresponding
        // pseudoprime list anyway).

        let vk_equals_two = !q_is_one || (vk == two || vk == minus_two);

        if check == LucasCheck::ExtraStrong && uk == zero && vk_equals_two {
            return Primality::ProbablyPrime;
        }

        // "Almost extra strong" check skips the `U_d` check.
        // Since we have `U_d` anyway, it does not improve performance,
        // so it is only here for testing purposes, since we have a corresponding pseudoprime list.
        if check == LucasCheck::AlmostExtraStrong && vk_equals_two {
            return Primality::ProbablyPrime;
        }
    }

    // Second sufficient condition requires further propagating `V_k` up to `V_{n+1}`.

    // Check if V_{2^t d} == 0 mod n for some 0 <= t < s.
    // (unless we're in Lucas-V mode, then we just propagate V_k)

    if check != LucasCheck::LucasV && vk == zero {
        return Primality::ProbablyPrime;
    }

    for _ in 1..s {
        // Optimization: V_k = ±2 is a fixed point for V_k' = V_k^2 - 2 Q^k with Q = 1,
        // so if V_k = ±2, we can stop: we will never find a future V_k == 0.
        if check != LucasCheck::LucasV && q_is_one && (vk == two || vk == minus_two) {
            return Primality::Composite;
        }

        // k' = 2k
        // V_{k'} = V_k^2 - 2 Q^k
        vk = vk * vk - qk - qk;

        if check != LucasCheck::LucasV && vk == zero {
            return Primality::ProbablyPrime;
        }

        if !q_is_one {
            qk = qk.square();
        }
    }

    if check == LucasCheck::LucasV {
        // At this point vk = V_{d * 2^(s-1)}.
        // Double the index again:
        vk = vk * vk - qk - qk; // now vk = V_{d * 2^s} = V_{n+1}

        // Lucas-V check[^Baillie2021]: if V_{n+1} == 2 Q, report `n` as prime.
        if vk == q + q {
            return Primality::ProbablyPrime;
        }
    }

    Primality::Composite
}

#[cfg(test)]
mod tests {

    use alloc::format;

    use crypto_bigint::{Uint, U128, U64};

    #[cfg(feature = "tests-exhaustive")]
    use num_prime::nt_funcs::is_prime64;

    use super::{
        decompose, lucas_test, AStarBase, BruteForceBase, LucasBase, LucasCheck, SelfridgeBase,
    };
    use crate::hazmat::{primes, pseudoprimes, Primality};

    #[test]
    fn bases_derived_traits() {
        assert_eq!(format!("{SelfridgeBase:?}"), "SelfridgeBase");
        assert_eq!(SelfridgeBase.clone(), SelfridgeBase);
        assert_eq!(format!("{AStarBase:?}"), "AStarBase");
        assert_eq!(AStarBase.clone(), AStarBase);
        assert_eq!(format!("{BruteForceBase:?}"), "BruteForceBase");
        assert_eq!(BruteForceBase.clone(), BruteForceBase);

        assert_eq!(format!("{:?}", LucasCheck::Strong), "Strong");
        assert_eq!(LucasCheck::Strong.clone(), LucasCheck::Strong);
    }

    #[test]
    fn base_for_square() {
        // We can't find a base with Jacobi symbol = -1 for a square,
        // check that it is handled properly.
        let num = U64::from(131u32).square();
        assert_eq!(SelfridgeBase.generate(&num), Err(Primality::Composite));
        assert_eq!(AStarBase.generate(&num), Err(Primality::Composite));
        assert_eq!(BruteForceBase.generate(&num), Err(Primality::Composite));
    }

    #[test]
    fn base_early_quit() {
        // 5 is flagged as prime at the base generation stage
        assert_eq!(
            BruteForceBase.generate(&U64::from(5u32)),
            Err(Primality::Prime)
        )
    }

    #[test]
    fn lucas_early_quit() {
        // If the number is even, no need to run the test.
        assert_eq!(
            lucas_test(&U64::from(6u32), SelfridgeBase, LucasCheck::Strong),
            Primality::Composite
        );
    }

    #[test]
    fn gcd_check() {
        // Test that `gcd(2QD, n) == 1` is checked after the base is found.
        // All the bases already produce D such that `gcd(D, n) == 1`.
        // We need a special test base generator to test the situation where
        // `gcd(n, Q) != 1` and `n > Q`, because it just doesn't seem to happen normally
        // for "production" bases.

        struct TestBase;

        impl LucasBase for TestBase {
            fn generate<const L: usize>(&self, _n: &Uint<L>) -> Result<(u32, i32), Primality> {
                Ok((5, 5))
            }
        }

        assert_eq!(
            lucas_test(&U64::from(15u32), TestBase, LucasCheck::Strong),
            Primality::Composite
        );
    }

    #[test]
    fn decomposition() {
        assert_eq!(decompose(&U128::MAX), (128, U128::ONE));
        assert_eq!(decompose(&U128::ONE), (1, U128::ONE));
        assert_eq!(decompose(&U128::from(7766015u32)), (15, U128::from(237u32)));
    }

    fn is_slpsp(num: u32) -> bool {
        pseudoprimes::STRONG_LUCAS.iter().any(|x| *x == num)
    }

    fn is_aeslpsp(num: u32) -> bool {
        pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS
            .iter()
            .any(|x| *x == num)
    }

    fn is_eslpsp(num: u32) -> bool {
        pseudoprimes::EXTRA_STRONG_LUCAS.iter().any(|x| *x == num)
    }

    fn is_vpsp(num: u32) -> bool {
        pseudoprimes::LUCAS_V.iter().any(|x| *x == num)
    }

    fn test_composites_selfridge(numbers: &[u32], expected_result: bool) {
        for num in numbers.iter() {
            let false_positive = is_slpsp(*num);
            let actual_expected_result = if false_positive {
                true
            } else {
                expected_result
            };

            // Test both single-limb and multi-limb, just in case.
            assert_eq!(
                lucas_test(&Uint::<1>::from(*num), SelfridgeBase, LucasCheck::Strong)
                    .is_probably_prime(),
                actual_expected_result
            );
            assert_eq!(
                lucas_test(&Uint::<2>::from(*num), SelfridgeBase, LucasCheck::Strong)
                    .is_probably_prime(),
                actual_expected_result
            );
        }
    }

    fn test_composites_a_star(numbers: &[u32], expected_result: bool) {
        for num in numbers.iter() {
            let false_positive = is_vpsp(*num);
            let actual_expected_result = if false_positive {
                true
            } else {
                expected_result
            };

            // Test both single-limb and multi-limb, just in case.
            assert_eq!(
                lucas_test(&Uint::<1>::from(*num), AStarBase, LucasCheck::LucasV)
                    .is_probably_prime(),
                actual_expected_result
            );
            assert_eq!(
                lucas_test(&Uint::<2>::from(*num), AStarBase, LucasCheck::LucasV)
                    .is_probably_prime(),
                actual_expected_result
            );
        }
    }

    fn test_composites_brute_force(numbers: &[u32], almost_extra: bool, expected_result: bool) {
        for num in numbers.iter() {
            let false_positive = if almost_extra {
                is_aeslpsp(*num)
            } else {
                is_eslpsp(*num)
            };
            let actual_expected_result = if false_positive {
                true
            } else {
                expected_result
            };
            let check = if almost_extra {
                LucasCheck::AlmostExtraStrong
            } else {
                LucasCheck::ExtraStrong
            };

            // Test both single-limb and multi-limb, just in case.
            assert_eq!(
                lucas_test(&Uint::<1>::from(*num), BruteForceBase, check).is_probably_prime(),
                actual_expected_result,
                "Brute force base, n = {num}, almost_extra = {almost_extra}",
            );
            assert_eq!(
                lucas_test(&Uint::<2>::from(*num), BruteForceBase, check).is_probably_prime(),
                actual_expected_result
            );
        }
    }

    #[test]
    fn strong_fibonacci_pseudoprimes() {
        // Can't use `test_composites()` since `STRONG_FIBONACCI` is `U64`.
        // Good thing we don't need to test for intersection
        // with `EXTRA_STRONG_LUCAS` or `STRONG_LUCAS` - there's none.
        for num in pseudoprimes::STRONG_FIBONACCI.iter() {
            assert!(!lucas_test(num, SelfridgeBase, LucasCheck::Strong).is_probably_prime());
            assert!(!lucas_test(num, BruteForceBase, LucasCheck::ExtraStrong).is_probably_prime());
        }
    }

    #[test]
    fn fibonacci_pseudoprimes() {
        test_composites_selfridge(pseudoprimes::FIBONACCI, false);
        test_composites_a_star(pseudoprimes::FIBONACCI, false);
        test_composites_brute_force(pseudoprimes::FIBONACCI, false, false);
        test_composites_brute_force(pseudoprimes::FIBONACCI, true, false);
    }

    #[test]
    fn bruckman_lucas_pseudoprimes() {
        test_composites_selfridge(pseudoprimes::BRUCKMAN_LUCAS, false);
        test_composites_a_star(pseudoprimes::BRUCKMAN_LUCAS, false);
        test_composites_brute_force(pseudoprimes::BRUCKMAN_LUCAS, false, false);
        test_composites_brute_force(pseudoprimes::BRUCKMAN_LUCAS, true, false);
    }

    #[test]
    fn almost_extra_strong_lucas_pseudoprimes() {
        test_composites_selfridge(pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS, false);
        test_composites_a_star(pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS, false);

        // Check for the difference between the almost extra strong and extra strong tests.
        test_composites_brute_force(pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS, false, false);
        test_composites_brute_force(pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS, true, true);
    }

    #[test]
    fn extra_strong_lucas_pseudoprimes() {
        test_composites_selfridge(pseudoprimes::EXTRA_STRONG_LUCAS, false);
        test_composites_a_star(pseudoprimes::EXTRA_STRONG_LUCAS, false);

        // These are the known false positives for the extra strong test
        // with brute force base selection.
        test_composites_brute_force(pseudoprimes::EXTRA_STRONG_LUCAS, false, true);
    }

    #[test]
    fn lucas_pseudoprimes() {
        test_composites_selfridge(pseudoprimes::LUCAS, false);
        test_composites_a_star(pseudoprimes::LUCAS, false);
        test_composites_brute_force(pseudoprimes::LUCAS, false, false);
        test_composites_brute_force(pseudoprimes::LUCAS, true, false);
    }

    #[test]
    fn strong_lucas_pseudoprimes() {
        // These are the known false positives for the strong test
        // with Selfridge base selection.
        test_composites_selfridge(pseudoprimes::STRONG_LUCAS, true);

        test_composites_a_star(pseudoprimes::STRONG_LUCAS, false);
        test_composites_brute_force(pseudoprimes::STRONG_LUCAS, false, false);
        test_composites_brute_force(pseudoprimes::STRONG_LUCAS, true, false);
    }

    #[test]
    fn strong_pseudoprimes_base_2() {
        // Cross-test against the pseudoprimes that circumvent the MR test base 2.
        // We expect the Lucas test to correctly classify them as composites.
        test_composites_selfridge(pseudoprimes::STRONG_BASE_2, false);
        test_composites_a_star(pseudoprimes::STRONG_BASE_2, false);
        test_composites_brute_force(pseudoprimes::STRONG_BASE_2, false, false);
        test_composites_brute_force(pseudoprimes::STRONG_BASE_2, true, false);
    }

    #[test]
    fn large_carmichael_number() {
        let p = pseudoprimes::LARGE_CARMICHAEL_NUMBER;
        assert!(!lucas_test(&p, SelfridgeBase, LucasCheck::Strong).is_probably_prime());
        assert!(!lucas_test(&p, AStarBase, LucasCheck::LucasV).is_probably_prime());
        assert!(!lucas_test(&p, BruteForceBase, LucasCheck::AlmostExtraStrong).is_probably_prime());
        assert!(!lucas_test(&p, BruteForceBase, LucasCheck::ExtraStrong).is_probably_prime());
    }

    fn test_large_primes<const L: usize>(nums: &[Uint<L>]) {
        for num in nums {
            assert!(lucas_test(num, SelfridgeBase, LucasCheck::Strong).is_probably_prime());
            assert!(lucas_test(num, AStarBase, LucasCheck::LucasV).is_probably_prime());
            assert!(
                lucas_test(num, BruteForceBase, LucasCheck::AlmostExtraStrong).is_probably_prime()
            );
            assert!(lucas_test(num, BruteForceBase, LucasCheck::ExtraStrong).is_probably_prime());
        }
    }

    #[test]
    fn large_primes() {
        test_large_primes(primes::PRIMES_128);
        test_large_primes(primes::PRIMES_256);
        test_large_primes(primes::PRIMES_384);
        test_large_primes(primes::PRIMES_512);
        test_large_primes(primes::PRIMES_1024);
    }

    #[test]
    fn test_lucas_v_pseudoprimes() {
        for num in pseudoprimes::LARGE_LUCAS_V {
            // These are false positives for Lucas-V test
            assert!(lucas_test(num, AStarBase, LucasCheck::LucasV).is_probably_prime());

            // These tests should work correctly
            assert!(!lucas_test(num, SelfridgeBase, LucasCheck::Strong).is_probably_prime());
            assert!(
                !lucas_test(num, BruteForceBase, LucasCheck::AlmostExtraStrong).is_probably_prime()
            );
            assert!(!lucas_test(num, BruteForceBase, LucasCheck::ExtraStrong).is_probably_prime());
        }
    }

    #[cfg(feature = "tests-exhaustive")]
    #[test]
    fn exhaustive() {
        // Test all the odd numbers up to the limit where we know the false positives,
        // and compare the results with the reference.
        for num in (3..pseudoprimes::EXHAUSTIVE_TEST_LIMIT).step_by(2) {
            let res_ref = is_prime64(num.into());

            let eslpsp = is_eslpsp(num);
            let aeslpsp = is_aeslpsp(num);
            let slpsp = is_slpsp(num);
            let vpsp = is_vpsp(num);

            let res = lucas_test(
                &Uint::<1>::from(num),
                BruteForceBase,
                LucasCheck::AlmostExtraStrong,
            )
            .is_probably_prime();
            let expected = aeslpsp || res_ref;
            assert_eq!(
                res, expected,
                "Brute force base, almost extra strong: n={num}, expected={expected}, actual={res}",
            );

            let res = lucas_test(
                &Uint::<1>::from(num),
                BruteForceBase,
                LucasCheck::ExtraStrong,
            )
            .is_probably_prime();
            let expected = eslpsp || res_ref;
            assert_eq!(
                res, expected,
                "Brute force base: n={num}, expected={expected}, actual={res}",
            );

            let res = lucas_test(&Uint::<1>::from(num), SelfridgeBase, LucasCheck::Strong)
                .is_probably_prime();
            let expected = slpsp || res_ref;
            assert_eq!(
                res, expected,
                "Selfridge base: n={num}, expected={expected}, actual={res}",
            );

            let res = lucas_test(&Uint::<1>::from(num), AStarBase, LucasCheck::LucasV)
                .is_probably_prime();
            let expected = vpsp || res_ref;

            assert_eq!(
                res, expected,
                "A* base, Lucas-V: n={num}, expected={expected}, actual={res}",
            );
        }
    }
}
