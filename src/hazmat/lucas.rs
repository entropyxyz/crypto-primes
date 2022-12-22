//! Lucas primality test.
use crypto_bigint::{
    modular::{
        runtime_mod::{DynResidue, DynResidueParams},
        AddResidue, InvResidue, MulResidue, SubResidue,
    },
    Integer, Limb, Uint,
};

use super::jacobi::{jacobi_symbol, JacobiSymbol};

/// Returns the number of least-significant bits that are zero
fn trailing_zeros<B: Clone + Integer + core::ops::ShrAssign<usize>>(n: &B) -> usize {
    let mut i = 0_usize;
    let mut t = *n;
    while t.is_even().into() {
        i += 1;
        t >>= 1_usize;
    }
    i
}

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
// can be found in [^1], section 7; for the brute force method it seems to be about the same).
const ATTEMPTS_BEFORE_SQRT: u32 = 30;

/// A method for selecting the base `(P, Q)` for the Lucas primality test.
pub trait LucasBase {
    /// Given an odd integer, returns `Ok((P, Q))` on success,
    /// or `Err(is_prime)` if the primality for the given integer was discovered
    /// during the search for the base.
    fn generate<const L: usize>(&self, n: &Uint<L>) -> Result<(u32, i32), bool>;
}

/// "Method A" for selecting the base given in Baillie & Wagstaff[^Baillie1980],
/// attributed to Selfridge.
///
/// [^Baillie1980]:
///   R. Baillie, S. S. Wagstaff, "Lucas pseudoprimes",
///   Math. Comp. 35 1391-1417 (1980),
///   DOI: [10.2307/2006406](https://dx.doi.org/10.2307/2006406),
///   <http://mpqs.free.fr/LucasPseudoprimes.pdf>
pub struct SelfridgeBase;

impl LucasBase for SelfridgeBase {
    fn generate<const L: usize>(&self, n: &Uint<L>) -> Result<(u32, i32), bool> {
        // Try D = 1 - 4Q = 5, -7, 9, -11, 13, ... until Jacobi(D, n) = -1.
        // Return P = 1, Q = (1 - D) / 4).

        let mut d = 5_i32;
        let n_is_small = n.bits_vartime() < (Limb::BIT_SIZE - 1);
        let small_n = n.as_words()[0] as u32;
        let mut attempts = 0;
        loop {
            if attempts >= MAX_ATTEMPTS {
                panic!("internal error: cannot find (D/n) = -1 for {:?}", n)
            }

            if attempts >= ATTEMPTS_BEFORE_SQRT {
                let sqrt_n = n.sqrt();
                if &sqrt_n.wrapping_mul(&sqrt_n) == n {
                    return Err(false);
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
                    return Err(false);
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

/// "Method C" for selecting the base given by Baillie[^Baillie].
///
/// [^Baillie]:
///   R. Baillie, Mathematica code for extra strong Lucas pseudoprimes,
///   <https://oeis.org/A217719/a217719.txt>
pub struct BruteForceBase;

impl LucasBase for BruteForceBase {
    fn generate<const L: usize>(&self, n: &Uint<L>) -> Result<(u32, i32), bool> {
        // Try P = 3, 4, 5, ... until (D/n) = -1, where D = P^2 - 4Q.
        // Return the found P, and Q = 1.

        let mut p = 3_u32;
        let mut attempts = 0;

        loop {
            if attempts >= MAX_ATTEMPTS {
                panic!("internal error: cannot find (D/n) = -1 for {:?}", n)
            }

            if attempts >= ATTEMPTS_BEFORE_SQRT {
                let sqrt_n = n.sqrt();
                if &sqrt_n.wrapping_mul(&sqrt_n) == n {
                    return Err(false);
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
                return Err(n == &Uint::<L>::from(p + 2));
            }

            attempts += 1;
            p += 1;
        }

        Ok((p, 1))
    }
}

/// Performs the strong Lucas primality test by Baillie & Wagstaff[^Baillie1980],
/// with the "extra strong" check for available bases
/// (defined by Mo[^Mo1993], also described by Grantham[^Grantham2001]).
///
/// In more detail, there are three types of checks that can be done,
/// given the Lucas sequence built from some base `(P, Q)`
/// (determined by the choice of [`LucasBase`]) up to the elements `V(d)`, `U(d)`
/// where `d * 2^s == n - (D/n)`, `d` odd, and `D = P^2 - 4Q`:
///
/// 1. Check that that any of `V(d*2^r) == 0` for `0 <= r < s`.
/// 2. Check that `U(d) == 0`.
/// 3. Check that `V(d) == ±2` (only valid when `Q == 1`, and not performed otherwise).
///
/// If 1 or 2 are true, `n` is a "strong Lucas probable prime"[^Baillie1980].
/// If the base is [`SelfridgeBase`], known false positives constitute OEIS:A217255[^A217255].
///
/// If 1, or 2 and 3 together are true (and the base was chosen such that `Q == 1`),
/// `n` is an "extra strong Lucas probable prime"[^Mo1993].
/// If the base is [`BruteForceBase`], known false positives constitute OEIS:A217719[^A217719].
///
/// If 1 or 3 are true (and the base was chosen such that `Q == 1`),
/// `n` is an "almost extra strong Lucas probable prime".
/// If the base is [`BruteForceBase`], some known false positives are listed by Jacobsen[^Jacobsen].
///
/// One can disable the check 2 by passing `check_u = false`.
/// This option is mainly intended for testing with known false positives;
/// make sure you know the risks if you do it.
///
/// [^Baillie1980]:
///   R. Baillie, S. S. Wagstaff, "Lucas pseudoprimes",
///   Math. Comp. 35 1391-1417 (1980),
///   DOI: [10.2307/2006406](https://dx.doi.org/10.2307/2006406),
///   <http://mpqs.free.fr/LucasPseudoprimes.pdf>
///
/// [^A217255]: <https://oeis.org/A217255>
///
/// [^A217719]: <https://oeis.org/A217719>
///
/// [^Crandall2005]:
///   R. Crandall, C. Pomerance, "Prime numbers: a computational perspective",
///   2nd ed., Springer (2005) (ISBN: 0-387-25282-7, 978-0387-25282-7)
///
/// [^Grantham2001]:
///   J. Grantham, "Frobenius pseudoprimes",
///   Math. Comp. 70 873-891 (2001),
///   DOI: [10.1090/S0025-5718-00-01197-2](https://dx.doi.org/10.1090/S0025-5718-00-01197-2)
///
/// [^Jacobsen]:
///   D. Jacobsen, "Pseudoprime Statistics, Tables, and Data",
///   <http://ntheory.org/pseudoprimes.html>
///
/// [^Mo1993]:
///   Zhaiyu Mo, "Diophantine equations, Lucas sequences and pseudoprimes",
///   graduate thesis, University of Calgary, Calgary, AB (1993)
///   DOI: [10.11575/PRISM/10820](https://dx.doi.org/10.11575/PRISM/10820)
pub fn is_strong_lucas_prime<const L: usize>(
    n: &Uint<L>,
    base: impl LucasBase,
    check_u: bool,
) -> bool {
    if n.is_even().into() {
        return false;
    }

    // Find the base for the Lucas sequence.
    let (p, q) = match base.generate(n) {
        Ok((p, q)) => (p, q),
        Err(is_prime) => return is_prime,
    };

    // If either is true, it allows us to optimize certain parts of the calculations.
    let p_is_one = p == 1;
    let q_is_one = q == 1;

    // The definitions can be found in:
    // - "strong pseudoprime": [^Baillie1980], Section 3.
    // - "extra strong pseudoprime" (esprp): [^Mo1993], Def. 7.4,
    //   and [^Grantham2001], definition after Thm 2.3.
    // If `Q == 1`, every n satisfying "extra strong" conditions
    // also satisfies the "strong" ones for the same `P`.

    // Both of the definitions require gcd(n, 2QD) == 1.
    // We know gcd(n, D) = 1 by construction of the base (D is chosen such that (D/n) != 0).
    // We know gcd(n, 2) = 1 because n is odd.
    // TODO (#1): make sure that the logic here is correct:
    // If the checks below succeed, gcd(n, Q) == 1
    // (proved in [^Baillie1980], right after the definition of "strong pseudoprime")

    // Find d and s, such that d is odd and d * 2^s = (n - (D/n)).
    let mut d = n.wrapping_add(&Uint::<L>::ONE);
    let s = trailing_zeros(&d);
    d >>= s;

    // Some constants in Montgomery form

    let params = DynResidueParams::<L>::new(*n);

    let zero_m = DynResidue::<L>::zero(params);
    let one_m = DynResidue::<L>::one(params);
    let two_m = one_m.add(&one_m);
    let minus_two_m = zero_m.sub(&two_m);

    // Convert Q to Montgomery form

    let q_m = if q_is_one {
        one_m
    } else {
        let abs_q_m = DynResidue::<L>::new(Uint::<L>::from(q.abs_diff(0)), params);
        if q < 0 {
            DynResidue::<L>::zero(params).sub(&abs_q_m)
        } else {
            abs_q_m
        }
    };

    // Convert P to Montgomery form

    let p_m = if p_is_one {
        one_m
    } else {
        DynResidue::<L>::new(Uint::<L>::from(p), params)
    };

    // Compute d-th element of Lucas sequence V_d(P, Q), where:
    //
    //  V_0 = 2
    //  V_1 = P
    //  V_k = P V_{k-1} - Q V_{k-2}.
    //
    // In general V(k) = α^k + β^k, where α and β are roots of x^2 - Px + Q.
    // [^Crandall2005], eq. (3.14) observe that for 0 <= j <= k,
    //
    //  V_{j+k} = V_j V_k - Q^j * V_{k-j}.
    //
    // So in particular, to quickly double the subscript:
    //
    //  V_{2k} = V_k^2 - 2 * Q^k
    //  V_{2k+1} = V_k V_{k+1} - Q^k
    //
    // We can therefore start with k=0 and build up to k=d in log2(d) steps.

    let mut vk_m = two_m; // keeps V_k
    let mut vk1_m = p_m; // keeps V_{k+1}
    let mut qk = one_m; // keeps Q^k
    let mut qk_times_p = if p_is_one { one_m } else { p_m }; // keeps P Q^{k}

    for i in (0..d.bits_vartime()).rev() {
        if d.bit_vartime(i) == 1 {
            // k' = 2k+1

            // V_k' = V_{2k+1} = V_k V_{k+1} - P Q^k
            vk_m = vk_m.mul(&vk1_m).sub(&qk_times_p);

            // V_{k'+1} = V_{2k+2} = V_{k+1}^2 - 2 Q^{k+1}
            let qk1 = qk.mul(&q_m); // Q^{k+1}
            let two_qk1 = if q_is_one { two_m } else { qk1.add(&qk1) }; // 2 Q^{k+1}
            vk1_m = vk1_m.square().sub(&two_qk1);
            qk = qk.mul(&qk1);
        } else {
            // k' = 2k

            // V_{k'+1} = V_{2k+1} = V_k V_{k+1} - P Q^k
            vk1_m = vk_m.mul(&vk1_m).sub(&qk_times_p);

            // V_k' = V_{2k} = V_k^2 - 2 Q^k
            let two_qk = if q_is_one { two_m } else { qk.add(&qk) }; // 2 Q^k
            vk_m = vk_m.square().sub(&two_qk);
            qk = qk.square();
        }

        if p_is_one {
            qk_times_p = qk;
        } else {
            qk_times_p = qk.mul(&p_m);
        }
    }

    // Now k=d, so vk = V_d, vk_1 = V_{d+1}.

    // Extra strong check (from [^Mo1993]): `V_d == ±2 mod n`.
    // Do it first since it is cheap.
    //
    // Note that it only applies if Q = 1, since it is a consequence
    // of a property of Lucas series: V_k^2 - 4 Q^k = D U_k^2 mod n.
    // If Q = 1 we can easily decompose the left side of the equation leading to the check above.
    let vk_equals_two = if q_is_one {
        vk_m == two_m || vk_m == minus_two_m
    } else {
        true
    };

    if vk_equals_two {
        // Strong check:`U_d == 0 mod n`.
        // As suggested by [^Jacobsen], apply Eq. (3.13) from [^Crandall2005]:
        //
        //  U_k = D^{-1} (2 V_{k+1} - P V_k)
        //
        // Some implementations just test for 2 V_{k+1} == P V_{k},
        // but we don't have any reference pseudoprime lists for this, so we are not doing it.
        if check_u {
            let d = (p * p) as i32 - 4 * q;
            let abs_d_m = DynResidue::<L>::new(Uint::<L>::from(d.abs_diff(0)), params);
            let d_m = if d < 0 {
                DynResidue::<L>::zero(params).sub(&abs_d_m)
            } else {
                abs_d_m
            };
            let inv_d_m = d_m.inv().unwrap();

            let vk_times_p = if p_is_one { vk_m } else { vk_m.mul(&p_m) };
            let uk_m = inv_d_m.mul(&(vk1_m.add(&vk1_m).sub(&vk_times_p)));

            if uk_m == zero_m {
                return true;
            }
        } else {
            // This is "almost extra strong check": we only checked for `V_d ` earlier.
            return true;
        }
    }

    // Check if V_{2^d t} == 0 mod n for some 0 <= t < s.

    if vk_m == zero_m {
        return true;
    }

    for _ in 1..s {
        // Optimization: V_k = ±2 is a fixed point for V_k' = V_k^2 - 2 Q^k with Q = 1,
        // so if V_k = ±2, we can stop: we will never find a future V_k == 0.
        if q_is_one && (vk_m == two_m || vk_m == minus_two_m) {
            return false;
        }

        // k' = 2k
        // V(k') = V(2k) = V(k)² - 2 * Q^k
        vk_m = vk_m.mul(&vk_m).sub(&qk).sub(&qk);

        if !q_is_one {
            qk = qk.square();
        }

        if vk_m == zero_m {
            return true;
        }
    }

    false
}

#[cfg(test)]
mod tests {
    use crypto_bigint::Uint;
    use number_theory::NumberTheory;

    use super::{is_strong_lucas_prime as is_prime, BruteForceBase, SelfridgeBase};
    use crate::hazmat::{primes, pseudoprimes};

    fn is_slpsp(num: u32) -> bool {
        pseudoprimes::STRONG_LUCAS
            .iter()
            .position(|x| *x == num)
            .is_some()
    }

    fn is_aeslpsp(num: u32) -> bool {
        pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS
            .iter()
            .position(|x| *x == num)
            .is_some()
    }

    fn is_eslpsp(num: u32) -> bool {
        pseudoprimes::EXTRA_STRONG_LUCAS
            .iter()
            .position(|x| *x == num)
            .is_some()
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
                is_prime(&Uint::<1>::from(*num), SelfridgeBase, true),
                actual_expected_result
            );
            assert_eq!(
                is_prime(&Uint::<2>::from(*num), SelfridgeBase, true),
                actual_expected_result
            );
        }
    }

    fn test_composites_brute_force(numbers: &[u32], check_u: bool, expected_result: bool) {
        for num in numbers.iter() {
            let false_positive = if check_u {
                is_eslpsp(*num)
            } else {
                is_aeslpsp(*num)
            };
            let actual_expected_result = if false_positive {
                true
            } else {
                expected_result
            };

            // Test both single-limb and multi-limb, just in case.
            assert_eq!(
                is_prime(&Uint::<1>::from(*num), BruteForceBase, check_u),
                actual_expected_result,
                "Brute force base, n = {}, check_u = {}",
                num,
                check_u,
            );
            assert_eq!(
                is_prime(&Uint::<2>::from(*num), BruteForceBase, check_u),
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
            assert!(!is_prime(num, SelfridgeBase, true));
            assert!(!is_prime(num, BruteForceBase, true));
        }
    }

    #[test]
    fn fibonacci_pseudoprimes() {
        test_composites_selfridge(pseudoprimes::FIBONACCI, false);
        test_composites_brute_force(pseudoprimes::FIBONACCI, false, false);
        test_composites_brute_force(pseudoprimes::FIBONACCI, true, false);
    }

    #[test]
    fn bruckman_lucas_pseudoprimes() {
        test_composites_selfridge(pseudoprimes::BRUCKMAN_LUCAS, false);
        test_composites_brute_force(pseudoprimes::BRUCKMAN_LUCAS, false, false);
        test_composites_brute_force(pseudoprimes::BRUCKMAN_LUCAS, true, false);
    }

    #[test]
    fn almost_extra_strong_lucas_pseudoprimes() {
        test_composites_selfridge(pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS, false);

        // Check for the difference between the almost extra strong and extra strong tests.
        test_composites_brute_force(pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS, true, false);
        test_composites_brute_force(pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS, false, true);
    }

    #[test]
    fn extra_strong_lucas_pseudoprimes() {
        test_composites_selfridge(pseudoprimes::EXTRA_STRONG_LUCAS, false);

        // These are the known false positives for the extra strong test
        // with brute force base selection.
        test_composites_brute_force(pseudoprimes::EXTRA_STRONG_LUCAS, true, true);
    }

    #[test]
    fn lucas_pseudoprimes() {
        test_composites_selfridge(pseudoprimes::LUCAS, false);
        test_composites_brute_force(pseudoprimes::LUCAS, false, false);
        test_composites_brute_force(pseudoprimes::LUCAS, true, false);
    }

    #[test]
    fn strong_lucas_pseudoprimes() {
        // These are the known false positives for the strong test
        // with Selfridge base selection.
        test_composites_selfridge(pseudoprimes::STRONG_LUCAS, true);

        test_composites_brute_force(pseudoprimes::STRONG_LUCAS, false, false);
        test_composites_brute_force(pseudoprimes::STRONG_LUCAS, true, false);
    }

    #[test]
    fn strong_pseudoprimes_base_2() {
        // Cross-test against the pseudoprimes that circumvent the MR test base 2.
        // We expect the Lucas test to correctly classify them as composites.
        test_composites_selfridge(pseudoprimes::STRONG_BASE_2, false);
        test_composites_brute_force(pseudoprimes::STRONG_BASE_2, false, false);
        test_composites_brute_force(pseudoprimes::STRONG_BASE_2, true, false);
    }

    #[test]
    fn large_carmichael_number() {
        let p = pseudoprimes::LARGE_CARMICHAEL_NUMBER;
        assert!(!is_prime(&p, SelfridgeBase, true));
        assert!(!is_prime(&p, BruteForceBase, false));
        assert!(!is_prime(&p, BruteForceBase, true));
    }

    fn test_large_primes<const L: usize>(nums: &[Uint<L>]) {
        for num in nums {
            assert!(is_prime(&num, SelfridgeBase, true));
            assert!(is_prime(&num, BruteForceBase, false));
            assert!(is_prime(&num, BruteForceBase, true));
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
    fn exhaustive() {
        // Test all the odd numbers up to the limit where we know the false positives,
        // and compare the results with the reference.
        for num in (3..pseudoprimes::EXHAUSTIVE_TEST_LIMIT).step_by(2) {
            let res_ref = num.is_prime();

            let eslpsp = is_eslpsp(num);
            let aeslpsp = is_aeslpsp(num);
            let slpsp = is_slpsp(num);

            let res = is_prime(&Uint::<1>::from(num), BruteForceBase, false);
            let expected = aeslpsp || res_ref;
            assert_eq!(
                res, expected,
                "Brute force base, almost extra strong: n={}, expected={}, actual={}",
                num, expected, res,
            );

            let res = is_prime(&Uint::<1>::from(num), BruteForceBase, true);
            let expected = eslpsp || res_ref;
            assert_eq!(
                res, expected,
                "Brute force base: n={}, expected={}, actual={}",
                num, expected, res,
            );

            let res = is_prime(&Uint::<1>::from(num), SelfridgeBase, true);
            let expected = slpsp || res_ref;
            assert_eq!(
                res, expected,
                "Selfridge base: n={}, expected={}, actual={}",
                num, expected, res,
            );
        }
    }
}
