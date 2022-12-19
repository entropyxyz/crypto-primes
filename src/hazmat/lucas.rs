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

pub trait LucasBase {
    fn generate<const L: usize>(n: &Uint<L>) -> Result<(u32, i32), bool>;
}

pub struct SelfridgeBase;

impl LucasBase for SelfridgeBase {
    fn generate<const L: usize>(n: &Uint<L>) -> Result<(u32, i32), bool> {
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
                // Modification of Method A by Baillie [^6]:
                // If `d == (+,-)n`, (e.g., `n` = 5 or 11) try the next `d` instead of quitting;
                // this small modification of Selfridge's method A
                // enables 5 and 11 to be classified as Lucas probable primes.
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

pub struct BruteForceBase;

impl LucasBase for BruteForceBase {
    fn generate<const L: usize>(n: &Uint<L>) -> Result<(u32, i32), bool> {
        // Baillie-OEIS "method C" for choosing D, P, Q (see [^2]):
        // Try increasing P ≥ 3 such that D = P² - 4 (so Q = 1) until Jacobi(D, n) = -1.
        // The search is expected to succeed for non-square n after just a few trials.
        // After more than expected failures, check whether n is square
        // (which would cause Jacobi(D, n) = 1 for all D not dividing n).

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
                // D = p²-4 = (p-2)(p+2).
                // If (D/n) == 0 then D shares a prime factor with n.
                // Since the loop proceeds in increasing p and starts with p-2==1,
                // the shared prime factor must be p+2.
                // If p+2 == n, then n is prime; otherwise p+2 is a proper factor of n.
                return Err(n == &Uint::<L>::from(p + 2));
            }

            attempts += 1;
            p += 1;
        }

        Ok((p, 1))
    }
}

/// Performs the extra strong Lucas primality test, as defined by Grantham[^4].
///
/// NOTE: despite the name, this is **not** a superset of "strong Lucas test"
/// defined by Baillie and Wagstaff[^1]; although it appears[^5] that it has a similar ratio
/// of pseudoprimes, and, like for the strong test, they have not been demonstrated to intersect
/// with Miller-Rabin pseudoprimes.
///
/// [^1]: R. Baillie, S. S. Wagstaff, "Lucas pseudoprimes",
///       Math. Comp. 35 1391-1417 (1980),
///       DOI: [10.2307/2006406](https://dx.doi.org/10.2307/2006406),
///       <http://mpqs.free.fr/LucasPseudoprimes.pdf>
///
/// [^2]: R. Baillie, Mathematica code for extra strong Lucas pseudoprimes,
///       <https://oeis.org/A217719/a217719.txt>
///
/// [^3]: R. Crandall, C. Pomerance, "Prime numbers: a computational perspective",
///       2nd ed., Springer (2005) (ISBN: 0-387-25282-7, 978-0387-25282-7)
///
/// [^4]: J. Grantham, "Frobenius pseudoprimes",
///       Math. Comp. 70 873-891 (2001),
///       DOI: [10.1090/S0025-5718-00-01197-2](https://dx.doi.org/10.1090/S0025-5718-00-01197-2)
///
/// [^5]: D. Jacobsen, "Pseudoprime Statistics, Tables, and Data",
///       <http://ntheory.org/pseudoprimes.html>
///
/// [^6]: R. Baillie, Mathematica program to generate terms,
///       <https://oeis.org/A217120/a217120_1.txt>
///
/// [^7]: Zhaiyu Mo, "Diophantine equations, lucas sequences and pseudoprimes",
///       graduate thesis, University of Calgary, Calgary, AB (1993)
///       DOI: [10.11575/PRISM/10820](https://dx.doi.org/10.11575/PRISM/10820)
pub fn is_strong_lucas_prime<B: LucasBase, const L: usize>(
    n: &Uint<L>,
    full_uk_check: bool,
) -> bool {
    let (p, q) = match B::generate(n) {
        Ok((p, q)) => (p, q),
        Err(is_prime) => return is_prime,
    };

    let p_is_one = p == 1;
    let q_is_one = q == 1;

    // The definition of "extra strong Lucas pseudoprime" in [^7], Def. 7.4
    // or [^4], after Thm 2.3 on p. 876 ([^4] only defines it for `Q == 1`)
    // (D, P, Q above have become Δ, b, 1):
    //
    // Let U_n = U_n(b, 1), V_n = V_n(b, 1), and Δ = b²-4.
    // An extra strong Lucas pseudoprime to base b is a composite n = 2^r s + Jacobi(Δ, n),
    // where s is odd and gcd(n, 2*Δ) = 1, such that either (i) U_s ≡ 0 mod n and V_s ≡ ±2 mod n,
    // or (ii) V_{2^t s} ≡ 0 mod n for some 0 ≤ t < r-1.
    //
    // We know gcd(n, Δ) = 1 or else we'd have found Jacobi(d, n) == 0 above.
    // We know gcd(n, 2) = 1 because n is odd.
    //
    // Arrange s = (n - Jacobi(Δ, n)) / 2^r = (n+1) / 2^r.
    let mut s = n.wrapping_add(&Uint::<L>::ONE);
    let r = trailing_zeros(&s);
    s >>= r;

    // Compute Lucas sequence V_s(1, b), where:
    //
    //  V(0) = 2
    //  V(1) = P
    //  V(k) = P V(k-1) - Q V(k-2).
    //
    // (Remember that due to method C above, P = 1, Q = b.)
    //
    // In general V(k) = α^k + β^k, where α and β are roots of x² - Px + Q.
    // [^3], eq. (3.14) observe that for 0 ≤ j ≤ k,
    //
    //  V(j+k) = V(j)V(k) - Q^j * V(k-j).
    //
    // So in particular, to quickly double the subscript:
    //
    //  V(2k) = V(k)² - 2 * Q^k
    //  V(2k+1) = V(k) V(k+1) - Q^k
    //
    // We can therefore start with k=0 and build up to k=s in log₂(s) steps.
    let params = DynResidueParams::<L>::new(*n);

    let zero_m = DynResidue::<L>::zero(params);
    let one_m = DynResidue::<L>::one(params);
    let two_m = DynResidue::<L>::new(Uint::<L>::from(2u32), params);
    let minus_two_m = zero_m.sub(&two_m);

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

    let p_m = if p_is_one {
        one_m
    } else {
        DynResidue::<L>::new(Uint::<L>::from(p), params)
    };

    let mut vk_m = two_m;
    let mut vk1_m = p_m;
    let mut qk = one_m;
    let mut qk_times_p = if p_is_one { one_m } else { p_m };

    for i in (0..s.bits_vartime()).rev() {
        if s.bit_vartime(i) == 1 {
            // k' = 2k+1

            // V(k') = V(2k+1) = V(k) V(k+1) - Q^k * P
            vk_m = vk_m.mul(&vk1_m).sub(&qk_times_p);

            // V(k'+1) = V(2k+2) = V(k+1)² - 2 * Q^(k+1)
            let qk1 = qk.mul(&q_m);

            let two_qk1 = if q_is_one { two_m } else { qk1.add(&qk1) };

            vk1_m = vk1_m.square().sub(&two_qk1);
            qk = qk.mul(&qk1);
        } else {
            // k' = 2k

            // V(k'+1) = V(2k+1) = V(k) V(k+1) - Q^k * P
            vk1_m = vk_m.mul(&vk1_m).sub(&qk_times_p);

            // V(k') = V(2k) = V(k)² - 2 * Q^k
            let two_qk1 = if q_is_one { two_m } else { qk.add(&qk) };

            vk_m = vk_m.square().sub(&two_qk1);
            qk = qk.square();
        }

        if p_is_one {
            qk_times_p = qk;
        } else {
            qk_times_p = qk.mul(&p_m);
        }
    }

    // Strong Lucas pseudorpime involves checking for `U(s) ≡ 0 mod n`.
    // Extra strong check (from [^7]) adds to that checking for `V(s) ≡ ±2 mod n`.
    // Since the latter does not require additional calculations, we check it first.
    let vk_equals_two = if q_is_one {
        vk_m == two_m || vk_m == minus_two_m
    } else {
        true
    };

    if vk_equals_two {
        // Check U(s) ≡ 0.
        // As suggested by [^5], apply eq. 3.13 from [^3]:
        //
        //  U(k) = D⁻¹ (2 V(k+1) - P V(k))
        if full_uk_check {
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
            // TODO: could compare 2 V(k+1) == P V(k)
            return true;
        }
    }

    // Now k=s, so vk = V(s).
    if vk_m == zero_m {
        return true;
    }

    // Check V(2^t s) ≡ 0 mod n for some 1 ≤ t ≤ r-1 (we checked that for `t == 0` just above).
    for _ in 1..r {
        // Optimization: V(k) = ±2 is a fixed point for V(k') = V(k)² - 2,
        // so if V(k) = 2, we can stop: we will never find a future V(k) == 0.
        if q_is_one && (vk_m == two_m || vk_m == minus_two_m) {
            return false;
        }

        // k' = 2k
        // V(k') = V(2k) = V(k)² - 2 * Q^k
        vk_m = vk_m.mul(&vk_m).sub(&qk).sub(&qk);
        qk = qk.square();

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
    use crate::hazmat::pseudoprimes;

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

    fn test_sequence_selfridge(numbers: &[u32], expected_result: bool) {
        for num in numbers.iter() {
            // Many of the lists have intersections with strong Lucas pseudoprimes,
            // and our test is expected to report them as primes.
            // Skipping those, since we will test them separately.
            let false_positive = is_slpsp(*num);
            let actual_expected_result = if false_positive {
                true
            } else {
                expected_result
            };

            // Test both single-limb and multi-limb, just in case.
            assert_eq!(
                is_prime::<SelfridgeBase, 1>(&Uint::<1>::from(*num), true),
                actual_expected_result
            );
            assert_eq!(
                is_prime::<SelfridgeBase, 2>(&Uint::<2>::from(*num), true),
                actual_expected_result
            );
        }
    }

    fn test_sequence_brute_force(numbers: &[u32], full_uk_check: bool, expected_result: bool) {
        for num in numbers.iter() {
            // Many of the lists have intersections with extra strong Lucas pseudoprimes,
            // and our test is expected to report them as primes.
            // Skipping those, since we will test them separately.
            let false_positive = if full_uk_check {
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
                is_prime::<BruteForceBase, 1>(&Uint::<1>::from(*num), full_uk_check),
                actual_expected_result,
                "Brute force base, n = {}, full_uk_check = {}",
                num,
                full_uk_check,
            );
            assert_eq!(
                is_prime::<BruteForceBase, 2>(&Uint::<2>::from(*num), full_uk_check),
                actual_expected_result
            );
        }
    }

    #[test]
    fn strong_fibonacci_pseudoprimes() {
        // Can't use `test_sequence()` since `STRONG_FIBONACCI` is `u64`,
        // and we need to be mindful of 32-bit systems which don't have `Uint::from(u64)`.
        // Good thing we don't need to test for intersection
        // with `EXTRA_STRONG_LUCAS` or `STRONG_LUCAS` - there's none.
        for num in pseudoprimes::STRONG_FIBONACCI.iter() {
            let num_hi = num >> 32;
            let num_lo = num & ((1 << 32) - 1);
            let uint2 = Uint::<2>::from(num_lo) | (Uint::<2>::from(num_hi) << 32);
            assert!(!is_prime::<SelfridgeBase, 2>(&uint2, true));
            assert!(!is_prime::<BruteForceBase, 2>(&uint2, true));
        }
    }

    #[test]
    fn fibonacci_pseudoprimes() {
        test_sequence_selfridge(pseudoprimes::FIBONACCI, false);
        test_sequence_brute_force(pseudoprimes::FIBONACCI, false, false);
        test_sequence_brute_force(pseudoprimes::FIBONACCI, true, false);
    }

    #[test]
    fn bruckman_lucas_pseudoprimes() {
        test_sequence_selfridge(pseudoprimes::BRUCKMAN_LUCAS, false);
        test_sequence_brute_force(pseudoprimes::BRUCKMAN_LUCAS, false, false);
        test_sequence_brute_force(pseudoprimes::BRUCKMAN_LUCAS, true, false);
    }

    #[test]
    fn almost_extra_strong_lucas_pseudoprimes() {
        test_sequence_selfridge(pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS, false);

        // Check specifically that we are running an extra strong test with brute force (P,Q),
        // and not an almost extra strong one (that is, we check U(s) == 0).
        // If that condition is not checked, this should fail.
        test_sequence_brute_force(pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS, true, false);
        test_sequence_brute_force(pseudoprimes::ALMOST_EXTRA_STRONG_LUCAS, false, true);
    }

    #[test]
    fn extra_strong_lucas_pseudoprimes() {
        test_sequence_selfridge(pseudoprimes::EXTRA_STRONG_LUCAS, false);

        // We expect our Lucas test to incorrectly classify them as primes.
        test_sequence_brute_force(pseudoprimes::EXTRA_STRONG_LUCAS, true, true);
    }

    #[test]
    fn lucas_pseudoprimes() {
        test_sequence_selfridge(pseudoprimes::LUCAS, false);
        test_sequence_brute_force(pseudoprimes::LUCAS, false, false);
        test_sequence_brute_force(pseudoprimes::LUCAS, true, false);
    }

    #[test]
    fn strong_lucas_pseudoprimes() {
        // We expect our Lucas test to incorrectly classify them as primes.
        test_sequence_selfridge(pseudoprimes::STRONG_LUCAS, true);

        test_sequence_brute_force(pseudoprimes::STRONG_LUCAS, false, false);
        test_sequence_brute_force(pseudoprimes::STRONG_LUCAS, true, false);
    }

    #[test]
    fn strong_pseudoprimes_base_2() {
        // Cross-test against the pseudoprimes that circumvent the MR test base 2.
        // We expect the Lucas test to correctly classify them as composites.
        test_sequence_selfridge(pseudoprimes::STRONG_BASE_2, false);

        test_sequence_brute_force(pseudoprimes::STRONG_BASE_2, false, false);
        test_sequence_brute_force(pseudoprimes::STRONG_BASE_2, true, false);
    }

    #[test]
    fn exhaustive() {
        // Test all the odd numbers up to the last extra strong Lucas pseudoprime (approximately),
        // and compare the results with the reference.
        let last = pseudoprimes::STRONG_LUCAS[pseudoprimes::STRONG_LUCAS.len() - 1];
        for num in (3..last).step_by(2) {
            let res_ref = num.is_prime();

            let eslpsp = is_eslpsp(num);
            let aeslpsp = is_aeslpsp(num);
            let slpsp = is_slpsp(num);

            let res = is_prime::<BruteForceBase, 1>(&Uint::<1>::from(num), false);
            let expected = aeslpsp || res_ref;
            assert_eq!(
                res, expected,
                "Brute force base, almost extra strong: n={}, expected={}, actual={}",
                num, expected, res,
            );

            let res = is_prime::<BruteForceBase, 1>(&Uint::<1>::from(num), true);
            let expected = eslpsp || res_ref;
            assert_eq!(
                res, expected,
                "Brute force base: n={}, expected={}, actual={}",
                num, expected, res,
            );

            let res = is_prime::<SelfridgeBase, 1>(&Uint::<1>::from(num), true);
            let expected = slpsp || res_ref;
            assert_eq!(
                res, expected,
                "Selfridge base: n={}, expected={}, actual={}",
                num, expected, res,
            );
        }
    }
}
