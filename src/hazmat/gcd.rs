use crypto_bigint::{Integer, Limb, NonZero, Word};

/// Calculates the greatest common divisor of `n` and `m`.
/// By definition, `gcd(0, m) == m`.
/// `n` must be non-zero.
pub(crate) fn gcd_vartime<T: Integer>(n: &T, m: Word) -> Word {
    // This is an internal function, and it will never be called with `m = 0`.
    // Allowing `m = 0` would require us to have the return type of `Uint<L>`
    // (since `gcd(n, 0) = n`).
    debug_assert!(m != 0);

    // This we can check since it doesn't affect the return type,
    // even though `n` will not be 0 either in the application.
    if n.is_zero().into() {
        return m;
    }

    // Normalize input: the resulting (a, b) are both small, a >= b, and b != 0.
    let (mut a, mut b): (Word, Word) = if n.bits() > Word::BITS {
        // `m` is non-zero, so we can unwrap.
        let r = n.rem_limb(NonZero::new(Limb::from(m)).expect("divisor should be non-zero here"));
        (m, r.0)
    } else {
        // In this branch `n` is `Word::BITS` bits or shorter,
        // so we can safely take the first limb.
        let n = n.as_ref()[0].0;
        if n > m {
            (n, m)
        } else {
            (m, n)
        }
    };

    // Euclidean algorithm.
    // Binary GCD algorithm could be used here,
    // but the performance impact of this code is negligible.
    loop {
        let r = a % b;
        if r == 0 {
            return b;
        }
        (a, b) = (b, r)
    }
}

#[cfg(test)]
mod tests {
    use crypto_bigint::{Word, U128};
    use num_bigint::BigUint;
    use num_integer::Integer;
    use proptest::prelude::*;

    use super::gcd_vartime;

    #[test]
    fn corner_cases() {
        assert_eq!(gcd_vartime(&U128::from(0u64), 5), 5);
        assert_eq!(gcd_vartime(&U128::from(1u64), 11 * 13 * 19), 1);
        assert_eq!(gcd_vartime(&U128::from(7u64 * 11 * 13), 1), 1);
        assert_eq!(gcd_vartime(&U128::from(7u64 * 11 * 13), 11 * 13 * 19), 11 * 13);
    }

    prop_compose! {
        fn uint()(bytes in any::<[u8; 16]>()) -> U128 {
            U128::from_le_slice(&bytes) | U128::ONE
        }
    }

    proptest! {
        #[test]
        fn fuzzy(m in any::<Word>(), n in uint()) {
            if m == 0 {
                return Ok(());
            }

            let m_bi = BigUint::from(m);
            let n_bi = BigUint::from_bytes_be(n.to_be_bytes().as_ref());
            let gcd_ref: Word = n_bi.gcd(&m_bi).try_into().unwrap();

            let gcd_test = gcd_vartime(&n, m);
            assert_eq!(gcd_test, gcd_ref);
        }
    }
}
