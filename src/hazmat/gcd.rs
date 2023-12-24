use crypto_bigint::{Limb, NonZero, Uint};

/// Calculates the greatest common divisor of `n` and `m`.
/// By definition, `gcd(0, m) == m`.
/// `n` must be non-zero.
pub(crate) fn gcd<const L: usize>(n: &Uint<L>, m: u32) -> u32 {
    // This is an internal function, and it will never be called with `m = 0`.
    // Allowing `m = 0` would require us to have the return type of `Uint<L>`
    // (since `gcd(n, 0) = n`).
    debug_assert!(m != 0);

    // This we can check since it doesn't affect the return type,
    // even though `n` will not be 0 either in the application.
    if n == &Uint::<L>::ZERO {
        return m;
    }

    // Normalize input: the resulting (a, b) are both small, a >= b, and b != 0.
    let (mut a, mut b): (u32, u32) = if n.bits() > u32::BITS {
        // `m` is non-zero, so we can unwrap.
        let (_quo, n) = n.div_rem_limb(NonZero::new(Limb::from(m)).unwrap());
        // `n` is a remainder of a division by `u32`, so it can be safely cast to `u32`.
        let b: u32 = n.0.try_into().unwrap();
        (m, b)
    } else {
        // In this branch `n` is 32 bits or shorter,
        // so we can safely take the first limb and cast it to u32.
        let n: u32 = n.as_words()[0].try_into().unwrap();
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
    use crypto_bigint::U128;
    use num_bigint::BigUint;
    use num_integer::Integer;
    use proptest::prelude::*;

    use super::gcd;

    #[test]
    fn corner_cases() {
        assert_eq!(gcd(&U128::from(0u64), 5), 5);
        assert_eq!(gcd(&U128::from(1u64), 11 * 13 * 19), 1);
        assert_eq!(gcd(&U128::from(7u64 * 11 * 13), 1), 1);
        assert_eq!(gcd(&U128::from(7u64 * 11 * 13), 11 * 13 * 19), 11 * 13);
    }

    prop_compose! {
        fn uint()(bytes in any::<[u8; 16]>()) -> U128 {
            U128::from_le_slice(&bytes) | U128::ONE
        }
    }

    proptest! {
        #[test]
        fn fuzzy(m in any::<u32>(), n in uint()) {
            if m == 0 {
                return Ok(());
            }

            let m_bi = BigUint::from(m);
            let n_bi = BigUint::from_bytes_be(n.to_be_bytes().as_ref());
            let gcd_ref: u32 = n_bi.gcd(&m_bi).try_into().unwrap();

            let gcd_test = gcd(&n, m);
            assert_eq!(gcd_test, gcd_ref);
        }
    }
}
