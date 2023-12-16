//! Jacobi symbol calculation.

use core::fmt::Display;

use crypto_bigint::{BoxedUint, Limb, NonZero, Uint, Word};

use crate::UintLike;

#[allow(missing_docs)]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum JacobiSymbol {
    Zero,
    One,
    MinusOne,
}

impl core::ops::Neg for JacobiSymbol {
    type Output = Self;
    fn neg(self) -> Self {
        match self {
            Self::Zero => Self::Zero,
            Self::One => Self::MinusOne,
            Self::MinusOne => Self::One,
        }
    }
}

// A helper trait to generalize some functions over Word and Uint.
pub(crate) trait SmallMod {
    fn mod8(&self) -> Word;
    fn mod4(&self) -> Word;
}

impl SmallMod for Word {
    fn mod8(&self) -> Word {
        self & 7
    }
    fn mod4(&self) -> Word {
        self & 3
    }
}

impl<const L: usize> SmallMod for Uint<L> {
    fn mod8(&self) -> Word {
        self.as_limbs()[0].0 & 7
    }
    fn mod4(&self) -> Word {
        self.as_limbs()[0].0 & 3
    }
}

impl SmallMod for BoxedUint {
    fn mod8(&self) -> Word {
        return self.as_limbs().get(0).unwrap().0 % 7;
    }

    fn mod4(&self) -> Word {
        return self.as_limbs().get(0).expect("Missing limbs").0 % 3;
    }
}

/// Transforms `(a/p)` -> `(r/p)` for odd `p`, where the resulting `r` is odd, and `a = r * 2^s`.
/// Takes a Jacobi symbol value, and returns `r` and the new Jacobi symbol,
/// negated if the transformation changes parity.
///
/// Note that the returned `r` is odd.
fn reduce_numerator<V: SmallMod>(j: JacobiSymbol, a: Word, p: &V) -> (JacobiSymbol, Word) {
    let p_mod_8 = p.mod8();
    let s = a.trailing_zeros();
    let j = if (s & 1) == 1 && (p_mod_8 == 3 || p_mod_8 == 5) {
        -j
    } else {
        j
    };
    (j, a >> s)
}

/// Transforms `(a/p)` -> `(p/a)` for odd and coprime `a` and `p`.
/// Takes a Jacobi symbol value, and returns the swapped pair and the new Jacobi symbol,
/// negated if the transformation changes parity.
fn swap<T: SmallMod, V: SmallMod>(j: JacobiSymbol, a: T, p: V) -> (JacobiSymbol, V, T) {
    let j = if a.mod4() == 1 || p.mod4() == 1 {
        j
    } else {
        -j
    };
    (j, p, a)
}

/// Returns the Jacobi symbol `(a/p)` given an odd `p`. Panics on even `p`.
pub(crate) fn jacobi_symbol<T: UintLike + Display + SmallMod>(a: i32, p_long: &T) -> JacobiSymbol {
    if p_long.is_even().into() {
        panic!("`p_long` must be an odd integer, but got {}", p_long);
    }

    let result = JacobiSymbol::One; // Keep track of all the sign flips here.

    // Deal with a negative `a` first:
    // (-a/n) = (-1/n) * (a/n)
    //        = (-1)^((n-1)/2) * (a/n)
    //        = (-1 if n = 3 mod 4 else 1) * (a/n)
    let (result, a_pos) = {
        let result = if a < 0 && p_long.mod4() == 3 {
            -result
        } else {
            result
        };
        (result, a.abs_diff(0))
    };

    // A degenerate case.
    if a_pos == 1 || p_long == &T::one_with_precision(p_long.bits_precision()) {
        return result;
    }

    let a_limb = Limb::from(a_pos);

    // Normalize input: at the end we want `a < p`, `p` odd, and both fitting into a `Word`.
    let (result, a, p): (JacobiSymbol, Word, Word) = if p_long.bits_vartime() <= Limb::BITS {
        let a = a_limb.0;
        let p = p_long.as_limbs()[0].0;
        (result, a % p, p)
    } else {
        let (result, a) = reduce_numerator(result, a_limb.0, p_long);
        if a == 1 {
            return result;
        }
        // NOTE: prior to the UintLike generics it was *p_long, which works because Uint implements
        // Copy. However, BoxedUint does not implement Copy (which it should not anyways since it
        // is heap-allocated), so explicit cloning is the next best thing
        let (result, a_long, p) = swap(result, a, p_long.clone());
        // Can unwrap here, since `p` is swapped with `a`,
        // and `a` would be odd after `reduce_numerator()`.
        let (_, a) = a_long.div_rem_limb(NonZero::new(Limb::from(p)).unwrap());
        (result, a.0, p)
    };

    let mut result = result;
    let mut a = a;
    let mut p = p;

    loop {
        if a == 0 {
            return JacobiSymbol::Zero;
        }

        // At this point `p` is odd (either coming from outside of the `loop`,
        // or from the previous iteration, where a previously reduced `a`
        // was swapped into its place), so we can call this.
        (result, a) = reduce_numerator(result, a, &p);

        if a == 1 {
            return result;
        }

        // At this point both `a` and `p` are odd: `p` was odd before,
        // and `a` is odd after `reduce_numerator()`.
        // Note that technically `swap()` only returns a valid `result` if `a` and `p` are coprime.
        // But if they are not, we will return `Zero` eventually,
        // which is not affected by any sign changes.
        (result, a, p) = swap(result, a, p);

        a %= p;
    }
}

#[cfg(test)]
mod tests {

    use alloc::format;

    use crypto_bigint::{Encoding, U128};
    use num_bigint::{BigInt, Sign};
    use num_modular::ModularSymbols;
    use proptest::prelude::*;

    use super::{jacobi_symbol, JacobiSymbol};

    #[test]
    fn jacobi_symbol_derived_traits() {
        assert_eq!(format!("{:?}", JacobiSymbol::One), "One");
        assert_eq!(JacobiSymbol::One, JacobiSymbol::One);
        assert!(JacobiSymbol::One != JacobiSymbol::MinusOne);
        assert_eq!(JacobiSymbol::One.clone(), JacobiSymbol::One);
    }

    #[test]
    fn jacobi_symbol_neg_zero() {
        // This does not happen during normal operation, since we return zero as soon as we get it.
        // So just covering it for the completness' sake.
        assert_eq!(-JacobiSymbol::Zero, JacobiSymbol::Zero);
    }

    #[test]
    #[should_panic(
        expected = "`p_long` must be an odd integer, but got 00000000000000000000000000000004"
    )]
    fn jacobi_symbol_p_is_even() {
        let _j = jacobi_symbol(1, &U128::from(4u32));
    }

    // Reference from `num-modular` - supports long `p`, but only positive `a`.
    fn jacobi_symbol_ref(a: i32, p: &U128) -> JacobiSymbol {
        let a_bi = BigInt::from(a);
        let p_bi = BigInt::from_bytes_be(Sign::Plus, p.to_be_bytes().as_ref());
        let j = a_bi.jacobi(&p_bi);
        if j == 1 {
            JacobiSymbol::One
        } else if j == -1 {
            JacobiSymbol::MinusOne
        } else {
            JacobiSymbol::Zero
        }
    }

    #[test]
    fn small_values() {
        // Test small values, using a reference implementation.
        for a in -31i32..31 {
            for p in (1u32..31).step_by(2) {
                let p_long = U128::from(p);
                let j_ref = jacobi_symbol_ref(a, &p_long);
                let j = jacobi_symbol(a, &p_long);
                assert_eq!(j, j_ref);
            }
        }
    }

    #[test]
    fn big_values() {
        // a = x, p = x * y, where x and y are big primes. Should give 0.
        let a = 2147483647i32; // 2^31 - 1, a prime
        let p = U128::from_be_hex("000000007ffffffeffffffe28000003b"); // (2^31 - 1) * (2^64 - 59)
        assert_eq!(jacobi_symbol(a, &p), JacobiSymbol::Zero);
        assert_eq!(jacobi_symbol_ref(a, &p), JacobiSymbol::Zero);

        // a = x^2 mod p, should give 1.
        let a = 659456i32; // Obtained from x = 2^70
        let p = U128::from_be_hex("ffffffffffffffffffffffffffffff5f"); // 2^128 - 161 - not a prime
        assert_eq!(jacobi_symbol(a, &p), JacobiSymbol::One);
        assert_eq!(jacobi_symbol_ref(a, &p), JacobiSymbol::One);

        let a = i32::MIN; // -2^31, check that no overflow occurs
        let p = U128::from_be_hex("000000007ffffffeffffffe28000003b"); // (2^31 - 1) * (2^64 - 59)
        assert_eq!(jacobi_symbol(a, &p), JacobiSymbol::One);
        assert_eq!(jacobi_symbol_ref(a, &p), JacobiSymbol::One);
    }

    prop_compose! {
        fn odd_uint()(bytes in any::<[u8; 16]>()) -> U128 {
            U128::from_le_slice(&bytes) | U128::ONE
        }
    }

    proptest! {
        #[test]
        fn fuzzy(a in any::<i32>(), p in odd_uint()) {
            let j_ref = jacobi_symbol_ref(a, &p);
            let j = jacobi_symbol(a, &p);
            assert_eq!(j, j_ref);
        }
    }
}
