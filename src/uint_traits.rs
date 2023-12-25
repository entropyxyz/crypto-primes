#![allow(missing_docs)]

use crypto_bigint::{subtle::CtOption, BoxedUint, Integer, Random, RandomMod, Uint};
use rand_core::CryptoRngCore;

// would be nice to have: *Assign traits; arithmetic traits for &self (BitAnd and Shr in particular);
pub trait UintLike: Integer + RandomMod {
    fn bit_vartime(&self, index: u32) -> bool;
    fn set_bit_vartime(&mut self, index: u32, value: bool);
    fn trailing_zeros_vartime(&self) -> u32;
    fn trailing_ones_vartime(&self) -> u32;
    fn sqrt_vartime(&self) -> Self;
    fn overflowing_shl_vartime(&self, shift: u32) -> CtOption<Self>;
    fn overflowing_shr_vartime(&self, shift: u32) -> CtOption<Self>;
    fn wrapping_shl_vartime(&self, shift: u32) -> Self;
    fn wrapping_shr_vartime(&self, shift: u32) -> Self;
    fn random_bits(rng: &mut impl CryptoRngCore, bit_length: u32, bits_precision: u32) -> Self;
}

// Uint<L> impls

impl<const L: usize> UintLike for Uint<L> {
    fn set_bit_vartime(&mut self, index: u32, value: bool) {
        if value {
            *self |= Uint::ONE << index
        } else {
            *self &= (Uint::ONE << index).not()
        }
    }

    fn trailing_zeros_vartime(&self) -> u32 {
        Self::trailing_zeros_vartime(self)
    }

    fn trailing_ones_vartime(&self) -> u32 {
        Self::trailing_ones_vartime(self)
    }

    fn bit_vartime(&self, index: u32) -> bool {
        Self::bit_vartime(self, index)
    }

    fn sqrt_vartime(&self) -> Self {
        Self::sqrt_vartime(self)
    }

    fn overflowing_shl_vartime(&self, shift: u32) -> CtOption<Self> {
        Self::overflowing_shl_vartime(self, shift).into()
    }

    fn overflowing_shr_vartime(&self, shift: u32) -> CtOption<Self> {
        Self::overflowing_shr_vartime(self, shift).into()
    }

    fn wrapping_shl_vartime(&self, shift: u32) -> Self {
        Self::wrapping_shl_vartime(self, shift)
    }

    fn wrapping_shr_vartime(&self, shift: u32) -> Self {
        Self::wrapping_shr_vartime(self, shift)
    }

    /// TODO: bits_precision is required because BoxedUint::random requires bits_precision
    ///   we can require the user to input the correct precision:
    ///   in this case the documentation will need to communicate to the users about this
    ///   requirement, and we could put an assert statement to check
    ///
    /// We can also accept whatever value and just ignore it.
    fn random_bits(rng: &mut impl CryptoRngCore, bit_length: u32, _bits_precision: u32) -> Self {
        if bit_length > Self::BITS {
            panic!("The requested bit length ({bit_length}) is larger than the chosen Uint size");
        }
        let random = Self::random(rng);
        random >> (Self::BITS - bit_length)
    }
}

impl UintLike for BoxedUint {
    fn set_bit_vartime(&mut self, index: u32, value: bool) {
        if value {
            *self |= Self::one() << index
        } else {
            *self &= (Self::one() << index).not()
        }
    }

    fn trailing_zeros_vartime(&self) -> u32 {
        Self::trailing_zeros(self)
    }

    fn trailing_ones_vartime(&self) -> u32 {
        Self::trailing_ones_vartime(self)
    }

    fn bit_vartime(&self, index: u32) -> bool {
        Self::bit_vartime(self, index)
    }

    fn sqrt_vartime(&self) -> Self {
        Self::sqrt_vartime(self)
    }

    fn overflowing_shl_vartime(&self, shift: u32) -> CtOption<Self> {
        let (res, is_some) = Self::overflowing_shl(self, shift);
        CtOption::new(res, is_some)
    }

    fn overflowing_shr_vartime(&self, shift: u32) -> CtOption<Self> {
        let (res, is_some) = Self::overflowing_shr(self, shift);
        CtOption::new(res, is_some)
    }

    fn wrapping_shl_vartime(&self, shift: u32) -> Self {
        Self::wrapping_shl_vartime(self, shift)
    }

    fn wrapping_shr_vartime(&self, shift: u32) -> Self {
        Self::wrapping_shr_vartime(self, shift)
    }

    fn random_bits(rng: &mut impl CryptoRngCore, bit_length: u32, bits_precision: u32) -> Self {
        if bit_length > bits_precision {
            panic!("The requested bit length ({bit_length}) is larger than the chosen Uint size");
        }
        let random = Self::random(rng, bits_precision);
        random >> (bits_precision - bit_length)
    }
}
