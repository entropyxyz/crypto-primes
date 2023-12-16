use crate::hazmat::{
    gcd,
    jacobi::{self, JacobiSymbol},
};
use core::ops::{Add, Mul, Neg, Sub};

use crypto_bigint::{
    modular::{BoxedResidue, BoxedResidueParams, DynResidue, DynResidueParams},
    subtle::CtOption,
    BoxedUint, ConstChoice, Integer, Limb, NonZero, PowBoundedExp, Random, RandomMod, Reciprocal,
    Uint, Word,
};
use rand_core::CryptoRngCore;

// would be nice to have: *Assign traits; arithmetic traits for &self (BitAnd and Shr in particular);
#[allow(missing_docs)]
pub trait UintLike: Integer + RandomMod {
    type Modular: UintModLike<Raw = Self>;

    // We can get by with non-small versions of jacobi_symbol and gcd, they don't have a big impact
    // on the performance.
    fn jacobi_symbol_small(lhs: i32, rhs: &Self) -> JacobiSymbol;
    fn gcd_small(&self, rhs: u32) -> u32;
    fn bit_vartime(&self, index: u32) -> bool;
    fn trailing_zeros(&self) -> u32;
    fn trailing_ones(&self) -> u32;
    fn wrapping_sub(&self, rhs: &Self) -> Self;
    fn wrapping_mul(&self, rhs: &Self) -> Self;
    fn sqrt_vartime(&self) -> Self;
    fn shr_vartime(&self, shift: u32) -> (Self, ConstChoice);
    fn shl_vartime(&self, shift: u32) -> (Self, ConstChoice);
    fn random_bits(rng: &mut impl CryptoRngCore, bit_length: u32, bits_precision: u32) -> Self;
    fn ct_div_rem_limb_with_reciprocal(&self, reciprocal: &Reciprocal) -> (Self, Limb);
    fn try_into_u32(&self) -> Option<u32>; // Will have to be implemented at Uint<L> level if we want to use TryFrom trait

    fn as_limbs(&self) -> &[Limb];
    fn as_words(&self) -> &[Word];
    fn div_rem_limb(&self, rhs: NonZero<Limb>) -> (Self, Limb);

    // This is necessary for making sure that every BoxedUint has the same bits_precision
    // TODO: remove them after BoxedUint can work with different bits_precision
    fn one_with_precision(bits_precision: u32) -> Self;
    fn zero_with_precision(bits_precision: u32) -> Self;
    fn widen(&self, bits_precision: u32) -> Self;
}

#[allow(missing_docs)]
pub trait UintModLike:
    core::fmt::Debug
    + Clone
    + Eq
    + Sized
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + Neg<Output = Self>
    + PowBoundedExp<Self::Raw>
{
    type Raw: UintLike<Modular = Self>;
    type Params: Clone + Eq + core::fmt::Debug;

    fn new_params(modulus: &Self::Raw) -> CtOption<Self::Params>;
    fn new(raw: &Self::Raw, params: &Self::Params) -> Self;

    fn zero(params: &Self::Params) -> Self;
    fn one(params: &Self::Params) -> Self;
    fn square(&self) -> Self;
    fn div_by_2(&self) -> Self;
}

/// Uint<L> impls
impl<const L: usize> UintLike for Uint<L> {
    type Modular = DynResidue<L>;

    fn as_words(&self) -> &[Word] {
        self.as_words()
    }

    fn jacobi_symbol_small(lhs: i32, rhs: &Self) -> JacobiSymbol {
        jacobi::jacobi_symbol(lhs, rhs)
    }

    fn gcd_small(&self, rhs: u32) -> u32 {
        gcd::gcd(self, rhs)
    }

    fn trailing_zeros(&self) -> u32 {
        Self::trailing_zeros(self)
    }

    fn trailing_ones(&self) -> u32 {
        Self::trailing_ones(self)
    }

    fn wrapping_sub(&self, rhs: &Self) -> Self {
        Self::wrapping_sub(self, rhs)
    }

    fn wrapping_mul(&self, rhs: &Self) -> Self {
        Self::wrapping_mul(self, rhs)
    }

    fn bit_vartime(&self, index: u32) -> bool {
        Self::bit_vartime(self, index)
    }

    fn sqrt_vartime(&self) -> Self {
        Self::sqrt_vartime(self)
    }

    fn shr_vartime(&self, shift: u32) -> (Self, ConstChoice) {
        Self::overflowing_shr_vartime(self, shift)
    }

    fn shl_vartime(&self, shift: u32) -> (Self, ConstChoice) {
        Self::overflowing_shl_vartime(self, shift)
    }

    fn as_limbs(&self) -> &[Limb] {
        self.as_limbs()
    }

    fn random_bits(rng: &mut impl CryptoRngCore, bit_length: u32, _bits_precision: u32) -> Self {
        let random = Self::random(rng) & Self::MAX >> (Self::BITS - bit_length);
        let random = random | Self::ONE << (bit_length - 1);
        return random;
    }

    fn ct_div_rem_limb_with_reciprocal(&self, reciprocal: &Reciprocal) -> (Self, Limb) {
        self.div_rem_limb_with_reciprocal(reciprocal)
    }

    fn try_into_u32(&self) -> Option<u32> {
        self.as_words()[0].try_into().ok()
    }

    fn div_rem_limb(&self, rhs: NonZero<Limb>) -> (Self, Limb) {
        self.div_rem_limb(rhs)
    }

    fn one_with_precision(_bits_precision: u32) -> Self {
        Self::ONE
    }

    fn zero_with_precision(_bits_precision: u32) -> Self {
        Self::ZERO
    }

    fn widen(&self, _bits_precision: u32) -> Self {
        *self
    }
}

impl<const L: usize> UintModLike for DynResidue<L> {
    type Raw = Uint<L>;
    type Params = DynResidueParams<L>;

    fn new_params(modulus: &Self::Raw) -> CtOption<Self::Params> {
        Self::Params::new(modulus)
    }

    fn new(value: &Self::Raw, params: &Self::Params) -> Self {
        Self::new(value, *params)
    }

    fn zero(params: &Self::Params) -> Self {
        Self::zero(*params)
    }

    fn one(precomputed: &Self::Params) -> Self {
        Self::one(*precomputed)
    }

    fn square(&self) -> Self {
        Self::square(self)
    }

    fn div_by_2(&self) -> Self {
        Self::div_by_2(self)
    }
}

impl UintLike for BoxedUint {
    type Modular = BoxedResidue;

    fn jacobi_symbol_small(lhs: i32, rhs: &Self) -> JacobiSymbol {
        jacobi::jacobi_symbol(lhs, rhs)
    }

    fn gcd_small(&self, rhs: u32) -> u32 {
        gcd::gcd(self, rhs)
    }

    /// TODO: BoxedUint does not implement bit_vartime
    /// TODO: this needs to be tested
    fn bit_vartime(&self, index: u32) -> bool {
        if index >= self.bits_precision() {
            return false;
        }
        (self.as_limbs()[(index / Limb::BITS) as usize].0 >> (index % Limb::BITS)) & 1 == 1
    }

    fn trailing_zeros(&self) -> u32 {
        self.trailing_zeros()
    }

    /// TODO: BoxedUint does not implement trailing_ones, but it should be doable
    /// TODO: this needs to be testesd, but this is not the implementation I would go with in
    /// crypto-bigint
    fn trailing_ones(&self) -> u32 {
        let limbs = self.as_limbs();

        let mut count = 0;
        let mut i = 0;
        let mut nonmax_limb_not_encountered = ConstChoice::TRUE;
        while i < limbs.len() {
            let l = limbs[i];
            let z = l.trailing_ones();
            let should_count: bool = nonmax_limb_not_encountered.into();
            if should_count {
                count += z;
            }
            let is_max = l.0 == Limb::MAX.0;
            if should_count && is_max {
                nonmax_limb_not_encountered = ConstChoice::TRUE;
            } else {
                nonmax_limb_not_encountered = ConstChoice::FALSE;
            }
            i += 1;
        }

        count
    }

    fn wrapping_sub(&self, rhs: &Self) -> Self {
        self.wrapping_sub(rhs)
    }

    fn wrapping_mul(&self, rhs: &Self) -> Self {
        self.wrapping_mul(rhs)
    }

    /// BoxedUint does not implement sqrt_vartime
    fn sqrt_vartime(&self) -> Self {
        self.sqrt_vartime()
    }

    /// TODO: BoxedUint::shr_vartime should behave similarly as Uint::shr_vartime; instead of
    /// returning Option, return (val, choice)
    fn shr_vartime(&self, shift: u32) -> (Self, ConstChoice) {
        let (val, overflow) = self.overflowing_shr(shift);
        if overflow.into() {
            (val, ConstChoice::TRUE)
        } else {
            (val, ConstChoice::FALSE)
        }
    }

    /// TODO: BoxedUint::shl_vartime should behave similarly as Uint::shr_vartime; instead of
    /// returning Option, return (val, choice)
    fn shl_vartime(&self, shift: u32) -> (Self, ConstChoice) {
        let (val, overflow) = self.overflowing_shl(shift);
        if overflow.into() {
            (val, ConstChoice::TRUE)
        } else {
            (val, ConstChoice::FALSE)
        }
    }

    fn random_bits(rng: &mut impl CryptoRngCore, bit_length: u32, bits_precision: u32) -> Self {
        let random = Self::random(rng, bits_precision)
            & Self::max(bits_precision)
                .overflowing_shr(bits_precision - bit_length)
                .0;
        let random = random
            | Self::one_with_precision(bits_precision)
                .overflowing_shl(bit_length - 1)
                .0;
        return random;
    }

    /// TODO: BoxedUint does not implement div_rem_limb_with_reciprocal
    /// TODO: BoxedUint does not implement shl_limb
    fn ct_div_rem_limb_with_reciprocal(&self, reciprocal: &Reciprocal) -> (Self, Limb) {
        self.div_rem_limb_with_reciprocal(reciprocal)
    }

    fn try_into_u32(&self) -> Option<u32> {
        self.as_words()[0].try_into().ok()
    }

    fn as_limbs(&self) -> &[Limb] {
        self.as_limbs()
    }

    fn as_words(&self) -> &[Word] {
        self.as_words()
    }

    /// TODO: BoxedUint does not implement div_rem_limb
    fn div_rem_limb(&self, rhs: NonZero<Limb>) -> (Self, Limb) {
        self.div_rem_limb(rhs)
    }

    fn one_with_precision(bits_precision: u32) -> Self {
        Self::one_with_precision(bits_precision)
    }

    fn zero_with_precision(bits_precision: u32) -> Self {
        Self::zero_with_precision(bits_precision)
    }

    fn widen(&self, bits_precision: u32) -> Self {
        self.widen(bits_precision)
    }
}

impl UintModLike for BoxedResidue {
    type Raw = BoxedUint;
    type Params = BoxedResidueParams;

    fn new_params(modulus: &Self::Raw) -> CtOption<Self::Params> {
        Self::Params::new(modulus.clone())
    }

    fn new(raw: &Self::Raw, params: &Self::Params) -> Self {
        Self::new(raw.widen(params.bits_precision()).clone(), params.clone())
    }

    fn zero(params: &Self::Params) -> Self {
        Self::zero(params.clone())
    }

    fn one(params: &Self::Params) -> Self {
        Self::one(params.clone())
    }

    fn square(&self) -> Self {
        self.square()
    }

    /// TODO: BoxedUint does not implement div_by_2
    fn div_by_2(&self) -> Self {
        self.div_by_2()
    }
}
