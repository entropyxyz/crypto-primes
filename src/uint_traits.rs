use crate::hazmat::{
    gcd::gcd,
    jacobi::{self, JacobiSymbol},
};
use core::ops::{Add, Mul, Neg, Sub};

use crypto_bigint::{
    modular::{DynResidue, DynResidueParams},
    subtle::CtOption,
    ConstChoice, Integer, Limb, NonZero, PowBoundedExp, RandomMod, Reciprocal, Uint, Word,
    Random
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
    fn random_bits(rng: &mut impl CryptoRngCore, bit_length: u32) -> Self;
    fn ct_div_rem_limb_with_reciprocal(&self, reciprocal: &Reciprocal) -> (Self, Limb);
    fn try_into_u32(&self) -> Option<u32>; // Will have to be implemented at Uint<L> level if we want to use TryFrom trait

    fn as_limbs(&self) -> &[Limb];
    fn as_words(&self) -> &[Word];
    fn div_rem_limb(&self, rhs: NonZero<Limb>) -> (Self, Limb);
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
        gcd(self, rhs)
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
        Self::shr_vartime(self, shift)
    }

    fn shl_vartime(&self, shift: u32) -> (Self, ConstChoice) {
        Self::shl_vartime(self, shift)
    }

    fn as_limbs(&self) -> &[Limb] {
        self.as_limbs()
    }

    fn random_bits(rng: &mut impl CryptoRngCore, bit_length: u32) -> Self {
        Self::random(rng).shl(Self::BITS - bit_length).0
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
