use core::ops::{Add, Mul, Neg, Sub};

use crypto_bigint::{
    modular::{DynResidue, DynResidueParams},
    subtle::{Choice, CtOption},
    CheckedAdd, Integer, Uint, Zero,
};

pub trait UintLike:
    Clone + core::fmt::Debug + Eq + From<u32> + Ord + for<'a> CheckedAdd<&'a Self> + Zero
{
    type Modular: UintModLike<Raw = Self>;

    // We can get by with non-small versions of jacobi_symbol and gcd, they don't have a big impact
    // on the performance.
    fn jacobi_symbol_small(lhs: i32, rhs: &Self) -> JacobiSymbol;
    fn gcd_small(&self, rhs: u32) -> u32;
    fn bits_vartime(&self) -> usize;
    fn bit_vartime(&self, index: usize) -> bool;
    fn trailing_ones(&self) -> usize;
    fn wrapping_mul(&self, rhs: &Self) -> Self;
    fn sqrt_vartime(&self) -> Self;
    fn is_even(&self) -> Choice;
    fn is_odd(&self) -> Choice;
    fn shr(&self, shift: usize) -> Self;
    fn one() -> Self;
}

pub trait UintModLike:
    Clone
    + Eq
    + Sized
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + Neg<Output = Self>
{
    type Raw: UintLike<Modular = Self>;
    type Params;

    fn new_params(modulus: &Self::Raw) -> CtOption<Self::Params>;
    fn new(raw: &Self::Raw, params: &Self::Params) -> Self;

    fn zero(params: &Self::Params) -> Self;
    fn one(params: &Self::Params) -> Self;
    fn square(&self) -> Self;
    fn div_by_2(&self) -> Self;
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum JacobiSymbol {
    Zero,
    One,
    MinusOne,
}

impl Neg for JacobiSymbol {
    type Output = Self;
    fn neg(self) -> Self {
        match self {
            Self::Zero => Self::Zero,
            Self::One => Self::MinusOne,
            Self::MinusOne => Self::One,
        }
    }
}

// Uint<L> impls

impl<const L: usize> UintLike for Uint<L> {
    type Modular = DynResidue<L>;

    fn jacobi_symbol_small(lhs: i32, rhs: &Self) -> JacobiSymbol {
        unimplemented!()
    }

    fn gcd_small(&self, rhs: u32) -> u32 {
        unimplemented!()
    }

    fn trailing_ones(&self) -> usize {
        Self::trailing_ones(self)
    }

    fn wrapping_mul(&self, rhs: &Self) -> Self {
        Self::wrapping_mul(self, rhs)
    }

    fn bits_vartime(&self) -> usize {
        Self::bits_vartime(self)
    }

    fn bit_vartime(&self, index: usize) -> bool {
        Self::bit_vartime(self, index)
    }

    fn sqrt_vartime(&self) -> Self {
        Self::sqrt_vartime(self)
    }

    fn is_even(&self) -> Choice {
        Integer::is_even(self)
    }

    fn is_odd(&self) -> Choice {
        Integer::is_odd(self)
    }

    fn shr(&self, shift: usize) -> Self {
        Self::shr(self, shift)
    }

    fn one() -> Self {
        Self::ONE
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
