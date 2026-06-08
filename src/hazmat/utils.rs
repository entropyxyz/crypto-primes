use crypto_bigint::{Unsigned, Word};

pub(crate) fn first_limb<T: Unsigned>(num: &T) -> Word {
    num.as_limbs().first().expect("a big integer has at least one limb").0
}

pub(crate) fn equals_primitive<T>(num: &T, primitive: Word) -> bool
where
    T: Unsigned,
{
    num.bits_vartime() <= Word::BITS && first_limb(num) == primitive
}
