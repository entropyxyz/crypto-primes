use rand_core::{CryptoRng, RngCore, TryCryptoRng, TryRngCore};

/// Adapter from [`CryptoRng`] to [`TryCryptoRng`]
///
/// This is pending the release of a fix availale in this PR:
/// <https://github.com/rust-random/rand/pull/1593>
#[doc(hidden)]
#[derive(Debug)]
pub struct MaybeRng<'r, R>(pub &'r mut R)
where
    R: ?Sized;

impl<R> TryRngCore for MaybeRng<'_, R>
where
    R: RngCore + ?Sized,
{
    type Error = core::convert::Infallible;

    #[inline]
    fn try_next_u32(&mut self) -> Result<u32, Self::Error> {
        Ok(self.0.next_u32())
    }
    #[inline]
    fn try_next_u64(&mut self) -> Result<u64, Self::Error> {
        Ok(self.0.next_u64())
    }
    #[inline]
    fn try_fill_bytes(&mut self, dst: &mut [u8]) -> Result<(), Self::Error> {
        self.0.fill_bytes(dst);
        Ok(())
    }
}

impl<R: CryptoRng + ?Sized> TryCryptoRng for MaybeRng<'_, R> {}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;
    #[test]
    fn test_rng() {
        let mut rng = OsRng.unwrap_err();

        let mut rng = MaybeRng(&mut rng);
        rng.try_next_u32().unwrap();
        rng.try_next_u64().unwrap();
        rng.try_fill_bytes(&mut []).unwrap();
    }
}
