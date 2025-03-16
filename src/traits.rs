use rand_core::CryptoRng;

/// A type producing sieves for random prime generation.
pub trait SieveFactory {
    /// The type of items returning by the sieves.
    type Item;

    /// The resulting sieve.
    type Sieve: Iterator<Item = Self::Item>;

    /// Makes a sieve given an RNG and the previous exhausted sieve (if any).
    ///
    /// Returning `None` signals that the prime generation should stop.
    fn make_sieve<R: CryptoRng + ?Sized>(
        &mut self,
        rng: &mut R,
        previous_sieve: Option<&Self::Sieve>,
    ) -> Option<Self::Sieve>;
}
