pub mod hazmat;
mod presets;

pub use presets::{
    is_prime, is_prime_with_rng, is_safe_prime, is_safe_prime_with_rng, prime, prime_with_rng,
    safe_prime, safe_prime_with_rng,
};
