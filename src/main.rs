use crypto_bigint::{Integer, Random, U128};
use rand_core::OsRng;

fn main() {
    let max = U128::MAX;
    let s = max.trailing_ones();
    println!("{}, {}", max, s);
    println!("{}", max.clone().shr(s).0);
}
