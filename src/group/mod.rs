use alga::general::AbstractGroup;
use alga::general::Operator;
use num::BigInt;
use num::BigUint;

pub mod class;
pub mod rsa;

/// Selection of an element in the generating set.
pub trait Generator<O: Operator>: AbstractGroup<O> {
  fn generator() -> Self;
}

/// Efficient computation of group inverses.
pub trait Inverse<O: Operator>: Pow<O> {
  fn efficient_inverse(&self, exp: &BigUint) -> Self;
  fn pow_signed(&self, exp: &BigInt) -> Self {
    match exp.to_biguint() {
      Some(value) => self.pow(&value),
      None => self.efficient_inverse(&(-exp).to_biguint().expect("negative BigInt expected"))
    }
  }
}

/// Efficient exponentiation in a group.
pub trait Pow<O: Operator>: AbstractGroup<O> {
  // TODO: Write default impl using repeated squaring.
  fn pow(&self, exp: &BigUint) -> Self;
}
