//! RSA (3x5) group using GMP integers in the `rug` crate.
use super::{ElemFrom, Group, UnknownOrderGroup};
use crate::util::{int, TypeRep};
use rug::Integer;
use std::str::FromStr;

#[allow(clippy::module_name_repetitions)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
/// RSA-3x5 group implementation.
/// **Note**: If you want to use `Rsa3x5` outside the context of this crate, be advised that
/// it treats `x` and `-x` as the same element for sound proofs-of-exponentiation.
/// See BBF (page 9).
pub enum Rsa3x5 {}

/// RSA-3x5 modulus.
const RSA3X5_MODULUS_DECIMAL: &str = "15";

lazy_static! {
    static ref RSA3X5_MODULUS: Integer = Integer::from_str(RSA3X5_MODULUS_DECIMAL).unwrap();
    static ref HALF_MODULUS: Integer = RSA3X5_MODULUS.clone() / 2;
}

#[allow(clippy::module_name_repetitions)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
/// An RSA 3x5 group element, directly wrapping a GMP integer from the `rug` crate.
pub struct Rsa3x5Elem(Integer);

impl TypeRep for Rsa3x5 {
    type Rep = Integer;
    fn rep() -> &'static Self::Rep {
        &RSA3X5_MODULUS
    }
}

impl Group for Rsa3x5 {
    type Elem = Rsa3x5Elem;

    fn id_(_: &Integer) -> Rsa3x5Elem {
        Self::elem(1)
    }

    fn op_(modulus: &Integer, a: &Rsa3x5Elem, b: &Rsa3x5Elem) -> Rsa3x5Elem {
        Self::elem(int(&a.0 * &b.0) % modulus)
    }

    fn exp_(modulus: &Integer, x: &Rsa3x5Elem, n: &Integer) -> Rsa3x5Elem {
        // A side-channel resistant impl is 40% slower; we'll consider it in the future if we need to.
        Self::elem(x.0.pow_mod_ref(n, modulus).unwrap())
    }

    fn inv_(modulus: &Integer, x: &Rsa3x5Elem) -> Rsa3x5Elem {
        Self::elem(x.0.invert_ref(modulus).unwrap())
    }
}

impl<T> ElemFrom<T> for Rsa3x5
    where
        Integer: From<T>,
{
    fn elem(t: T) -> Rsa3x5Elem {
        let modulus = Self::rep();
        let val = int(t) % modulus;
        if val > *HALF_MODULUS {
            Rsa3x5Elem(<(Integer, Integer)>::from((-val).div_rem_euc_ref(&modulus)).1)
        } else {
            Rsa3x5Elem(val)
        }
    }
}

impl UnknownOrderGroup for Rsa3x5 {
    fn unknown_order_elem_(_: &Integer) -> Rsa3x5Elem {
        Self::elem(2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rsa3x5_init() {
        let x = Rsa3x5::rep();
        println!("Rsa3x5::rep() - {:?}", x);
    }

    #[test]
    fn test_rsa3x5_op() {
        let a = Rsa3x5::op(&Rsa3x5::elem(2), &Rsa3x5::elem(3));
        println!("Rsa3x5:  2*3= {:?}", a);
        assert_eq!(a, Rsa3x5::elem(6));

        let b = Rsa3x5::op(&Rsa3x5::elem(-2), &Rsa3x5::elem(-3));
        println!("Rsa3x5: -2*(-3)= {:?}", b);
        assert_eq!(b, Rsa3x5::elem(6));

        let c = Rsa3x5::op(&Rsa3x5::elem(-2), &Rsa3x5::elem(3));
        println!("Rsa3x5: -2*3={:?}", c);
        assert_eq!(c, Rsa3x5::elem(-6));

        let d = Rsa3x5::op(&Rsa3x5::elem(4), &Rsa3x5::elem(7));
        println!("Rsa3x5:  4*7= {:?}", d);
        assert_eq!(d, Rsa3x5::elem(13));
    }

    /// Tests that `-x` and `x` are treated as the same element.
    #[test]
    fn test_rsa3x5_cosets() {
        assert_eq!(Rsa3x5::elem(3), Rsa3x5::elem(RSA3X5_MODULUS.clone() - 3));
        println!(
            "Rsa3x5::elem(3): {:?} \n\
             Rsa3x5::elem(RSA3X5_MODULUS.clone() - 3): {:?}\n\
             Rsa3x5::elem(-3): {:?}",
            Rsa3x5::elem(3),
            Rsa3x5::elem(RSA3X5_MODULUS.clone() - 3),
            Rsa3x5::elem(-3),
        );
        // TODO: Add a trickier coset test involving `op`.
    }

    #[test]
    fn test_rsa3x5_exp() {
        let a = Rsa3x5::exp(&Rsa3x5::elem(2), &int(3));
        assert!(a == Rsa3x5::elem(8));
        let b = Rsa3x5::exp(&Rsa3x5::elem(2), &int(4096));
        assert_eq!(b, Rsa3x5::elem(1));
        let c = Rsa3x5::exp(&Rsa3x5::elem(2), &RSA3X5_MODULUS);
        dbg!(c);
        let d = Rsa3x5::exp(&Rsa3x5::elem(2), &(RSA3X5_MODULUS.clone() * int(2)));
        dbg!(d);
    }

    #[test]
    fn test_rsa3x5_inv() {
        let x = Rsa3x5::elem(2);
        let inv = Rsa3x5::inv(&x);
        dbg!(inv.clone());
        assert_eq!(Rsa3x5::op(&x, &inv), Rsa3x5::id());
    }
}
