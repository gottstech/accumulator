//! RSA (100) group using GMP integers in the `rug` crate.
use super::{ElemFrom, Group, UnknownOrderGroup};
use crate::util::{int, TypeRep};
use rug::Integer;
use std::fmt;
use std::str::FromStr;

#[allow(clippy::module_name_repetitions)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
/// RSA-100 group implementation.
/// Modulus taken from [here](https://en.wikipedia.org/wiki/RSA_numbers#RSA-100).
/// **Note**: If you want to use `Rsa100` outside the context of this crate, be advised that
/// it treats `x` and `-x` as the same element for sound proofs-of-exponentiation.
/// See BBF (page 9).
pub enum Rsa100 {}

/// RSA-100 modulus.
const RSA100_MODULUS_DECIMAL: &str = "15226050279225333605356183781326374297180681149613\
     80688657908494580122963258952897654000350692006139
    ";

lazy_static! {
    static ref RSA100_MODULUS: Integer = Integer::from_str(RSA100_MODULUS_DECIMAL).unwrap();
    static ref HALF_MODULUS: Integer = RSA100_MODULUS.clone() / 2;
}

#[allow(clippy::module_name_repetitions)]
#[derive(Clone, PartialEq, Eq, Hash)]
/// An RSA 100 group element, directly wrapping a GMP integer from the `rug` crate.
pub struct Rsa100Elem(Integer);

// Manual debug impl required because want to be one line.
impl fmt::Debug for Rsa100Elem {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Rsa100Elem( {:?} )", self.0)
    }
}

impl TypeRep for Rsa100 {
    type Rep = Integer;
    fn rep() -> &'static Self::Rep {
        &RSA100_MODULUS
    }
}

impl Group for Rsa100 {
    type Elem = Rsa100Elem;

    fn id_(_: &Integer) -> Rsa100Elem {
        Self::elem(1)
    }

    fn op_(modulus: &Integer, a: &Rsa100Elem, b: &Rsa100Elem) -> Rsa100Elem {
        Self::elem(int(&a.0 * &b.0) % modulus)
    }

    fn exp_(modulus: &Integer, x: &Rsa100Elem, n: &Integer) -> Rsa100Elem {
        // A side-channel resistant impl is 40% slower; we'll consider it in the future if we need to.
        Self::elem(x.0.pow_mod_ref(n, modulus).unwrap())
    }

    fn inv_(modulus: &Integer, x: &Rsa100Elem) -> Rsa100Elem {
        Self::elem(x.0.invert_ref(modulus).unwrap())
    }
}

impl<T> ElemFrom<T> for Rsa100
where
    Integer: From<T>,
{
    fn elem(t: T) -> Rsa100Elem {
        let modulus = Self::rep();
        let val = int(t) % modulus;
        if val > *HALF_MODULUS {
            Rsa100Elem(<(Integer, Integer)>::from((-val).div_rem_euc_ref(&modulus)).1)
        } else {
            Rsa100Elem(val)
        }
    }
}

impl UnknownOrderGroup for Rsa100 {
    fn unknown_order_elem_(_: &Integer) -> Rsa100Elem {
        Self::elem(2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rsa100_init() {
        let x = Rsa100::rep();
        println!("Rsa100::rep() - {:?}", x);
    }

    #[test]
    fn test_rsa100_op() {
        let a = Rsa100::op(&Rsa100::elem(2), &Rsa100::elem(3));
        println!("Rsa100:  2*3= {:?}", a);
        assert_eq!(a, Rsa100::elem(6));

        let b = Rsa100::op(&Rsa100::elem(-2), &Rsa100::elem(-3));
        println!("Rsa100: -2*(-3)= {:?}", b);
        assert_eq!(b, Rsa100::elem(6));

        let c = Rsa100::op(&Rsa100::elem(-2), &Rsa100::elem(3));
        println!("Rsa100: -2*3={:?}", c);
        assert_eq!(c, Rsa100::elem(-6));
    }

    /// Tests that `-x` and `x` are treated as the same element.
    #[test]
    fn test_rsa100_cosets() {
        assert_eq!(Rsa100::elem(3), Rsa100::elem(RSA100_MODULUS.clone() - 3));
        println!(
            "Rsa100::elem(3): {:?} \n\
             Rsa100::elem(RSA100_MODULUS.clone() - 3): {:?}\n\
             Rsa100::elem(-3): {:?}",
            Rsa100::elem(3),
            Rsa100::elem(RSA100_MODULUS.clone() - 3),
            Rsa100::elem(-3),
        );
        // TODO: Add a trickier coset test involving `op`.
    }

    #[test]
    fn test_rsa100_exp() {
        let a = Rsa100::exp(&Rsa100::elem(2), &int(3));
        assert!(a == Rsa100::elem(8));
        let b = Rsa100::exp(&Rsa100::elem(2), &int(4096));
        assert_eq!(
            b,
            Rsa100::elem(
                Integer::parse(
                    "23853855460373584876145204365446184722448625612051\
                     4468122906297969891609543350919691048700369060845"
                )
                .unwrap()
            )
        );
        let c = Rsa100::exp(&Rsa100::elem(2), &RSA100_MODULUS);
        dbg!(c);
        let d = Rsa100::exp(&Rsa100::elem(2), &(RSA100_MODULUS.clone() * int(2)));
        dbg!(d);
    }

    #[test]
    fn test_rsa100_inv() {
        let x = Rsa100::elem(2);
        let inv = Rsa100::inv(&x);
        dbg!(inv.clone());
        assert_eq!(Rsa100::op(&x, &inv), Rsa100::id());
    }
}
