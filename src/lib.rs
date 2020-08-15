use simba::scalar::{ComplexField, RealField};
//use num_traits::cast::ToPrimitive;
use num_traits::{Num, NumCast,FromPrimitive, ToPrimitive};
#[derive(Copy, Clone, Debug)]
pub enum PolygammaError{
    NegInt,
    Misc
}

/// Cofficients in odd terms of asymptotic expansion starting at 1/x^3
static TRIGAMMA_ASYMPT_ODD : [f64; 5] = [
    1.0/6.0,
    -1.0/30.0,
    1.0/42.0,
    -1.0/30.0,
    5.0/66.0,
];

/// Evaluate \sum_{k=0}^n 1 / (z + k)^2
fn jump_sum<T: ComplexField>(z: T, n: i32) -> T{
    let mut s = T::zero();

    for k in 0..n+1{
        let x: T = (z + T::from_i32(k).unwrap()).recip();
        s += x*x;
    }

    s
}

/// Evaluate the asymptotic sum (y = 1/x)
fn asym_sum<T: Copy + Num + FromPrimitive>(y: T) -> T{
    let y2: T = y * y;
    let y3: T = y2 * y;
    let y4: T = y2 * y2;

    #[allow(non_snake_case)]
    let mut A : [T; 5] = [T::zero() ; 5];
    for (a, &b) in A.iter_mut().zip(TRIGAMMA_ASYMPT_ODD.iter()){
        *a = T::from_f64(b).unwrap();
    }
    //let A = TRIGAMMA_ASYMPT_ODD.iter().map(|&a| T::from_f64(a).unwrap()).co

    return
        y + T::from_f64(1.0 / 2.0).unwrap() * y2
          + y3 * (
                A[0] +    A[1] * y2   +    A[2] * y4
                     + A[3] * y3 * y3 + A[4] * y4 * y4
            )
}

pub fn trigamma<T: ComplexField>( z: T) -> Result<T, PolygammaError>
where T::RealField : RealField + FromPrimitive + NumCast
{
    let pi = T::from_real(T::RealField::pi());
    let yn = 5;
    let y : T::RealField = FromPrimitive::from_i32(yn).unwrap();

    let re_z = z.real();
    // Check that re_z is not an negative integer and reflect if negative
    if re_z.is_sign_negative() {
        let x : T = pi * ( (pi * z).sin().recip());
        if !x.is_finite(){
            return Err(PolygammaError::NegInt);
        }

        return trigamma( T::one() - z).map(|psi2| x*x - psi2);
    }
    // For small re_z, use the recurrence relation to evaluate for larger re_z
    if re_z < y {
        let dy : T::RealField = y - re_z;
        let yn: i32 = dy.ceil().to_i32().unwrap();

        return trigamma( z + T::from_i32(yn).unwrap() + T::one())
            .map(|psi2| psi2 + jump_sum(z, yn))
    }

    //For large re_z, use the asymptotic series
    let w = z.recip();
    let psi = asym_sum(w);

    Ok(psi)

}

#[cfg(test)]
mod tests {
    use core::f64::consts::PI;
    use super::trigamma;
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn psi_onehalf(){
        let psi = trigamma(0.5).unwrap();
        let y = PI * PI / 2.0;
        println!("psi(1/2) = {}\npi^2/2 = {}", psi, y);
    }
    #[test]
    fn psi_six(){
        let psi = trigamma(6.0).unwrap();
        let y = 0.18132295573711532536;
        println!("computer psi(6) = {}\nexpected psi(6) = {}", psi, y);
    }

    #[test]
    fn psi_three(){
        let psi = trigamma(3.0).unwrap();
        let y = 0.39493406684822643647;
        println!("computer psi(3) = {}\nexpected psi(3) = {}", psi, y);
    }
}
